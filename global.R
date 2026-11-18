## global.R -- standalone design functions for the snpGuide Shiny app
##
## GENERATED FROM THE snpguideR PACKAGE, version 0.8.0.
## Do not edit by hand: change the package sources and regenerate, otherwise
## the deployed app and the package will drift apart.
##
## Source of record: https://github.com/ybaglaenko/snpguideR
## Data source: NCBI dbSNP + E-utilities (the Ensembl REST API this app used
## previously is frozen and unreliable).

library(httr)
library(jsonlite)
library(stringr)
library(dplyr)
library(purrr)
library(tibble)
library(rlang)
library(methods)
library(Biostrings)

.data <- rlang::.data

## NCBI dbSNP and E-utilities helpers -----------------------------------------

#' Strip a leading "chr" from a chromosome name
#'
#' @param x Chromosome name, for example `"chr10"` or `"10"`.
#'
#' @return A character vector with any leading `"chr"` removed.
#' @noRd
normalize_chr <- function(x) {
  sub("^chr", "", as.character(x), ignore.case = TRUE)
}

#' Translate a human RefSeq chromosome accession to a chromosome name
#'
#' @param accession RefSeq accession such as `"NC_000010.11"`.
#'
#' @return Chromosome name such as `"10"`, `"X"` or `"MT"`, or `NA` when the
#'   accession is not a human chromosome.
#' @noRd
refseq_to_chr <- function(accession) {
  acc <- sub("\\..*$", "", as.character(accession))
  if (identical(acc, "NC_012920")) return("MT")

  m <- regmatches(acc, regexpr("^NC_0000([0-9]{2})$", acc))
  if (length(m) == 0L) return(NA_character_)

  n <- as.integer(sub("^NC_0000", "", m))
  if (is.na(n)) return(NA_character_)
  if (n == 23L) return("X")
  if (n == 24L) return("Y")
  if (n >= 1L && n <= 22L) return(as.character(n))
  NA_character_
}

#' Issue a GET request to NCBI with retries
#'
#' @param url Full request URL.
#' @param accept Media type passed to [httr::accept()].
#' @param timeout Seconds allowed for the whole request.
#' @param connect_timeout Seconds allowed to establish the connection. curl
#'   gives up after 10 by default, which is what produces "Connection timed out
#'   after 10002 milliseconds" against an unresponsive host.
#' @param retries Number of attempts made for transient failures.
#'
#' @return An [httr::response], or a character scalar holding the error message
#'   when no request could be made.
#' @noRd
ncbi_get <- function(url, accept, timeout, connect_timeout, retries) {
  tryCatch(
    httr::RETRY(
      "GET", url,
      httr::accept(accept),
      httr::timeout(timeout),
      httr::config(connecttimeout = connect_timeout),
      httr::user_agent("snpguideR"),
      times = retries,
      pause_base = 1,
      pause_cap = 8,
      terminate_on = c(400, 404),
      quiet = TRUE
    ),
    error = function(e) conditionMessage(e)
  )
}

#' Explain a failed NCBI request
#'
#' @param msg The curl error message.
#' @noRd
ncbi_failure_note <- function(msg) {
  if (grepl("[Tt]imed out|[Tt]imeout", msg)) {
    paste0(
      "Could not reach NCBI: ", msg,
      "\n  The connection was never established, so this is the network path ",
      "or the service, not your query.",
      "\n  Check with ncbi_available(). If your network is slow, raise ",
      "connect_timeout."
    )
  } else {
    paste0("Could not reach NCBI: ", msg)
  }
}

#' Check whether the NCBI services are responding
#'
#' Calls the E-utilities information endpoint. Use this to tell an outage apart
#' from a problem with a particular variant or query.
#'
#' @param timeout Numeric scalar. Seconds allowed for the request.
#' @param connect_timeout Numeric scalar. Seconds allowed to establish the
#'   connection.
#'
#' @return `TRUE` when NCBI answers, otherwise `FALSE`. Never raises an error.
#'
#' @examples
#' \donttest{
#' ncbi_available()
#' }
#' @export
ncbi_available <- function(timeout = 10, connect_timeout = 10) {
  res <- ncbi_get(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi?db=snp&retmode=json",
    "application/json", timeout, connect_timeout, retries = 1
  )
  if (is.character(res)) return(FALSE)
  httr::status_code(res) == 200
}

#' Look up a variant in dbSNP
#'
#' Retrieves the genomic placement and alleles of a reference SNP cluster
#' identifier from the NCBI dbSNP refsnp service.
#'
#' The function requires an internet connection. If dbSNP is unreachable, the
#' identifier is unknown, or the record carries no placement on the requested
#' assembly, a warning is issued and `NULL` is returned rather than failing.
#'
#' Positions reported by dbSNP use SPDI, which is 0-based. The `position`
#' returned here is converted to the 1-based convention used elsewhere in this
#' package and by genome browsers.
#'
#' @param rsid Character scalar. Reference SNP cluster identifier, with or
#'   without the `"rs"` prefix, for example `"rs61839660"` or `"61839660"`.
#' @param assembly Character scalar. Assembly to take the placement from,
#'   `"GRCh38"` or `"GRCh37"`.
#' @param timeout Numeric scalar. Request timeout in seconds.
#' @param connect_timeout Numeric scalar. Seconds allowed to establish the
#'   connection before giving up.
#' @param retries Integer scalar. Number of attempts made for transient
#'   failures, with exponential backoff. HTTP 400 and 404 are treated as final.
#' @param server Character scalar. Base URL of the NCBI variation service.
#'
#' @return A list with elements `name` (the canonical rsID), `requested`,
#'   `chromosome`, `accession` (the RefSeq sequence the position refers to),
#'   `position` (1-based), `assembly`, `ref`, `alleles` and `allele_string`
#'   (reference allele first, separated by `"/"`), or `NULL` on failure.
#'
#' @seealso [ncbi_region_sequence()] for the companion sequence call.
#'
#' @examples
#' \donttest{
#' dbsnp_variant("rs61839660")
#' }
#' @export
dbsnp_variant <- function(rsid,
                          assembly = c("GRCh38", "GRCh37"),
                          timeout = 20,
                          connect_timeout = 30,
                          retries = 3,
                          server = "https://api.ncbi.nlm.nih.gov") {
  rsid <- as.character(rsid)
  if (length(rsid) != 1L || is.na(rsid) || !nzchar(rsid)) {
    stop("`rsid` must be a single non-empty character string.", call. = FALSE)
  }
  assembly <- match.arg(assembly)

  num <- sub("^rs", "", trimws(rsid), ignore.case = TRUE)
  if (!grepl("^[0-9]+$", num)) {
    warning("'", rsid, "' is not a valid rsID.", call. = FALSE)
    return(NULL)
  }

  res <- ncbi_get(paste0(server, "/variation/v0/refsnp/", num),
                  "application/json", timeout, connect_timeout, retries)
  if (is.character(res)) {
    warning(ncbi_failure_note(res), call. = FALSE)
    return(NULL)
  }
  if (httr::status_code(res) != 200) {
    warning("dbSNP returned status ", httr::status_code(res), " for 'rs",
            num, "'.", call. = FALSE)
    return(NULL)
  }

  dat <- tryCatch(
    jsonlite::fromJSON(rawToChar(res$content), simplifyVector = FALSE),
    error = function(e) NULL
  )
  if (is.null(dat)) {
    warning("Could not parse the dbSNP response for 'rs", num, "'.",
            call. = FALSE)
    return(NULL)
  }

  canonical <- paste0("rs", if (is.null(dat$refsnp_id)) num else dat$refsnp_id)
  if (!identical(canonical, paste0("rs", num))) {
    message("rs", num, " has been merged into ", canonical, ".")
  }

  placements <- dat$primary_snapshot_data$placements_with_allele
  if (is.null(placements) || length(placements) == 0L) {
    warning(canonical, " has no placement data; it may be withdrawn or ",
            "unplaced.", call. = FALSE)
    return(NULL)
  }

  pick <- NULL
  for (pl in placements) {
    traits <- pl$placement_annot$seq_id_traits_by_assembly
    if (is.null(traits) || length(traits) == 0L) next
    asm <- traits[[1]]$assembly_name
    if (is.null(asm) || !startsWith(asm, assembly)) next
    if (!isTRUE(traits[[1]]$is_chromosome)) next
    pick <- pl
    pick_assembly <- asm
    break
  }

  if (is.null(pick)) {
    warning(canonical, " has no chromosomal placement on ", assembly, ".",
            call. = FALSE)
    return(NULL)
  }

  spdis <- lapply(pick$alleles, function(a) a$allele$spdi)
  spdis <- Filter(Negate(is.null), spdis)
  if (length(spdis) == 0L) {
    warning("No alleles reported for ", canonical, ".", call. = FALSE)
    return(NULL)
  }

  ref <- spdis[[1]]$deleted_sequence
  pos <- as.integer(spdis[[1]]$position) + 1L
  accession <- as.character(pick$seq_id)

  alts <- unique(vapply(spdis, function(s) as.character(s$inserted_sequence),
                        character(1)))
  alts <- setdiff(alts, ref)
  alleles <- c(ref, alts)

  list(
    name          = canonical,
    requested     = rsid,
    chromosome    = refseq_to_chr(accession),
    accession     = accession,
    position      = pos,
    assembly      = pick_assembly,
    ref           = ref,
    alleles       = alleles,
    allele_string = paste(alleles, collapse = "/")
  )
}

#' Retrieve a genomic region sequence from NCBI
#'
#' Fetches the reference sequence of an interval from an NCBI nucleotide record
#' using E-utilities.
#'
#' The function requires an internet connection and returns `NULL` with a
#' warning on failure.
#'
#' NCBI asks that unauthenticated clients make no more than three requests per
#' second. Supply an API key, or set the `NCBI_API_KEY` environment variable,
#' to raise that limit.
#'
#' @param accession Character scalar. Nucleotide accession, for example
#'   `"NC_000010.11"`.
#' @param start,end Integer scalars giving the 1-based inclusive interval.
#' @param strand Integer scalar, `1` for the plus strand and `2` for the minus
#'   strand, following the E-utilities convention.
#' @param api_key Character scalar. NCBI API key, taken from the `NCBI_API_KEY`
#'   environment variable by default. Use `""` for none.
#' @inheritParams dbsnp_variant
#'
#' @return A [Biostrings::DNAString], or `NULL` on failure.
#'
#' @examples
#' \donttest{
#' ncbi_region_sequence("NC_000010.11", 6052714, 6052754)
#' }
#' @export
ncbi_region_sequence <- function(accession, start, end, strand = 1,
                                 timeout = 20,
                                 connect_timeout = 30,
                                 retries = 3,
                                 api_key = Sys.getenv("NCBI_API_KEY")) {
  accession <- as.character(accession)
  start <- as.integer(start)
  end <- as.integer(end)

  if (length(accession) != 1L || !nzchar(accession)) {
    warning("`accession` must be a single non-empty character string.",
            call. = FALSE)
    return(NULL)
  }
  if (is.na(start) || is.na(end) || start < 1L || end < 1L || end < start) {
    warning("Invalid genomic interval requested.", call. = FALSE)
    return(NULL)
  }

  url <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
    "?db=nuccore&id=", accession,
    "&rettype=fasta&retmode=text",
    "&seq_start=", start, "&seq_stop=", end,
    "&strand=", as.integer(strand),
    "&tool=snpguideR"
  )
  if (nzchar(api_key)) url <- paste0(url, "&api_key=", api_key)

  res <- ncbi_get(url, "text/plain", timeout, connect_timeout, retries)
  if (is.character(res)) {
    warning(ncbi_failure_note(res), call. = FALSE)
    return(NULL)
  }
  if (httr::status_code(res) != 200) {
    warning("NCBI returned status ", httr::status_code(res), " for ",
            accession, ":", start, "-", end, ".", call. = FALSE)
    return(NULL)
  }

  txt <- httr::content(res, as = "text", encoding = "UTF-8")
  lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  lines <- lines[!grepl("^>", lines)]
  seq_txt <- toupper(gsub("[^A-Za-z]", "", paste(lines, collapse = "")))

  if (!nzchar(seq_txt)) {
    warning("NCBI returned an empty sequence for ", accession, ":", start,
            "-", end, ".", call. = FALSE)
    return(NULL)
  }

  tryCatch(Biostrings::DNAString(seq_txt), error = function(e) {
    warning("NCBI returned a sequence that is not DNA: ",
            conditionMessage(e), call. = FALSE)
    NULL
  })
}

## Editing window annotation --------------------------------------------------

#' Annotate the base editing window of a protospacer
#'
#' Extracts the editing window from a protospacer and records the positions of
#' the bases a deaminase could act on.
#'
#' @param protospacer Character scalar giving the protospacer sequence.
#' @param min_w,max_w Integer scalars giving the first and last protospacer
#'   position of the base editing window, counted from the
#'   protospacer-adjacent-motif-distal end.
#' @param bases Character vector of bases considered editable, typically
#'   `"A"` for adenine base editors and `"C"` for cytosine base editors.
#'
#' @return A one-row [tibble::tibble] with columns `editing_window`,
#'   `editable_bases`, `n_editable`, `editable_pos`, `editable_pos_by_base`,
#'   `n_A` and `n_C`. Positions are 1-based within the protospacer.
#'
#' @examples
#' annotate_window("ACGTACGTACGTACGTACGT", min_w = 4, max_w = 8)
#' @export
annotate_window <- function(protospacer, min_w, max_w, bases = c("A", "C")) {
  w <- stringr::str_sub(protospacer, start = min_w, end = max_w)

  pos_list <- purrr::map(bases, function(b) {
    hits <- gregexpr(b, w, fixed = TRUE)[[1]]
    if (length(hits) == 1L && hits[1] == -1L) {
      integer(0)
    } else {
      as.integer(hits) + (as.integer(min_w) - 1L)
    }
  })
  names(pos_list) <- bases

  base_counts <- purrr::map_int(pos_list, length)

  tibble::tibble(
    editing_window = w,
    editable_bases = paste(bases, collapse = ","),
    n_editable     = sum(base_counts),
    editable_pos   = paste(unlist(pos_list), collapse = ","),
    editable_pos_by_base = paste(
      purrr::map2_chr(names(pos_list), pos_list,
                      function(b, p) paste0(b, ":", paste(p, collapse = "|"))),
      collapse = ";"
    ),
    n_A = if ("A" %in% names(base_counts)) base_counts[["A"]] else 0L,
    n_C = if ("C" %in% names(base_counts)) base_counts[["C"]] else 0L
  )
}

## Protospacer design ---------------------------------------------------------

#' Match a protospacer-plus-motif pattern against a sequence
#'
#' @param subject A [Biostrings::DNAString].
#' @param guide_length Integer protospacer length.
#' @param pam Character protospacer adjacent motif, possibly containing IUPAC
#'   ambiguity codes.
#'
#' @return A [tibble::tibble] with columns `start`, `end` and `seq`, one row
#'   per match (matches may overlap), or `NULL` when there is no match.
#' @noRd
match_protospacers <- function(subject, guide_length, pam) {
  pattern <- Biostrings::DNAString(
    paste0(strrep("N", guide_length), pam)
  )
  hits <- Biostrings::matchPattern(pattern, subject, fixed = FALSE)
  if (length(hits) == 0L) return(NULL)

  df <- as.data.frame(hits)
  tibble::tibble(
    start = as.integer(df$start),
    end   = as.integer(df$end),
    seq   = as.character(df$seq)
  )
}

#' Design base editing protospacers from a supplied sequence
#'
#' Identifies protospacers in a reference sequence that place a variant inside
#' the base editing window. This is the computational core of [snpGuide()] and
#' runs offline, which makes it useful for testing and for variants that are
#' not in Ensembl.
#'
#' Two modes are available.
#'
#' * Transition mode (`force_protospacers = FALSE`) only returns guides for
#'   variants whose alleles contain a C/T or A/G transition, that is a change a
#'   cytosine or adenine base editor could install. Both the reference and each
#'   alternative allele sequence are searched, so guides are reported for
#'   installing and for reverting the variant.
#' * Forced mode (`force_protospacers = TRUE`) ignores the alleles and returns
#'   every protospacer that places the variant inside the editing window. The
#'   deaminase class is then assigned from the bases present in the window:
#'   a cytosine base editor when the window contains C, an adenine base editor
#'   when it contains A.
#'
#' @param sequence Reference sequence covering the variant, either a character
#'   string or a [Biostrings::DNAString].
#' @param snp_position Integer scalar. 1-based position of the variant within
#'   `sequence`.
#' @param allele_string Character scalar of alleles separated by `"/"`, as used
#'   by Ensembl, for example `"C/T"`. Required in transition mode and optional
#'   in forced mode.
#' @param snp_id Character scalar used to label the output, typically an rsID.
#' @param pam Character scalar giving the protospacer adjacent motif. IUPAC
#'   ambiguity codes are allowed, for example `"NGG"` or `"NGN"`.
#' @param guide_length Integer scalar. Protospacer length in nucleotides.
#' @param min_editing_window,max_editing_window Integer scalars giving the
#'   first and last protospacer position of the base editing window.
#' @param force_protospacers Logical scalar selecting forced mode.
#' @param force_edit_bases Character vector of bases treated as editable in
#'   forced mode.
#' @param search_both_strands Logical scalar. In forced mode, also search the
#'   reverse complement strand.
#' @param keep_only_window_has_editable Logical scalar. In forced mode, discard
#'   protospacers whose editing window contains no editable base.
#'
#' @return A [tibble::tibble] of candidate guides, or `NULL` when none are
#'   found. Transition mode returns the columns `snp`, `allele_string`, `mode`,
#'   `strand`, `protospacer`, `pam`, `snp_position`, `targeted_nuc`,
#'   `targeted_allele`, `enzyme`, `conversion`, `editing_window` and
#'   `bystander_likely`. Forced mode returns `snp`, `allele_string`, `mode`,
#'   `strand`, `protospacer`, `pam`, `snp_pos_in_protospacer`, `targeted_nuc`,
#'   `enzyme`, `conversion`, `editing_window`, `editable_bases`, `n_A`, `n_C`,
#'   `n_editable`, `editable_pos` and `editable_pos_by_base`.
#'
#' @seealso [snpGuide()] to design guides directly from an rsID.
#'
#' @examples
#' seq <- paste0(
#'   "TTGACCTAGGACCTTGGACA", "C", "GGTACCGGTTACCGGAACCT"
#' )
#' snpGuide_from_sequence(
#'   sequence = seq,
#'   snp_position = 21,
#'   allele_string = "C/T",
#'   snp_id = "example"
#' )
#' @export
snpGuide_from_sequence <- function(sequence,
                                   snp_position,
                                   allele_string = NA_character_,
                                   snp_id = NA_character_,
                                   pam = "NGN",
                                   guide_length = 20L,
                                   min_editing_window = 4L,
                                   max_editing_window = 8L,
                                   force_protospacers = FALSE,
                                   force_edit_bases = c("A", "C"),
                                   search_both_strands = TRUE,
                                   keep_only_window_has_editable = TRUE) {

  guide_length <- as.integer(guide_length)
  min_editing_window <- as.integer(min_editing_window)
  max_editing_window <- as.integer(max_editing_window)
  snp_position <- as.integer(snp_position)
  pam <- toupper(as.character(pam))

  if (is.na(min_editing_window) || is.na(max_editing_window) ||
      min_editing_window < 1L ||
      max_editing_window > guide_length ||
      min_editing_window > max_editing_window) {
    stop("Invalid editing window. Ensure 1 <= min <= max <= guide_length.",
         call. = FALSE)
  }

  reference_seq <- if (methods::is(sequence, "DNAString")) {
    sequence
  } else {
    Biostrings::DNAString(toupper(as.character(sequence)))
  }

  seq_len_nt <- length(reference_seq)
  if (is.na(snp_position) || snp_position < 1L || snp_position > seq_len_nt) {
    stop("`snp_position` must lie inside `sequence`.", call. = FALSE)
  }

  pam_width <- nchar(pam)
  rc_index <- function(i) seq_len_nt - i + 1L

  ## ---- forced mode ---------------------------------------------------------
  if (isTRUE(force_protospacers)) {

    build_forced_df <- function(subject, strand_label, snp_idx) {
      hits <- match_protospacers(subject, guide_length, pam)
      if (is.null(hits)) return(NULL)

      df <- dplyr::mutate(
        hits,
        strand = strand_label,
        protospacer = stringr::str_sub(.data$seq, 1L, guide_length),
        pam = stringr::str_sub(.data$seq, guide_length + 1L,
                               guide_length + pam_width),
        snp_pos_in_protospacer = snp_idx - .data$start + 1L
      )
      df <- dplyr::filter(
        df,
        .data$snp_pos_in_protospacer >= 1L,
        .data$snp_pos_in_protospacer <= guide_length,
        .data$snp_pos_in_protospacer >= min_editing_window,
        .data$snp_pos_in_protospacer <= max_editing_window
      )
      df <- dplyr::select(df, -"seq", -"start", -"end")

      if (nrow(df) == 0L) return(NULL)
      df
    }

    df_fwd <- build_forced_df(reference_seq, "forward", snp_position)

    df_rev <- NULL
    if (isTRUE(search_both_strands)) {
      df_rev <- build_forced_df(
        Biostrings::reverseComplement(reference_seq),
        "reverse",
        rc_index(snp_position)
      )
    }

    guides <- dplyr::bind_rows(df_fwd, df_rev)
    if (is.null(guides) || nrow(guides) == 0L) return(NULL)

    window_ann <- dplyr::bind_rows(purrr::map(
      guides$protospacer,
      function(p) annotate_window(
        protospacer = p,
        min_w = min_editing_window,
        max_w = max_editing_window,
        bases = force_edit_bases
      )
    ))

    guides <- dplyr::bind_cols(guides, window_ann)
    guides <- dplyr::mutate(
      guides,
      conversion = dplyr::case_when(
        .data$n_A > 0L & .data$n_C > 0L ~ "C>T;A>G",
        .data$n_C > 0L                  ~ "C>T",
        .data$n_A > 0L                  ~ "A>G",
        TRUE                            ~ NA_character_
      ),
      enzyme = dplyr::case_when(
        .data$n_A > 0L & .data$n_C > 0L ~ "CBE;ABE",
        .data$n_C > 0L                  ~ "CBE",
        .data$n_A > 0L                  ~ "ABE",
        TRUE                            ~ NA_character_
      ),
      snp = as.character(snp_id),
      allele_string = as.character(allele_string),
      mode = "forced_overlap_snp_in_window",
      targeted_nuc = stringr::str_sub(.data$protospacer,
                                      .data$snp_pos_in_protospacer,
                                      .data$snp_pos_in_protospacer)
    )

    if (isTRUE(keep_only_window_has_editable)) {
      guides <- dplyr::filter(guides, .data$n_editable > 0L)
    }
    if (nrow(guides) == 0L) return(NULL)

    guides <- dplyr::select(
      guides,
      "snp", "allele_string", "mode", "strand",
      "protospacer", "pam",
      "snp_pos_in_protospacer", "targeted_nuc",
      "enzyme", "conversion",
      "editing_window", "editable_bases", "n_A", "n_C", "n_editable",
      "editable_pos", "editable_pos_by_base"
    )

    return(guides)
  }

  ## ---- transition mode -----------------------------------------------------
  if (is.na(allele_string) || !nzchar(allele_string)) {
    stop("`allele_string` is required unless `force_protospacers = TRUE`.",
         call. = FALSE)
  }

  alleles <- unlist(stringr::str_split(allele_string, pattern = "/"))
  alleles <- toupper(trimws(alleles))
  if (length(alleles) < 2L) return(NULL)

  guides_design_func <- function(subject_seq, allele1, allele2, label) {

    if ((allele1 == "C" && allele2 == "T") ||
        (allele1 == "A" && allele2 == "G")) {
      subject <- subject_seq
      strand_label <- "forward"
      snp_idx <- snp_position
      conversion <- if (allele1 == "C") "C>T" else "A>G"

    } else if ((allele1 == "T" && allele2 == "C") ||
               (allele1 == "G" && allele2 == "A")) {
      subject <- Biostrings::reverseComplement(subject_seq)
      strand_label <- "reverse"
      snp_idx <- rc_index(snp_position)
      conversion <- if (allele1 == "T") "A>G" else "C>T"

    } else {
      return(NULL)
    }

    hits <- match_protospacers(subject, guide_length, pam)
    if (is.null(hits)) return(NULL)

    df <- dplyr::mutate(hits, snp_position = snp_idx - .data$start + 1L)
    df <- dplyr::filter(
      df,
      .data$snp_position >= min_editing_window,
      .data$snp_position <= max_editing_window
    )
    if (nrow(df) == 0L) return(NULL)

    df <- dplyr::mutate(
      df,
      strand = strand_label,
      targeted_allele = label,
      conversion = conversion,
      enzyme = dplyr::case_when(
        conversion == "C>T" ~ "CBE",
        conversion == "A>G" ~ "ABE",
        TRUE                ~ NA_character_
      ),
      protospacer = stringr::str_sub(.data$seq, 1L, guide_length),
      pam = stringr::str_sub(.data$seq, guide_length + 1L,
                             guide_length + pam_width),
      targeted_nuc = stringr::str_sub(.data$seq, .data$snp_position,
                                      .data$snp_position)
    )

    dplyr::select(df, -"seq", -"start", -"end")
  }

  if (!(all(c("C", "T") %in% alleles) || all(c("A", "G") %in% alleles))) {
    return(NULL)
  }

  guides <- guides_design_func(
    reference_seq, alleles[1], alleles[2],
    paste0("reference_", alleles[1])
  )

  seq_str <- as.character(reference_seq)

  for (i in 2:length(alleles)) {
    allele2 <- alleles[i]
    if (nchar(allele2) != 1L || !allele2 %in% c("A", "C", "G", "T")) next

    seq_alt <- Biostrings::DNAString(paste0(
      stringr::str_sub(seq_str, 1L, snp_position - 1L),
      allele2,
      stringr::str_sub(seq_str, snp_position + 1L, nchar(seq_str))
    ))

    guide_alt <- guides_design_func(
      seq_alt, allele2, alleles[1], paste0("alternative_", allele2)
    )
    guides <- dplyr::bind_rows(guides, guide_alt)
  }

  if (is.null(guides) || nrow(guides) == 0L) return(NULL)

  guides <- dplyr::mutate(
    guides,
    snp = as.character(snp_id),
    allele_string = as.character(allele_string),
    mode = "editable_transition",
    editing_window = stringr::str_sub(.data$protospacer,
                                      min_editing_window,
                                      max_editing_window),
    bystander_likely = stringr::str_count(.data$editing_window,
                                          pattern = .data$targeted_nuc) > 1L
  )

  dplyr::select(
    guides,
    "snp", "allele_string", "mode", "strand",
    "protospacer", "pam",
    "snp_position", "targeted_nuc", "targeted_allele",
    "enzyme", "conversion",
    "editing_window", "bystander_likely"
  )
}

#' Design base editing protospacers for a variant
#'
#' Retrieves the genomic position, alleles and flanking reference sequence of a
#' variant from NCBI, then designs base editing protospacers that place the
#' variant inside the editing window.
#'
#' The variant is resolved with [dbsnp_variant()] and the flanking sequence is
#' taken from the same RefSeq record that dbSNP placed it on, with
#' [ncbi_region_sequence()], so the coordinates and the sequence cannot
#' disagree about the assembly. The function requires an internet connection;
#' if the variant or its sequence cannot be retrieved a warning is issued and
#' `NULL` is returned. See [snpGuide_from_sequence()] for the design rules and
#' the two modes.
#'
#' @param rsid Character scalar. Reference SNP cluster identifier, for example
#'   `"rs61839660"`.
#' @param flank Integer scalar. Number of nucleotides retrieved either side of
#'   the variant.
#' @param assembly Character scalar. Assembly to design against, `"GRCh38"` or
#'   `"GRCh37"`.
#' @inheritParams snpGuide_from_sequence
#' @param ... Further arguments passed to [dbsnp_variant()] and
#'   [ncbi_region_sequence()], such as `timeout`, `connect_timeout`, `retries`
#'   or `api_key`.
#'
#' @return A [tibble::tibble] of candidate guides as described in
#'   [snpGuide_from_sequence()], or `NULL` when none are found.
#'
#' @seealso [snpGuide_batch()] for several variants,
#'   [run_snpGuide_app()] for the interactive application.
#'
#' @examples
#' \donttest{
#' # requires an internet connection
#' snpGuide("rs61839660", pam = "NGN")
#'
#' snpGuide("rs61839660", force_protospacers = TRUE)
#'
#' # design against GRCh37 coordinates instead
#' snpGuide("rs61839660", assembly = "GRCh37")
#' }
#' @export
snpGuide <- function(rsid,
                     pam = "NGN",
                     guide_length = 20L,
                     min_editing_window = 4L,
                     max_editing_window = 8L,
                     flank = 20L,
                     force_protospacers = FALSE,
                     force_edit_bases = c("A", "C"),
                     search_both_strands = TRUE,
                     keep_only_window_has_editable = TRUE,
                     assembly = c("GRCh38", "GRCh37"),
                     ...) {

  rsid <- as.character(rsid)
  if (length(rsid) != 1L || is.na(rsid) || !nzchar(rsid)) return(NULL)

  assembly <- match.arg(assembly)
  flank <- as.integer(flank)

  variant <- dbsnp_variant(rsid, assembly = assembly, ...)
  if (is.null(variant)) return(NULL)

  start <- variant$position - flank
  end <- variant$position + flank
  if (start < 1L) return(NULL)

  reference_seq <- ncbi_region_sequence(
    accession = variant$accession,
    start = start,
    end = end,
    strand = 1,
    ...
  )
  if (is.null(reference_seq)) return(NULL)

  if (length(reference_seq) != (end - start + 1L)) {
    warning("NCBI returned ", length(reference_seq), " nucleotides for ",
            variant$name, " where ", end - start + 1L, " were requested.",
            call. = FALSE)
    return(NULL)
  }

  snpGuide_from_sequence(
    sequence = reference_seq,
    snp_position = flank + 1L,
    allele_string = variant$allele_string,
    snp_id = variant$name,
    pam = pam,
    guide_length = guide_length,
    min_editing_window = min_editing_window,
    max_editing_window = max_editing_window,
    force_protospacers = force_protospacers,
    force_edit_bases = force_edit_bases,
    search_both_strands = search_both_strands,
    keep_only_window_has_editable = keep_only_window_has_editable
  )
}

## Batch design --------------------------------------------------------------

#' Design base editing protospacers for several variants
#'
#' Applies [snpGuide()] to a vector of rsIDs and row-binds the results.
#' Variants that return no guides, or that cannot be retrieved from Ensembl,
#' are skipped rather than aborting the run.
#'
#' @param rsids Character vector of reference SNP cluster identifiers.
#'   Duplicated and empty entries are dropped.
#' @param verbose Logical scalar. Report progress with [message()].
#' @param delay Numeric scalar. Seconds to pause between variants. NCBI asks
#'   unauthenticated clients to stay under three requests per second, and each
#'   variant costs two requests. Set to `0` if you supply an API key.
#' @param ... Further arguments passed to [snpGuide()].
#'
#' @return A [tibble::tibble] of guides for all variants, with zero rows when
#'   nothing was found.
#'
#' @examples
#' \donttest{
#' # requires an internet connection
#' snpGuide_batch(c("rs12722517", "rs7920946"))
#' }
#' @export
snpGuide_batch <- function(rsids, verbose = FALSE, delay = 0.34, ...) {
  rsids <- unique(trimws(as.character(rsids)))
  rsids <- rsids[!is.na(rsids) & nzchar(rsids)]
  if (length(rsids) == 0L) return(tibble::tibble())

  safe_snpGuide <- purrr::possibly(snpGuide, otherwise = NULL)

  out <- purrr::map(seq_along(rsids), function(i) {
    if (isTRUE(verbose)) {
      message("Processing ", rsids[i], " (", i, "/", length(rsids), ")")
    }
    out_i <- safe_snpGuide(rsid = rsids[i], ...)
    if (delay > 0 && i < length(rsids)) Sys.sleep(delay)
    out_i
  })

  dplyr::bind_rows(out)
}

#' Example rsIDs shipped with the package
#'
#' Returns the rsIDs used as the example batch input of the 'shiny'
#' application.
#'
#' @return A character vector of rsIDs.
#'
#' @examples
#' example_rsids()
#' @export
example_rsids <- function() {
  path <- system.file("extdata", "example_rsids.txt", package = "snpguideR")
  if (!nzchar(path)) return(character(0))
  readLines(path)
}
