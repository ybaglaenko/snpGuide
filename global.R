## global.R

## Load libraries
library(httr)
library(jsonlite)
library(stringr)
library(dplyr)
library(purrr)
library(Biostrings)
library(shiny)
library(DT)

# -----------------------------
# Ensembl REST helpers (no BSgenome)
# -----------------------------
normalize_chr <- function(x) {
  x <- as.character(x)
  sub("^chr", "", x, ignore.case = TRUE)
}

get_variant_mapping <- function(rsid) {
  url <- paste0(
    "https://rest.ensembl.org/variation/homo_sapiens/",
    rsid,
    "?content-type=application/json"
  )

  res <- httr::GET(
    url,
    httr::accept("application/json"),
    httr::timeout(20)
  )

  if (httr::status_code(res) != 200) return(NULL)

  dat <- tryCatch(jsonlite::fromJSON(rawToChar(res$content)), error = function(e) NULL)
  if (is.null(dat) || is.null(dat$mappings) || nrow(dat$mappings) < 1) return(NULL)

  dat
}

get_region_sequence <- function(chr, start, end, strand = 1) {
  chr <- normalize_chr(chr)
  start <- as.integer(start)
  end <- as.integer(end)
  strand <- as.integer(strand)

  if (is.na(start) || is.na(end) || start < 1 || end < 1 || end < start) return(NULL)

  region <- paste0(chr, ":", start, "..", end, ":", strand)
  url <- paste0(
    "https://rest.ensembl.org/sequence/region/human/",
    region,
    "?content-type=text/plain"
  )

  res <- httr::GET(
    url,
    httr::accept("text/plain"),
    httr::timeout(20)
  )

  if (httr::status_code(res) != 200) return(NULL)

  seq_txt <- httr::content(res, as = "text", encoding = "UTF-8")
  seq_txt <- gsub("\\s+", "", seq_txt)

  if (!nzchar(seq_txt)) return(NULL)

  Biostrings::DNAString(seq_txt)
}

# -----------------------------
# Helper: annotate editing window
# -----------------------------
annotate_window <- function(protospacer, min_w, max_w, bases = c("A", "C")) {
  w <- str_sub(protospacer, start = min_w, end = max_w)

  pos_list <- purrr::map(bases, function(b) {
    hits <- gregexpr(b, w, fixed = TRUE)[[1]]
    if (length(hits) == 1 && hits[1] == -1) integer(0) else (hits + (min_w - 1))
  })
  names(pos_list) <- bases

  base_counts <- purrr::map_int(pos_list, length)

  tibble(
    editing_window = w,
    editable_bases = paste(bases, collapse = ","),
    n_editable     = sum(base_counts),
    editable_pos   = paste(unlist(pos_list), collapse = ","),
    editable_pos_by_base = paste(
      purrr::map2_chr(names(pos_list), pos_list, ~ paste0(.x, ":", paste(.y, collapse = "|"))),
      collapse = ";"
    ),
    n_A = if ("A" %in% names(base_counts)) base_counts[["A"]] else 0L,
    n_C = if ("C" %in% names(base_counts)) base_counts[["C"]] else 0L
  )
}

# -----------------------------
# Main function: snpGuide
# -----------------------------
snpGuide <- function(rsid,
                     pam = "NGN",
                     guide_length = 20,
                     min_editing_window = 4,
                     max_editing_window = 8,
                     flank = 20,
                     force_protospacers = FALSE,
                     force_edit_bases = c("A", "C"),
                     search_both_strands = TRUE,
                     keep_only_window_has_editable = TRUE) {

  rsid <- as.character(rsid)
  if (!nzchar(rsid)) return(NULL)

  guide_length <- as.integer(guide_length)
  min_editing_window <- as.integer(min_editing_window)
  max_editing_window <- as.integer(max_editing_window)
  flank <- as.integer(flank)

  if (min_editing_window < 1 ||
      max_editing_window > guide_length ||
      min_editing_window > max_editing_window) {
    stop("Invalid editing window. Ensure 1 <= min <= max <= guide_length.")
  }

  dat <- get_variant_mapping(rsid)
  if (is.null(dat)) return(NULL)

  # Use first mapping
  mapping <- dat$mappings[1, , drop = FALSE]
  chr <- as.character(mapping$seq_region_name)
  pos <- as.integer(mapping$start)

  if (is.na(pos) || !nzchar(chr)) return(NULL)

  # Centered flanking region
  start <- pos - flank
  end <- pos + flank
  if (start < 1) return(NULL)

  reference_seq <- get_region_sequence(chr = chr, start = start, end = end, strand = 1)
  if (is.null(reference_seq)) return(NULL)

  snp_index <- flank + 1 # 1-based position within reference_seq

  # Pattern: N...N + PAM, with IUPAC allowed
  pam.char <- pam
  guide.pattern <- paste0(paste(rep("N", guide_length), collapse = ""), pam.char) %>% DNAString()

  # Alleles string
  alleles <- mapping$allele_string %>% str_split(pattern = "/") %>% unlist()
  if (length(alleles) < 2) return(NULL)

  # -----------------------------
  # Normal-mode helper
  # -----------------------------
  guides_design_func <- function(sequence, allele1, allele2, label) {
    if ((allele1 == "C" & allele2 == "T") | (allele1 == "A" & allele2 == "G")) {
      hits <- matchPattern(pattern = guide.pattern, subject = sequence, fixed = FALSE)
      strand_label <- "forward"
      conversion <- ifelse(allele1 == "C", "C>T", "A>G")

    } else if ((allele1 == "T" & allele2 == "C") | (allele1 == "G" & allele2 == "A")) {
      sequence_rc <- reverseComplement(sequence)
      hits <- matchPattern(pattern = guide.pattern, subject = sequence_rc, fixed = FALSE)
      strand_label <- "reverse"
      conversion <- ifelse(allele1 == "T", "A>G", "C>T")

    } else {
      hits <- NULL
    }

    if (length(hits) == 0) return(NULL)

    df <- as.data.frame(hits) %>%
      select(-width) %>%
      mutate(
        snp_position = snp_index - start + 1
      ) %>%
      filter(
        snp_position >= min_editing_window,
        snp_position <= max_editing_window
      ) %>%
      mutate(
        strand = strand_label,
        targeted_allele = label,
        conversion = conversion,
        enzyme = case_when(
          conversion == "C>T" ~ "CBE",
          conversion == "A>G" ~ "ABE",
          .default = NA_character_
        ),
        protospacer = str_sub(seq, 1, guide_length),
        pam = str_sub(seq, guide_length + 1, guide_length + nchar(pam.char)),
        targeted_nuc = str_sub(seq, snp_position, snp_position)
      ) %>%
      select(-seq)

    if (nrow(df) == 0) return(NULL)
    df
  }

  # ==========================================================
  # Forced mode: protospacers MUST overlap SNP AND SNP in window
  # Enzyme assigned based on A/C presence in editing window
  # ==========================================================
  if (isTRUE(force_protospacers)) {

    build_forced_df <- function(subject_seq, strand_label) {
      hits <- matchPattern(pattern = guide.pattern, subject = subject_seq, fixed = FALSE)
      if (length(hits) == 0) return(NULL)

      df <- as.data.frame(hits) %>%
        select(-width) %>%
        mutate(
          strand = strand_label,
          protospacer = str_sub(seq, 1, guide_length),
          pam = str_sub(seq, guide_length + 1, guide_length + nchar(pam.char))
        ) %>%
        select(-seq) %>%
        mutate(
          snp_pos_in_protospacer = snp_index - start + 1
        ) %>%
        filter(
          snp_pos_in_protospacer >= 1,
          snp_pos_in_protospacer <= guide_length,
          snp_pos_in_protospacer >= min_editing_window,
          snp_pos_in_protospacer <= max_editing_window
        )

      if (nrow(df) == 0) return(NULL)
      df
    }

    df_fwd <- build_forced_df(reference_seq, "forward")

    df_rev <- NULL
    if (isTRUE(search_both_strands)) {
      df_rev <- build_forced_df(reverseComplement(reference_seq), "reverse")
    }

    guides <- bind_rows(df_fwd, df_rev)
    if (is.null(guides) || nrow(guides) == 0) return(NULL)

    window_ann <- pmap_dfr(
      list(guides$protospacer),
      function(protospacer) annotate_window(
        protospacer = protospacer,
        min_w = min_editing_window,
        max_w = max_editing_window,
        bases = force_edit_bases
      )
    )

    guides <- bind_cols(guides, window_ann) %>%
      mutate(
        conversion = case_when(
          n_A > 0 & n_C > 0 ~ "C>T;A>G",
          n_C > 0           ~ "C>T",
          n_A > 0           ~ "A>G",
          TRUE              ~ NA_character_
        ),
        enzyme = case_when(
          n_A > 0 & n_C > 0 ~ "CBE;ABE",
          n_C > 0           ~ "CBE",
          n_A > 0           ~ "ABE",
          TRUE              ~ NA_character_
        ),
        snp = dat$name,
        allele_string = mapping$allele_string,
        mode = "forced_overlap_snp_in_window",
        targeted_nuc = str_sub(protospacer, snp_pos_in_protospacer, snp_pos_in_protospacer)
      )

    if (isTRUE(keep_only_window_has_editable)) {
      guides <- guides %>% filter(n_editable > 0)
    }

    guides <- guides %>%
      select(
        snp, allele_string, mode, strand,
        protospacer, pam,
        snp_pos_in_protospacer, targeted_nuc,
        enzyme, conversion,
        editing_window, editable_bases, n_A, n_C, n_editable, editable_pos, editable_pos_by_base
      )

    return(guides)
  }

  # ==========================================================
  # Normal mode: only if editable transition exists
  # ==========================================================
  if (all(c("C", "T") %in% alleles) || all(c("A", "G") %in% alleles)) {

    guides_ref <- guides_design_func(reference_seq, alleles[1], alleles[2], paste0("reference_", alleles[1]))
    guides <- guides_ref

    # Generate alternative-allele sequences at SNP site
    seq_str <- as.character(reference_seq)

    for (i in 2:length(alleles)) {
      allele2 <- alleles[i]

      seq_alt <- DNAString(paste0(
        str_sub(seq_str, 1, snp_index - 1),
        allele2,
        str_sub(seq_str, snp_index + 1, nchar(seq_str))
      ))

      guide_alt <- guides_design_func(seq_alt, allele2, alleles[1], paste0("alternative_", allele2))
      guides <- bind_rows(guides, guide_alt)
    }

    if (is.null(guides) || nrow(guides) == 0) return(NULL)

    guides <- guides %>%
      select(-start, -end) %>%
      mutate(
        snp = dat$name,
        editing_window = str_sub(protospacer, min_editing_window, max_editing_window),
        bystander_likely = str_count(editing_window, pattern = targeted_nuc) > 1
      )

    return(guides)
  }

  NULL
}
