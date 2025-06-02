## Load libraries
library(BSgenome.Hsapiens.UCSC.hg38)
library(httr)
library(jsonlite)
library(stringr)
library(dplyr)
library(purrr)
library(Biostrings) # for reading in fasta files
library(shiny)
library(DT)

genome <- BSgenome.Hsapiens.UCSC.hg38

snpGuide <- function(rsid,
    pam = "NGN", # desired PAM
    guide_length = 20, # desired guide length
    min_editing_window=4,
    max_editing_window=8,
    flank = 20 #desired flanking length
                    ){
    
## Pull location
url <- paste0("https://rest.ensembl.org/variation/homo_sapiens/", rsid, "?content-type=application/json")
res <- GET(url)
data <- fromJSON(rawToChar(res$content))

# Extract ±20 bp around the SNP
reference_seq <- getSeq(genome, 
              paste0("chr",data$mappings$seq_region_name), 
              start = data$mappings$start - flank, 
              end = data$mappings$start + flank)

# Establish PAM information
pam.char = pam
pam.DNAString = pam %>% DNAString()
pam_length = nchar(pam.char)

# Establish guide information
guide.DNAString = rep("N", guide_length) %>%
  paste0(., collapse = "") %>%
  paste0(., pam.char) %>%
  DNAString()

##Pull alleles
alleles <- data$mappings$allele_string %>% str_split(pattern = "/") %>% unlist 

guides_design_func <- function(sequence, allele1, allele2, label) { 
    ## Design guides to reference

if (allele1 == "C" & allele2 == "T" | allele1 == "A" & allele2 == "G") { 
    guides <- matchPattern(pattern = guide.DNAString, subject = sequence, fixed = FALSE) 
    CRISPR_target <- "forward"
    conversion = ifelse(allele1 == "C", "C>T", "A>G")
        } else if (allele1 == "T" & allele2 == "C" | allele1 == "G" & allele2 == "A") { 
    reverse_sequence <- reverseComplement(sequence) # target bottom strand
    guides <- matchPattern(pattern = guide.DNAString, subject = reverse_sequence, fixed = FALSE) 
    CRISPR_target <- "reverse"
    conversion = ifelse(alleles[1] == "T", "A>G", "C>T")  
    } else { 
    guides<- NULL}

    
if (length(guides) > 0) {
    guides %<>% as.data.frame %>% select(-width) %>% 
    mutate(snp_position = 22-start) %>% 
    filter(snp_position >= min_editing_window & 
    snp_position <= max_editing_window) %>% 
        mutate(strand = CRISPR_target, 
        targeted_allele = label, 
        conversion = conversion, 
        enzyme = case_when(conversion == "C>T" ~ "CBE",
                           conversion == "A>G" ~ "ABE", .default = NA), 
        protospacer = str_sub(seq, start = 1, end = 20), 
        pam = str_sub(seq, start = 21, end = 23), 
        targeted_nuc = str_sub(seq, start = snp_position, end = snp_position)) %>% 
    select(-seq) }
    
    return(guides)}

##identify if any of the reference or alternatives are editable 
## is an editable transition present ? 
if ((all(c("C", "T") %in% alleles) | all(c("A", "G") %in% alleles)) == TRUE) {

#Reference guides
guides_ref <- guides_design_func(reference_seq, alleles[1], alleles[2], paste0("reference_", alleles[1]))

guides <- guides_ref

### for all alternative alleles - generate the same conversions. 
for (i in 2:length(alleles)) { 
    allele2 <- alleles[i]

seq_alt <- DNAString(paste0(str_sub(as.character(reference_seq), 1, 20), 
                 allele2,
                 str_sub(as.character(reference_seq), 22, 41)))
    
guide_alt <- guides_design_func(seq_alt, allele2, alleles[1], paste0("alternative_",allele2))

guides <- rbind(guides,guide_alt)}

guides <- select(guides,-start,-end) %>% mutate(snp = data$name)

### bystander calculation

guides <- mutate(guides,editing_window = str_sub(protospacer, start = min_editing_window, end = max_editing_window), 
                  bystander_likely = str_count(editing_window, pattern = targeted_nuc) > 1)

} else {guides <- NULL} 
return(guides)}
