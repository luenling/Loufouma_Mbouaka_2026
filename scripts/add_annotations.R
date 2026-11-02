library(tidyverse)
conflicts_prefer(dplyr::select)

#library(ape)
# need to add some annotations to genes, to be able to map to A. cast. neff AMDB
setwd("~/Cemm_sync/P_2024_ACRS_AcanthamoebaRNASeq_LE/")
# read aca neff gff
aca_neff_gff <- ape::read.gff("References/AmoebaDB-67_AcastellaniiNeff/AmoebaDB-67_AcastellaniiNeff.gff")
nfowl_gff <- ape::read.gff("Other_Annos/AmoebaDB-68_NfowleriTy.gff")
ngrub_gff <- ape::read.gff("Other_Annos/AmoebaDB-68_NgruberiNEG-M.gff")

# split attributes mRNA column into key value columns and add them at the end ot get all gene/transcript annotations
# Descriptions and Note has %2C = , and %3B = ; , replace
# useless descriptions:   
not_useful <- c("unspecified product", "predicted protein", "Hypothetical protein", "hypothetical protein")
aca_neff_gff <- aca_neff_gff %>% filter(type == "mRNA") %>%
  separate_longer_delim(attributes,";") %>% 
  separate_wider_delim(attributes, "=", names=c("key", "value")) %>%
  mutate(value = str_replace_all(value, c("%2C" = ",", "%3B" = ";"))) %>%
  pivot_wider(names_from = key, values_from = value) %>% 
  mutate( description = if_else(description %in% not_useful, NA, description)) %>%
  select(ID, Note_aca_neff = Note, descr_aca_neff = description)
  
nfowl_gff <- nfowl_gff %>% filter(type == "mRNA") %>%
  separate_longer_delim(attributes,";") %>% 
  separate_wider_delim(attributes, "=", names=c("key", "value")) %>%
  mutate(value = str_replace_all(value, c("%2C" = ",", "%3B" = ";"))) %>%
  pivot_wider(names_from = key, values_from = value)  %>% 
  mutate( description = if_else(description %in% not_useful, NA, description)) %>%
  select(ID, Note_n_fowl = Note, descr_n_fowl = description)

# Problem with multiple isoforms having the same start and end
# transcript id unique - needs to be extracted beforehand
ngrub_gff <- ngrub_gff %>% filter(type == "mRNA") %>%
  mutate(geneid = str_extract(attributes, "(?<=ID=)[^;]*")) %>%
  separate_longer_delim(attributes,";") %>% 
  separate_wider_delim(attributes, "=", names=c("key", "value")) %>%
  mutate(value = str_replace_all(value, c("%2C" = ",", "%3B" = ";"))) %>%
  pivot_wider(names_from = key, values_from = value) %>%
  mutate( description = if_else(description %in% not_useful, NA, description)) %>%
  select(ID, Note_n_grub = Note, descr_n_grub = description)

## get diamond results for each mapping, replace -p1 at end of nfowl sseqid string, take only the blast result with the highest bitscore or if multiple, 
## the highest query coverage
conflicts_prefer(dplyr::slice)
blast_nfowl <- read_tsv("References/Acastellanii_genome_annotation/C3/diamond_others/diamond_blastp_Nfowleri_top3.tab", 
                        col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                                      "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen")) %>%
  mutate(sseqid = gsub("-p[0-9]+$", "", sseqid)) %>%
  mutate(qcov = (qend - qstart)/qlen, scov = (send - sstart)/slen) 
#%>% 
#  group_by(qseqid) %>% slice(which.max(bitscore)) %>% slice(which.max(qcov))
blast_ngrub <- read_tsv("References/Acastellanii_genome_annotation/C3/diamond_others/diamond_blastp_Ngruberi_top3.tab", 
                        col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                                      "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen")) %>%
  mutate(sseqid = gsub("-p[0-9]+$", "", sseqid)) %>%
  mutate(qcov = (qend - qstart)/qlen, scov = (send - sstart)/slen) 
blast_neff <- read_tsv("References/Acastellanii_genome_annotation/C3/diamond_others/diamond_blastp_neff_aDB_neff_top3.tab", 
                        col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                                      "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen")) %>%
  mutate(sseqid = gsub("-p[0-9]+$", "", sseqid)) %>%
  mutate(qcov = (qend - qstart)/qlen, scov = (send - sstart)/slen) 

# for start only take best bitscore, if same, take highest qcov, pident, and qlen

neff_annos <- blast_neff %>% group_by(qseqid) %>% arrange(-bitscore, -pident, -qcov, .by_group = TRUE) %>%
  filter(sseqid != "*") %>% 
  group_by(qseqid) %>% arrange(-bitscore, -pident, -qcov, .by_group = TRUE) %>%
  slice(1) %>% ungroup() %>%
  select(qseqid, sseqid, evalue, bitscore, pident, length, qcov) %>%
  rename_with( ~ paste0(.x, "_aca_neff", recycle0 = TRUE), 
               all_of(c("evalue", "bitscore", "pident", "length", "qcov"))) %>%
  left_join(aca_neff_gff, by = join_by(sseqid == ID), suffix = c("", "_aca_neff")) %>% 
  dplyr::rename( Aca_neff_ID = sseqid) 
nfowl_annos <- blast_nfowl %>% filter(sseqid != "*") %>% 
  group_by(qseqid) %>% arrange(-bitscore, -pident, -qcov, .by_group = TRUE) %>%
  slice(1) %>% ungroup() %>%
  select(qseqid, sseqid, evalue, bitscore, pident, length, qcov) %>%
  rename_with( ~ paste0(.x, "_n_fowl", recycle0 = TRUE), 
               all_of(c("evalue", "bitscore", "pident", "length", "qcov"))) %>%
  left_join(nfowl_gff, by = join_by(sseqid == ID), suffix = c("", "_n_fowl")) %>% 
  dplyr::rename( N_fowl_ID = sseqid) 
ngrub_annos <- blast_ngrub %>% group_by(qseqid) %>% arrange(-bitscore, -pident, -qcov, .by_group = TRUE) %>%
  filter(sseqid != "*") %>% 
  group_by(qseqid) %>% arrange(-bitscore, -pident, -qcov, .by_group = TRUE) %>%
  slice(1) %>% ungroup() %>%
  select(qseqid, sseqid, evalue, bitscore, pident, length, qcov) %>%
  rename_with( ~ paste0(.x, "_n_grub", recycle0 = TRUE), 
               all_of(c("evalue", "bitscore", "pident", "length", "qcov"))) %>%
  left_join(ngrub_gff, by = join_by(sseqid == ID)) %>% 
  dplyr::rename( N_grub_ID = sseqid) 

all_annos <- full_join(neff_annos, nfowl_annos, by="qseqid") %>% full_join(ngrub_annos) %>% 
  mutate(qseqid = gsub("-T[0-9]+$", "", qseqid)) %>% dplyr::rename(GeneID = qseqid)
write_tsv(all_annos, "References/Acastellanii_genome_annotation/C3/blasted_C3_vs_acaneff_nfowl_ngrub_annos.tsv")
library(openxlsx)
openxlsx::write.xlsx(all_annos, file = "blasted_C3_vs_acaneff_nfowl_ngrub_annos.xlsx",  asTable = TRUE)
openxlsx::write.xlsx(genes, file = "all_gene_annos_C3_annos.xlsx",  asTable = TRUE)
