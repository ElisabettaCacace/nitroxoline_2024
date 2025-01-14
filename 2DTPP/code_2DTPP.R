library(TPP2D)
library(tidyverse)
library(readxl)
library(ggplot2)
library(openxlsx)

inputdir = "/Users/elisabettacacace/Documents/EMBL/Nx/Nx_paper/REVISION/nitroxoline_2024/2DTPP/input/"
outputdir = "/Users/elisabettacacace/Documents/EMBL/Nx/Nx_paper/REVISION/nitroxoline_2024/2DTPP/output/"

Cmpd_raw_df <- read_xlsx(paste0(inputdir,'2020-07-27_results_2D_TPP-nitroxoline.xlsx'), sheet = 'pEC50') %>% 
  filter(!grepl("##", protein_id)) %>% 
  dplyr::select(representative = protein_id,
                clustername = gene_name,
                experiment,
                qupm,
                qssm,
                temperature,
                matches("signal_sum"),
                matches("norm_rel_fc_protein"),
                -matches("transformed"))  %>% 
  gather(key, value, matches("signal_sum"), matches("norm_rel_fc_protein")) %>% 
  mutate(conc = as.numeric(gsub("norm_rel_fc_protein_", "", gsub("signal_sum_", "", gsub("_unmodified", "", key)))),
         key = case_when(grepl("signal_sum_", key) ~ "raw_value",
                         grepl("rel_fc", key) ~ "rel_value")) %>% 
  spread(key, value) %>% 
  arrange(representative, temperature, conc) %>% 
  group_by(clustername, temperature, conc) %>% 
  filter(qupm == max(qupm), 
         qssm == max(qssm), 
         raw_value == max(raw_value)) %>% 
  filter(!duplicated(clustername)) %>% 
  ungroup %>% 
  mutate(log2_value = log2(raw_value),
         log_conc = log10(conc/1e6)) %>% 
  filter(qupm > 1) %>% 
  arrange(clustername, temperature) %>% 
  filter(!is.na(log2_value),
         !is.na(rel_value),
         !is.na(clustername))


# resolve ambiguous protein names
Cmpd_raw_fil_df <- resolveAmbiguousProteinNames(Cmpd_raw_df)

# recompute reporter ion signal from robust Isobarquant fold changes
Cmpd_df <- recomputeSignalFromRatios(Cmpd_raw_fil_df)


#Compute null and alternative model fits and extract parameters (this may run for 10 min or so)
Cmpd_params_df <- getModelParamsDf(Cmpd_df, maxit = 500, minObs = 24)
Cmpd_fstat_df <- computeFStatFromParams(Cmpd_params_df)


#Compute a bootstrapping derived null distribution / this may take a while / parallelize by choosing appropriate BiocParallel BPPARM option
BPPARAM <-  BiocParallel::SerialParam(progressbar = TRUE)


set.seed(12, kind = "L'Ecuyer-CMRG")
Cmpd_null_df <- bootstrapNullAlternativeModel(
  df = Cmpd_df,# %>% filter(!clustername %in% c("ATP5MG", "EIF1","RAD50", "TMED7")), 
  params_df = Cmpd_params_df, 
  maxit = 500, B = 20,
  minObs = 24,
  BPPARAM = BPPARAM,
  verbose = FALSE)

# Get FDR
Cmpd_fdr_df <- getFDR(Cmpd_fstat_df, Cmpd_null_df)
Cmpd_hits_df <- findHits(Cmpd_fdr_df, 0.1)


subsets <- list("fdr_tab" = Cmpd_fdr_df %>% filter(dataset == "true"), 
                "hits" = Cmpd_hits_df,
                "fil_df" = Cmpd_df)

write.table(Cmpd_fdr_df %>% filter(dataset == "true"),paste0(outputdir,"fdr_df.txt"), sep = "\t", quote = FALSE,row.names = FALSE)
write.table(Cmpd_hits_df,paste0(outputdir,"hits_df.txt"), sep = "\t", quote = FALSE,row.names = FALSE)

