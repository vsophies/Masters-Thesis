
####################Master Thesis - Supplementary Material: S1 R script#######
##############################################################################

# Overview ----

## Description of the data set ----

#Treatment = Bti, Ctrl   
#Ponds = 12 Ponds: Ctrl = 1,3,5,7,9,11; Bti = 2,4,6,8,10,12   
#Applications (Date (Doy,Week)) = 14.4 (105,16), 4.5 (125,18), 25.5.2020 (146,21) 
#Sampling_Date = 26 Sampling events, from 14.4 - 30.7   
#Doy = Day of the year 
#Week = Week (of the year)  
#OTU  = Operational taxonomic unit   
#Number_of_Individuals (Ind) = Number of Individuals per sample 
#OTU_sum = sum of OTUs per sample 
#Cumsum/OTU_cumsum = cumulative sum of OTUs detected over time (only first detection) 
#Ind_per_OTU = Individuals per Sample / Sum of OTUs per Sample 
#Start of emergence = Day of first detection of an OTU 
#Duration of Emergence = Days betw. first and last detection of an OTU 
#nrOTUs = non-responsive OTUs (no significant treatment effect observed(GLM)) 
#rOTUs = responsive OTUs (significant treatment effect observed (GLM)) 
  
# Analysis ----
## load data 
setwd("~/M_Sc_Ecotoxicology/Masters_Thesis/Data_Analysis/Supplementary_Material")
data_all_OTUs <- read.csv("S2.2_OTU_table.csv")
env <- read.csv("S2.3_Environmental_data.csv")
OTUs_complete <- read.csv("S2.4_OTU_identification.csv")

## install & load packages
if (!require("plyr")) install.packages("plyr"); library(plyr)
if (!require("bookdown")) install.packages("bookdown"); library(bookdown)
if (!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if (!require("proxy")) install.packages("proxy"); library(proxy)
if (!require("vegan")) install.packages("vegan"); library(vegan)
if (!require("gplots")) install.packages("gplots"); library(gplots)
if (!require("gplots")) install.packages("gplots"); library(ggplot2)
if (!require("ggpubr")) install.packages("ggpubr"); library(ggpubr)
if (!require("betapart")) install.packages("betapart"); library(betapart)
if (!require("qpcR")) install.packages("qpcR"); library(qpcR)
if (!require("seqinr")) install.packages("seqinr"); library(seqinr)
if (!require("mgcv")) install.packages("mgcv"); library(mgcv)
if (!require("data.table")) install.packages("data.table"); library(data.table)
if (!require("synchrony")) install.packages("synchrony"); library(synchrony)
if (!require("mvabund")) install.packages("mvabund"); library(mvabund)
if (!require("reshape2")) install.packages("reshape2"); library(reshape2)
if (!require("janitor")) install.packages("janitor"); library(janitor)
if (!require("ggforce")) install.packages("ggforce"); library(ggforce)
if (!require("knitr")) install.packages("knitr"); library(knitr)


# Bioinformatic analysis ----
  ## exploration of identified OTUs
OTUs <- OTUs_complete %>% 
  rename( OTU_ID = 1) %>%
  filter(!OTU_ID == "")

OTUs_Chiro <- OTUs %>%
  filter(Family == "Chironomidae") # only chiros = 185 OTUs

OTUs_Chiro %>%
  group_by(Result) %>%
  transmute(count = n()) %>%
  unique() %>%
  arrange(desc(count)) # species names assigned to more than one OTU

OTUs %>%
  filter(Family == "No match") # 113 OTUs not identified

OTUs %>%
  filter(Family != "Chironomidae" & Family != "No match") %>%
  filter(Phylum != "Arthropoda") # 14 OTUs belong to phyla other than Arthropoda

OTUs %>%
  filter(Family != "Chironomidae" & Family != "No match") %>%
  filter(Phylum == "Arthropoda") %>%
  filter(Class != "Insecta") # 11 OTUs belong to classes other than Insecta

OTUs %>%
  filter(Family != "Chironomidae" & Family != "No match") %>%
  filter(Phylum == "Arthropoda") %>%
  filter(Class == "Insecta") %>%
  filter(Order != "Diptera") # 10 OTUs belong to orders other than Diptera

OTUs %>%
  filter(Family != "Chironomidae" & Family != "No match") %>%
  filter(Phylum == "Arthropoda") %>%
  filter(Class == "Insecta") %>%
  filter(Order == "Diptera") %>%
  filter(Family != "Chironomidae") # 19 OTUs belong to families other than Chironomidae

  ## prepare data
Chiro_OTUs_vec <- OTUs_Chiro$OTU_ID
data <- data_all_OTUs %>% 
  filter(ID %in% Chiro_OTUs_vec )
data_flip <- as.data.frame(t(data)) 
data_flip <- tibble::rownames_to_column(data_flip, "ID") 
colnames(data_flip) <- data_flip[2,]
data_flip <- data_flip[-(1:2), ]
data_flip <- data_flip %>%
  mutate_at(vars(starts_with("OTU_")), funs(as.numeric)) %>% 
  filter(!ID == "sequ") %>% 
  mutate(ID = str_remove(ID, '_B.*'))

## Octave plot ----
abundance_repl <- data_flip %>%
  mutate(ID = str_remove(ID, '_B.*')) %>% 
  remove_rownames %>%
  column_to_rownames(var = "ID") 

abundance_colSums <- as.data.frame(abundance_repl) %>%
  colSums(.) %>%
  data.frame() %>%
  filter(.!= 0) %>%
  mutate(sum = case_when(
    . <= 1 ~ "1",
    . <= 2 ~ "2",
    . <= 4 ~ "4",
    . <= 8 ~ "8",
    . <= 16 ~ "16",
    . <= 32 ~ "32",
    . <= 64 ~ "64",
    . <= 128 ~ "128",
    . <= 256 ~ "256",
    . <= 512 ~ "512",
    . <= 1000 ~ "1k",
    . <= 2000 ~ "2k",
    . <= 4000 ~ "4k",
    . <= 8000 ~ "8k",
    . <= 16000 ~ "16k",
    . <= 32000 ~ "32k",
    . <= 64000 ~ "64k",
    . <= 128000 ~ "128k",
    . <= 256000 ~ "256k",
    . <= 512000~ "512k",
    . <= 1000000 ~ "1 mio",
    . <= 2000000 ~ "2 mio",
    . <= 3000000 ~ "3 mio"))

abundance_colSums$sum <- factor(abundance_colSums$sum, levels = c("1",
                                                                  "2",
                                                                  "4",
                                                                  "8",
                                                                  "16",
                                                                  "32",
                                                                  "64",
                                                                  "128",
                                                                  "256",
                                                                  "512",
                                                                  "1k",
                                                                  "2k",
                                                                  "4k",
                                                                  "8k",
                                                                  "16k",
                                                                  "32k",
                                                                  "64k",
                                                                  "128k",
                                                                  "256k",
                                                                  "512k",
                                                                  "1 mio",
                                                                  "2 mio",
                                                                  "3 mio"))


ggplot(abundance_colSums, aes(sum)) +
  geom_histogram(stat = "count") + 
  ylab("Number of OTUs") +
  xlab("Number of reads")


# Data analysis ----

  ## check number of replicates per sample
data_flip %>%
  group_by(str_extract(ID, "[^_]+"))  %>% 
  summarize(n = n()) %>% 
  mutate(n = as.numeric(n)) %>% 
  arrange(desc(n)) 

## Jaccard dissimilarity between replicates ----
  ## abundance matrix
abundance_repl <- data_flip %>%
  remove_rownames %>%
  column_to_rownames(var = "ID") 

  ## occurrence matrix
occurrence_repl <- data_flip %>%
  mutate_if(is.numeric, ~1 * (. > 0)) %>% 
  mutate(ID = str_remove(ID, '_B.*')) %>% 
  remove_rownames %>%
  column_to_rownames(var = "ID")

dist_ab <- as.matrix(proxy::dist(abundance_repl, by_rows = TRUE, method = "Bray")) # calculate Bray Curtis dissimilarity
dist_oc <- as.matrix(proxy::dist(occurrence_repl, by_rows = TRUE, method = "Jaccard")) # calculate Jaccard dissimilarity

  ## extract only comparisons between replicates
row.ind_ab <- grep("_t",rownames(dist_ab),value=TRUE)
col.ind_ab <- match(sub("_t","",row.ind_ab),colnames(dist_ab))
dist_repl_ab <- data.frame(ID=colnames(dist_ab)[col.ind_ab],
                           dist=diag(dist_ab[row.ind_ab,col.ind_ab])
)

row.ind_oc <- grep("_t",rownames(dist_oc),value=TRUE)
col.ind_oc <- match(sub("_t","",row.ind_oc),colnames(dist_oc))
dist_repl_oc <- data.frame(ID=colnames(dist_oc)[col.ind_oc],
                           dist=diag(dist_oc[row.ind_oc,col.ind_oc])
)

dist_repl_ab %>%
  arrange(desc(dist)) %>%
  rename(Abundance_matrix = dist)

dist_repl_oc %>%
  arrange(desc(dist)) %>%
  rename(Occurrence_matrix = dist)

ab_oc <- data.frame(dist_repl_ab, dist_repl_oc) %>%
  arrange(desc(dist)) %>%
  dplyr::select(!ID.1) %>%
  rename(Abundance = dist, Occurrence = dist.1) 

ab_oc # Jaccard dissimilarity between replicates (read abundance/occurrence matrix)

ab_oc_above <- ab_oc %>%
  filter(Abundance > 0.4 | Occurrence > 0.4) 

ab_oc_above # replicates with Jaccard dissimilarity > 0.4 in either read abundance or 
              # occurrence matrix


## Extract means of replicates 
  ## means of occurrence: 1, 0, 0.5 
means_unclean <- occurrence_repl %>%
  rownames_to_column("ID") %>% 
  filter(ID !="MT245" &  ID != "MT297_t" & ID !="MT47_t") %>% # remove because only one 
  group_by(ID = str_extract(ID, "[^_]+"))  %>%                  #replicate is used
  summarize_at(vars(starts_with("OTU_")), mean)

means_unclean2 <- means_unclean %>% 
  remove_rownames %>%
  column_to_rownames(var = "ID") %>%
  transmute(., replace(., . == 0.5, 0))

means_plusrare <- means_unclean2 %>%
  dplyr::select(which(!colSums(means_unclean2) %in% 0)) # remove columns that sum to zero

means_colsumiszero <- means_unclean2 %>%
  dplyr::select(which(colSums(means_unclean2) %in% 0)) # which OTUs are removed?

means_plusrare <- means_plusrare %>%
  rownames_to_column("Sample_ID") %>%
  filter(Sample_ID != "MT8" & Sample_ID != "MT43" & Sample_ID != "MT238" ) %>%
  column_to_rownames("Sample_ID") # remove samples with too high Jacc. diss.


means_colsumiszero %>% colnames() # OTUs removed after combining replicates


## Species accumulation curves ----
  ## prepare environmental data
env2 <- env %>% 
  group_by(Sample_ID, Sampling_Date, Pond) %>%
  dplyr::summarise(Number_of_Individuals = sum(Number_of_Individuals)) %>%
  mutate(Treatment = ifelse(Pond == 1 | Pond == 3 | Pond == 5 | Pond == 7 | Pond == 9 | 
                              Pond == 11, "Ctrl", "Bti")) 

env2$Sampling_Date <- as.Date(env2$Sampling_Date, format = "%d/%m/%Y")
env2 <- env2 %>%
  mutate(Pond = as.factor(Pond)) %>%
  mutate(Treatment = as.factor(Treatment)) %>%
  mutate(Doy = as.factor(yday(Sampling_Date))) %>%
  mutate(Week = as.factor(week(Sampling_Date)))

  ## subset community data for each Pond
env2 <- env2 %>%
  mutate(Pond = paste("Pond", Pond, sep="")) %>%
  mutate(Pond = as.factor(Pond))

env2 <- env2 %>%
  filter(Sample_ID != "MT8" & Sample_ID != "MT43" & Sample_ID != "MT238" )


list_df <- split(env2, env2$Pond) 

split_frames <- lapply(list_df, function(x)
  means_plusrare %>%
    rownames_to_column("ID") %>% 
    filter(ID %in% x$Sample_ID) %>%
    remove_rownames %>%
    column_to_rownames(var = "ID")
)

list2env(split_frames,envir=.GlobalEnv) 

  ## calculate species accumulation curves per Pond using specpool from 'vegan'
set.seed(10)
specpool(Pond1)
specpool(Pond2)
specpool(Pond3)
specpool(Pond4)
specpool(Pond5)
specpool(Pond6)
specpool(Pond7)
specpool(Pond8)
specpool(Pond9)
specpool(Pond10)
specpool(Pond11)
specpool(Pond12)

  ### calculate species accumulation curves per Treatment using specpool from 'vegan'

list_df_Treatment <- split(env2, env2$Treatment) # split environmental data for each Treatment

split_frames_Treatment <- lapply(list_df_Treatment, function(x)
  means_plusrare %>%
    rownames_to_column("ID") %>% 
    filter(ID %in% x$Sample_ID) %>%
    remove_rownames %>%
    column_to_rownames(var = "ID")
) 

list2env(split_frames_Treatment,envir=.GlobalEnv) 

set.seed(10)
specpool(Ctrl)
specpool(Bti)

## OTU occurrence ----
  ## Number of occurrences per OTU
OTU_occurrence <- means_plusrare %>%
  colSums(.) %>%
  as.data.frame() %>%
  rename(Total_occurrences = ".") %>%
  rownames_to_column("OTU_ID")

OTU_occurrence %>% 
  filter(Total_occurrences > 150) # 10 OTUs occur more than 150 times

OTU_occurrence %>% 
  filter(!Total_occurrences > 150) %>%
  filter(Total_occurrences > 50) # 34 OTUs occur between 50 and 150 times

OTU_occurrence %>% 
  filter(!Total_occurrences > 50) # 124 OTUs occur less than 50 times

OTU_occurrence %>%
  filter(Total_occurrences == 1 | Total_occurrences == 2) # 56 OTUs occur only 1 or 2 times in the whole data set and will be removed for further analysis.


####################################### check occurrence in ponds

means <- means_plusrare %>%
  dplyr::select(!
                  OTU_occurrence %>%
                  filter(Total_occurrences == 1 | Total_occurrences == 2) %>%
                  pull(OTU_ID)
  ) # remove rare OTUs


## Community composition ----
  ## Number of OTUs 

Ctrl_spec <- c(70, 72,88,84,92,80) 
Bti_spec <- c(77,79,91,83,72,66) 
No_spec_Pond <- melt(as.data.frame(cbind(Ctrl_spec, Bti_spec)))


  ### Descriptive statistics 
Final_cumsum_ds <- as.data.frame(matrix(c(mean(Ctrl_spec), mean(Bti_spec),
                                          sd(Ctrl_spec), sd(Bti_spec)), ncol = 2, byrow = TRUE)) %>% rename(Ctrl = V1,  Bti = V2) %>% mutate(Treatment = c("mean", "sd")) %>%
  column_to_rownames("Treatment")


  ### ANOVA 
aov_No_spec <- aov(value ~ variable, data = No_spec_Pond)
summary(aov_No_spec)
plot(aov_No_spec)

  ### Effect size
aov_fin_cumsum_all_effectsize <- effectsize::eta_squared(aov_No_spec, digits = 3)


  ## prepare data
vPond1 <- split_frames$Pond1 %>%
  dplyr::select(which(!colSums(.) %in% 0)) %>% 
  colnames(.)
vPond2 <- split_frames$Pond2 %>%
  dplyr::select(which(!colSums(.) %in% 0)) %>%  
  colnames(.)
vPond3 <- split_frames$Pond3 %>%
  dplyr::select(which(!colSums(.) %in% 0)) %>%  
  colnames(.)
vPond4 <- split_frames$Pond4 %>%
  dplyr::select(which(!colSums(.) %in% 0)) %>%  
  colnames(.)
vPond5 <- split_frames$Pond5 %>%
  dplyr::select(which(!colSums(.) %in% 0)) %>%  
  colnames(.)
vPond6 <- split_frames$Pond6 %>%
  dplyr::select(which(!colSums(.) %in% 0)) %>%  
  colnames(.)
vPond7 <- split_frames$Pond7 %>%
  dplyr::select(which(!colSums(.) %in% 0)) %>%  
  colnames(.)
vPond8 <- split_frames$Pond8 %>%
  dplyr::select(which(!colSums(.) %in% 0)) %>%  
  colnames(.)
vPond9 <- split_frames$Pond9 %>%
  dplyr::select(which(!colSums(.) %in% 0)) %>%  
  colnames(.)
vPond10 <- split_frames$Pond10 %>%
  dplyr::select(which(!colSums(.) %in% 0)) %>%  
  colnames(.)
vPond11 <- split_frames$Pond11 %>%
  dplyr::select(which(!colSums(.) %in% 0)) %>%  
  colnames(.)
vPond12 <- split_frames$Pond12 %>%
  dplyr::select(which(!colSums(.) %in% 0)) %>%  
  colnames(.)

Ctrl_OTUs <- unique(c(vPond1,vPond3,vPond5,vPond7,vPond9,vPond11))
Bti_OTUs <- unique(c(vPond2,vPond4,vPond6,vPond8,vPond10,vPond12))

OTU_list <- list(Control = Ctrl_OTUs, Bti = Bti_OTUs)

  ## Pond-specific OTU table
Pond_OTUs <- c(vPond1,vPond3,vPond5,vPond7,vPond9,vPond11,vPond2,vPond4,vPond6,vPond8,vPond10,vPond12)

OTUs_Chiro_without12 <- OTUs_Chiro %>%
  filter(OTU_ID %in% colnames(means))

Pond_OTUs_df <- OTUs_Chiro_without12 %>%
  dplyr::select(OTU_ID, Result) %>%
  mutate(Pond1 = if_else(OTU_ID %in% vPond1, 1, 0)) %>%
  mutate(Pond2 = if_else(OTU_ID %in% vPond2, 1, 0)) %>%
  mutate(Pond3 = if_else(OTU_ID %in% vPond3, 1, 0)) %>%
  mutate(Pond4 = if_else(OTU_ID %in% vPond4, 1, 0)) %>%
  mutate(Pond5 = if_else(OTU_ID %in% vPond5, 1, 0)) %>%
  mutate(Pond6 = if_else(OTU_ID %in% vPond6, 1, 0)) %>%
  mutate(Pond7 = if_else(OTU_ID %in% vPond7, 1, 0)) %>%
  mutate(Pond8 = if_else(OTU_ID %in% vPond8, 1, 0)) %>%
  mutate(Pond9 = if_else(OTU_ID %in% vPond9, 1, 0)) %>%
  mutate(Pond10 = if_else(OTU_ID %in% vPond10, 1, 0)) %>%
  mutate(Pond11 = if_else(OTU_ID %in% vPond11, 1, 0)) %>%
  mutate(Pond12 = if_else(OTU_ID %in% vPond12, 1, 0)) %>%
  mutate(Control = Pond1 + Pond3 + Pond5 + Pond7 + Pond9 + Pond11) %>%
  mutate(Bti = Pond2 + Pond4 + Pond6 + Pond8 + Pond10 + Pond12) %>%
  filter(!(Control == 0 & Bti == 0)) %>%
  bind_rows(summarise_all(.[1:14], ~if(is.numeric(.)) sum(.) else "Total"))

Ctrl_Bti_OTUs <- Pond_OTUs_df %>%
  dplyr::select(OTU_ID, Result, Control, Bti) %>% filter(Bti != "NA")

Ctrl_Bti_OTUs %>%
  filter(Control >= 0 & Bti == 0) # OTUs only in Ctrl

Ctrl_Bti_OTUs %>%
  filter(Bti >= 0 & Control == 0) # OTUs only in Bti

OTU_occ_Ponds <- Pond_OTUs_df %>%
  filter(!OTU_ID == "Total") %>%
  dplyr::select(-c(Result, Control, Bti)) %>%
  mutate_at(vars(-OTU_ID), as.numeric) # OTU occurrence per Pond

Pond_OTU_occ <- as.data.frame(t(OTU_occ_Ponds))
colnames(Pond_OTU_occ) <- Pond_OTU_occ[1,]
Pond_OTU_occ <- Pond_OTU_occ[-(1), ]
Pond_OTU_occ <- Pond_OTU_occ %>%
  mutate_if(is.character, as.numeric)

OTU_list_CtrlBti <- list(P1C = vPond1, P3C = vPond3, P5C = vPond5, P7C = vPond7, P9C = vPond9, P11C = vPond11,
                         P2B = vPond2, P4B = vPond4, P6B = vPond6, P8B = vPond8, P10B = vPond10, P12B = vPond12)

Pond <- colnames(OTU_occ_Ponds[2:13])

## OTUs occurring only 1 or 2 times: distribution in treatments?

OTUs_Chiro_with12 <- OTUs_Chiro %>%
  filter(!OTU_ID %in% colnames(means))

Pond_OTUs_df_test <- OTUs_Chiro_with12 %>%
  dplyr::select(OTU_ID, Result) %>%
  mutate(Pond1 = if_else(OTU_ID %in% vPond1, 1, 0)) %>%
  mutate(Pond2 = if_else(OTU_ID %in% vPond2, 1, 0)) %>%
  mutate(Pond3 = if_else(OTU_ID %in% vPond3, 1, 0)) %>%
  mutate(Pond4 = if_else(OTU_ID %in% vPond4, 1, 0)) %>%
  mutate(Pond5 = if_else(OTU_ID %in% vPond5, 1, 0)) %>%
  mutate(Pond6 = if_else(OTU_ID %in% vPond6, 1, 0)) %>%
  mutate(Pond7 = if_else(OTU_ID %in% vPond7, 1, 0)) %>%
  mutate(Pond8 = if_else(OTU_ID %in% vPond8, 1, 0)) %>%
  mutate(Pond9 = if_else(OTU_ID %in% vPond9, 1, 0)) %>%
  mutate(Pond10 = if_else(OTU_ID %in% vPond10, 1, 0)) %>%
  mutate(Pond11 = if_else(OTU_ID %in% vPond11, 1, 0)) %>%
  mutate(Pond12 = if_else(OTU_ID %in% vPond12, 1, 0)) %>%
  mutate(Control = Pond1 + Pond3 + Pond5 + Pond7 + Pond9 + Pond11) %>%
  mutate(Bti = Pond2 + Pond4 + Pond6 + Pond8 + Pond10 + Pond12) %>%
  filter(!(Control == 0 & Bti == 0)) %>%
  bind_rows(summarise_all(.[1:14], ~if(is.numeric(.)) sum(.) else "Total"))

Ctrl_Bti_OTUs_test <- Pond_OTUs_df_test %>%
  dplyr::select(OTU_ID, Control, Bti) %>% filter(Bti != "NA")

Ctrl_Bti_OTUs_test %>%
  filter(Control >= 0 & Bti == 0) # 25 OTUs only in Ctrl

Ctrl_Bti_OTUs_test %>%
  filter(Bti >= 0 & Control == 0) # 23 OTUs only in Bti

Ctrl_Bti_OTUs_test %>%
  mutate(All = Bti + Control) %>%
  filter(All < 3) %>%
  filter(Bti > 0 & Control > 0 ) # 8 OTUs in both treatments


## Community dynamics ----

### Multivariate GLM ----

means_abund <- mvabund(means)

  ## full model
set.seed(99)
mglm_full = manyglm(means_abund ~ Treatment + Pond + Doy + Treatment:Doy + Treatment:Pond, data = env2, family = binomial("cloglog"))
plot(mglm_full, which = 1:2)
mglm_aov <- anova(mglm_full, nBoot = 999, show.time = "all")
  ## fit single GLMs per OTU
set.seed(99)
SingleGLM_full <- anova(mglm_full, nBoot = 999, show.time = "all", p.uni = "unadjusted")
SingleGLM_full
pvals_full <- as.data.frame(SingleGLM_full$uni.p[2,])
pvalssign_full <- pvals_full %>%
  filter(pvals_full$`SingleGLM_full$uni.p[2, ]` < 0.05)
responsive_OTUs <- row.names(pvalssign_full)


  ## Deviance over time (resp. OTUs)
mod_pt_full <- NULL
for (i in levels(factor(env2$Sampling_Date))) {
  take_abu <- means_abund[env2$Sampling_Date == i, ]
  take_env <- env2[env2$Sampling_Date == i, ]
  
  mod_pt_full[[i]]$mod <- manyglm(take_abu ~ Treatment + Pond + Doy + Treatment:Doy + 
                                    Treatment:Pond, data = take_env, 
                                  family = binomial("cloglog"))
  mod_pt_full[[i]]$aov <- anova(mod_pt_full[[i]]$mod, nBoot = 100, 
                                p.uni = 'adjusted', test = 'LR', show.time = "none")
  mod_pt_full[[i]]$sum <- summary(mod_pt_full[[i]]$mod, nBoot = 100, 
                                  p.uni = 'adjusted', test = 'LR')
}

get_pvals_full <- function(x){
  comm <- c(community = x$aov$table[2, 4])
  spec <- x$aov$uni.p[2, ]
  c(comm, spec)
}


Sampling_Datelypvals_full <- ldply(mod_pt_full, get_pvals_full)

devs_full <- ldply(mod_pt_full, function(x) x$aov$uni.test[2, ])
devs_full <-  devs_full %>%
  dplyr::select(c(".id", all_of(responsive_OTUs)))

plotdf_full <- melt(devs_full, id.vars = '.id')
plotdf_full <- plotdf_full %>%
  rename(Sampling_Date = .id) %>%
  mutate(Sampling_Date = as.Date(Sampling_Date)) %>%
  rename(OTU_ID = variable)
 
plotdf_full # deviance per rOTU and Sampling Date 

### Sum of OTUs ----

  ## prepare data
OTUs_perMT <- means %>%
  transmute(OTU_sum = rowSums(.)) %>%
  rownames_to_column("Sample_ID")

GLAM_data <- merge(env2, OTUs_perMT, by = "Sample_ID")
GLAM_data$Sampling_Date <- as.Date(GLAM_data$Sampling_Date, format = "%d/%m/%y")

GLAM_data <- GLAM_data %>%
  mutate(Doy = yday(Sampling_Date)) %>% # Doy = Day of year
  mutate(Week = week(Sampling_Date)) %>%
  mutate(Ind_per_OTU = Number_of_Individuals / OTU_sum) %>%
  mutate(OTUs_per_Ind = OTU_sum / Number_of_Individuals)

GLAM_data %>%
  filter(Ind_per_OTU < 1) # No sample with more OTUs than individuals 

GLAM_data_Ctrl <- GLAM_data %>%
  filter(Treatment == "Ctrl")

GLAM_data_Bti <- GLAM_data %>%
  filter(Treatment == "Bti")

#### Descr. statistics
Sum_ds_all <- as.data.frame(matrix(c(mean(GLAM_data_Ctrl$OTU_sum), mean(GLAM_data_Bti$OTU_sum),
                                     sd(GLAM_data_Ctrl$OTU_sum), sd(GLAM_data_Bti$OTU_sum)), ncol = 2, byrow = TRUE)) %>% rename(Ctrl = V1,  Bti = V2) %>% mutate(Treatment = c("mean", "sd")) %>%
  column_to_rownames("Treatment")

#### Two-way ANOVA
aov_sum <- aov(OTU_sum ~ Treatment * as.factor(Week), data = GLAM_data)
summary(aov_sum)
plot(aov_sum)

aov_sum_effectsize_all <- effectsize::eta_squared(aov_sum)

#### Ds per week

  #### Mean
Sum_mean_Bti <- lapply(c(16:31), function(x)
  mean(GLAM_data_Bti[GLAM_data_Bti$Week == x, ]$OTU_sum)
)

Sum_mean_Ctrl <- lapply(c(16:31), function(x)
  mean(GLAM_data_Ctrl[GLAM_data_Ctrl$Week == x, ]$OTU_sum)
)

names(Sum_mean_Bti) <- c(16:31)
names(Sum_mean_Ctrl) <- c(16:31)

Sum_mean_Bti_df <- as.data.frame(do.call(rbind, Sum_mean_Bti))
Sum_mean_Ctrl_df <- as.data.frame(do.call(rbind, Sum_mean_Ctrl))

Sum_mean <- cbind(Sum_mean_Ctrl_df, mean_Bti = Sum_mean_Bti_df$V1) %>% rename(mean_Ctrl = V1)

  #### Standard deviation
Sum_sd_Bti <- lapply(c(16:31), function(x)
  sd(GLAM_data_Bti[GLAM_data_Bti$Week == x, ]$OTU_sum)
)
Sum_sd_Ctrl <- lapply(c(16:31), function(x)
  sd(GLAM_data_Ctrl[GLAM_data_Ctrl$Week == x, ]$OTU_sum)
)

names(Sum_sd_Bti) <- c(16:31)
names(Sum_sd_Ctrl) <- c(16:31)

Sum_sd_Bti_df <- as.data.frame(do.call(rbind, Sum_sd_Bti))
Sum_sd_Ctrl_df <- as.data.frame(do.call(rbind, Sum_sd_Ctrl))

Sum_sd <- cbind(Sum_sd_Ctrl_df, sd_Bti = Sum_sd_Bti_df$V1) %>% 
  rename(sd_Ctrl = V1)

Sum_ds <- cbind(Sum_mean, Sum_sd) %>% relocate(sd_Ctrl, .after = mean_Ctrl)

#### ANOVA per week
Sum_aov <- lapply(c(16:31), function(x) 
  aov(OTU_sum ~ Treatment, data = GLAM_data[GLAM_data$Week == x, ])) 

Sum_anovas <- lapply(c(16:31), function(x) 
  summary(aov(OTU_sum ~ Treatment, data = GLAM_data[GLAM_data$Week == x, ]))) 

plot_Sum_aov <- lapply(Sum_aov, function(x)
  plot(x)
) # hit Esc if you can't stop plotting ;)

Sum_aov_effectsize <- lapply(c(1:16), function(x)
  effectsize::eta_squared(Sum_aov[[x]])
)

Sum_aov_data <- lapply(c(1:16), function(x)
  Sum_aov_effectsize[[x]]["Eta2"]
) 


### Sum of nrOTUs ----
  ## prepare data
OTUs_perMT_nrOTUs <- means %>%
  dplyr::select(!all_of(responsive_OTUs)) %>%
  transmute(OTU_sum = rowSums(.)) %>%
  rownames_to_column("Sample_ID") 

GLAM_data_nrOTUs <- merge(env2, OTUs_perMT_nrOTUs, by = "Sample_ID")
GLAM_data_nrOTUs$Sampling_Date <- as.Date(GLAM_data_nrOTUs$Sampling_Date, format = "%d/%m/%y")

GLAM_data_nrOTUs <- GLAM_data_nrOTUs %>%
  mutate(Doy = yday(Sampling_Date)) %>% 
  mutate(Week = week(Sampling_Date)) %>%
  mutate(Ind_per_OTU = Number_of_Individuals / OTU_sum) %>%
  mutate(OTUs_per_Ind = OTU_sum / Number_of_Individuals)

GLAM_data_nrOTUs %>%
  filter(Ind_per_OTU < 1) # No sample with more OTUs than individuals 

GLAM_data_nrOTUs_Ctrl <- GLAM_data_nrOTUs %>%
  filter(Treatment == "Ctrl")

GLAM_data_nrOTUs_Bti <- GLAM_data_nrOTUs %>%
  filter(Treatment == "Bti")

#### Descr. statistics

Sum_ds_nrOTUs <- as.data.frame(matrix(c(mean(GLAM_data_nrOTUs_Ctrl$OTU_sum), mean(GLAM_data_nrOTUs_Bti$OTU_sum),
                                        sd(GLAM_data_nrOTUs_Ctrl$OTU_sum), sd(GLAM_data_nrOTUs_Bti$OTU_sum)), ncol = 2, byrow = TRUE)) %>% rename(Ctrl = V1,  Bti = V2) %>% mutate(Treatment = c("mean", "sd")) %>%
  column_to_rownames("Treatment")

#### ANOVA
aov_sum_nrOTUs <- aov(OTU_sum ~ Treatment * as.factor(Week), data = GLAM_data_nrOTUs)
summary(aov_sum_nrOTUs)
plot(aov_sum_nrOTUs)

aov_sum_effectsize_nrOTUs <- effectsize::eta_squared(aov_sum_nrOTUs)

#### Ds per week
  ## Mean
Sum_nrOTUs_mean_Bti <- lapply(c(16:31), function(x)
  mean(GLAM_data_nrOTUs_Bti[GLAM_data_nrOTUs_Bti$Week == x, ]$OTU_sum)
)
Sum_nrOTUs_mean_Ctrl <- lapply(c(16:31), function(x)
  mean(GLAM_data_nrOTUs_Ctrl[GLAM_data_nrOTUs_Ctrl$Week == x, ]$OTU_sum)
)

names(Sum_nrOTUs_mean_Bti) <- c(16:31)
names(Sum_nrOTUs_mean_Ctrl) <- c(16:31)

Sum_nrOTUs_mean_Bti_df <- as.data.frame(do.call(rbind, Sum_nrOTUs_mean_Bti))
Sum_nrOTUs_mean_Ctrl_df <- as.data.frame(do.call(rbind, Sum_nrOTUs_mean_Ctrl))

Sum_nrOTUs_mean <- cbind(Sum_nrOTUs_mean_Ctrl_df, mean_Bti = Sum_nrOTUs_mean_Bti_df$V1) %>% rename(mean_Ctrl = V1)

  ## Standard deviation
Sum_nrOTUs_sd_Bti <- lapply(c(16:31), function(x)
  sd(GLAM_data_nrOTUs_Bti[GLAM_data_nrOTUs_Bti$Week == x, ]$OTU_sum)
)
Sum_nrOTUs_sd_Ctrl <- lapply(c(16:31), function(x)
  sd(GLAM_data_nrOTUs_Ctrl[GLAM_data_nrOTUs_Ctrl$Week == x, ]$OTU_sum)
)

names(Sum_nrOTUs_sd_Bti) <- c(16:31)
names(Sum_nrOTUs_sd_Ctrl) <- c(16:31)

Sum_nrOTUs_sd_Bti_df <- as.data.frame(do.call(rbind, Sum_nrOTUs_sd_Bti))
Sum_nrOTUs_sd_Ctrl_df <- as.data.frame(do.call(rbind, Sum_nrOTUs_sd_Ctrl))

Sum_nrOTUs_sd <- cbind(Sum_nrOTUs_sd_Ctrl_df, sd_Bti = Sum_nrOTUs_sd_Bti_df$V1) %>% rename(sd_Ctrl = V1)

Sum_nrOTUs_ds <- cbind(Sum_nrOTUs_mean, Sum_nrOTUs_sd) %>% relocate(sd_Ctrl, .after = mean_Ctrl)

#### ANOVA per week

Sum_nrOTUs_aov <- lapply(c(16:31), function(x) 
  aov(OTU_sum ~ Treatment, data = GLAM_data_nrOTUs[GLAM_data_nrOTUs$Week == x, ])) 

Sum_nrOTUs_anovas <- lapply(c(16:31), function(x) 
  summary(aov(OTU_sum ~ Treatment, data = GLAM_data_nrOTUs[GLAM_data_nrOTUs$Week == x, ]))) 

plot_Sum_nrOTUs_aov <- lapply(Sum_nrOTUs_aov, function(x)
  plot(x)
) # hit Esc if you can't stop plotting ;)

Sum_nrOTUs_aov_effectsize <- lapply(c(1:16), function(x)
  effectsize::eta_squared(Sum_nrOTUs_aov[[x]])
)

#### One-way ANOVA
owANOVA_sum_nrOTUs <- aov(OTU_sum ~ Treatment, data = GLAM_data_nrOTUs)
summary(owANOVA_sum_nrOTUs)
plot(owANOVA_sum_nrOTUs)

### Sum of rOTUs ----
  ## prepare data
OTUs_perMT_rOTUs <- means %>%
  dplyr::select(all_of(responsive_OTUs)) %>%
  transmute(OTU_sum = rowSums(.)) %>%
  filter(OTU_sum != 0) %>%
  rownames_to_column("Sample_ID") 

GLAM_data_rOTUs <- merge(env2, OTUs_perMT_rOTUs, by = "Sample_ID")
GLAM_data_rOTUs$Sampling_Date <- as.Date(GLAM_data_rOTUs$Sampling_Date, format = "%d/%m/%y")

GLAM_data_rOTUs <- GLAM_data_rOTUs %>%
  mutate(Doy = yday(Sampling_Date)) %>% 
  mutate(Week = week(Sampling_Date)) %>%
  mutate(Ind_per_OTU = Number_of_Individuals / OTU_sum) %>%
  mutate(OTUs_per_Ind = OTU_sum / Number_of_Individuals) 

GLAM_data_rOTUs %>%
  filter(Ind_per_OTU < 1) # No samples with more OTUs than individuals

GLAM_data_rOTUs_Ctrl <- GLAM_data_rOTUs %>%
  filter(Treatment == "Ctrl")

GLAM_data_rOTUs_Bti <- GLAM_data_rOTUs %>%
  filter(Treatment == "Bti")


#### Descr. statistics
Sum_ds_rOTUs <- as.data.frame(matrix(c(mean(GLAM_data_rOTUs_Ctrl$OTU_sum), mean(GLAM_data_rOTUs_Bti$OTU_sum),
                                       sd(GLAM_data_rOTUs_Ctrl$OTU_sum), sd(GLAM_data_rOTUs_Bti$OTU_sum)), ncol = 2, byrow = TRUE)) %>% rename(Ctrl = V1,  Bti = V2) %>% mutate(Treatment = c("mean", "sd")) %>%
  column_to_rownames("Treatment")

#### ANOVA
aov_sum_rOTUs <- aov(OTU_sum ~ Treatment * as.factor(Week), data = GLAM_data_rOTUs)
summary(aov_sum_rOTUs)
plot(aov_sum_rOTUs)

aov_sum_effectsize_rOTUs <- effectsize::eta_squared(aov_sum_rOTUs)

#### Ds per week
  ## Mean
Sum_rOTUs_mean_Bti <- lapply(c(16:31), function(x)
  mean(GLAM_data_rOTUs_Bti[GLAM_data_rOTUs_Bti$Week == x, ]$OTU_sum)
)
Sum_rOTUs_mean_Ctrl <- lapply(c(16:31), function(x)
  mean(GLAM_data_rOTUs_Ctrl[GLAM_data_rOTUs_Ctrl$Week == x, ]$OTU_sum)
)

names(Sum_rOTUs_mean_Bti) <- c(16:31)
names(Sum_rOTUs_mean_Ctrl) <- c(16:31)

Sum_rOTUs_mean_Bti_df <- as.data.frame(do.call(rbind, Sum_rOTUs_mean_Bti))
Sum_rOTUs_mean_Ctrl_df <- as.data.frame(do.call(rbind, Sum_rOTUs_mean_Ctrl))

Sum_rOTUs_mean <- cbind(Sum_rOTUs_mean_Ctrl_df, mean_Bti = Sum_rOTUs_mean_Bti_df$V1) %>% rename(mean_Ctrl = V1)

  ## Standard deviation
Sum_rOTUs_sd_Bti <- lapply(c(16:31), function(x)
  sd(GLAM_data_rOTUs_Bti[GLAM_data_rOTUs_Bti$Week == x, ]$OTU_sum)
)
Sum_rOTUs_sd_Ctrl <- lapply(c(16:31), function(x)
  sd(GLAM_data_rOTUs_Ctrl[GLAM_data_rOTUs_Ctrl$Week == x, ]$OTU_sum)
)

names(Sum_rOTUs_sd_Bti) <- c(16:31)
names(Sum_rOTUs_sd_Ctrl) <- c(16:31)

Sum_rOTUs_sd_Bti_df <- as.data.frame(do.call(rbind, Sum_rOTUs_sd_Bti))
Sum_rOTUs_sd_Ctrl_df <- as.data.frame(do.call(rbind, Sum_rOTUs_sd_Ctrl))

Sum_rOTUs_sd <- cbind(Sum_rOTUs_sd_Ctrl_df, sd_Bti = Sum_rOTUs_sd_Bti_df$V1) %>% rename(sd_Ctrl = V1)

Sum_rOTUs_ds <- cbind(Sum_rOTUs_mean, Sum_rOTUs_sd) %>% relocate(sd_Ctrl, .after = mean_Ctrl)

#### ANOVA per week
Sum_rOTUs_aov <- lapply(c(16:31), function(x) 
  aov(OTU_sum ~ Treatment, data = GLAM_data_rOTUs[GLAM_data_rOTUs$Week == x, ])) 

Sum_rOTUs_anovas <- lapply(c(16:31), function(x) 
  summary(aov(OTU_sum ~ Treatment, data = GLAM_data_rOTUs[GLAM_data_rOTUs$Week == x, ]))) 

plot_Sum_rOTUs_aov <- lapply(Sum_rOTUs_aov, function(x)
  plot(x)
) # hit Esc if you can't stop plotting ;)

Sum_rOTUs_aov_effectsize <- lapply(c(1:16), function(x)
  effectsize::eta_squared(Sum_rOTUs_aov[[x]])
)

#### One-way ANOVA
owANOVA_sum_rOTUs <- aov(OTU_sum ~ Treatment, data = GLAM_data_rOTUs)
summary(owANOVA_sum_rOTUs)
plot(owANOVA_sum_rOTUs)

### Cumsum of OTUs ----
  ## prepare data
means_env <- means %>%
  rownames_to_column("Sample_ID")
means_env <- merge(means_env, env2, by.y = "Sample_ID")
means_env <- means_env %>%
  dplyr::select(-Number_of_Individuals) %>%
  column_to_rownames("Sample_ID") 
means_env <- melt(means_env, id.vars = c('Treatment', 'Doy', 'Week', 'Sampling_Date', 'Pond'))
means_env <- means_env %>%
  mutate(Week = as.numeric(as.character(Week)))%>%
  mutate(Doy = as.numeric(as.character(Doy)))

PondDate <- env2 %>%
  dplyr::select(Sampling_Date, Pond) %>%
  arrange(Pond, Sampling_Date) 

Cumsum_data <- means_env %>%
  arrange(variable, Pond, Sampling_Date) %>%
  filter(value != 0) %>%
  mutate(variable = as.character(variable)) %>%
  group_by(Pond) %>%
  filter(!duplicated(variable)) %>%
  arrange(Pond, Sampling_Date) %>%
  dplyr::mutate(OTU_sum = cumsum(value)) %>%
  filter(!duplicated(Sampling_Date, fromLast = TRUE)) %>%
  mutate(Doy = as.character(Doy)) %>%
  mutate(Doy = as.numeric(Doy)) %>%
  mutate(Week = as.character(Week)) %>%
  mutate(Week = as.numeric(Week))

Cumsum_data <- merge(PondDate, Cumsum_data, by = c("Sampling_Date", "Pond"), all = TRUE) %>%
  arrange(Pond, Sampling_Date) %>%
  mutate(OTU_sum = ifelse(Sampling_Date == "2020-04-14" & is.na(OTU_sum), 0, OTU_sum)) %>%
  fill(OTU_sum, .direction = "down") %>%
  mutate(Treatment = ifelse(Pond == "Pond1" | Pond == "Pond3" | Pond == "Pond5" | Pond == "Pond7" | 
                              Pond == "Pond9" | Pond == "Pond11", "Ctrl", "Bti")) %>%
  mutate(Treatment = as.factor(Treatment)) %>%
  dplyr::select(Sampling_Date, Treatment, Pond, OTU_sum) %>% 
  mutate(Doy = as.numeric(yday(Sampling_Date))) %>%
  mutate(Week = as.numeric(week(Sampling_Date)))

Cumsum_data$Treatment <- factor(Cumsum_data$Treatment, levels = c("Ctrl", "Bti"))

Cumsum_data_Ctrl <- Cumsum_data[Cumsum_data$Treatment == "Ctrl", ]
Cumsum_data_Bti <- Cumsum_data[Cumsum_data$Treatment == "Bti", ]

#### Descr. statistics
Cumsum_ds_all <- as.data.frame(matrix(c(mean(Cumsum_data_Ctrl$OTU_sum), mean(Cumsum_data_Bti$OTU_sum),
                                        sd(Cumsum_data_Ctrl$OTU_sum), sd(Cumsum_data_Bti$OTU_sum)), ncol = 2, byrow = TRUE)) %>% rename(Ctrl = V1,  Bti = V2) %>% mutate(Treatment = c("mean", "sd")) %>%
  column_to_rownames("Treatment")

#### Two-way ANOVA
aov_sum <- aov(OTU_sum ~ Treatment * as.factor(Week), data = Cumsum_data)
summary(aov_sum)
plot(aov_sum)

aov_Cumsum_effectsize_all <- effectsize::eta_squared(aov_sum)

#### Ds per week
  ## Mean
Cumsum_mean_Bti <- lapply(c(16:31), function(x)
  mean(Cumsum_data_Bti[Cumsum_data_Bti$Week == x, ]$OTU_sum)
)
Cumsum_mean_Ctrl <- lapply(c(16:31), function(x)
  mean(Cumsum_data_Ctrl[Cumsum_data_Ctrl$Week == x, ]$OTU_sum)
)

names(Cumsum_mean_Bti) <- c(16:31)
names(Cumsum_mean_Ctrl) <- c(16:31)

Cumsum_mean_Bti_df <- as.data.frame(do.call(rbind, Cumsum_mean_Bti))
Cumsum_mean_Ctrl_df <- as.data.frame(do.call(rbind, Cumsum_mean_Ctrl))

Cumsum_mean <- cbind(Cumsum_mean_Ctrl_df, mean_Bti = Cumsum_mean_Bti_df$V1) %>% rename(mean_Ctrl = V1)

  ## Standard deviation
Cumsum_sd_Bti <- lapply(c(16:31), function(x)
  sd(Cumsum_data_Bti[Cumsum_data_Bti$Week == x, ]$OTU_sum)
)
Cumsum_sd_Ctrl <- lapply(c(16:31), function(x)
  sd(Cumsum_data_Ctrl[Cumsum_data_Ctrl$Week == x, ]$OTU_sum)
)

names(Cumsum_sd_Bti) <- c(16:31)
names(Cumsum_sd_Ctrl) <- c(16:31)

Cumsum_sd_Bti_df <- as.data.frame(do.call(rbind, Cumsum_sd_Bti))
Cumsum_sd_Ctrl_df <- as.data.frame(do.call(rbind, Cumsum_sd_Ctrl))

Cumsum_sd <- cbind(Cumsum_sd_Ctrl_df, sd_Bti = Cumsum_sd_Bti_df$V1) %>% 
  rename(sd_Ctrl = V1)

Cumsum_ds <- cbind(Cumsum_mean, Cumsum_sd) %>% relocate(sd_Ctrl, .after = mean_Ctrl)

#### ANOVA per week
Cumsum_aov <- lapply(c(16:31), function(x) 
  aov(OTU_sum ~ Treatment, data = Cumsum_data[Cumsum_data$Week == x, ])) 

Cumsum_anovas <- lapply(c(16:31), function(x) 
  summary(aov(OTU_sum ~ Treatment, data = Cumsum_data[Cumsum_data$Week == x, ]))) 

plot_Cumsum_aov <- lapply(Cumsum_aov, function(x)
  plot(x)
) # hit Esc if you can't stop plotting ;)

Cumsum_aov_effectsize <- lapply(c(1:16), function(x)
  effectsize::eta_squared(Cumsum_aov[[x]])
)


### Cumsum of nrOTUs ----
  ## prepare data
Cumsum_data_nrOTUs <- means_env[!means_env$variable %in% responsive_OTUs, ] %>%
  arrange(variable, Pond, Sampling_Date) %>%
  filter(value != 0) %>%
  mutate(variable = as.character(variable)) %>%
  group_by(Pond) %>%
  filter(!duplicated(variable)) %>%
  arrange(Pond, Sampling_Date) %>%
  dplyr::mutate(OTU_sum = cumsum(value)) %>%
  filter(!duplicated(Sampling_Date, fromLast = TRUE)) %>%
  mutate(Doy = as.character(Doy)) %>%
  mutate(Doy = as.numeric(Doy)) %>%
  mutate(Week = as.character(Week)) %>%
  mutate(Week = as.numeric(Week))

Cumsum_data_nrOTUs <- merge(PondDate, Cumsum_data_nrOTUs, by = c("Sampling_Date", "Pond"), all = TRUE) %>%
  arrange(Pond, Sampling_Date) %>%
  mutate(OTU_sum = ifelse(Sampling_Date == "2020-04-14" & is.na(OTU_sum), 0, OTU_sum)) %>%
  fill(OTU_sum, .direction = "down") %>%
  mutate(Treatment = ifelse(Pond == "Pond1" | Pond == "Pond3" | Pond == "Pond5" | Pond == "Pond7" | 
                              Pond == "Pond9" | Pond == "Pond11", "Ctrl", "Bti")) %>%
  mutate(Treatment = as.factor(Treatment)) %>%
  dplyr::select(Sampling_Date, Treatment, Pond, OTU_sum) %>%
  mutate(Doy = as.numeric(yday(Sampling_Date))) %>%
  mutate(Week = as.numeric(week(Sampling_Date)))


Cumsum_data_nrOTUs_Ctrl <- Cumsum_data_nrOTUs[Cumsum_data_nrOTUs$Treatment == "Ctrl", ]
Cumsum_data_nrOTUs_Bti <- Cumsum_data_nrOTUs[Cumsum_data_nrOTUs$Treatment == "Bti", ]

#### Descr. statistics
Cumsum_ds_nrOTUs <- as.data.frame(matrix(c(mean(Cumsum_data_nrOTUs_Ctrl$OTU_sum), mean(Cumsum_data_nrOTUs_Bti$OTU_sum),
                                           sd(Cumsum_data_nrOTUs_Ctrl$OTU_sum), sd(Cumsum_data_nrOTUs_Bti$OTU_sum)), ncol = 2, byrow = TRUE)) %>% rename(Ctrl = V1,  Bti = V2) %>% mutate(Treatment = c("mean", "sd")) %>%
  column_to_rownames("Treatment")

#### ANOVA
aov_Cumsum_nrOTUs <- aov(OTU_sum ~ Treatment * as.factor(Week), data = Cumsum_data_nrOTUs)
summary(aov_Cumsum_nrOTUs)
plot(aov_Cumsum_nrOTUs)

aov_Cumsum_effectsize_nrOTUs <- effectsize::eta_squared(aov_Cumsum_nrOTUs)


#### Ds per week
  ## Mean
Cumsum_nrOTUs_mean_Bti <- lapply(c(16:31), function(x)
  mean(Cumsum_data_nrOTUs_Bti[Cumsum_data_nrOTUs_Bti$Week == x, ]$OTU_sum)
)
Cumsum_nrOTUs_mean_Ctrl <- lapply(c(16:31), function(x)
  mean(Cumsum_data_nrOTUs_Ctrl[Cumsum_data_nrOTUs_Ctrl$Week == x, ]$OTU_sum)
)

names(Cumsum_nrOTUs_mean_Bti) <- c(16:31)
names(Cumsum_nrOTUs_mean_Ctrl) <- c(16:31)

Cumsum_nrOTUs_mean_Bti_df <- as.data.frame(do.call(rbind, Cumsum_nrOTUs_mean_Bti))
Cumsum_nrOTUs_mean_Ctrl_df <- as.data.frame(do.call(rbind, Cumsum_nrOTUs_mean_Ctrl))

Cumsum_nrOTUs_mean <- cbind(Cumsum_nrOTUs_mean_Ctrl_df, mean_Bti = Cumsum_nrOTUs_mean_Bti_df$V1) %>% rename(mean_Ctrl = V1)

  ## Standard deviation
Cumsum_nrOTUs_sd_Bti <- lapply(c(16:31), function(x)
  sd(Cumsum_data_nrOTUs_Bti[Cumsum_data_nrOTUs_Bti$Week == x, ]$OTU_sum)
)
Cumsum_nrOTUs_sd_Ctrl <- lapply(c(16:31), function(x)
  sd(Cumsum_data_nrOTUs_Ctrl[Cumsum_data_nrOTUs_Ctrl$Week == x, ]$OTU_sum)
)

names(Cumsum_nrOTUs_sd_Bti) <- c(16:31)
names(Cumsum_nrOTUs_sd_Ctrl) <- c(16:31)

Cumsum_nrOTUs_sd_Bti_df <- as.data.frame(do.call(rbind, Cumsum_nrOTUs_sd_Bti))
Cumsum_nrOTUs_sd_Ctrl_df <- as.data.frame(do.call(rbind, Cumsum_nrOTUs_sd_Ctrl))

Cumsum_nrOTUs_sd <- cbind(Cumsum_nrOTUs_sd_Ctrl_df, sd_Bti = Cumsum_nrOTUs_sd_Bti_df$V1) %>% rename(sd_Ctrl = V1)

Cumsum_nrOTUs_ds <- cbind(Cumsum_nrOTUs_mean, Cumsum_nrOTUs_sd) %>% relocate(sd_Ctrl, .after = mean_Ctrl)

#### ANOVA per week
Cumsum_nrOTUs_aov <- lapply(c(16:31), function(x) 
  aov(OTU_sum ~ Treatment, data = Cumsum_data_nrOTUs[Cumsum_data_nrOTUs$Week == x, ])) 

Cumsum_nrOTUs_anovas <- lapply(c(16:31), function(x) 
  summary(aov(OTU_sum ~ Treatment, data = Cumsum_data_nrOTUs[Cumsum_data_nrOTUs$Week == x, ]))) 

plot_Cumsum_nrOTUs_aov <- lapply(Cumsum_nrOTUs_aov, function(x)
  plot(x)
) # hit Esc if you can't stop plotting ;)

Cumsum_nrOTUs_aov_effectsize <- lapply(c(1:16), function(x)
  effectsize::eta_squared(Cumsum_nrOTUs_aov[[x]])
)

## Final cumsum of non-responsive OTUs
Cumsum_data_nrOTUs_max <- Cumsum_data_nrOTUs %>%
  group_by(Pond) %>%
  summarise(No_of_OTUs = max(OTU_sum)) %>%
  mutate(Treatment = ifelse(Pond == "Pond1" | Pond == "Pond3" | Pond == "Pond5" | Pond == "Pond7" | 
                              Pond == "Pond9" | Pond == "Pond11", "Ctrl", "Bti")) %>%
  mutate(Treatment = as.factor(Treatment)) 

aov_max_nrOTUs <- summary(aov(No_of_OTUs ~ Treatment, data = Cumsum_data_nrOTUs_max))
plot(aov_max_nrOTUs)
max_nrOTUs_es <- effectsize::eta_squared(aov(No_of_OTUs ~ Treatment, data = Cumsum_data_nrOTUs_max), partial = TRUE)


### Cumsum of rOTUs ----
Cumsum_data_rOTUs <- means_env[means_env$variable %in% responsive_OTUs, ] %>%
  arrange(variable, Pond, Sampling_Date) %>%
  filter(value != 0) %>%
  mutate(variable = as.character(variable)) %>%
  group_by(Pond) %>%
  filter(!duplicated(variable)) %>%
  arrange(Pond, Sampling_Date) %>%
  dplyr::mutate(OTU_sum = cumsum(value)) %>%
  filter(!duplicated(Sampling_Date, fromLast = TRUE)) %>%
  mutate(Doy = as.character(Doy)) %>%
  mutate(Doy = as.numeric(Doy)) %>%
  mutate(Week = as.character(Week)) %>%
  mutate(Week = as.numeric(Week))

Cumsum_data_rOTUs <- merge(PondDate, Cumsum_data_rOTUs, by = c("Sampling_Date", "Pond"), all = TRUE) %>%
  arrange(Pond, Sampling_Date) %>%
  mutate(OTU_sum = ifelse(Sampling_Date == "2020-04-14" & is.na(OTU_sum), 0, OTU_sum)) %>%
  fill(OTU_sum, .direction = "down") %>%
  mutate(Treatment = ifelse(Pond == "Pond1" | Pond == "Pond3" | Pond == "Pond5" | Pond == "Pond7" | 
                              Pond == "Pond9" | Pond == "Pond11", "Ctrl", "Bti")) %>%
  mutate(Treatment = as.factor(Treatment)) %>%
  dplyr::select(Sampling_Date, Treatment, Pond, OTU_sum) %>%
  mutate(Doy = as.numeric(yday(Sampling_Date))) %>%
  mutate(Week = as.numeric(week(Sampling_Date))) %>%
  replace(is.na(.), 0)

Cumsum_data_rOTUs_Ctrl <- Cumsum_data_rOTUs[Cumsum_data_rOTUs$Treatment == "Ctrl", ]
Cumsum_data_rOTUs_Bti <- Cumsum_data_rOTUs[Cumsum_data_rOTUs$Treatment == "Bti", ]

#### Descr. statistics
Cumsum_ds_rOTUs <- as.data.frame(matrix(c(mean(Cumsum_data_rOTUs_Ctrl$OTU_sum), mean(Cumsum_data_rOTUs_Bti$OTU_sum),
                                          sd(Cumsum_data_rOTUs_Ctrl$OTU_sum), sd(Cumsum_data_rOTUs_Bti$OTU_sum)), ncol = 2, byrow = TRUE)) %>% rename(Ctrl = V1,  Bti = V2) %>% mutate(Treatment = c("mean", "sd")) %>%
  column_to_rownames("Treatment")


#### ANOVA
aov_Cumsum_rOTUs <- aov(OTU_sum ~ Treatment * as.factor(Week), data = Cumsum_data_rOTUs)
summary(aov_Cumsum_rOTUs)
plot(aov_Cumsum_rOTUs)

aov_Cumsum_effectsize_rOTUs <- effectsize::eta_squared(aov_Cumsum_rOTUs)

#### Ds per week
  ## Mean
Cumsum_rOTUs_mean_Bti <- lapply(c(16:31), function(x)
  mean(Cumsum_data_rOTUs_Bti[Cumsum_data_rOTUs_Bti$Week == x, ]$OTU_sum)
)
Cumsum_rOTUs_mean_Ctrl <- lapply(c(16:31), function(x)
  mean(Cumsum_data_rOTUs_Ctrl[Cumsum_data_rOTUs_Ctrl$Week == x, ]$OTU_sum)
)

names(Cumsum_rOTUs_mean_Bti) <- c(16:31)
names(Cumsum_rOTUs_mean_Ctrl) <- c(16:31)

Cumsum_rOTUs_mean_Bti_df <- as.data.frame(do.call(rbind, Cumsum_rOTUs_mean_Bti))
Cumsum_rOTUs_mean_Ctrl_df <- as.data.frame(do.call(rbind, Cumsum_rOTUs_mean_Ctrl))

Cumsum_rOTUs_mean <- cbind(Cumsum_rOTUs_mean_Ctrl_df, mean_Bti = Cumsum_rOTUs_mean_Bti_df$V1) %>% rename(mean_Ctrl = V1)

  ## Standard deviation
Cumsum_rOTUs_sd_Bti <- lapply(c(16:31), function(x)
  sd(Cumsum_data_rOTUs_Bti[Cumsum_data_rOTUs_Bti$Week == x, ]$OTU_sum)
)
Cumsum_rOTUs_sd_Ctrl <- lapply(c(16:31), function(x)
  sd(Cumsum_data_rOTUs_Ctrl[Cumsum_data_rOTUs_Ctrl$Week == x, ]$OTU_sum)
)

names(Cumsum_rOTUs_sd_Bti) <- c(16:31)
names(Cumsum_rOTUs_sd_Ctrl) <- c(16:31)

Cumsum_rOTUs_sd_Bti_df <- as.data.frame(do.call(rbind, Cumsum_rOTUs_sd_Bti))
Cumsum_rOTUs_sd_Ctrl_df <- as.data.frame(do.call(rbind, Cumsum_rOTUs_sd_Ctrl))

Cumsum_rOTUs_sd <- cbind(Cumsum_rOTUs_sd_Ctrl_df, sd_Bti = Cumsum_rOTUs_sd_Bti_df$V1) %>% rename(sd_Ctrl = V1)

Cumsum_rOTUs_ds <- cbind(Cumsum_rOTUs_mean, Cumsum_rOTUs_sd) %>% relocate(sd_Ctrl, .after = mean_Ctrl)

#### ANOVA per week
Cumsum_rOTUs_aov <- lapply(c(16:31), function(x) 
  aov(OTU_sum ~ Treatment, data = Cumsum_data_rOTUs[Cumsum_data_rOTUs$Week == x, ])) 

Cumsum_rOTUs_anovas <- lapply(c(16:31), function(x) 
  summary(aov(OTU_sum ~ Treatment, data = Cumsum_data_rOTUs[Cumsum_data_rOTUs$Week == x, ]))) 

plot_Cumsum_rOTUs_aov <- lapply(Cumsum_rOTUs_aov, function(x)
  plot(x)
) # hit Esc if you can't stop plotting ;)

Cumsum_rOTUs_aov_effectsize <- lapply(c(1:16), function(x)
  effectsize::eta_squared(Cumsum_rOTUs_aov[[x]])
)

#### One-way ANOVA
  ## Final cumsum of responsive OTUs
Cumsum_data_rOTUs_max <- Cumsum_data_rOTUs %>%
  group_by(Pond) %>%
  summarise(No_of_OTUs = max(OTU_sum)) %>%
  mutate(Treatment = ifelse(Pond == "Pond1" | Pond == "Pond3" | Pond == "Pond5" | Pond == "Pond7" | 
                              Pond == "Pond9" | Pond == "Pond11", "Ctrl", "Bti")) %>%
  mutate(Treatment = as.factor(Treatment)) 

aov_max_rOTUs <- summary(aov(No_of_OTUs ~ Treatment, data = Cumsum_data_rOTUs_max))
plot(aov_max_rOTUs)
max_rOTUs_es <- effectsize::eta_squared(aov(No_of_OTUs ~ Treatment, data = Cumsum_data_rOTUs_max), partial = TRUE)


### Comparison of slopes ----
  ## cumsum as linear model
lm <- lm(OTU_sum ~ Week * Treatment, data =  Cumsum_data[Cumsum_data$Week < 23, ])
anova(lm)
plot(lm)
vlm <- emmeans::emtrends(lm, specs = pairwise ~ Treatment, var = "Week")
svlm <- summary(vlm, infer = TRUE)
effectsize::eta_squared(lm, partial = FALSE)

lm_nrOTUs <- lm(OTU_sum ~ Week * Treatment, data =  Cumsum_data_nrOTUs[Cumsum_data_nrOTUs$Week < 23, ]) 
anova(lm_nrOTUs)
plot(lm_nrOTUs)
vlm_nrOTUs <- emmeans::emtrends(lm_nrOTUs, specs = pairwise ~ Treatment, var = "Week")
svlm_nrOTUs <- summary(vlm_nrOTUs, infer = TRUE)
effectsize::eta_squared(lm_nrOTUs, partial = FALSE)

lm_rOTUs <- lm(OTU_sum ~ Week * Treatment, data =  Cumsum_data_rOTUs[Cumsum_data_rOTUs$Week < 23, ]) 
anova(lm_rOTUs)
plot(lm_rOTUs)
vlm_rOTUs <- emmeans::emtrends(lm_rOTUs, specs = pairwise ~ Treatment, var = "Week")
svlm_rOTUs <- summary(vlm_rOTUs, infer = TRUE)
effectsize::eta_squared(lm_rOTUs, partial = FALSE)

### Number of individuals ----
#### Descr. statistics
NoInd_ds <- as.data.frame(matrix(c(mean(GLAM_data_Ctrl$Number_of_Individuals), mean(GLAM_data_Bti$Number_of_Individuals),
                                   sd(GLAM_data_Ctrl$Number_of_Individuals), sd(GLAM_data_Bti$Number_of_Individuals)), ncol = 2, byrow = TRUE)) %>% rename(Ctrl = V1,  Bti = V2) %>% mutate(Treatment = c("mean", "sd")) %>%
  column_to_rownames("Treatment")

#### ANOVA
aov_NoInd <- aov(Number_of_Individuals ~ Treatment * as.factor(Week), data = GLAM_data)
summary(aov_NoInd)
aov_NoInd_es <- effectsize::eta_squared(aov_NoInd, partial = TRUE)
plot(aov_NoInd)

#### Ds per week
  ## Mean
Ind_mean_Bti <- lapply(c(16:31), function(x)
  mean(GLAM_data_Bti[GLAM_data_Bti$Week == x, ]$Number_of_Individuals)
)
Ind_mean_Ctrl <- lapply(c(16:31), function(x)
  mean(GLAM_data_Ctrl[GLAM_data_Ctrl$Week == x, ]$Number_of_Individuals)
)

names(Ind_mean_Bti) <- c(16:31)
names(Ind_mean_Ctrl) <- c(16:31)

Ind_mean_Bti_df <- as.data.frame(do.call(rbind, Ind_mean_Bti))
Ind_mean_Ctrl_df <- as.data.frame(do.call(rbind, Ind_mean_Ctrl))

Ind_mean <- cbind(Ind_mean_Ctrl_df, mean_Bti = Ind_mean_Bti_df$V1) %>% rename(mean_Ctrl = V1)

  ## Standard deviation
Ind_sd_Bti <- lapply(c(16:31), function(x)
  sd(GLAM_data_Bti[GLAM_data_Bti$Week == x, ]$Number_of_Individuals)
)
Ind_sd_Ctrl <- lapply(c(16:31), function(x)
  sd(GLAM_data_Ctrl[GLAM_data_Ctrl$Week == x, ]$Number_of_Individuals)
)

names(Ind_sd_Bti) <- c(16:31)
names(Ind_sd_Ctrl) <- c(16:31)

Ind_sd_Bti_df <- as.data.frame(do.call(rbind, Ind_sd_Bti))
Ind_sd_Ctrl_df <- as.data.frame(do.call(rbind, Ind_sd_Ctrl))

Ind_sd <- cbind(Ind_sd_Ctrl_df, sd_Bti = Ind_sd_Bti_df$V1) %>% 
  rename(sd_Ctrl = V1)

Ind_ds <- cbind(Ind_mean, Ind_sd) %>% relocate(sd_Ctrl, .after = mean_Ctrl)

  ### ANOVA per week
Ind_aov <- lapply(c(16:31), function(x) 
  aov(Number_of_Individuals ~ Treatment, data = GLAM_data[GLAM_data$Week == x, ]))

Ind_anovas <- lapply(c(16:31), function(x) 
  summary(aov(Number_of_Individuals ~ Treatment, data = GLAM_data[GLAM_data$Week == x, ])))

names(Ind_anovas) <- c(16:31)

plot_aov_Ind <- lapply(Ind_aov, function(x)
 plot(x)
) 

Ind_aov_effectsize <- lapply(c(1:16), function(x)
  effectsize::eta_squared(Ind_aov[[x]])
)


### Individuals per OTU ----

#### Descr. statistics
IndOTU_ds <- as.data.frame(matrix(c(mean(GLAM_data_Ctrl$Ind_per_OTU), mean(GLAM_data_Bti$Ind_per_OTU),
                                    sd(GLAM_data_Ctrl$Ind_per_OTU), sd(GLAM_data_Bti$Ind_per_OTU)), ncol = 2, byrow = TRUE)) %>% rename(Ctrl = V1,  Bti = V2) %>% mutate(Treatment = c("mean", "sd")) %>%
  column_to_rownames("Treatment")

### ANOVA
aov_IndOTU <- aov(Ind_per_OTU ~ Treatment * as.factor(Week), data = GLAM_data)
summary(aov_IndOTU)
plot(aov_IndOTU)
aov_IndOTU_es <- effectsize::eta_squared(aov_IndOTU, partial = TRUE)

### Ds per week
  ## Mean
IndOTU_mean_Bti <- lapply(c(16:31), function(x)
  mean(GLAM_data_Bti[GLAM_data_Bti$Week == x, ]$Ind_per_OTU)
)
IndOTU_mean_Ctrl <- lapply(c(16:31), function(x)
  mean(GLAM_data_Ctrl[GLAM_data_Ctrl$Week == x, ]$Ind_per_OTU)
)

names(IndOTU_mean_Bti) <- c(16:31)
names(IndOTU_mean_Ctrl) <- c(16:31)

IndOTU_mean_Bti_df <- as.data.frame(do.call(rbind, IndOTU_mean_Bti))
IndOTU_mean_Ctrl_df <- as.data.frame(do.call(rbind, IndOTU_mean_Ctrl))

IndOTU_mean <- cbind(IndOTU_mean_Ctrl_df, mean_Bti = IndOTU_mean_Bti_df$V1) %>% rename(mean_Ctrl = V1)

  ## Standard deviation
IndOTU_sd_Bti <- lapply(c(16:31), function(x)
  sd(GLAM_data_Bti[GLAM_data_Bti$Week == x, ]$Ind_per_OTU)
)
IndOTU_sd_Ctrl <- lapply(c(16:31), function(x)
  sd(GLAM_data_Ctrl[GLAM_data_Ctrl$Week == x, ]$Ind_per_OTU)
)

names(IndOTU_sd_Bti) <- c(16:31)
names(IndOTU_sd_Ctrl) <- c(16:31)

IndOTU_sd_Bti_df <- as.data.frame(do.call(rbind, IndOTU_sd_Bti))
IndOTU_sd_Ctrl_df <- as.data.frame(do.call(rbind, IndOTU_sd_Ctrl))

IndOTU_sd <- cbind(IndOTU_sd_Ctrl_df, sd_Bti = IndOTU_sd_Bti_df$V1) %>% 
  rename(sd_Ctrl = V1)

IndOTU_ds <- cbind(IndOTU_mean, IndOTU_sd) %>% relocate(sd_Ctrl, .after = mean_Ctrl)

### ANOVA per week
IndOTU_aov <- lapply(c(16:31), function(x) 
  aov(Ind_per_OTU ~ Treatment, data = GLAM_data[GLAM_data$Week == x, ]))

IndOTU_anovas <- lapply(c(16:31), function(x) 
  summary(aov(Ind_per_OTU ~ Treatment, data = GLAM_data[GLAM_data$Week == x, ])))

names(IndOTU_anovas) <- c(16:31)

plot_aov_IndOTU <- lapply(IndOTU_aov, function(x)
  plot(x)
) 

IndOTU_aov_effectsize <- lapply(c(1:16), function(x)
  effectsize::eta_squared(aov_Ind[[x]])
)

### Duration of emergence ----
  ## preparation 
  ## Control data preparation 
OTUs_means <- as.vector(colnames(means))
means_env_Ctrl <- means_env %>%
  filter(Treatment == "Ctrl") 

Ctrl_OTUdynamic_list <- lapply(OTUs_means, function(x)
  
  means_env_Ctrl %>%
    filter(variable == x) %>%
    arrange(Pond, Sampling_Date) %>%
    mutate(value = as.numeric(value)) %>%
    group_by(Pond) %>%
    dplyr::mutate(cumsum = cumsum(value)) %>%
    filter(!duplicated(cumsum)) %>%
    filter(!cumsum == 0) %>%
    slice(which.min(cumsum), which.max(cumsum)) %>%
    ungroup() %>%
    mutate(Sampling_Date2 = Sampling_Date) %>%
    mutate(Doy2 = Doy) %>%
    mutate_at(c("Sampling_Date2", "Doy2", "cumsum"), funs(lead), n = 1) %>%
    filter(row_number() %% 2 != 0)  %>%
    mutate(Duration = as.Date(as.character(Sampling_Date2), format = "%Y-%m-%d")-
             as.Date(as.character(Sampling_Date), format = "%Y-%m-%d")) %>%
    dplyr::rename(Number_of_occurrences = cumsum) %>%
    dplyr::select(-value)
  
)

Ctrl_OTUdynamic <- do.call("rbind", Ctrl_OTUdynamic_list) 
Ctrl_OTUdynamic %>%
  summarise(count = n_distinct(variable))

  ## Bti data preparation
means_env_Bti <- means_env %>%
  filter(Treatment == "Bti") 

Bti_OTUdynamic_list <- lapply(OTUs_means, function(x)
  
  means_env_Bti %>%
    filter(variable == x) %>%
    arrange(Pond, Sampling_Date) %>%
    mutate(value = as.numeric(value)) %>%
    group_by(Pond) %>%
    dplyr::mutate(cumsum = cumsum(value)) %>%
    filter(!duplicated(cumsum)) %>%
    filter(!cumsum == 0) %>%
    slice(which.min(cumsum), which.max(cumsum)) %>%
    ungroup() %>%
    mutate(Sampling_Date2 = Sampling_Date) %>%
    mutate(Doy2 = Doy) %>%
    mutate_at(c("Sampling_Date2", "Doy2", "cumsum"), funs(lead), n = 1) %>%
    filter(row_number() %% 2 != 0)  %>%
    mutate(Duration = as.Date(as.character(Sampling_Date2), format = "%Y-%m-%d")-
             as.Date(as.character(Sampling_Date), format = "%Y-%m-%d"))  %>%
    dplyr::rename(Number_of_occurrences = cumsum) %>%
    dplyr::select(-value)
  
)


Bti_OTUdynamic <- do.call("rbind", Bti_OTUdynamic_list) 
Bti_OTUdynamic %>%
  summarise(count = n_distinct(variable))


OTU_dynamic <- rbind(Ctrl_OTUdynamic, Bti_OTUdynamic)

OTU_dynamic <- OTU_dynamic %>%
  mutate(Duration_num = as.numeric(Duration)) %>%
  mutate(variable = as.character(variable)) %>%
  mutate(Duration_num = ifelse(Duration_num == 0, 1, as.numeric(Duration_num)))

  ## OTU dynamic without OTUs occurring on first or last day of sampling
OTU_dynamic_complDuration <- OTU_dynamic[
  OTU_dynamic$Doy != 105 & 
    OTU_dynamic$Doy2 != 212, ] 

OTUs_cd <- OTU_dynamic_complDuration %>%
  group_by(variable) %>%
  summarise(count = n_distinct(Treatment)) %>%
  filter(!count == 1) %>%
  pull(variable) 

OTU_dynamic_cd <- OTU_dynamic_complDuration[
  OTU_dynamic_complDuration$variable %in% OTUs_cd, ]

OTU_dynamic_cd_Ctrl <- OTU_dynamic_cd %>%
  filter(Treatment == "Ctrl")

OTU_dynamic_cd_Bti <- OTU_dynamic_cd %>%
  filter(Treatment == "Bti")

#### Descr. statistics
Duration_ds <- as.data.frame(matrix(c(mean(OTU_dynamic_cd_Ctrl$Duration_num), mean(OTU_dynamic_cd_Bti$Duration_num),
                                      sd(OTU_dynamic_cd_Ctrl$Duration_num), sd(OTU_dynamic_cd_Bti$Duration_num)), ncol = 2, byrow = TRUE)) %>% rename(Ctrl = V1,  Bti = V2) %>% mutate(Treatment = c("mean", "sd")) %>%
  column_to_rownames("Treatment")

#### ANOVA
aov_Dur <- aov(Duration_num ~ Treatment, 
               data = OTU_dynamic_cd)
summary(aov_Dur)
plot(aov_Dur)

aov_Dur_es <- effectsize::eta_squared(aov_Dur, partial = TRUE)

#### Ds per OTU
  ## Mean
Dur_mean_Bti <- lapply(OTUs_cd, function(x)
  mean(OTU_dynamic_cd_Bti[OTU_dynamic_cd_Bti$variable == x, ]$Duration_num)
)
Dur_mean_Ctrl <- lapply(OTUs_cd, function(x)
  mean(OTU_dynamic_cd_Ctrl[OTU_dynamic_cd_Ctrl$variable == x, ]$Duration_num)
)

names(Dur_mean_Bti) <- OTUs_cd
names(Dur_mean_Ctrl) <- OTUs_cd

Dur_mean_Bti_df <- as.data.frame(do.call(rbind, Dur_mean_Bti))
Dur_mean_Ctrl_df <- as.data.frame(do.call(rbind, Dur_mean_Ctrl))

Dur_mean <- cbind(Dur_mean_Ctrl_df, mean_Bti = Dur_mean_Bti_df$V1) %>% rename(mean_Ctrl = V1)

  ## Standard deviation
Dur_sd_Bti <- lapply(OTUs_cd, function(x)
  sd(OTU_dynamic_cd_Bti[OTU_dynamic_cd_Bti$variable == x, ]$Duration_num)
)
Dur_sd_Ctrl <- lapply(OTUs_cd, function(x)
  sd(OTU_dynamic_cd_Ctrl[OTU_dynamic_cd_Ctrl$variable == x, ]$Duration_num)
)

names(Dur_sd_Bti) <- OTUs_cd
names(Dur_sd_Ctrl) <- OTUs_cd

Dur_sd_Bti_df <- as.data.frame(do.call(rbind, Dur_sd_Bti))
Dur_sd_Ctrl_df <- as.data.frame(do.call(rbind, Dur_sd_Ctrl))

Dur_sd <- cbind(Dur_sd_Ctrl_df, sd_Bti = Dur_sd_Bti_df$V1) %>% 
  rename(sd_Ctrl = V1)

Dur_ds <- cbind(Dur_mean, Dur_sd) %>% relocate(sd_Ctrl, .after = mean_Ctrl)

### ANOVA per OTU
OTUID_cd_aov <- Dur_ds %>% 
  na.omit() %>%
  rownames_to_column("OTUID_cd_aov") %>%
  filter(OTUID_cd_aov != "OTU_37" & OTUID_cd_aov != "OTU_357" & OTUID_cd_aov != "OTU_95") %>%
  pull(OTUID_cd_aov)

Duration_aovs <- lapply(OTUID_cd_aov, function(x) 
  aov(Duration_num ~ Treatment, 
      data = OTU_dynamic_cd[
        OTU_dynamic_cd$variable == x, ])
)

names(Duration_aovs) <- OTUID_cd_aov

Duration_anovas <- lapply(OTUID_cd_aov, function(x) 
  summary(aov(Duration_num ~ Treatment, 
              data = OTU_dynamic_cd[
                OTU_dynamic_cd$variable == x, ]))
) 

names(Duration_anovas) <- OTUID_cd_aov

plot_Duration_aovs <- lapply(Duration_aovs, function(x)
  plot(x)
)

Dur_aov_effectsize <- lapply(OTUID_cd_aov, function(x)
  effectsize::eta_squared(Duration_aovs[[x]])
)
names(Dur_aov_effectsize) <- OTUID_cd_aov

### Duration of emergence (nrOTUs) ----
#### Descr. statistics
OTU_dynamic_cd_nrOTUs <- OTU_dynamic_cd[!OTU_dynamic_cd$variable %in% responsive_OTUs, ]
OTU_dynamic_cd_nrOTUs_Ctrl <- OTU_dynamic_cd_nrOTUs %>% filter(Treatment == "Ctrl")
OTU_dynamic_cd_nrOTUs_Bti <- OTU_dynamic_cd_nrOTUs %>% filter(Treatment == "Bti")

Duration_ds_nrOTUs <- as.data.frame(matrix(c(mean(OTU_dynamic_cd_nrOTUs_Ctrl$Duration_num), mean(OTU_dynamic_cd_nrOTUs_Bti$Duration_num),
                                             sd(OTU_dynamic_cd_nrOTUs_Ctrl$Duration_num), sd(OTU_dynamic_cd_nrOTUs_Bti$Duration_num)), ncol = 2, byrow = TRUE)) %>% rename(Ctrl = V1,  Bti = V2) %>% mutate(Treatment = c("mean", "sd")) %>%
  column_to_rownames("Treatment")


#### ANOVA
aov_Dur_nrOTUs <- aov(Duration_num ~ Treatment, 
                      data = OTU_dynamic_cd_nrOTUs)
summary(aov_Dur_nrOTUs)
plot(aov_Dur_nrOTUs)
aov_Dur_es_nrOTUs <- effectsize::eta_squared(aov_Dur, partial = TRUE)

### Duration of emergence (rOTUs) ----
#### Descr. statistics
OTU_dynamic_cd_rOTUs <- OTU_dynamic_cd[OTU_dynamic_cd$variable %in% responsive_OTUs, ]
OTU_dynamic_cd_rOTUs_Ctrl <- OTU_dynamic_cd_rOTUs %>% filter(Treatment == "Ctrl")
OTU_dynamic_cd_rOTUs_Bti <- OTU_dynamic_cd_rOTUs %>% filter(Treatment == "Bti")


Duration_ds_rOTUs <- as.data.frame(matrix(c(mean(OTU_dynamic_cd_rOTUs_Ctrl$Duration_num), mean(OTU_dynamic_cd_rOTUs_Bti$Duration_num),
                                            sd(OTU_dynamic_cd_rOTUs_Ctrl$Duration_num), sd(OTU_dynamic_cd_rOTUs_Bti$Duration_num)), ncol = 2, byrow = TRUE)) %>% rename(Ctrl = V1,  Bti = V2) %>% mutate(Treatment = c("mean", "sd")) %>%
  column_to_rownames("Treatment")

#### ANOVA
aov_Dur_rOTUs <- aov(Duration_num ~ Treatment, 
                     data = OTU_dynamic_cd_rOTUs)
summary(aov_Dur_rOTUs)
plot(aov_Dur_es_rOTUs)
aov_Dur_es_rOTUs <- effectsize::eta_squared(aov_Dur, partial = TRUE)

### Start of emergence ----

#### Descr. statistics
OTU_dynamic_Ctrl <- OTU_dynamic %>% filter(Treatment == "Ctrl")
OTU_dynamic_Bti <- OTU_dynamic %>% filter(Treatment == "Bti")

StartEmergence_ds <- as.data.frame(matrix(c(mean(OTU_dynamic_Ctrl$Doy), mean(OTU_dynamic_Bti$Doy),
                                            sd(OTU_dynamic_Ctrl$Doy), sd(OTU_dynamic_Bti$Doy)), ncol = 2, byrow = TRUE)) %>% rename(Ctrl = V1,  Bti = V2) %>% mutate(Treatment = c("mean", "sd")) %>%
  column_to_rownames("Treatment")

#### ANOVA
aov_SE <- aov(Doy ~ Treatment, 
              data = OTU_dynamic)
summary(aov_SE)
plot(aov_SE)
aov_StartEm_es <- effectsize::eta_squared(aov_SE, partial = TRUE)

#### Ds per OTU
  ## Mean
StartEm_mean_Bti <- lapply(OTUs_means, function(x)
  mean(OTU_dynamic_Bti[OTU_dynamic_Bti$variable == x, ]$Doy)
)
StartEm_mean_Ctrl <- lapply(OTUs_means, function(x)
  mean(OTU_dynamic_Ctrl[OTU_dynamic_Ctrl$variable == x, ]$Doy)
)

names(StartEm_mean_Bti) <- OTUs_means
names(StartEm_mean_Ctrl) <- OTUs_means

StartEm_mean_Bti_df <- as.data.frame(do.call(rbind, StartEm_mean_Bti))
StartEm_mean_Ctrl_df <- as.data.frame(do.call(rbind, StartEm_mean_Ctrl))

StartEm_mean <- cbind(StartEm_mean_Ctrl_df, mean_Bti = StartEm_mean_Bti_df$V1) %>% rename(mean_Ctrl = V1)

  ## Standard deviation
StartEm_sd_Bti <- lapply(OTUs_means, function(x)
  sd(OTU_dynamic_Bti[OTU_dynamic_Bti$variable == x, ]$Doy)
)
StartEm_sd_Ctrl <- lapply(OTUs_means, function(x)
  sd(OTU_dynamic_Ctrl[OTU_dynamic_Ctrl$variable == x, ]$Doy)
)

names(StartEm_sd_Bti) <- OTUs_means
names(StartEm_sd_Ctrl) <- OTUs_means

StartEm_sd_Bti_df <- as.data.frame(do.call(rbind, StartEm_sd_Bti))
StartEm_sd_Ctrl_df <- as.data.frame(do.call(rbind, StartEm_sd_Ctrl))

StartEm_sd <- cbind(StartEm_sd_Ctrl_df, sd_Bti = StartEm_sd_Bti_df$V1) %>% 
  rename(sd_Ctrl = V1)

StartEm_ds <- cbind(StartEm_mean, StartEm_sd) %>% relocate(sd_Ctrl, .after = mean_Ctrl)

#### ANOVA per OTU
OTUID_SE <- StartEm_ds %>% 
  na.omit() %>%
  rownames_to_column("OTUID_SE") %>%
  pull(OTUID_SE)

StartEmergence_aovs <- lapply(OTUID_SE, function(x) 
  aov(Doy ~ Treatment, 
      data = OTU_dynamic[
        OTU_dynamic$variable == x, ])
)

names(StartEmergence_aovs) <- OTUID_SE

StartEmergence_anovas <- lapply(OTUID_SE, function(x) 
  summary(aov(Doy ~ Treatment, 
              data = OTU_dynamic[
                OTU_dynamic$variable == x, ]))
) 

names(StartEmergence_anovas) <- OTUID_SE

plot_SE_aovs <- lapply(StartEmergence_aovs, function(x)
  plot(x)
)

StartEm_aov_effectsize <- lapply(OTUID_SE, function(x)
  effectsize::eta_squared(StartEmergence_aovs[[x]])
)
names(StartEm_aov_effectsize) <- OTUID_SE

### Start of emergence (nrOTUs) ----
#### Descr. statistics
OTU_dynamic_nrOTUs <- OTU_dynamic[!OTU_dynamic$variable %in% responsive_OTUs, ]
OTU_dynamic_nrOTUs_Ctrl <- OTU_dynamic_nrOTUs %>% filter(Treatment == "Ctrl")
OTU_dynamic_nrOTUs_Bti <- OTU_dynamic_nrOTUs %>% filter(Treatment == "Bti")

StartEmergence_ds_nrOTUs <- as.data.frame(matrix(c(mean(OTU_dynamic_nrOTUs_Ctrl$Doy), mean(OTU_dynamic_nrOTUs_Bti$Doy),
                                                   sd(OTU_dynamic_nrOTUs_Ctrl$Doy), sd(OTU_dynamic_nrOTUs_Bti$Doy)), ncol = 2, byrow = TRUE)) %>% rename(Ctrl = V1,  Bti = V2) %>% mutate(Treatment = c("mean", "sd")) %>%
  column_to_rownames("Treatment")

#### ANOVA
aov_SE_nrOTUs <- aov(Doy ~ Treatment, 
                     data = OTU_dynamic_nrOTUs)
summary(aov_SE_nrOTUs)
plot(aov_SE_nrOTUs)
aov_StartEm_es_nrOTUs <- effectsize::eta_squared(aov_SE_nrOTUs, partial = TRUE)

### Start of emergence (rOTUs) ----
#### Descr. statistics
OTU_dynamic_rOTUs <- OTU_dynamic[OTU_dynamic$variable %in% responsive_OTUs, ]
OTU_dynamic_rOTUs_Ctrl <- OTU_dynamic_rOTUs %>% filter(Treatment == "Ctrl")
OTU_dynamic_rOTUs_Bti <- OTU_dynamic_rOTUs %>% filter(Treatment == "Bti")

StartEmergence_ds_rOTUs <- as.data.frame(matrix(c(mean(OTU_dynamic_rOTUs_Ctrl$Doy), mean(OTU_dynamic_rOTUs_Bti$Doy),
                                                  sd(OTU_dynamic_rOTUs_Ctrl$Doy), sd(OTU_dynamic_rOTUs_Bti$Doy)), ncol = 2, byrow = TRUE)) %>% rename(Ctrl = V1,  Bti = V2) %>% mutate(Treatment = c("mean", "sd")) %>%
  column_to_rownames("Treatment")

#### ANOVA
aov_SE_rOTUs <- aov(Doy ~ Treatment, 
                    data = OTU_dynamic_rOTUs)
summary(aov_SE_rOTUs)
plot(aov_SE_rOTUs)
aov_StartEm_es_rOTUs <- effectsize::eta_squared(aov_SE_rOTUs, partial = TRUE)


## Synchrony ----

  ## Data preparation
means_synchr_Ponds_list <- lapply(split_frames, function(x)
  
  x %>%
    rownames_to_column("Sample_ID") 
)

means_synchr_Ponds_list <- lapply(means_synchr_Ponds_list, function(x)
  merge(env2[, c("Sample_ID", "Doy")], x, by.y = "Sample_ID")
)

means_synchr_Ponds_list <- lapply(means_synchr_Ponds_list, function(x)
  x %>%
    dplyr::select(-Sample_ID) %>%
    mutate(Doy = as.character(Doy)) %>%
    mutate(Doy = as.numeric(Doy)) %>%
    arrange(Doy) %>%
    column_to_rownames("Doy")
)


  ## synchrony in each pond
Pond_synchrony_list_shuffle <- lapply(means_synchr_Ponds_list, function(x)
  
  community.sync(x, type = 1, nrands = 999, quiet = TRUE) # shuffle each column
)

Pond_synchrony_values_shuffle <- lapply(Pond, function(x)
  Pond_synchrony_list_shuffle[[x]][[1]]
)

Pond_synchrony_p_shuffle <- lapply(Pond, function(x)
  Pond_synchrony_list_shuffle[[x]][[4]]
)

Pond_synchr_df1 <- as.data.frame(do.call(rbind, Pond_synchrony_values_shuffle))

Pond_synchr <- Pond_synchr_df1 %>% 
  mutate(Pond = as.factor(Pond)) %>%
  dplyr::rename(Synchrony = V1) %>%
  mutate(Treatment = ifelse(Pond == "Pond1" | Pond == "Pond3" | Pond == "Pond5" | Pond == "Pond7" |  Pond == "Pond9" | Pond == "Pond11", "Ctrl", "Bti")) %>%
  mutate(Treatment = as.factor(Treatment)) 

Pond_synchr_Ctrl <- Pond_synchr[Pond_synchr$Treatment == "Ctrl", ]
Pond_synchr_Bti <- Pond_synchr[Pond_synchr$Treatment == "Bti", ]

### Descr. statistics
Sync_ds <- as.data.frame(matrix(c(mean(Pond_synchr_Ctrl$Synchrony), mean(Pond_synchr_Bti$Synchrony),
                                  sd(Pond_synchr_Ctrl$Synchrony), sd(Pond_synchr_Bti$Synchrony)), ncol = 2, byrow = TRUE)) %>% rename(Ctrl = V1,  Bti = V2) %>% mutate(Treatment = c("mean", "sd")) %>%
  column_to_rownames("Treatment")

### ANOVA
Synchrony_aov <- aov(Synchrony ~ Treatment, data = Pond_synchr)
summary(Synchrony_aov)
plot(Synchrony_aov)
aov_Sync_es <- effectsize::eta_squared(Synchrony_aov, partial = TRUE)

### Correlation (synchrony and asynchrony) ----
Pond_corr <- lapply(means_synchr_Ponds_list, function(x)
  
  as.data.frame(cor(x, method = "pearson"))
)

Pond_corr2 <- lapply(Pond_corr, function(x)
  
  x %>%
    rownames_to_column("Corr1")
)

Pond_corr_all <- bind_rows(Pond_corr2, .id = "Pond")
Pond_corr_all <- Pond_corr_all %>%
  mutate(Treatment = ifelse(Pond == "Pond1" | Pond == "Pond3" | Pond == "Pond5" | Pond == "Pond7" |Pond == "Pond9" | Pond == "Pond11", "Ctrl", "Bti")) %>%
  mutate(Treatment = as.factor(Treatment)) 

Pond_corr_all_long <- melt(Pond_corr_all, id.vars = c('Pond', 'Corr1', 'Treatment'))

Pond_corr_all_long <- Pond_corr_all_long %>%
  dplyr::rename(Corr2 = variable) %>%
  filter(!(Corr1 == "OTU_77" | Corr1 == "OTU_106" | Corr2 == "OTU_77" | Corr2 == "OTU_106"))

Pond_corr_all_long_combined <- Pond_corr_all_long %>%
  unite(Corr, Corr1, Corr2, sep = ";") %>%
  mutate(Pond = as.factor(Pond)) %>%
  mutate(value = ifelse(is.na(value), 0, value))

Pond_corr_all_Ctrl <- Pond_corr_all_long_combined[Pond_corr_all_long_combined$Treatment == "Ctrl", ]
Pond_corr_all_Bti <- Pond_corr_all_long_combined[Pond_corr_all_long_combined$Treatment == "Bti", ]

### Ds per OTU pair
  ## Mean
PC_mean_Ctrl <- lapply(unique(Pond_corr_all_Ctrl$Corr), function(x)
  mean(Pond_corr_all_Ctrl[Pond_corr_all_Ctrl$Corr == x, ]$value)
)
PC_mean_Bti <- lapply(unique(Pond_corr_all_Bti$Corr), function(x)
  mean(Pond_corr_all_Bti[Pond_corr_all_Bti$Corr == x, ]$value)
)

names(PC_mean_Ctrl) <- unique(Pond_corr_all_Ctrl$Corr)
names(PC_mean_Bti) <- unique(Pond_corr_all_Bti$Corr)

PC_mean_Ctrl_df <- as.data.frame(do.call(rbind, PC_mean_Ctrl))
PC_mean_Bti_df <- as.data.frame(do.call(rbind, PC_mean_Bti))

PC_mean <- cbind(PC_mean_Ctrl_df, mean_Bti = PC_mean_Bti_df$V1) %>% rename(mean_Ctrl = V1)

  ## Standard deviation
PC_sd_Ctrl <- lapply(unique(Pond_corr_all_Ctrl$Corr), function(x)
  sd(Pond_corr_all_Ctrl[Pond_corr_all_Bti$Corr == x, ]$value)
)
PC_sd_Bti <- lapply(unique(Pond_corr_all_Bti$Corr), function(x)
  sd(Pond_corr_all_Bti[Pond_corr_all_Bti$Corr == x, ]$value)
)

names(PC_sd_Ctrl) <- unique(Pond_corr_all_Ctrl$Corr)
names(PC_sd_Bti) <- unique(Pond_corr_all_Bti$Corr)

PC_sd_Bti_df <- as.data.frame(do.call(rbind, PC_sd_Bti))
PC_sd_Ctrl_df <- as.data.frame(do.call(rbind, PC_sd_Ctrl))

PC_sd <- cbind(PC_sd_Ctrl_df, sd_Bti = PC_sd_Bti_df$V1) %>% 
  rename(sd_Ctrl = V1)

PC_ds_raw <- cbind(PC_mean, PC_sd) %>% relocate(sd_Ctrl, .after = mean_Ctrl)
PC_ds <- PC_ds_raw %>%
  filter_all(any_vars(. != 0))

### ANOVA per OTU pair 
list_Corr <- split(Pond_corr_all_long_combined, Pond_corr_all_long_combined$Corr) 

aovs_Corr <- lapply(list_Corr, function(x)
  aov(value ~ Treatment, data = x)
)

anovas_Corr <- lapply(list_Corr, function(x)
  summary(aov(value ~ Treatment, data = x))
)

anovas_Corr_data <- lapply(names(anovas_Corr), function(x)
  anovas_Corr[[x]][[1]][c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")]
)

names(anovas_Corr_data) <- names(anovas_Corr)

Corr_data <- as.data.frame(do.call(rbind, anovas_Corr_data)) %>%
  rownames_to_column("OTU_pairs") %>%
  filter(!grepl("Residuals", OTU_pairs)) %>%
  mutate(OTU_pairs = gsub(".Treatment", "", OTU_pairs)) %>%
  na.omit()

OTU_pairs_noO <- c(rownames(PC_ds))

Corr_aov_effectsize <- lapply(OTU_pairs_noO, function(x)
  effectsize::eta_squared(aovs_Corr[[x]])
)
names(Corr_aov_effectsize) <- OTU_pairs_noO


Corr_sign <- Corr_data[Corr_data$`Pr(>F)` < 0.05, c("OTU_pairs")] # significant diff. in correlation

PC_ds_sign <- PC_ds %>% 
  rownames_to_column("OTU_pairs") %>%
  filter(OTU_pairs %in% Corr_sign) %>%
  dplyr::select(OTU_pairs, Ctrl = mean_Ctrl, Bti = mean_Bti)

PC_sign_corr_diff <- PC_ds_sign %>%
  mutate(Diff = case_when(  
    Ctrl >= 0 & Bti >= 0 & Ctrl > Bti ~ "Ctrl > Bti",
    Ctrl >= 0 & Bti >= 0 & Ctrl < Bti ~ "Bti > Ctrl",
    Ctrl <= 0 & Bti <= 0 & Ctrl > Bti ~ "Ctrl < Bti",
    Ctrl <= 0 & Bti <= 0 & Ctrl < Bti ~ "Bti < Ctrl",
    Ctrl >= 0 & Bti < 0  ~ "-Bti +Ctrl",
    Ctrl <= 0 & Bti > 0  ~ "+Bti -Ctrl"
  )) # categories

PC_sign_corr_diff %>%
  group_by(Diff) %>%
  summarise(count = n()) # count per category


## Check if there are correlations of OTUs identified as the same species
Corr_sign_comp <- Corr_sign %>%
  as.data.frame() %>%
  rename(Corr = 1) %>%
  mutate(Correl = Corr) %>%
  separate(Corr, into = c("Corr1", "Corr2"), sep = ";")
  
Corr1_df <- Corr_sign_comp[, c("Correl", "Corr1")] %>% rename(OTU_ID = Corr1)
Corr2_df <- Corr_sign_comp[, c("Correl",  "Corr2")]%>% rename(OTU_ID = Corr2)

Result_Ident <- OTUs[, c("OTU_ID", "Result")]

Corr1_merge <- merge(Corr1_df, Result_Ident, by = "OTU_ID")
Corr2_merge <- merge(Corr2_df, Result_Ident, by = "OTU_ID")

Corr_merge <- merge(Corr1_merge, Corr2_merge, by = "Correl")
Corr_merge_TF <- Corr_merge %>%
  mutate(Identical = if_else(OTU_ID.x == OTU_ID.y, TRUE, FALSE))

  ## Are there sign. different correlations of two OTUs of the same species?
Corr_merge_TF %>%
  filter(Identical == "TRUE") # No!
