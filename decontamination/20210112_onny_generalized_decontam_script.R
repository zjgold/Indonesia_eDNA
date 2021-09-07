# ---
# title: "Generalized Decontamination Script"
# author: "Zack Gold"
# date: "11/02/2020"
# output: html_document
# ---

# Load Libraries ----------------------------------------------------------------------------------------------------------------
library (tidyverse)
library (vegan)
library (proxy)
library(reshape2)
library(microDecon)
library(stringr)
library(knitr)
library(ranacapa)
library(dplyr)
library(tibble)
library(reshape2)
library(phyloseq)
library(ggplot2)
library(plotly)
library(optparse)
library(fitdistrplus)
library(rstan)
library(shinystan)
library(bayesplot)
library(broom)
library(analyze.stuff)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source("/Users/zackgold/Documents/UCLA_phd/Projects/anacapa/decontam/decontamination_utilities.R")
library(here)

#Test inputs ----------------------------------------------------------------------------------------------------------------
working_directory <- here("decontamination")
input_anacapa_path_1 <- here("data","CO1_ASV_raw_taxonomy_60_edited.txt")
input_anacapa_path_2 <- here("data","bio_12S_miu_ASV_raw_taxonomy_60_edited.txt")
input_anacapa_path_3 <- here("data","bio_12S_elas_ASV_raw_taxonomy_60_edited.txt")

input_meta_path <- "/Users/zackgold/Documents/UCLA_phd/Projects/Indonesia/Indo_moorea_oct2020/ONNY/onny_metadata_11022020.txt"
number_anacapa_tables <- "three"

read_type <- "merged_only"

low_dist_probability_cutoff <- 0.05
minimum_read_cutoff <- 5000

step3_runs <- 2
step3_thresh <- 0.7
step3_prop.thresh <- 5e-05
step3_regression <- 0

som_level <- "single"
som_chains<-10
som_iterations<-1000
som_filter_choice <- "som_min_threshold_1"

max_dist_probability_cutoff <- 0.9999
#---

#---

#Load Data ----------------------------------------------------------------------------------------------------------------
setwd(working_directory)
dir.create("Output_csv")
dir.create("Output_R")
dir.create("Output_plots")

metadata <- read.table(input_meta_path, header = TRUE, sep = "\t", stringsAsFactors = F)

anacapa_table_1 <- read.table(input_anacapa_path_1, header = 1, sep = "\t", stringsAsFactors = F)
anacapa_table_2 <- read.table(input_anacapa_path_2, header = 1, sep = "\t", stringsAsFactors = F)
anacapa_table_3 <- read.table(input_anacapa_path_3, header = 1, sep = "\t", stringsAsFactors = F)

#Generate Hash.key
anacapa_table_1 %>% 
  dplyr::select(seq_number, sum.taxonomy) -> hash.key_1

anacapa_table_2 %>% 
  dplyr::select(seq_number, sum.taxonomy) -> hash.key_2

anacapa_table_3 %>% 
  dplyr::select(seq_number, sum.taxonomy) -> hash.key_3

rbind(hash.key_1,hash.key_2,hash.key_3) -> hash.key

#Format Anacapa Table 1
anacapa_table_1 %>%
  filter(., str_detect(seq_number,"merged")) %>% 
  dplyr::select(seq_number) -> merged_hash_1
merged_hash_1$seq_number[[1]] %>% str_remove(.,"merged_") %>% str_sub(., end=-3) -> barcode_1
anacapa_table_1$Miseq_run <- barcode_1

###Fix Names
anacapa_table_1_names <- colnames(anacapa_table_1)
anacapa_table_1_names %>% as.data.frame() -> anacapa_table_1_names
colnames(anacapa_table_1_names) <- c("Seq_number")
left_join(anacapa_table_1_names, metadata) %>% dplyr::select(Seq_number, New_name) %>% distinct() %>% 
  mutate(final_names = coalesce(New_name, Seq_number)) %>%  dplyr::select(final_names) %>% as.list() -> anacapa_table_1_names
colnames(anacapa_table_1) <- anacapa_table_1_names$final_names

#Format Anacapa Table 2
anacapa_table_2 %>%
  filter(., str_detect(seq_number,"merged")) %>% 
  dplyr::select(seq_number) -> merged_hash_2
merged_hash_2$seq_number[[1]] %>% str_remove(.,"merged_") %>% str_sub(., end=-3) -> barcode_2
anacapa_table_2$Miseq_run <- barcode_2

####Fix Names
anacapa_table_2_names <- colnames(anacapa_table_2)
anacapa_table_2_names %>% as.data.frame() -> anacapa_table_2_names
colnames(anacapa_table_2_names) <- c("Seq_number")
left_join(anacapa_table_2_names, metadata) %>% dplyr::select(Seq_number, New_name) %>% 
  mutate(final_names = coalesce(New_name, Seq_number)) %>%  dplyr::select(final_names) %>% as.list() -> anacapa_table_2_names
colnames(anacapa_table_2) <- anacapa_table_2_names$final_names

#Format Anacapa Table 3
anacapa_table_3 %>%
  filter(., str_detect(seq_number,"merged")) %>% 
  dplyr::select(seq_number) -> merged_hash_3
merged_hash_3$seq_number[[1]] %>% str_remove(.,"merged_") %>% str_sub(., end=-3) -> barcode_3
anacapa_table_3$Miseq_run <- barcode_3

####Fix Names
anacapa_table_3_names <- colnames(anacapa_table_3)
anacapa_table_3_names %>% as.data.frame() -> anacapa_table_3_names
colnames(anacapa_table_3_names) <- c("Seq_number")
left_join(anacapa_table_3_names, metadata) %>% dplyr::select(Seq_number, New_name) %>% 
  mutate(final_names = coalesce(New_name, Seq_number)) %>%  dplyr::select(final_names) %>% as.list() -> anacapa_table_3_names
colnames(anacapa_table_3) <- anacapa_table_3_names$final_names

#Merge All Tables
ASV.table <- bind_rows(anacapa_table_1,anacapa_table_2,anacapa_table_3)

ASV.table %>% 
  replace(is.na(.), 0) ->ASV.table

tail(ASV.table)
head(ASV.table)


# Convert to Long Data

###Format for Long Data
ASV.table$seq_number <- factor(ASV.table$seq_number)
ASV.table$Miseq_run <- factor(ASV.table$Miseq_run)

columns <- colnames(ASV.table)
remove <- c("seq_number","sum.taxonomy","Miseq_run")

gathercols <-  columns[! columns %in% remove] 

ASV.table %>% 
  pivot_longer(., cols=gathercols, names_to="sample", values_to="reads") -> ASV.table_long

ASV.table_long$reads <- as.numeric(ASV.table_long$reads)

saveRDS(ASV.table_long,  file="ASV.table_long")

#---

#---

# Cleaning Process 0: Remove all Forward, Reverse, and Unmerged reads & Remove Singletons ----------------------------------------------------------------------------------------------------------------

#Filter Merged only reads

if(read_type =="merged_only"){
  #Filter out forward, reverse, and unmerged reads
  ASV.table_long %>% 
    filter(., str_detect(seq_number,"merged")) %>% 
    filter(., !str_detect(seq_number,"unmerged")) -> ASV.table_used
} else if ( read_type =="merged_and_forward"){
  #Remove Singletons
  ASV.table_long %>% 
    filter(., !str_detect(seq_number,"reverse")) %>% 
    filter(., !str_detect(seq_number,"unmerged")) -> ASV.table_used
} else {
    #Remove Singletons
    ASV.table_long %>%
      dplyr::group_by(seq_number) %>%
      mutate (TotalReadsperSample = sum(reads)) %>% 
      filter(., TotalReadsperSample > 1) %>% 
      dplyr::select(-TotalReadsperSample) -> ASV.table_used
}

#Calculate % ASVs Kept
ASV.table_long %>%  dim() -> all_dim
ASV.table_used %>%  dim() -> used_only_dim

paste0(round(used_only_dim[[1]]/all_dim[[1]]*100,2),"% ASVs retained")

ASV.table_used %>% 
  ungroup() %>% 
  filter(., sum.taxonomy !="" ) %>% 
  dplyr::group_by(Miseq_run) %>% 
  filter(., !str_detect(sample,"SB_")) %>% 
  dplyr::summarise(ASVs = n_distinct(seq_number), Tot_reads = sum(reads), n_distinct(sample))
#---

#---

# Cleaning Process 1: Estimation of *Tag-jumping* or sample *Cross-talk* ----------------------------------------------------------------------------------------------------------------

## Step 1: Nest the dataset by origin of ASVs

###Identify Positives, Negatives, and Samples

###Create list of control samples
metadata %>% 
  filter(Sample_Control=="Control") %>% 
  dplyr::select(New_name) %>% unique() -> controls
controls <- controls$New_name

metadata %>% 
  filter(Control_Type=="Pos") %>% 
  dplyr::select(New_name) %>% unique()-> pos_controls
pos_controls <- pos_controls$New_name

metadata %>% 
  filter(Control_Type=="Blank") %>% 
  dplyr::select(New_name) %>% unique() -> neg_controls
neg_controls <- neg_controls$New_name

###New column that labels each ASV as from Positive (control) or Sample
ASV.table_used %>% 
  mutate(source = case_when(sample %in% pos_controls~"Positives",
                            sample %in% neg_controls~"Blanks",
                            TRUE ~"Samples")) -> ASV.table_used


###Convert to tibble
ASV.table_used <- as_tibble(ASV.table_used)

###Remove empty sequences
ASV.table_used %>% 
  filter(reads != 0)  -> ASV.table_used

###Rename Columns and remove seq_number
ASV.table_used %>%
  mutate(sample = as.character(sample),
         nReads = reads)  -> ASV.table_used

###ASVs in Positive Controls
ASV.table_used %>% 
  filter (source == "Positives") %>%
  dplyr::group_by(seq_number) %>% 
  dplyr::summarise(tot = sum(reads)) %>% 
  arrange(desc(tot)) %>% 
  pull(seq_number) -> all.seqs.in.positives

hash.key %>% 
  filter(seq_number %in% all.seqs.in.positives) %>% as_tibble() -> pos.contam.species

write.csv(pos.contam.species, file="pos.contam.species.csv")

###ASVs in Negative Controls

ASV.table_used %>% 
  filter (source == "Blanks") %>%
  dplyr::group_by(seq_number) %>% 
  dplyr::summarise(tot = sum(reads)) %>% 
  arrange(desc(tot)) %>% 
  pull(seq_number) -> all.seqs.in.blanks

hash.key %>% 
  filter(seq_number %in% all.seqs.in.blanks) %>% as_tibble() -> blank.contam.species

write.csv(blank.contam.species, file="blank.contam.species.csv")

### Visualize Read Counts Across Samples for Barcode_1
ASV.table_used %>% 
  group_by(sample) %>%
  filter(., Miseq_run== barcode_1) %>% 
  mutate (TotalReadsperSample = sum(nReads)) %>%
  arrange(desc(TotalReadsperSample)) %>%
  ggplot(., aes(x=sample, y=TotalReadsperSample, color=source)) + 
  geom_point() +ggtitle("Read Count Across Samples") + 
  theme(axis.text.x = element_text(angle = 90)) -> plot_1

ggsave(plot=plot_1,"Output_plots/Sample_Read_Depth_Barcode_1.png", device = "png", width = 12, height = 8, units = "in")

### Visualize Read Counts Across Samples for Barcode_2

ASV.table_used %>% 
  group_by(sample) %>%
  filter(., Miseq_run== barcode_2) %>% 
  mutate (TotalReadsperSample = sum(nReads)) %>%
  arrange(desc(TotalReadsperSample)) %>%
  ggplot(., aes(x=sample, y=TotalReadsperSample, color=source)) + 
  geom_point() +ggtitle("Read Count Across Samples") + 
  theme(axis.text.x = element_text(angle = 90)) -> plot_2

ggsave(plot=plot_2,"Output_plots/Sample_Read_Depth_Barcode_2.png", device = "png", width = 12, height = 8, units = "in")

ASV.table_used %>% 
  group_by(sample) %>%
  filter(., Miseq_run== barcode_3) %>% 
  mutate (TotalReadsperSample = sum(nReads)) %>%
  arrange(desc(TotalReadsperSample)) %>%
  ggplot(., aes(x=sample, y=TotalReadsperSample, color=source)) + 
  geom_point() + ggtitle("Read Count Across Samples") + 
  theme(axis.text.x = element_text(angle = 90)) -> plot_3

ggsave(plot=plot_3,"Output_plots/Sample_Read_Depth_Barcode_3.png", device = "png", width = 12, height = 8, units = "in")


###Nesting the dataset
ASV.table_used %>% 
  dplyr::group_by(Miseq_run, source) %>% 
  nest() %>% 
  pivot_wider(names_from=source, values_from=data) -> ASV.nested


####Summary.file.1
ASV.nested %>% 
  ungroup() %>% 
  dplyr::transmute(.,Miseq_run,Summary = purrr::map(Samples, ~ how.many(ASVtable = ., round = 0)))  -> ASV.summary

ASV.summary$Summary

## Step 2: Model the composition of the positive controls of each run 

###Jumping vector

ASV.nested %>% 
  ungroup() %>% 
  mutate(contam.tibble = purrr::map(Positives, 
                                     function(.x){
                                       .x %>%
                                         ungroup() %>% 
                                         group_by(sample) %>%
                                         mutate (TotalReadsperSample = sum(nReads)) %>%
                                         mutate (proportion = nReads/TotalReadsperSample) %>%
                                         group_by(seq_number) %>%
                                         dplyr::summarise (vector_contamination = max(proportion))
                                     }) ) -> ASV.nested

###Vector Contamination in Barcode_1
ASV.nested$contam.tibble[[1]] %>% as.data.frame() %>% 
  ggplot(aes(x= vector_contamination))+
  geom_histogram() -> vc_plot_1 # Check how it looks
ggsave(plot=vc_plot_1,"Output_plots/Vector_Contamination_Barcode_1.png", device = "png", width = 12, height = 8, units = "in")

###Vector Contamination in Barcode_2

ASV.nested$contam.tibble[[2]] %>% as.data.frame() %>% 
  ggplot(aes(x= vector_contamination))+
  geom_histogram() -> vc_plot_2 # Check how it looks
ggsave(plot=vc_plot_2,"Output_plots/Vector_Contamination_Barcode_2.png", device = "png", width = 12, height = 8, units = "in")

ASV.nested$contam.tibble[[3]] %>% as.data.frame() %>% 
  ggplot(aes(x= vector_contamination))+
  geom_histogram() -> vc_plot_3 # Check how it looks
ggsave(plot=vc_plot_3,"Output_plots/Vector_Contamination_Barcode_3.png", device = "png", width = 12, height = 8, units = "in")


##Step 3: Substract the composition of the positive controls from the environment samples

ASV.nested %>% 
  ungroup() %>% 
  mutate(Step1.cleaned.tibble = map2(Samples, contam.tibble, function(.x,.y){ 
    .x %>%
      dplyr::group_by (sample) %>%
      mutate (TotalReadsperSample = sum (nReads)) %>%
      left_join(.y, by = "seq_number") %>%
      mutate (Updated_nReads = ifelse (!is.na(vector_contamination),  nReads - (ceiling(vector_contamination*TotalReadsperSample)), nReads)) %>%
      filter (Updated_nReads > 0) %>%
      ungroup() %>% 
      dplyr::select (sample, seq_number, nReads = Updated_nReads)
  })) -> ASV.nested

###Add this step to the summary table we were creating

####Summary.file.2
ASV.nested %>% 
  transmute(Miseq_run, Summary.1 = purrr::map(Step1.cleaned.tibble, ~ how.many(ASVtable = .,round = "1.Jump"))) %>% 
  left_join(ASV.summary) %>% #use left join when there are many miseq runs to join
  mutate(Summary   = map2(Summary, Summary.1, bind_rows)) %>%
  dplyr::select(-Summary.1) -> ASV.summary 

ASV.summary$Summary

#---

#---

# Cleaning Process 2: microDecon - Clearance of Negative Control Contaminants ----------------------------------------------------------------------------------------------------------------

## Anacapa_barcode_1

### Format for microDecon

####Arrange Blank Samples to be the first 3 columns
##### Identify blank samples
if (is.null(ASV.nested$Blanks[[1]])){
  ASV.nested$Step1.cleaned.tibble[[1]] %>% 
    arrange(desc(nReads)) %>% 
    slice(1) %>% 
    mutate(., nReads =1, sample="Fake_blank_1") -> test1
  tibble(sample=c("Fake_blank_1"), seq_number=c("fake_asv_1"),nReads =c(10000)) -> blank_ASVs_1
  full_join(test1,blank_ASVs_1) -> blank_ASVs_1
  #### Number of Blanks
  blank_ASVs_1$sample %>% unique() %>% length() -> blank_num
  #### Columns to Ignore for Grouping
  blank_ASVs_1$sample %>% unique() -> blank_samples
} else{
  ASV.nested$Blanks[[1]] %>% 
    dplyr::select(sample,seq_number, nReads) -> blank_ASVs_1
  #### Number of Blanks
  ASV.nested$Blanks[[1]]$sample %>% unique() %>% length() -> blank_num
  #### Columns to Ignore for Grouping
  ASV.nested$Blanks[[1]]$sample %>% unique() -> blank_samples
}

##### Add blank samples
rbind(blank_ASVs_1, ASV.nested$Step1.cleaned.tibble[[1]] 
      %>% dplyr::select(sample,seq_number,nReads)%>% arrange(sample)) -> Step.2_barcode1

#### Pivot to Wide format
Step.2_barcode1 %>% 
  pivot_wider(names_from=sample, values_from=nReads, values_fill = list(nReads =0)) -> Step.2_barcode1_wide

Step.2_barcode1_wide %>%
  left_join(hash.key)  %>% as.data.frame -> Step.2_barcode1_wide

### Identify Groups
#### Number of Samples
ASV.nested$Step1.cleaned.tibble[[1]]$sample %>% unique() %>% length() -> sample_num

append(blank_samples, c("sum.taxonomy","seq_number")) -> col_delete

#### Assign Groups -> should be replaced with pre-defined groups in the metadata file
if(som_level =="double"){
colnames(Step.2_barcode1_wide) %>%  as.data.frame() %>% 
  filter(., !(. %in% col_delete)) %>% 
  mutate(., New_name=.) %>% 
  left_join(metadata)  %>%
  dplyr::select(New_name,Site,Bio_rep) %>% 
  distinct() %>% 
  group_by(Site,Bio_rep) %>% 
    dplyr::count() -> sample_counter
} else{
  colnames(Step.2_barcode1_wide) %>%  as.data.frame() %>% 
    filter(., !(. %in% col_delete)) %>% 
    mutate(., New_name=.) %>% 
    left_join(metadata, by="New_name")  %>%
    dplyr::select(New_name,Site) %>% 
    distinct() %>% 
    group_by(Site) %>%
    dplyr::count() -> sample_counter
}

colnames(Step.2_barcode1_wide) %>%  as.data.frame() %>% 
  filter(., !(. %in% col_delete)) %>% 
  mutate(., New_name=.) %>% 
  left_join(metadata, by="New_name")  %>%
  dplyr::select(New_name,Site) %>% 
  distinct() %>%  dplyr::select(Site) %>% distinct() -> correct_order

correct_order %>% 
  left_join(sample_counter) -> sample_counter

Step.2_barcode1_decon <- decon(data=Step.2_barcode1_wide, numb.blanks = blank_num, numb.ind=sample_counter$n, taxa = T,
                        runs = step3_runs,
                        thresh = step3_thresh,
                        prop.thresh = step3_prop.thresh,
                        regression = step3_regression)

saveRDS(Step.2_barcode1_decon, file="Step.2_barcode1_decon.RDS")

### Convert to Long Data

columns <- colnames(Step.2_barcode1_decon$decon.table)
remove <- c("seq_number","sum.taxonomy")

gathercols <-  columns[! columns %in% remove] 


Step.2_barcode1_decon$decon.table %>% 
  pivot_longer(., cols=gathercols, names_to="sample", values_to="nReads") -> Step.2_barcode1_decon_clean

### Final Clean Up
Step.2_barcode1_decon_clean %>% as_tibble() %>% 
  dplyr::select(-sum.taxonomy) %>% 
  filter(., !sample %in% c("Mean.blank", "fake_asv_1")) %>% 
  filter(., !seq_number =="fake_asv_1") %>% 
  filter(., nReads >0) -> Step.2_barcode1_decon_clean_tibble


## Anacapa_barcode_2

### Format for microDecon

####Arrange Blank Samples to be the first 3 columns
##### Identify blank samples
if (is.null(ASV.nested$Blanks[[2]])){
  ASV.nested$Step1.cleaned.tibble[[2]] %>% 
    arrange(desc(nReads)) %>% 
    slice(1) %>% 
    mutate(., nReads =1, sample="Fake_blank_1") -> test1
  tibble(sample=c("Fake_blank_1"), seq_number=c("fake_asv_1"),nReads =c(10000)) -> blank_ASVs_2
  full_join(test1,blank_ASVs_2) -> blank_ASVs_2
  #### Number of Blanks
  blank_ASVs_2$sample %>% unique() %>% length() -> blank_num
  #### Columns to Ignore for Grouping
  blank_ASVs_2$sample %>% unique() -> blank_samples
} else{
ASV.nested$Blanks[[2]] %>% 
  dplyr::select(sample,seq_number, nReads) -> blank_ASVs_2
  #### Number of Blanks
  ASV.nested$Blanks[[2]]$sample %>% unique() %>% length() -> blank_num
  #### Columns to Ignore for Grouping
  ASV.nested$Blanks[[2]]$sample %>% unique() -> blank_samples
}
##### Add blank samples
rbind(blank_ASVs_2, ASV.nested$Step1.cleaned.tibble[[2]]%>% dplyr::select(sample,seq_number,nReads) %>% arrange(sample)) -> Step.2_barcode2

#### Pivot to Wide format
Step.2_barcode2 %>% 
  pivot_wider(names_from=sample, values_from=nReads, values_fill = list(nReads =0)) -> Step.2_barcode2_wide

Step.2_barcode2_wide %>%
  left_join(hash.key)  %>% as.data.frame -> Step.2_barcode2_wide

### Identify Groups

#### Number of Samples
ASV.nested$Step1.cleaned.tibble[[2]]$sample %>% unique() %>% length() -> sample_num

append(blank_samples, c("sum.taxonomy","seq_number")) -> col_delete

#### Assign Groups -> should be replaced with pre-defined groups in the metadata file
if(som_level =="double"){
  colnames(Step.2_barcode2_wide) %>%  as.data.frame() %>% 
    filter(., !(. %in% col_delete)) %>% 
    mutate(., New_name=.) %>% 
    left_join(metadata)  %>%
    dplyr::select(New_name,Site,Bio_rep) %>% 
    distinct() %>% 
    group_by(Site,Bio_rep) %>% 
    dplyr::count() -> sample_counter
} else{
  colnames(Step.2_barcode2_wide) %>%  as.data.frame() %>% 
    filter(., !(. %in% col_delete)) %>% 
    mutate(., New_name=.) %>% 
    left_join(metadata)  %>%
    dplyr::select(New_name,Site) %>% 
    distinct() %>% 
    group_by(Site) %>% 
    dplyr::count() -> sample_counter
}

colnames(Step.2_barcode2_wide) %>%  as.data.frame() %>% 
  filter(., !(. %in% col_delete)) %>% 
  mutate(., New_name=.) %>% 
  left_join(metadata, by="New_name")  %>%
  dplyr::select(New_name,Site) %>% 
  distinct() %>%  dplyr::select(Site) %>% distinct() -> correct_order

correct_order %>% 
  left_join(sample_counter) -> sample_counter

Step.2_barcode2_decon <- decon(data=Step.2_barcode2_wide, numb.blanks = blank_num, numb.ind=sample_counter$n, taxa = T,
                               runs = step3_runs,
                               thresh = step3_thresh,
                               prop.thresh = step3_prop.thresh,
                               regression = step3_regression)

saveRDS(Step.2_barcode2_decon, file="Step.2_barcode2_decon.RDS")

### Convert to Long Data

columns <- colnames(Step.2_barcode2_decon$decon.table)
remove <- c("seq_number","sum.taxonomy")

gathercols <-  columns[! columns %in% remove] 

Step.2_barcode2_decon$decon.table %>% 
  pivot_longer(., cols=gathercols, names_to="sample", values_to="nReads") -> Step.2_barcode2_decon_clean

### Final Clean Up
Step.2_barcode2_decon_clean %>% as_tibble() %>% 
  dplyr::select(-sum.taxonomy) %>% 
  filter(., !sample %in% c("Mean.blank", "fake_asv_1")) %>%
  filter(., !seq_number =="fake_asv_1") %>% 
  filter(., nReads >0) -> Step.2_barcode2_decon_clean_tibble


## Three
if (is.null(ASV.nested$Blanks[[3]])){
  ASV.nested$Step1.cleaned.tibble[[3]] %>% 
    arrange(desc(nReads)) %>% 
    slice(1) %>% 
    mutate(., nReads =1, sample="Fake_blank_1") -> test1
  tibble(sample=c("Fake_blank_1"), seq_number=c("fake_asv_1"),nReads =c(10000)) -> blank_ASVs_3
  full_join(test1,blank_ASVs_3) -> blank_ASVs_3
  #### Number of Blanks
  blank_ASVs_3$sample %>% unique() %>% length() -> blank_num
  #### Columns to Ignore for Grouping
  blank_ASVs_3$sample %>% unique() -> blank_samples
} else{
  ASV.nested$Blanks[[3]] %>% 
    dplyr::select(sample,seq_number, nReads) -> blank_ASVs_3
  #### Number of Blanks
  ASV.nested$Blanks[[3]]$sample %>% unique() %>% length() -> blank_num
  #### Columns to Ignore for Grouping
  ASV.nested$Blanks[[3]]$sample %>% unique() -> blank_samples
}
##### Add blank samples
rbind(blank_ASVs_3, ASV.nested$Step1.cleaned.tibble[[3]]%>% dplyr::select(sample,seq_number,nReads) %>%  arrange(sample)) -> Step.2_barcode3

#### Pivot to Wide format
Step.2_barcode3 %>% 
  pivot_wider(names_from=sample, values_from=nReads, values_fill = list(nReads =0)) -> Step.2_barcode3_wide

Step.2_barcode3_wide %>%
  left_join(hash.key)  %>% as.data.frame -> Step.2_barcode3_wide

### Identify Groups

#### Number of Samples
ASV.nested$Step1.cleaned.tibble[[3]]$sample %>% unique() %>% length() -> sample_num

append(blank_samples, c("sum.taxonomy","seq_number")) -> col_delete

#### Assign Groups -> should be replaced with pre-defined groups in the metadata file
if(som_level =="double"){
  colnames(Step.2_barcode3_wide) %>%  as.data.frame() %>% 
    filter(., !(. %in% col_delete)) %>% 
    mutate(., New_name=.) %>% 
    left_join(metadata)  %>%
    dplyr::select(New_name,Site,Bio_rep) %>% 
    distinct() %>% 
    group_by(Site,Bio_rep) %>% 
    dplyr::count() -> sample_counter
} else{
  colnames(Step.2_barcode3_wide) %>%  as.data.frame() %>% 
    filter(., !(. %in% col_delete)) %>% 
    mutate(., New_name=.) %>% 
    left_join(metadata)  %>%
    dplyr::select(New_name,Site) %>% 
    distinct() %>% 
    group_by(Site) %>% 
    dplyr::count() -> sample_counter
}

colnames(Step.2_barcode3_wide) %>%  as.data.frame() %>% 
  filter(., !(. %in% col_delete)) %>% 
  mutate(., New_name=.) %>% 
  left_join(metadata, by="New_name")  %>%
  dplyr::select(New_name,Site) %>% 
  distinct() %>%  dplyr::select(Site) %>% distinct() -> correct_order

correct_order %>% 
  left_join(sample_counter) -> sample_counter

Step.2_barcode3_decon <- decon(data=Step.2_barcode3_wide, numb.blanks = blank_num, numb.ind=sample_counter$n, taxa = T,
                               runs = step3_runs,
                               thresh = step3_thresh,
                               prop.thresh = step3_prop.thresh,
                               regression = step3_regression)

saveRDS(Step.2_barcode3_decon, file="Step.2_barcode3_decon.RDS")

### Convert to Long Data

columns <- colnames(Step.2_barcode3_decon$decon.table)
remove <- c("seq_number","sum.taxonomy")

gathercols <-  columns[! columns %in% remove] 

Step.2_barcode3_decon$decon.table %>% 
  pivot_longer(., cols=gathercols, names_to="sample", values_to="nReads") -> Step.2_barcode3_decon_clean

### Final Clean Up
Step.2_barcode3_decon_clean %>% as_tibble() %>% 
  dplyr::select(-sum.taxonomy) %>% 
  filter(., !sample %in% c("Mean.blank", "fake_asv_1")) %>%
  filter(., !seq_number =="fake_asv_1") %>% 
  filter(., nReads >0) -> Step.2_barcode3_decon_clean_tibble

## Add to ASV.nested
ASV.nested %>% 
  mutate(., Step2.tibble = list(Step.2_barcode1_decon_clean_tibble,Step.2_barcode2_decon_clean_tibble,Step.2_barcode3_decon_clean_tibble)) -> ASV.nested

ASV.nested %>% 
  transmute(Miseq_run, Summary.1 = purrr::map(Step2.tibble, ~ how.many(ASVtable = .,round = "2.microDecon"))) %>%
  left_join(ASV.summary) %>% 
  mutate(Summary = map2(Summary, Summary.1, bind_rows)) %>%
  dplyr::select(-Summary.1) -> ASV.summary

ASV.summary$Summary


#---

#---

# Cleaning Process 3: **Discarding PCR replicates with low number of reads** ----------------------------------------------------------------------------------------------------------------

###Pull out Sample Read Depth

ASV.nested %>% 
  dplyr::select(Miseq_run, Step2.tibble) %>% 
  unnest(Step2.tibble) %>% 
  group_by(Miseq_run,sample) %>%
  dplyr::summarise(tot = sum(nReads)) %>% 
  arrange(desc(tot)) %>% 
  unite(col = sample_name, c("Miseq_run","sample"))-> all.reps


all.reps %>%  
  pull(tot) -> reads.per.sample

names(reads.per.sample) <- all.reps %>% pull(sample_name)  


### Fit Normal Distribution

fit <- fitdist(reads.per.sample, "gamma", lower=c(0,0), start=list(scale=1,shape=1))

all.reps %>%  
  mutate(prob = pgamma(tot, shape = fit$estimate[[2]], scale = fit$estimate[[1]], lower.tail = TRUE,
                       log.p = FALSE)) -> all.reps

### Remove Outlier Samples
outliers <- all.reps %>% 
  filter(prob < low_dist_probability_cutoff  | tot < minimum_read_cutoff) # changed to 0.05 to save the two samples

ASV.nested %>% 
  mutate (Step2.tibble.edited = purrr::map(Step2.tibble, 
                                                   function(.x){
                                                     .x %>%
                                                       mutate(., Miseq_run= case_when(str_detect(seq_number,barcode_1) ~ barcode_1,
                                                                                      str_detect(seq_number,barcode_2) ~ barcode_2,
                                                                                      str_detect(seq_number,barcode_3) ~ barcode_3)) %>% 
                                                       unite(col = sample_name, c("Miseq_run","sample"), remove = FALSE) %>% 
                                                       dplyr::select(-Miseq_run)
                                                   }) ) -> ASV.nested


ASV.nested %>% 
  mutate(Step.3.low.read.depth = purrr::map (Step2.tibble.edited,
                                             ~ filter(.,!sample_name %in% outliers$sample_name) %>% ungroup)) -> ASV.nested

ASV.nested %>% 
  transmute(Miseq_run, Summary.1 = purrr::map(Step.3.low.read.depth, ~ how.many(ASVtable = .,round = "3.Low.Read.Depth"))) %>% 
  left_join(ASV.summary) %>% 
  mutate(Summary   = map2(Summary, Summary.1, bind_rows)) %>%
  dplyr::select(-Summary.1) -> ASV.summary 

ASV.summary$Summary
#---

#---

# Cleaning Process 4: **Dissimilarity between PCR (biological) replicates** ----------------------------------------------------------------------------------------------------------------

ASV.nested %>% 
  dplyr::select(Miseq_run,Step.3.low.read.depth) %>% 
  unnest(Step.3.low.read.depth) %>% 
  as.data.frame() %>% 
  ungroup() %>% 
  left_join(metadata, by=c("sample"="New_name"))-> cleaned.tibble.pre_occ


## How many samples, how many ASVs

cleaned.tibble.pre_occ %>% 
  dplyr::summarise(n_distinct(sample),
                   n_distinct(seq_number))


## eDNA Index
cleaned.tibble.pre_occ %>%
  dplyr::group_by (sample) %>%
  mutate (Tot = sum(nReads),
          Row.sums = nReads / Tot) %>% 
  dplyr::group_by (seq_number) %>%
  mutate (Colmax = max (Row.sums),
          Normalized.reads = Row.sums / Colmax) -> cleaned.tibble.pre_occ #transforms raw number of reads to eDNA index

#Calculate Vegan Distances
tibble_to_vegdist_all (cleaned.tibble.pre_occ) -> all.distances.full.pre


as.tibble(subset(melt(as.matrix(all.distances.full.pre)))) -> all.distances.melted.pre


  # Now, create a three variables for all distances, they could be PCR replicates, BIOL replicates, or from the same site
  all.distances.melted.pre %>%
    separate(Var1, into = c("Miseq_run", "New_name"), sep=":", remove=F) %>% 
    left_join(metadata) %>% 
    mutate(., Site1=Site,
           Bio_rep1=Bio_rep,
           Miseq_run1=Miseq_run) %>% 
    dplyr::select(Var1,Miseq_run1,Site1,Bio_rep1,Var2,value) %>% 
    separate(Var2, into = c("Miseq_run", "New_name"), sep=":", remove=F) %>% 
    left_join(metadata) %>% 
    mutate(., Site2=Site,
           Bio_rep2=Bio_rep,
           Miseq_run2=Miseq_run) %>% 
    dplyr::select(Var1,Var2,Site1,Miseq_run1,Bio_rep1,Miseq_run2,Site2,Bio_rep2,value) %>%
    unite( Site1, Bio_rep1, col= "station1", remove=F) %>% 
    unite( Site2, Bio_rep2, col= "station2", remove=F) %>% 
    unite( Miseq_run1, station1, col= "typers1", remove=F) %>% 
    unite( Miseq_run2, station2, col= "typers2", remove=F) %>% 
    unite( Miseq_run1, Site1, col= "type_site1", remove=F) %>% 
    unite( Miseq_run2, Site2, col= "type_site2", remove=F) %>%
    ungroup() %>% 
    mutate(Distance.type = case_when(type_site1 == type_site2 ~ "Biol.replicates",
                                     typers1 == typers2 ~ "PCR.replicates",
                                    
                                     Miseq_run1 == Miseq_run2 ~ "Same Barcode Different Site", 
                                     TRUE ~ "Different Barcode")) %>%
    dplyr::select(Sample1 = Var1, Sample2 = Var2 , value , Distance.type) %>%
    filter (Sample1 != Sample2) -> all.distances.to.plot.pre
  
  
  
  # Checking all went well
  sapply(all.distances.to.plot.pre, function(x) summary(is.na(x)))
  
  all.distances.to.plot.pre$Distance.type <- all.distances.to.plot.pre$Distance.type  %>% fct_relevel( "PCR.replicates", "Biol.replicates", "Same Barcode Different Site", "Different Barcode")
  
  plot_5 <- ggplot (all.distances.to.plot.pre , aes (fill = Distance.type, x = value,after_stat(density))) +
    geom_histogram(stat = 'bin', alpha = 0.9, binwidth = 0.05) + #xlim(0, 1) +
    facet_wrap( ~ Distance.type) +
    labs (x = "Pairwise Dissimilarity", y = "Density" ,
          fill = "Groups", title = "eDNA Pairwise Dissimilarity Between Samples", subtitle = "Pre Occupancy") +theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                  panel.background = element_blank())
  ggsave(plot=plot_5,"Output_plots/eDNA_Pairwise_Dissimilarity_Between_Samples_pre_occupancy.png", device = "png", width = 12, height = 8, units = "in")


# Instead of chosing based on the pairwise distances, we will use distance to centroid

  #Need to Remove Single Sampels
  cleaned.tibble.pre_occ  %>% 
    group_by(Miseq_run,Site) %>% 
    dplyr::summarise(cases = n_distinct(Bio_rep)) %>% 
    filter(., cases == 1) %>% 
    unite(col="location_id", Miseq_run,Site, remove = FALSE, sep=":") %>% ungroup() %>% 
    dplyr::select(location_id) -> singles
  #Currently this method does not filter out these single tech reps. There are only 3
  
  #Nest Data
  cleaned.tibble.pre_occ  %>% 
    left_join(hash.key) %>% 
    unite(col="location_id", Miseq_run,Site, remove = FALSE, sep=":") %>% 
    filter(., !location_id %in% singles$location_id) %>% 
    group_by(location_id) %>% nest() -> nested.cleaning
  
  nested.cleaning %>% 
    mutate(matrix = purrr::map(data, tibble_to_vegdist_bio_rep_single)) -> nested.cleaning
  
  # Convert
  nested.cleaning %>% mutate(ncomparisons = purrr::map(matrix, length)) -> nested.cleaning
  
  #Calculate distances to centroid
  nested.cleaning <- nested.cleaning %>% mutate(distances = map2(matrix, location_id, dist_to_centroid))
  nested.cleaning %>% 
    separate(location_id, c("Miseq_run","Site"), sep=":") %>% 
    unnest_longer(distances) %>% 
    dplyr::select(Miseq_run,Site, distances_id,distances) -> all_distances.groups #unnest data



#Calculate normal distribution of distances to centroid
normparams.step5 <- MASS::fitdistr(all_distances.groups$distances, "normal")$estimate

#Calculate Probability
probs.step5 <- pnorm(all_distances.groups$distances, normparams.step5[1], normparams.step5[2])

#Determine Outliers
outliers.step5 <- which(probs.step5>max_dist_probability_cutoff)

#Remove Outliers
discard.step5 <-names(all_distances.groups$distances[outliers.step5])

to_write_discarded.step5 <- tibble(distances_id = discard.step5,
                                     distance = all_distances.groups$distances[outliers.step5])

to_write_discarded.step5 <- to_write_discarded.step5 %>% left_join(tibble(distances_id = discard.step5,
                                                                              probs = probs.step5[outliers.step5]))
write_csv(to_write_discarded.step5 ,"step5.discared_samples.csv")

## Plot Final Replication Levels

  cleaned.tibble.pre_occ %>% 
    unite(col="distances_id", Miseq_run,Site,Bio_rep, remove = FALSE, sep=":") ->cleaned.tibble.pre_occ_for_plotting
  all_distances.groups %>%
    dplyr::select(-Miseq_run,-Site) -> all_distances.groups_for_plotting
  
  cleaned.tibble.pre_occ_for_plotting %>%
    filter(., !distances_id %in% discard.step5) %>% 
    dplyr::group_by(Miseq_run,Site,Bio_rep) %>% 
    dplyr::summarise(cases = n_distinct(distances_id)) %>% 
    ggplot()+
    geom_raster(aes(x= Site, y = Bio_rep, fill = cases))+
    geom_text(aes(x= Site, y = Bio_rep, label = cases),color="white") +
    ggtitle("Sample Replication Level") +facet_wrap(~Miseq_run) +
    theme(axis.text.x = element_text(angle = 90))-> plot_7
  
  ggsave(plot=plot_7,"Output_plots/Final_sample_replication_level.png", device = "png", width = 20, height = 8, units = "in")


##Remove Samples

  ASV.nested %>% 
    mutate (Step3.tibble_edited = purrr::map(Step.3.low.read.depth, 
                                             function(.x){
                                               .x %>%
                                                 mutate(., Miseq_run= str_detect(seq_number,barcode_1)) %>% 
                                                 mutate(distances_id =TRUE )
                                             }) ) -> ASV.nested
  
  
  #Filter Sample
  ASV.nested %>% 
    mutate(Step4.tibble = purrr::map (Step3.tibble_edited,  ~ filter(.,! distances_id %in% to_write_discarded.step5$distances_id))) -> ASV.nested


#Visualize Results of Clearence Process 4
ASV.nested %>% 
  transmute(Miseq_run, Summary.1 = purrr::map(Step4.tibble, ~ how.many(ASVtable = .,round ="4.Dissimilarity"))) %>% 
  left_join(ASV.summary) %>%
  mutate(Summary   = map2(Summary, Summary.1, bind_rows)) %>%
  dplyr::select(-Summary.1) -> ASV.summary


ASV.summary$Summary[[3]] %>%  View()

#---

# Cleaning Process 5: **Remove Santa Barbara Samples** ----------------------------------------------------------------------------------------------------------------

###Create list of control samples
metadata %>% 
  filter(Site=="SB") %>% 
  dplyr::select(New_name) %>% unique() -> sb_samples
sb_samples <- sb_samples$New_name


#Filter Sample
ASV.nested %>% 
  mutate(Step5.tibble = purrr::map (Step4.tibble,  ~ filter(.,! sample %in% sb_samples))) -> ASV.nested

##Save Data
saveRDS(ASV.nested,file="ASV.nested_final.RDS")
saveRDS(ASV.summary,file="ASV.summary_final.RDS")


ASV.nested <- readRDS(file="ASV.nested_final.RDS")
ASV.summary <- readRDS(file="ASV.summary_final.RDS")
#---

#---

#Generate Final Outputs ----------------------------------------------------------------------------------------------------------------

##Code for merging ASV tables

###Identify Unique Species for Sum.taxonomy merging

hash.key %>% 
  distinct(.,sum.taxonomy) -> hashes_unique

hashes_unique$number <- row.names(hashes_unique)
hashes_unique$number <- paste0("taxon_",hashes_unique$number)
row.names(hashes_unique)<-hashes_unique$number

hash.key %>% 
  left_join(hashes_unique, by="sum.taxonomy") -> hash.key.updated

#---

#---

### Pre Occupancy Sum by Taxonomy, All PCR Tech Reps Separate Samples

hash.key.updated$number %>% unique() -> total_taxa

ASV.nested$Step5.tibble[[1]] %>% 
  filter(.,! sample %in% sb_samples) %>% 
  mutate(miseq = ASV.nested$Miseq_run[[1]]) %>% 
  unite(miseq,sample, col="Sample") %>% 
  left_join(hash.key.updated) %>% 
  dplyr::group_by(number,Sample) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Sample, nReads) %>% #convert to wide data format
  replace(is.na(.), 0) -> barcode_1_pre_df
barcode_1_pre_df$number %>%  unique() -> barcode_1_pre_df_taxa

ASV.nested$Step5.tibble[[2]] %>% 
  filter(.,! sample %in% sb_samples) %>% 
  mutate(miseq = ASV.nested$Miseq_run[[2]]) %>% 
  unite(miseq,sample, col="Sample") %>% 
  left_join(hash.key.updated) %>% 
  dplyr::group_by(number,Sample) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Sample, nReads) %>% #convert to wide data format
  replace(is.na(.), 0) -> barcode_2_pre_df
barcode_2_pre_df$number %>%  unique() -> barcode_2_pre_df_taxa

ASV.nested$Step5.tibble[[3]] %>% 
  filter(.,! sample %in% sb_samples) %>% 
  mutate(miseq = ASV.nested$Miseq_run[[3]]) %>% 
  unite(miseq,sample, col="Sample") %>% 
  left_join(hash.key.updated) %>% 
  dplyr::group_by(number,Sample) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Sample, nReads) %>% #convert to wide data format
  replace(is.na(.), 0) -> barcode_3_pre_df
barcode_3_pre_df$number %>%  unique() -> barcode_3_pre_df_taxa

total_kept_taxa <- (append(append(barcode_1_pre_df_taxa,barcode_2_pre_df_taxa),barcode_3_pre_df_taxa)) %>% unique()

if (rlang::is_empty(setdiff(total_kept_taxa,barcode_1_pre_df_taxa))) {
  barcode_1_pre_df %>% ungroup() %>% 
    arrange(number) -> barcode_1_pre_df
} else {
  barcode_1_pre_df %>% ungroup() %>% 
    add_row(number=setdiff(total_kept_taxa,barcode_1_pre_df_taxa)) %>% 
    arrange(number) %>% 
    replace(is.na(.), 0) -> barcode_1_pre_df
}

if (rlang::is_empty(setdiff(total_kept_taxa,barcode_2_pre_df_taxa))) {
  barcode_2_pre_df %>% ungroup() %>% 
    arrange(number) -> barcode_2_pre_df
} else {
  barcode_2_pre_df %>% ungroup() %>% 
    add_row(number=setdiff(total_kept_taxa,barcode_2_pre_df_taxa)) %>% 
    arrange(number) %>% 
    replace(is.na(.), 0) -> barcode_2_pre_df
}


if (rlang::is_empty(setdiff(total_kept_taxa,barcode_3_pre_df_taxa))) {
  barcode_3_pre_df %>% ungroup() %>% 
    arrange(number) -> barcode_3_pre_df
} else {
  barcode_3_pre_df %>% ungroup() %>% 
    add_row(number=setdiff(total_kept_taxa,barcode_3_pre_df_taxa)) %>% 
    arrange(number) %>% 
    replace(is.na(.), 0) -> barcode_3_pre_df
}


barcode_1_pre_df <- as.data.frame(barcode_1_pre_df)
row.names(barcode_1_pre_df) <- barcode_1_pre_df$number
barcode_1_pre_df %>% ungroup() %>% dplyr::select(-number) -> barcode_1_pre_df

barcode_2_pre_df <- as.data.frame(barcode_2_pre_df)
row.names(barcode_2_pre_df) <- barcode_2_pre_df$number
barcode_2_pre_df %>% ungroup() %>% dplyr::select(-number) -> barcode_2_pre_df

barcode_3_pre_df <- as.data.frame(barcode_3_pre_df)
row.names(barcode_3_pre_df) <- barcode_3_pre_df$number
barcode_3_pre_df %>% ungroup() %>% dplyr::select(-number) -> barcode_3_pre_df

####First, we want to create proportions by dividing by the rowsums:
####we could do this with sweep() or mutate_all() or other ways, but using vegan:

barcode_1_prop <- decostand(barcode_1_pre_df, method = "total", MARGIN = 2)
barcode_2_prop <- decostand(barcode_2_pre_df, method = "total", MARGIN = 2)
barcode_3_prop <- decostand(barcode_3_pre_df, method = "total", MARGIN = 2)

####Second, we want to ask how the proprortion for each species has changed across columns (samples). 
####We do this by scaling everything to the max observed in each row. 

####to do this WITHIN a dataset, we could just do (again, using vegan):
barcode_1_prop_index <- decostand(barcode_1_prop, method = "max", MARGIN = 1)
barcode_2_prop_index <- decostand(barcode_2_prop, method = "max", MARGIN = 1)
barcode_3_prop_index <- decostand(barcode_3_prop, method = "max", MARGIN = 1)

####This gives us an index between 0 and 1 for each species in each dataset.  

####But if we want to combine datasets, this second step has to happen in the combined dataset, so it all gets scaled to 0-1.  
####easy enough:

combined_index <- decostand(cbind(barcode_1_prop,barcode_2_prop, barcode_3_prop), method = "max", MARGIN = 1)
####How both datasets are combined, on a common, comparable scale.

### Output Read Count Data
pre_results_reads = cbind(barcode_1_prop,barcode_2_prop,barcode_3_prop)

hash.key.updated.2 <- hash.key.updated[!duplicated(hash.key.updated$number), ]

pre_results_reads$number <- rownames(pre_results_reads)

pre_results_reads %>% 
  left_join(hash.key.updated.2, by="number") %>% 
  dplyr::select(-number,-seq_number) -> pre_results_reads

saveRDS(pre_results_reads,file="Output_R/pre_occupancy_results_sum.taxonomy_tech_reps_separate_read_counts.RDS")
write_csv(pre_results_reads ,"Output_csv/pre_occupancy_results_sum.taxonomy_tech_reps_separate_read_counts.csv")

### Output eDNA Index Data

combined_index$number <- rownames(combined_index)

combined_index %>% 
  left_join(hash.key.updated.2, by="number") %>% 
  dplyr::select(-number,-seq_number) -> combined_index

saveRDS(combined_index,file="Output_R/pre_occupancy_results_sum.taxonomy_tech_reps_separate_eDNA_index.RDS")
write_csv(combined_index ,"Output_csv/pre_occupancy_results_sum.taxonomy_tech_reps_separate_eDNA_index.csv")

#---

#---

#Code for Merging Tech Reps

###Pre Occupancy Sum by Taxonomy, Biological Replicates Separate

ASV.nested$Step5.tibble[[1]] %>%
  filter(.,! sample %in% sb_samples) %>% 
  mutate(Sample = sample) %>% 
  left_join(hash.key.updated) %>% 
  dplyr::group_by(number,Sample) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Sample, nReads) %>% #convert to wide data format
  replace(is.na(.), 0)  -> barcode_1_pre_df
barcode_1_pre_df$number %>%  unique() -> barcode_1_pre_df_taxa

metadata %>% 
  filter(., !(New_name %in% controls)) %>% 
  filter(.,!(New_name %in% colnames(barcode_1_pre_df))) %>% 
  filter(.,!(New_name %in% sb_samples)) %>% 
  pull(New_name) %>%  unique()-> columns2add

barcode_1_pre_df <- as.data.frame(barcode_1_pre_df)

barcode_1_pre_df %>% 
  tibble::add_column(!!!set_names(as.list(rep(NA, length(columns2add))),nm=columns2add)) %>% 
  replace(is.na(.), 0) %>% 
  dplyr::select(sort(tidyselect::peek_vars()))-> barcode_1_pre_df

row.names(barcode_1_pre_df) <- barcode_1_pre_df$number
barcode_1_pre_df %>% ungroup() -> barcode_1_pre_df

ASV.nested$Step5.tibble[[2]] %>% 
  filter(.,! sample %in% sb_samples) %>% 
  mutate(Sample = sample) %>% 
  left_join(hash.key.updated) %>% 
  dplyr::group_by(number,Sample) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Sample, nReads) %>% #convert to wide data format
  replace(is.na(.), 0)  -> barcode_2_pre_df
barcode_2_pre_df$number %>%  unique() -> barcode_2_pre_df_taxa

metadata %>% 
  filter(., !(New_name %in% controls)) %>% 
  filter(.,!(New_name %in% colnames(barcode_2_pre_df))) %>% 
  filter(.,!(New_name %in% sb_samples)) %>% 
  pull(New_name) %>%  unique()-> columns2add

barcode_2_pre_df <- as.data.frame(barcode_2_pre_df)

barcode_2_pre_df %>% 
  tibble::add_column(!!!set_names(as.list(rep(NA, length(columns2add))),nm=columns2add)) %>% 
  replace(is.na(.), 0) %>% 
  dplyr::select(sort(tidyselect::peek_vars())) -> barcode_2_pre_df

row.names(barcode_2_pre_df) <- barcode_2_pre_df$number
barcode_2_pre_df %>% ungroup() -> barcode_2_pre_df

ASV.nested$Step5.tibble[[3]] %>% 
  filter(.,! sample %in% sb_samples) %>% 
  mutate(Sample = sample) %>% 
  left_join(hash.key.updated) %>% 
  dplyr::group_by(number,Sample) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Sample, nReads) %>% #convert to wide data format
  replace(is.na(.), 0)  -> barcode_3_pre_df
barcode_3_pre_df$number %>%  unique() -> barcode_3_pre_df_taxa

metadata %>% 
  filter(., !(New_name %in% controls)) %>% 
  filter(.,!(New_name %in% colnames(barcode_3_pre_df))) %>%
  filter(.,!(New_name %in% sb_samples)) %>% 
  pull(New_name) %>%  unique()-> columns2add

barcode_3_pre_df <- as.data.frame(barcode_3_pre_df)

barcode_3_pre_df %>% 
  tibble::add_column(!!!set_names(as.list(rep(NA, length(columns2add))),nm=columns2add)) %>% 
  replace(is.na(.), 0) %>% 
  dplyr::select(sort(tidyselect::peek_vars())) -> barcode_3_pre_df

row.names(barcode_3_pre_df) <- barcode_3_pre_df$number
barcode_3_pre_df %>% ungroup() -> barcode_3_pre_df

total_kept_taxa <- append(append(barcode_1_pre_df_taxa,barcode_2_pre_df_taxa),barcode_3_pre_df_taxa) %>% unique()

if (rlang::is_empty(setdiff(total_kept_taxa,barcode_1_pre_df_taxa))) {
  barcode_1_pre_df %>% ungroup() %>% 
    arrange(number) -> barcode_1_pre_df
} else {
  barcode_1_pre_df %>% ungroup() %>% 
    add_row(number=setdiff(total_kept_taxa,barcode_1_pre_df_taxa)) %>% 
    arrange(number) %>% 
    replace(is.na(.), 0) -> barcode_1_pre_df
}

if (rlang::is_empty(setdiff(total_kept_taxa,barcode_2_pre_df_taxa))) {
  barcode_2_pre_df %>% ungroup() %>% 
    arrange(number) -> barcode_2_pre_df
} else {
  barcode_2_pre_df %>% ungroup() %>% 
    add_row(number=setdiff(total_kept_taxa,barcode_2_pre_df_taxa)) %>% 
    arrange(number) %>% 
    replace(is.na(.), 0) -> barcode_2_pre_df
}


if (rlang::is_empty(setdiff(total_kept_taxa,barcode_3_pre_df_taxa))) {
  barcode_3_pre_df %>% ungroup() %>% 
    arrange(number) -> barcode_3_pre_df
} else {
  barcode_3_pre_df %>% ungroup() %>% 
    add_row(number=setdiff(total_kept_taxa,barcode_3_pre_df_taxa)) %>% 
    arrange(number) %>% 
    replace(is.na(.), 0) -> barcode_3_pre_df
}

barcode_1_pre_df <- as.data.frame(barcode_1_pre_df)
row.names(barcode_1_pre_df) <- barcode_1_pre_df$number
barcode_1_pre_df %>% ungroup() %>% dplyr::select(-number) -> barcode_1_pre_df

barcode_2_pre_df <- as.data.frame(barcode_2_pre_df)
row.names(barcode_2_pre_df) <- barcode_2_pre_df$number
barcode_2_pre_df %>% ungroup() %>% dplyr::select(-number) -> barcode_2_pre_df

barcode_3_pre_df <- as.data.frame(barcode_3_pre_df)
row.names(barcode_3_pre_df) <- barcode_3_pre_df$number
barcode_3_pre_df %>% ungroup() %>% dplyr::select(-number) -> barcode_3_pre_df

barcode_1_pre_df+barcode_2_pre_df+barcode_3_pre_df -> combined_df

####First, we want to create proportions by dividing by the rowsums:
####we could do this with sweep() or mutate_all() or other ways, but using vegan:

combined_df_prop <- decostand(combined_df, method = "total", MARGIN = 2)

####Second, we want to ask how the proprortion for each species has changed across columns (samples). 
####We do this by scaling everything to the max observed in each row. 

####To do this WITHIN a dataset, we could just do (again, using vegan):
####eDNA Index between 0-1 made by straight adding reads across each decontaminated datatable
combined_df_index <- decostand(combined_df_prop, method = "max", MARGIN = 1)

#Output read counts
hash.key.updated.2 <- hash.key.updated[!duplicated(hash.key.updated$number), ]

combined_df$number <- rownames(combined_df)

combined_df %>% 
  left_join(hash.key.updated.2, by="number") %>% 
  dplyr::select(-number,-seq_number) -> combined_df

saveRDS(combined_df,file="Output_R/pre_occupancy_results_sum.taxonomy_tech_reps_summed_read_counts.RDS")
write_csv(combined_df,"Output_csv/pre_occupancy_results_sum.taxonomy_tech_reps_summed_read_counts.csv")

### Output eDNA Index
combined_df_index$number <- rownames(combined_df_index)

combined_df_index %>% 
  left_join(hash.key.updated.2, by="number") %>% 
  dplyr::select(-number,-seq_number) -> combined_df_index

saveRDS(combined_df_index,file="Output_R/pre_occupancy_results_sum.taxonomy_tech_reps_summed_eDNA_index.RDS")
write_csv(combined_df_index,"Output_csv/pre_occupancy_results_sum.taxonomy_tech_reps_summed_eDNA_index.csv")

#---

#---

#Code for Merging Sites ---

###Pre Occupancy Merge by Site
ASV.nested %>% 
  dplyr::select(Step5.tibble,Miseq_run) %>% 
  unnest(Step5.tibble) %>% 
  filter(.,! sample %in% sb_samples) %>% 
  left_join(hash.key.updated) %>% 
  left_join(metadata, by=c("sample"="New_name")) %>% 
  dplyr::group_by(sum.taxonomy,Site) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Site, nReads) %>% #convert to wide data format
  replace(is.na(.), 0)  -> pre_wide_reads

saveRDS(pre_wide_reads,file="Output_R/pre_occ_site_sum.taxonomy_reads_summed.RDS")
write_csv(pre_wide_reads ,"Output_csv/pre_occ_site_sum.taxonomy_reads_summed.csv")

ASV.nested %>% 
  dplyr::select(Step5.tibble,Miseq_run) %>% 
  unnest(Step5.tibble) %>% 
  filter(.,! sample %in% sb_samples) %>% 
  left_join(hash.key.updated) %>% 
  left_join(metadata, by=c("sample"="New_name")) %>% 
  ungroup() %>% 
  group_by(Site, sum.taxonomy) %>% 
  summarise(meanreads = mean(nReads)) %>% 
  dplyr::group_by(Site) %>%
  mutate (Tot = sum(meanreads),
          Row.sums = meanreads / Tot) %>% 
  dplyr::group_by (sum.taxonomy) %>%
  mutate (Colmax = max (Row.sums),
          Normalized.reads = Row.sums / Colmax) %>% 
  ungroup() %>% 
  dplyr::select(Site,Normalized.reads, sum.taxonomy) %>% 
  spread(., Site, Normalized.reads) %>% #convert to wide data format
  replace(is.na(.), 0)  -> pre_wide

saveRDS(pre_wide,file="Output_R/pre_occ_site_averaged_sum.taxonomy_e_index.RDS")
write_csv(pre_wide,"Output_csv/pre_occ_site_averaged_sum.taxonomy_e_index.csv")

metadata <- read.table(input_meta_path, header = TRUE, sep = "\t", stringsAsFactors = F)
metadata %>%  dplyr::select(-Seq_number) %>% distinct() -> metadata2

#### ASV Tables
###Pre Occupancy ASV, All PCR Tech Reps Separate Samples
ASV.nested %>% 
  dplyr::select(Step5.tibble) %>% 
  unnest(Step5.tibble) %>% 
  filter(.,! sample %in% sb_samples) %>% 
  ungroup() %>% 
  left_join(hash.key.updated) %>% 
  left_join(metadata2, by=c("sample"="New_name")) %>% 
  dplyr::select(sample,seq_number, nReads, sum.taxonomy) %>%
  spread(., sample, nReads) %>% #convert to wide data format
  replace(is.na(.), 0)  -> pre_wide_reads

saveRDS(pre_wide_reads,file="Output_R/pre_occ_site_asv_table.RDS")
write_csv(pre_wide_reads ,"Output_csv/pre_occ_site_asv_table.csv")

ASV.nested %>% 
  dplyr::select(Step5.tibble) %>% 
  unnest(Step5.tibble) %>% 
  ungroup() %>% 
  left_join(hash.key.updated) %>% 
  left_join(metadata2, by=c("sample"="New_name")) %>% 
  dplyr::select(sample,seq_number, nReads,sum.taxonomy) %>%
  dplyr::group_by(sample) %>%
  mutate (Tot = sum(nReads),
          Row.sums = nReads / Tot) %>% 
  dplyr::group_by (seq_number) %>%
  mutate (Colmax = max(Row.sums),
          Normalized.reads = Row.sums / Colmax) %>% 
  dplyr::select(sample,Normalized.reads, seq_number, sum.taxonomy) %>% 
  spread(., sample, Normalized.reads) %>% #convert to wide data format
  replace(is.na(.), 0)  -> pre_wide

saveRDS(pre_wide,file="Output_R/pre_occ_site_asv_table_edna_index.RDS")
write_csv(pre_wide ,"Output_csv/pre_occ_site_asv_table_edna_index.csv")


# Generate 12S Fish Only --------------------


##Code for merging ASV tables

###Identify Unique Species for Sum.taxonomy merging

hash.key %>% 
  filter(., str_detect(seq_number, "12S")) %>% 
  distinct(.,sum.taxonomy) -> hashes_unique_12S

hashes_unique_12S$number <- row.names(hashes_unique_12S)
hashes_unique_12S$number <- paste0("taxon_",hashes_unique_12S$number)
row.names(hashes_unique_12S)<-hashes_unique_12S$number

hash.key %>% 
  filter(., str_detect(seq_number, "12S")) %>% 
  left_join(hashes_unique_12S, by="sum.taxonomy") -> hash.key.updated_12S

saveRDS(hash.key.updated_12S, "hash.key.updated_12S")
#---

#---

### Pre Occupancy Sum by Taxonomy, All PCR Tech Reps Separate Samples

hash.key.updated_12S$number %>% unique() -> total_taxa

ASV.nested$Step5.tibble[[2]] %>% 
  filter(.,! sample %in% sb_samples) %>% 
  mutate(miseq = ASV.nested$Miseq_run[[2]]) %>% 
  unite(miseq,sample, col="Sample") %>% 
  left_join(hash.key.updated_12S) %>% 
  dplyr::group_by(number,Sample) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Sample, nReads) %>% #convert to wide data format
  replace(is.na(.), 0) -> barcode_2_pre_df
barcode_2_pre_df$number %>%  unique() -> barcode_2_pre_df_taxa

ASV.nested$Step5.tibble[[3]] %>% 
  filter(.,! sample %in% sb_samples) %>% 
  mutate(miseq = ASV.nested$Miseq_run[[3]]) %>% 
  unite(miseq,sample, col="Sample") %>% 
  left_join(hash.key.updated_12S) %>% 
  dplyr::group_by(number,Sample) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Sample, nReads) %>% #convert to wide data format
  replace(is.na(.), 0) -> barcode_3_pre_df
barcode_3_pre_df$number %>%  unique() -> barcode_3_pre_df_taxa

total_kept_taxa <- (append(barcode_2_pre_df_taxa,barcode_3_pre_df_taxa)) %>% unique()

if (rlang::is_empty(setdiff(total_kept_taxa,barcode_2_pre_df_taxa))) {
  barcode_2_pre_df %>% ungroup() %>% 
    arrange(number) -> barcode_2_pre_df
} else {
  barcode_2_pre_df %>% ungroup() %>% 
    add_row(number=setdiff(total_kept_taxa,barcode_2_pre_df_taxa)) %>% 
    arrange(number) %>% 
    replace(is.na(.), 0) -> barcode_2_pre_df
}


if (rlang::is_empty(setdiff(total_kept_taxa,barcode_3_pre_df_taxa))) {
  barcode_3_pre_df %>% ungroup() %>% 
    arrange(number) -> barcode_3_pre_df
} else {
  barcode_3_pre_df %>% ungroup() %>% 
    add_row(number=setdiff(total_kept_taxa,barcode_3_pre_df_taxa)) %>% 
    arrange(number) %>% 
    replace(is.na(.), 0) -> barcode_3_pre_df
}

barcode_2_pre_df <- as.data.frame(barcode_2_pre_df)
row.names(barcode_2_pre_df) <- barcode_2_pre_df$number
barcode_2_pre_df %>% ungroup() %>% dplyr::select(-number) -> barcode_2_pre_df

barcode_3_pre_df <- as.data.frame(barcode_3_pre_df)
row.names(barcode_3_pre_df) <- barcode_3_pre_df$number
barcode_3_pre_df %>% ungroup() %>% dplyr::select(-number) -> barcode_3_pre_df

####First, we want to create proportions by dividing by the rowsums:
####we could do this with sweep() or mutate_all() or other ways, but using vegan:

barcode_2_prop <- decostand(barcode_2_pre_df, method = "total", MARGIN = 2)
barcode_3_prop <- decostand(barcode_3_pre_df, method = "total", MARGIN = 2)

####Second, we want to ask how the proprortion for each species has changed across columns (samples). 
####We do this by scaling everything to the max observed in each row. 

####to do this WITHIN a dataset, we could just do (again, using vegan):
barcode_2_prop_index <- decostand(barcode_2_prop, method = "max", MARGIN = 1)
barcode_3_prop_index <- decostand(barcode_3_prop, method = "max", MARGIN = 1)

####This gives us an index between 0 and 1 for each species in each dataset.  

####But if we want to combine datasets, this second step has to happen in the combined dataset, so it all gets scaled to 0-1.  
####easy enough:

combined_index <- decostand(cbind(barcode_2_prop, barcode_3_prop), method = "max", MARGIN = 1)
####How both datasets are combined, on a common, comparable scale.

### Output Read Count Data
pre_results_reads = cbind(barcode_2_prop,barcode_3_prop)

hash.key.updated_12S.2 <- hash.key.updated_12S[!duplicated(hash.key.updated_12S$number), ]

pre_results_reads$number <- rownames(pre_results_reads)

pre_results_reads %>% 
  left_join(hash.key.updated_12S.2, by="number") %>% 
  dplyr::select(-number,-seq_number) -> pre_results_reads

saveRDS(pre_results_reads,file="Output_R/fish_12S_pre_occupancy_results_sum.taxonomy_tech_reps_separate_read_counts.RDS")
write_csv(pre_results_reads ,"Output_csv/fish_12S_pre_occupancy_results_sum.taxonomy_tech_reps_separate_read_counts.csv")

### Output eDNA Index Data

combined_index$number <- rownames(combined_index)

combined_index %>% 
  left_join(hash.key.updated_12S.2, by="number") %>% 
  dplyr::select(-number,-seq_number) -> combined_index

saveRDS(combined_index,file="Output_R/fish_12S_pre_occupancy_results_sum.taxonomy_tech_reps_separate_eDNA_index.RDS")
write_csv(combined_index ,"Output_csv/fish_12S_pre_occupancy_results_sum.taxonomy_tech_reps_separate_eDNA_index.csv")

#---

#---

#Code for Merging Tech Reps

###Pre Occupancy Sum by Taxonomy, Biological Replicates Separate

ASV.nested$Step5.tibble[[2]] %>% 
  filter(.,! sample %in% sb_samples) %>% 
  mutate(Sample = sample) %>% 
  left_join(hash.key.updated_12S) %>% 
  dplyr::group_by(number,Sample) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Sample, nReads) %>% #convert to wide data format
  replace(is.na(.), 0)  -> barcode_2_pre_df
barcode_2_pre_df$number %>%  unique() -> barcode_2_pre_df_taxa

metadata %>% 
  filter(., !(New_name %in% controls)) %>% 
  filter(.,!(New_name %in% colnames(barcode_2_pre_df))) %>% 
  filter(.,!(New_name %in% sb_samples)) %>% 
  pull(New_name) %>%  unique()-> columns2add

barcode_2_pre_df <- as.data.frame(barcode_2_pre_df)

barcode_2_pre_df %>% 
  tibble::add_column(!!!set_names(as.list(rep(NA, length(columns2add))),nm=columns2add)) %>% 
  replace(is.na(.), 0) %>% 
  dplyr::select(sort(tidyselect::peek_vars())) -> barcode_2_pre_df

row.names(barcode_2_pre_df) <- barcode_2_pre_df$number
barcode_2_pre_df %>% ungroup() -> barcode_2_pre_df

ASV.nested$Step5.tibble[[3]] %>% 
  filter(.,! sample %in% sb_samples) %>% 
  mutate(Sample = sample) %>% 
  left_join(hash.key.updated_12S) %>% 
  dplyr::group_by(number,Sample) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Sample, nReads) %>% #convert to wide data format
  replace(is.na(.), 0)  -> barcode_3_pre_df
barcode_3_pre_df$number %>%  unique() -> barcode_3_pre_df_taxa

metadata %>% 
  filter(., !(New_name %in% controls)) %>% 
  filter(.,!(New_name %in% colnames(barcode_3_pre_df))) %>% 
  filter(.,!(New_name %in% sb_samples)) %>% 
  pull(New_name) %>%  unique()-> columns2add

barcode_3_pre_df <- as.data.frame(barcode_3_pre_df)

barcode_3_pre_df %>% 
  tibble::add_column(!!!set_names(as.list(rep(NA, length(columns2add))),nm=columns2add)) %>% 
  replace(is.na(.), 0) %>% 
  dplyr::select(sort(tidyselect::peek_vars())) -> barcode_3_pre_df

row.names(barcode_3_pre_df) <- barcode_3_pre_df$number
barcode_3_pre_df %>% ungroup() -> barcode_3_pre_df

total_kept_taxa <- append(barcode_2_pre_df_taxa,barcode_3_pre_df_taxa) %>% unique()

if (rlang::is_empty(setdiff(total_kept_taxa,barcode_2_pre_df_taxa))) {
  barcode_2_pre_df %>% ungroup() %>% 
    arrange(number) -> barcode_2_pre_df
} else {
  barcode_2_pre_df %>% ungroup() %>% 
    add_row(number=setdiff(total_kept_taxa,barcode_2_pre_df_taxa)) %>% 
    arrange(number) %>% 
    replace(is.na(.), 0) -> barcode_2_pre_df
}


if (rlang::is_empty(setdiff(total_kept_taxa,barcode_3_pre_df_taxa))) {
  barcode_3_pre_df %>% ungroup() %>% 
    arrange(number) -> barcode_3_pre_df
} else {
  barcode_3_pre_df %>% ungroup() %>% 
    add_row(number=setdiff(total_kept_taxa,barcode_3_pre_df_taxa)) %>% 
    arrange(number) %>% 
    replace(is.na(.), 0) -> barcode_3_pre_df
}

barcode_2_pre_df <- as.data.frame(barcode_2_pre_df)
row.names(barcode_2_pre_df) <- barcode_2_pre_df$number
barcode_2_pre_df %>% ungroup() %>% dplyr::select(-number) -> barcode_2_pre_df

barcode_3_pre_df <- as.data.frame(barcode_3_pre_df)
row.names(barcode_3_pre_df) <- barcode_3_pre_df$number
barcode_3_pre_df %>% ungroup() %>% dplyr::select(-number) -> barcode_3_pre_df

barcode_2_pre_df+barcode_3_pre_df -> combined_df

####First, we want to create proportions by dividing by the rowsums:
####we could do this with sweep() or mutate_all() or other ways, but using vegan:

combined_df_prop <- decostand(combined_df, method = "total", MARGIN = 2)

####Second, we want to ask how the proprortion for each species has changed across columns (samples). 
####We do this by scaling everything to the max observed in each row. 

####To do this WITHIN a dataset, we could just do (again, using vegan):
####eDNA Index between 0-1 made by straight adding reads across each decontaminated datatable
combined_df_index <- decostand(combined_df_prop, method = "max", MARGIN = 1)

#Output read counts
hash.key.updated_12S.2 <- hash.key.updated_12S[!duplicated(hash.key.updated_12S$number), ]

combined_df$number <- rownames(combined_df)

combined_df %>% 
  left_join(hash.key.updated_12S.2, by="number") %>% 
  dplyr::select(-number,-seq_number) -> combined_df

saveRDS(combined_df,file="Output_R/fish_12S_pre_occupancy_results_sum.taxonomy_tech_reps_summed_read_counts.RDS")
write_csv(combined_df,"Output_csv/fish_12S_pre_occupancy_results_sum.taxonomy_tech_reps_summed_read_counts.csv")

### Output eDNA Index
combined_df_index$number <- rownames(combined_df_index)

combined_df_index %>% 
  left_join(hash.key.updated_12S.2, by="number") %>% 
  dplyr::select(-number,-seq_number) -> combined_df_index

saveRDS(combined_df_index,file="Output_R/fish_12S_pre_occupancy_results_sum.taxonomy_tech_reps_summed_eDNA_index.RDS")
write_csv(combined_df_index,"Output_csv/fish_12S_pre_occupancy_results_sum.taxonomy_tech_reps_summed_eDNA_index.csv")


#---

#---

#Code for Merging Sites

###Pre Occupancy Merge by Site
ASV.nested %>% 
  dplyr::select(Step5.tibble,Miseq_run) %>% 
  unnest(Step5.tibble) %>% 
  filter(.,! sample %in% sb_samples) %>% 
  filter(., str_detect(seq_number, "12S")) %>% 
  left_join(hash.key.updated_12S) %>% 
  left_join(metadata, by=c("sample"="New_name")) %>% 
  dplyr::group_by(sum.taxonomy,Site) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Site, nReads) %>% #convert to wide data format
  replace(is.na(.), 0)  -> pre_wide_reads

saveRDS(pre_wide_reads,file="Output_R/fish_12S_pre_occ_site_sum.taxonomy_reads_summed.RDS")
write_csv(pre_wide_reads ,"Output_csv/fish_12S_pre_occ_site_sum.taxonomy_reads_summed.csv")

ASV.nested %>% 
  dplyr::select(Step5.tibble,Miseq_run) %>% 
  unnest(Step5.tibble) %>% 
  filter(.,! sample %in% sb_samples) %>% 
  filter(., str_detect(seq_number, "12S")) %>% 
  left_join(hash.key.updated_12S) %>% 
  left_join(metadata, by=c("sample"="New_name")) %>% 
  ungroup() %>% 
  group_by(Site, sum.taxonomy) %>% 
  summarise(meanreads = mean(nReads)) %>% 
  dplyr::group_by(Site) %>%
  mutate (Tot = sum(meanreads),
          Row.sums = meanreads / Tot) %>% 
  dplyr::group_by (sum.taxonomy) %>%
  mutate (Colmax = max (Row.sums),
          Normalized.reads = Row.sums / Colmax) %>% 
  ungroup() %>% 
  dplyr::select(Site,Normalized.reads, sum.taxonomy) %>% 
  spread(., Site, Normalized.reads) %>% #convert to wide data format
  replace(is.na(.), 0)  -> pre_wide

saveRDS(pre_wide,file="Output_R/fish_12S_pre_occ_site_averaged_sum.taxonomy_e_index.RDS")
write_csv(pre_wide,"Output_csv/fish_12S_pre_occ_site_averaged_sum.taxonomy_e_index.csv")

metadata <- read.table(input_meta_path, header = TRUE, sep = "\t", stringsAsFactors = F)
metadata %>%  dplyr::select(-Seq_number) %>% distinct() -> metadata2

###Pre Occupancy Merge by Site
ASV.nested %>% 
  dplyr::select(Step5.tibble) %>% 
  unnest(Step5.tibble) %>% 
  filter(., str_detect(seq_number, "12S")) %>% 
  ungroup() %>% 
  left_join(metadata2, by=c("sample"="New_name")) %>% 
  dplyr::group_by(sum.taxonomy,Site) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Site, nReads) %>% #convert to wide data format
  replace(is.na(.), 0)  -> pre_wide_reads

saveRDS(pre_wide_reads,file="Output_R/fish_12S_pre_occ_site_sum.taxonomy_reads_summed.RDS")
write_csv(pre_wide_reads ,"Output_csv/fish_12S_pre_occ_site_sum.taxonomy_reads_summed.csv")

ASV.nested %>% 
  dplyr::select(Step5.tibble) %>% 
  unnest(Step5.tibble) %>% 
  filter(., str_detect(seq_number, "12S")) %>% 
  ungroup() %>% 
  left_join(metadata2, by=c("sample"="New_name")) %>% 
  group_by(Site, sum.taxonomy) %>% 
  summarise(meanreads = mean(nReads)) %>% 
  dplyr::group_by(Site) %>%
  mutate (Tot = sum(meanreads),
          Row.sums = meanreads / Tot) %>% 
  dplyr::group_by (sum.taxonomy) %>%
  mutate (Colmax = max(Row.sums),
          Normalized.reads = Row.sums / Colmax) %>% 
  dplyr::select(Site,Normalized.reads, sum.taxonomy) %>% 
  spread(., Site, Normalized.reads) %>% #convert to wide data format
  replace(is.na(.), 0)  -> pre_wide

saveRDS(pre_wide,file="Output_R/fish_12S_pre_occ_site_averaged_sum.taxonomy_e_index.RDS")
write_csv(pre_wide,"Output_csv/fish_12S_pre_occ_site_averaged_sum.taxonomy_e_index.csv")


#### ASV Tables
###Pre Occupancy ASV, All PCR Tech Reps Separate Samples
ASV.nested %>% 
  dplyr::select(Step5.tibble) %>% 
  unnest(Step5.tibble) %>% 
  filter(.,! sample %in% sb_samples) %>% 
  filter(., str_detect(seq_number, "12S")) %>% 
  ungroup() %>% 
  left_join(hash.key.updated_12S) %>% 
  left_join(metadata2, by=c("sample"="New_name")) %>% 
  dplyr::select(sample,seq_number, nReads, sum.taxonomy) %>%
  spread(., sample, nReads) %>% #convert to wide data format
  replace(is.na(.), 0)  -> pre_wide_reads

saveRDS(pre_wide_reads,file="Output_R/fish_12S_pre_occ_site_asv_table.RDS")
write_csv(pre_wide_reads ,"Output_csv/fish_12S_pre_occ_site_asv_table.csv")

ASV.nested %>% 
  dplyr::select(Step5.tibble) %>% 
  unnest(Step5.tibble) %>% 
  filter(., str_detect(seq_number, "12S")) %>% 
  ungroup() %>% 
  left_join(hash.key.updated_12S) %>% 
  left_join(metadata2, by=c("sample"="New_name")) %>% 
  dplyr::select(sample,seq_number, nReads,sum.taxonomy) %>%
  dplyr::group_by(sample) %>%
  mutate (Tot = sum(nReads),
          Row.sums = nReads / Tot) %>% 
  dplyr::group_by (seq_number) %>%
  mutate (Colmax = max(Row.sums),
          Normalized.reads = Row.sums / Colmax) %>% 
  dplyr::select(sample,Normalized.reads, seq_number, sum.taxonomy) %>% 
  spread(., sample, Normalized.reads) %>% #convert to wide data format
  replace(is.na(.), 0)  -> pre_wide

saveRDS(pre_wide,file="Output_R/fish_12S_pre_occ_site_asv_table_edna_index.RDS")
write_csv(pre_wide ,"Output_csv/fish_12S_pre_occ_site_asv_table_edna_index.csv")

###Pre Occupancy Merge by Site
ASV.nested %>% 
  dplyr::select(Step5.tibble) %>% 
  unnest(Step5.tibble) %>% 
  filter(., str_detect(seq_number, "12S")) %>% 
  ungroup() %>% 
  left_join(metadata2, by=c("sample"="New_name")) %>% 
  dplyr::group_by(seq_number,Site) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Site, nReads) %>% #convert to wide data format
  replace(is.na(.), 0) %>% 
  left_join(hash.key.updated_12S)-> pre_wide_reads

saveRDS(pre_wide_reads,file="Output_R/fish_12S_pre_occ_site_ASV_reads_summed.RDS")
write_csv(pre_wide_reads ,"Output_csv/fish_12S_pre_occ_site_ASV_reads_summed.csv")


###Pre Occupancy Merge by Site_Viz
ASV.nested %>% 
  dplyr::select(Step5.tibble) %>% 
  unnest(Step5.tibble) %>% 
  filter(., str_detect(seq_number, "12S")) %>% 
  ungroup() %>% 
  left_join(metadata2, by=c("sample"="New_name")) %>% 
  dplyr::group_by(seq_number,Site_Viz) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Site_Viz, nReads) %>% #convert to wide data format
  replace(is.na(.), 0) %>% 
  left_join(hash.key.updated_12S)-> pre_wide_reads

saveRDS(pre_wide_reads,file="Output_R/fish_12S_pre_occ_site_ASV_reads_summed_4_viz.RDS")
write_csv(pre_wide_reads ,"Output_csv/fish_12S_pre_occ_site_ASV_reads_summed_4_viz.RDS")

