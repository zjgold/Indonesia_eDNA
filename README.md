# Environmental DNA in a Global Biodiversity Hotspot: Lessons from Coral Reef Fish Diversity Across the Indonesian Archipelago

Pre-print is Available [here](https://www.biorxiv.org/content/10.1101/2021.02.19.432056v1)


Onny Marwayana<sup>1</sup>, Zachary Gold<sup>1</sup>,Christopher Meyer<sup>2</sup>,Paul H. Barber<sup>1</sup>


<sup>1</sup>Department of Ecology and Evolutionary Biology, UCLA, Los Angeles

<sup>2</sup>Smithonsian Institute, Washington D.C,



## Abstract
Indonesia is the heart of the Coral Triangle, the world’s most diverse marine ecosystem. Preserving the biological and economic value of this marine biodiversity requires efficient and economical ecosystem monitoring, yet our understanding of marine biodiversity in this region remains limited. Towards this end, this study uses environmental DNA (eDNA) to survey fish communities across a well-documented biodiversity gradient in Indonesia. A total of  6,608,693 sequence reads of MiFish 12S rRNA from 39 sites spanning 7 regions of Indonesia revealed 1,099 fish Amplified Sequence Variants (ASVs), 80.4% of which could be identified to species through the inclusion of new reference sequences from Mo’orea BIOCODE. Patterns of regional fish diversity inferred from eDNA broadly conformed to expectations, with the highest fish biodiversity in Raja Ampat and lowest in Western Indonesia. Similarly, zeta diversity analysis showed greater community turnover in higher diversity reefs of Eastern Indonesia, and greater community similarity in low diversity regions of Western Indonesia. However, several results highlight challenges for eDNA in megadiverse ecosystems. Despite a two-fold difference in fish diversity between Eastern and Western Indonesia, mean ASVs recovered per one-liter seawater was relatively similar across 7 regions of Indonesia. Moreover, although ASV recovery from individual seawater samples saturated, ASV recovery did not saturate at the site or region level, indicating that sampling/sequencing efforts employed in lower diversity temperate marine ecosystems are insufficient for biodiversity hotspots like the Coral Triangle. Despite these limitations, 36.3% to 84.1% (mean 57.1%) of fish species detected by eDNA at 8 sites within Raja Ampat were not observed during intensive visual surveys. Taxa missed include pelagic (tuna, jacks, scads, mackerels), nocturnal (soldierfish, lanternfish), and crevice dwelling species (eels, blennies, gobies) that are difficult to document in visual surveys. These results demonstrate the added value of eDNA in biodiversity hotspots like the Coral Triangle and the need for further research to understand how best to sample eDNA in high diversity regions like the Coral Triangle to deliver on the promise of eDNA as a tool to monitor marine biodiversity effectively and efficiently.   

## Description
This page is dedicated to hosting code generated for the Accepted Environmental DNA Manuscript. 

Included on this page is
1. **Scripts used to Conduct Analyses**
    1. *analysis/20210420_onny_zjg_analysis.Rmd* This script does most of the analyses in the paper. Inputs are available on [Dryad](https://doi.org/10.5068/D1H963)
    2. *analysis/20210415_onny_visual_eDNA_comparison.Rmd* This script compares visual survey data to eDNA results.
    3. *Anacapa_commands/* Includes the scripts used to run the Anacapa Toolkit
    4. *decontamination/20210112_onny_generalized_decontam_script.R* script to run decontamination
2. **Data**
    1. *elas_bio_12S_taxonomy_tables* Anacapa Output of MiFish 12S Elasmobranch data
    2. *miu_bio_12S_taxonomy_tables* Anacapa Output of MiFish 12S Elasmobranch data
    3. *20201228_raja_ampat_data_visual_census_FIX_zjg.txt* Visual census data
    4. *global_database_taxonomy_20201204.csv* CRUX taxonomy database

The CRUX-12S reference database was generated using standard CRUX parameters using EMBL and NCBI GenBank sequences downloaded in October 2019. Additional sequences were added as described in the manuscript.
