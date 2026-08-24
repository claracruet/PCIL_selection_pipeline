library(devtools)
library(dplyr)
#devtools::install_github("awkena/panGenomeBreedr")
library(panGenomeBreedr)
library(ggplot2)
options(timeout = 300)


##########################################################
#PCIL Selection Guide — Gene (PCV) Branch: PCV selection #
##########################################################

# using pg_query_db to get type of variants that we have in our gene
pg_ann_region<-  panGenomeBreedr::fetch_table_region(
    table_name = c("annotations"),
    chrom = "Chr03",
    gene_name = "Sobic.003G260300", connect_db_mode = 'online')

head(pg_ann_region[1:10,1:10])

# Displaying annotation impacts
table(pg_ann_region$impact)

# Filtering to MODERATE
pg_ann_region_mod<- pg_ann_region[pg_ann_region$impact %in% c("MODERATE"),]

# Looking at the HIGH impact
head(pg_ann_region_mod[1:10,1:10])

# Now we are going to extract the genotypic information to get the maf to assist in our decision
# Extracting genotypes at high impact variants, we are usign pg_query_genotypes
variant_geno<- fetch_genotypes_by_id(variant_ids = pg_ann_region_mod$variant_id, connect_db_mode = 'online')

# looking at our genotypes
head(variant_geno[1:10,1:10])

# we are going to create an object with the main information to assist our decision for PCV
# Joining for a short summary
pg_variant_summary<- left_join(pg_ann_region_mod[c("variant_id","annotation","impact","gene_name","variant_id")], 
                               variant_geno[c("variant_id","chrom","pos","ref","alt", "minor_allele","minor_allele_freq")],
                               by=c("variant_id"))

# The variants can be duplicated because annotations happens by transcript, therefore I will run the unique function to keep them only once
pg_variant_summary<- unique(pg_variant_summary)

# looking at our result
head(pg_variant_summary[1:10,1:10])

# From here I am going to pursue variant "INDEL_Chr03_66131272 "


#######################################################################################################
# Part 2: Identifying PCIL parents segregating for our PCV and their corresponding families and lines #
#######################################################################################################

# Variant selection
selection=c("INDEL_Chr03_66131272")

# Now we are going to source the "select_pcil_fam_by_variant_pangb.R" from gist
source("https://gist.githubusercontent.com/claracruet/7fb4425da272d985d747eb032550c80f/raw/select_pcil_fam_by_variant_pangb.R")

# Runing select_pcil_fam_by_variant_pangb()

# the only parameter needed is your selection, it will use the pcil metadata to know which libraries from our database are PCIL parents, you do not have to update the link
results<-select_pcil_families_by_variant(selection=selection)

# Results output a list of three elements
names(results)

# geno_pi, provides the expected genotype for each of the PCIL parents, it includes both recurrent and donors
head(results$geno_pi)

#  pciL_family_summary, provides a summary of all of the lines and the parents by variants
head(results$pcil_family_summary)

# pcil_summary, provides the full of lines by families that were identified to hypothesized to be segregating
head(results$pcil_summary)

## At this point, you have a list of lines that are hypothesized to be segregating for your PCV based on their parental genotypes
## We will use the "sample_id" column in the results$pcil_summary later to select lines with introgressions that are hypothesized to be carrying the PCV
## Is important that we will need to use KASP markers to validate the actual prescence of our PCV in the introgressed regions


#############################################
# Part 3: Identifying PCIL (+) for our PCVs #
#############################################

# we are going to first source the genomic infomation for the PCILs that we will be using
source("https://gist.githubusercontent.com/claracruet/b6ade06ffa38c1e6bb97c813621632ea/raw/load_pcil_data.R")

# loading up the PCIL data
# Run the loading function and save the PCIL data list into an object.
pcil_data <- load_pcil_data()

# loading the pcil_pos function
# Load the function that identifies PCIL-positive lines.
source("https://gist.githubusercontent.com/claracruet/189e3a4a2aabf0527ef0845832597439/raw/select_pcil_positive.R")

#creating selection dataframe
input_pcil<- data.frame(Region=selection, Chr="Chr03", pos=66131272)

# running pcil_pos
pcil_positives<- select_pcil_positive(pcil_data = pcil_data, 
                                      input = input_pcil, 
                                      type = "pos",  
                                      sel = 5,
                                      window = 1000, 
                                      available_ids = results$pcil_summary[c("sample_id","selection")],
                                      )

# select_pcil_positive, returns a list
names(pcil_positives)

# pcil_postive, has all of the lines that are segregating among the families hypothesized to be segregating for your PCV
head(pcil_positives$pcil_positive)

# you can check how many you have by
nrow(pcil_positives$pcil_positive)

# if you want to understand what families are represented, you can use the metadata for it
left_join(pcil_positives$pcil_positive,pcil_data$metadata, by = "SampleID") %>% 
  count(Region, Family)

# best_lines, are the best lines recommended according to you "sel", in this case the top 5 lines for each region, if available.
# the criteria for selection can be found here:
head(pcil_positives$best_lines)

##############################################
# Part 4: Let's plot them to see them better #
##############################################

# sourcing function to plot all of the PCIL (+)
source("https://gist.githubusercontent.com/claracruet/d4b8c7dfb22d50c5e31f9a1c3f1ffb94/raw/plot_all_pcil_positive.R")

# running the plotting function
all_pcil_pos_plot<-plot_all_pcil_positive(pcil_pos_pcv = pcil_positives)
all_pcil_pos_plot

# If you ran the "sel" option, you can use the plot_pcil_pos_best function to view them
source("https://gist.githubusercontent.com/claracruet/7613608cf517186d49bd234668d61e02/raw/plot_best_pcil_positive.R")

# running the plotting function
best_pcil_pos_plot<- plot_best_pcil_positive(pcil_data = pcil_data, pcil_pos_pcv = pcil_positives)
best_pcil_pos_plot


####################################################################################
# Part 5: Selecting PCIL (+) only within a user defined phentoype, condition, etc..#
####################################################################################

# Example 1: You would like to only select PCIL lines for a specific family, maybe you already know the phenotype for the parents

# let's say that you want to select only IRAT204 clans 
# and families for which  donor for which you also have phentoypes
irat_families<-results$pcil_summary[results$pcil_summary$clan=="IRAT204",]
head(unique(irat_families$family))

# lets choose SC1345
# You can gather the samples IDs directly from the results$pcil_summary
samples_to_keep<- results$pcil_summary$sample_id[results$pcil_summary$clan=="IRAT204" & results$pcil_summary$family=="SC1345"]


# Now I'm going to use the global_available_ids for filtering, this is going to be applied first and it will be across all variants
# running pcil_pos
pcil_positives_fam<- select_pcil_positive(pcil_data = pcil_data, donor_thresh = 0.8, sel = 5, type = "pos", 
                                             input = input_pcil, window = 1000, 
                                             available_ids = results$pcil_summary[c("sample_id","selection")], 
                                             global_available_ids = samples_to_keep) # global ids, require specific sample ids to keep

# total number of PCIL now is 5, instead of 113
nrow(pcil_positives_fam$pcil_positive)

# running the plotting function
best_pcil_pos_plot_fam<- plot_best_pcil_positive(pcil_data = pcil_data, pcil_pos_pcv = pcil_positives_fam)
best_pcil_pos_plot_fam


# Example 2: You would like to only select PCIL lines within a specific height range

# loading the phenotype data
pheno_pcils<- read.csv("~/Downloads/BLUES_for_DTF_and_Height_PCILs.csv")


# looking at range
range(pheno_pcils$Height_BLUE)

# getting list of lines with height less than 200
filter_pheno<- pheno_pcils[pheno_pcils$Height_BLUE<=200,]

# lookg at filter pheno
head(filter_pheno)

# Now I'm going to use the global_available_ids for filtering, this is going to be applied first and it will be across all variants
# running pcil_pos
pcil_positives_height<- select_pcil_positive(pcil_data = pcil_data, donor_thresh = 0.8, sel = 5, type = "pos", 
                                      input = input_pcil, window = 1000, 
                                      available_ids = results$pcil_summary[c("sample_id","selection")], 
                                      global_available_ids = filter_pheno$sample_id) # global ids, require specific sample ids to keep

# total number of PCIL now is 41, instead of 113
nrow(pcil_positives_height$pcil_positive)

# running the plotting function
best_pcil_pos_plot_height<- plot_best_pcil_positive(pcil_data = pcil_data, pcil_pos_pcv = pcil_positives_height)
best_pcil_pos_plot_height


# Example 3: You would like to only select PCIL lines that have a specific number of seeds available

# For this we are going to use the CAL inventory. Please note that this inventory is updated weekly and the seed availability might change quickly
cal_inventory<- read.csv("Downloads/Seed_Inventory_for_BMS_Import_dataframe.csv")

# we are also going to use our seed master list, again this will be updated as we advance generations
pcil_seed_masterlist<- read.csv("~/Downloads/PCIL_genotype_to_seed_metadata_CLEAN_and_SEED.csv")

# we are going to use the seed master list to extract current seed weights from our inventory
pcil_seed_current_weight<- cal_inventory[cal_inventory$LOT.UID %in% pcil_seed_masterlist$CURRENT_LOT.UID,]

# as we normally want number of seeds we are going to use a conservative estamate based on TKW of, and we are going to substract the 100 seeds for stock
pcil_seed_current_weight$seed_number<- (pcil_seed_current_weight$AMOUNT*28.6)-100

# we can look at the ranges
range(pcil_seed_current_weight$seed_number)

# now let's say that i want to use only lines that have over 200 seeds because of my experimental design
pcil_enough_seed<- pcil_seed_current_weight[pcil_seed_current_weight$seed_number>=200,]

# now we need to connect back our seed id to our sample ID to be able to filter the genomic data
seed_samples_to_keep<- pcil_seed_masterlist$sample_id[pcil_seed_masterlist$CURRENT_LOT.UID %in% pcil_enough_seed$LOT.UID]

# Now I'm going to use the global_available_ids for filtering, this is going to be applied first and it will be across all variants
# running pcil_pos
pcil_positives_seed<- select_pcil_positive(pcil_data = pcil_data, donor_thresh = 0.8, sel = 5, type = "pos", 
                                             input = input_pcil, window = 1000, 
                                             available_ids = results$pcil_summary[c("sample_id","selection")], 
                                             global_available_ids = seed_samples_to_keep) # global ids, require specific sample ids to keep

# total number of PCIL now is 33, instead of 113
nrow(pcil_positives_seed$pcil_positive)

# running the plotting function
best_pcil_pos_plot_seed<- plot_best_pcil_positive(pcil_data = pcil_data, pcil_pos_pcv = pcil_positives_seed)
best_pcil_pos_plot_seed


######################################################
# Part 6: Selecting PCIL (-) for your PCIL (+) lines #
######################################################

#NOTE YOU SHOULD ALWAYS RESTRICT YOUR PCIL (-) selection by the same criterias you have restricted your PCIL (+)
# You can pass 'available_ids' and 'global_available_ids' the same way you did for the PCIL (+)

# For the PCIL (-), I'm going to continue using our 'pcil_positives_seed' example as it is the most common scenario to filter by

# loading the PCIL (-)
source("https://gist.githubusercontent.com/claracruet/3f758a2f7d74a7d2f8278309b9500f67/raw/select_pcil_negative.R")

# running pcil_neg
pcil_negatives_seed<- select_pcil_negative(pcil_data = pcil_data, 
                                            n_neg = 3, 
                                            # any other information such as phenotypes for your final selection
                                            pcil_positive_df = pcil_positives$best_lines, 
                                            # for your PCV
                                            regions =   pcil_positives$regions, 
                                            available_ids = results$pcil_summary[c("sample_id","selection")], 
                                            ) 


# select_pcil_positive, returns a list
names(pcil_negatives_seed)

# we have now a list of two
# 'pairs_best', provides the number one PCIL (-) for each PCIL positive.
head(pcil_negatives_seed$pairs_best)

# we have now a list of two
# 'pairs_extended', provides the top ranked PCIL (-) for the 'n_neg' you provided.
head(pcil_negatives_seed$pairs_extended)

# Sourcing function to plot pcil pairs
source("https://gist.githubusercontent.com/claracruet/88850837f726cdc1a797993e93261847/raw/569350602db187b7c37a087a00c20f7064f51bde/plot_pcil_pairs.R")

# loading pcil pair plotting function
plot_pcil_pairs_negatives<- plot_pcil_pairs(pcil_neg_sel = pcil_negatives_seed$pairs_extended,  # you must indicate if you want 'pair_best' or 'pairs_extended'
                pcil_data = pcil_data, pcil_pos = pcil_positives
                )
names(plot_pcil_pairs_negatives)

# Plotting my number one ranked PCIL (+)
plot_pcil_pairs_negatives$INDEL_Chr03_66131272_GMS_MN2025_125057
