# PCIL Selection Guide — Seed availability filter example

This guide shows how to move from a **gene interest** to PCIL (+) and PCIL (-) lines. This pipeline will select all lines that have an introgression covering a gene despite the prescence or abscence of specific PCVs.



# 1.  Load PCIL genomic data

```r
# we are going to first source the genomic infomation for the PCILs that we will be using
source("https://gist.githubusercontent.com/claracruet/b6ade06ffa38c1e6bb97c813621632ea/raw/load_pcil_data.R")

# loading up the PCIL data
# Run the loading function and save the PCIL data list into an object.
pcil_data <- load_pcil_data()

```
---
# 2. Extracting sampleIDs based on seed availability

Google Drive file ID for the CAL seed inventory. The inventory contains the current seed amount associated with each seed LOT.UID.
```r
# Reading inventory
file_id <- "1Yqdb-rlUgAA5pLR631rt-AAPEGvm035z"

# Create a direct-download URL from the Google Drive file ID.
url <- paste0(
  "https://drive.google.com/uc?export=download&id=",
  file_id
)

# Read the current CAL seed inventory.
cal_inventory <- read.csv(url)

# Inspect the inventory.
head(cal_inventory)
dim(cal_inventory)
```
<img width="1251" height="273" alt="image" src="https://github.com/user-attachments/assets/7faccb1e-8340-4bcf-bf1b-4c3f969750b4" />

The CAL inventory is organized by seed LOT.UID, while the genomic PCIL datasets are organized by sample_id.
```r
# This metadata connects PCIL genomic sample IDs to their corresponding seed lots across generations.
seed_metadata <- read.csv(
  "https://raw.githubusercontent.com/claracruet/PCIL_selection_pipeline/refs/heads/main/PCIL_genotype_to_seed_metadata_CLEAN_and_SEED%20(1).csv"
)

# Inspect the seed metadata.
head(seed_metadata)

```
<img width="1280" height="349" alt="image" src="https://github.com/user-attachments/assets/e4746910-2ea3-485c-9e20-8b221cccf4b9" />

Here, CURRENT_LOT.UID is used to identify the seed lot representing the current generation of each PCIL line.
If seed availability from a specific generation is needed,
- CURRENT_LOT.UID can be replaced with the corresponding
- generation-specific LOT.UID column (e.g., BC1F3_LOT.UID , BC1F4_LOT.UID ,BC1F5_LOT.UID,  etc.).

```r
#  Filtering inventory to PCILs
pcil_inventory <- cal_inventory[cal_inventory$LOT.UID %in% seed_metadata$CURRENT_LOT.UID,]

# Inspect the PCIL-specific inventory.
head(pcil_inventory)
dim(pcil_inventory)
```
<img width="1271" height="272" alt="image" src="https://github.com/user-attachments/assets/e8da84c9-9f4e-4c39-afa7-b111e910bd56" />

Seed inventory is recorded as seed weight (g), so seed number is estimated using a conservative thousand-kernel weight (TKW)of approximately 35 g:

``` text
35 g / 1000 seeds ≈ 0.035 g/seed
1000 seeds / 35 g ≈ 28.6 seeds/g

Therefore:
estimated seeds = seed weight (g) × 28.6

An additional 100 seeds are reserved as breeding stock and
removed from the estimated number available for experiments.
usable seeds = (seed weight × 28.6) - 100
```
``` r
# calculating usable seeds
pcil_inventory$est_seed_number <-  (pcil_inventory$AMOUNT * 28.6) - 100

# establishing the seed amount needed
seed_n=300

# Retain only PCIL lots with enough estimated usable seed
# for the experiment.
pcil_inventory_enough <- pcil_inventory[pcil_inventory$est_seed_number >= seed_n,]
```
Identify the genomic PCIL samples associated with the seed lots that passed the seed-availability filter.
```r
# Filtering by lot.uid with enough seed
seed_metadata_enough <- seed_metadata[seed_metadata$CURRENT_LOT.UID %in% pcil_inventory_enough$LOT.UID,]

# Keep the genomic sample ID and seed lot ID needed for mapping.
seed_metadata_enough_samples <- seed_metadata_enough[c("sample_id", "CURRENT_LOT.UID")]
```

 The resulting table connects each qualifying genomic PCIL line to its seed ID, pedigree, generation, and estimated number of usable seeds.
```r
# Combine genomic sample IDs with seed inventory information.
inventory_and_sampleid <- left_join(
  seed_metadata_enough_samples,
  pcil_inventory_enough[ c("SEED_ID", "LOT.UID", "GID", "PEDIGREE", "GENERATION", "est_seed_number" )],
  by = c("CURRENT_LOT.UID" = "LOT.UID")
)

# Inspect the final filtered PCIL inventory.
head(inventory_and_sampleid)

# now i will extract my sample_ids to use
sample_ids<- inventory_and_sampleid$sample_id
```
<img width="814" height="134" alt="image" src="https://github.com/user-attachments/assets/7a3ce8ee-1d55-4e05-b38d-c430036b26df" />

---

# 3. Load PCIL (+) selection fucntion

```r
# loading the pcil_pos function
# Load the function that identifies PCIL-positive lines.
source("https://gist.githubusercontent.com/claracruet/189e3a4a2aabf0527ef0845832597439/raw/select_pcil_positive.R")
```

---
# 4. Prepare the region input


```r
#creating selection 
input_pcil<- "Sobic.003G260300"

```
---

# 5. Select PCIL (+)

Search for PCIL (+) with introgression for your region

```r
# running pcil_pos
# running pcil_pos
pcil_positives<- select_pcil_positive(pcil_data = pcil_data, 
                                      input = input_pcil, 
                                      type = "gene",  
                                      sel = 5)

# select_pcil_positive, returns a list
names(pcil_positives)
```
<img width="313" height="32" alt="image" src="https://github.com/user-attachments/assets/1a6ac3a1-0c00-41e1-ad30-4618bd061f76" />

Inspect all PCIL (+):

```r
# pcil_postive, has all of the lines have an introgression in your region
head(pcil_positives$pcil_positive)

# you can check how many you have by
nrow(pcil_positives$pcil_positive)
```
<img width="1023" height="194" alt="image" src="https://github.com/user-attachments/assets/c0e54567-8d99-4169-ad66-5f1747e03b3f" />


Check family representation:

```r
# if you want to understand what families are represented, you can use the metadata for it
left_join(pcil_positives$pcil_positive,pcil_data$metadata, by = "SampleID") %>% 
  count(Region, Family)
```
<img width="555" height="447" alt="image" src="https://github.com/user-attachments/assets/97f929ae-af9a-4c31-a190-1e0e4d8ba880" />

Inspect the selected PCIL (+):

```r
# best_lines, are the best lines recommended according to you "sel", in this case the top 5 lines for each region, if available.
# the criteria for selection can be found here:
head(pcil_positives$best_lines)
```
<img width="1075" height="117" alt="image" src="https://github.com/user-attachments/assets/5798e11e-f455-47a2-a9af-1bb084608f4a" />


[PCIL (+) selection documentation](https://github.com/claracruet/PCIL_selection_pipeline/blob/main/select_pcil_positive_ReadMe.md)

---

# 6. Visualize PCIL (+)

### All PCIL (+)

```r
# sourcing function to plot all of the PCIL (+)
source("https://gist.githubusercontent.com/claracruet/d4b8c7dfb22d50c5e31f9a1c3f1ffb94/raw/plot_all_pcil_positive.R")

# running the plotting function
all_pcil_pos_plot<-plot_all_pcil_positive(pcil_pos_pcv = pcil_positives)
all_pcil_pos_plot
```
<img width="3110" height="1952" alt="image" src="https://github.com/user-attachments/assets/190d9a68-9002-4200-9e57-c06ccd5c982b" />


### Selected PCIL (+)

```r
# If you ran the "sel" option, you can use the plot_pcil_pos_best function to view them
source("https://gist.githubusercontent.com/claracruet/7613608cf517186d49bd234668d61e02/raw/plot_best_pcil_positive.R")

# running the plotting function
best_pcil_pos_plot<- plot_best_pcil_positive(pcil_data = pcil_data, pcil_pos_pcv = pcil_positives)
best_pcil_pos_plot
```
<img width="3110" height="1952" alt="image" src="https://github.com/user-attachments/assets/58076b76-1d25-4e3a-9951-f750733b4654" />


[PCIL plotting documentation](https://github.com/claracruet/PCIL_selection_pipeline/blob/main/plotting_functions_ReadMe.md)

---


# 8. Select PCIL (-)

Select PCIL (-) controls for the final PCIL (+) lines. IMPORTANT, you must also filter on the PCIL (-), the filtering options are not passed

```r
# loading the PCIL (-)
source("https://gist.githubusercontent.com/claracruet/3f758a2f7d74a7d2f8278309b9500f67/raw/select_pcil_negative.R")

# running pcil_neg
pcil_negatives<- select_pcil_negative(pcil_data = pcil_data, 
                                      n_neg = 3, 
                                      pcil_positive_df = pcil_positives$best_lines, 
                                      regions =   pcil_positives$regions, 
                                      global_available_ids = sample_ids)
                                      
) 

```
<img width="554" height="428" alt="image" src="https://github.com/user-attachments/assets/493a1946-113e-44ee-b449-0a11249e4700" />


Initial candidates: Total number of PCIL (-)
Using subset PCIL (+), best PCIL (+) from the previous step
<img width="525" height="78" alt="image" src="https://github.com/user-attachments/assets/64928a2c-022a-42cd-8e07-ed502f629578" />


This shows you the selection process for each PCIL (+), it shows the PCIL (+) , then the number of candidates within the same family, then the ones that have the closest IBS, then it shows you the reccomended and it's IBS distance. The lowest the number the more similar to the PCIL (+).


Best PCIL (-) match:

```r
# select_pcil_positive, returns a list
names(pcil_negatives)
```
<img width="287" height="39" alt="image" src="https://github.com/user-attachments/assets/baff55e0-c904-4242-a987-791af3158dd2" />


```
# we have now a list of two
# 'pairs_best', provides the number one PCIL (-) for each PCIL positive.
head(pcil_negatives$pairs_best)
```
<img width="1090" height="113" alt="image" src="https://github.com/user-attachments/assets/fd9dd21e-1819-4b89-b5bc-b2e83e9d8610" />


Additional ranked PCIL (-) candidates:

```r
# we have now a list of two
# 'pairs_extended', provides the top ranked PCIL (-) for the 'n_neg' you provided.
head(pcil_negatives$pairs_extended)
```
<img width="1090" height="238" alt="image" src="https://github.com/user-attachments/assets/77e09f60-96e8-4327-8a76-3d9869def2d9" />


[PCIL (-) selection documentation](https://github.com/claracruet/PCIL_selection_pipeline/blob/main/select_pcil_negative_ReadMe.md)

---

# 9. Visualize PCIL (+) / PCIL (-) pairs

```r
# Sourcing function to plot pcil pairs
source("https://gist.githubusercontent.com/claracruet/88850837f726cdc1a797993e93261847/raw/plot_pcil_pairs.R")

# loading pcil pair plotting function
plot_pcil_pairs_negatives<- plot_pcil_pairs(pcil_neg_sel = pcil_negatives$pairs_extended,  # you must indicate if you want 'pair_best' or 'pairs_extended'
                                            pcil_data = pcil_data, pcil_pos = pcil_positives
)
names(plot_pcil_pairs_negatives)

```
<img width="1025" height="41" alt="image" src="https://github.com/user-attachments/assets/34d8ae9d-8bd8-4fff-b427-617036c5c949" />


You will obtain one plot per each variant and per each PCIL (+).

```
# Plotting my number one ranked PCIL (+)
plot_pcil_pairs_negatives$Sobic.003G260300_25ALM_BC1F3s1_0948
```
<img width="3110" height="1952" alt="image" src="https://github.com/user-attachments/assets/d2921b23-3cd4-4eea-b631-4e8e07b7216e" />


Use:

```r
pcil_negatives$pairs_best
```

instead if only the best PCIL (-) pair should be plotted.

---
