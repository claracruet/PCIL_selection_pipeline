# PCIL Selection Guide — Phenotype filter example

This example shows how to restrict PCIL (+) and PCIL (-) selection to lines meeting a **phenotypic criterion** using `global_available_ids`.

The genomic target can still be a gene, position, or region. In this example, a gene is used as the target.

```text
Phenotype data
 ↓
Filter lines by phenotype
 ↓
Map phenotype IDs to genomic SampleIDs
 ↓
Use as global_available_ids
 ↓
Select PCIL (+)
 ↓
Select PCIL (-) using the same filter
```

# 1.  Load PCIL genomic data

```r
# we are going to first source the genomic infomation for the PCILs that we will be using
source("https://gist.githubusercontent.com/claracruet/b6ade06ffa38c1e6bb97c813621632ea/raw/load_pcil_data.R")

# loading up the PCIL data
# Run the loading function and save the PCIL data list into an object.
pcil_data <- load_pcil_data()

```
---
# 2. Extracting sampleIDs for a particular phenotype

```r
# reading pheno
phenotype <- read.csv("https://raw.githubusercontent.com/claracruet/PCIL_selection_pipeline/refs/heads/main/phenotype_example.csv")

# Looking at the range for Days to flowering time
range(phenotype$DTF_BLUE)

# Keeping lines that flower between 60 and 70 days
pheno_for_pcil <- phenotype[phenotype$DTF_BLUE %in% 60:70,]
head(pheno_for_pcil)
```
<img width="421" height="127" alt="image" src="https://github.com/user-attachments/assets/6a0bb7fa-2d68-41ba-96d7-aba3a68c5e04" />


```r
# In this case I have my phenotype by BC1F5 seed ID, so I need to use the seed metadata to recover my sampleIDs
# NOTE: If you have questions on using the seed metadata please contact Clara Cruet-Burgos
seed_metadata<-read.csv("https://raw.githubusercontent.com/claracruet/PCIL_selection_pipeline/refs/heads/main/PCIL_genotype_to_seed_metadata_CLEAN_and_SEED%20(1).csv")

# getting samples for my seed list
seed_to_sample_list<- seed_metadata[seed_metadata$BC1F5_SEED_ID %in% pheno_for_pcil$SEED.ID,]
seed_to_sample_list<-seed_to_sample_list[c("BC1F5_SEED_ID", "sample_id")]

# lets do a join
pheno_to_pcil_sample<- left_join(pheno_for_pcil, seed_to_sample_list, by=c("SEED.ID"="BC1F5_SEED_ID"))

head(pheno_to_pcil_sample)
```
<img width="557" height="126" alt="image" src="https://github.com/user-attachments/assets/8007c170-3f3f-4a34-8f72-847b84483aa6" />

```r
 # extract unique genomic SampleIDs associated with the phenotype filter
sample_ids <- unique(na.omit(pheno_to_pcil_sample$sample_id))
```
---

# 3. Load PCIL (+) selection function

```r
# loading the pcil_pos function
# Load the function that identifies PCIL-positive lines.
source("https://gist.githubusercontent.com/claracruet/189e3a4a2aabf0527ef0845832597439/raw/select_pcil_positive.R")
```

---
# 4. Prepare the gene input

```r
#creating selection 
input_pcil<- "Sobic.003G260300"

```
---

# 5. Select PCIL (+)

Identify PCIL (+) lines carrying an introgression that fully spans the target gene, restricting the search to lines that meet the selected phenotypic criterion.

Here we will use the `global_available_ids` to restrict the search to lines meeting the phenotype criterion.
```r
# running pcil_pos
pcil_positives<- select_pcil_positive(pcil_data = pcil_data, 
                                      input = input_pcil, 
                                      type = "gene",  
                                      sel = 5, 
                                      global_available_ids = sample_ids)

# select_pcil_positive, returns a list
names(pcil_positives)
```
<img width="313" height="32" alt="image" src="https://github.com/user-attachments/assets/1a6ac3a1-0c00-41e1-ad30-4618bd061f76" />

Inspect all PCIL (+):

```r
# pcil_positive contains all PCIL (+) lines with an introgression fully spanning the target gene within the
# phenotype-filtered population
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
# best_lines contains the preferred PCIL (+) lines selected
# according to sel; here, up to 5 lines are returned for the target:
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


# 7. Select PCIL (-)

Select PCIL (-) controls for the final PCIL (+) lines. 

Important: Population filters are not automatically carried from `select_pcil_positive()` into `select_pcil_negative()`. Apply the same `global_available_ids` restriction again during PCIL (-) selection.

```r
# loading the PCIL (-)
source("https://gist.githubusercontent.com/claracruet/3f758a2f7d74a7d2f8278309b9500f67/raw/select_pcil_negative.R")

# running pcil_neg
pcil_negatives<- select_pcil_negative(pcil_data = pcil_data, 
                                      n_neg = 3, 
                                      pcil_positive_df = pcil_positives$best_lines, 
                                      regions =   pcil_positives$regions, 
                                      global_available_ids = sample_ids)
                                       

```
<img width="554" height="428" alt="image" src="https://github.com/user-attachments/assets/493a1946-113e-44ee-b449-0a11249e4700" />


Initial candidates: Total number of PCIL (-)
Using subset PCIL (+), best PCIL (+) from the previous step
<img width="525" height="78" alt="image" src="https://github.com/user-attachments/assets/64928a2c-022a-42cd-8e07-ed502f629578" />


The console output summarizes the PCIL (-) matching process for each PCIL (+). It reports the candidate population, whether same-family candidates are available, the closest candidates based on IBS distance, and the recommended PCIL (-) match. Lower IBS distance indicates greater genome-wide genetic similarity to the focal PCIL (+).



Best PCIL (-) match:

```r
# select_pcil_negative() returns a list
names(pcil_negatives)
```
<img width="287" height="39" alt="image" src="https://github.com/user-attachments/assets/baff55e0-c904-4242-a987-791af3158dd2" />


```r
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

# 8. Visualize PCIL (+) / PCIL (-) pairs

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


You will obtain one plot for each target gene × PCIL (+) combination.

```r
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
