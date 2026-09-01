# PCIL Selection Guide — Gene (PCV) Branch 

This guide shows how to move from a **PCV of interest** to PCIL (+) and PCIL (-) lines. If you have not identified a PCV yet and you have a candidate gene, please refer to the [PCV_from_gene_selection.md](https://github.com/claracruet/PCIL_selection_pipeline/blob/main/PCV_from_gene_selection.md)

```text
Target PCV
 ↓
Identify segregating PCIL families
 ↓
Select PCIL (+)
 ↓
Select PCIL (-)
 ↓
Visualize PCIL (+) / PCIL (-) pairs
```

---

# 1. Identify PCIL families segregating for the PCV

```r
# Variant selection
selection=c("INDEL_Chr03_66131272")

# Now we are going to source the "select_pcil_fam_by_variant_pangb.R" from gist
source("https://gist.githubusercontent.com/claracruet/7fb4425da272d985d747eb032550c80f/raw/select_pcil_fam_by_variant_pangb.R")

# Runing select_pcil_fam_by_variant_pangb()
# the only parameter needed is your selection, it will use the pcil metadata to know which libraries from our database are PCIL parents, you do not have to update the link
results<-select_pcil_families_by_variant(selection=selection)

# Results output a list of three elements
names(results)
```
<img width="432" height="39" alt="image" src="https://github.com/user-attachments/assets/21702c6c-1f4e-468c-9fae-31eb752b3d47" />


Check the outputs:

```r
# Parent genotypes
head(results$geno_pi)

# Family summary
head(results$pcil_family_summary)

# PCIL lines from families expected to segregate
head(results$pcil_summary)
```
<img width="407" height="177" alt="image" src="https://github.com/user-attachments/assets/737ca4aa-77e0-4f4e-8c8b-ce5d1d360a42" />
<img width="535" height="100" alt="image" src="https://github.com/user-attachments/assets/74d17cb2-14ce-4917-8421-c9ee975852e1" />
<img width="607" height="126" alt="image" src="https://github.com/user-attachments/assets/ee82f1d4-2d69-4f5c-b989-4b0fb7b43b21" />


For PCV-based selection, `results$pcil_summary` will be used to restrict the downstream search to lines from families expected to segregate for the PCV.

> PCV presence in individual PCIL lines should ultimately be validated with a diagnostic marker such as KASP.

[Family selection documentation](https://github.com/claracruet/PCIL_selection_pipeline/blob/main/select_pcil_families_by_variant_ReadMe.md)

---

# 2. Load PCIL genomic data

```r
# we are going to first source the genomic infomation for the PCILs that we will be using
source("https://gist.githubusercontent.com/claracruet/b6ade06ffa38c1e6bb97c813621632ea/raw/load_pcil_data.R")

# loading up the PCIL data
# Run the loading function and save the PCIL data list into an object.
pcil_data <- load_pcil_data()

```

[PCIL data documentation](https://github.com/claracruet/PCIL_selection_pipeline/blob/main/load_pcil_data_ReadMe.md)

---

# 3. Load PCIL (+) selection function

```r
# loading the pcil_pos function
# Load the function that identifies PCIL-positive lines.
source("https://gist.githubusercontent.com/claracruet/189e3a4a2aabf0527ef0845832597439/raw/select_pcil_positive.R")
```

[PCIL positive documentation](https://github.com/claracruet/PCIL_selection_pipeline/blob/main/select_pcil_positive_ReadMe.md)

---
# 4. Prepare the PCV input

Get the position of the selected PCV:

```r
#creating selection dataframe
input_pcil<- data.frame(Region=selection, Chr="Chr03", pos=66131272)

```
<img width="270" height="45" alt="image" src="https://github.com/user-attachments/assets/6e0f630b-2dc9-4e13-90b7-9229078852b4" />

---

# 5. Select PCIL (+)

Search for PCIL (+) within families expected to segregate for the PCV:

```r
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
```
<img width="313" height="32" alt="image" src="https://github.com/user-attachments/assets/1a6ac3a1-0c00-41e1-ad30-4618bd061f76" />

Inspect all PCIL (+):

```r
# pcil_positive contains PCIL (+) candidates from families
# hypothesized to segregate for the selected PCV
head(pcil_positives$pcil_positive)

# you can check how many you have by
nrow(pcil_positives$pcil_positive)
```
<img width="823" height="271" alt="image" src="https://github.com/user-attachments/assets/23549418-8788-45f9-b5e1-925eef3d80df" />


Check family representation:

```r
# if you want to understand what families are represented, you can use the metadata for it
left_join(pcil_positives$pcil_positive,pcil_data$metadata, by = "SampleID") %>% 
  count(Region, Family)
```
<img width="551" height="442" alt="image" src="https://github.com/user-attachments/assets/58d285ac-3e29-4347-83fd-63ea6808f000" />


Inspect the selected PCIL (+):

```r
# best_lines, are the best lines recommended according to you "sel", in this case the top 5 lines for each region, if available.
# the criteria for selection can be found here:
head(pcil_positives$best_lines)
```
<img width="891" height="205" alt="image" src="https://github.com/user-attachments/assets/0f3f2f9a-666e-4a46-abd1-bbae0321dd21" />

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
<img width="2624" height="1638" alt="image" src="https://github.com/user-attachments/assets/5c739585-b4fc-4dd2-8aa0-5ebf4d8d3467" />

### Selected PCIL (+)

```r
# If you ran the "sel" option, you can use the plot_pcil_pos_best function to view them
source("https://gist.githubusercontent.com/claracruet/7613608cf517186d49bd234668d61e02/raw/plot_best_pcil_positive.R")

# running the plotting function
best_pcil_pos_plot<- plot_best_pcil_positive(pcil_data = pcil_data, pcil_pos_pcv = pcil_positives)
best_pcil_pos_plot
```
<img width="2624" height="1638" alt="image" src="https://github.com/user-attachments/assets/80f8d4b4-98fa-4a48-ab5f-b429fbc280f0" />


[PCIL plotting documentation](https://github.com/claracruet/PCIL_selection_pipeline/blob/main/plotting_functions_ReadMe.md)

---


# 7. Select PCIL (-)

Select PCIL (-) controls for the final PCIL (+) lines.

```r
# loading the PCIL (-)
source("https://gist.githubusercontent.com/claracruet/3f758a2f7d74a7d2f8278309b9500f67/raw/select_pcil_negative.R")

# running pcil_neg
pcil_negatives<- select_pcil_negative(pcil_data = pcil_data, 
                                            n_neg = 3, 
                                            # any other information such as phenotypes for your final selection
                                            pcil_positive_df = pcil_positives$best_lines, 
                                            # for your PCV
                                            regions =   pcil_positives$regions, 
                                            available_ids = results$pcil_summary[c("sample_id","selection")], 
                                            ) 
```
<img width="792" height="429" alt="image" src="https://github.com/user-attachments/assets/fd8a87ab-ae6b-4870-8e4e-51ca5df833f6" />
Initial candidates: Total number of PCIL (-)
Using subset PCIL (+), best PCIL (+) from the previous step

<img width="605" height="81" alt="image" src="https://github.com/user-attachments/assets/66969f89-e0c8-4b4a-b7f4-4eb26ee51bc5" />
The console output summarizes the PCIL (-) selection process for each PCIL (+), including the number of candidates within the same family, the candidates retained based on IBS similarity, and the recommended PCIL (-) match with its IBS distance. Lower IBS distance indicates greater genetic similarity to the PCIL (+).


Best PCIL (-) match:

```r
# select_pcil_negative() returns a list
names(pcil_negatives)
```

<img width="265" height="36" alt="image" src="https://github.com/user-attachments/assets/ac05747c-81db-46ba-a3a5-4e80702f4cb2" />


```
# we have now a list of two
# 'pairs_best', provides the number one PCIL (-) for each PCIL positive.
head(pcil_negatives$pairs_best)
```
<img width="821" height="173" alt="image" src="https://github.com/user-attachments/assets/cc6359c0-9d71-460c-8ba5-97ed99c87ede" />

Additional ranked PCIL (-) candidates:

```r
# we have now a list of two
# 'pairs_extended', provides the top ranked PCIL (-) for the 'n_neg' you provided.
head(pcil_negatives$pairs_extended)
```
<img width="813" height="199" alt="image" src="https://github.com/user-attachments/assets/c8be6ae6-0c9d-401f-9c0e-e54c71d30ddf" />

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
<img width="543" height="57" alt="image" src="https://github.com/user-attachments/assets/cad4eec8-2588-4d5f-b097-275e666682be" />
One plot is generated for each target × PCIL (+) combination.

```
# Plotting my number one ranked PCIL (+)
plot_pcil_pairs_negatives$INDEL_Chr03_66131272_GMS_MN2025_125057
```
<img width="2623" height="1638" alt="image" src="https://github.com/user-attachments/assets/7aa040c8-4847-4c13-b0b7-433a8d84d47f" />

Use:

```r
pcil_negatives$pairs_best
```

instead if only the best PCIL (-) pair should be plotted.

[PCIL plotting documentation](https://github.com/claracruet/PCIL_selection_pipeline/blob/main/plotting_functions_ReadMe.md)

---

# Final outputs

At the end of the workflow, the primary objects are:

```text
results$pcil_summary
        ↓
Families expected to segregate for the PCV

pcil_positives$pcil_positive
        ↓
PCIL (+) candidates

pcil_positives$best_lines
        ↓
Preferred PCIL (+)

pcil_negatives$pairs_best
        ↓
Best PCIL (+) / PCIL (-) pairs

pcil_negatives$pairs_extended
        ↓
Alternative PCIL (-) matches
```

The final PCIL (+)/PCIL (-) pairs can then be moved forward for **PCV validation and experimental testing**.
