# PCIL Selection Guide — Family filter example

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
# 2. Extracting sampleIDs for a family

```r
# Extracting samples in family
# family selection
family<-"IRAT204/SC1074"

# extracting lines 
sample_ids<- pcil_data$metadata[pcil_data$metadata$Family==family,]
```
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

Search for PCIL (+) that have introgression for the gene

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
# pcil_postive, has all of the lines that have an introgression in your gene
head(pcil_positives$pcil_positive)

# you can check how many you have by
nrow(pcil_positives$pcil_positive)
```
<img width="991" height="173" alt="image" src="https://github.com/user-attachments/assets/53e0b5ea-09c9-4eb4-a239-4c333bc0a274" />


Check family representation:

```r
# if you want to understand what families are represented, you can use the metadata for it
left_join(pcil_positives$pcil_positive,pcil_data$metadata, by = "SampleID") %>% 
  count(Region, Family)
```
<img width="546" height="60" alt="image" src="https://github.com/user-attachments/assets/71224517-cc48-452e-a71d-bb93593dd179" />

Inspect the selected PCIL (+):

```r
# best_lines, are the best lines recommended according to you "sel", in this case the top 5 lines for each region, if available.
# the criteria for selection can be found here:
head(pcil_positives$best_lines)
```
<img width="1060" height="97" alt="image" src="https://github.com/user-attachments/assets/799c8f36-3987-476e-af1d-d60d7e43885d" />


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
<img width="3110" height="1952" alt="image" src="https://github.com/user-attachments/assets/d735f8c6-503b-4410-b652-642fcc14caee" />


### Selected PCIL (+)

```r
# If you ran the "sel" option, you can use the plot_pcil_pos_best function to view them
source("https://gist.githubusercontent.com/claracruet/7613608cf517186d49bd234668d61e02/raw/plot_best_pcil_positive.R")

# running the plotting function
best_pcil_pos_plot<- plot_best_pcil_positive(pcil_data = pcil_data, pcil_pos_pcv = pcil_positives)
best_pcil_pos_plot
```
<img width="3110" height="1952" alt="image" src="https://github.com/user-attachments/assets/2d3feca9-afa2-4287-89ee-9272a25937fb" />


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
<img width="539" height="354" alt="image" src="https://github.com/user-attachments/assets/b90d4e3f-de11-4de6-abad-bd4502d49dd4" />

Initial candidates: Total number of PCIL (-)
Using subset PCIL (+), best PCIL (+) from the previous step
<img width="539" height="87" alt="image" src="https://github.com/user-attachments/assets/7fa16983-2759-4721-88b0-01c59c6c9823" />


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
<img width="1092" height="112" alt="image" src="https://github.com/user-attachments/assets/af087e3f-86ab-45fc-8970-dabff2cbdd82" />

Additional ranked PCIL (-) candidates:

```r
# we have now a list of two
# 'pairs_extended', provides the top ranked PCIL (-) for the 'n_neg' you provided.
head(pcil_negatives$pairs_extended)
```
<img width="1061" height="100" alt="image" src="https://github.com/user-attachments/assets/7944b0eb-a4ba-44ef-8454-d3e844b4bf8b" />


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
plot_pcil_pairs_negatives$Sobic.003G260300_25ALM_BC1F3s1_0272
```
<img width="3110" height="1952" alt="image" src="https://github.com/user-attachments/assets/ba814ab7-5936-4c16-a305-fb9e58d2cdd5" />


Use:

```r
pcil_negatives$pairs_best
```

instead if only the best PCIL (-) pair should be plotted.

---
