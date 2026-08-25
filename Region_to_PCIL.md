# PCIL Selection Guide — Gene Branch

This guide shows how to move from a **gene interest** to PCIL (+) and PCIL (-) lines. This pipeline will select all lines that have an introgression covering a gene despite the prescence or abscence of specific PCVs.

```text
Target Gene
 ↓
Select PCIL (+)
 ↓
Select PCIL (-)
 ↓
Visualize PCIL (+) / PCIL (-) pairs
```

---

# 1.  Load PCIL genomic data

```r
# we are going to first source the genomic infomation for the PCILs that we will be using
source("https://gist.githubusercontent.com/claracruet/b6ade06ffa38c1e6bb97c813621632ea/raw/load_pcil_data.R")

# loading up the PCIL data
# Run the loading function and save the PCIL data list into an object.
pcil_data <- load_pcil_data()

```

[PCIL data documentation](https://github.com/claracruet/PCIL_selection_pipeline/blob/main/load_pcil_data_ReadMe.md)

---

# 3. Load PCIL (+) selection fucntion

```r
# loading the pcil_pos function
# Load the function that identifies PCIL-positive lines.
source("https://gist.githubusercontent.com/claracruet/189e3a4a2aabf0527ef0845832597439/raw/select_pcil_positive.R")
```

[PCIL positive documentation](https://github.com/claracruet/PCIL_selection_pipeline/blob/main/select_pcil_positive_ReadMe.md)

---
# 4. Prepare the region input

```r
#creating selection 
input_pcil<- data.frame(Region="QTL_Chr03", Chr="Chr03", Start=66940361, End=67050361)

```
---

# 5. Select PCIL (+)

Search for PCIL (+) that have an introgression including our region

```r
# running pcil_pos
pcil_positives<- select_pcil_positive(pcil_data = pcil_data, 
                                      input = input_pcil, 
                                      type = "region",  
                                      sel = 5)

# select_pcil_positive, returns a list
names(pcil_positives)
```
<img width="313" height="32" alt="image" src="https://github.com/user-attachments/assets/1a6ac3a1-0c00-41e1-ad30-4618bd061f76" />

Inspect all PCIL (+):

```r
# pcil_postive, has all of the lines that have an introgression for your region
head(pcil_positives$pcil_positive)

# you can check how many you have by
nrow(pcil_positives$pcil_positive)
```
<img width="950" height="172" alt="image" src="https://github.com/user-attachments/assets/8357a020-38f9-4627-9726-009a97470443" />


Check family representation:

```r
# if you want to understand what families are represented, you can use the metadata for it
left_join(pcil_positives$pcil_positive,pcil_data$metadata, by = "SampleID") %>% 
  count(Region, Family)
```
<img width="562" height="755" alt="image" src="https://github.com/user-attachments/assets/ca795d95-28b0-4a49-ba04-2932ca4e15c6" />

Inspect the selected PCIL (+):

```r
# best_lines, are the best lines recommended according to you "sel", in this case the top 5 lines for each region, if available.
# the criteria for selection can be found here:
head(pcil_positives$best_lines)
```
<img width="1048" height="116" alt="image" src="https://github.com/user-attachments/assets/9eec9aa0-5812-4d9e-97de-e4cacfcfe279" />


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
<img width="3110" height="1952" alt="image" src="https://github.com/user-attachments/assets/ae111cd6-6b9f-419f-aa39-8a1f03c59435" />


### Selected PCIL (+)

```r
# If you ran the "sel" option, you can use the plot_pcil_pos_best function to view them
source("https://gist.githubusercontent.com/claracruet/7613608cf517186d49bd234668d61e02/raw/plot_best_pcil_positive.R")

# running the plotting function
best_pcil_pos_plot<- plot_best_pcil_positive(pcil_data = pcil_data, pcil_pos_pcv = pcil_positives)
best_pcil_pos_plot
```
<img width="3110" height="1952" alt="image" src="https://github.com/user-attachments/assets/74bb1f0d-5cb5-40cf-915b-18938d05febd" />


[PCIL plotting documentation](https://github.com/claracruet/PCIL_selection_pipeline/blob/main/plotting_functions_ReadMe.md)

---


# 8. Select PCIL (-)

Select PCIL (-) controls for the final PCIL (+) lines.

```r
# loading the PCIL (-)
source("https://gist.githubusercontent.com/claracruet/3f758a2f7d74a7d2f8278309b9500f67/raw/select_pcil_negative.R")

# running pcil_neg
pcil_negatives<- select_pcil_negative(pcil_data = pcil_data, 
                                           n_neg = 3, 
                                           # any other information such as phenotypes for your final selection
                                           pcil_positive_df = pcil_positives$best_lines, 
                                           regions =   pcil_positives$regions
) 

```
<img width="618" height="426" alt="image" src="https://github.com/user-attachments/assets/90bbe491-df2a-47d8-aae1-e1d4be32b52a" />

Initial candidates: Total number of PCIL (-)
Using subset PCIL (+), best PCIL (+) from the previous step
<img width="533" height="79" alt="image" src="https://github.com/user-attachments/assets/fd8d83b8-6f79-4299-8943-43862a5fcdcd" />


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
<img width="1036" height="114" alt="image" src="https://github.com/user-attachments/assets/ab671bd1-48e4-44b1-8446-3643bcd1bd95" />


Additional ranked PCIL (-) candidates:

```r
# we have now a list of two
# 'pairs_extended', provides the top ranked PCIL (-) for the 'n_neg' you provided.
head(pcil_negatives$pairs_extended)
```
<img width="1064" height="127" alt="image" src="https://github.com/user-attachments/assets/40229efb-5aa4-445f-9f5d-7fe58cca56b9" />


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
<img width="1046" height="41" alt="image" src="https://github.com/user-attachments/assets/f6bf4ad2-c4da-4641-b8df-6aa2e07fb44f" />

You will obtain one plot per each variant and per each PCIL (+).

```
# Plotting my number one ranked PCIL (+)
plot_pcil_pairs_negatives$Sobic.003G260300_GMS_MN2025_127047
```
<img width="3110" height="1952" alt="image" src="https://github.com/user-attachments/assets/41edf3b8-224c-4c98-98a2-2a047090a9d1" />


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

