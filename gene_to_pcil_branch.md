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
# 4. Prepare the PCV input

Get the position of the selected PCV:

```r
#creating selection 
input_pcil<- "Sobic.003G260300"

```
---

# 5. Select PCIL (+)

Search for PCIL (+) within families expected to segregate for the PCV:

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
# pcil_postive, has all of the lines that are segregating among the families hypothesized to be segregating for your PCV
head(pcil_positives$pcil_positive)

# you can check how many you have by
nrow(pcil_positives$pcil_positive)
```
<img width="1021" height="161" alt="image" src="https://github.com/user-attachments/assets/586e10fd-52b6-4a04-bc9b-be44d98b3d4b" />


Check family representation:

```r
# if you want to understand what families are represented, you can use the metadata for it
left_join(pcil_positives$pcil_positive,pcil_data$metadata, by = "SampleID") %>% 
  count(Region, Family)
```
<img width="522" height="775" alt="image" src="https://github.com/user-attachments/assets/32a0f05d-6f5f-4875-a766-20d45201db07" />

Inspect the selected PCIL (+):

```r
# best_lines, are the best lines recommended according to you "sel", in this case the top 5 lines for each region, if available.
# the criteria for selection can be found here:
head(pcil_positives$best_lines)
```
<img width="1087" height="120" alt="image" src="https://github.com/user-attachments/assets/2f561ca8-297d-492d-a2c1-c9b0bb304c73" />


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
<img width="3110" height="1952" alt="image" src="https://github.com/user-attachments/assets/becc9152-3147-4fa5-98b9-59530efd571c" />


### Selected PCIL (+)

```r
# If you ran the "sel" option, you can use the plot_pcil_pos_best function to view them
source("https://gist.githubusercontent.com/claracruet/7613608cf517186d49bd234668d61e02/raw/plot_best_pcil_positive.R")

# running the plotting function
best_pcil_pos_plot<- plot_best_pcil_positive(pcil_data = pcil_data, pcil_pos_pcv = pcil_positives)
best_pcil_pos_plot
```
<img width="3110" height="1952" alt="image" src="https://github.com/user-attachments/assets/5c855949-ea18-40f3-bcf4-72159aebb5ab" />


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
                                           # for your PCV
                                           regions =   pcil_positives$regions
) 

```
<img width="540" height="430" alt="image" src="https://github.com/user-attachments/assets/382ae74c-811d-45e3-9e74-eef862ae3238" />

Initial candidates: Total number of PCIL (-)
Using subset PCIL (+), best PCIL (+) from the previous step

<img width="526" height="74" alt="image" src="https://github.com/user-attachments/assets/dbab2aea-bd0f-4d57-abdc-f849be084c3e" />

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
<img width="1091" height="240" alt="image" src="https://github.com/user-attachments/assets/ffa80b94-3336-499e-8d82-4d8f51b2c3a5" />


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
<img width="543" height="66" alt="image" src="https://github.com/user-attachments/assets/92eb1578-4f51-4981-b3c2-e9f72a1ac6a7" />

You will obtain one plot per each variant and per each PCIL (+).

```
# Plotting my number one ranked PCIL (+)
plot_pcil_pairs_negatives$Sobic.003G260300_GMS_MN2025_127047
```
<img width="3110" height="1952" alt="image" src="https://github.com/user-attachments/assets/6eac2061-c87f-497e-972a-4f9a8a14c136" />


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
