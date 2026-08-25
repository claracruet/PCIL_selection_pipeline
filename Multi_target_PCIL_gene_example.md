# Multi-target PCIL selection

This example shows how to identify **PCIL (+) and PCIL (-) lines that can be used across multiple genomic targets**.

The targets are first evaluated independently to identify the complete set of PCIL (+) candidates. After the targets are evaluated independently, lines that are PCIL (+) for more than one target can optionally be identified and used to restrict the population before ranking preferred multi-target PCIL (+) lines.

The same approach can be applied to PCIL (-) selection by first generating a broad set of candidate controls and then identifying PCIL (-) lines that can serve as controls across multiple targets.

> **Note:** Multiple targets must be of the same `type`. For example, multiple genes can be analyzed together, but genes and regions should be evaluated in separate runs.

Important: Multi-target selection does not require the same PCIL to be used for every target. Each target is first evaluated independently. The shared-line filtering shown here is an additional strategy for prioritizing PCILs that can be informative across multiple targets.

## Load the PCIL data

First, source `load_pcil_data()` and load the genomic information used throughout the PCIL selection pipeline.

```r
# source the function used to load the PCIL genomic data
source("https://gist.githubusercontent.com/claracruet/b6ade06ffa38c1e6bb97c813621632ea/raw/load_pcil_data.R")

# load the PCIL genomic datasets
pcil_data <- load_pcil_data()
```

## Identify PCIL (+) across multiple targets

Source `select_pcil_positive()` and define the targets to evaluate.

In this example, two genes are evaluated simultaneously.

```r
# source the PCIL (+) selection function
source("https://gist.githubusercontent.com/claracruet/189e3a4a2aabf0527ef0845832597439/raw/select_pcil_positive.R")

# define multiple gene targets
input_pcil <- c(
  "Sobic.003G260300",
  "Sobic.007G200100"
)

# identify all PCIL (+) independently for each target
pcil_positives <- select_pcil_positive(
  pcil_data = pcil_data,
  input = input_pcil,
  type = "gene"
)
```

At this stage, each gene is evaluated independently and `pcil_positives$pcil_positive` contains the PCIL (+) lines identified for each target.

### Identify PCIL (+) shared across targets

To identify lines that are PCIL (+) for more than one target, count the number of distinct `Region` values associated with each `SampleID`.

```r
# identify PCIL (+) lines containing more than one target
multi_target_summary_pos <- pcil_positives$pcil_positive %>%
  dplyr::group_by(SampleID) %>%
  dplyr::summarise(
    n_targets = dplyr::n_distinct(Region),
    targets = paste(unique(Region), collapse = ", "),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_targets > 1)

# inspect the multi-target PCIL (+) lines
multi_target_summary_pos
```

<img width="499" height="250" alt="image" src="https://github.com/user-attachments/assets/d29a7ee3-1e42-4d14-b1da-53ee7020e32f" />


`n_targets` indicates the number of distinct targets for which each line was identified as PCIL (+).

### Select the best multi-target PCIL (+)

The PCIL (+) selection can now be repeated while restricting the population to the lines identified above.

`global_available_ids` is used because the same set of SampleIDs is being made available across all targets.

```r
# re-run PCIL (+) selection using only multi-target lines
pcil_positives_multi <- select_pcil_positive(
  pcil_data = pcil_data,
  input = input_pcil,
  type = "gene",
  sel = 5,
  global_available_ids = unique(multi_target_summary_pos$SampleID)
)

names(pcil_positives_multi)
```
<img width="311" height="40" alt="image" src="https://github.com/user-attachments/assets/7c74e2bc-60f3-4015-b3bf-27603c479cfb" />


`pcil_positive` contains all PCIL (+) lines identified after applying the multi-target population restriction.

```r
head(pcil_positives_multi$pcil_positive)

nrow(pcil_positives_multi$pcil_positive)
```

<img width="751" height="285" alt="image" src="https://github.com/user-attachments/assets/6a010187-2a5b-48ac-8dd1-0f57e234fb8c" />


The PCIL metadata can also be used to determine which families are represented for each target.

```r
left_join(pcil_positives_multi$pcil_positive, pcil_data$metadata,by = "SampleID") %>%
  count(Region, Family)
```

<img width="627" height="441" alt="image" src="https://github.com/user-attachments/assets/79174eac-8583-4c13-9d37-b9c6261f4fc1" />


```r
head(pcil_positives_multi$best_lines)
```
<img width="779" height="240" alt="image" src="https://github.com/user-attachments/assets/7aac223e-974a-4c44-b570-09d4f0c4fa6f" />


## Visualize the multi-target PCIL (+)

```r
# source the plotting function
source("https://gist.githubusercontent.com/claracruet/d4b8c7dfb22d50c5e31f9a1c3f1ffb94/raw/plot_all_pcil_positive.R")

# generate PCIL (+) plots
all_pcil_pos_plot <- plot_all_pcil_positive(
  pcil_pos_pcv = pcil_positives_multi
)

all_pcil_pos_plot
```
<img width="1704" height="1049" alt="image" src="https://github.com/user-attachments/assets/baacc32c-d48a-4b90-936d-6754821fc331" />
<img width="1704" height="1049" alt="image" src="https://github.com/user-attachments/assets/fd3bb892-f7a2-482c-a53b-de6f82079b8b" />


If `sel` was used, the preferred PCIL (+) lines can also be visualized genome-wide using `plot_best_pcil_positive()`.

```r
# source the plotting function
source("https://gist.githubusercontent.com/claracruet/7613608cf517186d49bd234668d61e02/raw/plot_best_pcil_positive.R")

# plot the selected PCIL (+)
best_pcil_pos_plot <- plot_best_pcil_positive(
  pcil_data = pcil_data,
  pcil_pos_pcv = pcil_positives_multi
)

best_pcil_pos_plot
```

<img width="1704" height="1049" alt="image" src="https://github.com/user-attachments/assets/09318da2-1bf2-4ccd-9349-d7846275790f" />
<img width="1704" height="1049" alt="image" src="https://github.com/user-attachments/assets/91a01dba-b613-42d2-a81a-d8de157db253" />

In this example, 25ALM_BC1F3s1_0948 is the highest-ranked PCIL (+) for both genes.

# Multi-target PCIL (-) selection

The same general approach can be used to identify **PCIL (-) lines that can serve as controls across multiple targets**.

First, a broad PCIL (-) candidate pool is generated for the selected PCIL (+) lines. PCIL (-) candidates represented across multiple targets are then identified, and the selection is repeated using only those shared candidates.

## Generate a broad PCIL (-) candidate pool

Source `select_pcil_negative()` and initially request a large number of PCIL (-) candidates.

```r
# source the PCIL (-) selection function
source("https://gist.githubusercontent.com/claracruet/3f758a2f7d74a7d2f8278309b9500f67/raw/select_pcil_negative.R")

# generate a broad PCIL (-) candidate pool
pcil_negatives <- select_pcil_negative(
  pcil_data = pcil_data,
  pcil_positive_df = pcil_positives_multi$best_lines,
  regions = pcil_positives_multi$regions,
  n_neg = 100
)
```

Here, `n_neg = 100` is used to generate a broad ranked PCIL (-) candidate pool for each target × PCIL (+) combination. This first run is used to identify negative lines that are candidates for more than one target rather than to make the final PCIL (-) selection.

### Identify PCIL (-) shared across targets

The extended candidate list can be grouped by `SampleID_Negative` to determine how many distinct targets each negative candidate represents.

```r
# identify PCIL (-) candidates represented across multiple targets
multi_target_summary_neg <- pcil_negatives$pairs_extended %>%
  dplyr::group_by(SampleID_Negative) %>%
  dplyr::summarise(
    n_targets = dplyr::n_distinct(Region),
    targets = paste(unique(Region), collapse = ", "),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_targets > 1)

# inspect the multi-target PCIL (-) candidates
multi_target_summary_neg
```

<img width="492" height="240" alt="image" src="https://github.com/user-attachments/assets/e38b51a9-65e0-47e6-94e3-c71d0702a1c1" />


As with the PCIL (+) selection, `n_targets` indicates the number of targets for which the line was identified as a candidate PCIL (-).

## Select the final multi-target PCIL (-)

PCIL (-) selection can now be repeated while restricting the available population to the shared multi-target negative candidates.

```r
# re-run PCIL (-) selection using only multi-target candidates
pcil_negatives_multi <- select_pcil_negative(
  pcil_data = pcil_data,
  pcil_positive_df = pcil_positives_multi$best_lines,
  regions = pcil_positives_multi$regions,
  n_neg = 3,
  global_available_ids = unique(multi_target_summary_neg$SampleID_Negative)
)
```
<img width="641" height="854" alt="image" src="https://github.com/user-attachments/assets/6d895c22-d69d-4d67-9e16-8bbf49891b2e" />

Here, `global_available_ids` restricts the PCIL (-) candidate population to lines that were represented across multiple targets, while `n_neg = 3` returns up to three ranked PCIL (-) candidates for each PCIL (+).

## Inspect the final PCIL (-) selection

```r
head(pcil_negatives_multi$pairs_best)

head(pcil_negatives_multi$pairs_extended)
```

<img width="818" height="489" alt="image" src="https://github.com/user-attachments/assets/2d9923a6-5d79-4a37-936f-3f9a139e06ae" />


## Visualize the multi-target PCIL (+)/PCIL (-) pairs

The final PCIL (+)/PCIL (-) combinations can be visualized using `plot_pcil_pairs()`.

```r
# source the PCIL pair plotting function
source("https://gist.githubusercontent.com/claracruet/88850837f726cdc1a797993e93261847/raw/plot_pcil_pairs.R")

# generate PCIL (+)/PCIL (-) pair plots
plot_pcil_pairs_negatives <- plot_pcil_pairs(
  pcil_neg_sel = pcil_negatives_multi$pairs_extended,
  pcil_data = pcil_data,
  pcil_pos = pcil_positives_multi
)

# inspect the available plots
names(plot_pcil_pairs_negatives)
```

<img width="803" height="89" alt="image" src="https://github.com/user-attachments/assets/91029370-5059-407d-b2c0-7bb843d52dfd" />


Each plot represents a specific **target × PCIL (+)** combination. This makes it possible to inspect the same focal PCIL (+) independently for each target.

In this example, `25ALM_BC1F3s1_0948` is the highest-ranked PCIL (+) for both targets, while `25ALM_BC1F3s1_0504` is the highest-ranked shared PCIL (-) candidate for both genes.

```r
plot_pcil_pairs_negatives$Sobic.003G260300_25ALM_BC1F3s1_0948
plot_pcil_pairs_negatives$Sobic.007G200100_25ALM_BC1F3s1_0948
```

<img width="1704" height="1049" alt="image" src="https://github.com/user-attachments/assets/5910e9c5-07aa-46bb-a7de-162a5e2ac8d8" />
<img width="1704" height="1049" alt="image" src="https://github.com/user-attachments/assets/b32de2bc-11ea-47ef-a221-fae96497aaae" />
