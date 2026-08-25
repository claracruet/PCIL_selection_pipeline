# Fine-mapping PCIL selection

This example shows how `select_pcil_positive()` can be used to identify **PCIL (+) lines with informative introgression boundaries around a target genomic position**.

The same position is evaluated using progressively larger windows. Increasing the `window` progressively requires the PCIL (+) introgression to span a larger interval around the target. Comparing which lines remain PCIL (+) across increasing window sizes provides a simple way to examine differences in introgression boundaries around the focal position.

In this example, the target is **SNP_Chr06_3939483** at **Chr06:3,939,483**, and PCIL (+) candidates are selected using **1 kb, 5 kb, 10 kb, 25 kb, 50 kb, and 100 kb windows**.

---

## Load the PCIL data and functions

First, load the PCIL genomic data and the functions required for PCIL selection.

```r
# source the PCIL data loading function
source("https://gist.githubusercontent.com/claracruet/b6ade06ffa38c1e6bb97c813621632ea/raw/load_pcil_data.R")

# source the PCIL (+) selection function
source("https://gist.githubusercontent.com/claracruet/189e3a4a2aabf0527ef0845832597439/raw/select_pcil_positive.R")

# source the PCIL (-) selection function
source("https://gist.githubusercontent.com/claracruet/3f758a2f7d74a7d2f8278309b9500f67/raw/select_pcil_negative.R")

# load PCIL genomic data
pcil_data <- load_pcil_data()
```

---

## Define the target position

For a position-based analysis, the input requires a target name (`Region`), chromosome (`Chr`), and genomic position (`pos`).

```r
region_sandeep <- data.frame(
  Region = "SNP_Chr06_3939483",
  Chr = "Chr06",
  pos = 3939483
)
```

---

## Identify PCIL (+) at different window sizes

First, run the target without an additional window. In this example, this run is used only to recover the original target coordinates in pcil_pos_no_window$regions. These coordinates are later assigned to the combined fine-mapping object so that the original focal position, rather than one of the expanded windows, is used for visualization.

```r
pcil_pos_no_window <- select_pcil_positive(
  pcil_data = pcil_data,
  input = region_sandeep,
  type = "pos"
)
```

The target is then evaluated using progressively larger windows.

### 1 kb window

```r
pcil_pos_best_available_1kb <- select_pcil_positive(
  pcil_data = pcil_data,
  input = region_sandeep,
  type = "pos",
  window = 1000,
  sel = 3
)
```

### 5 kb window

```r
pcil_pos_best_available_5kb <- select_pcil_positive(
  pcil_data = pcil_data,
  input = region_sandeep,
  type = "pos",
  window = 5000,
  sel = 3
)
```

### 10 kb window

```r
pcil_pos_best_available_10kb <- select_pcil_positive(
  pcil_data = pcil_data,
  input = region_sandeep,
  type = "pos",
  window = 10000,
  sel = 3
)
```

### 25 kb window

```r
pcil_pos_best_available_25kb <- select_pcil_positive(
  pcil_data = pcil_data,
  input = region_sandeep,
  type = "pos",
  window = 25000,
  sel = 3
)
```

### 50 kb window

```r
pcil_pos_best_available_50kb <- select_pcil_positive(
  pcil_data = pcil_data,
  input = region_sandeep,
  type = "pos",
  window = 50000,
  sel = 3
)
```

### 100 kb window

```r
pcil_pos_best_available_100kb <- select_pcil_positive(
  pcil_data = pcil_data,
  input = region_sandeep,
  type = "pos",
  window = 100000,
  sel = 3
)
```

As the window increases, PCIL (+) lines must carry a donor introgression spanning a progressively larger interval around the target. Different `sel` values can be used to retain the desired number of candidates at each window. Therefore, PCIL (+) lines retained only at smaller windows can be particularly useful for fine mapping because their introgression boundaries occur closer to the focal position. Lines retained across larger windows carry introgressions extending farther from the target.

---

## Combine the fine-mapping candidates

The PCIL (+) results from all six windows can be combined into a single object with the same general structure returned by `select_pcil_positive()`.

```r
top_fine_map <- list(
  
  pcil_positive = dplyr::bind_rows(
    pcil_pos_best_available_1kb$pcil_positive,
    pcil_pos_best_available_5kb$pcil_positive,
    pcil_pos_best_available_10kb$pcil_positive,
    pcil_pos_best_available_25kb$pcil_positive,
    pcil_pos_best_available_50kb$pcil_positive,
    pcil_pos_best_available_100kb$pcil_positive
  ),
  
  best_lines = dplyr::bind_rows(
    pcil_pos_best_available_1kb$best_lines,
    pcil_pos_best_available_5kb$best_lines,
    pcil_pos_best_available_10kb$best_lines,
    pcil_pos_best_available_25kb$best_lines,
    pcil_pos_best_available_50kb$best_lines,
    pcil_pos_best_available_100kb$best_lines
  ),
  
  regions = pcil_pos_no_window$regions
)
```

The resulting object contains:

- `pcil_positive` — all PCIL (+) records identified across the different windows.
- `best_lines` — the preferred PCIL (+) candidates selected across the different windows.
- `regions` — the original target position used to visualize the combined fine-mapping candidates.

Because the same PCIL can satisfy more than one window, a line may occur more than once in the combined object.

---

## Visualize the fine-mapping candidates

The combined object can be passed directly to `plot_all_pcil_positive()` to visualize the introgression boundaries of the PCIL (+) candidates around the original target.

```r
# source the PCIL (+) plotting function
source("https://gist.githubusercontent.com/claracruet/d4b8c7dfb22d50c5e31f9a1c3f1ffb94/raw/c556cf16337fe8d0e5692ec78a712389e264771e/plot_all_pcil_positive.R")

# plot the combined fine-mapping PCIL (+)
all_pcil_pos_plot <- plot_all_pcil_positive(
  pcil_pos_pcv = top_fine_map
)

all_pcil_pos_plot
```
<img width="3110" height="1952" alt="image" src="https://github.com/user-attachments/assets/e97bed77-2d4d-432f-8a6e-d868b3982c53" />


The red target marker identifies the genomic position being fine mapped, while the horizontal segments represent the donor introgressions carried by the PCIL (+) lines. Lines with introgression boundaries close to the target are particularly informative because they can help define a smaller donor interval surrounding the focal position. The selected PCIL (+) lines can then be paired with genetically similar PCIL (-) controls using the standard `select_pcil_negative()` workflow.

---

