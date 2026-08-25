## Multiple target selection edits in functions and how to do it!

The pipeline supports the simultaneous evaluation of **multiple targets of the same type**. A single workflow can evaluate multiple PCVs, multiple genes, or multiple genomic regions.

Only **one target type can be used per run**. For example, multiple PCVs can be analyzed together, but PCVs and genes should be analyzed in separate runs.

Each target is tracked independently through the pipeline using its `Region` identifier.

### How multiple targets move through the pipeline

#### `select_pcil_families_by_variant()`

This function is specific to the **PCV branch** and can evaluate multiple PCVs simultaneously.

Provide the PCV IDs as a vector:

```r
selection <- c(
  "PCV_1",
  "PCV_2",
  "PCV_3"
)

results <- select_pcil_families_by_variant(
  selection = selection
)
```

The function evaluates parental genotypes independently for each PCV and identifies the PCIL families and lines hypothesized to segregate for each variant.

The target identity is retained in:

```r
results$pcil_summary$selection
```

This is important because the families available for one PCV may differ from those available for another PCV.

---

#### `select_pcil_positive()`

`select_pcil_positive()` can evaluate multiple targets in the same run.

The input format depends on the selected `type`.

**Multiple PCVs**

Provide one row per PCV:

```r
input_pcil <- data.frame(
  Region = c("PCV_1", "PCV_2", "PCV_3"),
  Chr = c("Chr01", "Chr03", "Chr07"),
  pos = c(1234567, 2345678, 3456789)
)

pcil_positives <- select_pcil_positive(
  pcil_data = pcil_data,
  input = input_pcil,
  type = "pos",
  available_ids = results$pcil_summary[
    c("sample_id", "selection")
  ],
  sel = 5
)
```

For PCVs, `available_ids` can contain a different set of SampleIDs for each PCV because the `selection` column is matched to the corresponding `Region`.

**Multiple genes**

Provide multiple gene IDs:

```r
genes <- c(
  "Sobic.001G000100",
  "Sobic.003G000200",
  "Sobic.007G000300"
)

pcil_positives <- select_pcil_positive(
  pcil_data = pcil_data,
  input = genes,
  type = "gene",
  sel = 5
)
```

**Multiple regions**

Provide one row per genomic region:

```r
regions <- data.frame(
  Region = c("Region_1", "Region_2", "Region_3"),
  Chr = c("Chr01", "Chr03", "Chr07"),
  Start = c(1000000, 2000000, 3000000),
  End = c(1500000, 2500000, 3500000)
)

pcil_positives <- select_pcil_positive(
  pcil_data = pcil_data,
  input = regions,
  type = "region",
  sel = 5
)
```

For all three target types, PCIL (+) identification and ranking are performed **independently for each target**.

The output retains the target identity in:

```r
pcil_positives$pcil_positive$Region
pcil_positives$best_lines$Region
```

Therefore, if:

```r
sel = 5
```

the function attempts to select up to five preferred PCIL (+) lines **for each target**, rather than five lines across the entire analysis.

---

#### `plot_all_pcil_positive()`

No special changes are required for multiple targets.

```r
all_pcil_pos_plot <- plot_all_pcil_positive(
  pcil_pos_pcv = pcil_positives
)
```

The function automatically generates a separate plot for each `Region`.

For example:

```text
all_pcil_pos_plot
│
├── PCV_1
├── PCV_2
└── PCV_3
```

The same behavior applies to multiple genes or multiple genomic regions.

---

#### `plot_best_pcil_positive()`

The selected PCIL (+) lines for multiple targets can also be plotted in a single call:

```r
best_pcil_pos_plot <- plot_best_pcil_positive(
  pcil_data = pcil_data,
  pcil_pos_pcv = pcil_positives
)
```

A separate genome-wide plot is generated for each target.

---

#### `select_pcil_negative()`

Multiple targets do not need to be separated before PCIL (-) selection.

The `regions` and PCIL (+) information can be passed directly from the multi-target `select_pcil_positive()` output:

```r
pcil_negatives <- select_pcil_negative(
  pcil_data = pcil_data,
  pcil_positive_df = pcil_positives$best_lines,
  regions = pcil_positives$regions,
  n_neg = 3
)
```

PCIL (-) matching is performed independently for each:

```text
Target × PCIL (+)
```

combination.

For PCV analyses using target-specific family restrictions, the same `available_ids` can be passed to PCIL (-) selection:

```r
pcil_negatives <- select_pcil_negative(
  pcil_data = pcil_data,
  pcil_positive_df = pcil_positives$best_lines,
  regions = pcil_positives$regions,
  available_ids = results$pcil_summary[
    c("sample_id", "selection")
  ],
  n_neg = 3
)
```

This preserves the target-specific population restriction for each PCV.

---

#### `plot_pcil_pairs()`

PCIL (+)/PCIL (-) pairs for all targets can be visualized in a single call:

```r
pcil_pair_plots <- plot_pcil_pairs(
  pcil_data = pcil_data,
  pcil_neg_sel = pcil_negatives$pairs_extended
)
```

The function generates a separate plot for each:

```text
Target × PCIL (+)
```

combination.

---

### Multi-target workflow

```text
Multiple targets of ONE type
          ↓
PCV | PCV | PCV
or
Gene | Gene | Gene
or
Region | Region | Region
          ↓
select_pcil_positive()
          ↓
PCIL (+) identified and ranked
independently for each target
          ↓
select_pcil_negative()
          ↓
PCIL (-) matched independently
for each Target × PCIL (+)
          ↓
Plotting functions
          ↓
Separate visualizations
for each target
```

