# PCIL Selection Guide — Gene (PCV) Branch

This guide shows how to move from a **gene of interest** to PCIL (+) and PCIL (-) lines for a selected **putative causal variant (PCV)**.

```text
Gene
 ↓
Identify candidate variants
 ↓
Select PCV
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

# 1. Identify candidate variants

Start with a gene of interest and retrieve its annotated variants.

```r
pg_ann_region <- pg_query_db(
  table_name = c("annotations"),
  chrom = "Chr03",
  gene_name = "Sobic.003G260300"
)

table(pg_ann_region$impact)
```

Filter variants based on the criteria relevant to your analysis.

For example, to retain `MODERATE` impact variants:

```r
pg_ann_region_mod <- pg_ann_region[
  pg_ann_region$impact %in% c("MODERATE"),
]
```

Retrieve their genotypes:

```r
variant_geno <- pg_query_genotypes(
  variant_ids = pg_ann_region_mod$variant_id
)
```

Create a short variant summary:

```r
pg_variant_summary <- left_join(
  pg_ann_region_mod[
    c("variant_id", "annotation", "impact", "gene_name", "variant_id")
  ],
  variant_geno[
    c(
      "variant_id",
      "chrom",
      "pos",
      "ref",
      "alt",
      "minor_allele",
      "minor_allele_freq"
    )
  ]
)

pg_variant_summary <- unique(pg_variant_summary)

pg_variant_summary
```

Select the PCV(s) you want to pursue:

```r
selection <- c("INDEL_Chr03_66131272")
```

---

# 2. Identify PCIL families segregating for the PCV

```r
results <- select_pcil_families_by_variant(
  selection = selection
)
```

Check the outputs:

```r
# Parent genotypes
head(results$geno_pi)

# Family summary
head(results$pcil_family_summary)

# PCIL lines from families expected to segregate
head(results$pcil_summary)
```

For PCV-based selection, `results$pcil_summary` will be used to restrict the downstream search to lines from families expected to segregate for the PCV.

> PCV presence in individual PCIL lines should ultimately be validated with a diagnostic marker such as KASP.

[Family selection documentation](LINK_TO_FAMILY_README)

---

# 3. Load PCIL genomic data

```r
pcil_data <- load_pcil_data()
```

[PCIL data documentation](LINK_TO_PCIL_DATA_README)

---

# 4. Prepare the PCV input

Get the position of the selected PCV:

```r
variant_geno_sel <- variant_geno[
  variant_geno$variant_id %in% selection,
]

input_pcil <- variant_geno_sel[
  c("variant_id", "chrom", "pos")
]

names(input_pcil) <- c(
  "Region",
  "Chr",
  "pos"
)

input_pcil
```

---

# 5. Select PCIL (+)

Search for PCIL (+) within families expected to segregate for the PCV:

```r
pcil_positives <- select_pcil_positive(
  pcil_data = pcil_data,
  donor_thresh = 0.8,
  sel = 5,
  type = "pos",
  input = input_pcil,
  window = 1000,
  available_ids = results$pcil_summary[
    c("sample_id", "selection")
  ]
)
```

Inspect all PCIL (+):

```r
head(pcil_positives$pcil_positive)

nrow(pcil_positives$pcil_positive)
```

Check family representation:

```r
left_join(
  pcil_positives$pcil_positive,
  pcil_data$metadata,
  by = "SampleID"
) %>%
  count(Region, Family)
```

Inspect the selected PCIL (+):

```r
head(pcil_positives$best_lines)
```

[PCIL (+) selection documentation](LINK_TO_PCIL_POSITIVE_README)

---

# 6. Visualize PCIL (+)

### All PCIL (+)

```r
all_pcil_pos_plot <- plot_all_pcil_positive(
  pcil_pos_pcv = pcil_positives
)

all_pcil_pos_plot
```

### Selected PCIL (+)

```r
best_pcil_pos_plot <- plot_best_pcil_positive(
  pcil_data = pcil_data,
  pcil_pos_pcv = pcil_positives
)

best_pcil_pos_plot
```

[PCIL plotting documentation](LINK_TO_PLOTTING_README)

---

# 7. Optional — Restrict the PCIL population

You can restrict the PCIL search using any external information that can be converted into a vector of `SampleID`s.

For example:

```text
Family
Phenotype
Seed availability
Other experimental requirements
```

Create the SampleID vector:

```r
samples_to_keep <- ...
```

Then rerun PCIL (+) selection using:

```r
pcil_positives_filtered <- select_pcil_positive(
  pcil_data = pcil_data,
  donor_thresh = 0.8,
  sel = 5,
  type = "pos",
  input = input_pcil,
  window = 1000,
  available_ids = results$pcil_summary[
    c("sample_id", "selection")
  ],
  global_available_ids = samples_to_keep
)
```

For example, `samples_to_keep` could contain only lines with sufficient seed for the planned experiment.

> Use the same experimental restrictions when selecting the corresponding PCIL (-).

---

# 8. Select PCIL (-)

Select PCIL (-) controls for the final PCIL (+) lines.

```r
pcil_negatives <- select_pcil_negative(
  pcil_data = pcil_data,
  n_neg = 3,
  pcil_positive_df = pcil_positives_filtered$best_lines,
  regions = pcil_positives_filtered$regions,
  available_ids = results$pcil_summary[
    c("sample_id", "selection")
  ],
  global_available_ids = samples_to_keep
)
```

Best PCIL (-) match:

```r
head(pcil_negatives$pairs_best)
```

Additional ranked PCIL (-) candidates:

```r
head(pcil_negatives$pairs_extended)
```

[PCIL (-) selection documentation](LINK_TO_PCIL_NEGATIVE_README)

---

# 9. Visualize PCIL (+) / PCIL (-) pairs

```r
pcil_pair_plots <- plot_pcil_pairs(
  pcil_neg_sel = pcil_negatives$pairs_extended,
  pcil_data = pcil_data
)

pcil_pair_plots
```

Use:

```r
pcil_negatives$pairs_best
```

instead if only the best PCIL (-) pair should be plotted.

[PCIL plotting documentation](LINK_TO_PLOTTING_README)

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
