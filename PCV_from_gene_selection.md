# PCIL Selection Guide — Gene (PCV) Branch: PCV selection

This guide shows how to move from a **gene of interest** to a selected **putative causal variant (PCV)**.

```text
Gene
 ↓
Identify candidate variants
 ↓
Select PCV
```

---

# 1. Identify candidate variants

Start with a gene of interest and retrieve its annotated variants.

```r
# using pg_query_db to get type of variants that we have in our gene
pg_ann_region<-  panGenomeBreedr::fetch_table_region(
    table_name = c("annotations"),
    chrom = "Chr03",
    gene_name = "Sobic.003G260300", connect_db_mode = 'online')


head(pg_ann_region[1:10,1:10])
```

<img width="999" height="108" alt="image" src="https://github.com/user-attachments/assets/d8e79efa-aab9-4b4c-b8b6-4d73098bcd4a" />

```
table(pg_ann_region$impact)
```
<img width="268" height="83" alt="image" src="https://github.com/user-attachments/assets/5bd5a119-02e7-4de7-a2c7-66eea1376e11" />


Filter variants based on the criteria relevant to your analysis.

For example, to retain `MODERATE` impact variants:

```r
# Filtering to MODERATE
pg_ann_region_mod<- pg_ann_region[pg_ann_region$impact %in% c("MODERATE"),]

head(pg_ann_region_mod[1:10,1:10])
```
<img width="1037" height="116" alt="image" src="https://github.com/user-attachments/assets/50e154ae-ef81-48e4-8ecb-b289e12dad17" />


Retrieve their genotypes:

```r
# Now we are going to extract the genotypic information to get the maf to assist in our decision
# # retrieve genotype information for the retained variants
variant_geno<- panGenomeBreedr::fetch_genotypes_by_id(variant_ids = pg_ann_region_mod$variant_id, connect_db_mode = 'online')

head(variant_geno[1:10,1:10])
```
<img width="754" height="111" alt="image" src="https://github.com/user-attachments/assets/77c2e370-2b2d-4e45-a62b-7dd08730b223" />

Create a short variant summary to assist in the decision (joining annotations and genotype information):

```r
# we are going to create an object with the main information to assist our decision for PCV
# Joining for a short summary
pg_variant_summary<- left_join(pg_ann_region_mod[c("variant_id","annotation","impact","gene_name")], 
                               variant_geno[c("variant_id","chrom","pos","ref","alt", "minor_allele","minor_allele_freq")],
 by=c("variant_id"))

# The variants can be duplicated because annotations happens by transcript, therefore I will run the unique function to keep them only once
pg_variant_summary<- unique(pg_variant_summary)

# looking at our result
head(pg_variant_summary[1:10,1:10])
```
<img width="1742" height="256" alt="image" src="https://github.com/user-attachments/assets/10668097-5827-4d55-ba17-a751e5f3a34e" />

## Select the PCV

Use the annotation, predicted impact, allele frequency, and other available information to select the PCV to pursue.

In this example:

```r
# select the PCV to pursue
selection <- c("INDEL_Chr03_66131272")
```

The selected PCV can now be passed to `select_pcil_families_by_variant()` to identify PCIL families and lines hypothesized to segregate for the variant based on parental genotypes.

[Continue to PCV family and PCIL pair selection](https://github.com/claracruet/PCIL_selection_pipeline/blob/main/pcv_to_pcil_pairs_pipeline_ReadMe.md)
