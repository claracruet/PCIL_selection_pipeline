# PCIL Family Selection by Variant

`select_pcil_families_by_variant()` identifies PCIL families expected to segregate for one or more selected variants.

The function uses parental genotypes together with the PCIL breeding design to determine which donor × recurrent-parent combinations are informative for each variant. The resulting set of PCIL families can then be used to restrict downstream PCIL (+) searches to families in which the selected variant is expected to segregate.

The function is run using:

```r
pcil_families <- select_pcil_families_by_variant(
  selection = selection
)
```

This function is particularly useful when the target of interest is a **specific polymorphism or set of polymorphisms** and the objective is to identify the PCIL families in which those variants are expected to segregate.

---

## Selection logic

PCIL families were generated using different combinations of recurrent and donor parents. Consequently, not every PCIL family is informative for every variant.

For a PCIL family to be considered potentially segregating for a selected variant:

1. The recurrent parent must have a known homozygous genotype at the variant.
2. A potential donor must carry a known homozygous genotype different from the recurrent parent.
3. That donor must have been used to generate a PCIL family within that recurrent-parent background.

Importantly, a family identified by this function is **expected to segregate for the variant based on parental genotypes and pedigree**. This does not mean that every PCIL line within that family carries the variant or the corresponding donor introgression.

Individual PCIL (+) carriers are identified in subsequent steps of the pipeline.

---

## Input

### `selection`

Character vector containing one or more variant IDs to evaluate.

For example:

```r
selection <- c(
  "SNP_Chr03_66131272",
  "INDEL_Chr03_66131500"
)
```

Each variant is evaluated independently. Therefore, different variants can result in different sets of potentially segregating PCIL families.

---

### `metadata_pcil_current`

Location of the current PCIL metadata file.

By default:

```r
metadata_pcil_current =
  "https://github.com/claracruet/File_sharing/releases/download/pcil_metadata/metadata_picls_summer2026.txt"
```

The metadata describe the PCIL breeding structure and provide the recurrent- and donor-parent library identifiers required to connect parental genotypes with individual PCIL families.

The function uses:

- `clan_rs_Lib` to identify recurrent-parent libraries
- `fam_rs_Lib` to identify family/donor-parent libraries

These parental libraries are used to retrieve genotype information for the selected variants.

---

## Parental genotype retrieval

The function identifies all parental libraries represented in the current PCIL population and retrieves their genotypes at the variants provided in `selection`.

Variant genotypes are obtained using:

```r
pg_query_genotypes(
  variant_ids = selection,
  meta_data = c(
    "variant_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "variant_type"
  )
)
```

The raw genotype calls are converted into interpretable genotype classes:

| Genotype | Classification |
| -------- | -------------- |
| `0\|0` | `Reference` |
| `1\|1` | `Alternate` |
| `0\|1` or `1\|0` | `Heterozygous` |
| `./.` | Missing (`NA`) |

The resulting genotype information is then connected to the parental accession metadata.

---

## Identification of potentially segregating families

For each selected variant, the function evaluates each recurrent parent represented in the PCIL population.

The genotype of the recurrent parent is compared against the genotypes of the other PCIL parental accessions.

Potential donors are defined as parental accessions that:

- have a non-missing genotype,
- are not heterozygous, and
- carry a genotype different from the recurrent parent.

For example:

```text
Variant: PCV_1

Recurrent parent     Reference

Potential donor A    Alternate      → informative
Potential donor B    Reference      → same as recurrent parent
Potential donor C    Heterozygous   → excluded
Potential donor D    Missing        → excluded
```

The PCIL metadata are then used to determine whether each potential donor was actually crossed with that recurrent parent.

Only existing PCIL families generated from an informative recurrent-parent × donor combination are retained.

This process is repeated independently for every variant supplied through `selection`.

---

## Output

When potentially segregating PCIL families are identified, the function returns a list containing three objects:

```text
pcil_families
│
├── pcil_family_summary
├── geno_pi
└── pcil_summary
```

---

## `geno_pi`

Contains the parental genotype information used to determine which PCIL families are potentially segregating for each selected variant.

The data are organized by PI identity and retain the library and sample identifiers associated with the parental genotype records.

| Column | Description |
| ------ | ----------- |
| `pinumber` | PI identifier associated with the parental accession. |
| `lib` | Library identifier used to connect the parental accession with the PCIL metadata and genotype data. |
| `sample` | Sample identifier associated with the parental genotype record. |
| `<variant ID>` | Genotype of the parental accession at the selected variant. Values are represented as `Reference`, `Alternate`, `Heterozygous`, or missing. One column is generated for each variant supplied through `selection`. |

For multiple variants, the structure is conceptually:

```text
pinumber   lib   sample   Variant_A   Variant_B   Variant_C
...
```

This table provides the parental genotype information underlying the family-selection process.

---

## `pcil_family_summary`

Provides a summary of the PCIL material predicted to segregate for each selected variant.

Each row represents a **variant × recurrent-parent clan combination** for which at least one potentially segregating PCIL family was identified.

| Column | Description |
| ------ | ----------- |
| `clan` | Recurrent-parent clan associated with the selected PCIL families. |
| `families` | Number of distinct PCIL families within the clan predicted to segregate for the selected variant. |
| `lines` | Total number of PCIL lines belonging to those selected families. |
| `rp_genotype` | Genotype classification of the recurrent parent at the selected variant (`Reference` or `Alternate`). |
| `recurrent_allele` | Actual allele carried by the recurrent parent at the selected variant. This is obtained from the reference or alternate allele associated with the variant according to `rp_genotype`. |
| `selection` | Variant ID for which the family selection was performed. |

This table provides a rapid summary of **how much potentially informative PCIL material exists for each variant and recurrent-parent background**.

---

## `pcil_summary`

Contains the individual PCIL lines belonging to the families predicted to segregate for each selected variant.

Each row represents an individual PCIL line belonging to an informative recurrent-parent × donor family for a particular variant.

| Column | Description |
| ------ | ----------- |
| `clan` | Recurrent-parent clan of the PCIL line. |
| `clan_pi` | PI identifier associated with the recurrent parent. |
| `family` | PCIL family identifier. |
| `family_pi` | PI identifier associated with the donor/family parent. |
| `lib_id` | Library identifier associated with the PCIL line. |
| `sample_id` | Unique sample identifier for the individual PCIL line. |
| `selection` | Variant for which the PCIL family was identified as potentially segregating. |

The `selection` column is particularly important when multiple variants are evaluated because the same PCIL family or line may be associated with more than one selected variant.

This output defines the individual PCIL lines belonging to the families that are expected to segregate for each variant.

---

## Relationship to PCIL (+) selection

`select_pcil_families_by_variant()` and `select_pcil_positive()` answer two different questions.

```text
select_pcil_families_by_variant()
              ↓
Which PCIL families are expected
to segregate for this variant?
              ↓
       Candidate population
              ↓
select_pcil_positive()
              ↓
Which individual PCIL lines contain
a donor introgression spanning
the target region?
```

Family selection therefore provides a **parental-genotype- and pedigree-based restriction of the PCIL population**.

PCIL (+) selection subsequently evaluates the introgression structure of individual lines within that population.

This distinction is important:

> Membership in a potentially segregating family does not by itself classify an individual line as PCIL (+).

---

## Using the output for variant-specific PCIL selection

The `pcil_summary` output can be used to create a variant-specific list of PCIL lines for downstream selection.

For example:

```r
available_ids <- pcil_families$pcil_summary[
  c("sample_id", "selection")
]
```

This preserves the relationship between each individual PCIL line and the variant for which its family is expected to segregate.

The resulting information can be supplied to downstream PCIL (+) selection so that each variant is evaluated within the appropriate set of PCIL families.

---

## Example

```r
# Variants of interest
selection <- c(
  "SNP_Chr03_66131272",
  "INDEL_Chr03_66131500"
)

# Identify PCIL families expected to segregate
pcil_families <- select_pcil_families_by_variant(
  selection = selection
)

# View family-level summary
pcil_families$pcil_family_summary

# View parental genotype information
pcil_families$geno_pi

# View individual PCIL lines from potentially segregating families
pcil_families$pcil_summary
```

[View `select_pcil_families_by_variant()` source](https://gist.githubusercontent.com/claracruet/7fb4425da272d985d747eb032550c80f/raw/2cdfd5661511cd9574d5e01b5e47ec9a853c83b0/select_pcil_fam_by_variant_pangb.R)
