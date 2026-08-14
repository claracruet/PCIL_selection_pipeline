# PCIL Data

The PCIL selection pipeline relies on a standardized set of genomic and metadata resources that describe the structure of the **Pangenome Characterized Introgression Line (PCIL)** population.

These datasets are assembled using:

```r
pcil_data <- load_pcil_data()
```

`load_pcil_data()` collects the core PCIL datasets into a single R object that can be used consistently throughout the PCIL (+) and PCIL (−) selection workflow.

The function does not generate the underlying genomic analyses. Instead, it loads and harmonizes outputs generated from the PCIL DArTseq genotyping data, the PCIL introgression-mapping pipeline, PCIL metadata, and genome annotation resources.
The introgression mapping was conducted by Dr. Gaia Cortinovis as part of the BMG Allele Mining project and is currently being developed for publication as a manuscript.

The resulting object contains:

```text
pcil_data
│
├── introgressions
├── genomewide_introgressions
├── metadata
├── inbreeding_coefficient
├── IBS_dis
└── gene_regions
```

---

## `introgressions`

Contains the donor-derived introgression blocks identified across the PCIL population.

Each row represents one introgression block detected in an individual PCIL line. The main columns used by the PCIL selection and visualization functions are:

| Column            | Description                                                                                                                                                |
| ----------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `SampleID`        | Unique identifier for the PCIL line carrying the introgression block. This identifier is used to connect the line across the different PCIL datasets.      |
| `ChrLabel`        | Chromosome on which the introgression block was detected, represented as `Chr01`–`Chr10`.                                                                  |
| `block_start_bp`  | Genomic start coordinate of the introgression block in base pairs.                                                                                         |
| `block_end_bp`    | Genomic end coordinate of the introgression block in base pairs.                                                                                           |
| `block_len_Mb`    | Length of the introgression block in megabases (Mb).                                                                                                       |
| `mean_donor_frac` | Mean donor-allele fraction across the windows composing the introgression block. Higher values indicate stronger evidence that the block is donor-derived. |
| `startMb`         | Genomic start coordinate of the introgression block in Mb. Used primarily for genome-wide visualization.                                                   |
| `endMb`           | Genomic end coordinate of the introgression block in Mb. Used primarily for genome-wide visualization.                                                     |
| `Clan`            | PCIL clan/recurrent-parent background associated with the line. Used primarily for labeling and visualization.                                             |
| `Family`          | PCIL family associated with the line. Used primarily for labeling and visualization.                                                                       |

These blocks were reconstructed from DArTseq genotypes using the PCIL breeding design and parental genotype information.

The introgression data are used to determine whether a PCIL line contains a donor segment that fully spans a gene, variant, or genomic region of interest. Introgression size and donor fraction are also used during PCIL (+) selection.

**Used by:**

* `select_pcil_positive()`
* `select_pcil_negative()`
* PCIL plotting functions

[Detailed introgression-mapping pipeline](https://github.com/CropAdaptationLab/Lab-Notebooks/tree/62f4389149293c69dd295452915ea30990e2bcc4/members/Gaia/Introgression_Mapping_PCILs)

---

## `genomewide_introgressions`

Provides genome-wide introgression summaries for each PCIL line.

Unlike `introgressions`, where each row represents an individual introgression block, this dataset summarizes the overall donor-introgression burden of each line.

| Column         | Description                                                                                 |
| -------------- | ------------------------------------------------------------------------------------------- |
| `SampleID`     | Unique identifier for the PCIL line.                                                        |
| `total_blocks` | Total number of donor-derived introgression blocks detected across the genome.              |
| `total_Mb`     | Total genomic length, in Mb, represented by donor-derived introgressions across the genome. |

These values provide a measure of **recurrent-parent genome recovery**.

A lower `total_Mb` indicates that less donor-derived genome remains across the PCIL line, while a lower `total_blocks` indicates that the remaining donor genome is distributed across fewer separate introgression segments.

Lines with lower `total_Mb` and fewer `total_blocks` are generally preferred during PCIL selection.

**Used by:**

* PCIL (+) ranking
* PCIL (−) ranking

[Detailed introgression-mapping pipeline](https://github.com/CropAdaptationLab/Lab-Notebooks/tree/62f4389149293c69dd295452915ea30990e2bcc4/members/Gaia/Introgression_Mapping_PCILs)

---

## `metadata`

Contains the PCIL breeding structure associated with each genomic SampleID.

The original PCIL metadata contain additional pedigree and breeding information. `load_pcil_data()` retains the information required for downstream PCIL selection and combines the clan and family identifiers into a single family field.

| Column     | Description                                                                                                                           |
| ---------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| `SampleID` | Unique identifier for the PCIL line and the key used to connect metadata with the genomic datasets.                                   |
| `Family`   | Combined PCIL family identifier represented as `clan/family`. This identifies the breeding background associated with each PCIL line. |

Within `load_pcil_data()`, the family identifier is generated as:

```r
Family <- paste(clan, family, sep = "/")
```

This information is particularly important during PCIL (−) selection, where candidate negative controls from the same family as a PCIL (+) are preferred when available.

**Used by:**

* family-level summaries
* PCIL (+) visualization
* PCIL (−) matching

---

## `inbreeding_coefficient`

Contains the genome-wide inbreeding coefficient (`F`) for each PCIL line.

`F` was calculated from the merged PCIL DArTseq genotype dataset using PLINK. The original PLINK output contains several heterozygosity statistics, but `load_pcil_data()` retains the two columns required by the PCIL selection pipeline.

| Column     | Description                                                                                                                                                                     |
| ---------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `SampleID` | Unique identifier for the PCIL line. This corresponds to the `IID` column in the original PLINK output and is renamed to `SampleID` when the data are loaded.                   |
| `F`        | Genome-wide inbreeding coefficient estimated by PLINK. Higher values indicate greater homozygosity relative to expectation and are used to favor more highly inbred PCIL lines. |

Higher `F` values are preferred when other genomic characteristics are comparable.

`F` is used as part of the ranking process for both PCIL (+) and PCIL (−) selection.

**Used by:**

* `select_pcil_positive()`
* `select_pcil_negative()`

[Inbreeding calculation workflow](https://github.com/CropAdaptationLab/Lab-Notebooks/tree/62f4389149293c69dd295452915ea30990e2bcc4/members/Clara/PCIL_Selection_Script/Development_files)

---

## `IBS_dis`

Contains the pairwise genetic distance between PCIL lines.

Distances were calculated from the PCIL DArTseq genotype data as:

```text
1 - IBS
```

Therefore:

```text
lower distance = greater genome-wide genetic similarity
```

The dataset is stored as a square distance matrix with an additional `SampleID` column.

| Column               | Description                                                                                                                                             |
| -------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `SampleID`           | Identifier of the PCIL line represented by each row of the distance matrix.                                                                             |
| `<SampleID>` columns | Each additional column corresponds to another PCIL SampleID. The value at the intersection of two samples is their pairwise `1 - IBS` genetic distance. |

For example, conceptually:

```text
SampleID   PCIL_A   PCIL_B   PCIL_C
PCIL_A       0       0.08     0.15
PCIL_B      0.08      0       0.12
PCIL_C      0.15     0.12      0
```

A value of `0` indicates no IBS distance between a sample and itself, while smaller values between different samples indicate greater genome-wide genetic similarity.

The IBS matrix is primarily used during PCIL (−) selection to identify negative lines that are genetically similar to a selected PCIL (+).

**Used by:**

* `select_pcil_negative()`

[IBS calculation workflow](https://github.com/CropAdaptationLab/Lab-Notebooks/tree/62f4389149293c69dd295452915ea30990e2bcc4/members/Clara/PCIL_Selection_Script/Development_files)

---

## `gene_regions`

Contains genomic coordinates for genes in the BTx623 v5.1 reference annotation.

The full annotation file contains additional gene information, but the PCIL selection pipeline primarily uses the following columns:

| Column   | Description                                                                                                                                             |
| -------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `GeneID` | BTx623 v5.1 gene identifier. The `.v5.1` suffix is removed when the data are loaded so that standard gene IDs can be provided directly to the pipeline. |
| `seqid`  | Chromosome containing the gene, represented using the chromosome identifiers used by the PCIL genomic datasets.                                         |
| `start`  | Start coordinate of the gene in base pairs.                                                                                                             |
| `end`    | End coordinate of the gene in base pairs.                                                                                                               |

Conceptually:

```text
GeneID
  ↓
Chromosome
Start
End
```

When:

```r
select_pcil_positive(
  type = "gene"
)
```

is used, `select_pcil_positive()` searches `gene_regions` for the requested `GeneID` and internally converts the gene into the standardized format:

```text
Region | Chr | Start | End
```

This allows users to provide gene IDs directly rather than manually supplying genomic coordinates.

**Used by:**

* `select_pcil_positive(type = "gene")`

---

For users of the PCIL selection pipeline, these datasets are loaded automatically through `load_pcil_data()`.

The detailed generation and QC procedures for each dataset are documented separately so that the main PCIL selection guide can remain focused on how to use the data rather than how each genomic resource was generated.

[View `load_pcil_data()` source](https://gist.githubusercontent.com/claracruet/b6ade06ffa38c1e6bb97c813621632ea/raw/3604e4db55aef3d90f86b78af4541bd0fb26cb32/load_pcil_data.R)
