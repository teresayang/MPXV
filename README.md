---
output:
  html_document: default
  pdf_document: default
editor_options: 
  markdown: 
    wrap: 72
---

## Variable Number Tandem Repeats Analysis

*2022-12-05*

## Description

Since May 2022, confirmed cases monkeypox (mpox) have been reported in
many countries. WHO declares the ongoing monkeypox outbreak a Public
Health Emergency of Interntional Concern. The genome of monkeypox
contains a large number of variable number tandem repeats (VNTR). By
analyzing VNTR, researchers will have better understand of the ongoing
monkeypox outbreak. Therefore, package **`MPXV`** provides the functions
not only to estimate the copy number of VNTRs but also to track the
origin of a given MPXV by using VNTRs.

## Installation

To install **`MPXV`** from
[**Github**](https://github.com/teresayang/MPXV_VNTR.git):

```{r Installation from GitHub, eval = FALSE}
if(!require(devtools)) install.packages("devtools")
remotes::install_git("https://github.com/teresayang/MPXV_VNTR.git")
```

To load the installed **`MPXV`** in R:

```{r Load MPXV, eval = FALSE}
library(MPXV)
```

## Output

### VNTR caller:

The function **`VNTRcaller`** estimates the copy number of VNTRs and the
output of the function **`VNTRcaller`** are described below.

For each query sequence of the repetitive unit, there are two output
files provided (.csv file and .fas file). In addition, a summary table
for all query sequence of the repetitive units is exported to a .csv
file.

The table (.csv) includes the following informations: 
1. ID: name of MPXV sequences 
2. r: copy of tandem repeats 
3. match: number of matches
4. mismatch: number of mismatches 
5. indel: number of insertions and
deletions (indels) 
6. score: alignment score of a VNTR region for a query strain 
7. start_pos: starting position of the VNTR region for a query strain

|        ID        |  r  | match | mismatch | indel | score | start_pos |
|:----------------:|:---:|:-----:|:--------:|:-----:|:-----:|:---------:|
|     ON563414     | 28  |  28   |    0     |   0   |  56   |  133095   |
| EPI_ISL_13308158 | 19  |  19   |    0     |   0   |  38   |  132741   |
| EPI_ISL_13889908 | 30  |  30   |    0     |   0   |  60   |  133200   |
| EPI_ISL_14699911 | 32  |  32   |    0     |   0   |  64   |  133102   |
| EPI_ISL_14166709 | 20  |  20   |    0     |   0   |  40   |  133080   |

The VNTR region sequence for query strains are shown in FASTA format
(.fas file).

    >ON563414
    TTTTTTTTTTTTTTTTTTTTTTTTTTTT
    >EPI_ISL_13308158
    TTTTTTTTTTTTTTTTTTT
    >EPI_ISL_13889908
    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    >EPI_ISL_14699911
    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    >EPI_ISL_14166709
    TTTTTTTTTTTTTTTTTTTT

Summary table for all query sequence of the repetitive units are showed
as below. Each row represents a MPXV strain and the first column is the
names of MPXV sequences. Each of the remaining columns represents a
sequence of the repetitive units. Each cell represents a copy of tandem
repeats.

|        ID        | nt133095_T | nt150554_TATGATGGA | nt173267_AT | nt179074_ATATACATT |
|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|
|     ON563414     |     28     |         7          |     24      |         16         |
| EPI_ISL_13308158 |     19     |         7          |     20      |         15         |
| EPI_ISL_13889908 |     30     |         9          |     24      |         7          |
| EPI_ISL_14699911 |     32     |         7          |     25      |         7          |
| EPI_ISL_14166709 |     20     |         7          |     24      |         16         |

### VNTR tracker:

In the VNTR tracker output, each row represents a strain in our
database, and column variables contain accession ID, country name,
collection date, clade, lineage, and three distances (see Yang et al.
(2022)) between a query MPXV strain and a strains in our database. The
rank of strains are ordered by ploymorphism information content distance
(column: "Distance_PIC") (see Yang et al. (2022)) from low to high.

|    Accession     |   Country    | Collection.date | clade | lineage  | Distance_PIC | Distance_L  | Distance_entropy |
|:-------:|:-------:|:-------:|:-------:|:-------:|:-------:|:-------:|:-------:|
|     MN346690     | Cote dIvoire |    2017/3/13    |  IIa  | unassign |      0       |      0      |        0         |
|     MN346692     | Cote dIvoire |    2017/3/5     |  IIa  | unassign | 0.004872201  | 0.007389163 |   0.004750482    |
|     MN346694     | Cote dIvoire |    2017/3/28    |  IIa  | unassign | 0.004872201  | 0.007389163 |   0.004750482    |
|     MN346696     | Cote dIvoire |    2017/4/3     |  IIa  | unassign | 0.019488804  | 0.02955665  |   0.019001927    |
|     KP849470     | Cote dIvoire |      1971       |  IIa  | unassign | 0.030844508  | 0.073219884 |   0.030573661    |
|     MN346702     | Cote dIvoire |    2018/5/11    |  IIa  | unassign | 0.048722011  | 0.073891626 |   0.047504818    |
| EPI_ISL_13544250 |    Canada    |    2022/5/30    |  IIb  |  B.1.4   | 0.048765167  | 0.07320153  |   0.049289862    |
|     MN346698     | Cote dIvoire |    2017/1/24    |  IIa  | unassign | 0.060627581  | 0.100761308 |   0.059279532    |
|     MN346699     | Cote dIvoire |    2017/1/23    |  IIa  | unassign | 0.060627581  | 0.100761308 |   0.059279532    |
|     MN346700     | Cote dIvoire |    2017/1/31    |  IIa  | unassign | 0.060627581  | 0.100761308 |   0.059279532    |

## Reference

Yang, H.-C. et al (2022) **Monkeypox genome contains variable number
tandem repeats enabling accurate virus tracking.**
