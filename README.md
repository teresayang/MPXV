## README
Description

*2022-12-01*

## Installation

To install **`MPXV`** from [**Github**](https://github.com/teresayang/MPXV_VNTR.git):

```{r Installation from GitHub, eval = FALSE}
if(!require(devtools)) install.packages("devtools")
remotes::install_git("https://github.com/teresayang/MPXV_VNTR.git")
```

To load the installed **`MPXV`** in R:

```{r Load MPXV, eval = FALSE}
library(MPXV)
```

## Usage



## Output

### VNTR caller: ###

| ID       | r  | match | mismatch | indel | score | start_pos  |
|:--------:|:---:|:-----:|:--------:|:-----:|:-----:|:----------:|
| MT903346 | 11 | 11    | 1        | 0     | 17    | 135074     |
| MT903347 | 11 | 11    | 1        | 0     | 17    | 135075     |
| MT903348 | 11 | 11    | 1        | 0     | 17    | 135075     |
| MN346690 | 11 | 11    | 1        | 0     | 17    | 134676     |
| MN346692 | 11 | 11    | 1        | 0     | 17    | 134636     |

| ID       | nt133095_T | nt150554_TATGATGGA | nt173267_AT | nt179074_ATATACATT  |
|:--------:|:----------:|:------------------:|:-----------:|:-------------------:|
| MT903346 | 11         | 3                  | 4           | 26                  |
| MT903347 | 11         | 3                  | 4           | 26                  |
| MT903348 | 11         | 3                  | 4           | 26                  |
| MN346690 | 11         | 5                  | 4           | 14                  |
| MN346692 | 11         | 5                  | 4           | 15                  |

```
>MT903346
TTTTTTTTCTTT
>MT903347
TTTTTTTTCTTT
>MT903348
TTTTTTTTCTTT
>MN346690
TTTTTTTTCTTT
>MN346692
TTTTTTTTCTTT
```

### VNTR tracker: ###

| Accession        | Country      | Collection.date | clade | lineage  | Distance_PIC | Distance_L  | Distance_entropy  |
|:----------------:|:------------:|:---------------:|:-----:|:--------:|:------------:|:-----------:|:-----------------:|
| MN346690         | Cote dIvoire | 2017/3/13       | IIa   | unassign | 0            | 0           | 0                 |
| MN346692         | Cote dIvoire | 2017/3/5        | IIa   | unassign | 0.004872201  | 0.007389163 | 0.004750482       |
| MN346694         | Cote dIvoire | 2017/3/28       | IIa   | unassign | 0.004872201  | 0.007389163 | 0.004750482       |
| MN346696         | Cote dIvoire | 2017/4/3        | IIa   | unassign | 0.019488804  | 0.02955665  | 0.019001927       |
| KP849470         | Cote dIvoire | 1971            | IIa   | unassign | 0.030844508  | 0.073219884 | 0.030573661       |
| MN346702         | Cote dIvoire | 2018/5/11       | IIa   | unassign | 0.048722011  | 0.073891626 | 0.047504818       |
| EPI_ISL_13544250 | Canada       | 2022/5/30       | IIb   | B.1.4    | 0.048765167  | 0.07320153  | 0.049289862       |
| MN346698         | Cote dIvoire | 2017/1/24       | IIa   | unassign | 0.060627581  | 0.100761308 | 0.059279532       |
| MN346699         | Cote dIvoire | 2017/1/23       | IIa   | unassign | 0.060627581  | 0.100761308 | 0.059279532       |
| MN346700         | Cote dIvoire | 2017/1/31       | IIa   | unassign | 0.060627581  | 0.100761308 | 0.059279532       |
