## README
Genotyping Variable Number Tandem Repeats (VNTR) for the genome sequence of monkeypox virus (MPXV)

*2022-10-26*

## Introduction


## Installation

## Usage

## Arguments

`data` sequences from a file in FASTA format

`VNTR`

`match_s` matching weight

`mismatch_s` mismatching penalty

`regionStart`

`regionEnd`

`baseonly` logical. If TRUE, only uses the letters in base alphabet i.e. A,C,G,T.

`VNTRoutput` logical. If TRUE, export output to .csv file.

`finder`
`brief`


## Value

The output of `VNTR` is a matrix, which rows are strains and columns contain the following items:
* `ID` sequence name
* `r` the number of tandem repeats
* `match` the number of matches
* `mismatch` the number of mismatches
* `indel` the number of indels
* `score` alignment score for VNTR region
* `start_pos` start position of the VNTR region 

## Examples

```{Read sequences from a file in FASTA format}
data <- read.fasta("example.fasta", as.string = T)

length(data)
```

```{r example}
out <- VNTR(data, STR=STR, match_s=match_s, mismatch_s=mismatch_s, 
            regionStart=regionStart, regionEnd=regionEnd,baseonly = baseonly,VNTRoutput=VNTRoutput,finder=finder,brief=brief)

```
