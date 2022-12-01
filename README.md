## README
Description

*2022-12-01*

## Installation

To install **`MPXV`** from [**Github**](https://github.com/teresayang/MPXV_VNTR.git):

```{r Installation from GitHub, eval = FALSE}
if(!require(devtools)) install.packages("devtools")
remotes::install_git("https://github.com/teresayang/MPXV_VNTR.git",
                     credentials=git2r::cred_user_pass("teresayang", "password"))
```

To load the installed **`MPXV`** in R:

```{r Load MPXV, eval = FALSE}
library(scDEseq)
```
