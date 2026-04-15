# NETDOM
Implementation of NETDOM in R

## Details about `statFun` options
The `statFun` argument determines how location-specific test statistics ("T(v)") are computed. We've provided several built-in options (e.g., `statFun='lm'` or `statFun='gam.deltaRsq'`), and there is also an option to specify a custom one (`statFun='custom'`).
The code for NETDOM is largely derived from the work of Sarah Weinstein on her method NEST ([@smweinst](https://github.com/smweinst/NEST)).

Example code for implementing NETDOM in R.
================

#### Step 1: Install the package

``` r
if (!require(devtools)){
  install.packages("devtools")
}
devtools::install_github("Nhillman19/NETDOM")
library(NETDOM)
```

#### Step 2: Prepare the data

Import necessary libraries for NETDOM:

Parallel processing can be used to distribute permutations across multiple cores and thereby speed up computation.

``` r
library(parallel)
library(NETDOM)
```

First load preprocessed brain imaging data and phenotype (e.g., demographic or diagnosis) information:

- Brain imaging data (e.g., freesurfer or cifti data) is loaded
- Network labels (network_labels) are loaded from yeo_7 network and
  filtered to remove unspecified indices (-1).
- A list (net) of a vector (net_7) is created for the specific default
  network, indicating vertex membership.

``` r
cifti_data <- load_dataset() # (replace with code to load your data)
network_labels <- load_yeo_7() # (replace with code to load your network labels)
idx <- which(network_labels != -1) # identify which locations should be ignored (e.g., medial wall)
network_label <- network_label[idx] # remove labels outside idx
X <- cifti_data[, idx] # subset image (X) locations to idx
net_7 <- as.integer(network_label == 7) # create a binary vector with 1's at locations corresponding to network 7 and 0's at locations outside network 7

# List networks of interest (could be more than one)
net <- list(
    Default <- net_7 
)
```

- Phenotype of interest (y) and covariates (Z)

``` r
y <- load_y() # replace with code to load your phenotype (e.g., age)
Z <- load_Z() # replace with code to load other covariates/confounders (e.g., sex)
```

After loading the data, we recommend checking that the dimensions of your input data are correct:
- dimension of `X` should be N x P (number of participants x number of image locations)
- dimension of phenotype vector `y` should be N
- dimensions of covariate/confounder vector or matrix `Z` should be `N` (x however many covariates you're adjusting for)

#### Step 3: Define the following dictionary of arguments passed to the NETDOM method. The args can be defined as follows, assuming vertex-wise linear models will be fit to estimate local brain-phenotype associations (i.e., specifying statFun=“lm” in step 4.).
- `X`: N x P matrix of P imaging features (e.g., vertices) for N
  participants
- `y`: N-dimensional vector of the phenotype of interest (i.e., testing
  enrichment of X-y associations).
- `Z`: Optional. Specify one or more covariates (matrix with N rows and q
  columns for q covariates). Default is NULL (no covariates to be
  included).
- `FL`: Optional. Default is FALSE, set to TRUE to use the Freedman-Lane
  procedure to account for dependence between covariates in permutation.
- `n.perm`: Optional. Default is 999, with the smallest possible p-value
  of 1/1000.

``` r
args.lm <- list(
    X = X, # brain measurements (dimension N subjects x P image locations)
    y = y, # phenotype of interest (dimension N)
    Z = Z, # covariates (dimension (N x # number of covariates)
    type = "coef", # what type of test statistic to extract from linear regression model (note: if using a different type of model for statistic, this may be different. see source code)
    n.perm = 999 # how many permutations to use to obtain null distribution
)
```
Note: non-linear regression-based statistics can also be used. This example is just for the regression-based statistic.

#### Step 4: Apply NETDOM to test enrichment of brain-phenotype associations in specified networks.
Arguments for NETDOM function: 
- `statFun`: specify the method to get vertex-level test statistics (e.g., "lm.fast"). Must correspond to a statFun R script (e.g., statFun.lm.R or statFun.gam.mvwald.R, or another one customized by the user)
- `args`: arguments needed for whatever was specified as statFun
- `net.maps`: list of binary vector(s) indicating locations inside (1) or outside (0) network(s) of interest.
- `direction`: Specifies which direction (e.g., "right" or "left") to test for when conducting a one-sided enrichment test
- `n.cores`: specify the number of CPU cores to be employed for parallel processing tasks within the function
- `seed`: random seed for reproducible permutation. Default is NULL, but we recommend setting one.
- `what.to.return`: specify what values to return. "everything" will include p-value, enrichment score, null distribution, etc. If left unspecified, the default is to return only the p-value for each network.
``` r
out <- NETDOM(statFun = "lm.fast",
            args = args.lm, # arguments defined above (specific to statFun="lm.fast" setting)
            net.maps = net, 
            direction = "right",
            n.cores = 1, 
            seed = 10, 
            what.to.return = "everything")

```

## Citation
If you use our method in any project or publication, please use the following citation [Testing for Network Specificity in Brain-Behavior Associations Using Ordinal Dominance Curves](http://dx.doi.org/10.1002/hbm.70493).

```
@article{Hillman2026NETDOM,
  title = {Testing for Network Specificity in Brain-Behavior Associations Using Ordinal Dominance Curves},
	author = {Hillman, Noah and Weinstein, Sarah M. and Bagautdinova, Joëlle and Sun, Kevin Y. and Cieslak, Matthew and Salo, Taylor and Fan, Yong and Keller, Arielle S. and Alexander-Bloch, Aaron F. and Vandekar, Simon N. and Raznahan, Armin and Satterthwaite, Theodore D. and Shou, Haochang and Shinohara, Russell T.},
	journal = {Human Brain Mapping},
	note = {e70493 HBM-25-0951},
	volume = {47},
	number = {5},
	pages = {e70493},
	year = {2026}}
```

