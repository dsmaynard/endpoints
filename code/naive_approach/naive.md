Fitting endpoints: naive approach
================
Daniel Maynard, Zachary Miller and Stefano Allesina

Here we show how to fit the endpoints or predict out of fit the endpoints using the naive approach presented in the Materials and Methods. This method is similar to a linear regression. The advantage is that it is quite fast and allows one to understand the mechanics of the model very well. The downside is that it can return nonsensical solutions (e.g., negative abundances) if the data deviate from what expected.

All the code needed for the analysis is in the file `naive.R`. The libraries `MASS` and `tidyverse` need to be installed for the code to run.

``` r
source("naive.R")
```

### Step 1: label the data

First, call the function `prepare_data` using one of the available data sets, or any other data set organized in the same manner. The data file should be a `csv` with one column for each species, and one row for each recorded endpoint. The header should specify the name of the species. For example:

``` r
dt <- read_csv("../../data/Kuebbing_plants/natives.csv")
dt %>% sample_n(10) # show 10 endpoints sampled at random
```

    ## # A tibble: 10 x 4
    ##       as    fa    la    po
    ##    <dbl> <dbl> <dbl> <dbl>
    ##  1 1.80   0    0      0   
    ##  2 0      5.47 0.780  0   
    ##  3 0      0    1.44   0   
    ##  4 1.19   0    0      2.23
    ##  5 0      2.99 0      2.04
    ##  6 0.487  0    0      3.07
    ##  7 0.472  1.52 0      2.39
    ##  8 0.602  1.94 0.182  2.75
    ##  9 0      3.22 0      2.53
    ## 10 2.21   0    0.288  0

The function `prepare_data` simply adds a column containing a label for the community. For example:

``` r
dt <- prepare_data("../../data/Kuebbing_plants/natives.csv")
dt %>% sample_n(10) # show 10 endpoints sampled at random
```

    ## # A tibble: 10 x 5
    ##       as    fa    la    po community  
    ##    <dbl> <dbl> <dbl> <dbl> <chr>      
    ##  1 0      4.85 1.35   0    fa-la      
    ##  2 1.51   0    0.309  2.90 as-la-po   
    ##  3 0.602  1.94 0.182  2.75 as-fa-la-po
    ##  4 3.20   0    0.591  0    as-la      
    ##  5 1.13   3.89 0      1.74 as-fa-po   
    ##  6 2.74   1.72 0.621  0    as-fa-la   
    ##  7 0      0    0.752  5.23 la-po      
    ##  8 1.80   0    0      0    as         
    ##  9 2.88   3.13 0      0    as-fa      
    ## 10 0      0    0.226  2.58 la-po

### Step 2: fitting all endpoints

To use all of the available data to parameterize the matrix *B*, simply call:

``` r
out <- fit_endpoints(dt)
```

The function returns a list, containing the best-fitting matrix *B*:

``` r
out$B
```

    ##            [,1]        [,2]       [,3]       [,4]
    ## [1,] -0.3218798 -0.08661629 -0.3816623 -0.1730324
    ## [2,] -0.2149084 -0.14837835 -0.3080815 -0.1697552
    ## [3,] -0.2454120 -0.06011190 -0.7612364 -0.1560388
    ## [4,] -0.2066253 -0.09150894 -0.2194261 -0.2343866

And a tibble containing, for each community, the observed vs. fitted values:

``` r
out$results
```

    ## # A tibble: 289 x 4
    ##    community species observed predicted
    ##    <chr>     <chr>      <dbl>     <dbl>
    ##  1 as        as          1.47      3.11
    ##  2 as        as          1.80      3.11
    ##  3 as        as          1.90      3.11
    ##  4 as        as          2.06      3.11
    ##  5 as        as          2.27      3.11
    ##  6 as        as          2.31      3.11
    ##  7 as        as          2.36      3.11
    ##  8 as        as          2.46      3.11
    ##  9 as        as          2.81      3.11
    ## 10 as        as          3.18      3.11
    ## # ... with 279 more rows

The results can be visualized using:

``` r
stats_and_plot_results(out$results)
```

    ## [1] "Correlation"
    ## [1] 0.8600766

![](naive_files/figure-markdown_github/unnamed-chunk-7-1.png)

### Step 3: predicting out of fit

The function `predict_out_of_fit` implements the out-of-fit prediction presented in the manuscript. In turn, each endpoint (and all its replicates) is removed, and the remaining data is used to predict the coexistence abundance of the remaining species. Negative abundances mean lack of coexistence for the community. The program sets all abundances to -1 in this case.

For example:

``` r
loo <- predict_out_of_fit(dt)
stats_and_plot_results(loo)
```

    ## [1] "Correlation"
    ## [1] 0.7607995

![](naive_files/figure-markdown_github/unnamed-chunk-8-1.png)

### Poor performance of naive approach

While the results are good for the data by Kuebbing *et al.* (2015), presented above, the other data sets show that this approach cannot successfully predict the endpoints. For example:

``` r
# in-fit Rakowski and Cardinale (2016)
dtrc <- prepare_data("../../data/Cardinale_daphnia/cd_algae.csv")
out <- fit_endpoints(dtrc)
stats_and_plot_results(out$results) 
```

    ## [1] "Correlation"
    ## [1] 0.7291944

![](naive_files/figure-markdown_github/unnamed-chunk-9-1.png)

Where several communities that were observed experimentally are predicted to collapse. Similar problems are found when using the data from Pennekamp et al. (2018):

``` r
# in-fit Pennekamp et al. 2018
dtp <- prepare_data("../../data/Pennekamp_protists/temp_17.csv")
out <- fit_endpoints(dtp)
stats_and_plot_results(out$results) 
```

    ## [1] "Correlation"
    ## [1] 0.2783566

![](naive_files/figure-markdown_github/unnamed-chunk-10-1.png)
