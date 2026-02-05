# Get cost data from DAEDALUS scenarios

Get cost data from DAEDALUS scenarios

## Usage

``` r
get_cost_data(dt, disease_tags = "default", format = c("long", "wide"))
```

## Arguments

- dt:

  A `<data.table>` resulting from
  [`run_scenarios()`](https://jameel-institute.github.io/daedalus.compare/reference/run_scenarios.md).

- disease_tags:

  A character vector giving names for replicates within each scenario.

- format:

  A string for whether the data should be returned in `"long"` or
  `"wide"` format. Returning a 'wide' data.frame is only recommended
  when there is a small number of `disease_tags`, typically a triplet of
  "low", "medium", and "high" risk scenarios.

## Value

A `<data.frame>` summarising epidemic costs over the timeframe of the
modelled epidemics.

A `<data.frame>` of the costs for each model scenario.
