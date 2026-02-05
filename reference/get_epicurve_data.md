# Get epidemiological curves from model data

Get epidemiological curves from model data

## Usage

``` r
get_epicurve_data(dt, disease_tags = "default", format = c("long", "wide"))
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

A `<data.frame>` of individuals in each epidemiological compartments
over the timeframe of the modelled epidemics.
