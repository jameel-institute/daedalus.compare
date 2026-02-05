# Get summary data from DAEDALUS scenarios

Get summary data from DAEDALUS scenarios

## Usage

``` r
get_summary_data(dt, disease_tags, format = c("long", "wide"), ...)
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

- ...:

  Additional arguments passed to
  [`daedalus::get_epidemic_summary()`](https://jameel-institute.github.io/daedalus/reference/epi_output_helpers.html).
  If passing `groups`, only `format = "long"` is supported.

## Value

A `<data.frame>` summarising epidemiological outcomes: cumulative
infections, deaths, and hospitalisations over the timeframe of the
modelled epidemics.
