# Get economic cost details from an output data.table

Get economic cost details from an output data.table

## Usage

``` r
get_econ_cost_data(dt)
```

## Arguments

- dt:

  A `<data.table>` resulting from
  [`run_scenarios()`](https://jameel-institute.github.io/daedalus.compare/reference/run_scenarios.md).

## Value

A `<data.frame>` summarising economic costs over the timeframe of the
modelled epidemics, broken down into costs due to restrictions and
illness- related absences.
