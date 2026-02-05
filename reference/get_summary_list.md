# Get epidemic summary data from a list of model outputs

Get epidemic summary data from a list of model outputs

## Usage

``` r
get_summary_list(l, names, ...)
```

## Arguments

- l:

  A list of `<daedalus_output>` objects; each object will be passed to
  [`daedalus::get_costs()`](https://jameel-institute.github.io/daedalus/reference/daedalus_costs.html).

- names:

  A character vector, intended to apply to elements of `l`.

- ...:

  Additional arguments passed to
  [`daedalus::get_epidemic_summary()`](https://jameel-institute.github.io/daedalus/reference/epi_output_helpers.html).

## Value

A `<data.table>`.
