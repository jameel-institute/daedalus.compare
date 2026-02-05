# Get incidence data from a list of model outputs

Get incidence data from a list of model outputs

## Usage

``` r
get_epidata_list(l, names)
```

## Arguments

- l:

  A list of `<daedalus_output>` objects; each object will be passed to
  [`daedalus::get_costs()`](https://jameel-institute.github.io/daedalus/reference/daedalus_costs.html).

- names:

  A character vector, intended to apply to elements of `l`.

## Value

A `<data.table>`.
