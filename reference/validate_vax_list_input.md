# Check multi-vaccination inputs

A helper function to check whether a list of vaccination objects has
been passed as input.

## Usage

``` r
validate_vax_list_input(x)
```

## Arguments

- x:

  A list to be checked as a list of `<daedalus_vaccination>`. Must have
  a minimum length of 2 elements.

## Value

Either `x`, if it is a list of `<daedalus_vaccination>`, or an error
side-effect if not.
