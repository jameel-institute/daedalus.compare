# Check multi-NPI inputs

A helper function to check whether a list of NPI objects has been passed
as input. This could live in daedalus similar to
`daedalus:::validate_infection_list_input()`, but lives here instead.

## Usage

``` r
validate_npi_list_input(x, timed_only = TRUE)
```

## Arguments

- x:

  A list to be checked as a list of `<daedalus_npi>`. Must have a
  minimum length of 2 elements.

- timed_only:

  A boolean defaulting to `TRUE` for whether only time-limited NPIs are
  allowed. This is typically the case for a real-time modelling
  exercise.

## Value

Either `x`, if it is a list of `<daedalus_npi>`, or an error side-effect
if not.
