# Generate multiple `<daedalus_infection>`s from parameter distributions

Generate multiple `<daedalus_infection>`s from parameter distributions

## Usage

``` r
make_infection_samples(
  name,
  param_distributions,
  param_ranges = NULL,
  samples = 100
)
```

## Arguments

- name:

  An infection name from among
  [`daedalus.data::epidemic_names`](https://jameel-institute.github.io/daedalus.data/reference/epidemic_data.html).

- param_distributions:

  A named list of `<distribution>` class objects provided by
  distributional, with names corresponding to
  [`daedalus.data::infection_parameter_names`](https://jameel-institute.github.io/daedalus.data/reference/epidemic_data.html).
  These are used to generate `samples` draws from each distribution for
  the corresponding infection parameters. Arguments which vary by age
  are supported, with the drawn and scaled value treated as the median
  of the profile vector. See **Examples**.

- param_ranges:

  An optional named list of two-element vectors, giving the ranges to
  which samples drawn using `param_distributions` should be rescaled.
  This allows easy rescaling of values drawn from bounded distributions
  (e.g. the Beta distribution). If passed, the list must have names
  corresponding to `param_distribution` names. When empty, the generated
  values are used without scaling. See **Examples**.

- samples:

  The number of samples to generate.

## Value

A list of `<daedalus_infection>` objects.

## Examples

``` r
# make 10 infection objects varying in R0, with a skewed distribution
# scaled between 1.0 and 2.0
make_infection_samples(
  "influenza_2009",
  samples = 3,
  list(
    r0 = distributional::dist_beta(2, 5)
  ),
  list(
    r0 = c(0.1, 0.2)
  )
)
#> [[1]]
#> <daedalus_infection>
#> • Epidemic name: influenza_2009
#> • R0: 0.2
#> • sigma: 0.909
#> • p_sigma: 0.669
#> • epsilon: 0.58
#> • rho: 0.003
#> • eta: 0.003, 0.003, 0.001, and 0.006
#> • hfr: 0.04, 0.063, 0.041, and 0.634
#> • gamma_Ia: 0.4
#> • gamma_Is: 0.4
#> • gamma_H_recovery: 0.2
#> • gamma_H_death: 0.2
#> 
#> [[2]]
#> <daedalus_infection>
#> • Epidemic name: influenza_2009
#> • R0: 0.1
#> • sigma: 0.909
#> • p_sigma: 0.669
#> • epsilon: 0.58
#> • rho: 0.003
#> • eta: 0.003, 0.003, 0.001, and 0.006
#> • hfr: 0.04, 0.063, 0.041, and 0.634
#> • gamma_Ia: 0.4
#> • gamma_Is: 0.4
#> • gamma_H_recovery: 0.2
#> • gamma_H_death: 0.2
#> 
#> [[3]]
#> <daedalus_infection>
#> • Epidemic name: influenza_2009
#> • R0: 0.182723943621094
#> • sigma: 0.909
#> • p_sigma: 0.669
#> • epsilon: 0.58
#> • rho: 0.003
#> • eta: 0.003, 0.003, 0.001, and 0.006
#> • hfr: 0.04, 0.063, 0.041, and 0.634
#> • gamma_Ia: 0.4
#> • gamma_Is: 0.4
#> • gamma_H_recovery: 0.2
#> • gamma_H_death: 0.2
#> 
```
