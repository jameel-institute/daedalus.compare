# Changelog

## daedalus.compare 0.0.5

This patch version allows users to run multiple vaccination scenarios by
passing a list of `<daedalus_vaccination>` objects to
[`run_scenarios()`](https://jameel-institute.github.io/daedalus.compare/reference/run_scenarios.md).

- Some reorganisation of data processing functions to handle scenarios
  that vary across vaccination strategies;

- Reorganised function documentation.

## daedalus.compare 0.0.4

This patch version updates {daedalus.compare} to work with {daedalus}
v0.3.0.

### Breaking changes

- [`run_scenarios()`](https://jameel-institute.github.io/daedalus.compare/reference/run_scenarios.md)
  only accepts `<daedalus_npi>`, a list of `<daedalus_npi>` or `NULL` as
  inputs to argument `response_strategy()`; in addition, the NPIs must
  be time-limited and reactive NPIs are not allowed.

### Other changes

- Tests and documentation have been changed to reflect the change in
  accepted NPI input types.

- A helper function to check that the NPI input is acceptable.

## daedalus.compare 0.0.3

Restores ability to get epidemic summary disaggregated by age and other
groups by passing arguments to
[`daedalus::get_epidemic_summary()`](https://jameel-institute.github.io/daedalus/reference/epi_output_helpers.html)
via `...` in
[`get_summary_data()`](https://jameel-institute.github.io/daedalus.compare/reference/get_summary_data.md).

## daedalus.compare 0.0.2

Updates for compatibility with *daedalus* \>= v0.2.16 which
re-introduces list-infection inputs.

- Input checking in
  [`run_scenarios()`](https://jameel-institute.github.io/daedalus.compare/reference/run_scenarios.md)
  removed to rely on errors bubbled up from *daedalus* functions.

## daedalus.compare 0.0.1

This is the first version of *daedalus.compare* and adds:

1.  A scenario-runner function
    [`run_scenarios()`](https://jameel-institute.github.io/daedalus.compare/reference/run_scenarios.md)
    that wraps around `daedalus::daedalus_rtm()` and allows running
    multiple scenarios of pandemic mitigation responses with uncertainty
    in the infection parameters.

2.  Multiple functions downstream of
    [`run_scenarios()`](https://jameel-institute.github.io/daedalus.compare/reference/run_scenarios.md)
    that help get data from scenarios in formats suitable for plotting
    or summarisation:
    [`get_epicurve_data()`](https://jameel-institute.github.io/daedalus.compare/reference/get_epicurve_data.md),
    [`get_summary_data()`](https://jameel-institute.github.io/daedalus.compare/reference/get_summary_data.md),
    [`get_cost_data()`](https://jameel-institute.github.io/daedalus.compare/reference/get_cost_data.md),
    and
    [`get_econ_cost_data()`](https://jameel-institute.github.io/daedalus.compare/reference/get_econ_cost_data.md).

3.  Helper functions and constants:
    [`ci()`](https://jameel-institute.github.io/daedalus.compare/reference/ci.md),
    [`make_infection_samples()`](https://jameel-institute.github.io/daedalus.compare/reference/make_infection_samples.md),
    and `NAMES_VECTOR_INF_PARAMS`.

4.  Basic tests and function documentation.

- This project now includes a
  [`NEWS.md`](https://r-pkgs.org/other-markdown.html#sec-news) file to
  inform users about changes and new features.
