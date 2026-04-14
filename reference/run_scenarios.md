# Run multiple DAEDALUS scenarios

Run multiple DAEDALUS scenarios

## Usage

``` r
run_scenarios(
  country,
  infection,
  response_strategy = NULL,
  vaccination_strategy = NULL,
  time_end = 100,
  initial_state_manual = NULL
)
```

## Arguments

- country:

  A country or territory object of class `<daedalus_country>`, **or** a
  country or territory name from those included in the package; see
  [daedalus.data::country_names](https://jameel-institute.github.io/daedalus.data/reference/country_names_codes.html),
  **or** a country ISO2 or ISO3 code; see
  [daedalus.data::country_codes_iso2c](https://jameel-institute.github.io/daedalus.data/reference/country_names_codes.html)
  and
  [daedalus.data::country_codes_iso3c](https://jameel-institute.github.io/daedalus.data/reference/country_names_codes.html).
  Country-specific data such as the community and workplace contacts,
  the demography, and the distribution of the workforce into economic
  sectors is automatically accessed from package data for the relevant
  country name if it is passed as a string. To override package defaults
  for country characteristics, pass a `<daedalus_country>` object
  instead. See
  [`daedalus_country()`](https://jameel-institute.github.io/daedalus/reference/class_country.html)
  for more.

- infection:

  A list of `<daedalus_infection>` objects. Must have a minimum length
  of 2.

- response_strategy:

  A time-limited response strategy specified as a single
  `<daedalus_npi>` created using
  [`daedalus::daedalus_timed_npi()`](https://jameel-institute.github.io/daedalus/reference/class_npi.html),
  or a list of such objects, or `NULL` for no response. Lists may
  include `NULL`, which is useful when comparing response scenarios
  against the no-response counterfactual. If a list is passed, all
  elements must be named, or the list may have no names, in which case
  synthetic names will be assigned and a message printed to screen.
  Lists with some elements named and some unnamed are not accepted.

- vaccination_strategy:

  A vaccination strategy specified as a single `<daedalus_vaccination>`
  created using
  [`daedalus::daedalus_vaccination()`](https://jameel-institute.github.io/daedalus/reference/class_vaccination.html),
  or a list of such objects, or `NULL` for no response. Lists may
  include `NULL`, which is useful when comparing vaccination scenarios
  against the no-vaccination counterfactual. If a list is passed, all
  elements must be named, or the list may have no names, in which case
  synthetic names will be assigned and a message printed to screen.
  Lists with some elements named and some unnamed are not accepted.

- time_end:

  A vector of integer-ish numbers giving the durations over which to run
  scenarios. Each scenario is run for each `time_end`. The intention is
  to be able to run response scenarios for different durations, to be
  able to generate a time-series of how different costs accumulate.

- initial_state_manual:

  An optional **named** list with the names `p_infectious`,
  `p_asymptomatic`, and `p_immune`. `p_infectious` and `p_asymptomatic`
  give the proportion of infectious and symptomatic individuals in each
  age group and economic sector. Defaults to `1e-6` and `0.0`
  respectively. `p_immune` may be a single number in the range
  `0.0 <= p_immune <= 1.0` or a 4-element vector in that range (the
  number of age groups in the model), for the proportion of individuals
  in the population or in each age group that have some pre-existing
  immunity to infection (reduced susceptibility). See **Details** for
  more.

## Value

A `<data.table>`, with data held in list-columns.
