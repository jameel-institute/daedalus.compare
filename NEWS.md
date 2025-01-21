# daedalus.compare 0.0.1

This is the first version of _daedalus.compare_ and adds:

1. A scenario-runner function `run_scenarios()` that wraps around `daedalus::daedalus_rtm()` and allows running multiple scenarios of pandemic mitigation responses with uncertainty in the infection parameters.

2. Multiple functions downstream of `run_scenarios()` that help get data from scenarios in formats suitable for plotting or summarisation:
`get_epicurve_data()`, `get_summary_data()`, `get_cost_data()`, and `get_econ_cost_data()`.

3. Helper functions and constants: `ci()`, `make_infection_samples()`, and `NAMES_VECTOR_INF_PARAMS`.

4. Basic tests and function documentation.

* This project now includes a
   [`NEWS.md`](https://r-pkgs.org/other-markdown.html#sec-news) file to inform
   users about changes and new features.
