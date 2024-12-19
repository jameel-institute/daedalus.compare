# TODO: move setup to separate file
gbr <- daedalus::daedalus_country("GBR")
infection <- daedalus::daedalus_infection("influenza_2009")

low_severity <- daedalus::get_data(daedalus_infection(""))
high_severity <- daedalus::get_data(daedalus_infection("sars_cov_1"), "gamma_H")

r0_samples <- round(withr::with_seed(1, rnorm(100, r0, 0.2)), 2)
r0_interval <- ci(r0_samples)
r0_limits <- c(
  lower = r0 - r0_interval,
  upper = r0 + r0_interval
)

# estimate final sizes in an unmitigated pandemic
lapply(r0_limits, finalsize::final_size)

disease_x <- lapply(
  r0_limits, function(x) {
    daedalus_infection("influenza_2009", r0 = x)
  }
)

disease_tags <- c("lower", "upper")
names(disease_x) <- disease_tags
