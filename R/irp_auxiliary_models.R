#### Models ####


#' Model: Enthalpy of combustion per mol C versus electrons transferred during combustion per mol C
#' 
#' @export
irp_make_m1 <- function(irp_db_obigt_organic, irp_mcmc_settings) {
  
  res_prior <-
    c(
      brms::prior_string("normal(0.0, 1.0)", class = "b"),
      brms::prior_string("normal(-2.0, 0.5)", class = "Intercept"), # ---note: exp(-10) is near 0, a value which is plausible based on prior knowledge that no electron transfer means no reaction means a small absolute enthalpy of combustion
      brms::prior_string(paste0("gamma(1.0, ", 1.0/500.0, ")"), class = "shape", lb = 0.0)
    )
  
  # define scaling variables
  res_config <- 
    tibble::tibble(
      y_min = 
        irp_db_obigt_organic |> 
        dplyr::filter(! is.na(enthalpy_of_combustion_per_C) & as.numeric(enthalpy_of_combustion_per_C) < 0) |>
        dplyr::pull(enthalpy_of_combustion_per_C) |>
        min(),
      x_max = 
        irp_db_obigt_organic |> 
        dplyr::filter(! is.na(enthalpy_of_combustion_per_C) & as.numeric(enthalpy_of_combustion_per_C) < 0) |>
        dplyr::pull(E_per_C) |>
        max()
    )
  
  # fit the model
  res_model <-  
    brms::brm(
      enthalpy_of_combustion_per_C ~ log(E_per_C), 
      data = 
        irp_db_obigt_organic |>
        dplyr::mutate(
          enthalpy_of_combustion_per_C = enthalpy_of_combustion_per_C / res_config$y_min,
          E_per_C = E_per_C / res_config$x_max,
          dplyr::across(
            dplyr::all_of(c("enthalpy_of_combustion_per_C", "E_per_C")),
            as.numeric
          )
        ) |>
        dplyr::filter(enthalpy_of_combustion_per_C > 0),
      family = Gamma(link = "log"),
      chains = irp_mcmc_settings$chains,
      iter = irp_mcmc_settings$iter,
      warmup = irp_mcmc_settings$warmup,
      prior = res_prior,
      cores = 1L
    )
  
  list(
    res_prior = res_prior,
    res_config = res_config,
    res_model = res_model
  )
  
}


#' Model: Entropy of formation per mol C versus sum of standard atomic entropies per mol C
#' 
#' @export
irp_make_m2 <- function(irp_db_obigt_organic, irp_d_battley_1999_table_1, irp_mcmc_settings) {
 
  # merge data
  res_data <- 
    dplyr::bind_rows(
      irp_db_obigt_organic |>
        dplyr::filter(state == "cr"),
      irp_d_battley_1999_table_1
    ) |>
    dplyr::filter(! duplicated(formula), fromLast = TRUE) |>
    dplyr::filter(!is.na(entropy_of_formation_per_C) & as.numeric(entropy_of_formation_per_C) < 0 & ! is.na(sum_of_element_standard_entropies_per_C) & ! is.infinite(sum_of_element_standard_entropies_per_C))
  
  # priors
  res_prior <-
    c(
      brms::prior_string("normal(0.0, 1.0)", class = "b"),
      brms::prior_string("normal(-2.0, 0.5)", class = "Intercept"), # ---note: exp(-10) is near 0, a value which is plausible based on prior knowledge that no electron transfer means no reaction means a small absolute enthalpy of combustion
      brms::prior_string(paste0("gamma(1.0, ", 1.0/500.0, ")"), class = "shape", lb = 0.0)
    )
  
  # define scaling variables
  res_config <- 
    tibble::tibble(
      y_min = 
        res_data |> 
        dplyr::pull(entropy_of_formation_per_C) |>
        min(),
      x_max = 
        res_data |> 
        dplyr::pull(sum_of_element_standard_entropies_per_C) |>
        max()
    )
  
  # fit the model
  res_model <-  
    brms::brm(
      entropy_of_formation_per_C ~ 1 + log(sum_of_element_standard_entropies_per_C),
      data = 
        res_data |>
        dplyr::mutate(
          entropy_of_formation_per_C = entropy_of_formation_per_C / res_config$y_min,
          sum_of_element_standard_entropies_per_C = sum_of_element_standard_entropies_per_C / res_config$x_max,
          dplyr::across(
            dplyr::all_of(c("entropy_of_formation_per_C", "sum_of_element_standard_entropies_per_C")),
            as.numeric
          )
        ),
      family = Gamma(link = "log"),
      chains = irp_mcmc_settings$chains,
      iter = irp_mcmc_settings$iter,
      warmup = irp_mcmc_settings$warmup,
      prior = res_prior,
      cores = 1L
    )
  
  list(
    res_data = res_data,
    res_prior = res_prior,
    res_config = res_config,
    res_model = res_model
  )
  
}


#' Model: Volume fraction of different pore classes in peat
#' 
#' @export
irp_make_m3 <- function(irp_liu2019_pmird, irp_wang2015b, irp_mcmc_settings) {
  
  res_data <- 
    irp_liu2019_pmird |>
    dplyr::mutate(
      non_macroporosity = porosity - macroporosity,
      volume_fraction_solids = 1 - porosity
    ) |>
    dplyr::filter(! is.na(bulk_density) & ! is.na(porosity) & ! is.na(macroporosity)) |>
    dplyr::bind_rows(
      irp_wang2015b |>
        dplyr::rename(porosity = "porosity_helium") |>
        dplyr::mutate(
          dplyr::across(where(is.numeric), as.numeric),
          non_macroporosity = porosity - 0.001,
          volume_fraction_solids = 1 - porosity,
          macroporosity = 0 + 0.001,
          sample_type = "shale"
        )
    )
  
  res_prior <- 
    c(
      brms::prior_string(paste0("normal(", ohenery::inv_smax(mu = c(0.001, 0.001, 1-0.001*2))[[2]],", 2)"), dpar = "munonmacroporosity", class = "Intercept"),
      brms::prior_string(paste0("normal(", ohenery::inv_smax(mu = c(0.001, 0.001, 1-0.001*2))[[3]],", 2)"), dpar = "mumacroporosity", class = "Intercept"),
      brms::prior_string("normal(0, 2)", dpar = "munonmacroporosity", class = "b"),
      brms::prior_string("normal(0, 2)", dpar = "mumacroporosity", class = "b"),
      brms::prior_string("gamma(1.0, 1.0/200.0)", class = "phi")
    )
  
  # compute the Dirichlet model
  bind <- function(...) cbind(...)
  
  res_model <- 
    brms::brm(
      bind(volume_fraction_solids, non_macroporosity, macroporosity) ~ log(bulk_density) + bulk_density,
      data = 
        res_data |>
        dplyr::select(sample_type, bulk_density, porosity, macroporosity, non_macroporosity, volume_fraction_solids), 
      family = brms::dirichlet(),
      chains = irp_mcmc_settings$chains,
      iter = irp_mcmc_settings$iter,
      warmup = irp_mcmc_settings$warmup,
      prior = res_prior,
      cores = 1L,
      control = 
        list(
          adapt_delta = irp_mcmc_settings$adapt_delta, 
          max_treedepth = irp_mcmc_settings$max_treedepth
        )
    )
  
  list(
    res_data = res_data,
    res_prior = res_prior,
    res_config = tibble::tibble(),
    res_model = res_model
  )
  
}


#' Model: Saturated hydraulic conductivity
#' 
#' @export
irp_make_m4 <- function(irp_liu2019_pmird, irp_mcmc_settings) {
  
  res_data <- 
    irp_liu2019_pmird |>
    dplyr::filter(! is.na(bulk_density) & ! is.na(hydraulic_conductivity))
  
  res_prior <- 
    c(
      brms::prior_string(paste0("normal(", brms::logit_scaled(2000/3000),", 0.8)"), class = "Intercept"),
      brms::prior_string(paste0("normal(", 0.0,", 2.0)"), class = "b"),
      brms::prior_string(paste0("normal(", 0.0,", 1.0)"), class = "Intercept", dpar = "phi"),
      brms::prior_string(paste0("normal(", 0.0,", 2.0)"), class = "b", dpar = "phi")
    ) 
  
  # define scaling variables
  res_config <- 
    tibble::tibble(
      y_scale = 3000.0
    )
  
  res_model <- 
    brms::brm(
      brms::bf(
        hydraulic_conductivity ~ log(bulk_density) + bulk_density, 
        phi ~ log(bulk_density) + bulk_density
      ),
      data = 
        res_data |>
        dplyr::select(sample_type, bulk_density, hydraulic_conductivity) |>
        dplyr::mutate(hydraulic_conductivity = hydraulic_conductivity/res_config$y_scale), 
      family = brms::Beta(link = "logit", link_phi = "log"),
      chains = irp_mcmc_settings$chains,
      iter = irp_mcmc_settings$iter,
      warmup = irp_mcmc_settings$warmup,
      prior = res_prior,
      cores = 1L,
      control = 
        list(
          adapt_delta = irp_mcmc_settings$adapt_delta, 
          max_treedepth = irp_mcmc_settings$max_treedepth
        )
    )
  
  list(
    res_data = res_data,
    res_prior = res_prior,
    res_config = res_config,
    res_model = res_model
  )
  
}


#' Model: Specific heat capacity
#' 
#' @export
irp_make_m5 <- function(irp_d_gnatowski2022, irp_d_cp_air, irp_mcmc_settings) {
  
  res_data <- 
    dplyr::bind_rows(
      irp_d_gnatowski2022,
      irp_d_cp_air |>
        dplyr::mutate(
          temperature = temperature - 273.15
        ) %>%
        dplyr::filter(
          temperature >= min(irp_d_gnatowski2022$temperature, na.rm = TRUE) & temperature <= max(irp_d_gnatowski2022$temperature, na.rm = TRUE)
        )
    )
  
  res_prior <- 
    c(
      brms::prior_string(paste0("normal(", log(irp_d_cp_air$specific_heat_capacity[irp_d_cp_air$temperature == 270]),", 0.1)"), class = "Intercept"),
      brms::prior_string(paste0("normal(", 0,", 0.2)"), class = "sd", group = "sample_label", lb = 0),
      brms::prior_string(paste0("normal(", 0,", 1.0)"), class = "b"),
      brms::prior_string(paste0("gamma(1.0, 1.0/200.0)"), class = "shape", lb = 0.0)
    )
  
  # define scaling variables
  res_config <- 
    tibble::tibble(
      temperature_scale = 
        attr(scale(res_data$temperature, center = FALSE, scale = TRUE), "scaled:scale"),
      N_scale = max(res_data$N)
    )
  
  res_model <- 
    brms::brm(
      brms::bf(
        specific_heat_capacity | mi(specific_heat_capacity_error) ~ N * temperature + (1|sample_label)
      ),
      data = 
        res_data |>
        dplyr::select(specific_heat_capacity, specific_heat_capacity_error, N, temperature, sample_label) |>
        dplyr::mutate(
          N = N/res_config$N_scale,
          temperature = temperature/res_config$temperature_scale
        ), 
      family = Gamma(link = "log"),
      chains = irp_mcmc_settings$chains,
      iter = irp_mcmc_settings$iter,
      warmup = irp_mcmc_settings$warmup,
      prior = res_prior,
      cores = irp_mcmc_settings$chains,
      backend = "cmdstanr",
      control = 
        list(
          adapt_delta = irp_mcmc_settings$adapt_delta, 
          max_treedepth = irp_mcmc_settings$max_treedepth
        )
    )
  
  list(
    res_data = res_data,
    res_prior = res_prior,
    res_config = res_config,
    res_model = res_model
  )
  
}


#' Model: Dry thermal conductivity
#' 
#' @export
irp_make_m6 <- function(irp_d_oconnor2022, irp_mcmc_settings) {
  
  res_data <- 
    irp_d_oconnor2022 |>
    dplyr::filter(!is.na(bulk_density) & !is.na(dry_thermal_conductivity))
  
  res_prior <- 
    c(
      brms::prior_string(paste0("normal(", log(.9894 * 2.414 * 10^(-4) * 100),", 0.5)"), class = "Intercept"), #---note: thermal conductivity of dry air at 1 atm (@Hilsenrath.1955, with multiplication factor 2.414 * 10^(-4) from @Hilsenrath.1955 and additional multiplication factor 100 to normalize per m)
      brms::prior_string(paste0("normal(", 0.0,", 1.0)"), class = "b"),
      #brms::prior_string(paste0("gamma(1.0, 1.0/200.0)"), class = "shape"),
      brms::prior_string(paste0("normal(", 0.0,", 1.0)"), class = "Intercept", dpar = "shape"),
      brms::prior_string(paste0("normal(", 0.0,", 1.0)"), class = "b", dpar = "shape")
    )
  
  # define scaling variables
  res_config <- tibble::tibble()
  
  res_model <- 
    brms::brm(
      brms::bf(
        dry_thermal_conductivity ~ bulk_density + log(bulk_density),
        shape ~ log(bulk_density) + bulk_density
      ),
      data = res_data, 
      family = Gamma(link = "log"),
      chains = irp_mcmc_settings$chains,
      iter = irp_mcmc_settings$iter,
      warmup = irp_mcmc_settings$warmup,
      prior = res_prior,
      cores = irp_mcmc_settings$chains,
      backend = "cmdstanr",
      control = 
        list(
          adapt_delta = irp_mcmc_settings$adapt_delta, 
          max_treedepth = irp_mcmc_settings$max_treedepth
        )
    )
  
  list(
    res_data = res_data,
    res_prior = res_prior,
    res_config = res_config,
    res_model = res_model
  )
  
}


#### Helper functions ####

#' Takes a data frame of molar elemental compositions of samples and estimates the standard enthalpy of combustion per mol C
#' 
#' @param x A data frame with columns representing chemical elements (and named this way),
#' with units defined by elco (e.g. mol_C/mol_sample). Columns must match arguments in
#' `elco_get_electrons_tansferred_combustion()`.
#' 
#' @export
irp_estimate_dhc0_per_C <- function(x, irp_m1) {
  
  stopifnot(is.data.frame(x))
  
  # add predictors
  res <- 
    x |>
    dplyr::mutate(
      E_per_C = 
        elco_get_electrons_tansferred_combustion(C = C, H = H, N = N, O = O, P = P, S = S) |>
        magrittr::divide_by(C) |>
        magrittr::divide_by(irp_m1$res_config$x_max) |>
        as.numeric()
    )
  
  # make predictions
  res <- 
    res |>
    dplyr::mutate(
      # standard enthalpy of combustion per C
      dhc0_per_C = 
        irp_m1$res_model |>
        brms::posterior_predict(
          type = "response",
          newdata = res
        ) |>
        posterior::rvar() |>
        magrittr::multiply_by(as.numeric(irp_m1$res_config$y_min))
    )
  
  # add units
  posterior::draws_of(res$dhc0_per_C) <- units::set_units(posterior::draws_of(res$dhc0_per_C), "J/mol_C", mode = "standard") 
  
  res |>
    dplyr::select(-"E_per_C")
  
}


#' Takes a data frame of molar elemental compositions of samples and estimates the standard entropy of formation per mol C
#' 
#' @param x A data frame with columns representing chemical elements (and named this way),
#' with units defined by elco (e.g. mol_C/mol_sample). Columns must match arguments in
#' `elco_get_electrons_tansferred_combustion()`.
#' 
#' @export
irp_estimate_dsf0_per_C <- function(x, irp_m2) {
  
  stopifnot(is.data.frame(x))
  
  # standard entropies of elements
  irp_standard_entropies_elements <- 
    irp_make_standard_entropies_elements()
  
  # extract chemical elements
  res_cf <- 
    x |>
    dplyr::select(dplyr::any_of(irp_standard_entropies_elements$chemical_element))
  
  irp_standard_entropies_elements <- 
    irp_standard_entropies_elements |>
    dplyr::filter(chemical_element %in% colnames(res_cf)) |>
    dplyr::slice(match(chemical_element, colnames(res_cf)))
  
  
  # add predictors
  res <- 
    x |>
    dplyr::mutate(
      sum_of_element_standard_entropies =
        units::set_units(as.vector(as.matrix(res_cf) %*% irp_standard_entropies_elements$standard_entropy), "J/(K * mol_sample)", mode = "standard"),
      sum_of_element_standard_entropies_per_C = sum_of_element_standard_entropies/(res_cf$C / units::set_units(1, "mol_sample", mode = "standard")),
      sum_of_element_standard_entropies_per_C =
        sum_of_element_standard_entropies_per_C |>
        magrittr::divide_by(irp_m2$res_config$x_max) |>
        as.numeric()
    ) |>
    dplyr::select(-"sum_of_element_standard_entropies")
  
  # make predictions
  res <- 
    res |>
    dplyr::mutate(
      # standard enthalpy of combustion per C
      dsf0_per_C = 
        irp_m2$res_model |>
        brms::posterior_predict(
          type = "response",
          newdata = res
        ) |>
        posterior::rvar() |>
        magrittr::multiply_by(as.numeric(irp_m2$res_config$y_min))
    )
  
  # add units
  posterior::draws_of(res$dsf0_per_C) <- units::set_units(posterior::draws_of(res$dsf0_per_C), "J/K/mol_C", mode = "standard") 
  
  res |>
    dplyr::select(-"sum_of_element_standard_entropies_per_C")
  
}


#' Takes a data frame of molar elemental compositions of samples and estimates the standard enthalpy of formation per mol C
#' 
#' @export
irp_estimate_dhf0_per_C <- function(x, irp_m1, irp_d_compounds_standard_enthalpies_of_formation) {
  
  # estimate dhc0
  res <- 
    x |>
    irp_estimate_dhc0_per_C(irp_m1 = irp_m1)
  
  # assign combustion product standard enthalpies of formation to chemical elements and invert data frame
  irp_d_compounds_standard_enthalpies_of_formation <- 
    irp_d_compounds_standard_enthalpies_of_formation |>
    dplyr::mutate(
      chemical_element =
        dplyr::case_when(
          chemical_compound == "CO2" ~ "C",
          chemical_compound == "H2O" ~ "H",
          chemical_compound == "N2" ~ "N",
          chemical_compound == "SO2" ~ "S",
          chemical_compound == "O2" ~ "O",
          chemical_compound == "P4O10" ~ "P",
          chemical_compound == "K2O" ~ "K",
          chemical_compound == "MgO" ~ "Mg",
          chemical_compound == "CaO" ~ "Ca",
          chemical_compound == "Fe2O3" ~ "Fe",
          TRUE ~ NA_character_
        ),
      conversion_factor =
        dplyr::case_when(
          chemical_element == "C" ~ 1,
          chemical_element == "H" ~ 1/2,
          chemical_element == "N" ~ 1/2,
          chemical_element == "S" ~ 1,
          chemical_element == "O" ~ 1/2,
          chemical_element == "P" ~ 1/4,
          chemical_element == "K" ~ 1/2,
          chemical_element == "Mg" ~ 1,
          chemical_element == "Ca" ~ 1,
          chemical_element == "Fe" ~ 1/2,
          TRUE ~ NA_real_
        ) |>
        units::set_units("1", mode = "standard")
    )
  
  res_dhf0_combustion_products <- 
    irp_d_compounds_standard_enthalpies_of_formation |>
    dplyr::select(chemical_element, standard_enthalpy_of_formation) |>
    dplyr::select(standard_enthalpy_of_formation) |> 
    t() |>
    tibble::as_tibble() |>
    setNames(nm = irp_d_compounds_standard_enthalpies_of_formation$chemical_element) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::everything(),
        function(.x) .x * units::set_units(1, paste0("J/mol_", dplyr::cur_column()), mode = "standard")
      )
    )
  
  # get factors for conversion between amount of an element and amount of a combustion product
  res_chemical_elements_conversion_factors <- 
    irp_d_compounds_standard_enthalpies_of_formation |>
    dplyr::select(chemical_element, conversion_factor) |>
    dplyr::select(conversion_factor) |> 
    t() |>
    tibble::as_tibble() |>
    setNames(nm = irp_d_compounds_standard_enthalpies_of_formation$chemical_element) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::everything(),
        function(.x) units::set_units(.x, "1", mode = "standard")
      )
    )
  
  # extract chemical elements
  res_cf <- 
    x |>
    dplyr::select(dplyr::any_of(colnames(res_chemical_elements_conversion_factors))) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::everything(),
        function(.x) elco::elco_convert(.x, to = "mol")
      )
    )
  res_dhf0_combustion_products <- 
    res_dhf0_combustion_products |>
    dplyr::select(dplyr::all_of(colnames(res_cf)))
  res_dhf0_combustion_products <- 
    res_dhf0_combustion_products |>
    dplyr::select(match(colnames(res_cf), colnames(res_dhf0_combustion_products))) |>
    dplyr::slice(rep(1, nrow(res)))
  res_chemical_elements_conversion_factors <- 
    res_chemical_elements_conversion_factors |>
    dplyr::select(dplyr::all_of(colnames(res_cf)))
  res_chemical_elements_conversion_factors <-
    res_chemical_elements_conversion_factors |>
    dplyr::select(match(colnames(res_cf), colnames(res_chemical_elements_conversion_factors))) |>
    dplyr::slice(rep(1, nrow(res)))
  
  # compute sum of enthalpies of formation of combustion products
  res1 <- 
    (res_cf * res_dhf0_combustion_products * res_chemical_elements_conversion_factors) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      sum_dhf0 = sum(dplyr::c_across(dplyr::everything()), na.rm = TRUE)
    ) |>
    dplyr::ungroup()
  
  # eq (10) in @Popovic.2019
  new_unit <- as.character(units((res1$sum_dhf0[[1]] - units::set_units(1, as.character(units(posterior::draws_of(res$dhc0_per_C[[1]])[[1]] * res_cf$C[[1]])), mode = "standard")) / res_cf$C[[1]]))
  posterior::draws_of(res$dhc0_per_C) <- units::drop_units(posterior::draws_of(res$dhc0_per_C))
  
  res <- 
    res |>
    dplyr::mutate(
      dhf0_per_C = (units::drop_units(res1$sum_dhf0) - dhc0_per_C * units::drop_units(res_cf$C)) / units::drop_units(res_cf$C),
      dhf0_per_C = 
        {
          posterior::draws_of(dhf0_per_C) <- units::set_units(posterior::draws_of(dhf0_per_C), new_unit, mode = "standard") 
          dhf0_per_C
        }
    ) |>
    dplyr::select(-dhc0_per_C)
  
  res
  
}


#' Takes a data frame of molar elemental compositions of samples and estimates the standard Gibbs free energy of formation per mol C
#' 
#' @export
irp_estimate_dgf0_per_C <- function(x, irp_m1, irp_m2, irp_d_compounds_standard_enthalpies_of_formation) {
  
  # get enthalpy and entropy of formation
  res <- 
    x |>
    irp_estimate_dsf0_per_C(irp_m2 = irp_m2) |>
    irp_estimate_dhf0_per_C(
      irp_m1 = irp_m1, 
      irp_d_compounds_standard_enthalpies_of_formation = irp_d_compounds_standard_enthalpies_of_formation
    )
  
  # compute Gibbs free energy
  res <- 
    res |>
    dplyr::mutate(
      temperature =
        {
          res <- posterior::rvar(298.15)
          posterior::draws_of(res) <- units::set_units(posterior::draws_of(res), "K", mode = "standard")
          res
        },
      dgf0_per_C =
        dhf0_per_C - dsf0_per_C * temperature
    )
  
  x$dgf0_per_C <- res$dgf0_per_C
  x
  
}



#' Helper functions to make predictions with irp_m3
#' 
#' @param x Data frame with a column `bulk_density`.
#' 
#' @export
irp_predict_with_m3 <- function(x, irp_m3) {

  # make predictions
  res <- 
    x |>
    dplyr::bind_cols(
      {
        res <- 
          irp_m3$res_model |>
          brms::posterior_predict(
            type = "response",
            newdata = x
          ) |>
          posterior::rvar() |>
          tibble::as_tibble()
        colnames(res) <- paste0(colnames(res), "_1")
        res
      }
    ) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(c("volume_fraction_solids_1", "non_macroporosity_1", "macroporosity_1")),
        function(.x) irp_set_units_rvar(.x, "L/L", mode = "standard")
      )
    )
  
  res
  
}


#' Helper functions to make predictions with irp_m3
#' 
#' @param x Data frame with a column `bulk_density`.
#' 
#' @export
irp_predict_with_m4 <- function(x, irp_m4) {
  
  # make predictions
  res <- 
    x |>
    dplyr::mutate(
      saturated_hydraulic_conductivity_1 = 
        irp_m4$res_model |>
        brms::posterior_predict(
          type = "response",
          newdata = x
        ) |>
        posterior::rvar() |>
        magrittr::multiply_by(as.numeric(irp_m4$res_config$y_scale)) |>
        irp_set_units_rvar("cm/h", mode = "standard")
    )
  
  res
  
}


#' Makes predictions with model `irp_m5`
#' 
#' @param x A data frame with columns `N` (g/g) and `temperature` (in °C).
#' 
#' @param ... Arguments passed to `brms::posterior_predict()`.
#'
#' @export
irp_predict_with_m5 <- function(x, irp_m5, ...) {
  
  stopifnot(is.data.frame(x))
  
  # add predictors
  res <- 
    x |>
    dplyr::mutate(
      N =
        N |>
        magrittr::divide_by(irp_m5$res_config$N_scale),
      temperature =
        temperature |>
        magrittr::divide_by(irp_m5$res_config$temperature_scale)
    )
  
  # make predictions
  res <- 
    res |>
    dplyr::mutate(
      # standard enthalpy of combustion per C
      specific_heat_capacity_1 = 
        irp_m5$res_model |>
        brms::posterior_predict(
          type = "response",
          newdata = res,
          ...
        ) |>
        posterior::rvar() |>
        irp_set_units_rvar("J/K/g", mode = "standard")
    )
  
  x$specific_heat_capacity_1 <- res$specific_heat_capacity_1
  
  x
  
}


#' Makes predictions with model `irp_m6`
#' 
#' @param x A data frame with column `bulk_density` (g/cm$^{-3}$).
#' 
#' @param ... Arguments passed to `brms::posterior_predict()`.
#'
#' @export
irp_predict_with_m6 <- function(x, irp_m6, ...) {
  
  stopifnot(is.data.frame(x))
  
  # add predictors
  res <- x
  
  # make predictions
  res <- 
    res |>
    dplyr::mutate(
      # standard enthalpy of combustion per C
      dry_thermal_conductivity_1 = 
        irp_m6$res_model |>
        brms::posterior_predict(
          type = "response",
          newdata = res,
          ...
        ) |>
        posterior::rvar() |>
        irp_set_units_rvar("W/m/K", mode = "standard")
    )
  
  x$dry_thermal_conductivity_1 <- res$dry_thermal_conductivity_1
  
  x
  
}
