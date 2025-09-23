#### target variables ####

#' Defines the target variables
#' 
#' @export
irp_get_target_variables <- function() {
  
  c(
    "carbon_content",
    "nitrogen_content",
    "oxygen_content",
    "hydrogen_content",
    "phosphorus_content",
    "sulfur_content",
    "potassium_content",
    "titanium_content",
    "silicon_content",
    "calcium_content",
    "d13C",
    "d15N",
    "nosc",
    "dgf0",
    "loss_on_ignition",
    "bulk_density",
    "C_to_N",
    "O_to_C",
    "H_to_C"
  )
  
}


#### irp_d_model_info ####

#' Defines modeling options for each target variable
#' 
#' @param x One of `irp_target_variables`.
#' 
#' @export
irp_make_d_model_info <- function(n_preprocessing) {
  
  tibble::tibble(
    target_variable = irp_get_target_variables(),
    id_preprocessing = list(seq_len(n_preprocessing)),
    target_variable_label =
      dplyr::case_when(
        target_variable == "carbon_content" ~ "C",
        target_variable == "nitrogen_content" ~ "N",
        target_variable == "oxygen_content" ~ "O",
        target_variable == "hydrogen_content" ~ "H",
        target_variable == "phosphorus_content" ~ "P",
        target_variable == "sulfur_content" ~ "S",
        target_variable == "potassium_content" ~ "K",
        target_variable == "titanium_content" ~ "Ti",
        target_variable == "silicon_content" ~ "Si",
        target_variable == "calcium_content" ~ "Ca",
        target_variable == "d13C" ~ "d13C",
        target_variable == "d15N" ~ "d15N",
        target_variable == "nosc" ~ "nosc",
        target_variable == "dgf0" ~ "dgf0",
        target_variable == "loss_on_ignition" ~ "loss_on_ignition",
        target_variable == "bulk_density" ~ "bulk_density",
        target_variable == "C_to_N" ~ "C_to_N",
        target_variable == "O_to_C" ~ "O_to_C",
        target_variable == "H_to_C" ~ "H_to_C"
      )
  ) |>
    tidyr::unnest("id_preprocessing") |>
    dplyr::mutate(
      id_model = paste0(target_variable, "_", id_preprocessing)
    ) |>
    dplyr::relocate(id_model, .before = dplyr::everything())
  
}



#### Helper functions for `irp_d_model_info` #####

#' Adds seeds for random number generation to `irp_d_model_info`
#' 
#' @export
irp_add_rng_seed_to_d_model_info <- function(d_model_info) {
  
  d_model_info |>
    dplyr::mutate(
      rng_seed =
        dplyr::case_when(
          id_model == "bulk_density_1" ~ 324324,
          id_model == "bulk_density_2" ~ 324324,
          id_model == "bulk_density_3" ~ 777878,
          id_model == "O_to_C_2" ~ 533363680,
          id_model == "O_to_C_3" ~ 356324,
          TRUE ~ NA_real_
        )
    )
  
}



#' Adds information on which samples to use for training and which for testing to `irp_d_model_info`
#' 
#' @export
irp_add_data_partition_to_d_model_info <- function(d_model_info, irp_data_model_preprocessed) {
  
  d_model_info |>
    dplyr::mutate(
      data_partition =
        purrr::map(seq_along(id_model), function(i) {
          ir::ir_sample_prospectr(
            irp_data_model_preprocessed[[1L]] |> #---note: use the same preprocessed data for all models to allow model comparison
              dplyr::filter(id_measurement %in% id_measurement_all[[i]]) |>
              dplyr::select(id_measurement, spectra),
            sampling_function = prospectr::kenStone,
            metric = "euclid",
            k = n_training_sample[[i]],
            return_prospectr_output = FALSE
          ) |>
            ir::ir_drop_spectra()
        }) 
    )
  
}


#' Adds the training sample size to `irp_d_model_info`
#' 
#' @export
irp_add_n_training_sample_to_d_model_info <- function(d_model_info, irp_n_training_sample_max) {
  
  d_model_info |>
    dplyr::mutate(
      n_training_sample =
        purrr::map_int(seq_along(id_model), function(i) {
          min(floor(length(id_measurement_all[[i]]) * 0.8), irp_n_training_sample_max)
        }),
      n_training_sample = 
        dplyr::case_when(
          target_variable %in% c("C_to_N", "sulfur_content", "potassium_content", "d13C", "d15N", "carbon_content", "nitrogen_content", "phosphorus_content", "calcium_content", "bulk_density", "titanium_content") ~ 200L,
          TRUE ~ n_training_sample
        )
    )

}

#' Adds `id_measurement_all` to `irp_d_model_info`
#' 
#' @export
irp_add_id_measurement_all_to_d_model_info <- function(d_model_info, irp_data_model_preprocessed) {
  
  d_model_info |>
    dplyr::mutate(
      id_measurement_all =
        purrr::map(seq_along(target_variable), function(i) {
          irp_data_model_preprocessed[[id_preprocessing[[i]]]] |>
            dplyr::filter(eval(parse(text = filter_criteria[[i]]))) |>
            dplyr::pull(id_measurement)
        })
    )
  
}


#' Adds filter criteria to `irp_d_model_info`
#' 
#' @export
irp_add_filter_criteria_to_d_model_info <- function(d_model_info) {
  
  d_model_info |>
    dplyr::mutate(
      filter_criteria =
        dplyr::case_when(
          target_variable == "carbon_content" ~ '! is.na(C)',
          target_variable == "nitrogen_content" ~ '! is.na(N)',
          target_variable == "oxygen_content" ~ '! is.na(O) & sample_type != "dom"',
          target_variable == "hydrogen_content" ~ '! is.na(H) & sample_type != "dom"',
          target_variable == "phosphorus_content" ~ '! is.na(P)',
          target_variable == "sulfur_content" ~ '! is.na(S)',
          target_variable == "potassium_content" ~ '! is.na(K)',
          target_variable == "titanium_content" ~ '! is.na(Ti)',
          target_variable == "silicon_content" ~ '! is.na(Si)',
          target_variable == "calcium_content" ~ '! is.na(Ca)',
          target_variable == "d13C" ~ '! is.na(d13C)',
          target_variable == "d15N" ~ '! is.na(d15N)',
          target_variable == "nosc" ~ '! (is.na(C) | is.na(H) | is.na(O)) & sample_type != "dom"',
          target_variable == "dgf0" ~ '! (is.na(C) | is.na(H) | is.na(O)) & sample_type != "dom"',
          target_variable == "loss_on_ignition" ~ '! is.na(loss_on_ignition)',
          target_variable == "bulk_density" ~ '! is.na(bulk_density)',
          target_variable == "C_to_N" ~ '! is.na(N) & ! is.na(C) & (Ca < 5000 | is.na(Ca))',
          target_variable == "O_to_C" ~ '! is.na(O) & ! is.na(C) & sample_type != "dom"',
          target_variable == "H_to_C" ~ '! is.na(H) & ! is.na(C) & sample_type != "dom"'
        ) |>
        paste0(' & (mir_mode != "atr_ftir" | is.na(mir_mode)) & (! is_baseline_corrected | stringr::str_detect(core_label, "^peatbog")) & sample_type != "dom" & id_dataset != 16') 
    )
  
}


#' Adds information on variables with measurement error
#' 
#' 
#' @export
irp_add_y_has_error_to_d_model_info <- function(d_model_info) {
  
  d_model_info |>
    dplyr::mutate(
      y_has_error =
        target_variable %in% c("dgf0")
    )
  
}


#' Adds values for the target variable to `irp_d_model_info`.
#' 
#' @export
irp_add_y_to_d_model_info <- function(d_model_info, irp_data_model_preprocessed, irp_isotope_standards_isotope_fraction, irp_pmird_units, irp_m1, irp_m2, irp_d_compounds_standard_enthalpies_of_formation) {
  
  d_model_info |>
    dplyr::mutate(
      y =
        purrr::map(seq_along(target_variable), function(i) {
          switch(
            target_variable[[i]],
            "carbon_content" = ,
            "nitrogen_content" = ,
            "oxygen_content" = ,
            "hydrogen_content" = ,
            "loss_on_ignition" = {
              irp_data_model_preprocessed[[id_preprocessing[[i]]]] |>
                dplyr::filter(id_measurement %in% id_measurement_all[[i]]) |>
                dplyr::select(dplyr::all_of(target_variable_label[[i]])) |>
                setNames(nm = "y") |>
                dplyr::mutate(
                  y =
                    dplyr::case_when(
                      y <= 0.0 ~ 0.00001,
                      y >= 1.0 ~ 1 - 0.00001,
                      TRUE ~ y
                    )
                )
            },
            "bulk_density" = {
              irp_data_model_preprocessed[[id_preprocessing[[i]]]] |>
                dplyr::filter(id_measurement %in% id_measurement_all[[i]]) |>
                dplyr::select(dplyr::all_of(target_variable_label[[i]])) |>
                setNames(nm = "y")
            },
            "calcium_content" = ,
            "phosphorus_content" = ,
            "sulfur_content" = ,
            "potassium_content" = ,
            "titanium_content" = ,
            "silicon_content" = {
              irp_data_model_preprocessed[[id_preprocessing[[i]]]] |>
                dplyr::filter(id_measurement %in% id_measurement_all[[i]]) |>
                dplyr::select(dplyr::all_of(target_variable_label[[i]])) |>
                dplyr::mutate(
                  dplyr::across(
                    dplyr::all_of(target_variable_label[[i]]),
                    function(.x) {
                      ifelse(.x <= 0.0, 0.01, .x) |>
                      units::set_units("ug/g", mode = "standard") |>
                        units::set_units("g/g", mode = "standard") |>
                        units::drop_units()
                    }
                  )
                ) |>
                setNames(nm = "y")
            },
            "d13C" = {
              irp_data_model_preprocessed[[id_preprocessing[[i]]]] |>
                dplyr::filter(id_measurement %in% id_measurement_all[[i]]) |>
                dplyr::select(dplyr::all_of(target_variable_label[[i]])) |>
                #dplyr::mutate(
                #  d13C =
                #    d13C |>
                #    irp_delta_to_atom_fraction(r = irp_isotope_standards_isotope_fraction$r_13C)
                #) |>
                setNames(nm = "y")
            },
            "d15N" = {
              irp_data_model_preprocessed[[id_preprocessing[[i]]]] |>
                dplyr::filter(id_measurement %in% id_measurement_all[[i]]) |>
                dplyr::select(dplyr::all_of(target_variable_label[[i]])) |>
                #dplyr::mutate(
                #  d15N =
                #    d15N |>
                #    irp_delta_to_atom_fraction(r = irp_isotope_standards_isotope_fraction$r_15N)
                #) |>
                setNames(nm = "y")
            },
            "nosc" = {
              irp_data_model_preprocessed[[id_preprocessing[[i]]]] |>
                dplyr::filter(id_measurement %in% id_measurement_all[[i]]) |>
                dplyr::mutate(
                  dplyr::across(
                    dplyr::all_of(c("C", "H", "N", "O")), 
                    function(.x2) {
                      quantities::set_quantities(
                        .x2, 
                        errors = 
                          dplyr::cur_data() |> 
                          dplyr::pull(paste0(dplyr::cur_column(), "_err")),
                        unit = paste0("g_", dplyr::cur_column(), "/g_sample"), 
                        mode = "standard"
                      ) |>
                        magrittr::multiply_by(quantities::set_quantities(1, errors = 0, unit = "g_sample")) |>
                        elco::elco_convert(to = "mol")
                    }
                  ),
                  nosc_masiello2008 = 
                    elco::elco_nosc(C = C, H = H, N = N, O = O) |>
                    as.numeric(),
                ) |>
                dplyr::select(nosc_masiello2008) |>
                setNames(nm = "y")
              
            },
            "dgf0" = {
              
              irp_standard_entropies_elements <- 
                irp_make_standard_entropies_elements()
              
              res <- 
                irp_data_model_preprocessed[[id_preprocessing[[i]]]] |>
                dplyr::filter(id_measurement %in% id_measurement_all[[i]]) |>
                dplyr::mutate(
                  dplyr::across(
                    dplyr::any_of(c(irp_standard_entropies_elements$chemical_element, "Fe")), 
                    function(.x2) {
                      quantities::set_quantities(
                        ifelse(is.na(.x2), 0.0, .x2), 
                        errors = {
                          res <- 
                            dplyr::cur_data() |> 
                            dplyr::pull(paste0(dplyr::cur_column(), "_err")) 
                          ifelse(is.na(res), 0, res)
                        },
                        unit = 
                          irp_pmird_units |> 
                          dplyr::filter(attribute_name == dplyr::cur_column()) |> 
                          dplyr::pull(udunits_unit) |>
                          stringr::str_replace(pattern = "g/", replacement = paste0("g_", dplyr::cur_column(), "/")) |>
                          stringr::str_replace(pattern = "/g$", replacement = "/g_sample"), 
                        mode = "standard"
                      ) |>
                        magrittr::multiply_by(quantities::set_quantities(1, errors = 0, unit = "g_sample")) |>
                        elco::elco_convert(to = "mol") |>
                        errors::drop_errors()
                    }
                  )
                )
              
              res |>
                dplyr::mutate(
                  y = 
                    irp_estimate_dgf0_per_C(
                      x = as.data.frame(res),
                      irp_m1 = irp_m1, 
                      irp_m2 = irp_m2, 
                      irp_d_compounds_standard_enthalpies_of_formation = irp_d_compounds_standard_enthalpies_of_formation
                    ) |>
                    dplyr::pull(dgf0_per_C),
                ) |>
                dplyr::select(y) |>
                dplyr::mutate(
                  y = 
                    y |> 
                    irp_drop_units_rvar(),
                  y_err = posterior::sd(y),
                  y = mean(y)
                )
            },
            "C_to_N" = {
              irp_data_model_preprocessed[[id_preprocessing[[i]]]] |>
                dplyr::filter(id_measurement %in% id_measurement_all[[i]]) |>
                dplyr::mutate(
                  y = C/N,
                  #y = y / (1 + y)
                ) |>
                dplyr::select(y)
            },
            "O_to_C" = {
              irp_data_model_preprocessed[[id_preprocessing[[i]]]] |>
                dplyr::filter(id_measurement %in% id_measurement_all[[i]]) |>
                dplyr::mutate(y = O/C) |>
                dplyr::select(y)
            },
            "H_to_C" = {
              irp_data_model_preprocessed[[id_preprocessing[[i]]]] |>
                dplyr::filter(id_measurement %in% id_measurement_all[[i]]) |>
                dplyr::mutate(y = H/C) |>
                dplyr::select(y)
            },
          )
        })
        
    )
  
}


#' Adds the likelihood to use in models
#' 
#' @export
irp_add_likelihood_to_d_model_info <- function(d_model_info) {
  
  d_model_info |>
    dplyr::mutate(
      likelihood_name =
        dplyr::case_when(
          target_variable %in% c(paste0(c("carbon", "hydrogen", "oxygen", "nitrogen", "sulfur", "phosphorus", "potassium", "titanium", "silicon", "calcium"), "_content"), "nosc", "loss_on_ignition", "C_to_N") ~ "Beta",
          target_variable %in% c("bulk_density", "O_to_C", "H_to_C") ~ "Gamma",
          target_variable %in% c("dgf0", "d13C", "d15N") ~ "Gaussian",
          TRUE ~ NA_character_
        ),
      likelihood =
        dplyr::case_when(
          likelihood_name == "Beta" ~ list(brms::Beta(link = "logit", link_phi = "log")),
          likelihood_name == "Gamma" ~ list(Gamma(link = "log")),
          likelihood_name == "Gaussian" ~ list(gaussian(link = "identity")),
          TRUE ~ list(NULL)
        )
    )
  
}




#' Adds y scaling factors to `d_model_info`
#' 
#' @export
irp_add_y_scale_to_d_model_info <- function(d_model_info, irp_isotope_standards_isotope_fraction, irp_mass_density_fe) {
  
  d_model_info |>
    dplyr::mutate(
      y_center = 
        dplyr::case_when(
          target_variable %in% c("nosc") ~ -4, #---note: add 4 to scale values in [0, 1]
          target_variable %in% c("dgf0", "d13C", "d15N") ~ purrr::map_dbl(y, function(.y1) mean(.y1$y)),
          target_variable %in% c("carbon_content", "hydrogen_content", "oxygen_content", "nitrogen_content", "sulfur_content", "phosphorus_content", "potassium_content", "titanium_content", "calcium_content", "bulk_density", "O_to_C", "H_to_C", "C_to_N") ~ 0,
          #target_variable == "d15N" ~ irp_delta_to_atom_fraction(x = -500, r = irp_isotope_standards_isotope_fraction$r_15N), #---note: assuming that no peat sample has a d15N value smaller than -500.
          #target_variable == "d13C" ~ irp_delta_to_atom_fraction(x = -500, r = irp_isotope_standards_isotope_fraction$r_13C), #---note: assuming that no peat sample has a d13C value smaller than -500.
          TRUE ~ 0.0
        ),
      y_scale =
        dplyr::case_when(
          target_variable %in% c("nosc") ~ 8, #---note: divide by 8 to scale values in [0, 1]
          target_variable %in% c("dgf0", "d13C", "d15N") ~ purrr::map_dbl(y, function(.y1) sd(.y1$y)),
          target_variable %in% c("carbon_content", "hydrogen_content", "oxygen_content") ~ 1.0,
          target_variable %in% c("nitrogen_content", "sulfur_content", "calcium_content") ~ 0.1, #---note: assuming that no peat/DOM N, S content is ever > 10 mass-%. This will set the range [0, 10] mass-% to [0, 100] mass-% to facilitate MCMC sampling and incorporate additional prior information
          target_variable %in% c("phosphorus_content", "potassium_content", "titanium_content") ~ 0.05, #---note: assuming that no peat/DOM P, K content is ever > 5 mass-%. This will set the range [0, 5] mass-% to [0, 100] mass-% to facilitate MCMC sampling and incorporate additional prior information
          target_variable == "bulk_density" ~ 1.0,# irp_mass_density_fe, #---note: This is the mass density of pure Fe which can safely be assumed as upper bound. ---todo: I need a reliable source for this digit 
          target_variable == "O_to_C" ~ (4 * 16 )/ (1 * 12), #---note: Assumed maximum possible O/C mass ratio. Assumed based on a maximum molar ratio of 4/1 (4 O atoms per C atom) and then converted to a mass ratio.
          target_variable == "H_to_C" ~ (4 * 1)/(12), #---note: Assumed maximum possible H/C mass ratio. Assumed based on a maximum molar ratio of 4/1 (4 H atoms per C atom) and then converted to a mass ratio.
          target_variable == "C_to_N" ~ 0.7/0.001, # 0.05/0.3, # 0.7/0.001, #---note: Assumed maximum possible C/N mass ratio
          #target_variable == "d15N" ~ irp_delta_to_atom_fraction(x =  1000, r = irp_isotope_standards_isotope_fraction$r_15N), #---note: assuming that no peat sample has a d15N value larger than 1000.
          #target_variable == "d13C" ~ irp_delta_to_atom_fraction(x =  1000, r = irp_isotope_standards_isotope_fraction$r_13C), #---note: assuming that no peat sample has a d13C value larger than 1000.
          TRUE ~ 1.0
        )
    )
  
}


#' Adds `x_train` (training data as matrix, scaled) to `d_model_info`
#' 
#' @export
irp_add_x_train_to_d_model_info <- function(d_model_info, irp_data_model_preprocessed) {
  
  d_model_info |>
    dplyr::mutate(
      x_train =
        purrr::map(seq_along(target_variable), function(i) {
          irp_data_model_preprocessed[[id_preprocessing[[i]]]] |>
            dplyr::filter(id_measurement %in% id_measurement_all[[i]]) |>
            ir::ir_flatten() |>
            dplyr::select(-1) |>
            t() |>
            scale(center = TRUE, scale = TRUE)
        })
    )
  
}

#' Adds predictor scaling factors to `d_model_info`
#' 
#' @export
irp_add_predictor_scale_to_d_model_info <- function(d_model_info) {
  
  d_model_info |>
    dplyr::mutate(
      x_center = 
        x_train |>
        purrr::map(attr, which = "scaled:center"),
      x_scale = 
        x_train |>
        purrr::map(attr, which = "scaled:scale")
    )
  
}


#' Adds the training prediction domain to `d_model_info`
#' 
#' @export
irp_add_training_prediction_domain_to_d_model_info <- function(d_model_info, irp_data_model_preprocessed) {
  
  d_model_info |>
    dplyr::mutate(
      prediction_domain =
        purrr::map(seq_along(target_variable), function(i) {
          .x1 <- x_train[[i]]
          .y1 <- id_measurement_all[[i]]
            
            .x1 |>
              tibble::as_tibble() |>
              purrr::map(function(.x2) {
                tibble::tibble(
                  ymin = min(.x2), 
                  ymax = max(.x2), 
                  ymin_id = .y1[[which.min(.x2)[[1]]]], 
                  ymax_id = .y1[[which.max(.x2)[[1]]]]
                )
              }) |>
              dplyr::bind_rows() |>
              dplyr::mutate(
                x = 
                  irp_data_model_preprocessed[[id_preprocessing[[i]]]]$spectra[[1]]$x |> 
                  stringr::str_extract_all(pattern = "\\d+") |> 
                  as.numeric()
              ) |>
              dplyr::relocate(x, .before = "ymin") |>
              irpeat::new_irp_prediction_domain()
        })
      
    )
  
}


#' Adds the number of components to compute during dimension reduction
#' 
#' @export
irp_add_dimension_reduction_ncomps_to_d_model_info <- function(d_model_info) {
  
  d_model_info |>
    dplyr::mutate(
      ncomps = 5L
    )
  
}


#' Adds information on priors to use in the models
#' 
#' @export
irp_add_priors_to_d_model_info <- function(d_model_info, irp_isotope_standards_isotope_fraction) {
  
  d_model_info |>
    dplyr::mutate(
      prior_phi = 
        dplyr::case_when(
          target_variable == "dgf0" ~ "gamma(prior_phi_shape, prior_phi_scale)",
          TRUE ~ "gamma(prior_phi_shape, prior_phi_scale)"
        ),
      prior_phi_shape = 
        dplyr::case_when(
          target_variable %in% c("d13C", "d15N") ~ 5.0,
          target_variable %in% c("dgf0") ~ 10.0,
          #target_variable == "carbon_content" ~ 50.0,
          TRUE ~ 5
        ),
      prior_phi_scale = 
        dplyr::case_when(
          target_variable %in% c("d13C", "d15N") ~ 5.0/0.5, #---note: zero-avoiding prior
          target_variable %in% c("dgf0") ~ 10.0/0.4,
          #target_variable == "carbon_content" ~ 50.0 / 300.0,
          target_variable %in% c("titanium_content", "potassium_content", "calcium_content", "silicon_content", "sulfur_content") ~ 5.0 / 100.0,
          target_variable %in% c("carbon_content", "nitrogen_content", "C_to_N") ~ 5.0 / 100.0,
          TRUE ~ 5.0/300
        ),
      prior_phi_class =
        purrr::map_chr(seq_along(target_variable), function(i) {
          switch(
            likelihood_name[[i]],
            "Gaussian" = "sigma",
            "Gamma" = "shape",
            "Beta" = "phi",
            stop("Wrong `likelihood_name` while defining `prior_phi_class`!")
          )
        }),
      prior_intercept = "normal(prior_intercept_mu, prior_intercept_sd)",
      prior_intercept_mu = 
        dplyr::case_when(
          target_variable == "carbon_content" ~ brms::logit_scaled(0.47), # see Loisel.2014
          target_variable == "nitrogen_content" ~ -2, # see Loisel.2014
          target_variable == "oxygen_content" ~ brms::logit_scaled(0.34), # see Moore.2018a
          target_variable == "hydrogen_content" ~ brms::logit_scaled(0.052), # see Moore.2018a
          target_variable == "phosphorus_content" ~ brms::logit_scaled(0.0005/y_scale), #---note The average value (0.001) is from Moore.2018a
          target_variable == "sulfur_content" ~ brms::logit_scaled(0.002/y_scale), #---note: See Moore.2018a
          target_variable == "potassium_content" ~ brms::logit_scaled(0.001/y_scale), #---note The average value (0.001) is from Moore.2018a
          target_variable == "titanium_content" ~ brms::logit_scaled(0.001/y_scale), #---note: assuming the same content as for P
          target_variable == "silicon_content" ~ brms::logit_scaled(0.05/y_scale),
          target_variable == "calcium_content" ~ brms::logit_scaled(0.002/y_scale), #---note: assuming the same content as for P
          target_variable == "d13C" ~ 0.0,# brms::logit_scaled((irp_delta_to_atom_fraction(x = -25, r = irp_isotope_standards_isotope_fraction$r_13C) - y_center)/y_scale), # Approximate average d13C value for terrestrial C3 plants from Kohn.2010.
          target_variable == "d15N" ~ 0.0,# brms::logit_scaled((irp_delta_to_atom_fraction(x = 0, r = irp_isotope_standards_isotope_fraction$r_15N) - y_center)/y_scale), # Isotope value of the AIR standard
          target_variable == "nosc" ~ brms::logit_scaled((-0.4 + 4.0)/y_scale), # Corresponds roughly to the average value for deeper peat given in Moore.2018
          target_variable == "dgf0" ~ 0.0,
          target_variable == "loss_on_ignition" ~ brms::logit_scaled(0.97/y_scale), # see Loisel.2014
          target_variable == "bulk_density" ~ log(0.076/y_scale), # see Loisel.2014
          target_variable == "C_to_N" ~ brms::logit_scaled(55/y_scale), # brms::logit_scaled((55) / (1 + (55))), # log((1/55)/y_scale), # brms::logit_scaled(55/y_scale), # see Loisel.2014
          target_variable == "O_to_C" ~ units::set_units(0.65, "mol_O/mol_C") |> units::set_units("g_O/g_C") |> as.numeric() |> magrittr::divide_by(y_scale) |> log(), # units::set_units(0.65, "mol_O/mol_C") |> units::set_units("g_O/g_C") |> as.numeric() |> magrittr::divide_by(y_scale) |> brms::logit_scaled(), # see values cited for various peat materials in Moore.2018a (similar to bog peat)
          target_variable == "H_to_C" ~ units::set_units(1.4, "mol_H/mol_C") |> units::set_units("g_H/g_C") |> as.numeric() |> magrittr::divide_by(y_scale) |> log(), # units::set_units(1.4, "mol_H/mol_C") |> units::set_units("g_H/g_C") |> as.numeric() |> magrittr::divide_by(y_scale) |> brms::logit_scaled(), # see values cited for various peat materials in Moore.2018a and Leifeld.2020
          TRUE ~ NA_real_
        ),
      prior_intercept_sd = 
        dplyr::case_when(
          target_variable == "carbon_content" ~ 0.5, 
          target_variable == "nitrogen_content" ~ 0.2,
          target_variable == "oxygen_content" ~ 0.2,
          target_variable == "hydrogen_content" ~ 0.2,
          target_variable == "phosphorus_content" ~ 0.2,
          target_variable == "sulfur_content" ~ 0.5,
          target_variable == "potassium_content" ~ 0.1,
          target_variable == "titanium_content" ~ 0.5,
          target_variable == "silicon_content" ~ 0.1,
          target_variable == "calcium_content" ~ 0.2,
          target_variable == "d13C" ~ 0.5,
          target_variable == "d15N" ~ 0.5,
          target_variable == "nosc" ~ 0.2,
          target_variable == "dgf0" ~ 0.05,
          target_variable == "loss_on_ignition" ~ 0.5,
          target_variable == "bulk_density" ~ 0.2,
          target_variable == "C_to_N" ~ 0.2,
          target_variable == "O_to_C" ~ 0.2,
          target_variable == "H_to_C" ~ 0.2,
          TRUE ~ NA_real_
        ),
      prior_b_par_ratio =
        purrr::map(seq_along(target_variable), function(i) {
          p <- ncol(x_train[[i]])
          p0 <- 
            dplyr::case_when(
              target_variable[[i]] == "carbon_content" ~ 5,
              TRUE ~ 5
            )
          p0/p
        }),
      prior_b_df = 
        dplyr::case_when(
          target_variable %in% c("oxygen_content", "hydrogen_content", "phosphorus_content", "titanium_content", "C_to_N", "O_to_C", "dgf0_1") ~ 4,
          TRUE ~ 3
        ),
      prior_b = paste0("horseshoe(df = ", prior_b_df, ", par_ratio = ", prior_b_par_ratio, ", df_global = 1, autoscale = TRUE)"),
      prior_all =
        purrr::map(seq_along(target_variable), function(i) {
          rbind(
            brms::prior_string(prior_intercept[[i]], class = "Intercept"),
            brms::prior_string(prior_phi[[i]], class = prior_phi_class[[i]], lb = 0.0),
            brms::prior_string(prior_b[[i]], class = "b")
          ) 
        }),
      prior_stanvars =
        purrr::map(seq_along(target_variable), function(i) {
          brms::stanvar(prior_phi_shape[[i]], name = "prior_phi_shape") +
            brms::stanvar(prior_phi_scale[[i]], name = "prior_phi_scale") +
            brms::stanvar(prior_intercept_mu[[i]], name = "prior_intercept_mu") +
            brms::stanvar(prior_intercept_sd[[i]], name = "prior_intercept_sd")
        })
    )
  
}


#' Adds the brms model formula for the models
#' 
#' @export
irp_add_model_formula_to_d_model_info <- function(d_model_info) {
  
  d_model_info |>
    dplyr::mutate(
      brms_formula = 
        purrr::map(seq_along(target_variable), function(i) {
        if(y_has_error[[i]]) { #---todo: need still an own formula for EAC/EDC
          brms::bf(y | mi(y_err) ~ .)
        } else {
          brms::bf(y ~ .)
        }
      })
    )
  
}


#' Adds `id_model` to `irp_d_model_info`
#' 
#' @export
irp_add_id_model_to_d_model_info <- function(d_model_info) {
  
  d_model_info |>
    dplyr::group_split(target_variable) |>
    purrr::map(function(.x) {
      .x |>
        dplyr::mutate(
          id_model = paste0(target_variable, "_", seq_along(target_variable))
        )
    }) |>
    dplyr::bind_rows() |>
    dplyr::relocate(id_model, .before = dplyr::everything())
  
}



#' Add validation fold id to `irp_d_model_info`
#' 
#' @export
irp_add_id_validation_fold_to_d_model_info <- function(irp_d_model_info) {
  
  irp_d_model_info |>
    dplyr::left_join(
      irp_d_model_info |>
        dplyr::group_by(target_variable, uses_only_high_quality_mirs, validation_mode) |>
        dplyr::summarise(
          .groups = "drop"
        ) |>
        dplyr::mutate(
          id_validation_fold = seq_along(target_variable)
        ),
      by = c("target_variable", "uses_only_high_quality_mirs", "validation_mode")
    )
  
}




#### Dimension reduction ####


#' Settings for dimension reduction methods
#' 
#' @export
irp_get_dimension_reduction_model_settings_1 <- function() {
  
  list(
    "plsr" =
      list(
        scale = FALSE, 
        center = FALSE,
        validation = "none"
      ),
    "ispca" =
      list(
        normalize = FALSE,
        preprocess = FALSE,
        permtest = TRUE,
        alpha = 0.01,
        perms = 1e3,
        sup.only = TRUE
      )
  )
  
}


#' Helper that prepares data for dimension reduction models
irp_prepare_data_for_dimension_reduction_models <- function(x, irp_data_model_preprocessed) {
  
  x <- 
    x |>
    irp_add_x_train_to_d_model_info(irp_data_model_preprocessed = irp_data_model_preprocessed)
  
  tibble::tibble(
    y = 
      x$y[[1]]$y |>
      scale(center = x$y_center[[1]], scale = x$y_scale[[1]]) |>
      x$likelihood[[1]]$linkfun() |>
      as.numeric(),
    x = 
      x$x_train[[1]] |>
      tibble::as_tibble() |>
      setNames(nm = paste0("V", irp_data_model_preprocessed[[x$id_preprocessing[[1]]]]$spectra[[1]]$x))
  )
  
}


#' Computes a PLSR model for a target variable
#' 
#' @param x One row in `irp_d_model_info`.
#' 
#' @export
irp_compute_plsr_1 <- function(x, irp_data_model_preprocessed, irp_dimension_reduction_model_settings_1) {
  
  stopifnot(nrow(x) == 1 && x$do_dimension_reduction && x$dimension_reduction_method == "plsr")
  
  res_data <- 
    irp_prepare_data_for_dimension_reduction_models(
      x = x,
      irp_data_model_preprocessed = irp_data_model_preprocessed
    )
  
  res_settings <- irp_dimension_reduction_model_settings_1$plsr
  
  pls::plsr(
    y ~ as.matrix(x),
    data = res_data,
    ncomp = x$ncomps[[1]],
    scale = res_settings$scale, 
    center = res_settings$center,
    validation = res_settings$validation
  )
  
}


#' Computes an ISPCA model for a target variable
#' 
#' @param x One row in `irp_d_model_info`.
#' 
#' @export
irp_compute_ispca_1 <- function(x, irp_data_model_preprocessed, irp_dimension_reduction_model_settings_1) {
  
  stopifnot(nrow(x) == 1 && x$do_dimension_reduction && x$dimension_reduction_method == "ispca")
  
  res_data <- 
    irp_prepare_data_for_dimension_reduction_models(
      x = x,
      irp_data_model_preprocessed = irp_data_model_preprocessed
    )
  
  res_settings <- irp_dimension_reduction_model_settings_1$ispca
  
  dimreduce::ispca(
    x = res_data$x,
    y =  res_data$y,
    nctot = x$ncomps[[1]],
    ncsup = x$ncomps[[1]],
    normalize = res_settings$normalize,
    preprocess = res_settings$preprocess,
    permtest = res_settings$permtest,
    alpha = res_settings$alpha,
    perms = res_settings$perms,
    sup.only = res_settings$sup.only
  )
  
}


#' Helper function that decides what dimension reduction method to use, based on `irp_d_model_info$dimension_reduction_method`
#' 
#' @param x A row in `irp_d_model_info` where `do_dimension_reduction` 
#' is `TRUE`.
#' 
#' @export
irp_compute_dimension_reduction_model_1 <- function(x, irp_data_model_preprocessed, irp_dimension_reduction_model_settings_1) {
  
  stopifnot(nrow(x) == 1 && x$do_dimension_reduction)
  
  res <- 
    switch(
      x$dimension_reduction_method,
      "plsr" = 
        irp_compute_plsr_1(
          x = x, 
          irp_data_model_preprocessed = irp_data_model_preprocessed, 
          irp_dimension_reduction_model_settings_1 = irp_dimension_reduction_model_settings_1
        ),
      "ispca" = 
        irp_compute_ispca_1(
          x = x, 
          irp_data_model_preprocessed = irp_data_model_preprocessed, 
          irp_dimension_reduction_model_settings_1 = irp_dimension_reduction_model_settings_1
        ),
      stop("Unknown `x$dimension_reduction_method`.")
    )
  
  res
  
}


#' Helper that extracts scores of dimension reduction models as data frame to use as inputs for a brms model
#' 
#' @param ncomps Integer. Number of components for which to predict scores.
#' 
#' @param irp_dimension_reduction_model_1 An object of class `mvr` or `ispca`.
#' 
#' @param newdata Either a data frame as returned by 
#' `irp_prepare_data_for_dimension_reduction_models()`, or `NULL`, in which case
#' scores will be extracted for the training data used to compute the dimension
#' reduction model.
#'
#' @export
irp_predict_scores_dimension_reduction_model_1 <- function(irp_dimension_reduction_model_1, newdata = NULL, ncomps) {
  
  res_model <- irp_dimension_reduction_model_1
  
  if(inherits(res_model, "mvr")) {
    
    res_scores_olddata_sd <- 
      stats::predict(
        res_model,
        ncomp = seq_len(ncomps),
        type = "scores"
      ) |>
      tibble::as_tibble() |>
      purrr::map_dbl(stats::sd)
    
    res_scores_newdata <- 
      stats::predict(
        res_model,
        ncomp = seq_len(ncomps),
        type = "scores", 
        newdata = newdata 
      ) |>
      tibble::as_tibble()
    
  } else if(inherits(res_model, "ispca")) {
    
    res_scores_olddata_sd <- 
      res_model$z |>
      tibble::as_tibble() |>
      dplyr::slice(seq_len(ncomps)) |>
      purrr::map_dbl(stats::sd)
    
    # need to do this manually because predict.dimred does not handle newdata = NULL
    res_scores_newdata <- 
      if(is.null(newdata)) {
        res_model$z
      } else {
        stats::predict(
          res_model,
          xnew = newdata$x 
        ) 
      }   
    
    res_scores_newdata <-
      res_scores_newdata |>
      tibble::as_tibble() |>
      dplyr::select(seq_len(ncomps))
    
  }
  
  res <- 
    res_scores_newdata |>
    setNames(nm = paste0("V", seq_len(ncol(res_scores_newdata)))) |>
    purrr::map2(res_scores_olddata_sd, function(.x, .y) {
      .x/.y
    }) |>
    dplyr::bind_cols() |>
    as.matrix()
  
  res
  
}




#### brms models ####


#' Helper that prepares data for the brms models using scores from dimension reduction models as predictors
#' 
#' @export
irp_prepare_data_for_brms_models_1 <- function(x) {
  
  d_model_info <- x
  
  res <- 
    tibble::tibble(
      y = 
        d_model_info$y[[1]] |>
        dplyr::filter(d_model_info$id_measurement_all[[1]] %in% (d_model_info$data_partition[[1]] |> dplyr::filter(for_prospectr_model) |> dplyr::pull(id_measurement))) |>
        dplyr::pull(y) |>
        scale(center = d_model_info$y_center[[1]], scale = d_model_info$y_scale[[1]]) |>
        as.numeric(),
      x = 
        d_model_info$x_train[[1]] |>
        tibble::as_tibble() |>
        dplyr::filter(d_model_info$id_measurement_all[[1]] %in% (d_model_info$data_partition[[1]] |> dplyr::filter(for_prospectr_model) |> dplyr::pull(id_measurement))) |>
        as.matrix()
    )
  
  if(d_model_info$y_has_error[[1]]) {
    res <- 
      res |>
      dplyr::mutate(
        y_err =
          d_model_info$y[[1]] |>
          dplyr::filter(d_model_info$id_measurement_all[[1]] %in% (d_model_info$data_partition[[1]] |> dplyr::filter(for_prospectr_model) |> dplyr::pull(id_measurement))) |>
          dplyr::pull(y_err) |>
          scale(center = 0.0, scale = d_model_info$y_scale[[1]]) |>
          as.numeric()
      ) |>
      dplyr::relocate(
        y_err, .after = "y"
      )
  }
  
  res
  
}


#' Compiles a brmsfit model (to avoid recompilation with different datasets)
#' 
#' @export
irp_make_brms_compiled_model <- function(irp_brms_compiled_model_metadata) {
  
  res_d_model_info <- 
    irp_brms_compiled_model_metadata
  
  res_d_model_info$brms_model <- 
    purrr::map(seq_len(nrow(res_d_model_info)), function(i) {
      
      brms::brm(
        formula = res_d_model_info$brms_formula[[i]],
        data = 
          tibble::tibble(
            y = 0.5, 
            y_err = 
              if(res_d_model_info$y_has_error[[i]]) {
                0.1
              } else {
                NULL
              }, 
            x = 3
          ), 
        family = res_d_model_info$likelihood[[i]],
        prior = res_d_model_info$prior_all[[i]],
        stanvars = res_d_model_info$prior_stanvars[[i]],
        chains = 0,
        save_pars = brms::save_pars(all = TRUE),
        backend = "cmdstanr",
        save_warmup = TRUE,
        sig_figs = 14L
      )
      
    })
  
  res_d_model_info |>
    dplyr::select(! dplyr::starts_with("prior_"))
    
}


#' Helper to fit brms models
#' 
#' @param x A row in `irp_d_model_info_byrow`
#' 
#' @export
irp_make_fit_1 <- function(x, irp_brms_compiled_model, irp_data_for_brms_models_1, irp_isotope_standards_isotope_fraction, irp_mcmc_settings, irp_fit_1_pathfinder) {
  
  stopifnot(nrow(x) == 1)
  
  res_brms_compiled <- 
    irp_brms_compiled_model |>
    dplyr::filter(likelihood_name == x$likelihood_name & y_has_error == x$y_has_error) |>
    dplyr::pull(brms_model)
  
  update(
    res_brms_compiled[[1]], 
    newdata = irp_data_for_brms_models_1,
    stanvars = x$prior_stanvars[[1]],
    prior = x$prior_all[[1]],
    iter = irp_mcmc_settings$iter,
    warmup = irp_mcmc_settings$warmup,
    chains = irp_mcmc_settings$chains,
    cores = 4L,
    control = 
      list(
        max_treedepth = irp_mcmc_settings$max_treedepth,
        adapt_delta = irp_mcmc_settings$adapt_delta
      ),
    init = {
      d_pathfinder <- 
        as.data.frame(irp_fit_1_pathfinder) |>
        dplyr::select(! dplyr::all_of(c("lprior", "lp__", "lp_approx__")))
      d_cuts_width <- nrow(d_pathfinder) / irp_mcmc_settings$chains
      d_cuts <- 
        tibble::tibble(
          start = seq(1, by = d_cuts_width, length.out = irp_mcmc_settings$chains),
          end = start + d_cuts_width - 1
        )
      purrr::map(seq_len(irp_mcmc_settings$chains), function(i) {
        res <- 
          tibble::tibble(
          value = apply(d_pathfinder[d_cuts$start[[i]]:d_cuts$end[[i]], ], 2, mean),
          name = 
            names(value) |>
            stringr::str_remove(pattern = "\\[\\d+\\]$")
        ) |>
          tidyr::nest(data = c(value))
        res1 <- purrr::map(res$data, function(.x) unname(.x$value))
        names(res1) <- res$name
        res1
      })
    },
    save_pars = brms::save_pars(all = TRUE),
    backend = "cmdstanr",
    save_warmup = TRUE,
    sig_figs = 14L
  )
  
}


#' Helper to fit brms pathfinder models
#' 
#' @param x A row in `irp_d_model_info_byrow`
#' 
#' @export
irp_make_fit_1_pathfinder <- function(x, irp_brms_compiled_model, irp_data_for_brms_models_1, irp_isotope_standards_isotope_fraction, irp_mcmc_settings) {
  
  stopifnot(nrow(x) == 1)
  #switch(
  #  x$id_model,
  #  "bulk_density_1" = 324324,
  #  "O_to_C_2" = 533363680,
  #  NA
  #) |>
  #  set.seed()
  
  res_brms_compiled <- 
    irp_brms_compiled_model |>
    dplyr::filter(likelihood_name == x$likelihood_name & y_has_error == x$y_has_error) |>
    dplyr::pull(brms_model)
  
  update(
    res_brms_compiled[[1]], 
    newdata = irp_data_for_brms_models_1,
    stanvars = x$prior_stanvars[[1]],
    algorithm = "pathfinder",
    seed = x$rng_seed[[1]],
    init =
      if(x$likelihood_name[[1]] == "Gamma") {
        NULL
      } else {
        purrr::map(seq_len(8L), function(i) {
          res <- list(Intercept = rnorm(1L, x$prior_intercept_mu[[1]], x$prior_intercept_sd[[1]]))
          res_phi_name <- 
            switch(
              x$likelihood_name[[1]],
              "Beta" = "phi",
              "Gamma" = "shape",
              "Gaussian" = "sigma"
            )
          res_phi <- 
            switch(
              x$likelihood_name[[1]],
              "Beta" = ,
              "Gamma" = ,
              "Gaussian" = rgamma(1L, x$prior_phi_shape[[1]], x$prior_phi_scale[[1]])
            )
          res[[res_phi_name]] <- res_phi
          res
        }) 
      },
    num_paths = 8L,
    history_size = 100,
    max_lbfgs_iters = 500,
    single_path_draws = 40,
    draws = 160,
    psis_resample = FALSE,
    sig_figs = 10
  )
  
}



#### Model validation ####

#### new:start ####

#' Compares ELPD between compatible models
#' 
#' @export
irp_make_fit_1_map_elpd_compare <- function(irp_d_model_info_enriched_1, irp_fit_1_map_elpd) {
  
  irp_d_model_info_enriched_1 |>
    dplyr::mutate(
      elpd = irp_fit_1_map_elpd
    ) |>
    dplyr::group_by(target_variable) |>
    dplyr::summarise(
      elpd_compare = 
        {
          res <-
            loo::loo_compare(x = elpd)
          res |>
            as.data.frame() |>
            dplyr::mutate(
              id_model = 
                rownames(res) |>
                stringr::str_remove(pattern = "^irp_fit_1_elpd_")
            )|>
            list()
        },
      .groups = "drop"
    ) |>
    tidyr::unnest(c("elpd_compare")) |>
    dplyr::relocate("id_model", .before = dplyr::everything())
  
}

#' Collects measured values and predictions for all models
#' 
#' @export
irp_fit_1_make_evaluation_1 <- function(irp_d_model_info_enriched_2, irp_fit_1, irp_isotope_standards_isotope_fraction) {
  
  d_model_info <- irp_d_model_info_enriched_2
  
  # prepare all data for prediction
  res_newdata <- 
    tibble::tibble(
      y = 
        d_model_info$y[[1]] |>
        dplyr::pull(y) |>
        scale(center = d_model_info$y_center[[1]], scale = d_model_info$y_scale[[1]]) |>
        as.numeric(),
      y_err =
        if(d_model_info$y_has_error) {
          d_model_info$y[[1]] |>
            dplyr::pull(y_err) |>
            scale(center = 0.0, scale = d_model_info$y_scale[[1]]) |>
            as.numeric()
        } else {
          NULL
        },
      x = 
        d_model_info$x_train[[1]] |>
        tibble::as_tibble() |>
        as.matrix()
    )
  
  # predictions
  res <- 
    tibble::tibble(
      id_model = d_model_info$id_model,
      target_variable = d_model_info$target_variable,
      id_measurement = d_model_info$id_measurement_all[[1]],
      is_training_data = d_model_info$data_partition[[1]]$for_prospectr_model,
      y = res_newdata$y,
      y_err = 
        if(is.null(res_newdata$y_err)) {
          NA_real_
        } else {
          res_newdata$y_err * d_model_info$y_scale[[1]]
        },
      yhat = 
        brms::posterior_predict(
          irp_fit_1,
          re_formula = NULL,
          newdata = res_newdata
        ) |>
        posterior::rvar()
    ) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(c("y", "yhat")),
        function(.x) {
          .x |>
            magrittr::multiply_by(d_model_info$y_scale[[1]]) |>
            magrittr::add(d_model_info$y_center[[1]])
        }
      ),
      unit_for_modeling =
        irp_get_units_target_variables(target_variable, what = "unit_for_modeling"),
      unit_for_plotting = 
        irp_get_units_target_variables(target_variable, what = "unit_for_plotting"),
      y = 
        y |>
        units::set_units(unit_for_modeling[[1]], mode = "standard") |>
        units::set_units(unit_for_plotting[[1]], mode = "standard") |>
        units::drop_units(),
      y_err = 
        y_err |>
        units::set_units(unit_for_modeling[[1]], mode = "standard") |>
        units::set_units(unit_for_plotting[[1]], mode = "standard") |>
        units::drop_units(),
      yhat = 
        yhat |>
        irp_set_units_rvar(unit_for_modeling[[1]], mode = "standard") |>
        irp_set_units_rvar(unit_for_plotting[[1]], mode = "standard") |>
        irp_drop_units_rvar(),
      y_mcse_mean = posterior::mcse_mean(yhat),
      y_mcse_sd = posterior::mcse_sd(yhat),
      y_mcse_lower = posterior::mcse_quantile(yhat, probs = 0.025),
      y_mcse_upper = posterior::mcse_quantile(yhat, probs = 0.975)
    ) |>
    dplyr::select(-target_variable)
  
  # compute rmse
  res_rmse <- 
    res |>
    dplyr::group_by(is_training_data) |>
    dplyr::summarise(
      rmse = irp_rmse_rvar(yhat = yhat, y = y),
      bias = posterior::rvar_mean(y - yhat),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      rmse_train_minus_test = rmse[is_training_data] - rmse[! is_training_data]
    )
  
  res <- 
    list(
      yhat = 
        res |>
        dplyr::mutate(
          yhat_err = posterior::sd(yhat),
          yhat = mean(yhat)
        ),
      rmse = 
        res_rmse |>
        dplyr::mutate(
          dplyr::across(
            dplyr::where(posterior::is_rvar),
            function(x) {
              purrr::map(x, function(.x) {
                res <- 
                  tibble::tibble(
                  mean = mean(.x),
                  lower = posterior::quantile2(.x, probs = 0.025),
                  upper = posterior::quantile2(.x, probs = 0.975)
                )
                res |>
                  setNames(nm = paste0(dplyr::cur_column(), "_", colnames(res)))
              })
            } 
          )
        ) |>
        tidyr::unnest(c("rmse", "bias", "rmse_train_minus_test")),
      num_divergent = 
        irp_fit_1$fit |>
        rstan::get_num_divergent(),
      max_rhat =
        irp_fit_1 |>
        brms::rhat() |>
        max()
    )
  
}


#' Computes the training and testing prediction domains for a model
#' 
#' @export
irp_make_fit_1_pd <- function(irp_d_model_info_enriched_2, irp_data_model_preprocessed) {
  
  res_ir <- 
    irp_data_model_preprocessed[[irp_d_model_info_enriched_2$id_preprocessing[[1]]]] |>
    ir::ir_scale(
      center = irp_d_model_info_enriched_2$x_center[[1]],
      scale = irp_d_model_info_enriched_2$x_scale[[1]]
    )
  
  list(
    train = 
      res_ir |>
      dplyr::filter(id_measurement %in% (irp_d_model_info_enriched_2$data_partition[[1]] |> dplyr::filter(for_prospectr_model) |> dplyr::pull(id_measurement))) |>
      irpeat::irp_as_irp_prediction_domain(),
    test = 
      res_ir |>
      dplyr::filter(id_measurement %in% (irp_d_model_info_enriched_2$data_partition[[1]] |> dplyr::filter(for_prospectr_test) |> dplyr::pull(id_measurement))) |>
      irpeat::irp_as_irp_prediction_domain(),
    all = 
      res_ir |>
      dplyr::filter( (mir_mode != "atr_ftir" | is.na(mir_mode)) & (! is_baseline_corrected | stringr::str_detect(core_label, "^peatbog")) & sample_type != "dom") |>
      dplyr::filter(id_dataset != 16) |> #---note: seem to be ATR spectra
      irpeat::irp_as_irp_prediction_domain()
  )
  
}


#' Extracts model coefficients (slopes)
#' 
#' @export
irp_make_fit_1_slopes <- function(irp_fit_1, irp_d_model_info_enriched_2, irp_data_model_preprocessed) {
  
  res_model <- irp_fit_1
  res_model_info <- irp_d_model_info_enriched_2
  
  res <- 
    tibble::tibble(
      x = 
        irp_data_model_preprocessed[[res_model_info$id_preprocessing]]$spectra[[1]]$x,
      slope =
        res_model |>
        as.data.frame() |>
        dplyr::select(dplyr::starts_with("b_x")) |>
        as.matrix() |>
        posterior::rvar()
    )
  
  res_file <- paste0("targets_rvars/irp_fit_1_", irp_d_model_info_enriched_2$id_model, ".rds")
  saveRDS_rvars(res, res_file)
  res_file
  
}



#### new:end ####


#' Estimates depth ranges for which selected peat properties are autocorrelated with depth
#' 
#' @param irp_data_model_preprocessed One element of `irp_data_model_preprocessed`.
#' 
#' @export
irp_estimate_acf_depth_range <- function(irp_data_model_preprocessed) {
  
  # settings
  irp_settings_cv_depth_acf_analysis <- 
    list(
      layer_thickness_max = 5,
      layer_gap_max = 10,
      core_n_samples_min = 30,
      acf_confidence_level = 0.95
    )
  
  # define the target variables
  cv_fold_temporal_target_variables <- rlang::quos(C, N, d13C, d15N)
  
  res <- 
    purrr::map(seq_along(cv_fold_temporal_target_variables), function(i) {
      
      irp_data_model_preprocessed |>
        dplyr::filter(sample_type == "peat" & ! is.na(sample_depth_upper) & ! is.na(sample_depth_lower) & (sample_depth_lower - sample_depth_upper) <= irp_settings_cv_depth_acf_analysis$layer_thickness_max & ! is.na(!!cv_fold_temporal_target_variables[[i]])) |>
        dplyr::arrange(id_dataset, core_label, sample_depth_upper) |>
        dplyr::mutate(
          .core_label = paste0(id_dataset, "_", core_label),
          is_new_block = 
            c(
              TRUE,
              purrr::map_lgl(seq_along(id_dataset)[-1], function(i) {
                if(.core_label[[i]] == .core_label[[i-1]] & (sample_depth_upper[[i]] - sample_depth_lower[[i-1]]) <= irp_settings_cv_depth_acf_analysis$layer_gap_max) {
                  FALSE
                } else {
                  TRUE
                }
              })
            ),
          id_block = cumsum(is_new_block)
        ) |>
        dplyr::group_split(id_block) |>
        purrr::map(function(.x) {
          if(nrow(.x) < irp_settings_cv_depth_acf_analysis$core_n_samples_min) {
            tibble::tibble()
          } else {
            .x |>
              dplyr::mutate(
                sample_depth_difference = c(sample_depth_upper[-1] - sample_depth_lower[-length(sample_depth_lower)], NA_real_),
                sample_thickness = sample_depth_lower - sample_depth_upper,
              ) |>
              dplyr::group_by(id_dataset, core_label, id_block) |>
              dplyr::summarise(
                n_layers = length(sample_depth_lower),
                lag_distance_median = median(sample_depth_difference + sample_thickness, na.rm = TRUE),
                block_core_depth_min = min(sample_depth_upper),
                block_core_depth_max = min(sample_depth_lower),
                acf = 
                  stats::acf(!!cv_fold_temporal_target_variables[[i]], plot = FALSE) |>
                  list(),
                # compute confidence interval for a white noise time series
                acf_limit_upper = 
                  qnorm((1 + irp_settings_cv_depth_acf_analysis$acf_confidence_level) / 2) / sqrt(acf[[1]]$n.used),
                acf_n_relevant_lags = 
                  which(acf[[1]]$acf <= acf_limit_upper)[[1]] - 1L,
                .groups = "drop"
              )
          }
        }) |>
        dplyr::bind_rows() |>
        dplyr::mutate(
          variable =
            cv_fold_temporal_target_variables[[i]] |>
            rlang::as_name(),
          depth_range = lag_distance_median * acf_n_relevant_lags
        ) |>
        dplyr::relocate(variable, .before = dplyr::everything())
      
    }) |>
    dplyr::bind_rows()
  
  list(
    res_settings = irp_settings_cv_depth_acf_analysis,
    res = res
  )
  
}

#' Defines settings for model validation
#' 
#' @export
irp_get_cv_settings <- function() {
  
  list(
    cv_group_depth_range = 50, #---note: defined based on irp_acf_depth_range
    cv_group_spatial_buffer_radius = 500,
    cv_separate_test_fraction = 0.2,
    cv_K = 10L
  )
  
}


#' Defines temporal groups to be used to define validation folds
#' 
#' @param x One element of `irp_data_model_preprocessed`.
#' 
#' @export
irp_add_cv_group_depth <- function(x, irp_cv_settings) {
  
  res_data <- x
  
  # define the depth groups
  res_cv_group_depth <- 
    tibble::tibble(
      group_depth_upper = seq(0, max(res_data$sample_depth_lower, na.rm = TRUE), by = irp_cv_settings$cv_group_depth_range),
      group_depth_lower = group_depth_upper + irp_cv_settings$cv_group_depth_range,
      id_group = seq_along(group_depth_upper) + 2L
    )
  
  res_data |>
    dplyr::mutate(
      id_cv_group_depth =
        purrr::map_int(seq_along(id_dataset), function(i) {
          if(sample_type[[i]] != "peat") {
            1L
          } else if (is.na(sample_depth_upper[[i]]) || is.na(sample_depth_lower[[i]])) {
            2L
          } else {
            index <- 
              res_cv_group_depth$group_depth_upper <= sample_depth_upper[[i]] & res_cv_group_depth$group_depth_lower >= sample_depth_lower[[i]]
            if(sum(index) == 0) {
              index <- 
                res_cv_group_depth$group_depth_upper <= sample_depth_upper[[i]]
            }
            res_cv_group_depth$id_group[rev(which(index))[[1]]]
          }
        })
    )
  
}


#' Defines spatial groups as all objects within a certain buffer to be used to define validation folds
#' 
#' @export
irp_add_cv_group_spatial <- function(x, irp_cv_settings) {
  
  res_sf <- 
    x |>
    dplyr::filter(! is.na(sampling_longitude)) |>
    dplyr::mutate(
      index_spatial = paste0(sampling_longitude, "_", sampling_latitude)
    ) |>
    dplyr::select(index_spatial, sampling_longitude, sampling_latitude) |>
    sf::st_as_sf(coords = c("sampling_longitude", "sampling_latitude"), crs = 4326)
  
  # define buffers
  res_buffer <- 
    res_sf |>
    sf::st_buffer(dist = irp_cv_settings$cv_group_spatial_buffer_radius) |> 
    sf::st_union() |> 
    sf::st_cast('POLYGON') |>
    sf::st_as_sf() |>
    sf::st_make_valid() |>
    dplyr::mutate(
      id_cv_group_spatial = seq_along(x)
    )  

  res_sf <- 
    res_sf |>
    sf::st_join(
      res_buffer, 
      join = sf::st_within
    )
  
  res <- x
  res$id_cv_group_spatial <- NA_integer_
  res$id_cv_group_spatial[! is.na(res$sampling_longitude)] <- res_sf$id_cv_group_spatial
  res$id_cv_group_spatial[is.na(res$sampling_longitude)] <- 
    max(res_sf$id_cv_group_spatial) + 
    res |> 
    dplyr::filter(is.na(sampling_longitude)) |>
    dplyr::mutate(
      id_dataset = 
        as.factor(id_dataset) |>
        as.integer()
    ) |>
    dplyr::pull(id_dataset)
  
  res
  
}


#' Defines groups to be used for validation
#' 
#' Combines `irp_add_cv_group_depth()` and `irp_add_cv_group_spatial()`.
#' 
#' @export
irp_make_cv_groups <- function(x, irp_cv_settings) {
  
  x |>
    irp_add_cv_group_depth(irp_cv_settings) |>
    irp_add_cv_group_spatial(irp_cv_settings = irp_cv_settings) |>
    dplyr::mutate(
      id_cv_group_all = 
        paste0(id_cv_group_depth, "_", id_cv_group_spatial) |>
        as.factor() |>
        as.integer()
    ) |>
    dplyr::select(dplyr::all_of(c("id_dataset", "id_sample", "id_measurement")) | dplyr::starts_with("id_cv_group_"))
    
}


#' Defines validation folds for entries in `irp_d_model_info`
#' 
#' @export
irp_make_validation_folds <- function(irp_d_model_info, irp_cv_groups, irp_cv_settings) {
  
  # define groups for cv folds
  irp_d_model_info <- 
    irp_d_model_info |>
    irp_add_id_validation_fold_to_d_model_info()
  
  x <- 
    irp_d_model_info |>
    dplyr::filter(! duplicated(id_validation_fold)) |>
    dplyr::select(id_validation_fold, validation_mode, id_measurement_all)
  
  res <- 
    x |>
    dplyr::mutate(
      info_validation_fold =
        purrr::map(seq_along(id_validation_fold), function(i) {
          
          res_cv_groups <- 
            irp_cv_groups |>
            dplyr::filter(id_measurement %in% x$id_measurement_all[[i]])
          
          switch(
            x$validation_mode[[i]],
            "separate_test_dataset" = {
              
              # define folds (by convention, fold 2 is the test fold)
              res <- 
                tibble::tibble(
                  id_cv_group_all = 
                    sample(
                      unique(res_cv_groups$id_cv_group_all), 
                      size = length(unique(res_cv_groups$id_cv_group_all)), 
                      replace = FALSE
                    )
                ) |>
                dplyr::left_join(
                  res_cv_groups |>
                    dplyr::group_by(id_cv_group_all) |>
                    dplyr::summarise(
                      n_sample = nrow(dplyr::pick(dplyr::everything())),
                      .groups = "drop"
                    ),
                  by = "id_cv_group_all"
                ) |>
                dplyr::mutate(
                  is_training_fold =
                    cumsum(n_sample) <= ((1.0 - irp_cv_settings$cv_separate_test_fraction) * sum(n_sample)),
                  id_validation_fold =
                    dplyr::case_when(
                      is_training_fold ~ 1L,
                      TRUE ~ 2L
                    )
                ) |>
                dplyr::select(id_cv_group_all, id_validation_fold)
              
              dplyr::left_join(
                res_cv_groups,
                res,
                by = "id_cv_group_all"
              ) |>
                dplyr::select(id_measurement, id_cv_group_all, id_validation_fold)
              
            },
            "cv" = {
              res_cv_groups |>
                dplyr::mutate(
                  id_validation_fold = loo::kfold_split_grouped(K = irp_cv_settings$cv_K, x = id_cv_group_all)
                ) |>
                dplyr::select(id_measurement, id_cv_group_all, id_validation_fold)
            },
            stop("Unknown `x$validation_mode`.")
          )
          
        })
    ) |>
    dplyr::select(id_validation_fold, info_validation_fold)
  
  dplyr::left_join(
    irp_d_model_info,
    res,
    by = c("id_validation_fold")
  ) |>
    dplyr::pull(info_validation_fold)
  
}


#' Creates a data frame akin to `irp_d_model_info` for the model validation
#' 
#' Creates a row for each validation step. Sets `id_measurement_all` to 
#' the training subset of `id_measurement_all` for the validation step. Adds 
#' `id_measurement_test` with the testing subset of `id_measurement_all`.
#' 
#' @param irp_d_model_info One row of `irp_d_model_info`.
#' 
#' @param irp_validation_folds One element of `irp_validation_folds`.
#' 
#' @export
irp_make_d_model_info_validation <- function(irp_d_model_info, irp_validation_folds) {
  
  stopifnot(nrow(irp_d_model_info) == 1L)
  
  purrr::map(sort(unique(irp_validation_folds$id_validation_fold)), function(i) {
    if(irp_d_model_info$validation_mode == "separate_test_dataset" && i == 2) {
      tibble::tibble()
    } else if (irp_d_model_info$validation_mode == "separate_test_dataset" && i == 1) {
      irp_d_model_info |>
        dplyr::mutate(
          id_validation_step = i,
          id_measurement_all = 
            irp_validation_folds |>
            dplyr::filter(id_validation_fold == i) |>
            dplyr::pull(id_measurement) |>
            list(),
          id_measurement_test =
            irp_validation_folds |>
            dplyr::filter(id_validation_fold != i) |>
            dplyr::pull(id_measurement) |>
            list(),
          id_model = paste0(id_model, "_", id_validation_step),
          y =
            purrr::map(y, function(.y) {
              .y |>
                dplyr::slice(
                  irp_validation_folds |>
                    dplyr::mutate(
                      index = seq_along(id_validation_fold)
                    ) |>
                    dplyr::filter(id_validation_fold == i) |>
                    dplyr::pull(index)
                )
            })
        )
    } else {
      irp_d_model_info |>
        dplyr::mutate(
          id_validation_step = i,
          id_measurement_all = 
            irp_validation_folds |>
            dplyr::filter(id_validation_fold != i) |>
            dplyr::pull(id_measurement) |>
            list(),
          id_measurement_test =
            irp_validation_folds |>
            dplyr::filter(id_validation_fold == i) |>
            dplyr::pull(id_measurement) |>
            list(),
          id_model = paste0(id_model, "_", id_validation_step),
          y =
            purrr::map(y, function(.y) {
              .y |>
                dplyr::slice(
                  irp_validation_folds |>
                    dplyr::mutate(
                      index = seq_along(id_validation_fold)
                    ) |>
                    dplyr::filter(id_validation_fold != i) |>
                    dplyr::pull(index)
                )
            })
        )
    }
  }) |>
    dplyr::bind_rows()
  
}


#' Computes the ELPD for all validation folds for a model and combines them
#' 
#' @export
irp_make_validation_elpd <- function(irp_d_model_info_enriched_2, irp_fit_1) {
  
  id_measurement_test <- 
    irp_d_model_info_enriched_2$data_partition[[1]] |>
    dplyr::filter(for_prospectr_test) |>
    dplyr::pull(id_measurement)
  
  newdata <- 
    tibble::tibble(
      y = 
        irp_d_model_info_enriched_2$y[[1]]$y[irp_d_model_info_enriched_2$id_measurement_all[[1]] %in% id_measurement_test] |>
        scale(
          center = irp_d_model_info_enriched_2$y_center[[1]], 
          scale = irp_d_model_info_enriched_2$y_scale[[1]]
        ) |>
        as.numeric(),
      x = 
        irp_d_model_info_enriched_2$x_train[[1]] |>
        tibble::as_tibble() |>
        dplyr::filter(irp_d_model_info_enriched_2$id_measurement_all[[1]] %in% id_measurement_test) |>
        as.matrix(),
      y_err = 
        if(irp_d_model_info_enriched_2$y_has_error) {
          irp_d_model_info_enriched_2$y[[1]]$y_err[irp_d_model_info_enriched_2$id_measurement_all[[1]] %in% id_measurement_test] |> 
            scale(
              center = 0, 
              scale = irp_d_model_info_enriched_2$y_scale[[1]]
              ) |>
            as.numeric()
        } else {
          NULL
        }
    )
  
  brms::log_lik(
    object = irp_fit_1,
    newdata = newdata,
    re_formula = NULL
  ) |>
    loo::elpd()
  
}



#' Computes the RMSE for all validation folds for a model and combines them
#' 
#' @param d_model_info One row of `irp_d_model_info`.
#' 
#' @param d_model_info_validation The rows in `irp_d_model_info_validation` 
#' corresponding to `d_model_info`.
#' 
#' @param irp_fit_1_validation_config The elements of `irp_fit_1_validation_config` 
#' corresponding to `d_model_info`.
#' 
#' @param irp_fit_1_validation The elements of `irp_fit_1_validation` 
#' corresponding to `d_model_info`. 
#' 
#' @export
irp_make_validation_summary_1 <- function(d_model_info_validation, d_model_info, irp_fit_1_validation_config, irp_fit_1_validation, irp_data_model_preprocessed, irp_fit_1, irp_isotope_standards_isotope_fraction, file = "targets_rvars/irp_validation_rmse.rds") {
  
  res_test <- 
    purrr::map(seq_len(nrow(d_model_info_validation)), function(i) {
      
      x <- 
        irp_data_model_preprocessed[[d_model_info_validation$id_preprocessing[[i]]]] |>
        dplyr::filter(id_measurement %in% d_model_info_validation$id_measurement_test[[i]]) |>
        ir::ir_scale(
          center = irp_fit_1_validation_config[[i]]$irp_preprocess$scale_center, 
          scale = irp_fit_1_validation_config[[i]]$irp_preprocess$scale_scale
        )
      
      newdata <- 
        irp_predict_for_eb1079_helper_1(
          x = x, 
          config = irp_fit_1_validation_config[[i]]
        ) |>
        dplyr::mutate(
          y = 
            d_model_info$y[[1]] |>
            dplyr::filter(d_model_info$id_measurement_all[[1]] %in% d_model_info_validation$id_measurement_test[[i]]) |>
            dplyr::mutate(
              y = 
                y |>
                scale(center = d_model_info$y_center[[1]], scale = d_model_info$y_scale[[1]]) |>
                as.numeric()
            ) |>
            dplyr::pull(y)
        )
      
      if(d_model_info$y_has_error[[1]]) {
        newdata <- 
          newdata |>
          dplyr::mutate(
            y_err =
              d_model_info$y[[1]] |>
              dplyr::filter(d_model_info$id_measurement_all[[1]] %in% d_model_info_validation$id_measurement_test[[i]]) |>
              dplyr::mutate(
                y_err = 
                  y_err |>
                  scale(center = FALSE, scale = d_model_info$y_scale[[1]]) |>
                  as.numeric()
              ) |>
              dplyr::pull(y_err)
          ) |>
          dplyr::relocate(
            y_err, .after = "y"
          )
      }
      
      res <- 
        tibble::tibble(
          id_model = d_model_info$id_model,
          y = newdata$y,
          yhat = 
            brms::posterior_predict(
              irp_fit_1_validation[[i]],
              re_formula = NULL,
              newdata = newdata
            ) |>
            posterior::rvar()
        ) |>
        dplyr::mutate(
          dplyr::across(
            dplyr::all_of(c("y", "yhat")),
            function(.x) {
              .x |>
                magrittr::multiply_by(d_model_info$y_scale[[1]]) |>
                magrittr::add(d_model_info$y_center[[1]])
            }
          )
        )
      
      # convert atom fraction to delta values
      if(stringr::str_detect(d_model_info$target_variable, pattern = "^(d13C|d15N)")) {
        isotope <- stringr::str_extract(d_model_info$target_variable, pattern = "^(d13C|d15N)")
        r <- 
          switch(
            isotope,
            "d13C" = irp_isotope_standards_isotope_fraction$r_13C,
            "d15N" = irp_isotope_standards_isotope_fraction$r_15N
          )
        res <- 
          res |>
          dplyr::mutate(
            dplyr::across(
              dplyr::all_of(c("y", "yhat")),
              function(.x) {
                if(is.numeric(.x)) {
                  irp_atom_fraction_to_delta(.x, r = r)  
                } else {
                  posterior::rfun(irp_atom_fraction_to_delta, rvar_args = "x")(.x, r = r)
                }
              }
            )
          )
      }
      
      res
      
    }) |>
    irp_do_call("rbind") |>
    dplyr::mutate(
      rmse_test = irp_rmse_rvar(yhat = yhat, y = y),
      y_minus_yhat_test = posterior::rvar_mean(y - yhat)
    ) |>
    dplyr::select(id_model, rmse_test, y_minus_yhat_test) |>
    dplyr::slice(1)
  
  res_train <- 
    tibble::tibble(
      y = d_model_info$y[[1]]$y,
      yhat = 
        brms::posterior_predict(
          irp_fit_1,
          re_formula = NULL,
          newdata = NULL
        ) |>
        posterior::rvar()
    ) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(c("yhat")),
        function(.x) {
          .x |>
            magrittr::multiply_by(d_model_info$y_scale[[1]]) |>
            magrittr::add(d_model_info$y_center[[1]])
        }
      )
    )
  
  # convert atom fraction to delta values
  if(stringr::str_detect(d_model_info$target_variable, pattern = "^(d13C|d15N)")) {
    isotope <- stringr::str_extract(d_model_info$target_variable, pattern = "^(d13C|d15N)")
    r <- 
      switch(
        isotope,
        "d13C" = irp_isotope_standards_isotope_fraction$r_13C,
        "d15N" = irp_isotope_standards_isotope_fraction$r_15N
      )
    res_train <- 
      res_train |>
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(c("y", "yhat")),
          function(.x) {
            if(is.numeric(.x)) {
              irp_atom_fraction_to_delta(.x, r = r)  
            } else {
              posterior::rfun(irp_atom_fraction_to_delta, rvar_args = "x")(.x, r = r)
            }
          }
        )
      )
  }
  
  # MCSE
  res_train <- 
    res_train |>
    dplyr::mutate(
      y_mcse_mean = posterior::mcse_mean(yhat),
      y_mcse_sd = posterior::mcse_sd(yhat),
      y_mcse_lower = posterior::mcse_quantile(yhat, probs = 0.025),
      y_mcse_upper = posterior::mcse_quantile(yhat, probs = 0.975)
    )
  
  res <- 
    res_test |>
    dplyr::mutate(
      rmse_train_all = irp_rmse_rvar(yhat = res_train$yhat, y = res_train$y),
      y_minus_yhat_train_all = posterior::rvar_mean(res_train$y - res_train$yhat),
      rmse_train_id_test =
        if(d_model_info$validation_mode == "cv") {
          rmse_train_all
        } else {
          res_train |>
            dplyr::filter(d_model_info$id_measurement_all[[1]] %in% d_model_info_validation$id_measurement_test[[1]]) |>
            dplyr::mutate(
              rmse = irp_rmse_rvar(yhat = res_train$yhat, y = res_train$y)
            ) |>
            dplyr::slice(1) |>
            dplyr::pull(rmse)
        },
      y_minus_yhat_train_id_test =
        if(d_model_info$validation_mode == "cv") {
          y_minus_yhat_train_all
        } else {
          res_train |>
            dplyr::filter(d_model_info$id_measurement_all[[1]] %in% d_model_info_validation$id_measurement_test[[1]]) |>
            dplyr::mutate(
              y_minus_yhat = posterior::rvar_mean(y - yhat)
            ) |>
            dplyr::slice(1) |>
            dplyr::pull(y_minus_yhat)
        },
      y_mcse_mean = list(res_train$y_mcse_mean),
      y_mcse_sd = list(res_train$y_mcse_sd),
      y_mcse_lower = list(res_train$y_mcse_lower),
      y_mcse_upper = list(res_train$y_mcse_upper),
      num_divergent = 
        irp_fit_1 |>
        brms::nuts_params(pars = c("divergent__")) |>
        dplyr::pull(Value) |>
        sum(),
      max_rhat =
        irp_fit_1 |>
        brms::rhat() |>
        max()
    )
  
  saveRDS_rvars(res, file)
  
  file
  
}


#' Compares models by their ELPD
#' 
#' @export
irp_make_validation_elpd_compare <- function(irp_fit_1_validation_elpd, irp_d_model_info) {
  
  res <- 
    irp_d_model_info |>
    dplyr::mutate(
      elpd = irp_fit_1_validation_elpd
    ) |>
    dplyr::group_by(target_variable, do_dimension_reduction, uses_only_high_quality_mirs) |>
    dplyr::summarise(
      elpd_compare = {
        
        # compute Delta elpd
        res <- 
          loo::loo_compare(elpd) |>
          as.data.frame()
        res$model_label <- rownames(res)
        
        # find best model
        res$is_best_model <- c(TRUE, rep(FALSE, nrow(res) - 1))
        for(i in seq_len(nrow(res))[-2]) {
          if(abs(res$elpd_diff[[i]]) <= 4 || (abs(res$elpd_diff[[i]]) - 2 * res$se_diff[[i]]) <= 0) {
            res$is_best_model[[i]] <- TRUE
          }
        }
        
        list(res)
      },
      .groups = "drop"
    ) |>
    tidyr::unnest("elpd_compare")
  
  res <-
    res |>
    dplyr::left_join(
      irp_d_model_info |>
        dplyr::mutate(
          model_label = names(irp_fit_1_validation_elpd)
        ),
      by = c("model_label", "target_variable", "do_dimension_reduction", "uses_only_high_quality_mirs")
    )
  
  res
  
}





