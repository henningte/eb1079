#' Fills gaps in the pmird database
#' 
#' @export
irp_make_pmird_gap_filled <- function(irp_pmird_mirs, irp_pmird_gap_filled_helper_1_combined, irp_m1, irp_m2, irp_d_compounds_standard_enthalpies_of_formation, irp_pmird_units) {
  
  index <- ! purrr::map_lgl(irp_pmird_gap_filled_helper_1_combined, is.null)
  target_variables <- stringr::str_extract(names(index), pattern = paste0("(", paste(irp_get_target_variables(), collapse = "|"), ")_\\d{1}$"))[index]
  target_variables_transfer_function <- c("volume_fraction_solids_1", "non_macroporosity_1", "macroporosity_1", "saturated_hydraulic_conductivity_1", "specific_heat_capacity_1", "dry_thermal_conductivity_1")
  target_variables <- c(target_variables_transfer_function, target_variables)
  
  target_variables_ordered <- 
    c(target_variables, target_variables_transfer_function)
  
  
  irp_pmird_gap_filled_helper_1_combined <- irp_pmird_gap_filled_helper_1_combined[index]
  irp_pmird_gap_filled_helper_1_combined <- purrr::map(irp_pmird_gap_filled_helper_1_combined, readRDS_rvars)
  irp_pmird_gap_filled_helper_1_combined <- 
    purrr::map(irp_pmird_gap_filled_helper_1_combined, function(.x) {
      .x$meta <- NULL
      .x
    }) |>
    purrr::transpose() |>
    purrr::map(function(.x) {
      .x |>
        unname() |>
        irp_do_call("cbind")
    })
  
  
  spectral_preprocessing_clip_range <- readRDS(system.file("extdata", "model_carbon_content_1_config.rds", package = "irpeatmodels"))$irp_preprocess$clip_range[[1]]
  res_spectra <- 
    irp_pmird_mirs |>
    dplyr::filter(x_variable_min <= spectral_preprocessing_clip_range$start & x_variable_max >= spectral_preprocessing_clip_range$end) |>
    dplyr::filter((mir_mode == "absorbance_ftir" | is.na(mir_mode)) & sample_type %in% c("peat", "vegetation", "litter") & id_dataset != 16)
  
  res_yhat_auxilliarly <- 
    res_spectra |>
    irpeat::irp_predict(
      variable = c("volume_fraction_solids_1", "non_macroporosity_1", "macroporosity_1", "saturated_hydraulic_conductivity_1", "dry_thermal_conductivity_1"), 
      bulk_density = res_spectra$bulk_density,
      do_summary = FALSE,
      return_as_list = FALSE
    ) |>
    irpeat::irp_predict(
      variable = c("specific_heat_capacity_1"), 
      nitrogen_content = res_spectra$N,
      do_summary = FALSE,
      return_as_list = FALSE
    ) |>
    # NOSC
    dplyr::bind_cols(
      res_spectra |>
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
          nosc_2 = 
            elco::elco_nosc(C = C, H = H, N = N, O = O) |>
            as.numeric()
        ) |>
        dplyr::select(nosc_2)
    ) |>
    # dgf0
    dplyr::bind_cols(
      {
        irp_standard_entropies_elements <- 
          irp_make_standard_entropies_elements()
        
        res <- 
          res_spectra |>
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
            dgf0_3 = 
              irp_estimate_dgf0_per_C(
                x = as.data.frame(res),
                irp_m1 = irp_m1, 
                irp_m2 = irp_m2, 
                irp_d_compounds_standard_enthalpies_of_formation = irp_d_compounds_standard_enthalpies_of_formation
              ) |>
              dplyr::pull(dgf0_per_C),
          ) |>
          dplyr::select(dgf0_3) |>
          dplyr::mutate(
            dgf0_3 =
              dplyr::case_when(
                ! is.na(res_spectra$C) & ! is.na(res_spectra$H) & ! is.na(res_spectra$O) & ! is.na(res_spectra$N) ~ dgf0_3,
                TRUE ~ NA
              )
          )
      }
    ) |>
    # element ratios
    dplyr::bind_cols(
      res_spectra |>
        dplyr::mutate(
          C_to_N_3 = C/N,
          O_to_C_3 = O/C,
          H_to_C_3 = H/C
        ) |>
        dplyr::select(dplyr::all_of(c("C_to_N_3", "O_to_C_3", "H_to_C_3")))
    ) |>
    dplyr::select(dplyr::any_of(target_variables_ordered))
  
  
  res <- 
    list(
      meta = 
        res_spectra |>
        dplyr::select(id_dataset, id_sample, id_measurement),
      yhat = 
        irp_pmird_gap_filled_helper_1_combined$yhat |>
        dplyr::select(dplyr::all_of(target_variables_ordered)),
      yhat_auxilliary =
        res_yhat_auxilliarly |>
        dplyr::bind_rows(
          irp_pmird_gap_filled_helper_1_combined$yhat |>
          dplyr::select(dplyr::starts_with(target_variables_ordered) & ! dplyr::ends_with("_in_pd")) |>
            dplyr::slice(0)
        ) |>
        dplyr::select(dplyr::all_of(target_variables_ordered)),
      is_in_training_pd =
        irp_pmird_gap_filled_helper_1_combined$is_in_training_pd |>
        dplyr::select(dplyr::all_of(target_variables_ordered)),
      is_in_testing_pd =
        irp_pmird_gap_filled_helper_1_combined$is_in_testing_pd |>
        dplyr::select(dplyr::all_of(target_variables_ordered))
    )
  
  res_file <- "targets_rvars/irp_pmird_gap_filled.rds"
  saveRDS_rvars(res, res_file)
  res_file
  
}



#' Helper function that makes the predictions for gap filling and checks prediction domains
#' 
#' @export
irp_make_pmird_gap_filled_helper_1 <- function(irp_d_model_info_enriched_2, irp_fit_1_irpeat_model_config, irp_id_model_best, irp_data_model_preprocessed) {
  
  if(! irp_d_model_info_enriched_2$id_model %in% irp_id_model_best) {
    return(NULL)
  }
  
  config <- irp_fit_1_irpeat_model_config
  model <- readRDS(system.file("extdata", paste0("model_", stringr::str_replace(config$target_variable, pattern = "_\\d{1}$", replacement = "_1"), ".rds"), package = "irpeatmodels"))
  
  # preprocessing
  res_spectra <- 
    irp_data_model_preprocessed[[irp_d_model_info_enriched_2$id_preprocessing]] |>
    dplyr::filter(x_variable_min <= config$irp_preprocess$clip_range[[1]]$start & x_variable_max >= config$irp_preprocess$clip_range[[1]]$end) |>
    dplyr::filter((mir_mode == "absorbance_ftir" | is.na(mir_mode)) & sample_type %in% c("peat", "vegetation", "litter") & id_dataset != 16) |>
    ir::ir_scale(center = config$irp_preprocess$scale_center, scale = config$irp_preprocess$scale_scale)
  
  # check prediction domain
  res_pd_train <- 
    irpeat::irp_is_in_prediction_domain(
      x = res_spectra, 
      prediction_domain = config$prediction_domain$train
    ) |>
    dplyr::select(is_in_prediction_domain) |>
    setNames(nm = paste0(config$target_variable))
  
  res_pd_test <- 
    irpeat::irp_is_in_prediction_domain(
      x = res_spectra, 
      prediction_domain = config$prediction_domain$test
    ) |>
    dplyr::select(is_in_prediction_domain) |>
    setNames(nm = paste0(config$target_variable))
    
  # predictions
  newdata <-
    tibble::tibble(
      x =
        res_spectra |>
        ir::ir_flatten() |>
        dplyr::select(-1) |>
        t() |>
        tibble::as_tibble(.name_repair = "minimal") |>
        setNames(nm = paste0("V", seq_along(res_spectra$spectra[[1]]$x))) |>
        as.matrix(),
      y_err =
        if(stringr::str_detect(config$target_variable, pattern = "(dgf0_)")) {
          0.00001
        } else {
          NULL
        }
    )
  
  yhat <-
    brms::posterior_predict(object = model, newdata = newdata) |>
    as.data.frame()
  
  yhat <-
    irpeat:::irp_summarize_predictions(
      x = yhat * config$model_scale$y_scale + config$model_scale$y_center,
      x_unit = config$unit,
      do_summary = FALSE,
      return_as_list = FALSE,
      summary_function_mean = mean,
      summary_function_sd = stats::sd
    )
  
  res <- 
    list(
      meta = 
        res_spectra |>
        dplyr::select(id_dataset, id_sample, id_measurement),
      yhat = 
        tibble::tibble(
          yhat = yhat
        ) |>
        setNames(nm = config$target_variable),
      is_in_training_pd = res_pd_train,
      is_in_testing_pd = res_pd_test
    )
  
  switch(
    stringr::str_remove(config$target_variable, pattern = "_\\d{1}$"),
    "bulk_density" = {
      res_porosity <- 
        irpeat:::irp_porosity_1(
          x = res_spectra, 
          do_summary = FALSE, 
          return_as_list = FALSE, 
          bulk_density = yhat
        )
      res_hydraulic_conductivity <- 
        irpeat::irp_saturated_hydraulic_conductivity_1(
          x = res_spectra, 
          do_summary = FALSE, 
          return_as_list = FALSE, 
          bulk_density = yhat
        )
      res_dry_thermal_conductivity_1 <- 
        irpeat::irp_dry_thermal_conductivity_1(
          x = res_spectra, 
          do_summary = FALSE, 
          return_as_list = FALSE, 
          bulk_density = yhat
        )
      
      res$yhat <-
        res$yhat |>
        dplyr::bind_cols(
          res_porosity |>
            dplyr::select("macroporosity_1", "non_macroporosity_1", "volume_fraction_solids_1"),
          res_hydraulic_conductivity |>
            dplyr::select("saturated_hydraulic_conductivity_1"),
          res_dry_thermal_conductivity_1 |>
            dplyr::select("dry_thermal_conductivity_1")
        )
      
      res$is_in_training_pd <- 
        res$is_in_training_pd |>
        dplyr::bind_cols(
          purrr::map(c("macroporosity_1", "non_macroporosity_1", "volume_fraction_solids_1", "saturated_hydraulic_conductivity_1", "dry_thermal_conductivity_1"), function(.x) {
            res$is_in_training_pd |>
              setNames(nm = .x)
          }) |>
            dplyr::bind_cols()
        )
      
      res$is_in_testing_pd <- 
        res$is_in_testing_pd |>
        dplyr::bind_cols(
          purrr::map(c("macroporosity_1", "non_macroporosity_1", "volume_fraction_solids_1", "saturated_hydraulic_conductivity_1", "dry_thermal_conductivity_1"), function(.x) {
            res$is_in_testing_pd |>
              setNames(nm = .x)
          }) |>
            dplyr::bind_cols()
        )
      
    },
    "nitrogen_content" = {
      res_specific_heat_capacity <- 
        irpeat::irp_specific_heat_capacity_1(
          x = res_spectra, 
          do_summary = FALSE, 
          return_as_list = FALSE, 
          nitrogen_content = yhat
        )
      
      res$yhat <-
        res$yhat |>
        dplyr::bind_cols(
          res_specific_heat_capacity |>
            dplyr::select("specific_heat_capacity_1")
        )
      
      res$is_in_training_pd <- 
        res$is_in_training_pd |>
        dplyr::bind_cols(
          purrr::map(c("specific_heat_capacity_1"), function(.x) {
            res$is_in_training_pd |>
              setNames(nm = .x)
          }) |>
            dplyr::bind_cols()
        )
      
      res$is_in_testing_pd <- 
        res$is_in_testing_pd |>
        dplyr::bind_cols(
          purrr::map(c("specific_heat_capacity_1"), function(.x) {
            res$is_in_testing_pd |>
              setNames(nm = .x)
          }) |>
            dplyr::bind_cols()
        )
    }
  )
  
  res_file <- paste0("targets_rvars/irp_pmird_gap_filled_", config$target_variable, ".rds")
  saveRDS_rvars(res, res_file)
  res_file
  
}


