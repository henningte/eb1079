# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(magrittr)
library(quantities)
library(ir)
library(crew)
library(ggplot2)

crew_sequential <-
  crew::crew_controller_local(
    name = "crew_sequential",
    workers = 1L,
    tasks_max = 1L
  )

crew_parallel_1 <-
  crew::crew_controller_local(
    name = "crew_parallel_1",
    workers = 4L
  )

crew_parallel_2 <-
  crew::crew_controller_local(
    name = "crew_parallel_2",
    workers = 8L
  )

# Set target options:
tar_option_set(
  packages = c("tibble", "CHNOSZ", "elco", "ggplot2", "dm", "ir", "cmdstanr", "brms", "posterior", "pls", "dimreduce", "sf"),
  format = "rds",
  controller = crew_controller_group(crew_sequential, crew_parallel_1, crew_parallel_2)
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()

# to avoid deletion of aux files
options(tinytex.clean = FALSE)

# Replace the target list below with your own:
list(
  #### data ####
  tar_target(
    irp_db_obigt_organic,
    command = irp_make_db_obigt_organic()
  ),
  tar_target(
    irp_d_battley_1999_table_1,
    command = irp_make_d_battley_1999_table_1()
  ),
  tar_target(
    irp_d_compounds_standard_enthalpies_of_formation,
    command = irp_make_d_compounds_standard_enthalpies_of_formation()
  ),
  tar_target(
    irp_pmird_mirs,
    command = irp_get_pmird_mirs()
  ),
  tar_target(
    irp_liu2019_pmird,
    command = irp_get_liu2019_pmird()
  ),
  tar_target(
    irp_pmird_units,
    command = irp_get_pmird_units()
  ),
  tar_target(
    irp_wang2015b,
    command = irp_get_wang2015b()
  ),
  tar_target(
    irp_d_gnatowski2022,
    command = irp_make_d_gnatowski2022()
  ),
  tar_target(
    irp_d_cp_air,
    command = irp_make_d_cp_air()
  ),
  tar_target(
    irp_d_oconnor2022,
    command = irp_make_d_oconnor2022()
  ),
  tar_target(
    irp_isotope_standards_isotope_fraction,
    command =
      list(
        r_13C = 0.0112372,
        r_15N = 0.0036765
      )
  ),
  tar_target(
    irp_mass_density_fe,
    command = 7.874 # from Rumble.2016
  ),
  #### settings ####
  tar_target(
    irp_mcmc_settings,
    command = irp_get_mcmc_settings()
  ),
  #### auxiliary models ####
  tar_target(
    irp_m1,
    command = 
      irp_make_m1(
        irp_db_obigt_organic = irp_db_obigt_organic, 
        irp_mcmc_settings = irp_mcmc_settings
      )
  ),
  tar_target(
    irp_m2,
    command = 
      irp_make_m2(
        irp_db_obigt_organic = irp_db_obigt_organic,
        irp_d_battley_1999_table_1 = irp_d_battley_1999_table_1,
        irp_mcmc_settings = irp_mcmc_settings
      )
  ),
  tar_target(
    irp_m3,
    command =
      irp_make_m3(
        irp_liu2019_pmird = irp_liu2019_pmird, 
        irp_wang2015b = irp_wang2015b, 
        irp_mcmc_settings = irp_mcmc_settings
      )
  ),
  tar_target(
    irp_m4,
    command =
      irp_make_m4(
        irp_liu2019_pmird = irp_liu2019_pmird, 
        irp_mcmc_settings = irp_mcmc_settings
      )
  ),
  tar_target(
    irp_m5,
    command =
      irp_make_m5(
        irp_d_gnatowski2022 = irp_d_gnatowski2022, 
        irp_d_cp_air = irp_d_cp_air,
        irp_mcmc_settings = irp_mcmc_settings
      )
  ),
  tar_target(
    irp_m6,
    command =
      irp_make_m6(
        irp_d_oconnor2022 = irp_d_oconnor2022, 
        irp_mcmc_settings = irp_mcmc_settings
      )
  ),
  #### spectral preprocessing ####
  tar_target(
    irp_mir_preprocessing_settings,
    command = irp_get_mir_preprocessing_settings()
  ),
  tar_target(
    irp_data_model_preprocessed,
    command = 
      irp_make_data_model_preprocessed(
        irp_pmird_mirs = irp_pmird_mirs, 
        irp_mir_preprocessing_settings = irp_mir_preprocessing_settings
      )
  ),
  #### spectral prediction models ####
  tar_target(
    irp_threshold_model_validation_sample_number,
    command = 500
  ),
  tar_target(
    irp_target_variables,
    command = irp_get_target_variables()
  ),
  tar_target(
    irp_d_model_info,
    pattern = map(irp_target_variables),
    command = 
      irp_make_d_model_info(
        x = irp_target_variables, 
        irp_threshold_model_validation_sample_number = irp_threshold_model_validation_sample_number, 
        irp_data_model_preprocessed = irp_data_model_preprocessed, 
        irp_isotope_standards_isotope_fraction = irp_isotope_standards_isotope_fraction, 
        irp_pmird_units = irp_pmird_units, 
        irp_m1 = irp_m1, 
        irp_m2 = irp_m2, 
        irp_d_compounds_standard_enthalpies_of_formation = irp_d_compounds_standard_enthalpies_of_formation,
        irp_mass_density_fe
      ),
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_1")
      )
  ),
  tar_target(
    irp_d_model_info_id_model,
    command = irp_d_model_info$id_model
  ),
  tar_target(
    irp_dimension_reduction_model_settings_1,
    command = irp_get_dimension_reduction_model_settings_1()
  ),
  tar_target(
    irp_dimension_reduction_model_1,
    pattern = map(irp_d_model_info_id_model),
    iteration = "list",
    command = 
      irp_compute_dimension_reduction_model_1(
        x = 
          irp_d_model_info |>
          dplyr::filter(id_model == irp_d_model_info_id_model), 
        irp_data_model_preprocessed = irp_data_model_preprocessed, 
        irp_dimension_reduction_model_settings_1 = irp_dimension_reduction_model_settings_1
      ),
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_1")
      )
  ),
  tar_target(
    irp_data_for_brms_models_1,
    pattern = map(irp_d_model_info_id_model, irp_dimension_reduction_model_1),
    iteration = "list",
    command =
      irp_prepare_data_for_brms_models_1(
        x = 
          irp_d_model_info |>
          dplyr::filter(id_model == irp_d_model_info_id_model),
        newdata = NULL, 
        irp_dimension_reduction_model_1 = irp_dimension_reduction_model_1, 
        irp_data_model_preprocessed = irp_data_model_preprocessed
      ),
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_1")
      )
  ),
  tar_target(
    irp_brms_compiled_model,
    command =
      irp_make_brms_compiled_model(
        irp_d_model_info = irp_d_model_info, 
        irp_data_for_brms_models_1 = irp_data_for_brms_models_1,
        irp_isotope_standards_isotope_fraction = irp_isotope_standards_isotope_fraction
      )
  ),
  tar_target(
    irp_fit_1,
    pattern = map(irp_d_model_info_id_model, irp_data_for_brms_models_1),
    iteration = "list",
    command = 
      irp_make_fit_1(
        x = 
          irp_d_model_info |>
          dplyr::filter(id_model == irp_d_model_info_id_model),
        irp_brms_compiled_model = irp_brms_compiled_model, 
        irp_data_for_brms_models_1 = irp_data_for_brms_models_1, 
        irp_isotope_standards_isotope_fraction = irp_isotope_standards_isotope_fraction, 
        irp_mcmc_settings = irp_mcmc_settings
      ),
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_1")
      )
  ),
  tar_target(
    irp_fit_1_config,
    pattern = map(irp_d_model_info_id_model, irp_dimension_reduction_model_1),
    iteration = "list",
    command = 
      irp_make_mir_preprocessing_config(
        irp_mir_preprocessing_settings = irp_mir_preprocessing_settings, 
        irp_d_model_info = 
          irp_d_model_info |>
          dplyr::filter(id_model == irp_d_model_info_id_model), 
        irp_dimension_reduction_model_1 = irp_dimension_reduction_model_1
      ),
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_2")
      )
  ),
  #### validation ####
  tar_target(
    irp_acf_depth_range,
    command = 
      irp_estimate_acf_depth_range(
        irp_data_model_preprocessed = irp_data_model_preprocessed[[1]]
      )
  ),
  tar_target(
    irp_cv_settings,
    command =
      irp_get_cv_settings()
  ),
  tar_target(
    irp_cv_groups,
    command = 
      irp_make_cv_groups(
        x = irp_data_model_preprocessed[[1]], 
        irp_cv_settings = irp_cv_settings
      )
  ),
  tar_target(
    irp_validation_folds,
    command = 
      irp_make_validation_folds(
        irp_d_model_info = irp_d_model_info, 
        irp_cv_groups = irp_cv_groups, 
        irp_cv_settings = irp_cv_settings
      ),
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_1")
      )
  ),
  tar_target(
    irp_d_model_info_validation,
    pattern = map(irp_d_model_info_id_model, irp_validation_folds),
    command = 
      irp_make_d_model_info_validation(
        irp_d_model_info = 
          irp_d_model_info |>
          dplyr::filter(id_model == irp_d_model_info_id_model), 
        irp_validation_folds = irp_validation_folds[[1]]
      ),
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_1")
      )
  ),
  tar_target(
    irp_d_model_info_validation_id_model,
    command = irp_d_model_info_validation$id_model
  ),
  tar_target(
    irp_dimension_reduction_model_1_validation,
    pattern = map(irp_d_model_info_id_model),
    iteration = "list",
    command = {
      index_validation <- 
        which(stringr::str_detect(irp_d_model_info_validation$id_model, pattern = paste0("^", irp_d_model_info_id_model)))
      purrr::map(index_validation, function(i) {
        irp_compute_dimension_reduction_model_1(
          x = 
            irp_d_model_info_validation |>
            dplyr::slice(i), 
          irp_data_model_preprocessed = irp_data_model_preprocessed, 
          irp_dimension_reduction_model_settings_1 = irp_dimension_reduction_model_settings_1
        )
      })
    },
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_2")
      )
  ),
  tar_target(
    irp_data_for_brms_models_1_validation,
    pattern = map(irp_d_model_info_id_model, irp_dimension_reduction_model_1_validation),
    iteration = "list",
    command = {
      index_validation <- 
        which(stringr::str_detect(irp_d_model_info_validation$id_model, pattern = paste0("^", irp_d_model_info_id_model)))
      purrr::map(seq_along(index_validation), function(i) {
        irp_prepare_data_for_brms_models_1(
          x = 
            irp_d_model_info_validation |>
            dplyr::slice(index_validation[[i]]), 
          newdata = NULL, 
          irp_dimension_reduction_model_1 = irp_dimension_reduction_model_1_validation[[i]], 
          irp_data_model_preprocessed = irp_data_model_preprocessed
        ) 
      })
    },
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_2")
      )
  ),
  tar_target(
    irp_fit_1_validation,
    pattern = map(irp_d_model_info_id_model, irp_data_for_brms_models_1_validation),
    iteration = "list",
    command = {
      index_validation <- 
        which(stringr::str_detect(irp_d_model_info_validation$id_model, pattern = paste0("^", irp_d_model_info_id_model)))
      purrr::map(seq_along(index_validation), function(i) {
        irp_make_fit_1(
          x = 
            irp_d_model_info_validation |>
            dplyr::slice(index_validation[[i]]), 
          irp_brms_compiled_model = irp_brms_compiled_model, 
          irp_data_for_brms_models_1 = irp_data_for_brms_models_1_validation[[i]], 
          irp_isotope_standards_isotope_fraction = irp_isotope_standards_isotope_fraction, 
          irp_mcmc_settings = irp_mcmc_settings
        )
      })
    },
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_2")
      )
  ),
  tar_target(
    irp_fit_1_validation_config,
    pattern = map(irp_d_model_info_id_model), #, irp_dimension_reduction_model_1_validation
    iteration = "list",
    command = {
      index_validation <- 
        which(stringr::str_detect(irp_d_model_info_validation_id_model, pattern = paste0("^", irp_d_model_info_id_model)))
      purrr::map(seq_along(index_validation), function(i) {
        irp_make_mir_preprocessing_config(
          irp_mir_preprocessing_settings = irp_mir_preprocessing_settings, 
          irp_d_model_info = 
            irp_d_model_info_validation |>
            dplyr::slice(index_validation[[i]]), 
          irp_dimension_reduction_model_1 = irp_dimension_reduction_model_1_validation[[i]]
        )
      })
    },
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_2")
      )
  ),
  tar_target(
    irp_fit_1_validation_elpd,
    pattern = map(irp_d_model_info_id_model, irp_d_model_info_validation, irp_fit_1_validation_config, irp_fit_1_validation, irp_dimension_reduction_model_1_validation),
    iteration = "list",
    command = 
      irp_make_validation_elpd(
        d_model_info_validation = irp_d_model_info_validation, 
        d_model_info = 
          irp_d_model_info |>
          dplyr::filter(id_model == irp_d_model_info_id_model), 
        irp_fit_1_validation_config = 
          purrr::map(seq_along(irp_fit_1_validation_config), function(i) {
            irp_fit_1_validation_config[[i]]$dimension_reduction_model <- irp_dimension_reduction_model_1_validation
            irp_fit_1_validation_config
          }), 
        irp_fit_1_validation = irp_fit_1_validation, 
        irp_data_model_preprocessed = irp_data_model_preprocessed
      ),
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_1")
      )
  ),
  tar_target(
    irp_validation_summary_1,
    pattern = map(irp_d_model_info_id_model, irp_d_model_info_validation, irp_fit_1_validation_config, irp_fit_1_validation, irp_fit_1, irp_dimension_reduction_model_1_validation),
    command = 
      irp_make_validation_summary_1(
        d_model_info_validation = irp_d_model_info_validation, 
        d_model_info = 
          irp_d_model_info |>
          dplyr::filter(id_model == irp_d_model_info_id_model), 
        irp_fit_1_validation_config = 
          purrr::map(seq_along(irp_fit_1_validation_config), function(i) {
            irp_fit_1_validation_config[[i]]$dimension_reduction_model <- irp_dimension_reduction_model_1_validation
            irp_fit_1_validation_config
          }),
        irp_fit_1_validation = irp_fit_1_validation, 
        irp_data_model_preprocessed = irp_data_model_preprocessed,
        irp_fit_1 = irp_fit_1,
        irp_isotope_standards_isotope_fraction = irp_isotope_standards_isotope_fraction,
        file = paste0("targets_rvars/irp_validation_summary_1_", irp_d_model_info_id_model, ".rds")
      ),
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_1")
      ),
    format = "file"
  ),
  tar_target(
    irp_fit_1_validation_elpd_diff,
    command =
      irp_make_validation_elpd_compare(
        irp_fit_1_validation_elpd = irp_fit_1_validation_elpd, 
        irp_d_model_info = irp_d_model_info
      )
  ),
  tar_target(
    irp_fit_1_y_yhat,
    pattern = map(irp_d_model_info_id_model, irp_d_model_info_validation, irp_fit_1_validation_config, irp_fit_1_validation, irp_fit_1, irp_validation_folds),
    command = 
      irp_fit_1_make_y_yhat(
        d_model_info_validation = irp_d_model_info_validation, 
        d_model_info = 
          irp_d_model_info |>
          dplyr::filter(id_model == irp_d_model_info_id_model), 
        irp_fit_1_validation_config = irp_fit_1_validation_config, 
        irp_fit_1_validation = irp_fit_1_validation, 
        irp_data_model_preprocessed = irp_data_model_preprocessed,
        irp_fit_1 = irp_fit_1,
        irp_isotope_standards_isotope_fraction = irp_isotope_standards_isotope_fraction,
        irp_validation_folds = irp_validation_folds[[1]],
        file = paste0("targets_rvars/irp_fit_1_y_yhat_", irp_d_model_info_id_model, ".rds")
      ),
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_1")
      ),
    format = "file"
  ),
  #### plots ####
  tar_target(
    irp_plot_1,
    command = 
      irp_make_plot_1(
        irp_db_obigt_organic = irp_db_obigt_organic, 
        irp_m1 = irp_m1, 
        file_plot = "figures/irp_plot_1.pdf"
      ),
    format = "file"
  ),
  tar_target(
    irp_plot_2,
    command = 
      irp_make_plot_2(
        irp_m1 = irp_m1, 
        file_plot = "figures/irp_plot_2.pdf"
      ),
    format = "file"
  ),
  tar_target(
    irp_plot_3,
    command = 
      irp_make_plot_3(
        irp_m2 = irp_m2, 
        file_plot = "figures/irp_plot_3.pdf"
      ),
    format = "file"
  ),
  tar_target(
    irp_plot_4,
    command = 
      irp_make_plot_4(
        irp_m2 = irp_m2, 
        file_plot = "figures/irp_plot_4.pdf"
      ),
    format = "file"
  ),
  tar_target(
    irp_plot_5,
    command = 
      irp_make_plot_5(
        irp_db_obigt_organic = irp_db_obigt_organic,
        irp_m1 = irp_m1, 
        irp_m2 = irp_m2, 
        irp_d_compounds_standard_enthalpies_of_formation = irp_d_compounds_standard_enthalpies_of_formation, 
        file_plot = "figures/irp_plot_5.pdf"
      ),
    format = "file"
  ),
  tar_target(
    irp_plot_6,
    command = 
      irp_make_plot_6(
        irp_m3 = irp_m3, 
        file_plot = "figures/irp_plot_6.pdf"
      ),
    format = "file"
  ),
  tar_target(
    irp_plot_7,
    command = 
      irp_make_plot_7(
        irp_m3 = irp_m3, 
        file_plot = "figures/irp_plot_7.pdf"
      ),
    format = "file"
  ),
  tar_target(
    irp_plot_8,
    command = 
      irp_make_plot_8(
        irp_m4 = irp_m4, 
        file_plot = "figures/irp_plot_8.pdf"
      ),
    format = "file"
  ),
  tar_target(
    irp_plot_9,
    command = 
      irp_make_plot_9(
        irp_m4 = irp_m4, 
        file_plot = "figures/irp_plot_9.pdf"
      ),
    format = "file"
  ),
  tar_target(
    irp_plot_10,
    command = 
      irp_make_plot_10(
        irp_m5 = irp_m5, 
        file_plot = "figures/irp_plot_10.pdf"
      ),
    format = "file"
  ),
  tar_target(
    irp_plot_11,
    command = 
      irp_make_plot_11(
        irp_m5 = irp_m5, 
        file_plot = "figures/irp_plot_11.pdf"
      ),
    format = "file"
  ),
  tar_target(
    irp_plot_12,
    command = 
      irp_make_plot_12(
        irp_m6 = irp_m6, 
        file_plot = "figures/irp_plot_12.pdf"
      ),
    format = "file"
  ),
  tar_target(
    irp_plot_13,
    command = 
      irp_make_plot_13(
        irp_m6 = irp_m6, 
        file_plot = "figures/irp_plot_13.pdf"
      ),
    format = "file"
  )
  
)
