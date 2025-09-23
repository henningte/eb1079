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


#### for static branching ####

# define the models to compute
irp_d_model_info <- irp_make_d_model_info(n_preprocessing = 3L)


# sub-pipelines
irp_fit_1_map <- 
  tarchetypes::tar_map(
    values = 
      list(
        cur_id_model = irp_d_model_info$id_model, 
        cur_target_variable = irp_d_model_info$target_variable, 
        cur_id_preprocessing = irp_d_model_info$id_preprocessing
      ),
    names = dplyr::all_of("cur_id_model"),
    tar_target(
      irp_d_model_info_enriched_2,
      command = {
        irp_d_model_info_enriched_1 |>
          dplyr::filter(id_model == cur_id_model)
      },
      resources =
        tar_resources(
          crew = tar_resources_crew(controller = "crew_parallel_2")
        )
    ),
    tar_target(
      irp_data_for_brms_models_1,
      command = irp_prepare_data_for_brms_models_1(irp_d_model_info_enriched_2),
      resources =
        tar_resources(
          crew = tar_resources_crew(controller = "crew_parallel_2")
        )
    ),
    tar_target(
      irp_fit_1_pathfinder,
      command = 
        irp_make_fit_1_pathfinder(
          x = irp_d_model_info_enriched_2, 
          irp_brms_compiled_model = irp_brms_compiled_model, 
          irp_data_for_brms_models_1 = irp_data_for_brms_models_1, 
          irp_isotope_standards_isotope_fraction = irp_isotope_standards_isotope_fraction, 
          irp_mcmc_settings = irp_mcmc_settings_horseshoe
        ),
      resources =
        tar_resources(
          crew = tar_resources_crew(controller = "crew_parallel_2")
        )
    ),
    tar_target(
      irp_fit_1,
      command = 
        irp_make_fit_1(
          x = irp_d_model_info_enriched_2, 
          irp_brms_compiled_model = irp_brms_compiled_model, 
          irp_data_for_brms_models_1 = irp_data_for_brms_models_1, 
          irp_isotope_standards_isotope_fraction = irp_isotope_standards_isotope_fraction, 
          irp_mcmc_settings = irp_mcmc_settings_horseshoe,
          irp_fit_1_pathfinder = irp_fit_1_pathfinder
        ),
      resources =
        tar_resources(
          crew = tar_resources_crew(controller = "crew_parallel_2")
        )
    ),
    tar_target( #5
      irp_fit_1_elpd,
      command = 
        irp_make_validation_elpd(
          irp_d_model_info_enriched_2, 
          irp_fit_1
        ),
      resources =
        tar_resources(
          crew = tar_resources_crew(controller = "crew_parallel_2")
        )
    ),
    tar_target(
      irp_fit_1_evaluation_1,
      command = 
        irp_fit_1_make_evaluation_1(
          irp_d_model_info_enriched_2 = irp_d_model_info_enriched_2, 
          irp_fit_1 = irp_fit_1, 
          irp_isotope_standards_isotope_fraction = irp_isotope_standards_isotope_fraction
        ),
      resources =
        tar_resources(
          crew = tar_resources_crew(controller = "crew_parallel_2")
        )
    ),
    tar_target(
      irp_fit_1_pd,
      command =
        irp_make_fit_1_pd(
          irp_d_model_info_enriched_2 = irp_d_model_info_enriched_2,
          irp_data_model_preprocessed = irp_data_model_preprocessed
        ),
      resources =
        tar_resources(
          crew = tar_resources_crew(controller = "crew_parallel_2")
        )
    ),
    tar_target(
      irp_fit_1_irpeat_model_config,
      command = 
        irp_make_mir_model_config(
          irp_mir_preprocessing_settings = irp_mir_preprocessing_settings, 
          irp_d_model_info = irp_d_model_info_enriched_2, 
          prediction_domain = irp_fit_1_pd
        ),
      resources =
        tar_resources(
          crew = tar_resources_crew(controller = "crew_parallel_2")
        )
    ),
    tar_target(
      irp_pmird_gap_filled_helper_1,
      command = 
        irp_make_pmird_gap_filled_helper_1(
          irp_d_model_info_enriched_2 = irp_d_model_info_enriched_2,
          irp_fit_1_irpeat_model_config = irp_fit_1_irpeat_model_config, 
          irp_id_model_best = irp_id_model_best, 
          irp_data_model_preprocessed = irp_data_model_preprocessed
        ),
      resources =
        tar_resources(
          crew = tar_resources_crew(controller = "crew_parallel_1")
        )
    ),
    tar_target( #10
      irp_fit_1_slopes,
      command = 
        if(! irp_d_model_info_enriched_2$id_model %in% irp_id_model_best) {
          NULL
        } else {
          irp_make_fit_1_slopes(
            irp_fit_1 = irp_fit_1, 
            irp_d_model_info_enriched_2 = irp_d_model_info_enriched_2, 
            irp_data_model_preprocessed = irp_data_model_preprocessed
          )
        },
      format = "file",
      resources =
        tar_resources(
          crew = tar_resources_crew(controller = "crew_parallel_2")
        )
    )
  )

# targets
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
    irp_data_evaluation_literature_values_file,
    command = "data/raw_data/literature_review_predictive_accuracy.csv",
    format = "file"
  ),
  tar_target(
    irp_data_evaluation_literature_values,
    command = irp_make_data_evaluation_literature_values(irp_data_evaluation_literature_values_file = irp_data_evaluation_literature_values_file)
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
  tar_target(
    irp_data_wittington2024,
    command = irp_make_data_wittington2024()
  ),
  #### settings ####
  tar_target(
    irp_mcmc_settings,
    command = irp_get_mcmc_settings()
  ),
  tar_target(
    irp_mcmc_settings_horseshoe,
    command = irp_get_mcmc_settings_horseshoe()
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
    irp_n_training_sample_max,
    command = 100L
  ),
  tar_target(
    irp_target_variables,
    command = irp_get_target_variables()
  ),
  tar_target(
    irp_d_model_info_enriched_1,
    command = {
      irp_d_model_info |>
        irp_add_filter_criteria_to_d_model_info() |>
        irp_add_y_has_error_to_d_model_info() |>
        irp_add_likelihood_to_d_model_info() |>
        irp_add_model_formula_to_d_model_info() |>
        irp_add_id_measurement_all_to_d_model_info(
          irp_data_model_preprocessed = irp_data_model_preprocessed
        ) |>
        irp_add_n_training_sample_to_d_model_info(
          irp_n_training_sample_max = irp_n_training_sample_max
        ) |>
        irp_add_data_partition_to_d_model_info(
          irp_data_model_preprocessed = irp_data_model_preprocessed
        ) |>
        irp_add_y_to_d_model_info(
          irp_data_model_preprocessed = irp_data_model_preprocessed, 
          irp_isotope_standards_isotope_fraction = irp_isotope_standards_isotope_fraction, 
          irp_pmird_units = irp_pmird_units, 
          irp_m1 = irp_m1, 
          irp_m2 = irp_m2, 
          irp_d_compounds_standard_enthalpies_of_formation = irp_d_compounds_standard_enthalpies_of_formation
        ) |>
        irp_add_y_scale_to_d_model_info(
          irp_isotope_standards_isotope_fraction = irp_isotope_standards_isotope_fraction, 
          irp_mass_density_fe = irp_mass_density_fe
        ) |>
        irp_add_x_train_to_d_model_info(
          irp_data_model_preprocessed = irp_data_model_preprocessed
        ) |>
        irp_add_predictor_scale_to_d_model_info() |>
        irp_add_priors_to_d_model_info(
          irp_isotope_standards_isotope_fraction = irp_isotope_standards_isotope_fraction
        ) |>
        irp_add_rng_seed_to_d_model_info()
    }
  ),
  irp_fit_1_map,
  tar_target(
    irp_brms_compiled_model_metadata,
    command = 
      irp_d_model_info_enriched_1 |>
      dplyr::filter(! duplicated(paste0(likelihood_name, "_", y_has_error))) |>
      dplyr::select(likelihood_name, y_has_error, brms_formula, likelihood, prior_all, prior_stanvars) |>
      dplyr::arrange(likelihood_name, y_has_error) |>
      dplyr::mutate(
        prior_all =
          purrr::map(prior_all, function(.x) {
            .x |>
              dplyr::mutate(
                prior =
                  dplyr::case_when(
                    ! stringr::str_detect(prior, "^horseshoe") ~ prior,
                    TRUE ~ stringr::str_replace(prior, "^horseshoe\\(df = \\d+", "horseshoe(df = 1")
                  )
              )
          })
      )
  ),
  tar_target(
    irp_brms_compiled_model,
    command =
      irp_make_brms_compiled_model(
        irp_brms_compiled_model_metadata = irp_brms_compiled_model_metadata
      )
  ),
  tar_combine(
    irp_fit_1_map_elpd,
    irp_fit_1_map[[5]],
    command = c(list(!!!.x))
  ),
  tar_target(
    irp_fit_1_map_elpd_compare,
    command =
      irp_make_fit_1_map_elpd_compare(
        irp_d_model_info_enriched_1 = irp_d_model_info_enriched_1,
        irp_fit_1_map_elpd = irp_fit_1_map_elpd
      )
  ),
  tar_combine(
    irp_fit_1_map_evaluation_1,
    irp_fit_1_map[[6]],
    command = c(list(!!!.x))
  ),
  tar_combine(
    irp_fit_1_map_pd,
    irp_fit_1_map[[7]],
    command = c(list(!!!.x))
  ),
  tar_combine(
    irp_fit_1_map_slopes,
    irp_fit_1_map[[10]],
    command = {
      res <- c(list(!!!.x))
      res[purrr::map_lgl(res, length) == 1L]
    } 
  ),
  #### for irpeat ####
  tar_target(
    irp_id_model_best,
    command = {
      irp_fit_1_map_elpd_compare |>
        dplyr::left_join(
          tibble::tibble(
            id_model =  
              purrr::map_chr(irp_fit_1_map_evaluation_1, function(.x) {
                .x$yhat$id_model[[1]]
              }),
            num_divergent =
              purrr::map_int(irp_fit_1_map_evaluation_1, function(.x) {
                .x$num_divergent
              })
          ),
          by = "id_model"
        ) |>
        dplyr::filter(num_divergent == 0) |>
        dplyr::arrange(target_variable, dplyr::desc(elpd_diff)) |>
        dplyr::filter(! duplicated(target_variable)) |>
        dplyr::pull(id_model)
    }
  ),
  #### gap filling ####
  tar_combine(
    irp_pmird_gap_filled_helper_1_combined,
    irp_fit_1_map[[9]],
    command = c(list(!!!.x))
  ),
  tar_target(
    irp_pmird_gap_filled,
    command = 
      irp_make_pmird_gap_filled(
        irp_pmird_mirs = irp_pmird_mirs,
        irp_pmird_gap_filled_helper_1_combined = irp_pmird_gap_filled_helper_1_combined, 
        irp_m1 = irp_m1, 
        irp_m2 = irp_m2, 
        irp_d_compounds_standard_enthalpies_of_formation = irp_d_compounds_standard_enthalpies_of_formation, 
        irp_pmird_units = irp_pmird_units
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
  ),
  tar_target(
    irp_plot_14_1,
    command = 
      irp_make_plot_14(
        irp_fit_1_map_evaluation_1 = irp_fit_1_map_evaluation_1, 
        irp_d_model_info_enriched_1 = irp_d_model_info_enriched_1, 
        irp_fit_1_map_elpd_compare = irp_fit_1_map_elpd_compare,
        irp_pmird_mirs = irp_pmird_mirs, 
        variable_color = sym("is_training_data"), 
        variable_color_legend_title = "",
        highlight_outliers = TRUE,
        file_plot = "figures/irp_plot_14_1.pdf"
      ),
    format = "file"
  ),
  tar_target(
    irp_plot_14_2,
    command = 
      irp_make_plot_14(
        irp_fit_1_map_evaluation_1 = irp_fit_1_map_evaluation_1, 
        irp_d_model_info_enriched_1 = irp_d_model_info_enriched_1, 
        irp_fit_1_map_elpd_compare = irp_fit_1_map_elpd_compare,
        irp_pmird_mirs = irp_pmird_mirs, 
        variable_color = sym("Ca"), 
        variable_color_legend_title = "Calcium content (&mu;g g<sup>-1</sup>)",
        highlight_outliers = TRUE,
        file_plot = "figures/irp_plot_14_2.pdf"
      ),
    format = "file"
  ),
  tar_target(
    irp_plot_14_3,
    command = 
      irp_make_plot_14(
        irp_fit_1_map_evaluation_1 = irp_fit_1_map_evaluation_1, 
        irp_d_model_info_enriched_1 = irp_d_model_info_enriched_1, 
        irp_fit_1_map_elpd_compare = irp_fit_1_map_elpd_compare,
        irp_pmird_mirs = irp_pmird_mirs, 
        variable_color = sym("N"), 
        variable_color_legend_title = "Nitrogen content (g g<sup>-1</sup>)",
        highlight_outliers = TRUE,
        file_plot = "figures/irp_plot_14_3.pdf"
      ),
    format = "file"
  ),
  tar_target(
    irp_plot_14_4,
    command = 
      irp_make_plot_14(
        irp_fit_1_map_evaluation_1 = irp_fit_1_map_evaluation_1, 
        irp_d_model_info_enriched_1 = irp_d_model_info_enriched_1, 
        irp_fit_1_map_elpd_compare = irp_fit_1_map_elpd_compare,
        irp_pmird_mirs = irp_pmird_mirs, 
        variable_color = sym("Si"), 
        variable_color_legend_title = "Silicon content (&mu;g g<sup>-1</sup>)",
        highlight_outliers = TRUE,
        file_plot = "figures/irp_plot_14_4.pdf"
      ),
    format = "file"
  ),
  tar_target(
    irp_plot_15,
    command = 
      irp_make_plot_15(
        irp_data_evaluation_literature_values = irp_data_evaluation_literature_values,
        irp_fit_1_map_evaluation_1 = irp_fit_1_map_evaluation_1, 
        irp_d_model_info_enriched_1 = irp_d_model_info_enriched_1, 
        irp_fit_1_map_elpd_compare = irp_fit_1_map_elpd_compare, 
        file_plot = "figures/irp_plot_15.pdf"
      ),
    format = "file"
  ),
  tar_target(
    irp_plot_16_1,
    command =
      irp_make_plot_16(
        irp_fit_1_map_evaluation_1 = irp_fit_1_map_evaluation_1, 
        irp_d_model_info_enriched_1 = irp_d_model_info_enriched_1, 
        irp_fit_1_map_elpd_compare = irp_fit_1_map_elpd_compare, 
        irp_pmird_mirs = 
          irp_pmird_mirs |>
          dplyr::mutate(
            Ca = Ca/1000000
          ),
        x_variable = sym("Ca"), 
        x_variable_axis_title = "Ca content g g<sup>-1</sup>", 
        highlight_outliers = TRUE, 
        file_plot = "figures/irp_plot_16_1.pdf"
      ),
    format = "file"
  ),
  tar_target(
    irp_plot_16_2,
    command =
      irp_make_plot_16(
        irp_fit_1_map_evaluation_1 = irp_fit_1_map_evaluation_1, 
        irp_d_model_info_enriched_1 = irp_d_model_info_enriched_1, 
        irp_fit_1_map_elpd_compare = irp_fit_1_map_elpd_compare, 
        irp_pmird_mirs = irp_pmird_mirs,
        x_variable = sym("N"), 
        x_variable_axis_title = "N content g g<sup>-1</sup>", 
        highlight_outliers = TRUE, 
        file_plot = "figures/irp_plot_16_2.pdf"
        ),
    format = "file"
  ),
  tar_target(
    irp_plot_16_3,
    command =
      irp_make_plot_16(
        irp_fit_1_map_evaluation_1 = irp_fit_1_map_evaluation_1, 
        irp_d_model_info_enriched_1 = irp_d_model_info_enriched_1, 
        irp_fit_1_map_elpd_compare = irp_fit_1_map_elpd_compare, 
        irp_pmird_mirs =
          irp_pmird_mirs |>
          dplyr::mutate(
            Si = Si/1000000
          ),
        x_variable = sym("Si"), 
        x_variable_axis_title = "Si content g g<sup>-1</sup>", 
        highlight_outliers = TRUE, 
        file_plot = "figures/irp_plot_16_3.pdf"
      ),
    format = "file"
  ),
  tar_target(
    irp_plot_17,
    command = 
      irp_make_plot_17(
        irp_data_model_preprocessed = irp_data_model_preprocessed, 
        file_plot = "figures/irp_plot_17.pdf"
      ),
    format = "file"
  ),
  tar_target(
    irp_plot_18,
    command =
      irp_make_plot_18(
        irp_fit_1_map_evaluation_1 = irp_fit_1_map_evaluation_1, 
        irp_d_model_info_enriched_1 = irp_d_model_info_enriched_1, 
        irp_fit_1_map_elpd_compare = irp_fit_1_map_elpd_compare, 
        irp_data_model_preprocessed = irp_data_model_preprocessed, 
        file_plot = "figures/irp_plot_18.pdf"
      )
  ),
  tar_target(
    irp_plot_19,
    command =
      irp_make_plot_19(
        irp_fit_1_map_evaluation_1 = irp_fit_1_map_evaluation_1, 
        irp_d_model_info_enriched_1 = irp_d_model_info_enriched_1, 
        irp_fit_1_map_elpd_compare = irp_fit_1_map_elpd_compare, 
        irp_fit_1_map_pd = irp_fit_1_map_pd, 
        irp_id_model_best = irp_id_model_best,
        file_plot = "figures/irp_plot_19.pdf"
      ),
    format = "file"
  ),
  tar_target(
    irp_plot_20,
    command =
      irp_make_plot_20(
        irp_fit_1_map_slopes = irp_fit_1_map_slopes,
        file_plot = "figures/irp_plot_20.pdf"
      ),
    format = "file"
  ),
  tar_target(
    irp_plot_21,
    command =
      irp_make_plot_21(
        irp_data_wittington2024 = irp_data_wittington2024,
        file_plot = "figures/irp_plot_21.pdf"
      ),
    format = "file"
  ),
  tar_target(
    irp_plot_22,
    command =
      irp_make_plot_22(
        irp_data_wittington2024 = irp_data_wittington2024,
        file_plot = "figures/irp_plot_22.pdf"
      ),
    format = "file"
  ),
  #### tables ####
  tar_target(
    irp_table_1,
    command = 
      irp_make_table_1(
        irp_fit_1_map_evaluation_1 = irp_fit_1_map_evaluation_1, 
        irp_d_model_info_enriched_1 = irp_d_model_info_enriched_1, 
        irp_fit_1_map_elpd_compare = irp_fit_1_map_elpd_compare,
        irp_id_model_best = irp_id_model_best
      )
  ),
  tar_target(
    irp_table_2,
    command = 
      irp_make_table_2(
        irp_pmird_mirs = irp_pmird_mirs, 
        irp_pmird_gap_filled = readRDS_rvars(irp_pmird_gap_filled) 
      )
  ),
  tar_target(
    irp_table_3,
    command = 
      irp_make_table_3(
        irp_fit_1_map_slopes = irp_fit_1_map_slopes
      )
  ),
  #### report ####
  tar_target(
    irp_references_file,
    command = "references.bib",
    format = "file"
  ),
  tar_target(
    irp_pmird_citations,
    command = irp_get_pmird_citations(irp_d_model_info_enriched_1 = irp_d_model_info_enriched_1)
  ),
  tar_render(
    irp_supporting_info,
    path = "irp-supporting-info.Rmd"
  ),
  tar_render(
    irp_paper,
    path = "irp-paper.Rmd"
  )
)
