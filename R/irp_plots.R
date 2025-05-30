#' Plots enthalpy of combustion versus electrons transferred during combustion from the OBIGT database along with predictions by `irp_m1`
#' 
#' @export
irp_make_plot_1 <- function(irp_db_obigt_organic, irp_m1, file_plot = "figures/irp_plot_1.pdf") {
  
  res <- 
    tibble::tibble(
      E_per_C = 
        seq(min(irp_db_obigt_organic$E_per_C), max(irp_db_obigt_organic$E_per_C), length.out = 20)
    )
  
  res <- 
    res |>
    dplyr::mutate(
      yhat =
        irp_m1$res_model |>
        brms::posterior_predict(
          type = "response",
          newdata = 
            res |>
            dplyr::mutate(
              E_per_C = 
                E_per_C |>
                magrittr::divide_by(irp_m1$res_config$x_max) |>
                as.numeric()
            )
        ) |>
        posterior::rvar() |>
        magrittr::multiply_by(as.numeric(irp_m1$res_config$y_min))
    )
  
  res_plot <- 
    irp_db_obigt_organic |>
    ggplot(aes(x = as.numeric(E_per_C))) +
    geom_point(aes(y = as.numeric(enthalpy_of_combustion_per_C))) +
    geom_path(res, mapping = aes(y = mean(yhat)), color = "grey50") +
    labs(
      y = "&Delta;H<sub>C</sub><sup>0</sup> (kJ mol<sub>C</sub><sup>-1</sup>)",
      x = "Electrons transferred during combustion (mol<sub>e<sup>-</sup></sub> mol<sub>C</sub><sup>-1</sup>)"
    ) +
    scale_y_continuous(labels = function(.x) .x/1000) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = ggtext::element_markdown(),
      axis.title.x = ggtext::element_markdown(size = 11),
      axis.title.y = ggtext::element_markdown(size = 11)
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 5, height = 5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#' Plots enthalpy of combustion predicted by `irp_m1` versus enthalpy of combustion from the OBIGT database
#' 
#' @export
irp_make_plot_2 <- function(irp_m1, file_plot = "figures/irp_plot_2.pdf") {
  
  res <- 
    irp_m1$res_model$data |>
    dplyr::rename(
      y = "enthalpy_of_combustion_per_C"
    ) |>
    dplyr::mutate(
      y = y * irp_m1$res_config$y_min,
      yhat = 
        irp_m1$res_model |>
        brms::posterior_predict(
          type = "response"
        ) |>
        posterior::rvar() |>
        magrittr::multiply_by(as.numeric(irp_m1$res_config$y_min))
    )
   
  res_plot <- 
    res |>
    ggplot(aes(y = as.numeric(y), xdist = yhat)) +
    ggdist::stat_pointinterval(.width = 0.95, interval_colour = "grey", linewidth = 0.5) +
    geom_abline(intercept = 0, slope = 1, color = "grey50") +
    labs(
      y = "OBIGT &Delta;H<sub>C</sub><sup>0</sup> (kJ mol<sub>C</sub><sup>-1</sup>)",
      x = "Predicted &Delta;H<sub>C</sub><sup>0</sup> (kJ mol<sub>C</sub><sup>-1</sup>)"
    ) +
    scale_y_continuous(labels = function(.x) .x/1000) +
    scale_x_continuous(labels = function(.x) .x/1000) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = ggtext::element_markdown(),
      axis.title.x = ggtext::element_markdown(size = 11),
      axis.title.y = ggtext::element_markdown(size = 11)
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 5, height = 5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}



#' Plots enthalpy of combustion versus electrons transferred during combustion from the OBIGT database along with predictions by `irp_m1`
#' 
#' @export
irp_make_plot_3 <- function(irp_m2, file_plot = "figures/irp_plot_3.pdf") {
  
  res_data <- irp_m2$res_data
  
  res <- 
    tibble::tibble(
      sum_of_element_standard_entropies_per_C = 
        seq(min(res_data$sum_of_element_standard_entropies_per_C), max(res_data$sum_of_element_standard_entropies_per_C), length.out = 20)
    )
  
  res <- 
    res |>
    dplyr::mutate(
      yhat =
        irp_m2$res_model |>
        brms::posterior_predict(
          type = "response",
          newdata = 
            res |>
            dplyr::mutate(
              sum_of_element_standard_entropies_per_C = 
                sum_of_element_standard_entropies_per_C |>
                magrittr::divide_by(irp_m2$res_config$x_max) |>
                as.numeric()
            )
        ) |>
        posterior::rvar() |>
        magrittr::multiply_by(as.numeric(irp_m2$res_config$y_min))
    )
  
  res_plot <- 
    res_data |>
    dplyr::mutate(
      data_type = 
        dplyr::case_when(
          is.na(file) ~ "Battley (1999)",
          TRUE ~ "OBIGT database"
        )
    ) |>
    ggplot(aes(x = as.numeric(sum_of_element_standard_entropies_per_C))) +
    geom_point(aes(y = as.numeric(entropy_of_formation_per_C), color = data_type)) +
    geom_path(res, mapping = aes(y = mean(yhat)), color = "grey50") +
    labs(
      y = "&Delta;S<sub>f</sub><sup>0</sup> (J K<sup>-1</sup> mol<sub>C</sub><sup>-1</sup>)",
      x = "Sum of standard atomic entropies (J K<sup>-1</sup> mol<sub>C</sub><sup>-1</sup>)"
    ) +
    guides(
      color = guide_legend(title = "Data source", override.aes = list(size = 3))
    ) +
    scale_color_manual(values = c("grey", "grey5")) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = ggtext::element_markdown(),
      axis.title.x = ggtext::element_markdown(size = 11),
      axis.title.y = ggtext::element_markdown(size = 11)
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 5, height = 5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#' Plots entropy of formation predicted by `irp_m2` versus entropy of formation from the OBIGT database and Battley (1999)
#' 
#' @export
irp_make_plot_4 <- function(irp_m2, file_plot = "figures/irp_plot_4.pdf") {
  
  res <- 
    irp_m2$res_model$data |>
    dplyr::rename(
      y = "entropy_of_formation_per_C"
    ) |>
    dplyr::mutate(
      y = y * irp_m2$res_config$y_min,
      yhat = 
        irp_m2$res_model |>
        brms::posterior_predict(
          type = "response"
        ) |>
        posterior::rvar() |>
        magrittr::multiply_by(as.numeric(irp_m2$res_config$y_min))
    )
  
  res_plot <- 
    res |>
    dplyr::mutate(
      data_type = 
        dplyr::case_when(
          is.na(irp_m2$res_data$file) ~ "Battley (1999)",
          TRUE ~ "OBIGT database"
        )
    ) |>
    ggplot(aes(y = as.numeric(y), xdist = yhat)) +
    ggdist::stat_pointinterval(aes(color = data_type), .width = 0.95, interval_colour = "grey", linewidth = 0.5) +
    geom_abline(intercept = 0, slope = 1, color = "grey50") +
    labs(
      y = "OBIGT or Battley (1999) &Delta;S<sub>f</sub><sup>0</sup> (J K<sup>-1</sup> mol<sub>C</sub><sup>-1</sup>)",
      x = "Predicted &Delta;S<sub>f</sub><sup>0</sup> (J K<sup>-1</sup> mol<sub>C</sub><sup>-1</sup>)"
    ) +
    guides(
      color = guide_legend(title = "Data source", override.aes = list(size = 3))
    ) +
    scale_color_manual(values = c("grey", "grey5")) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = ggtext::element_markdown(),
      axis.title.x = ggtext::element_markdown(size = 11),
      axis.title.y = ggtext::element_markdown(size = 11)
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 5, height = 5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}



#' Plots predicted standard Gibbs free energy of formation versus values in the OBIGT database
#' 
#' @export
irp_make_plot_5 <- function(irp_db_obigt_organic, irp_m1, irp_m2, irp_d_compounds_standard_enthalpies_of_formation, file_plot = "figures/irp_plot_5.pdf") {
  
  res_cf <- 
    irp_db_obigt_organic |>
    dplyr::filter(state == "cr") |>
    dplyr::pull(formula) |>
    irp_cf_to_df() |>
    dplyr::mutate(
      dplyr::across(
        dplyr::everything(),
        function(.x) {
          units::set_units(.x, paste0("mol_", dplyr::cur_column()), mode = "standard")
        }
      )
    )
  
  res <- 
    res_cf |>
    irp_estimate_dgf0_per_C(
      irp_m1 = irp_m1, 
      irp_m2 = irp_m2, 
      irp_d_compounds_standard_enthalpies_of_formation = irp_d_compounds_standard_enthalpies_of_formation
    )
  
  res_plot <- 
    irp_db_obigt_organic |>
    dplyr::filter(state == "cr") |>
    dplyr::mutate(
      dgf0_per_C = 
        res$dgf0_per_C |>
        irp_set_units_rvar("kJ/mol_C", mode = "standard") |>
        irp_drop_units_rvar(),
      y = 
        G |>
        units::set_units(paste0(unique(E_units), "/mol_sample"), mode = "standard") |>
        units::set_units("kJ/mol_sample", mode = "standard") |>
        magrittr::divide_by(res_cf$C) |>
        as.numeric()
    ) |>
    ggplot(aes(y = y, xdist = dgf0_per_C)) +
    ggdist::stat_pointinterval(.width = 0.95, interval_colour = "grey", linewidth = 0.5) +
    geom_abline(intercept = 0, slope = 1, color = "grey50") +
    labs(
      y = "OBIGT &Delta;G<sub>f</sub><sup>0</sup> (kJ mol<sub>C</sub><sup>-1</sup>)",
      x = "Predicted &Delta;G<sub>f</sub><sup>0</sup> (kJ mol<sub>C</sub><sup>-1</sup>)"
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = ggtext::element_markdown(),
      axis.title.x = ggtext::element_markdown(size = 11),
      axis.title.y = ggtext::element_markdown(size = 11)
    ) +
    coord_fixed()
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 5, height = 5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#' Plots volume fractions predicted by `irp_m3` versus measured volume fractions from the training data 
#' 
#' @export
irp_make_plot_6 <- function(irp_m3, file_plot = "figures/irp_plot_6.pdf") {
  
  res <- 
    irp_m3$res_data |>
    irp_predict_with_m3(irp_m3 = irp_m3) |>
    dplyr::mutate(
      all_pores = 1 - volume_fraction_solids,
      all_pores_1 = irp_set_units_rvar(posterior::rvar(1), "L/L", mode = "standard") - volume_fraction_solids_1
    ) |>
    tidyr::nest(
      macroporosity = c(macroporosity, macroporosity_1),
      non_macroporosity = c(non_macroporosity, non_macroporosity_1),
      volume_fraction_solids = c(volume_fraction_solids, volume_fraction_solids_1),
      all_pores = c(all_pores, all_pores_1)
    ) |>
    tidyr::pivot_longer(
      dplyr::all_of(c("all_pores", "macroporosity", "non_macroporosity", "volume_fraction_solids")),
      names_to = "variable",
      values_to = "value"
    ) |>
    dplyr::mutate(
      value =
        purrr::map(value, function(.x) {
          .x |>
            setNames(nm = c("y", "yhat"))
        })
    ) |>
    tidyr::unnest(value)
  
  res_plot <- 
    res |>
    dplyr::mutate(
      data_type = 
        dplyr::case_when(
          is.na(id_measurement) ~ "Organic shale",
          TRUE ~ "Peat"
        ),
      variable =
        dplyr::case_when(
          variable == "all_pores" ~ "All pores",
          variable == "macroporosity" ~ "Macropores",
          variable == "non_macroporosity" ~ "Non-macropores",
          variable == "volume_fraction_solids" ~ "Solids"
        )
    ) |>
    ggplot(aes(y = as.numeric(y), xdist = yhat)) +
    ggdist::stat_pointinterval(aes(fill = data_type), .width = 0.95, interval_colour = "grey", linewidth = 0.2, shape = 21, point_size = 1.5) +
    geom_abline(intercept = 0, slope = 1, color = "grey50") +
    labs(
      y = "Volume fraction (L L<sup>-1</sup><sub>sample</sub>)",
      x = "Predicted volume fraction (L L<sup>-1</sup><sub>sample</sub>)"
    ) +
    guides(
      fill = guide_legend(title = "Sample type", override.aes = list(size = 3))
    ) +
    scale_fill_manual(values = c("white", "grey5")) +
    facet_wrap(~ variable) +
    coord_fixed() +
    theme_classic() +
    theme(
      legend.position = "bottom",
      axis.title.x = ggtext::element_markdown(size = 14),
      axis.title.y = ggtext::element_markdown(size = 14),
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 13),
      strip.text.y = ggtext::element_markdown(size = 13),
      legend.text = ggtext::element_markdown(size = 13),
      legend.title = ggtext::element_markdown(size = 13)
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 5, height = 5.5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#' Plots volume fractions predicted by `irp_m3` versus bulk density
#' 
#' @export
irp_make_plot_7 <- function(irp_m3, file_plot = "figures/irp_plot_7.pdf") {
  
  res <- 
    tibble::tibble(
      bulk_density = seq(from = 0.00001, to = 3.5, by = 0.05)
    ) |>
    irp_predict_with_m3(irp_m3 = irp_m3) |>
    dplyr::mutate(
      all_pores_1 = irp_set_units_rvar(posterior::rvar(1), "L/L", mode = "standard") - volume_fraction_solids_1
    ) |>
    tidyr::pivot_longer(
      dplyr::all_of(paste0(c("all_pores", "macroporosity", "non_macroporosity", "volume_fraction_solids"), "_1")),
      names_to = "variable",
      values_to = "value"
    ) |>
    dplyr::mutate(
      variable =
        dplyr::case_when(
          variable == "all_pores_1" ~ "All pores",
          variable == "macroporosity_1" ~ "Macropores",
          variable == "non_macroporosity_1" ~ "Non-macropores",
          variable == "volume_fraction_solids_1" ~ "Solids"
        )
    )
  
  res_data <- 
    irp_m3$res_data |>
    dplyr::mutate(
      all_pores = 1 - volume_fraction_solids,
    ) |>
    tidyr::pivot_longer(
      dplyr::all_of(c("all_pores", "macroporosity", "non_macroporosity", "volume_fraction_solids")),
      names_to = "variable",
      values_to = "value"
    ) |>
    dplyr::mutate(
      data_type = 
        dplyr::case_when(
          is.na(id_measurement) ~ "Organic shale",
          TRUE ~ "Peat"
        ),
      variable =
        dplyr::case_when(
          variable == "all_pores" ~ "All pores",
          variable == "macroporosity" ~ "Macropores",
          variable == "non_macroporosity" ~ "Non-macropores",
          variable == "volume_fraction_solids" ~ "Solids"
        )
    )
  
  res_plot <- 
    res |>
    ggplot(aes(x = bulk_density)) +
    ggdist::stat_lineribbon(aes(ydist = value), linewidth = 0.5) +
    scale_fill_brewer() +
    guides(
      fill = guide_legend(title = "Significance level", override.aes = list(size = 3, color = NA))
    ) +
    ggnewscale::new_scale_fill() +
    geom_point(data = res_data, aes(y = value, fill = data_type), shape = 21) +
    scale_fill_manual(values = c("white", "grey5")) +
    guides(
      fill = guide_legend(title = "Sample type", override.aes = list(size = 3))
    ) +
    labs(
      y = "Measured or predicted volume fraction (L L<sup>-1</sup><sub>sample</sub>)",
      x = "Bulk density (g cm<sup>-3</sup>)"
    ) +
    facet_wrap(~ variable) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      axis.title.x = ggtext::element_markdown(size = 14),
      axis.title.y = ggtext::element_markdown(size = 14),
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 13),
      strip.text.y = ggtext::element_markdown(size = 13),
      legend.text = ggtext::element_markdown(size = 13),
      legend.title = ggtext::element_markdown(size = 13)
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 7.5, height = 5.5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}



#' Plots saturated hydraulic conductivity predicted by `irp_m4` versus measured saturated hydraulic conductivity
#' 
#' @export
irp_make_plot_8 <- function(irp_m4, file_plot = "figures/irp_plot_8.pdf") {
  
  res <- 
    irp_m4$res_data |>
    irp_predict_with_m4(irp_m4 = irp_m4) |>
    dplyr::mutate(
      saturated_hydraulic_conductivity = hydraulic_conductivity,
      saturated_hydraulic_conductivity_log10 = log10(saturated_hydraulic_conductivity),
      saturated_hydraulic_conductivity_1_log10 = log10(saturated_hydraulic_conductivity_1)
    ) |>
    tidyr::nest(
      saturated_hydraulic_conductivity = c(saturated_hydraulic_conductivity, saturated_hydraulic_conductivity_1),
      saturated_hydraulic_conductivity_log10 = c(saturated_hydraulic_conductivity_log10, saturated_hydraulic_conductivity_1_log10)
    ) |>
    tidyr::pivot_longer(
      dplyr::all_of(c("saturated_hydraulic_conductivity", "saturated_hydraulic_conductivity_log10")),
      names_to = "variable",
      values_to = "value"
    ) |>
    dplyr::mutate(
      value =
        purrr::map(value, function(.x) {
          .x |>
            setNames(nm = c("y", "yhat"))
        })
    ) |>
    tidyr::unnest(value)
  
  res_plot <- 
    res |>
    dplyr::mutate(
      variable =
        dplyr::case_when(
          variable == "saturated_hydraulic_conductivity" ~ "k<sub>s</sub> (cm/h)",
          variable == "saturated_hydraulic_conductivity_log10" ~ "log<sub>10</sub>(k<sub>s</sub>)"
        )
    ) |>
    ggplot(aes(y = as.numeric(y), xdist = yhat)) +
    ggdist::stat_pointinterval(.width = 0.95, interval_colour = "grey", linewidth = 0.2, point_size = 1.5) +
    geom_abline(intercept = 0, slope = 1, color = "grey50") +
    labs(
      y = "Hydraulic conductivity (cm h<sup>-1</sup> or log<sub>10</sub>-scaled)",
      x = "Predicted hydraulic conductivity (cm h<sup>-1</sup> or log<sub>10</sub>-scaled)"
    ) +
    facet_wrap(~ variable, scales = "free") +
    theme_classic() +
    theme(
      legend.position = "bottom",
      axis.title.x = ggtext::element_markdown(size = 14),
      axis.title.y = ggtext::element_markdown(size = 14),
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 13),
      strip.text.y = ggtext::element_markdown(size = 13),
      legend.text = ggtext::element_markdown(size = 13),
      legend.title = ggtext::element_markdown(size = 13)
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 7.5, height = 5.5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}



#' Plots saturated hydraulic conductivity predicted by `irp_m4` versus bulk density
#' 
#' @export
irp_make_plot_9 <- function(irp_m4, file_plot = "figures/irp_plot_9.pdf") {
  
  res <- 
    tibble::tibble(
      bulk_density = seq(from = 0.00001, to = 0.9, by = 0.05)
    ) |>
    irp_predict_with_m4(irp_m4 = irp_m4) |>
    dplyr::mutate(
      saturated_hydraulic_conductivity_1_log10 = log10(saturated_hydraulic_conductivity_1)
    ) |>
    tidyr::pivot_longer(
      dplyr::all_of(paste0(c("saturated_hydraulic_conductivity", "saturated_hydraulic_conductivity"), "_1", c("", "_log10"))),
      names_to = "variable",
      values_to = "value"
    ) |>
    dplyr::mutate(
      variable =
        dplyr::case_when(
          variable == "saturated_hydraulic_conductivity_1" ~ "k<sub>s</sub> (cm/h)",
          variable == "saturated_hydraulic_conductivity_1_log10" ~ "log<sub>10</sub>(k<sub>s</sub>)"
        )
    )
  
  res_data <- 
    irp_m4$res_data |>
    dplyr::mutate(
      saturated_hydraulic_conductivity = hydraulic_conductivity,
      saturated_hydraulic_conductivity_log10 = log10(saturated_hydraulic_conductivity),
    ) |>
    tidyr::pivot_longer(
      dplyr::all_of(c("saturated_hydraulic_conductivity", "saturated_hydraulic_conductivity_log10")),
      names_to = "variable",
      values_to = "value"
    ) |>
    dplyr::mutate(
      variable =
        dplyr::case_when(
          variable == "saturated_hydraulic_conductivity" ~ "k<sub>s</sub> (cm/h)",
          variable == "saturated_hydraulic_conductivity_log10" ~ "log<sub>10</sub>(k<sub>s</sub>)"
        )
    )
  
  res_plot <- 
    res |>
    ggplot(aes(x = bulk_density)) +
    ggdist::stat_lineribbon(aes(ydist = value), linewidth = 0.5) +
    scale_fill_brewer() +
    guides(
      fill = guide_legend(title = "Significance level", override.aes = list(size = 3, color = NA))
    ) +
    geom_point(data = res_data, aes(y = value)) +
    labs(
      y = "Measured or predicted hydraulic conductivity<br>(cm h<sup>-1</sup> or log<sub>10</sub>-scaled)",
      x = "Bulk density (g cm<sup>-3</sup>)"
    ) +
    facet_wrap(~ variable, scales = "free_y") +
    theme_classic() +
    theme(
      legend.position = "bottom",
      axis.title.x = ggtext::element_markdown(size = 14),
      axis.title.y = ggtext::element_markdown(size = 14),
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 13),
      strip.text.y = ggtext::element_markdown(size = 13),
      legend.text = ggtext::element_markdown(size = 13),
      legend.title = ggtext::element_markdown(size = 13)
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 8, height = 5.5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}



#' Plots specific heat capacity predicted by `irp_m5` versus measured specific heat capacity 
#' 
#' @export
irp_make_plot_10 <- function(irp_m5, file_plot = "figures/irp_plot_10.pdf") {
  
  res <- 
    irp_m5$res_data |>
    irp_predict_with_m5(
      irp_m5 = irp_m5, 
      re_formula = NULL,
      allow_new_levels = TRUE,
      sample_new_levels = "uncertainty"
    ) |>
    dplyr::rename(
      y = "specific_heat_capacity",
      yhat = "specific_heat_capacity_1"
    )
  
  res_plot <- 
    res |>
    dplyr::mutate(
      data_type = 
        dplyr::case_when(
          sample_type == "air" ~ "Air",
          sample_type == "peat" ~ "Peat"
        ),
      yhat = irp_drop_units_rvar(yhat)
    ) |>
    ggplot(aes(y = as.numeric(y), xdist = yhat)) +
    geom_errorbar(
      aes(
        ymin = as.numeric(y - specific_heat_capacity_error),
        ymax = as.numeric(y + specific_heat_capacity_error),
        x = mean(yhat)
      ),
      color = "grey",
      width = 0.0
    ) +
    ggdist::stat_pointinterval(aes(fill = data_type), .width = 0.95, interval_colour = "grey", linewidth = 0.5, shape = 21) +
    geom_abline(intercept = 0, slope = 1, color = "grey50") +
    labs(
      y = "Measured specific heat capacity (J K<sup>-1</sup> g<sup>-1</sup>)",
      x = "Predicted specific heat capacity (J K<sup>-1</sup> g<sup>-1</sup>)"
    ) +
    guides(
      fill = guide_legend(title = "Sample type", override.aes = list(size = 3))
    ) +
    scale_fill_manual(values = c("white", "grey5")) +
    coord_fixed() +
    theme_classic() +
    theme(
      legend.position = "bottom",
      axis.title.x = ggtext::element_markdown(size = 14),
      axis.title.y = ggtext::element_markdown(size = 14),
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 13),
      strip.text.y = ggtext::element_markdown(size = 13),
      legend.text = ggtext::element_markdown(size = 13),
      legend.title = ggtext::element_markdown(size = 13)
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 5, height = 5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#' Plots specific heat capacity predicted by `irp_m5` or measured at 0°C versus N content 
#' 
#' @export
irp_make_plot_11 <- function(irp_m5, file_plot = "figures/irp_plot_11.pdf") {
  
  res <- 
    tibble::tibble(
      N = seq(from = 0, to = 0.04, by = 0.001),
      temperature = 0,
      specific_heat_capacity_error = mean(irp_m5$res_data$specific_heat_capacity_error)
    ) |>
    irp_predict_with_m5(
      irp_m5 = irp_m5,
      re_formula = NULL,
      allow_new_levels = TRUE,
      sample_new_levels = "uncertainty"
    ) |>
    dplyr::mutate(
      specific_heat_capacity_1 = irp_drop_units_rvar(specific_heat_capacity_1)
    )
  
  res_data <- 
    irp_m5$res_data |>
    dplyr::filter(temperature == 0)
  
  res_plot <- 
    res |>
    ggplot(aes(x = N)) +
    ggdist::stat_lineribbon(aes(ydist = specific_heat_capacity_1), linewidth = 0.5) +
    scale_fill_brewer() +
    guides(
      fill = guide_legend(title = "Significance level", override.aes = list(size = 3, color = NA))
    ) +
    geom_errorbar(
      data = res_data, 
      aes(
        ymin = as.numeric(specific_heat_capacity - specific_heat_capacity_error),
        ymax = as.numeric(specific_heat_capacity + specific_heat_capacity_error)
      ),
      color = "grey",
      width = 0.0
    ) +
    geom_point(data = res_data, aes(y = specific_heat_capacity)) +
    labs(
      y = "Measured or predicted specific heat capacity (J K<sup>-1</sup> g<sup>-1</sup>)",
      x = "N (g g<sup>-1</sup><sub>sample</sub>)"
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      axis.title.x = ggtext::element_markdown(size = 14),
      axis.title.y = ggtext::element_markdown(size = 14),
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 13),
      strip.text.y = ggtext::element_markdown(size = 13),
      legend.text = ggtext::element_markdown(size = 13),
      legend.title = ggtext::element_markdown(size = 13)
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 6, height = 6, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#' Plots dry thermal conductivity predicted by `irp_m5` versus measured dry thermal conductivity
#' 
#' @export
irp_make_plot_12 <- function(irp_m6, file_plot = "figures/irp_plot_12.pdf") {
  
  res <- 
    irp_m6$res_data |>
    irp_predict_with_m6(
      irp_m6 = irp_m6
    ) |>
    dplyr::rename(
      y = "dry_thermal_conductivity",
      yhat = "dry_thermal_conductivity_1"
    )
  
  res_plot <- 
    res |>
    dplyr::mutate(
      yhat = irp_drop_units_rvar(yhat)
    ) |>
    ggplot(aes(y = as.numeric(y), xdist = yhat)) +
    ggdist::stat_pointinterval(.width = 0.95, interval_colour = "grey", linewidth = 0.2) +
    geom_abline(intercept = 0, slope = 1, color = "grey50") +
    labs(
      y = "Measured dry thermal conductivity (W K<sup>-1</sup> m<sup>-1</sup>)",
      x = "Predicted specific heat capacity (W K<sup>-1</sup> m<sup>-1</sup>)"
    ) +
    coord_fixed() +
    theme_classic() +
    theme(
      legend.position = "bottom",
      axis.title.x = ggtext::element_markdown(size = 14),
      axis.title.y = ggtext::element_markdown(size = 14),
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 13),
      strip.text.y = ggtext::element_markdown(size = 13),
      legend.text = ggtext::element_markdown(size = 13),
      legend.title = ggtext::element_markdown(size = 13)
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 5, height = 5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#' Plots dry thermal conductivity predicted by `irp_m6` or measured versus bulk density
#' 
#' @export
irp_make_plot_13 <- function(irp_m6, file_plot = "figures/irp_plot_13.pdf") {
  
  res <- 
    tibble::tibble(
      bulk_density = seq(from = 0.00001, to = 2, by = 0.01)
    ) |>
    irp_predict_with_m6(
      irp_m6 = irp_m6,
    ) |>
    dplyr::mutate(
      dry_thermal_conductivity_1 = irp_drop_units_rvar(dry_thermal_conductivity_1)
    )
  
  res_data <- irp_m6$res_data
  
  res_plot <- 
    res |>
    ggplot(aes(x = bulk_density)) +
    ggdist::stat_lineribbon(aes(ydist = dry_thermal_conductivity_1), linewidth = 0.5) +
    scale_fill_brewer() +
    guides(
      fill = guide_legend(title = "Significance level", override.aes = list(size = 3, color = NA))
    ) +
    geom_point(data = res_data, aes(y = dry_thermal_conductivity)) +
    labs(
      y = "Measured or predicted dry thermal conductivity (W K<sup>-1</sup> m<sup>-1</sup>)",
      x = "bulk density (g cm<sup>-3</sup>)"
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      axis.title.x = ggtext::element_markdown(size = 14),
      axis.title.y = ggtext::element_markdown(size = 14),
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 13),
      strip.text.y = ggtext::element_markdown(size = 13),
      legend.text = ggtext::element_markdown(size = 13),
      legend.title = ggtext::element_markdown(size = 13)
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 6, height = 6.2, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}



#' Plots measured versus fittted/predicted values for the best models
#' 
#' @export
irp_make_plot_14 <- function(irp_fit_1_validation_elpd_diff, irp_fit_1_y_yhat) {
  
  res <- 
    irp_fit_1_validation_elpd_diff |>
    dplyr::arrange(dplyr::desc(is_best_model), dplyr::desc(ncomps)) |>
    dplyr::filter(! duplicated(target_variable))
  
  res <- 
    irp_fit_1_y_yhat |>
    dplyr::filter(id_model %in% res$id_model) |>
    dplyr::mutate(
      id_measurement = as.numeric(id_measurement)
    ) |>
    dplyr::left_join(
      irp_pmird_mirs |>
        dplyr::select(Si, S, Ca, Ti, P, C, N, bulk_density, id_measurement) |>
        dplyr::mutate(
          dplyr::across(
            dplyr::where(is.numeric), as.numeric
          )
        ),
      by = "id_measurement"
    )
  
  # irp_fit_1_y_yhat <- purrr::map(tar_read(irp_fit_1_y_yhat), readRDS) |> irp_do_call("rbind")
  
  res |>
    ggplot(aes(y = y, x = yhat_train)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "grey50") +
    facet_wrap(~ id_model, scales = "free") +
    labs(
      y = "Measured values",
      x = "Predicted values"
    ) +
    theme_classic() +
    theme(
      strip.background.x = element_blank(),
      strip.background.y = element_blank()
    )
  
}
















