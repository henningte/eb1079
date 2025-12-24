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
      legend.direction = "vertical",
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
    width = 7.5, height = 6, 
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
          variable == "saturated_hydraulic_conductivity" ~ "K<sub>s</sub> (cm/h)",
          variable == "saturated_hydraulic_conductivity_log10" ~ "log<sub>10</sub>(K<sub>s</sub>)"
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
          variable == "saturated_hydraulic_conductivity_1" ~ "K<sub>s</sub> (cm/h)",
          variable == "saturated_hydraulic_conductivity_1_log10" ~ "log<sub>10</sub>(K<sub>s</sub>)"
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
          variable == "saturated_hydraulic_conductivity" ~ "K<sub>s</sub> (cm/h)",
          variable == "saturated_hydraulic_conductivity_log10" ~ "log<sub>10</sub>(K<sub>s</sub>)"
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
      x = "Predicted dry thermal conductivity (W K<sup>-1</sup> m<sup>-1</sup>)"
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



#' Plots measured versus fitted/predicted values for the best models
#' 
#' @param highlight_outliers Logical value. If `TRUE`, selected outliers are 
#' highlighted with a different point shape.
#' 
#' @export
irp_make_plot_14 <- function(irp_fit_1_map_evaluation_1, irp_d_model_info_enriched_1, irp_fit_1_map_elpd_compare, irp_pmird_mirs, variable_color = sym("is_training_data"), variable_color_legend_title = "", highlight_outliers = FALSE, file_plot = "figures/irp_plot_14.pdf") {
  
  id_model_best <- 
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
    dplyr::arrange(id_model, dplyr::desc(elpd_diff)) |>
    dplyr::filter(! duplicated(target_variable)) |>
    dplyr::pull(id_model)
    
  cur_variable_color <- 
    if(variable_color == "is_training_data") {
      sym("id_measurement")
    } else {
      variable_color
    }
  
  res <- 
    irp_d_model_info_enriched_1 |>
    dplyr::select(id_model, target_variable) |>
    dplyr::mutate(
      validation =
        irp_fit_1_map_evaluation_1 |>
        purrr::map(function(.x) {
          .x$yhat |>
            dplyr::mutate( #---note: needed to avoid vctrs error
              dplyr::across(
                dplyr::where(is.integer),
                as.integer
              )
            ) |>
            dplyr::select(-id_model)
        })
    ) |>
    dplyr::filter(id_model %in% id_model_best) |>
    tidyr::unnest("validation") |>
    dplyr::left_join(
      irp_pmird_mirs |>
        dplyr::select(id_measurement, dplyr::all_of(cur_variable_color), taxon_rank_value),
      by = "id_measurement"
    ) |>
    dplyr::mutate(
      point_label =
        dplyr::case_when(
          target_variable == "carbon_content" & y > 0.6 & abs(y - yhat) > 0.1 ~ id_measurement,
          target_variable == "carbon_content" & y < 0.1 & abs(y - yhat) > 0.2 ~ id_measurement,
          target_variable == "potassium_content" & taxon_rank_value == "Juncus effusus" & y > 17000 ~ id_measurement,
          target_variable == "phosphorus_content" & taxon_rank_value == "Juncus effusus" & y > 2000 ~ id_measurement,
          target_variable == "silicon_content" & y > 0.15 & abs(y - yhat) > 0.1 ~ id_measurement,
          target_variable == "calcium_content" & y > 0.04 & abs(y - yhat) > 0.01 ~ id_measurement,
          target_variable == "d13C" & y > -20 & abs(y - yhat) > 5 & taxon_rank_value == "Juncus effusus" ~ id_measurement,
          TRUE ~ NA_integer_
          )
    )
  
  res <- 
    res |>
    dplyr::mutate(
      is_training_data =
        dplyr::case_when(
          is_training_data ~ "Training data",
          TRUE ~ "Test data"
        ) |>
        factor(levels = c("Training data", "Test data")),
      target_variable_name_short = 
        irp_get_names_target_variables(target_variable, what = "target_variable_name_html_short"),
      unit_for_plotting_html_no_element_subscripts = 
        irp_get_units_target_variables(target_variable, what = "unit_for_plotting_html_no_element_subscripts"),
      facet_label =
        paste0(target_variable_name_short, " (", unit_for_plotting_html_no_element_subscripts, ")"),
      facet_label =
        factor(
          facet_label,
          levels =
            tibble::tibble(
              target_variable = unique(target_variable),
              facet_label = unique(facet_label),
              index = irp_get_order_target_variables(target_variable)
            ) |>
            dplyr::arrange(index) |>
            dplyr::pull(facet_label)
        )
    ) |>
    dplyr::arrange(facet_label, is_training_data)
  
  res_plot <- 
    res |>
    ggplot(aes(y = y, x = yhat)) +
    geom_errorbar(aes(ymin = y - y_err, ymax = y + y_err), width = 0, color = "grey")
    
  res_plot <- 
    if(highlight_outliers) {
      res_plot + 
        geom_point(aes(fill = eval(variable_color), shape = is.na(point_label), color = is.na(point_label))) +
        scale_shape_manual(values = c(23, 21), guide = "none") +
        scale_color_manual(values = c("red", "black"), guide = "none")
    } else {
      res_plot + 
        geom_point(aes(fill = eval(variable_color)), shape = 21)
    }
    
  res_plot <- 
    res_plot +
    facet_wrap(~ facet_label, scales = "free", ncol = 4L) +
    geom_abline(intercept = 0, slope = 1, color = "grey50") +
    #geom_text(
    #  data =
    #    res |>
    #    dplyr::filter(! is.na(point_label)),
    #  aes(label = point_label)
    #) +
    labs(
      y = "Measured values",
      x = "Predicted values"
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
  
  if(as.character(variable_color) == "is_training_data") {
    res_plot <- 
      res_plot + 
      # ggnewscale::new_scale_fill() +
      scale_fill_manual(values = c("white", "black")) +
      guides(
        fill = guide_legend(title = variable_color_legend_title, override.aes = list(shape = 21, size = 3))
      )
  } else if(is.character(res[[as.character(variable_color)]])) {
    res_plot <- 
      res_plot +
      guides(
        fill = guide_legend(title = variable_color_legend_title, override.aes = list(size = 3))
      )
  } else {
    res_plot <- 
      res_plot +
      guides(
        fill = guide_colorbar(title = variable_color_legend_title)
      )
  }
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 9, height = 12.5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#' Comparison of model RMSE to other studies
#' 
#' @param color_by
#' 
#' @export
irp_make_plot_15 <- function(irp_data_evaluation_literature_values, irp_fit_1_map_evaluation_1, irp_d_model_info_enriched_1, irp_fit_1_map_elpd_compare, file_plot = "figures/irp_plot_15.pdf") {
  
  id_model_best <- 
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
    dplyr::arrange(id_model, dplyr::desc(elpd_diff)) |>
    dplyr::filter(! duplicated(target_variable)) |>
    dplyr::pull(id_model)
  
  res1 <- 
    irp_d_model_info_enriched_1 |>
    dplyr::select(id_model, target_variable) |>
    dplyr::mutate(
      validation =
        irp_fit_1_map_evaluation_1 |>
        purrr::map(function(.x) {
          .x$rmse |>
            dplyr::filter(! is_training_data) |>
            dplyr::mutate(
              value_type = "RMSE",
              value_mean = rmse_mean,
              value_lower = rmse_lower,
              value_upper = rmse_upper
            ) |>
            dplyr::select(dplyr::starts_with("value_")) |>
            dplyr::bind_rows(
              tibble::tibble(
                value_type = "y_range",
                value_lower = min(.x$yhat$y),
                value_upper = max(.x$yhat$y)
              )
            )
        })
    ) |>
    dplyr::filter(id_model %in% id_model_best) |>
    tidyr::unnest("validation") |>
    dplyr::mutate(
      study_type = "This study"
    ) |>
    dplyr::select(-id_model)
  
  res2 <- 
    irp_data_evaluation_literature_values |>
    dplyr::filter(target_variable != "pH" & summary_statistics_label == "RMSE") |>
    dplyr::mutate(
      dplyr::across(
        dplyr::where(is.numeric),
        as.numeric
      )
    ) |>
    tidyr::nest(RMSE = c("summary_statistics_lower", "summary_statistics_upper"), y_range = c("y_lower", "y_upper")) |>
    dplyr::rename(
      study_type = "reference_pretty"
    ) |>
    dplyr::select(- summary_statistics_label) |>
    dplyr::mutate(
      RMSE = 
        purrr::map(RMSE, function(.x) {
          .x |>
            setNames(nm = c("value_lower", "value_upper"))
        }),
      y_range =
        purrr::map(y_range, function(.x) {
          .x |>
            setNames(nm = c("value_lower", "value_upper"))
        })
    ) |>
    tidyr::pivot_longer(
      dplyr::all_of(c("RMSE", "y_range")),
      names_to = "value_type",
      values_to = "value"
    ) |>
    tidyr::unnest("value") |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(c("value_lower", "value_upper")),
        as.numeric
      )
    ) |>
    dplyr::select(target_variable, study_type, value_type, value_lower, value_upper)
   
  res <-  
    dplyr::bind_rows(
      res1 |>
        dplyr::filter(target_variable %in% res2$target_variable), 
      res2
    ) |>
    dplyr::mutate(
      target_variable_name_short = irp_get_names_target_variables(target_variable, what = "target_variable_name_html_short"),
      unit_for_plotting_html_no_element_subscripts = irp_get_units_target_variables(target_variable, what = "unit_for_plotting_html_no_element_subscripts"),
      facet_label =
        paste0(target_variable_name_short, " (", unit_for_plotting_html_no_element_subscripts, ")"),
      facet_label =
        factor(
          facet_label,
          levels =
            tibble::tibble(
              target_variable = unique(target_variable),
              facet_label = unique(facet_label),
              index = irp_get_order_target_variables(target_variable)
            ) |>
            dplyr::arrange(index) |>
            dplyr::pull(facet_label)
        )
    )
  
  
  p1 <- 
    res |>
    dplyr::filter(value_type == "RMSE") |>
    ggplot(aes(x = study_type)) +
    geom_errorbar(aes(ymin = value_lower, ymax = value_upper), width = 0.4, color = "grey10") +
    geom_point(aes(y = value_mean)) +
    coord_flip() +
    facet_grid(value_type ~ facet_label, scales = "free") +
    labs(
      x = "Study",
      y = ""
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
      legend.title = ggtext::element_markdown(size = 13),
      panel.spacing = unit(1, "lines")
    )
  
  p2 <- 
    res |>
    dplyr::filter(value_type == "y_range") |>
    dplyr::mutate(
      value_type = "Value range"
    ) |>
    ggplot(aes(x = study_type)) +
    geom_errorbar(aes(ymin = value_lower, ymax = value_upper), width = 0.4, color = "grey10") +
    geom_point(aes(y = value_mean)) +
    coord_flip() +
    facet_grid(value_type ~ facet_label, scales = "free") +
    labs(
      x = "Study",
      y = "Value"
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      axis.title.x = ggtext::element_markdown(size = 14),
      axis.title.y = ggtext::element_markdown(size = 14),
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = element_blank(),
      strip.text.y = ggtext::element_markdown(size = 13),
      legend.text = ggtext::element_markdown(size = 13),
      legend.title = ggtext::element_markdown(size = 13),
      panel.spacing = unit(1, "lines")
    )
  
  res_plot <- 
    list(p1, p2) |>
    patchwork::wrap_plots(nrow = 2L, byrow = TRUE) +
    patchwork::plot_annotation(
      tag_levels = c('a'),
      tag_prefix = '(',
      tag_suffix = ')'
    ) +
    patchwork::plot_layout(guides = "collect", heights = c(1, 0.9)) &
    theme(legend.position = "bottom")
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 12, height = 6.5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#' Plots residuals versus a specified target variable
#' 
#' @param highlight_outliers Logical value. If `TRUE`, selected outliers are 
#' highlighted with a different point shape.
#' 
#' @export
irp_make_plot_16 <- function(irp_fit_1_map_evaluation_1, irp_d_model_info_enriched_1, irp_fit_1_map_elpd_compare, irp_pmird_mirs, x_variable = sym("N"), x_variable_axis_title = "N content g g<sup>-1</sup>", highlight_outliers = TRUE, file_plot = "figures/irp_plot_16.pdf") {
  
  id_model_best <- 
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
    dplyr::arrange(id_model, dplyr::desc(elpd_diff)) |>
    dplyr::filter(! duplicated(target_variable)) |>
    dplyr::pull(id_model)
  
  cur_x_variable <- 
    if(x_variable == "is_training_data") {
      sym("id_measurement")
    } else {
      x_variable
    }
  
  res <- 
    irp_d_model_info_enriched_1 |>
    dplyr::select(id_model, target_variable) |>
    dplyr::mutate(
      validation =
        irp_fit_1_map_evaluation_1 |>
        purrr::map(function(.x) {
          .x$yhat |>
            dplyr::mutate( #---note: needed to avoid vctrs error
              dplyr::across(
                dplyr::where(is.integer),
                as.integer
              )
            ) |>
            dplyr::select(-id_model)
        })
    ) |>
    dplyr::filter(id_model %in% id_model_best) |>
    tidyr::unnest("validation") |>
    dplyr::left_join(
      irp_pmird_mirs |>
        dplyr::select(id_measurement, dplyr::all_of(cur_x_variable), taxon_rank_value),
      by = "id_measurement"
    ) |>
    dplyr::mutate(
      point_label =
        dplyr::case_when(
          target_variable == "carbon_content" & y > 0.6 & abs(y - yhat) > 0.1 ~ id_measurement,
          target_variable == "carbon_content" & y < 0.1 & abs(y - yhat) > 0.2 ~ id_measurement,
          target_variable == "potassium_content" & taxon_rank_value == "Juncus effusus" & y > 17000 ~ id_measurement,
          target_variable == "phosphorus_content" & taxon_rank_value == "Juncus effusus" & y > 2000 ~ id_measurement,
          target_variable == "silicon_content" & y > 0.15 & abs(y - yhat) > 0.1 ~ id_measurement,
          target_variable == "calcium_content" & y > 0.04 & abs(y - yhat) > 0.02 ~ id_measurement,
          target_variable == "d13C" & y > -20 & abs(y - yhat) > 5 & taxon_rank_value == "Juncus effusus" ~ id_measurement,
          TRUE ~ NA_integer_
        )
    )
  
  res <- 
    res |>
    dplyr::mutate(
      is_training_data =
        dplyr::case_when(
          is_training_data ~ "Training data",
          TRUE ~ "Test data"
        ) |>
        factor(levels = c("Training data", "Test data")),
      target_variable_name_short = 
        irp_get_names_target_variables(target_variable, what = "target_variable_name_html_short"),
      unit_for_plotting_html_no_element_subscripts = 
        irp_get_units_target_variables(target_variable, what = "unit_for_plotting_html_no_element_subscripts"),
      facet_label =
        paste0(target_variable_name_short, " (", unit_for_plotting_html_no_element_subscripts, ")"),
      facet_label =
        factor(
          facet_label,
          levels =
            tibble::tibble(
              target_variable = unique(target_variable),
              facet_label = unique(facet_label),
              index = irp_get_order_target_variables(target_variable)
            ) |>
            dplyr::arrange(index) |>
            dplyr::pull(facet_label)
        )
    ) |>
    dplyr::arrange(facet_label, is_training_data)
  
  res_plot <- 
    res |>
    ggplot(aes(y = y - yhat, x = eval(x_variable))) +
    #geom_point() +
    facet_wrap(~ facet_label, scales = "free", ncol = 4L)
  
  res_plot <- 
    if(highlight_outliers) {
      res_plot + 
        geom_point(aes(fill = is_training_data, shape = is.na(point_label))) +
        scale_fill_manual(values = c("white", "black")) +
        guides(
          fill = guide_legend(title = "", override.aes = list(shape = 21, size = 3))
        ) +
        scale_shape_manual(values = c(23, 21), guide = "none")
    } else {
      res_plot + 
        scale_fill_manual(values = c("white", "black")) +
        guides(
          fill = guide_legend(title = "", override.aes = list(shape = 21, size = 3))
        ) +
        geom_point(aes(fill = is_training_data, shape = 21))
    }
  
  res_plot <-
    res_plot +
    guides(
      fill = guide_legend(title = "", override.aes = list(shape = 21, size = 3))
    ) +
    labs(
      y = "Measured - predicted values",
      x = x_variable_axis_title
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
    ) +
    geom_smooth(method = "loess", color = "grey50") +
    geom_hline(yintercept = 0)
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 9, height = 12.5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#' Plots MIRS from smaples with contrasting Ca content
#' 
#' @export
irp_make_plot_17 <- function(irp_data_model_preprocessed, file_plot = "figures/irp_plot_17.pdf") {
  
  res <- 
    dplyr::bind_rows(
      irp_data_model_preprocessed[[1]] |> 
        dplyr::filter(Ca > 15000 & ! is.na(Ca)),
      irp_data_model_preprocessed[[1]] |> 
        dplyr::filter(Ca < 500) |> 
        dplyr::slice(1:5)
    ) |>
    dplyr::mutate(
      Ca = Ca/1000000
    ) |>
    ir::ir_bc(method = "rubberband", do_impute = TRUE) |>
    ir::ir_normalise(method = "area")
 
  res_plot <- 
    plot(res) + 
    geom_path(aes(color = Ca)) +
    geom_vline(xintercept = c(871, 1415, 1650, 1750), color = "grey50") +
    labs(
      y = "Intensity (-)",
      x = "Wavenumber cm<sup>-1</sup>"
    ) +
    guides(
      color = guide_colorbar(title = "Ca (g g<sup>-1</sup>)")
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
    width = 7, height = 4.5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}



#' Plots outlier spectra
#' 
#' @export
#' Plots residuals versus a specified target variable
#' 
#' @param highlight_outliers Logical value. If `TRUE`, selected outliers are 
#' highlighted with a different point shape.
#' 
#' @export
irp_make_plot_18 <- function(irp_fit_1_map_evaluation_1, irp_d_model_info_enriched_1, irp_fit_1_map_elpd_compare, irp_data_model_preprocessed, file_plot = "figures/irp_plot_18.pdf") {
  
  id_model_best <- 
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
    dplyr::arrange(id_model, dplyr::desc(elpd_diff)) |>
    dplyr::filter(! duplicated(target_variable)) |>
    dplyr::pull(id_model)
  
  res <- 
    irp_d_model_info_enriched_1 |>
    dplyr::select(id_model, target_variable) |>
    dplyr::mutate(
      validation =
        irp_fit_1_map_evaluation_1 |>
        purrr::map(function(.x) {
          .x$yhat |>
            dplyr::mutate( #---note: needed to avoid vctrs error
              dplyr::across(
                dplyr::where(is.integer),
                as.integer
              )
            ) |>
            dplyr::select(-id_model)
        })
    ) |>
    dplyr::filter(id_model %in% id_model_best) |>
    tidyr::unnest("validation") |>
    dplyr::left_join(
      irp_data_model_preprocessed[[1]] |>
        dplyr::select(id_measurement, spectra, id_dataset, core_label, sample_type, taxon_rank_value),
      by = "id_measurement"
    ) |>
    dplyr::mutate(
      point_label =
        dplyr::case_when(
          target_variable == "carbon_content" & y > 0.6 & abs(y - yhat) > 0.1 ~ id_measurement,
          target_variable == "carbon_content" & y < 0.1 & abs(y - yhat) > 0.2 ~ id_measurement,
          target_variable == "potassium_content" & taxon_rank_value == "Juncus effusus" & y > 17000 ~ id_measurement,
          target_variable == "phosphorus_content" & taxon_rank_value == "Juncus effusus" & y > 2000 ~ id_measurement,
          target_variable == "silicon_content" & y > 0.15 & abs(y - yhat) > 0.1 ~ id_measurement,
          target_variable == "calcium_content" & y > 0.04 & abs(y - yhat) > 0.01 ~ id_measurement,
          target_variable == "d13C" & y > -20 & abs(y - yhat) > 5 & taxon_rank_value == "Juncus effusus" ~ id_measurement,
          TRUE ~ NA_integer_
        )
    ) |>
    dplyr::filter(! is.na(point_label)) |>
    dplyr::select(-y) |>
    ir::ir_as_ir() |>
    ir::ir_bc(method = "rubberband", do_impute = TRUE) |>
    ir::ir_normalise(method = "area")
  
  res <- 
    res |>
    dplyr::mutate(
      target_variable_name_short = 
        irp_get_names_target_variables(target_variable, what = "target_variable_name_html_short"),
      unit_for_plotting_html_no_element_subscripts = 
        irp_get_units_target_variables(target_variable, what = "unit_for_plotting_html_no_element_subscripts"),
      facet_label =
        paste0(target_variable_name_short),
      facet_label =
        factor(
          facet_label,
          levels =
            tibble::tibble(
              target_variable = unique(target_variable),
              facet_label = unique(facet_label),
              index = irp_get_order_target_variables(target_variable)
            ) |>
            dplyr::arrange(index) |>
            dplyr::pull(facet_label)
        )
    )
  
  res_plot <- 
    plot(res) +
    facet_wrap(~ facet_label) +
    labs(
      y = "Intensity (-)",
      x = "Wavenumber cm<sup>-1</sup>"
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
    width = 9, height = 5.5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  list(
    file_plot = file_plot,
    res =
      res |>
      dplyr::select(target_variable, id_dataset, core_label, sample_type, taxon_rank_value, id_measurement)
  )
  
}


#' Plots training and testing prediction domains
#' 
#' @export
irp_make_plot_19 <- function(irp_fit_1_map_evaluation_1, irp_d_model_info_enriched_1, irp_fit_1_map_elpd_compare, irp_fit_1_map_pd, irp_id_model_best, scale_prediction_domains_uniformly = FALSE, irp_data_model_preprocessed, file_plot = "figures/irp_plot_19.pdf") {
  
  id_model_best <- irp_id_model_best
  
  res <- 
    irp_d_model_info_enriched_1 |>
    dplyr::select(id_model, target_variable) |>
    dplyr::mutate(
      prediction_domain =
        if(scale_prediction_domains_uniformly) {
          purrr::map(seq_along(target_variable), function(i) {
            res <- 
              irp_data_model_preprocessed[[1]] |>
              ir::ir_scale(
                center = irp_d_model_info_enriched_1$x_center[[1]],
                scale = irp_d_model_info_enriched_1$x_scale[[1]]
              )
            
            dplyr::bind_rows(
              res |>
                dplyr::filter(id_measurement %in% (irp_d_model_info_enriched_1$data_partition[[i]] |> dplyr::filter(for_prospectr_model) |> dplyr::pull(id_measurement))) |>
                irpeat::irp_as_irp_prediction_domain() |>
                dplyr::mutate(
                  prediction_domain_type = "Training prediction domain"
                ),
              res |>
                dplyr::filter(id_measurement %in% (irp_d_model_info_enriched_1$data_partition[[i]] |> dplyr::filter(for_prospectr_test) |> dplyr::pull(id_measurement))) |>
                irpeat::irp_as_irp_prediction_domain() |>
                dplyr::mutate(
                  prediction_domain_type = "Testing prediction domain"
                ),
              res |>
                dplyr::filter( (mir_mode != "atr_ftir" | is.na(mir_mode)) & (! is_baseline_corrected | stringr::str_detect(core_label, "^peatbog")) & sample_type != "dom") |>
                dplyr::filter(id_dataset != 16) |> #---note: seem to be ATR spectra
                irpeat::irp_as_irp_prediction_domain() |>
                dplyr::mutate(
                  prediction_domain_type = "All spectra"
                )
            )
          })
        } else {
          purrr::map(irp_fit_1_map_pd, function(.x) {
            dplyr::bind_rows(
              .x$train |>
                dplyr::mutate(
                  prediction_domain_type = "Training prediction domain"
                ),
              .x$test |>
                dplyr::mutate(
                  prediction_domain_type = "Testing prediction domain"
                ),
              .x$all |>
                dplyr::mutate(
                  prediction_domain_type = "All spectra"
                )
            )
          })
        }
    ) |>
    dplyr::filter(id_model %in% id_model_best) |>
    tidyr::unnest(prediction_domain)
  
  res <- 
    res |>
    dplyr::mutate(
      prediction_domain_type =
        factor(
          prediction_domain_type, 
          levels = rev(c("Training prediction domain", "Testing prediction domain", "All spectra"))
        ), 
      target_variable_name_short = 
        irp_get_names_target_variables(target_variable, what = "target_variable_name_html_short"),
      facet_label = target_variable_name_short,
      facet_label =
        factor(
          facet_label,
          levels =
            tibble::tibble(
              target_variable = unique(target_variable),
              facet_label = unique(facet_label),
              index = irp_get_order_target_variables(target_variable)
            ) |>
            dplyr::arrange(index) |>
            dplyr::pull(facet_label)
        )
    ) |>
    dplyr::arrange(facet_label, prediction_domain_type, x)
  
  res_plot <- 
    res |>
    ggplot(aes(x = x)) +
    geom_ribbon(
      aes(ymin = ymin, ymax = ymax, fill = prediction_domain_type, group = paste0(id_model, "_", prediction_domain_type)), 
      color = NA, alpha = 1.0
    ) +
    scale_fill_manual(values = c("grey30", "steelblue", "darksalmon")) +
    geom_ribbon(
      data = res |> dplyr::filter(prediction_domain_type == "Testing prediction domain"),
      aes(ymin = ymin, ymax = ymax, group = paste0(id_model, "_", prediction_domain_type)), 
      fill = "steelblue", color = NA, alpha = 1.0
    ) +
    labs(
      y = "Scaled intensity (-)",
      x = "Wavenumber (cm<sup>-1</sup>)"
    ) +
    guides(
      fill = guide_legend(title = "")
    ) +
    facet_wrap(~ facet_label, scales = "free_y", ncol = 4L) +
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
    width = 12.5, height = 10, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#' Plots model coefficients
#' 
#' @export
irp_make_plot_20 <- function(irp_fit_1_map_slopes, file_plot = "figures/irp_plot_20.pdf") {
  
  res <- 
    purrr::map(irp_fit_1_map_slopes, function(.x) {
      readRDS_rvars(.x) |>
        dplyr::mutate(
          id_model = 
            .x |>
            stringr::str_remove(pattern = "^targets_rvars/irp_fit_1_") |>
            stringr::str_remove(pattern = "\\.rds$"),
          target_variable =
            id_model |>
            stringr::str_remove(pattern = "_\\d{1}$"),
          x_group = x <= 2330,
          slope = slope / max(c(abs(min(slope)), max(slope))) 
        )
    }) |>
    dplyr::bind_rows()
  
  # needs: irp_data_model_preprocessed
  # res_spectra <- 
  #   irp_data_model_preprocessed[[1]] |>
  #   dplyr::filter(id_measurement %in% c(4561, 1580, 1917)) |>
  #   dplyr::select(id_measurement, spectra) |>
  #   dplyr::arrange(id_measurement) |>
  #   dplyr::mutate(
  #     id_measurement = seq_along(id_measurement)
  #   )
  # 
  # res_offset <- 
  #   res |>
  #   dplyr::group_by(target_variable) |>
  #   dplyr::summarize(
  #     offset = 0.0 #max(max(posterior::quantile2(slope, probs = 0.95)) - 0.3, 0.05)
  #   )
  # 
  # res_spectra <-
  #   res_spectra |>
  #   dplyr::mutate(
  #     res_offset = list(res_offset)
  #   ) |>
  #   tidyr::unnest("spectra") |>
  #   tidyr::unnest("res_offset") |>
  #   dplyr::group_split(target_variable, id_measurement) |>
  #   purrr::map(function(.x) {
  #     .x |>
  #       dplyr::mutate(
  #         y = y - min(y),
  #         y = y / max(y) * 0.05,
  #         y = y  + 0.07 * id_measurement + offset
  #       )
  #   }) |>
  #   dplyr::bind_rows()
  
  res_plot <- 
    res |>
    dplyr::mutate(
      target_variable = 
        irp_get_names_target_variables(
          x = target_variable,
          what = "target_variable_name_html_short"
        )
    ) |>
    ggplot() +
    geom_hline(yintercept = 0, color = "grey50") +
    ggdist::stat_lineribbon(
      aes(ydist = slope, x = x, group = x_group), 
      linewidth = 0.2
    ) +
    scale_fill_brewer() +
    guides(
      fill = guide_legend(title = "Confidence level", override.aes = list(color = NA))
    ) +
    # geom_path(
    #   data = res_spectra,
    #   mapping = aes(y = y, x = x, group = id_measurement)
    # ) +
    labs(
      y = "Coefficient (-)", 
      x = "Wavenumber (cm<sup>-1</sup>)"
    ) +
    facet_wrap(~ target_variable, scales = "free_y", ncol = 3L) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      axis.title.x = ggtext::element_markdown(size = 17),
      axis.title.y = ggtext::element_markdown(size = 17),
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 15),
      strip.text.y = ggtext::element_markdown(size = 15),
      legend.text = ggtext::element_markdown(size = 15),
      legend.title = ggtext::element_markdown(size = 15)
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 13, height = 18, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#' Predicted total porosity versus measurements from Whittington (2024)
#' 
#' @export
irp_make_plot_21 <- function(irp_data_wittington2024, file_plot = "figures/irp_plot_21.pdf") {
  
  res_plot <- 
    irp_data_wittington2024 |>
    dplyr::mutate(
      core_label = stringr::str_remove(core_label, pattern = "^core"),
      ytype =
        dplyr::case_when(
          ytype == "" ~ "Measurements from Whittington (2024)",
          ytype == "rep" ~ "Predictions by our model"
        )
    ) |>
    dplyr::filter(! is.na(total_porosity_1_mean)) |>
    ggplot(aes(y = as.numeric(total_porosity_1_mean), x = sample_depth)) +
    geom_ribbon(
      aes(
        ymin = as.numeric(total_porosity_1_lower), 
        ymax = as.numeric(total_porosity_1_upper),
        group = ytype
      ), 
      color = NA, fill = "grey", alpha = 0.3
    ) +
    geom_path(aes(color = ytype)) +
    scale_color_manual(values = c("black", "coral")) +
    facet_wrap(~ core_label) +
    coord_flip() +
    scale_x_reverse() +
    labs(
      x = "Sample depth (cm)",
      y = "Total porosity (L L<sup>-1</sup>)"
    ) +
    guides(color = guide_legend(title = "", override.aes = list(linewidth = 3))) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      axis.title.x = ggtext::element_markdown(size = 17),
      axis.title.y = ggtext::element_markdown(size = 17),
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 15),
      strip.text.y = ggtext::element_markdown(size = 15),
      legend.text = ggtext::element_markdown(size = 15),
      legend.title = ggtext::element_markdown(size = 15),
      legend.direction = "vertical"
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 6.5, height = 6.5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}

#' Predicted saturated hydraulic conductivity versus measurements from Whittington (2024) and versus model 4 predictions from Morris et al. (2022)
#' 
#' @export
irp_make_plot_22 <- function(irp_data_wittington2024, file_plot = "figures/irp_plot_22.pdf") {
  
  res_plot <- 
    irp_data_wittington2024 |>
    dplyr::filter(ytype != "morris2022_model4") |>
    dplyr::mutate(
      core_label = stringr::str_remove(core_label, pattern = "^core"),
      ytype =
        dplyr::case_when(
          ytype == "" ~ "Measurements from Whittington (2024)",
          ytype == "rep" ~ "Predictions by our model",
          ytype == "morris2022_model1" ~ "Predictions by model 1 from Morris et al. (2022)"
        )
    ) |>
    dplyr::filter(! is.na(saturated_hydraulic_conductivity_1_mean)) |>
    ggplot(aes(y = as.numeric(saturated_hydraulic_conductivity_1_mean), x = sample_depth)) +
    geom_ribbon(
      aes(
        ymin = as.numeric(saturated_hydraulic_conductivity_1_lower), 
        ymax = as.numeric(saturated_hydraulic_conductivity_1_upper),
        group = ytype
      ), 
      color = NA, fill = "grey", alpha = 0.3
    ) +
    geom_path(aes(color = ytype)) +
    scale_color_manual(values = c("black", "lightsteelblue", "coral")) +
    facet_wrap(~ core_label) +
    coord_flip() +
    scale_x_reverse() +
    scale_y_log10() +
    labs(
      x = "Sample depth (cm)",
      y = "Hydraulic conductivity (m s<sup>-1</sup>)"
    ) +
    guides(color = guide_legend(title = "", override.aes = list(linewidth = 3))) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      axis.title.x = ggtext::element_markdown(size = 17),
      axis.title.y = ggtext::element_markdown(size = 17),
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 15),
      strip.text.y = ggtext::element_markdown(size = 15),
      legend.text = ggtext::element_markdown(size = 15),
      legend.title = ggtext::element_markdown(size = 15),
      legend.direction = "vertical"
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 6.5, height = 6.5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#### Tables ####

#' Model evaluation table for the manuscript
#' 
#' @export
irp_make_table_1 <- function(irp_fit_1_map_evaluation_1, irp_d_model_info_enriched_1, irp_fit_1_map_elpd_compare, irp_id_model_best) {
  
  irp_fit_1_map_elpd_compare |>
    dplyr::filter(id_model %in% irp_id_model_best) |>
    dplyr::select(-elpd, -se_elpd) |>
    dplyr::left_join(
      irp_d_model_info_enriched_1 |>
        dplyr::select(id_model, id_preprocessing, n_training_sample, id_measurement_all) |>
        dplyr::mutate(
          n_testing_sample = 
            purrr::map2_int(n_training_sample, id_measurement_all, function(.x, .y) {
              length(.y) - .x
            })
        ) |>
        dplyr::select(-id_measurement_all) |>
        dplyr::bind_cols(
          purrr::map(irp_fit_1_map_evaluation_1, function(.x) {
            .x$rmse |>
              dplyr::mutate(
                target_variable = 
                  .x$yhat$id_model[[1]] |>
                  stringr::str_remove(pattern = "_\\d+$"),
                digits = irp_get_rounding_digits_target_variables(target_variable),
                unit = irp_get_units_target_variables_latex(target_variable, what = "unit_latex_no_element_subscripts"),
                num_divergent = .x$num_divergent,
                dplyr::across(
                  dplyr::starts_with(c("rmse_", "bias_")),
                  function(.x1) round(.x1, digits = digits[[1]])
                ),
                rmse = paste0(rmse_mean, " (", rmse_lower,", ", rmse_upper, ")"),
                bias = paste0(bias_mean, " (", bias_lower,", ", bias_upper, ")"),
                rmse_train_minus_test = paste0(rmse_train_minus_test_mean, " (", rmse_train_minus_test_lower,", ", rmse_train_minus_test_upper, ")"),
                is_training_data = ifelse(is_training_data, "train", "test")
              ) |>
              dplyr::select(target_variable, unit, num_divergent, is_training_data, rmse, bias, rmse_train_minus_test) |>
              tidyr::nest(res = c("rmse", "bias", "rmse_train_minus_test")) |>
              tidyr::pivot_wider(
                values_from = "res",
                names_from = "is_training_data" 
              ) |>
              tidyr::unnest(c("train", "test"), names_sep = "_")
          }) |>
            dplyr::bind_rows(),
          # mcse
          purrr::map(irp_fit_1_map_evaluation_1, function(.x) {
            .x$yhat |>
              dplyr::summarise(
                dplyr::across(
                  dplyr::starts_with("y_mcse"),
                  max
                ),
                .groups = "drop"
              ) |>
              dplyr::mutate(
                target_variable = 
                  .x$yhat$id_model[[1]] |>
                  stringr::str_remove(pattern = "_\\d+$"),
                digits = irp_get_rounding_digits_target_variables(target_variable),
                y_range = 
                  .x$yhat$y |>
                  range() |>
                  round(digits = digits) |>
                  paste0(collapse = " to "),
                dplyr::across(
                  dplyr::starts_with(c("y_mcse")),
                  function(.x1) round(.x1, digits = 4L)
                )
              )
          }) |>
            dplyr::bind_rows() |>
            dplyr::select(-target_variable, -digits)
        ),
      by = c("id_model", "target_variable")
    ) |>
    dplyr::mutate(
      index_target_variable = irp_get_order_target_variables(target_variable),
      target_variable_pretty = irp_get_names_target_variables(target_variable, what = "target_variable_name_latex_short"),
      dplyr::across(
        dplyr::all_of(c("elpd_diff", "se_diff")),
        function(.x1) round(.x1, digits = 1L)
      )
    ) |>
    dplyr::arrange(index_target_variable, dplyr::desc(elpd_diff)) |>
    dplyr::select(-index_target_variable) |>
    dplyr::select(! dplyr::starts_with("test_rmse_train_minus_")) |>
    dplyr::relocate(elpd_diff, se_diff, .after = "n_testing_sample") |>
    dplyr::relocate(target_variable_pretty, .after = "target_variable")
  
}


#' Summary of gap filling results
#' 
#' @export
irp_make_table_2 <- function(irp_pmird_mirs, irp_pmird_gap_filled) {
  
  res1 <- 
    irp_pmird_mirs |>
    dplyr::mutate(
      to_keep = (mir_mode == "absorbance_ftir" | is.na(mir_mode)) & sample_type %in% c("peat", "vegetation", "litter")
    )
  
  res2 <- res1 |> dplyr::filter(to_keep)
  
  # summarize
  res <- 
    tibble::tibble(
      target_variable_irpeat = colnames(irp_pmird_gap_filled$yhat_auxilliary),
      target_variable =
        target_variable_irpeat |>
        stringr::str_remove(pattern = "_\\d{1}$"),
      target_variable_order = 
        irp_get_order_target_variables(target_variable),
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
          target_variable == "H_to_C" ~ "H_to_C",
          target_variable == "volume_fraction_solids" ~ NA_character_,
          target_variable == "non_macroporosity" ~ NA_character_,
          target_variable == "macroporosity" ~ "macroporosity",
          target_variable == "saturated_hydraulic_conductivity" ~ "hydraulic_conductivity",
          target_variable == "specific_heat_capacity" ~ NA_character_,
          target_variable == "dry_thermal_conductivity" ~ NA_character_
        ),
      target_variable_pretty = 
        irp_get_names_target_variables(
          target_variable, 
          what = "target_variable_name_latex_short"
        ),
      n_total = nrow(res2),
      index_measured =
        purrr::map(target_variable_label, function(.x) {
          res <- rep(FALSE, nrow(res2))
          #if(! is.na(.x)) {
          #  index <- 
          #    switch(
          #    .x,
          #    "nosc" =,
          #    "dgf0" = ! is.na(res2$C) & ! is.na(res2$H) & ! is.na(res2$O) & ! is.na(res2$N), 
          #    "C_to_N" = ! is.na(res2$C) & ! is.na(res2$N),
          #    "O_to_C" = ! is.na(res2$O) & ! is.na(res2$C),
          #    "H_to_C" = ! is.na(res2$H) & ! is.na(res2$C),
          #    ! is.na(res2[[.x]])
          #  )
          #  res[index] <- TRUE
          #}
          if(! is.na(.x)) {
            res[! is.na(res2[[.x]])] <- TRUE
          }
          res
        }),
      index_auxilliary =
        purrr::map(target_variable_irpeat, function(.x) {
          ! is.na(irp_pmird_gap_filled$yhat_auxilliary[[.x]])
        }),
      n_measured =
        purrr::map_int(index_measured, sum),
      n_fillable = n_total - n_measured,
      n_filled_auxilliary =
        purrr::map_int(seq_along(target_variable_irpeat), function(i) {
          sum(index_auxilliary[[i]] & ! index_measured[[i]])
        }),
      n_filled_train =
        purrr::map_int(seq_along(target_variable_irpeat), function(i) {
          sum(irp_pmird_gap_filled$is_in_training_pd[[target_variable_irpeat[[i]]]] & ! index_measured[[i]] & ! index_auxilliary[[i]])
        }),
      n_filled_test =
        purrr::map_int(seq_along(target_variable_irpeat), function(i) {
          sum(irp_pmird_gap_filled$is_in_testing_pd[[target_variable_irpeat[[i]]]] & ! index_measured[[i]] & ! index_auxilliary[[i]])
        }),
      n_filled = 
        purrr::map_int(seq_along(target_variable_irpeat), function(i) {
          sum((irp_pmird_gap_filled$is_in_training_pd[[target_variable_irpeat[[i]]]] | irp_pmird_gap_filled$is_in_testing_pd[[target_variable_irpeat[[i]]]] | index_auxilliary[[i]]) & ! index_measured[[i]])
        }),
      f_fillable_filled = n_filled/n_fillable,
      n_fillable_not_filled =
        purrr::map_int(target_variable_irpeat, function(.x) {
          sum(! (irp_pmird_gap_filled$is_in_training_pd[[.x]] | irp_pmird_gap_filled$is_in_testing_pd[[.x]]))
        })
    ) |>
    dplyr::select(-index_measured, -index_auxilliary) |>
    dplyr::mutate(
      target_variable_pretty =
        dplyr::case_when(
          ! is.na(target_variable_pretty) ~ target_variable_pretty,
          target_variable == "volume_fraction_solids" ~ "Volume fraction of solids",
          target_variable == "non_macroporosity" ~ "Non-macroporosity",
          target_variable == "macroporosity" ~ "Macroporosity",
          target_variable == "saturated_hydraulic_conductivity" ~ "Saturated hydraulic conductivity",
          target_variable == "specific_heat_capacity" ~ "Specific heat capacity",
          target_variable == "dry_thermal_conductivity" ~ "Dry thermal conductivity"
        )
    ) |>
    dplyr::arrange(target_variable_order) |>
    dplyr::select(-target_variable_order)
  
  attr(res, "n_total") <- nrow(res2)
  
  res
  
}



#' Table of model coefficients
#' 
#' @export
#' Plots training and testing prediction domains
#' 
#' @export
irp_make_table_3 <- function(irp_fit_1_map_slopes) {
  
    purrr::map(irp_fit_1_map_slopes, function(.x) {
      readRDS_rvars(.x) |>
        dplyr::mutate(
          id_model = 
            .x |>
            stringr::str_remove(pattern = "^targets_rvars/irp_fit_1_") |>
            stringr::str_remove(pattern = "\\.rds$"),
          target_variable =
            id_model |>
            stringr::str_remove(pattern = "_\\d{1}$"),
          pr_larger_than_0 = posterior::Pr(slope > 0),
          pr_smaller_than_0 = posterior::Pr(slope < 0)
        )
    }) |>
    dplyr::bind_rows() |>
    dplyr::filter(pr_larger_than_0 >= 0.9 | pr_smaller_than_0 >= 0.9) |>
    dplyr::mutate(
      slope = paste0(round(median(slope), 2), " (", posterior::quantile2(slope, probs = c(0.05)) |> round(2), ", ", posterior::quantile2(slope, probs = c(0.95)) |> round(2), ")")
    )
  
}





