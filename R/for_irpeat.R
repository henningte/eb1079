#' Creates a list with preprocessing configurations that is compatible with irpeat
#' 
#' @export
irp_make_mir_model_config <- function(irp_mir_preprocessing_settings, irp_d_model_info, prediction_domain) {
    
    res <- 
      irp_mir_preprocessing_settings |>
      dplyr::slice(irp_d_model_info$id_preprocessing[[1]]) |>
      dplyr::select(-id_preprocessing) |>
      dplyr::mutate(
        do_scale = TRUE,
        do_bc = TRUE
      ) |>
      unclass()
    
    res$scale_center <- irp_d_model_info$x_center[[1]]
    res$scale_scale <- irp_d_model_info$x_scale[[1]]
    
    list(
      target_variable = irp_d_model_info$id_model[[1]],
      likelihood = irp_d_model_info$likelihood[[1]],
      unit = 
        irp_d_model_info$target_variable |> 
        irp_get_units_for_irpeat(),
      model_scale =
        list(
          y_center = irp_d_model_info$y_center[[1]],
          y_scale = irp_d_model_info$y_scale[[1]]
        ),
      irp_preprocess = res,
      prediction_domain = prediction_domain[c(1,2)]
    )
  
}



#' Preprocesses spectra according to a preprocessing config object
#' 
#' @param x An ir object to be preprocessed.
#' 
#' @param config A list with arguments compatible with `irpeat::irp_preprocess()`.
#' 
#' @export
irp_preprocess_eb1079 <- function (x, config) {
  
  res <- x
  
  res |>
    #ir::ir_interpolate_region(range = tibble::tibble(start = c(645, 2230), end = c(695, 2410))) |>
    ir::ir_interpolate_region(range = tibble::tibble(start = c(650), end = c(695))) |>
    irpeat::irp_preprocess(
      do_interpolate = config$irp_preprocess$do_interpolate, 
      interpolate_start = config$irp_preprocess$interpolate_start[[1]], 
      interpolate_dw = config$irp_preprocess$interpolate_dw, 
      do_clip = config$irp_preprocess$do_clip[[1]], 
      clip_range = config$irp_preprocess$clip_range[[1]], 
      do_interpolate_region = FALSE, 
      interpolate_region_range = config$irp_preprocess$interpolate_region_range[[1]], 
      do_bc = config$irp_preprocess$do_bc, 
      bc_method = config$irp_preprocess$bc_method, 
      bc_cutoff = config$irp_preprocess$bc_cutoff, 
      bc_do_impute = config$irp_preprocess$bc_do_impute, 
      do_smooth = config$irp_preprocess$do_smooth, 
      smooth_method = config$irp_preprocess$smooth_method, 
      smooth_p = config$irp_preprocess$smooth_p, 
      smooth_n = config$irp_preprocess$smooth_n, 
      smooth_m = config$irp_preprocess$smooth_m, 
      smooth_ts = config$irp_preprocess$smooth_ts, 
      smooth_k = config$irp_preprocess$smooth_k, 
      do_normalise = config$irp_preprocess$do_normalise, 
      normalise_method = config$irp_preprocess$normalise_method, 
      do_bin = config$irp_preprocess$do_bin, 
      bin_width = config$irp_preprocess$bin_width, 
      bin_new_x_type = config$irp_preprocess$bin_new_x_type, 
      do_scale = FALSE, #config$irp_preprocess$do_scale, 
      scale_center = config$irp_preprocess$scale_center, 
      scale_scale = config$irp_preprocess$scale_scale, 
      do_return_as_ir = TRUE
    ) |> 
    ir::ir_clip(range = data.frame(start = c(650, 2400), end = c(2250, 4000))) |>
    ir::ir_scale(center = config$irp_preprocess$scale_center, scale = config$irp_preprocess$scale_scale)
  
}


#' Helper that prepares the data for the brms models
#' 
#'@export
irp_predict_for_eb1079_helper_1 <- function(x, config) {
  
  newdata <- 
    tibble::tibble(
      x = 
        x |>
        ir::ir_flatten() |>
        dplyr::select(-1) |>
        t() |>
        tibble::as_tibble() |>
        setNames(nm = paste0("V", x$spectra[[1]]$x)) |>
        as.matrix()
    )
  
  newdata
  
}


#' Makes predictions with a model
#' 
#' @param x An ir object
#' 
#' @param config A list with arguments compatible with `irpeat::irp_preprocess()`.
#' 
#' @export
irp_predict_for_eb1079 <- function(x, model, config) {
  
  x_or <- x
  x <- irp_preprocess_eb1079(x = x, config = config)
  x_in_pd <- irpeat::irp_is_in_prediction_domain(x = x, prediction_domain = config$prediction_domain$train)
  
  newdata <- 
    irp_predict_for_eb1079_helper_1(
      x = x, 
      config = config
    )
  
  res <- 
    tidybayes::predicted_rvars(
      newdata = newdata,
      object = model,
      value = "yhat"
    )
  
  x_or[[config$target_variable]] <- res$yhat
  x_or[[paste0(config$target_variable, "_in_pd")]] <- x_in_pd$is_in_prediction_domain
  
  x_or
  
}


