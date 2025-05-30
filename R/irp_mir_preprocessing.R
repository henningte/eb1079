#' Defines settings for the spectral preprocessing for this project
#' 
#' @export
irp_get_mir_preprocessing_settings <- function() {
  
  spectral_preprocessing_n <- 2L
  
  d_spectral_preprocessing_settings <- 
    tibble::tibble(
      do_interpolate = TRUE,
      interpolate_start = rep(list(NULL), spectral_preprocessing_n),
      interpolate_dw = 1,
      do_clip = TRUE,
      clip_range = list(tibble::tibble(start = 650, end = 4000)),
      do_interpolate_region = FALSE,
      interpolate_region_range = rep(list(NULL), spectral_preprocessing_n),
      do_bc = TRUE,
      bc_method = "rubberband",
      bc_cutoff = 5, 
      bc_do_impute = TRUE,
      do_smooth = c(FALSE),
      smooth_method = "sg",
      smooth_p = 3,
      smooth_n = 31,
      smooth_m = c(0, 1),
      smooth_ts = 1,
      smooth_k = NA_real_,
      do_normalise = TRUE,
      normalise_method = "snv",
      do_bin = TRUE,
      bin_width = 10,
      bin_new_x_type = "mean",
      do_scale = FALSE, # ---note: don't scale now because this can only be done with the final data subset
      scale_center = NA,
      scale_scale = NA,
      id_preprocessing = seq_along(scale_scale)
    )  
  
  d_spectral_preprocessing_settings
  
}



#' Preprocess MIRS for model development
#' 
#' @export
irp_make_data_model_preprocessed <- function(irp_pmird_mirs, irp_mir_preprocessing_settings) {
  
  res <- irp_pmird_mirs
  
  # load the reference spectra for atmospheric correction
  d_bg_co2 <- pmird::mir_quality_config$mir_co2_reference_spectrum_raw
  d_bg_wv <- pmird::mir_quality_config$mir_water_vapor_reference_spectrum_raw
  
  # First: get a baseline which can be subtracted from all spectra even with negative CO2 peaks. I have to do a SG smoothing and regional interpolation here to avoid negative CO2 peaks elsewhere to corrupt the baseline. I also have to interpolate linearly the CO2 peak around 670 cm$^{-1}$ since this peak corrupts the baseline and does not get smoothed out completely by the SG smoothing. 
  spectral_preprocessing_clip_range <- irp_mir_preprocessing_settings$clip_range[[1]]
  res_bl <- 
    res |>
    dplyr::filter(x_variable_min <= spectral_preprocessing_clip_range$start & x_variable_max >= spectral_preprocessing_clip_range$end) |>
    ir::ir_interpolate(start = NULL, dw = 1) |>
    ir::ir_smooth(method = "sg", n = 91) |>
    ir::ir_clip(range = spectral_preprocessing_clip_range) |>
    ir::ir_interpolate_region(range = tibble::tibble(start = c(645, 2230), end = c(695, 2410))) |>
    ir::ir_bc(method = "rubberband", do_impute = TRUE, return_bl = TRUE)
  
  ## Second: perform the actual correction
  
  # clipping and baseline correction
  res_preprocessed <-
    res |>
    dplyr::filter(x_variable_min <= spectral_preprocessing_clip_range$start & x_variable_max >= spectral_preprocessing_clip_range$end) |>
    ir::ir_interpolate(start = NULL, dw = 1) |>
    ir::ir_clip(range = spectral_preprocessing_clip_range) |>
    ir::ir_subtract(res_bl)
  
  # correct co2 artifacts
  res_preprocessed <- 
    res_preprocessed |>
    ir::ir_correct_atmosphere(
      ref = 
        d_bg_co2 |> 
        dplyr::slice(match(res_preprocessed$measurement_device, d_bg_co2$measurement_device_target)) |>
        ir::ir_as_ir(),
      wn1 = 2361, wn2 = 2349, 
      do_interpolate = FALSE
    ) |>
    ir::ir_interpolate_region(range = tibble::tibble(start = c(650), end = c(695))) #---note: This is necessary because for some spectra, the CO2 peak around 2300 cm$^{-1}$ was removed manually, but not the peak at ~670 cm$^{-1}$. This would otherwise cause an uplift.
  
  # correct water vapor artifacts
  res_preprocessed <- 
    res_preprocessed |>
    ir::ir_correct_atmosphere(
      ref = 
        d_bg_wv |> 
        dplyr::slice(match(res_preprocessed$measurement_device, d_bg_wv$measurement_device_target)) |>
        ir::ir_as_ir(),
      wn1 = 3902, wn2 = 3912, 
      do_interpolate = FALSE
    )
  
  # do a final baseline correction to avoid negative values
  res_preprocessed <- 
    res_preprocessed |>
    ir::ir_bc(method = "rubberband", do_impute = TRUE)
  
  d_spectral_preprocessing_settings <- irp_mir_preprocessing_settings
  
  d_spectral_preprocessing_settings |>
    dplyr::mutate(
      spectra =
        purrr::map(id_preprocessing, function(i) {
          
          # do the preprocessing with `irpeat::ir_preprocess()`
          res <- 
            res_preprocessed |>
            irpeat::irp_preprocess(
              do_interpolate = TRUE,
              interpolate_start = interpolate_start[[i]][[1]],
              interpolate_dw = interpolate_dw[[i]],
              do_clip = TRUE,
              clip_range = clip_range[[i]],
              do_interpolate_region = FALSE,
              interpolate_region_range = interpolate_region_range[[i]][[1]],
              do_bc = do_bc[[i]],
              bc_method = bc_method[[i]],
              bc_cutoff = bc_cutoff[[i]],
              bc_do_impute = bc_do_impute[[i]], 
              do_smooth = do_smooth[[i]],
              smooth_method = smooth_method[[i]],
              smooth_p = smooth_p[[i]],
              smooth_n = smooth_n[[i]],
              smooth_m = smooth_m[[i]],
              smooth_ts = smooth_ts[[i]],
              smooth_k = smooth_k[[i]],
              do_normalise = do_normalise[[i]],
              normalise_method = normalise_method[[i]],
              do_bin = do_bin[[i]],
              bin_width = bin_width[[i]],
              bin_new_x_type = bin_new_x_type[[i]],
              do_scale = do_scale[[i]],
              scale_center = scale_center[[i]],
              scale_scale = scale_scale[[i]],
              do_return_as_ir = TRUE
            ) |>
            ir::ir_clip(range = data.frame(start = c(650, 2400), end = c(2250, 4000)))
          
        })
    ) |>
    dplyr::pull(spectra)
  
}

















