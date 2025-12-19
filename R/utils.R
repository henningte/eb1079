#' Creates a connection to the pmird database
#' 
#' @export
irp_get_pmird_connection <- function() {
  
  RMariaDB::dbConnect(
    drv = RMariaDB::MariaDB(),
    dbname = "pmird",
    username = "root", # ---todo: adjust or get from config file
    password = "coucou",
    host = "mariadb"
  )
  
}

#' Pipe-friendly `do.call()`
#' 
#' @param f Character value. Name of function to apply
#' 
#' @export
irp_do_call <- function(x, f) {
  do.call(f, x)
}


#' Saving and loading tibbles with rvars columns
#'
#' See https://github.com/stan-dev/posterior/issues/307
#'
#' @param object Object to save.
#'
#' @param file Character value. Path to the file.
#'
#' @export
saveRDS_rvars <- function(object, file) {
  
  saveRDS(
    object = object,
    file = file,
    refhook = \(x) if (any(c("vec_proxy", "vec_proxy_equal") %in% names(x))) ""
  )
  
}

#' @export
readRDS_rvars <- function(file) {
  
  readRDS(
    file = file,
    refhook = \(x) new.env()
  )
  
}



#' Defines settings for MCMC sampling
#' 
#' @export
irp_get_mcmc_settings <- function() {
  
  list(
    iter = 5000L,
    warmup = 3000L,
    chains = 4L,
    adapt_delta = 0.9,
    max_treedepth = 12
  )
  
}

irp_get_mcmc_settings_horseshoe <- function() {
  
  list(
    iter = 5000L,
    warmup = 3000L,
    chains = 4L,
    adapt_delta = 0.99,
    max_treedepth = 14
  )
  
}

#' Sets units of the array of an rvar object
#' 
#' @param ... Additional arguments passed to `units::set_units()`.
#' 
#' @export
irp_set_units_rvar <- function(x, ...) {
  posterior::draws_of(x) <- units::set_units(posterior::draws_of(x), ...)
  x
}

#' Drops units of the array of an rvar object
#' 
#' @export
irp_drop_units_rvar <- function(x) {
  posterior::draws_of(x) <- units::drop_units(posterior::draws_of(x))
  x
}


#' Computes the RMSE with rvar objects
#' @export
irp_rmse_rvar <- function(yhat, y) {
  
  sqrt(posterior::rvar_mean((yhat - y)^2)) 
  
}

#' Computes the RMSE with rvar objects
#' @export
irp_rmse <- function(yhat, y) {
  
  sqrt(mean((yhat - y)^2)) 
  
}

#### Target variables ####

#' Produces a data frame with general information on the target variables
#' 
#' @export
irp_get_tv_table <- function(irp_target_variables) {
  
  tibble::tibble(
    variable = irp_target_variables,
    definition =
      dplyr::case_when(
        variable == "carbon_content" ~ "Mass content of C in 1 g bulk peat.",
        variable == "nitrogen_content" ~ "Mass content of N in 1 g bulk peat.",
        variable == "oxygen_content" ~ "Mass content of O in 1 g bulk peat.",
        variable == "hydrogen_content" ~ "Mass content of H in 1 g bulk peat.",
        variable == "phosphorus_content" ~"Mass content of P in 1 g bulk peat.",
        variable == "sulfur_content" ~ "Mass content of S in 1 g bulk peat.",
        variable == "potassium_content" ~ "Mass content of K in 1 g bulk peat.",
        variable == "titanium_content" ~ "Mass content of Ti in 1 g bulk peat.",
        variable == "silicon_content" ~ "Mass content of Si in 1 g bulk peat.",
        variable == "calcium_content" ~ "Mass content of Ca in 1 g bulk peat.",
        variable == "d13C" ~ "$\\delta^{13}$C value of bulk peat relative to the Vienna Pee Dee Bee standard.",
        variable == "d15N" ~ "$\\delta^{15}$N value of bulk peat relative to the Air N$_2$ standard.",
        variable == "nosc" ~ "Nominal oxidation state of carbon as defined in (ref:Masiello2008-textual)",
        variable == "dgf0" ~ "Standard free Gibbs energy of formation (25°C, 1 bar).",
        variable == "loss_on_ignition" ~ "Fraction of initial mass lost during combustion of the dried sample at 400°C.",
        variable == "bulk_density" ~ "Mass of the dried sample divided by its volume.",
        variable == "C_to_N" ~ "The mass ratio of a samples' C and N content.",
        variable == "O_to_C" ~ "The mass ratio of a samples' O and C content.",
        variable == "H_to_C" ~ "The mass ratio of a samples' H and C content."
      ),
    significance =
      dplyr::case_when(
        variable == "carbon_content" ~ "Estimating peat C stocks.",
        variable == "hydrogen_content" ~ "Estimating peat H stocks.",
        variable == "nitrogen_content" ~ "Estimating peat N stocks. Quantifying nutrient limitation. Quantifying peat decomposability.",
        variable == "oxygen_content" ~ "Estimating peat O stocks.",
        variable == "sulfur_content" ~ "Estimating peat S stocks.",
        variable == "phosphorus_content" ~ "Estimating peat P stocks. Quantifying nutrient limitations.",
        variable == "potassium_content" ~ "Estimating peat K stocks. Quantifying nutrient limitations.",
        variable == "calcium_content" ~ "Estimating peat Ca stocks. Quantifying minerotrophy.",
        variable == "titanium_content" ~ "Estimating peat Ti stocks. Quantifying mineral dust inputs and decomposition mass losses.",
        variable == "silicon_content" ~ "Estimating peat mineral inputs. Si controls the Fe and P cycle (ref:)",
        variable == "d13C" ~ "Estimating degree of decomposition, moisture conditions during photosynthesis.",
        variable == "d15N" ~ "Estimating degree of decomposition, moisture conditions during photosynthesis.",
        variable == "nosc" ~ "Estimating degree of decomposition. Computation of the oxidative ratio (ref:Masiello2008).",
        variable == "dgf0" ~ "Estimating thermodynamic feasibility of reactions.",
        variable == "loss_on_ignition" ~ "Estimating organic matter pools and mineral pools. Quantifying degree of decomposition.",
        variable == "bulk_density" ~ "Estimating peat hydraulic properties (ref:Liu2019). Quantifying peat element stocks. Quantifying peat degree of decomposition.",
        variable == "C_to_N" ~ "Estimating the degree of decomposition, estimating the abundance of organic matter fractions from Van Krevelen diagrams.",
        variable == "O_to_C" ~ "Estimating the degree of decomposition, estimating the abundance of organic matter fractions from Van Krevelen diagrams.",
        variable == "H_to_C" ~ "Estimating the degree of decomposition, estimating the abundance of organic matter fractions from Van Krevelen diagrams."
      )
  )
  
}


#' Defines units for the target variables
#' 
#' @param what Character value. One of `unit_for_modeling`, `unit_for_plotting`,
#' `unit_for_plotting_html`, `unit_for_plotting_html_no_element_subscripts`,
#' `unit_latex_no_element_subscripts`.
#' 
#' @export
irp_get_units_target_variables <- function(x, what) {
  
  res <- 
    tibble::tibble(
      target_variable = x
    ) |>
    dplyr::mutate(
      unit_for_modeling =
        dplyr::case_when(
          target_variable == "carbon_content" ~ "g_C/g_sample",
          target_variable == "nitrogen_content" ~ "g_N/g_sample",
          target_variable == "oxygen_content" ~ "g_O/g_sample",
          target_variable == "hydrogen_content" ~ "g_H/g_sample",
          target_variable == "phosphorus_content" ~ "g_P/g_sample",
          target_variable == "sulfur_content" ~ "g_S/g_sample",
          target_variable == "potassium_content" ~ "g_K/g_sample",
          target_variable == "titanium_content" ~ "g_Ti/g_sample",
          target_variable == "silicon_content" ~ "g_Si/g_sample",
          target_variable == "calcium_content" ~ "g_Ca/g_sample",
          target_variable == "d13C" ~ "1",
          target_variable == "d15N" ~ "1",
          target_variable == "nosc" ~ "1",
          target_variable == "dgf0" ~ "J/mol_C",
          target_variable == "loss_on_ignition" ~ "g/g_sample",
          target_variable == "bulk_density" ~ "g/cm^3",
          target_variable == "C_to_N" ~ "g_C/g_N",
          target_variable == "O_to_C" ~ "g_O/g_C",
          target_variable == "H_to_C" ~ "g_H/g_C" 
        ),
      unit_for_plotting =
        dplyr::case_when(
          target_variable == "carbon_content" ~ "g_C/g_sample",
          target_variable == "nitrogen_content" ~ "g_N/g_sample",
          target_variable == "oxygen_content" ~ "g_O/g_sample",
          target_variable == "hydrogen_content" ~ "g_H/g_sample",
          target_variable == "phosphorus_content" ~ "ug_P/g_sample",
          target_variable == "sulfur_content" ~ "ug_S/g_sample",
          target_variable == "potassium_content" ~ "ug_K/g_sample",
          target_variable == "titanium_content" ~ "ug_Ti/g_sample",
          target_variable == "silicon_content" ~ "g_Si/g_sample",
          target_variable == "calcium_content" ~ "g_Ca/g_sample",
          target_variable == "d13C" ~ "1",
          target_variable == "d15N" ~ "1",
          target_variable == "nosc" ~ "1",
          target_variable == "dgf0" ~ "kJ/mol_C",
          target_variable == "loss_on_ignition" ~ "g/g_sample",
          target_variable == "bulk_density" ~ "g/cm^3",
          target_variable == "C_to_N" ~ "g_C/g_N",
          target_variable == "O_to_C" ~ "g_O/g_C",
          target_variable == "H_to_C" ~ "g_H/g_C" 
        ),
      unit_for_plotting_html =
        dplyr::case_when(
          target_variable == "carbon_content" ~ "g<sub>C</sub> g<sup>-1</sup><sub>sample</sub>",
          target_variable == "nitrogen_content" ~ "g<sub>N</sub> g<sup>-1</sup><sub>sample</sub>",
          target_variable == "oxygen_content" ~ "g<sub>O</sub> g<sup>-1</sup><sub>sample</sub>",
          target_variable == "hydrogen_content" ~ "g<sub>H</sub> g<sup>-1</sup><sub>sample</sub>",
          target_variable == "phosphorus_content" ~ "&mu;g<sub>P</sub> g<sup>-1</sup><sub>sample</sub>",
          target_variable == "sulfur_content" ~ "&mu;g<sub>S</sub> g<sup>-1</sup><sub>sample</sub>",
          target_variable == "potassium_content" ~ "&mu;g<sub>K</sub> g<sup>-1</sup><sub>sample</sub>",
          target_variable == "titanium_content" ~ "&mu;g<sub>Ti</sub> g<sup>-1</sup><sub>sample</sub>",
          target_variable == "silicon_content" ~ "g<sub>Si</sub> g<sup>-1</sup><sub>sample</sub>",
          target_variable == "calcium_content" ~ "g<sub>Ca</sub> g<sup>-1</sup><sub>sample</sub>",
          target_variable == "d13C" ~ "&permil;",
          target_variable == "d15N" ~ "&permil;",
          target_variable == "nosc" ~ "-",
          target_variable == "dgf0" ~ "kJ mol<sup>-1</sup><sub>C</sub>",
          target_variable == "loss_on_ignition" ~ "g g<sup>-1</sup><sub>sample</sub>",
          target_variable == "bulk_density" ~ "g cm<sup>-3</sup>",
          target_variable == "C_to_N" ~ "g<sub>C</sub> g<sup>-1</sup><sub>N</sub>",
          target_variable == "O_to_C" ~ "g<sub>O</sub> g<sup>-1</sup><sub>C</sub>",
          target_variable == "H_to_C" ~ "g<sub>H</sub> g<sup>-1</sup><sub>C</sub>" 
        ),
      unit_for_plotting_html_no_element_subscripts =
        dplyr::case_when(
          target_variable == "carbon_content" ~ "g g<sup>-1</sup><sub>sample</sub>",
          target_variable == "nitrogen_content" ~ "g g<sup>-1</sup><sub>sample</sub>",
          target_variable == "oxygen_content" ~ "g g<sup>-1</sup><sub>sample</sub>",
          target_variable == "hydrogen_content" ~ "g g<sup>-1</sup><sub>sample</sub>",
          target_variable == "phosphorus_content" ~ "&mu;g g<sup>-1</sup><sub>sample</sub>",
          target_variable == "sulfur_content" ~ "&mu;g g<sup>-1</sup><sub>sample</sub>",
          target_variable == "potassium_content" ~ "&mu;g g<sup>-1</sup><sub>sample</sub>",
          target_variable == "titanium_content" ~ "&mu;g g<sup>-1</sup><sub>sample</sub>",
          target_variable == "silicon_content" ~ "g g<sup>-1</sup><sub>sample</sub>",
          target_variable == "calcium_content" ~ "g g<sup>-1</sup><sub>sample</sub>",
          target_variable == "d13C" ~ "&permil;",
          target_variable == "d15N" ~ "&permil;",
          target_variable == "nosc" ~ "-",
          target_variable == "dgf0" ~ "kJ mol<sub>C</sub><sup>-1</sup>",
          target_variable == "loss_on_ignition" ~ "g g<sup>-1</sup><sub>sample</sub>",
          target_variable == "bulk_density" ~ "g cm<sup>-3</sup>",
          target_variable == "C_to_N" ~ "g g<sup>-1</sup>",
          target_variable == "O_to_C" ~ "g g<sup>-1</sup>",
          target_variable == "H_to_C" ~ "g g<sup>-1</sup>" 
        )
    )
  
  res[[what]]
  
}


#' Latex units
#' 
#' @export
irp_get_units_target_variables_latex <- function(x, what) {
  
  res <- 
    tibble::tibble(
      target_variable = x
    ) |>
    dplyr::mutate(
      unit_latex_no_element_subscripts =
        dplyr::case_when(
          target_variable == "carbon_content" ~ "g g$^{-1}_\\text{sample}$",
          target_variable == "nitrogen_content" ~ "g g$^{-1}_\\text{sample}$",
          target_variable == "oxygen_content" ~ "g g$^{-1}_\\text{sample}$",
          target_variable == "hydrogen_content" ~ "g g$^{-1}_\\text{sample}$",
          target_variable == "phosphorus_content" ~ "$\\mu$g g$^{-1}_\\text{sample}$",
          target_variable == "sulfur_content" ~ "$\\mu$g g$^{-1}_\\text{sample}$",
          target_variable == "potassium_content" ~ "$\\mu$g g$^{-1}_\\text{sample}$",
          target_variable == "titanium_content" ~ "$\\mu$g g$^{-1}_\\text{sample}$",
          target_variable == "silicon_content" ~ "g g$^{-1}_\\text{sample}$",
          target_variable == "calcium_content" ~ "g g$^{-1}_\\text{sample}$",
          target_variable == "d13C" ~ "$\\text{\\textperthousand}$",
          target_variable == "d15N" ~ "$\\text{\\textperthousand}$",
          target_variable == "nosc" ~ "-",
          target_variable == "dgf0" ~ "kJ mol$^{-1}_\\text{C}$",
          target_variable == "loss_on_ignition" ~ "g g$^{-1}_\\text{sample}$",
          target_variable == "bulk_density" ~ "g$_\\text{sample}$ cm$^{-3}_\\text{sample}$",
          target_variable == "C_to_N" ~ "g g$^{-1}$",
          target_variable == "O_to_C" ~ "g g$^{-1}$",
          target_variable == "H_to_C" ~ "g g$^{-1}$" 
        )
    )
  
  res[[what]]
  
}

#' Unit for irpeat
#' 
#' Same as`unit_for_modeling`, but with generic units and not elco units.
#' 
#' @export
irp_get_units_for_irpeat <- function(x) {
  
  res <- 
    tibble::tibble(
      target_variable = x
    ) |>
    dplyr::mutate(
      unit_for_irpeat = 
        dplyr::case_when(
          target_variable == "carbon_content" ~ "g/g",
          target_variable == "nitrogen_content" ~ "g/g",
          target_variable == "oxygen_content" ~ "g/g",
          target_variable == "hydrogen_content" ~ "g/g",
          target_variable == "phosphorus_content" ~ "g/g",
          target_variable == "sulfur_content" ~ "g/g",
          target_variable == "potassium_content" ~ "g/g",
          target_variable == "titanium_content" ~ "g/g",
          target_variable == "silicon_content" ~ "g/g",
          target_variable == "calcium_content" ~ "g/g",
          target_variable == "d13C" ~ "1",
          target_variable == "d15N" ~ "1",
          target_variable == "nosc" ~ "1",
          target_variable == "dgf0" ~ "J/mol",
          target_variable == "loss_on_ignition" ~ "g/g",
          target_variable == "bulk_density" ~ "g/cm^3",
          target_variable == "C_to_N" ~ "g/g",
          target_variable == "O_to_C" ~ "g/g",
          target_variable == "H_to_C" ~ "g/g" 
        )
    )
  
  what <- "unit_for_irpeat"
  res[[what]]
  
}


#' Target variable names for plotting and printing
#'
#' @param what Character value. One of `target_variable_name_html_short`, 
#' `target_variable_name_latex_short`, `target_variable_name_latex_long`.
#'
#' @export
irp_get_names_target_variables <- function(x, what) {
  
  res <- 
    tibble::tibble(
      target_variable = x
    ) |>
    dplyr::mutate(
      target_variable_name_html_short =
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
          target_variable == "d13C" ~ "&delta;<sup>13</sup>C",
          target_variable == "d15N" ~ "&delta;<sup>15</sup>N",
          target_variable == "nosc" ~ "NOSC",
          target_variable == "dgf0" ~ "&Delta;G<sup>0</sup><sub>f</sub>",
          target_variable == "loss_on_ignition" ~ "LOI",
          target_variable == "bulk_density" ~ "BD",
          target_variable == "C_to_N" ~ "C/N",
          target_variable == "O_to_C" ~ "O/C",
          target_variable == "H_to_C" ~ "H/C" 
        ),
      target_variable_name_latex_short =
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
          target_variable == "d13C" ~ "$\\delta^{13}$C",
          target_variable == "d15N" ~ "$\\delta^{15}$N",
          target_variable == "nosc" ~ "NOSC",
          target_variable == "dgf0" ~ "$\\Delta\\text{G}_\\text{f}^0$",
          target_variable == "loss_on_ignition" ~ "LOI",
          target_variable == "bulk_density" ~ "BD",
          target_variable == "C_to_N" ~ "C/N",
          target_variable == "O_to_C" ~ "O/C",
          target_variable == "H_to_C" ~ "H/C" 
        ),
      target_variable_name_latex_long =
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
          target_variable == "d13C" ~ "$\\delta^{13}$C",
          target_variable == "d15N" ~ "$\\delta^{15}$N",
          target_variable == "nosc" ~ "NOSC",
          target_variable == "dgf0" ~ "$\\Delta\\text{G}_\\text{f}^0$",
          target_variable == "loss_on_ignition" ~ "Loss on ignition",
          target_variable == "bulk_density" ~ "Bulk density",
          target_variable == "C_to_N" ~ "C/N",
          target_variable == "O_to_C" ~ "O/C",
          target_variable == "H_to_C" ~ "H/C" 
        )
    )
  
  res[[what]]
  
}


#' Returns the order of target variables
#' 
#' @export
irp_get_order_target_variables <- function(x) {
  
  res <- 
    tibble::tibble(
      target_variable = x
    ) |>
    dplyr::mutate(
      index =
        dplyr::case_when(
          target_variable == "carbon_content" ~ 1,
          target_variable == "nitrogen_content" ~ 3,
          target_variable == "oxygen_content" ~ 4,
          target_variable == "hydrogen_content" ~ 2,
          target_variable == "phosphorus_content" ~ 6,
          target_variable == "sulfur_content" ~ 5,
          target_variable == "potassium_content" ~ 7,
          target_variable == "titanium_content" ~ 10,
          target_variable == "silicon_content" ~ 8,
          target_variable == "calcium_content" ~ 9,
          target_variable == "d13C" ~ 11,
          target_variable == "d15N" ~ 12,
          target_variable == "nosc" ~ 13,
          target_variable == "dgf0" ~ 14,
          target_variable == "loss_on_ignition" ~ 19,
          target_variable == "bulk_density" ~ 18,
          target_variable == "C_to_N" ~ 15,
          target_variable == "O_to_C" ~ 16,
          target_variable == "H_to_C" ~ 17,
          target_variable == "macroporosity" ~ 20,
          target_variable == "non_macroporosity" ~ 21,
          target_variable == "volume_fraction_solids" ~ 22,
          target_variable == "saturated_hydraulic_conductivity" ~ 23,
          target_variable == "specific_heat_capacity" ~ 24,
          target_variable == "dry_thermal_conductivity" ~ 25
        )
    ) |>
    dplyr::pull(index)
  
  res
  
}

#' Defines the number of digits for rounding
#' 
#' @export
irp_get_rounding_digits_target_variables <- function(x) {
  
  res <- 
    tibble::tibble(
      target_variable = x
    ) |>
    dplyr::mutate(
      index =
        dplyr::case_when(
          target_variable == "carbon_content" ~ 2,
          target_variable == "nitrogen_content" ~ 3,
          target_variable == "oxygen_content" ~ 2,
          target_variable == "hydrogen_content" ~ 2,
          target_variable == "phosphorus_content" ~ 0,
          target_variable == "sulfur_content" ~ 0,
          target_variable == "potassium_content" ~ 0,
          target_variable == "titanium_content" ~ 0,
          target_variable == "silicon_content" ~ 2,
          target_variable == "calcium_content" ~ 3,
          target_variable == "d13C" ~ 1,
          target_variable == "d15N" ~ 1,
          target_variable == "nosc" ~ 1,
          target_variable == "dgf0" ~ 1,
          target_variable == "loss_on_ignition" ~ 2,
          target_variable == "bulk_density" ~ 2,
          target_variable == "C_to_N" ~ 1,
          target_variable == "O_to_C" ~ 3,
          target_variable == "H_to_C" ~ 3 
        )
    ) |>
    dplyr::pull(index)
  
  res
  
}



#### Isotope values ####

#' Converts delta values to corresponding atom percentages
#'
#' @param x A numeric vector with values representing isotope ratios in the
#' delta notation for a specific isotope.
#'
#' @param r A numeric value representing the delta value for the corresponding
#' isotope standard.
#'
#' @return `x` with values converted to atom percentages (i.e. a value in
#' \\[0,1\\]).
#'
#' @seealso `irp_atom_fraction_to_delta()`
#'
#' @examples
#' # convert a d13C value to atom percent
#' irp_delta_to_atom_fraction(x = -25, r = 0.0112372)
#'
#' @export
irp_delta_to_atom_fraction <- function(x, r) {
  
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(r) && length(r) == 1L)
  
  (r * (x/1000 + 1))/(1 + r * (x/1000 + 1))
}

#' Converts atom percentages to corresponding delta values
#'
#' @param x A numeric vector with values representing the atom percentage of a
#' specific isotope.
#'
#' @param r A numeric value representing the delta value for the corresponding
#' isotope standard.
#'
#' @return `x` with values converted to delta values.
#'
#' @seealso `irp_delta_to_atom_fraction()`
#'
#' @examples
#' # convert a 13C atom percentage to delta values
#' irp_atom_fraction_to_delta(
#'   x = irp_delta_to_atom_fraction(x = -25, r = 0.0112372),
#'   r = 0.0112372
#' )
#'
#' @export
irp_atom_fraction_to_delta <- function(x, r) {
  
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(r) && length(r) == 1L)
  
  x <- x*100
  -1000*(1-x/((100-x)*r))
}


#### Helper functions for the paper ####

#' Formats author list
#' 
#' @export
irp_make_author_list <- function(x) {
  
  purrr::map_chr(x, function(.x) {
    paste0(.x$name, "$^{\\text{", .x$affil, "}}$")
  }) |>
    paste(collapse = "  \n")
  
}

#' Formats affiliation list
#' 
#' @export
irp_make_affiliation_list <- function(x) {
  
  purrr::map_chr(x, function(.x) {
    paste0("$^{", .x$number, "}$ ", .x$text)
  }) |>
    paste(collapse = "  \n")
  
}


