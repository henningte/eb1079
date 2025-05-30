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
    iter = 4000L,
    warmup = 2000L,
    chains = 4L,
    adapt_delta = 0.9,
    max_treedepth = 12
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
        variable == "d13C" ~ "$\\delta^{13}$C value of bulk peat relative to Vienna Pee Dee Bee standard.",
        variable == "d15N" ~ "$\\delta^{15}$N value of bulk peat relative to Air standard.",
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
        variable == "silocon_content" ~ "Estimating peat mineral inputs.",
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


