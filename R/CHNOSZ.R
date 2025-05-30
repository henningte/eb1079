#### Modified functions ####

#' Modified version of the `OBIGT()` function from the 'CHNOSZ' package to load data from the OBIGT database
#'
#' This is a modified version of the `OBIGT()` function from the 'CHNOSZ'
#' package to load data from the OBIGT database. The only modification is that
#' file names of the `csv` files are added for each row as a new column
#' `file`. This allows more efficient filtering for organic compounds.
#'
#' @inheritParams CHNOSZ::OBIGT
#'
#' @source [`CHNOSZ::OBIGT()`]; See:
#' https://cran.r-project.org/web/packages/CHNOSZ/index.html. (code licensed
#' under the GPL-3 license).
#'
#' @examples
#' irp_get_OBIGT()
#'
#' @export
irpp_get_OBIGT <- function (no.organics = FALSE) {
  if (!"thermo" %in% ls(CHNOSZ)) 
    stop("The CHNOSZ environment doesn't have a \"thermo\" object. Try running reset()")
  sources_aq <- paste0(c("H2O", "inorganic", "organic"), "_aq")
  sources_cr <- paste0(c("Berman", "inorganic", "organic"), 
                       "_cr")
  sources_liq <- paste0(c("organic"), "_liq")
  sources_gas <- paste0(c("inorganic", "organic"), "_gas")
  sources <- c(sources_aq, sources_cr, sources_gas, sources_liq)
  OBIGTdir <- system.file("extdata/OBIGT/", package = "CHNOSZ")
  datalist <- 
    lapply(sources, function(source) {
      file <- paste0(OBIGTdir, "/", source, ".csv")
      dat <- read.csv(file, as.is = TRUE)
      dat$file <- source
      if (no.organics & grepl("^organic", source)) 
        dat <- subset(dat, formula == "CH4")
      dat
    })
  OBIGT <- do.call(rbind, datalist)
  refs <- read.csv(file.path(OBIGTdir, "refs.csv"), as.is = TRUE)
  thermo <- get("thermo", CHNOSZ)
  thermo$OBIGT <- OBIGT
  thermo$refs <- refs
  assign("thermo", thermo, CHNOSZ)
  message(paste("OBIGT: loading", ifelse(no.organics,
                                         "inorganic", "default"), "database with",
                nrow(thermo$OBIGT[thermo$OBIGT$state == "aq", ]),
                "aqueous,", nrow(thermo$OBIGT), "total species"))
  idup <- duplicated(paste(thermo$OBIGT$name, thermo$OBIGT$state))
  if (any(idup))
    warning("OBIGT: duplicated species: ", paste(thermo$OBIGT$name[idup],
                                                 "(", thermo$OBIGT$state[idup], ")", sep = "",
                                                 collapse = " "))
}



#' Modified version of the `add.OBIGT()` function from the 'CHNOSZ' package to add data to the OBIGT database
#'
#' This is a modified version of the `add.OBIGT()` function from the
#' 'CHNOSZ' package to add data to the OBIGT database. The only modification is
#' that the additional column `file` in the OBIGT database (See
#' [`irpp_get_OBIGT()`]) is also considered when adding new data.
#'
#' @inheritParams CHNOSZ::add.OBIGT
#'
#' @source [`CHNOSZ::add.OBIGT()`]; See:
#' https://cran.r-project.org/web/packages/CHNOSZ/index.html. (code licensed
#' under the GPL-3 license).
#'
#' @export
irpp_add_OBIGT <-
  function (file, species = NULL, force = TRUE) {
    thermo <- get("thermo", CHNOSZ)
    to1 <- thermo$OBIGT
    id1 <- paste(to1$name, to1$state)
    sysfiles <- dir(system.file("extdata/OBIGT/", package = "CHNOSZ"))
    sysnosuffix <- sapply(strsplit(sysfiles, "\\."), "[", 1)
    isys <- match(file, sysnosuffix)
    if (!is.na(isys))
      file <- system.file(paste0("extdata/OBIGT/", sysfiles[isys]),
                          package = "CHNOSZ")
    to2 <- read.csv(file, as.is = TRUE)
    if (!"E_units" %in% colnames(to2))
      to2 <- data.frame(to2[, 1:7], E_units = "cal",
                        to2[, 8:ncol(to2)], stringsAsFactors = FALSE)
    Etxt <- paste(unique(to2$E_units), collapse = " and ")
    if (!is.null(species)) {
      idat <- match(species, to2$name)
      ina <- is.na(idat)
      if (!any(ina))
        to2 <- to2[idat, ]
      else stop(paste("file", file, "doesn't have",
                      paste(species[ina], collapse = ", ")))
    }
    id2 <- paste(to2$name, to2$state)
    tr <- tryCatch(rbind(to1, to2), error = identity)
    if (inherits(tr, "error"))
      stop(paste(file, "is not compatible with thermo$OBIGT data table."))
    does.exist <- id2 %in% id1
    ispecies.exist <- na.omit(match(id2, id1))
    nexist <- sum(does.exist)
    inew <- numeric()
    if (force) {
      if (nexist > 0) {
        to1[ispecies.exist, ] <- to2[does.exist, ]
        to2 <- to2[!does.exist, ]
        inew <- c(inew, ispecies.exist)
      }
    }
    else {
      to2 <- to2[!does.exist, ]
      nexist <- 0
    }
    if (nrow(to2) > 0) {
      to1 <- rbind(to1, to2)
      inew <- c(inew, (length(id1) + 1):nrow(to1))
    }
    thermo$OBIGT <- to1
    rownames(thermo$OBIGT) <- 1:nrow(thermo$OBIGT)
    assign("thermo", thermo, CHNOSZ)
    message("add.OBIGT: read ", length(does.exist), " rows; made ",
            nexist, " replacements, ", nrow(to2), " additions [energy units: ",
            Etxt, "]")
    return(invisible(inew))
  }


#### Own functions ####



#' Parses chemical formulas given as string into a data frame
#'
#' `irp_cf_to_df` takes a string of chemical formulas and parses this into a
#' data frame where each row corresponds to a chemical formula, each column to a
#' chemical element and each cell value to an atom count per molecule of a
#' compound.
#'
#' @param x A string
#'
#' @return A data frame with a row for each element in `x`, a column for each
#' unique chemical element in \code{x}, and cell values representing atom counts
#' for each compound and chemcial element.
#'
#' @export
irp_cf_to_df <- function(x) {
  
  # get unique elements
  x_el <-
    x |>
    stringr::str_extract_all(pattern = "[A-Z]{1}[a-z]*") |>
    unlist() |>
    unique()
  
  # get counts
  purrr::map_dfc(x_el, function(i) {
    
    res <-
      x |>
      stringr::str_extract(pattern = paste0(i, "[0-9]+\\.*[0-9]*")) |>
      stringr::str_remove(pattern = i) |>
      as.numeric()
    
    res[is.na(res)] <- 1
    
    res <-
      tibble::tibble(
        i = ifelse(stringr::str_detect(x, pattern = i), res, 0)
      )
    
    colnames(res) <- i
    res
    
  })
  
}



#' Batch computing enthalpies of combustion with the 'OBIGT' database
#'
#' `irpp_get_enthalpy_of_combustion` is a wrapper for [`CHNOSZ::subcrt()`] to
#' support batch computing of enthalpies of combustion for multiple compounds.
#'
#' @param x A character vector where each element is a chemical formula for a
#' compound for which to compute combustion enthalpies.
#'
#' @param temperature,pressure See [`CHNOSZ::subcrt()`] (argument `T` and `P`,
#' respectively).
#'
#' @seealso
#' [`CHNOSZ::subcrt()`]
#'
#' @export
irpp_get_enthalpy_of_combustion <- function(x, state, temperature = 298.15, pressure = 1) {
  
  purrr::map_dbl(seq_along(x), function(i) {
    res <- 
      tryCatch(
        suppressMessages(CHNOSZ::subcrt(x[[i]], c(-1), state = state[[i]], T = temperature, P = pressure, autobalance = TRUE)),
        error = function(cond) list(out = list(H = NA_real_))
      )
    res$out$H
  })
  
}


#' Computes the number of electrons transferred during combustion of a compound
#'
#' Computes the number of electrons transferred during combustion of a compound
#' ($E$) following \insertCite{Battley.1998;textual}{ir} (see also
#' \insertCite{Popovic.2019;textual}{ir}).
#'
#' @param C,H,O,N,P,S A units vector representing the mols of the respective 
#' elements in the compounds to combust. The unit must be defined based on the 
#' elco package, for example, C needs unit `mol_C`.
#'
#' @param coef_S A numeric value representing the number of electrons removed
#' from S during combustion. `coef_S = 4` corresponds to combustion to SO$_2$.
#' `coef_S = 6` corresponds to combustion to SO$_3$ (or H$_2$SO$_4$). See
#' \insertCite{Popovic.2019;textual}{ir} for details.
#'
#' @return A units vector with the amount of electrons transferred during
#' combustion of each compound \[mol$_\text{e}^{-1}$ mol$^{-1}$\].
#'
#' @source
#' \insertAllCited{}
#'
#' @examples
#' # with units objects
#' library(elco)
#' elco_get_electrons_tansferred_combustion(
#'   C = 0.04 |>
#'     quantities::set_quantities(unit = "mol_C", errors = 0),
#'   H = 0.05 |>
#'     quantities::set_quantities(unit = "mol_H", errors = 0),
#'   O = 0.02 |>
#'     quantities::set_quantities(unit = "mol_O", errors = 0),
#'   N = 0.0015 |>
#'     quantities::set_quantities(unit = "mol_N", errors = 0),
#'   P = 1.3e-5 |>
#'     quantities::set_quantities(unit = "mol_P", errors = 0),
#'   S = 1e-4 |>
#'     quantities::set_quantities(unit = "mol_S", errors = 0)
#' )
#' 
#' elco_get_electrons_tansferred_combustion(
#'   C = 0.04 |>
#'     units::set_units("mol_C"),
#'   H = 0.05 |>
#'     units::set_units("mol_H"),
#'   O = 0.02 |>
#'     units::set_units("mol_O"),
#'   N = 0.0015 |>
#'     units::set_units("mol_N"),
#'   P = 1.3e-5 |>
#'     units::set_units("mol_P"),
#'   S = 1e-4 |>
#'     units::set_units("mol_S")
#' )
#'
#' @export
elco_get_electrons_tansferred_combustion <-
  function(C = 0, H = 0, O = 0, N = 0, P = 0, S = 0, coef_S = 4) {
    
    stopifnot(length(coef_S) == 1 && is.numeric(coef_S))
    
    coefs <- list(C = 4, H = 1, O = -2, N = -0, P = 5, S = coef_S)
    for(i in seq_along(coefs)) {
      coefs[[i]] <- quantities::set_quantities(coefs[[i]], unit = paste0("mol/mol_", names(coefs)[[i]]), errors = 0, mode = "standard")
    }
    
    all_elements <- list(C = C, H = H, O = O, N = N, P = P, S = S)
    cond <- purrr::map_lgl(all_elements, inherits, what = "units")
    if(! all(cond)) {
      rlang::abort('All of `C,H,N,O,S,P` must be given as `units` objects.')
    } 
    cond <- purrr::map_lgl(all_elements, inherits, what = "errors")
    if(! any(cond)) {
      coefs <- purrr::map(coefs, errors::drop_errors)
    }
    
    # check that all are convertible to mol
    all_elements <- purrr::map(all_elements, elco::elco_convert, to = "mol")
    
    purrr::map2(coefs, all_elements, magrittr::multiply_by) |>
      purrr::reduce(.f = magrittr::add)
    
  }



#' Computes the nominal oxidation state of carbon (NOSC) for a compound
#'
#' `irpp_get_nosc` computes the nominal oxidation state of carbon (NOSC) for a
#' compound using the equation stated in @LaRowe.2011.
#'
#' @param Z Net charge of the compound.
#' @param C,H,N,O,P,S Mols of the respective element in the compound.
#'
#' @examples
#' # providing values as quantities/elco objects
#' elco_get_nosc_larowe2011(
#'   Z = quantities::set_quantities(0, unit = "mol", errors = 0),
#'   C = 1 |>
#'     quantities::set_quantities(unit = "mol_C", errors = 0),
#'   H = 1.77 |>
#'     quantities::set_quantities(unit = "mol_H", errors = 0),
#'   O = 0.49 |>
#'     quantities::set_quantities(unit = "mol_O", errors = 0),
#'   N = 0.24 |>
#'     quantities::set_quantities(unit = "mol_N", errors = 0),
#'   P = 0 |>
#'     quantities::set_quantities(unit = "mol_P", errors = 0),
#'   S = 0 |>
#'     quantities::set_quantities(unit = "mol_S", errors = 0)
#' )
#' 
#' elco_get_nosc_larowe2011(
#'   Z = units::set_units(0, "mol"),
#'   C = 1 |>
#'     units::set_units("mol_C"),
#'   H = 1.77 |>
#'     units::set_units("mol_H"),
#'   O = 0.49 |>
#'     units::set_units("mol_O"),
#'   N = 0.24 |>
#'     units::set_units("mol_N"),
#'   P = 0 |>
#'     units::set_units("mol_P"),
#'   S = 0 |>
#'     units::set_units("mol_S")
#' ) 
#'
#' @export
elco_get_nosc_larowe2011 <- function(Z, C, H = 0, N = 0, O = 0, P = 0, S = 0) {
  
  coefs <- list(Z = -1, C = 4, H = 1, N = -3, O = -2, P = +5, S = -2, multiplication_factor = -1)
  for(i in seq_along(coefs)) {
    unit <- 
      if(names(coefs)[[i]] %in% c("Z", "multiplication_factor")) {
        "1"
      } else {
        paste0("mol/mol_", names(coefs)[[i]])
      }
    coefs[[i]] <- quantities::set_quantities(coefs[[i]], unit = unit, errors = 0, mode = "standard")
  }
  
  all_elements <-
    tibble::tibble(
      Z = Z,
      C = C,
      H = H,
      N = N,
      O = O,
      P = P,
      S = S
    )
  
  cond <- purrr::map_lgl(all_elements, inherits, what = "units")
  if(! all(cond)) {
    rlang::abort('All of `Z,C,H,N,O,S,P` must be given as `units` objects.')
  } 
  cond <- purrr::map_lgl(all_elements, inherits, what = "errors")
  if(! any(cond)) {
    coefs <- purrr::map(coefs, errors::drop_errors)
  }
  
  # check that all are convertible to mol
  purrr::map(all_elements[, -1], elco::elco_convert, to = "mol")
  units::set_units(Z, "mol", mode = "standard")
  
  purrr::map2(all_elements, coefs[-length(coefs)], magrittr::multiply_by) |>
    purrr::reduce(.f = magrittr::add) |>
    magrittr::multiply_by(coefs$multiplication_factor) |>
    magrittr::divide_by(elco::elco_make_units_generic(all_elements$C)) |>
    magrittr::add(elco::elco_make_units_generic(coefs$C))
  
}




