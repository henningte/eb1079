#' Prepares data from the OBIGT database
#' 
#' @export
irp_make_db_obigt_organic <- function() {
  
  # define units to use
  CHNOSZ::T.units("K")
  CHNOSZ::P.units("bar")
  CHNOSZ::E.units("J")
  
  CHNOSZ::OBIGT(no.organics = FALSE)
  
  # load the OBIGT database
  irpp_get_OBIGT()
  
  ## convert units to cal/mol since this seems the unit assumed by add.OBIGT
  #read.csv("data/raw_data/CHNOSZ/to_add_to_OBIGT_J.csv", as.is = TRUE) |>
  #  dplyr::mutate(
  #    H =
  #      H |>
  #      units::set_units(value = "J/mol") |>
  #      units::set_units(value = "cal/mol") |>
  #      units::drop_units(),
  #    E_units = "cal"
  #  ) |>
  #  write.csv("data/raw_data/CHNOSZ/to_add_to_OBIGT_cal.csv", row.names = FALSE)
  
  # add
  #to_add_to_OBIGT <- irpp_add_OBIGT("data/raw_data/CHNOSZ/to_add_to_OBIGT_J.csv")
  CHNOSZ::mod.OBIGT("tetraphosphorus decaoxide", formula = "P4O10", state = "cr", ref1 = "NIST", G = NA, H = -3009.94, S = 228.82, E_units = "J", z.T = 0)
  
  # generate the target subset of the OBIGT database
  db_obigt_organic <- 
    thermo()$OBIGT |>
    dplyr::filter((stringr::str_detect(file, pattern = "(organic|biotic)") & !stringr::str_detect(file, pattern = "inorganic")) | name %in% c("CH4", "carbon dioxide", "graphite")) |> # discard inorganic species, except for methane
    dplyr::filter(state != "aq") |> # discard dissolved species
    dplyr::filter(! stringr::str_detect(formula, pattern = "\\:")) |> # discard chemical formulas that are difficult to parse (these are only 3 entries and therefore it is unlikely that these would have a large impact)
    dplyr::filter(! is.na(H)) |> # we need species with a value for the enthalpy of formation to compute the enthalpy of combustion
    dplyr::filter(! stringr::str_detect(formula, pattern = "[A-BD-GI-MQRT-Z]")) |> # only retain compounds that contain only C, H, N, O, S, P
    dplyr::filter(! stringr::str_detect(name, pattern = "\\[.?(-|>)+")) |> # discard compounds with non-standard chemical formula (13 entries)
    dplyr::filter(! stringr::str_detect(name, pattern = "\\[")) |> # the amino acids have strange chemical formulas
    dplyr::filter(! stringr::str_detect(formula, pattern = "(C.*C|O.*O|H.*H)")) # remove all compounds where the chemical formula is not a simple stoichiometric formula, but has e.g. multiple entries for "C" (note: this are only 97 formulas and is not likely to cause large bias)
  
  # Define the basis species
  CHNOSZ::basis(c("CO2", "H2O", "N2", "SO2", "O2", "P4O10", "K2O", "MgO", "CaO", "Fe2O3"), state = c("gas", "liq", "gas", "gas", "gas", "cr", "cr", "cr", "cr", "cr"), add = TRUE) # states from @Battley.1999
  
  # extract stoichiometry
  db_obigt_organic_cf <- 
    irp_cf_to_df(db_obigt_organic$formula) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::everything(),
        function(.x) {
          units::set_units(.x, paste0("mol_", dplyr::cur_column()), mode = "standard")
        }
      )
    )
  
  irp_standard_entropies_elements <- 
    irp_make_standard_entropies_elements() |>
    dplyr::filter(chemical_element %in% colnames(db_obigt_organic_cf)) |>
    dplyr::slice(match(chemical_element, colnames(db_obigt_organic_cf)))
  
  # compute enthalpy of combustion
  db_obigt_organic <-
    db_obigt_organic |>
    # dplyr::filter(!ref1 %in% c("HOKR98", "HOKR98.1", "RH98")) |>
    dplyr::mutate(
      enthalpy_of_combustion = 
        irpp_get_enthalpy_of_combustion(
          x = name, 
          state = state, 
          temperature = 298.15, 
          pressure = 1
        ), #---todo: this are standard conditions often used for thermodynamic calculations (https://en.wikipedia.org/wiki/Standard_state). However, these are not IUPAC standard conditions, where T = 273.15 K. I have to check this.
      enthalpy_of_combustion_per_C = enthalpy_of_combustion / db_obigt_organic_cf$C,
      E = 
        with(
          db_obigt_organic_cf, 
          elco_get_electrons_tansferred_combustion(C = C, H = H, N = N, O = O, P = P, S = S)
        ),
      E_per_C = E/units::set_units(db_obigt_organic_cf$C, mol_C),
      nosc = 
        with(
          db_obigt_organic_cf, 
          elco_get_nosc_larowe2011(
            Z = units::set_units(0, "mol", mode = "standard"), 
            C = C, H = H, N = N, O = O, P = P, S = S
          )
        ),
      sum_of_element_standard_entropies =
        units::set_units(as.vector(as.matrix(db_obigt_organic_cf) %*% irp_standard_entropies_elements$standard_entropy), "J/(K * mol)", mode = "standard"),
      sum_of_element_standard_entropies_per_C = sum_of_element_standard_entropies/(db_obigt_organic_cf$C / units::set_units(1, "mol", mode = "standard")),
      entropy_of_formation =
        units::set_units(S, paste0(unique(E_units), "/K/mol"), mode = "standard") |>
        units::set_units("J/K/mol", mode = "standard") |>
        magrittr::subtract(sum_of_element_standard_entropies),
      entropy_of_formation_per_C = entropy_of_formation/(db_obigt_organic_cf$C / units::set_units(1, "mol", mode = "standard"))
    ) |>
    dplyr::filter(as.numeric(db_obigt_organic_cf$C) > 0 & ! is.na(enthalpy_of_combustion))
  
  db_obigt_organic
  
}



#' Prepares standard entropies for chemical elements from @Battley.1999
#'
#' Data for standard entropies for chemical elements extracted and reformatted
#' from @Battley.1999.
#'
#' @format A data frame with 11 rows and 2 columns. Each row represents one
#' chemical element:
#' \describe{
#'   \item{`chemical_element`}{A character vector with symbols for chemical
#'     elements.}
#'   \item{`standard_entropy`}{A numeric vector with standard entropy values for
#'     each chemical element \[J K$^{-1}$ mol$^{-1}$\].}
#' }
#'
#' @source
#' Data extracted from @Battley.1999.
#'
#' @export
irp_make_standard_entropies_elements <- function() {
  
  tibble::tibble(
    chemical_element = c("C", "H", "N", "O", "S", "P", "Zn", "H2O", "Ca", "Mg", "K"),
    standard_entropy = 
      c(5.74, 65.34, 95.81, 102.57, 31.8, 41.09, 41.72, 69.95, 41.59, 32.67, 64.63) |>
      units::set_units("J/(K * mol)", mode = "standard")
  )
  
} 


#' Data from table 1 from @Battley.1999
#'
#' Data extracted and reformatted from table 1 from @Battley.1999. Data
#' were extracted using [`tabulizer::extract_tables()`].
#'
#' @format A data frame with 23 rows and 6 columns. Each row is one chemical
#' compound or microbial biomass sample:
#' \describe{
#'   \item{`name`}{A character vector with names of the analyzed chemical
#'     compound or biomass sample.}
#'   \item{`H2O`}{A numeric vector with the mols of H2O per mol of the compound
#'     or per mol C of the biomass sample \[mol mol$^{-1}$\].}
#'   \item{`formula`}{A character vector with the chemical formula of the
#'     chemical compound or microbial biomass normalized per mol of C.}
#'   \item{`entropy_of_formation`}{A numeric vector with entropies of formation
#'     derived from measurements \[J K$^{-1}$ mol$^{-1}$\].}
#'   \item{`sum_of_element_standard_entropies`}{A numeric vector with the sum of
#'     the standard entropies of the chemical elements in the compound or sample
#'     multiplied by their respective amounts \[mol\] in one mol (for chemcial
#'     compounds, or mol$_\text{C}$ in case of biomass samples)
#'     \[J K$^{-1}$ mol$^{-1}$\].}
#'   \item{`standard_entropy`}{A numeric vector with standard entropies for one
#'     mol of the the respective compounds or one mol$_\text{C}$ of the
#'     respective biomass samples, derived from experimental measurements.}
#' }
#'
#' @source
#' Data are extracted from @Battley.1999.
#' 
#' @export
irp_make_d_battley_1999_table_1 <- function() {
  
  res <- 
    tabulapdf::extract_tables("private/papers/Battley - 1999 - An empirical method for estimating the entropy of .pdf")[[1]] |>
    as.data.frame() |>
    dplyr::slice(-c(1:2))
  
  res$`1 2`[[20]] <- paste0(res$`1 2`[[20]], res$`1 2`[[21]])
  res$`1 2`[[23]] <- paste0(res$`1 2`[[23]], res$`1 2`[[24]])
  res$`1 2`[[25]] <- paste0(res$`1 2`[[25]], res$`1 2`[[26]])
  
  res <- 
  res |>
    dplyr::filter(`4` != "") |>
    dplyr::mutate(
      name =
        `1 2` |> 
        stringr::str_remove(pattern = "^[a-z]{1}"),
      H2O =
        dplyr::case_when(
          stringr::str_detect(name, "H2O") ~ stringr::str_extract(name, pattern = " \\d+\\.{1}\\d*H2O$"),
          TRUE ~ "0"
        ) |>
        stringr::str_remove(pattern = "H2O") |>
        as.numeric(),
      formula =
        name |>
        stringr::str_remove(pattern = " \\d+\\.{1}\\d*H2O$") |>
        stringr::str_split(pattern = " C") |>
        purrr::map_chr(function(x) paste0("C", x[[2]])),
      entropy_of_formation =
        `5` |>
        stringr::str_remove(pattern = "ÿ") |>
        as.numeric() |>
        magrittr::multiply_by(-1) |>
        units::set_units("J/K/mol", mode = "standard"),
      sum_of_element_standard_entropies = 
        as.numeric(`6`) |>
        units::set_units("J/K/mol", mode = "standard"),
      standard_entropy = 
        as.numeric(`4`) |>
        units::set_units("J/K/mol", mode = "standard"),
      name =
        name |>
        stringr::str_remove(pattern = " \\d+\\.{1}\\d*H2O$") |>
        stringr::str_split(pattern = " C") |>
        purrr::map_chr(function(x) x[[1]])
    ) |>
    dplyr::select(-c(1:12))
  
  # normalize by C content
  d_elements_standard_entropies <- irp_make_standard_entropies_elements()
  res_cf <-
    irp_cf_to_df(res$formula) |>
    dplyr::mutate(
      H2O = res$H2O,
      O = O + H2O,
      H = H + 2 * H2O
    ) |>
    dplyr::select(dplyr::any_of(d_elements_standard_entropies$chemical_element))  |>
    dplyr::select(-H2O)
  
  # entropy_of_formation_per_C
  res <-
    res |>
    dplyr::mutate(
      entropy_of_formation_per_C = entropy_of_formation/units::set_units(res_cf$C, "mol_C/mol", mode = "standard"),
      sum_of_element_standard_entropies_per_C = sum_of_element_standard_entropies/units::set_units(res_cf$C, "mol_C/mol", mode = "standard")
    )
  
  res
  
}


#' Generates a table with standard enthalpies of formation for species created when combusting OM
#' 
#' @export
irp_make_d_compounds_standard_enthalpies_of_formation <- function() {
  
  # define units to use
  CHNOSZ::T.units("K")
  CHNOSZ::P.units("bar")
  CHNOSZ::E.units("J")
  
  CHNOSZ::OBIGT(no.organics = FALSE)
  
  # load the OBIGT database
  irpp_get_OBIGT()
  
  # Define the basis species
  CHNOSZ::basis(c("CO2", "H2O", "N2", "SO2", "O2", "O10P4", "K2O", "MgO", "CaO", "Fe2O3"), state = c("gas", "liq", "gas", "gas", "gas", "cr", "cr", "cr", "cr", "cr"), add = TRUE, delete = TRUE) # states from @Battley.1999
  
  d_compounds_standard_enthalpies_of_formation <-
    tibble::tibble(
      chemical_compound = c("CO2", "H2O", "N2", "SO2", "O2", "K2O", "MgO", "CaO", "Fe2O3"),
      state = c("gas", "liq", "gas", "gas", "gas", "cr", "cr", "cr", "cr"),
      standard_enthalpy_of_formation =
        purrr::map2_dbl(chemical_compound, state, function(.x, .y) {
          res <- CHNOSZ::subcrt(.x, T = 298.15, P = 1, state = .y, property = "H")
          res$out[[1]]$H
        }) |>
        quantities::set_quantities(unit = "J/mol", errors = 0)
    ) |>
    dplyr::bind_rows(
      tibble::tibble(
        chemical_compound = "P4O10",
        state = "cr",
        standard_enthalpy_of_formation = quantities::set_quantities(-3009.94, unit = "kJ/mol", errors = 0) #---note: source: https://webbook.nist.gov/cgi/cbook.cgi?ID=C16752606&Mask=6F
      )
    )
  
}  


#' Imports data from Wang.2015b
#' 
#' @export
irp_get_wang2015b <- function() {
  readRDS("data/raw_data/d9.rds")
}



#### Data from the pmird database ####

#' Gets MIRS from the pmird database
#' 
#' @export
irp_get_pmird_mirs <- function() {
  
  # get connection to database
  con <- irp_get_pmird_connection()
  on.exit(RMariaDB::dbDisconnect(con))
  
  dm_pmird <-
    pmird::pm_get_dm(con, learn_keys = TRUE)
  
  # get spectra
  res <-
    dm_pmird |>
    dm::dm_zoom_to(data) |>
    dm::left_join(mir_metadata, by = "id_measurement") |>
    dm::left_join(data_to_samples, by = "id_measurement") |>
    dm::left_join(samples, by = "id_sample") |>
    dm::pull_tbl() |>
    tibble::as_tibble()

  # load spectra
  res <-
    res |>
    dplyr::filter(!is.na(mirs_file)) |>
    dplyr::mutate(
      sample_id = as.character(id_measurement)
    ) |>
    pmird::pm_load_spectra(.col_sample_id = "sample_id", directory = "data/raw_data/pmird/")
  
  # remove corrupted spectra
  res <- 
    res |>
    dplyr::group_split(id_dataset) |>
    purrr::map(function(.x) {
      .x |>
        dplyr::mutate(
          id_measurement_for_dataset = seq_along(id_measurement)
        )
    }) |>
    dplyr::bind_rows() |>
    dplyr::filter(! (id_dataset == 16 & ! id_measurement_for_dataset %in% 1:49))
  
  # remove NA in intensities, correct measurement device name, recompute wavenumber range
  res <- 
    res |>
    dplyr::mutate(
      spectra =
        purrr::map(spectra, function(.x) {
          .x |>
            dplyr::filter(! is.na(y))
        }),
      measurement_device = 
        measurement_device |> 
        stringr::str_replace_all(pattern = "\\n", replacement = " ")
    ) |>
    dplyr::select(! dplyr::all_of(c("x_variable_min", "x_variable_max"))) |>
    range(
      .dimension = "x", 
      .col_names = c("x_variable_min", "x_variable_max")
    )
 
  res
  
}


#' Gets data from Liu.2019 from the pmird database
#'
#'@export
irp_get_liu2019_pmird <- function() {
  
  con <- irp_get_pmird_connection()
  on.exit(RMariaDB::dbDisconnect(con))
  
  dm_pmird <-
    pmird::pm_get_dm(con, learn_keys = TRUE)
  
  res <-
    dm_pmird |>
    dm::dm_zoom_to(data) |>
    dm::left_join(mir_metadata, by = "id_measurement") |>
    dm::left_join(data_to_samples, by = "id_measurement") |>
    dm::left_join(samples, by = "id_sample") |>
    dm::pull_tbl() |>
    tibble::as_tibble() |>
    dplyr::filter(id_dataset == 25)
  
  res
  
}


#' Gets pmird units for attributes
#' 
#' @export
irp_get_pmird_units <- function() {
  
  con <- irp_get_pmird_connection()
  on.exit(RMariaDB::dbDisconnect(con))
  
  dm_pmird <-
    pmird::pm_get_dm(con, learn_keys = TRUE)
  
  pmird_units <- 
    dm_pmird %>%
    dm::dm_zoom_to(attributes) %>%
    dm::left_join(measurement_scales, by = "id_measurement_scale") %>%
    dm::left_join(measurement_scales_ratio, by = "id_measurement_scale") %>%
    dm::left_join(units, by = "id_unit") %>%
    dm::pull_tbl() %>%
    tibble::as_tibble()
  
}


#' Gets references to cite from the pmird database
#' 
#' @export
irp_get_pmird_citations <- function(irp_d_model_info_enriched_1) {
  
  con <- irp_get_pmird_connection()
  on.exit(RMariaDB::dbDisconnect(con))
  
  id_measurement <- 
    irp_d_model_info_enriched_1$id_measurement_all |>
    unlist() |>
    unique()
  
  pmird::pm_get_citations(con = con, x = id_measurement) |>
    dplyr::pull(BIBTEXKEY) |>
    unique()
  
}


#### Misc ####


#' Data from Gnatowski.2022
#' 
#' @export
irp_make_d_gnatowski2022 <- function() {
  
  res1 <- 
    read.csv("data/raw_data/specific_heat_capacity/Gnatowski.2022-Tab1.csv") |>
    dplyr::mutate(
      sample_label = paste0("P", id_sample),
      dplyr::across(dplyr::any_of(c("C", "N", "ash_content")), function(.x) .x/100)
    )
  
  res2 <- 
    read.csv("data/raw_data/specific_heat_capacity/applsci-1523717-supplementary.csv", sep = ";") |>
    dplyr::slice(-1) |>
    dplyr::select(1, 7:9) |>
    setNames(nm = c("sample_label", "temperature", "specific_heat_capacity", "specific_heat_capacity_error")) |>
    dplyr::mutate(
      dplyr::across(-1, as.numeric)  
    )
  
  res <- 
    res1 |>
    dplyr::left_join(
      res2, 
      by = "sample_label"
    )
  
  res
  
}


#' Data on specific heat capacity from Hilsenrath.1955
#' 
#' @export
irp_make_d_cp_air <- function() {
  
  read.csv("data/raw_data/specific_heat_capacity/specific_heat_capacity_air_nbscircular564.csv") %>%
    dplyr::mutate(
      specific_heat_capacity = cp_R * 0.287041,
      specific_heat_capacity_error = specific_heat_capacity * 0.01, #---note: See Fig. 2b in Hilsenrath.1955: Errors do rarely exceed 1%, so I use 1% as conservative error. Note that I specify errors here only because this makes modeling much easier with brms
      N = 0,
      C = 0,
      ash_content = 0,
      sample_label = "air"
    )
  
}



#' Data from Oconnor.2022 (received via email since supplement is not complete)
#' 
#' @export
irp_make_d_oconnor2022 <- function() {
  
  read.csv("data/raw_data/dry_thermal_conductivity/Rho-b_regressions.csv") |>
    dplyr::select(c(10, 11)) |>
    setNames(nm = c("bulk_density", "dry_thermal_conductivity")) |>
    dplyr::filter(! is.na(bulk_density))
  
}


#' Summary table of predictive accuracy from other models
#' 
#' @export
irp_make_data_evaluation_literature_values <- function(irp_data_evaluation_literature_values_file) {
  
  read.csv(irp_data_evaluation_literature_values_file) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(c("summary_statistics_lower", "summary_statistics_upper", "y_lower", "y_upper")),
        function(.x) {
          res <- 
            if(dplyr::cur_column() %in% c("summary_statistics_lower", "summary_statistics_upper")) {
              purrr::map2_dbl(.x, summary_statistics_function_to_unit, function(.x1, .y) {
                switch(
                  .y,
                  "expm1" = expm1(.x1),
                  .x1
                )
              })
            } else {
              .x
            } 
          
          res |>
            magrittr::divide_by(summary_statistics_factor_to_unit) |>
            units::mixed_units(summary_statistics_unit, mode = "standard") |>
            units::set_units(summary_statistics_unit_target, mode = "standard")
        }
      )
    )
  
}


#' Peat hydrophysical data from Whittington (2024)
#' 
#' @export
irp_make_data_wittington2024 <- function() {
  
  # extract data
  d83_1 <- 
    dplyr::bind_rows(
      readRDS("data/raw_data/caldat-wittington2024/Whittington.2024-Fig1a")$processed_data
    ) |>
    dplyr::select(id, x, y) |>
    dplyr::rename(
      sample_depth = "y",
      bulk_density = "x",
      core_label = "id"
    ) |>
    dplyr::mutate(
      spectra =
        purrr::map(seq_along(sample_depth), function(i) {
          tibble::tibble(x = numeric(), y = numeric())
        })
    ) |>
    ir::ir_as_ir() %>%
    irpeat::irp_saturated_hydraulic_conductivity_1(
      do_summary = FALSE, 
      bulk_density = .$bulk_density, 
      summary_function_mean = median
    ) %>%
    irpeat::irp_volume_fraction_solids_1(
      do_summary = FALSE, 
      bulk_density = .$bulk_density
    ) |>
    dplyr::mutate(
      total_porosity_1 = irpeat:::irp_set_units_rvar(posterior::rvar(1), L/L) - volume_fraction_solids_1,
      total_porosity_1_mean = 
        total_porosity_1 |>
        mean() |>
        units::set_units("L/L", mode = "standard"),
      total_porosity_1_lower = 
        total_porosity_1 |>
        posterior::quantile2(probs = c(0.025)) |>
        units::set_units("L/L", mode = "standard"),
      total_porosity_1_upper = 
        total_porosity_1 |>
        posterior::quantile2(probs = c(0.975)) |>
        units::set_units("L/L", mode = "standard"),
      saturated_hydraulic_conductivity_1 = 
        saturated_hydraulic_conductivity_1 |>
        irpeat:::irp_set_units_rvar("m/s", mode = "standard"), 
      saturated_hydraulic_conductivity_1_mean =
        saturated_hydraulic_conductivity_1 |>
        mean() |>
        units::set_units("m/s", mode = "standard"),
      saturated_hydraulic_conductivity_1_lower =
        saturated_hydraulic_conductivity_1 |>
        posterior::quantile2(probs = c(0.025)) |>
        units::set_units("m/s", mode = "standard"),
      saturated_hydraulic_conductivity_1_upper =
        saturated_hydraulic_conductivity_1 |>
        posterior::quantile2(probs = c(0.975)) |>
        units::set_units("m/s", mode = "standard")
    )
  
  d83_2 <- 
    dplyr::bind_rows(
      readRDS("data/raw_data/caldat-wittington2024/Whittington.2024-Fig1c")$processed_data
    ) |>
    dplyr::select(id, x, y) |>
    dplyr::rename(
      sample_depth = "y",
      total_porosity_1_mean = "x",
      core_label = "id"
    ) |>
    dplyr::mutate(
      total_porosity_1_mean = units::set_units(total_porosity_1_mean, L/L)
    )
  
  d83_3 <- 
    dplyr::bind_rows(
      readRDS("data/raw_data/caldat-wittington2024/Whittington.2024-Fig1d")$processed_data
    ) |>
    dplyr::select(id, x, y) |>
    dplyr::rename(
      sample_depth = "y",
      saturated_hydraulic_conductivity_1_mean = "x",
      core_label = "id"
    ) |>
    dplyr::mutate(
      core_label = stringr::str_remove(core_label, pattern = "_"),
      saturated_hydraulic_conductivity_1_mean = units::set_units(saturated_hydraulic_conductivity_1_mean, m/s)
    )
  
  # predictions from model 4 from Morris et al. 2022
  d83_4 <- 
    dplyr::bind_rows(
      list.files("data/raw_data/caldat-wittington2024/", pattern = "Fig3-core", full.names = TRUE) |>
        purrr::map(function(.x) {
          readRDS(.x)$processed_data |>
            dplyr::mutate(
              core_label = stringr::str_extract(.x, pattern = "core[A-Z]{1}")
            )
        })
    ) |>
    dplyr::select(x, y, core_label, id) |>
    dplyr::rename(
      sample_depth = "y",
      saturated_hydraulic_conductivity_1_mean = "x",
      id_model_morris2022 = "id"
    ) |>
    dplyr::mutate(
      saturated_hydraulic_conductivity_1_mean = units::set_units(saturated_hydraulic_conductivity_1_mean, m/s)
    )
  
  d83 <- 
    dplyr::bind_rows(
      d83_1 |>
        ir::ir_drop_spectra() |>
        dplyr::mutate(
          ytype = "rep"
        ),
      d83_2 |>
        dplyr::mutate(
          ytype = ""
        ),
      d83_3 |>
        dplyr::mutate(
          ytype = ""
        ),
      d83_4 |>
        dplyr::mutate(
          ytype = paste0("morris2022_", id_model_morris2022)
        )
    )
  
  d83 |>
    dplyr::select(-saturated_hydraulic_conductivity_1, -total_porosity_1)
  
}


