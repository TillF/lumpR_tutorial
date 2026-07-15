# Tutorial workflow for using lumpR for generating geodata input files for the hydrological model WASA-SED (alternative dataset with snow and reservoir parameters)
# This workflow was provided by Alban Doko. A more verbose version can be found in isabena_tutorial.R

# =============================================================================
# Recommended use: run section by section.


# =============================================================================
# 0. PACKAGES
# =============================================================================

required_cran <- c("DBI", "odbc", "rgrass", "devtools")

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

if (!requireNamespace("lumpR", quietly = TRUE)) {
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)
  devtools::install_github(
    "tpilz/lumpR",
    ref = "wo_rgeos_maptools"
  )
}

library(lumpR)
library(rgrass)
library(DBI)
library(odbc)


# SETTINGS ####

# switch to specified working directory (this is usually the home directory of this very script)
if (!require("rstudioapi", quietly = TRUE)) install.packages("rstudioapi")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#setwd("/home/myname/somefolder") #manual specification. Use "/" instead of "\" in Windows

## GRASS-GIS-SETTINGS ####

# GRASS engine
addon_path <- "/home/alban/.grass8/addons/bin"  # normally not required: path to your locally installed GRASS add-ons, e.g. "/home/myuser/.grass8/addons/bin". Must only be given if necessary, see ?lump_grass_prep
gisBase <- "/usr/lib/grass83" #path to GRASS installation, e.g. "/usr/lib/grass83" (use / instead of / under windows, e.g. "d:/programme/GRASS7.0" )
gisDbase <- "/home/alban/grassdata" # path to 'grassdata' directory containing the location specified above and all corresp. data. 
grass_location <- "seman" # GRASS location, i.e. GRASS structure were data are stored
grass_mapset <- "alban" # corresponding GRASS mapset

initGRASS(gisBase=gisBase, 
          home=getwd(), # The directory in which to create the .gisrc file
          location=location,
          mapset=mapset,
          gisDbase=gisDbase,  
          override=TRUE)

#test if GRASS connections working
lumpR:::test_grass()

# GRASS maps (input)
vname <- "St_Miras" #name of GRASS point vector layer containing subbasin coordinates
outlet_row <- 1 #row of outlet point of entire watershed

dem <- "dem" # DEM raster - MANDATORY
lcov <- "vegetation" "vegetation" # land / vegetation cover raster map in GRASS location (e.g. derived from CORINE) - MANDATORY
soil <- "soils" # soil raster map in GRASS location - MANDATORY (The package [SoilDataPrep](https://github.com/TillF/SoilDataPrep/tree/master) is designed to help you with the generation of soil data, especially in conjunction with lumpR) 

soil_depth = NULL # soil depth raster map
watermask = NULL # water mask raster map in GRASS location (1=water, 0=no water)
imperviousmask = NULL # impervious surface areas raster map in GRASS location (1=impervious, 0=permeable)
river = NULL # river vector map

# GRASS maps (output)
subbas <- "subbasin" # subbasin raster map generated - MANDATORY
stream_pref <- "stream_accum" # prefix of calculated stream segments raster and vector maps
drainp_processed <- "drain_points" # prefix of drainage point vector files (given points snapped to river and internally calculated points, respectively) - MANDATORY
eha <- "eha" # elementary hillslope areas raster map - MANDATORY
flowdir <- "flowdir" # flow direction raster map - MANDATORY
flowacc <- "flowacc" # flow accumulation raster map - MANDATORY
stream <- "stream" # stream segments raster map (based on eha calculation; much finer than 'stream_pref' to delineate hillslopes) - MANDATORY
stream_horton <- "stream_horton" # Horton stream order raster map (based on 'stream' above) - MANDATORY
elevriv <- "elevriv" # elevation relative to next river cell raster map - MANDATORY
distriv <- "distriv" # distance to next river cell raster map - MANDATORY
svc <- "svc" # soil vegetation components raster map - MANDATORY
lu <- "lu" # landscape units raster map - MANDATORY

# soil and vegetation properties
# path to prepared vegetation parameter table 'vegetation.dat' in WASA format - MANDATORY
#'vegetation.dat' can be e.g. be generated from CORINE data and the mapping table provided in 2bprepared/vegetation_from_corine
veg_path <- file.path(dat_dir, "2bprepared", "from_AlbansScripts")
soil_path <- file.path(dat_dir, "2bprepared", "from_SoilDataPrep")


# Text output names
svc_file <- "soil_vegetation_components.dat" # name of file containing svc parameters - MANDATORY
catena_out <- "rstats.txt"
catena_head_out <- "rstats_head.txt"
sub_ofile <- "sub_stats.txt" # Name of subbasin statistics file containing subbasin parameters
lu_ofile <- "lu_stats.txt" # Name of file containing subbasins and the corresponding LUs with their fraction of area in the subbasin
lupar_ofile <- "lu_pars.txt" # Name of file containing LUs and related parameters

snow_eha_file <- "snow_eha.txt" #optional file containing additional properties at the EHA-scale. Can be generated e.g. with prepare_snow_input()
snow_eha_head_file <- "snow_eha_head.txt" #corresponding header file

## lumpR- PARAMETERS ####
# STRONGLY CASE-STUDY SPECIFIC! Try out what suits your needs / data

### MISCELLANEOUS PARAMETERS ####
# parameters influencing some outputs but not directly discretisation complexity

# parameters to calc_subbas. For details, see ?calc_subbas
thresh_stream <- 1000   # Threshold for derivation of stream from flow accumulation map (in cells). Needs to be set only if river is not set - MANDATORY
snap_dist <- 500 # maximum distance for snapping of drain_points to stream cells in the units of your GRASS location (usually meters) - MANDATORY
rm_spurious <- 0.01 # threshold for small spurious subbasins created within calc_subbas() to be removed.

# parameters to lump_grass_prep.  For details, see ?lump_grass_prep  
sizefilter <- 40   # minimum size of EHAs (in map units, usually m²) to be preserved, smaller EHAs (artefacts) are removed; parameter for GRASS function r.reclass.area - MANDATORY
growrad <- 50 # growing radius (in raster cells) to remove artefacts in EHA data; parameter for GRASS function r.grow, see see ?lump_grass_prep - MANDATORY

# parameters to area2catena  For details, see ?area2catena
min_cell_in_slope <- 30 # minimum number of cells a hillslope area must have, all smaller ones are skipped
min_catena_length <- 3 # minimum number of sampling points (cells) a catena should have. If there are less, the catena is not saved
max_riv_dist <- 15 # maximum distance to river [in cells]: if the closest cell of an EHA is farther than max_riv_dist, the EHA is skipped, otherwise all distances within the EHA are redurced by the distance of the closest cell to river

### LANDSCAPE DISCRETISATION PARAMETERS ####
thresh_sub <- NULL # Parameter for GRASS function r.watershed defining the minimum size of an exterior watershed basin in number of grid cells. If NULL only the given drainage points are used for subbasin delineation
eha_thres <- 300  # parameter for GRASS function r.watershed. This is a crucial parameter affecting the size of delineated hillslopes - MANDATORY

# Landscape classification
supp_quant <- character(0) 
supp_quant_no_classes <- numeric(0) # named vector with names of GRASS raster maps containing *quantitative* attributes for LU deviation and orresponding number of classes, that should be created from the quantitative attributes

supp_qual <- c("svc", "soils")
supp_qual_no_classes <- c(1, 3)

# shape, extent, elevation weighting, supplemental fields, slope width
no_classes <- c(
  -3,
  3,
  10,
  supp_quant_no_classes,
  supp_qual_no_classes,
  1
)

no_TCs <- 3

# Runtime
keep_temp <- FALSE
overwrite <- TRUE
silent <- FALSE
plot_catena <- TRUE
plot_profclass <- TRUE
grass_files <- TRUE
ncores <- 4

# Database
dsn <- "MyMariaDBDataSource"
dbname <- dsn

driver <- "MariaDB"
server <- "localhost"
port <- 3306
database <- "mydatabase"
user <- "myuser"
password <- "mypassword"
odbc_file <- path.expand("~/.odbc.ini")

# Database cleaning
apply_small_area_filter <- TRUE

# 0.001 = 0.1%. This is less aggressive than 0.01 = 1%.
small_area_threshold <- 0.001

# Since watermask and imperviousmask are NULL, these are FALSE automatically.
apply_remove_water_svc <- !is.null(watermask)
apply_remove_impervious_svc <- !is.null(imperviousmask)

apply_compute_rocky_frac <- TRUE

# Reservoir processing
#
# "none":
#   No reservoir parameter files are imported. Reservoir database tables remain
#   empty and db_wasa_input() writes the standard empty reservoir outputs.
#
# "reuse_existing":
#   Do not run reservoir_outlet(), reservoir_strategic(), or
#   reservoir_lumped(). Import reservoir parameter files already in work_dir.
#
# "generate":
#   Regenerate reservoir files from the GRASS vectors defined below.
reservoir_mode <- "none"

strategic_reservoir_vector <- "reservoir_vec"
strategic_outlet_vector <- "res_outlets"
small_reservoir_vector <- "res_small"
small_reservoir_class_vector <- "res_small_2"

# Neutral LU parameters used by the original tutorial.
# Replace these later only when calibration/field information supports it.
lu_soil_depth <- -1
lu_allu_depth <- -1
lu_riverbed_depth <- 2000
lu_kf_bedrock <- -9999
lu_gw_dist <- -9999
lu_frgw_delay <- -9999
lu_sdr <- 1

# -----------------------------------------------------------------------------
# MODEL MODE AND OPTIONAL PARAMETER INPUTS
# -----------------------------------------------------------------------------

# FALSE:
#   Prepare a hydrology + snow model. MUSLE-C and coarse_fraction may remain
#   undefined because sediment simulation is not enabled.
#
# TRUE:
#   Prepare a sediment model. The script requires valid MUSLE-C values and
#   checks the exported sediment parameters.
enable_sediment <- FALSE

# Optional MUSLE-C values by vegetation ID.
# Fill these only when enable_sediment <- TRUE.
# The same value is written to all four seasonal MUSLE-C columns unless you
# replace this vector with a more detailed seasonal table.
musle_c_by_vegetation <- c(
  `104` = NA_real_,
  `105` = NA_real_,
  `106` = NA_real_,
  `107` = NA_real_,
  `108` = NA_real_,
  `109` = NA_real_,
  `110` = NA_real_,
  `112` = NA_real_,
  `113` = NA_real_
)

# Optional prepared rainy_season.dat.
#
# NULL:
#   Keep the empty file generated by db_wasa_input() and print a warning.
#
# Character path:
#   Copy and validate that prepared file after the WASA export.
rainy_season_source <- NULL

# Set TRUE when a non-empty rainy_season.dat is mandatory for the intended run.
require_rainy_season <- FALSE


# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================

assert_files_exist <- function(paths) {
  missing <- paths[!file.exists(paths)]

  if (length(missing) > 0) {
    stop(
      "Required file(s) not found:\n",
      paste(missing, collapse = "\n")
    )
  }
}


configure_odbc_dsn <- function() {
  dsn_lines <- c(
    paste0("[", dsn, "]"),
    "Description = lumpR analysis database",
    paste0("Driver = ", driver),
    paste0("Server = ", server),
    paste0("Port = ", port),
    paste0("Database = ", database),
    paste0("User = ", user),
    paste0("Password = ", password)
  )

  existing <- if (file.exists(odbc_file)) {
    readLines(odbc_file, warn = FALSE)
  } else {
    character(0)
  }

  start <- grep(
    paste0("^\\[", dsn, "\\]$"),
    existing
  )

  if (length(start) > 0) {
    section_starts <- grep("^\\[.*\\]$", existing)
    next_section <- section_starts[section_starts > start[1]]

    end <- if (length(next_section) > 0) {
      next_section[1] - 1
    } else {
      length(existing)
    }

    existing <- existing[-seq.int(start[1], end)]
  }

  writeLines(
    c(existing, "", dsn_lines),
    odbc_file
  )

  message("ODBC DSN updated: ", dsn)
}


reset_database_tables <- function() {
  con <- DBI::dbConnect(
    odbc::odbc(),
    dsn = dsn,
    uid = user,
    pwd = password
  )

  on.exit(DBI::dbDisconnect(con), add = TRUE)

  DBI::dbExecute(
    con,
    "SET FOREIGN_KEY_CHECKS = 0;"
  )

  tabs <- DBI::dbListTables(con)

  for (tbl in tabs) {
    sql <- paste0(
      "DROP TABLE IF EXISTS `",
      gsub("`", "``", tbl, fixed = TRUE),
      "`;"
    )

    DBI::dbExecute(con, sql)
  }

  DBI::dbExecute(
    con,
    "SET FOREIGN_KEY_CHECKS = 1;"
  )

  message(
    "Database reset complete. Tables removed: ",
    length(tabs)
  )
}


read_table_checked <- function(
    path,
    header = TRUE,
    sep = "",
    skip = 0
) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }

  read.table(
    path,
    header = header,
    sep = sep,
    skip = skip,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    fill = TRUE
  )
}


check_vegetation_overlap <- function() {
  execGRASS(
    "g.region",
    raster = subbas,
    flags = "a"
  )

  try(
    execGRASS("r.mask", flags = "r"),
    silent = TRUE
  )

  on.exit(
    try(
      execGRASS("r.mask", flags = "r"),
      silent = TRUE
    ),
    add = TRUE
  )

  execGRASS(
    "r.mask",
    raster = subbas
  )

  raw <- execGRASS(
    "r.stats",
    input = paste(subbas, lcov, sep = ","),
    flags = c("c", "n"),
    separator = "comma",
    intern = TRUE
  )

  if (length(raw) == 0) {
    stop("No vegetation/subbasin overlap was returned by r.stats.")
  }

  parts <- do.call(
    rbind,
    strsplit(raw, ",", fixed = TRUE)
  )

  result <- data.frame(
    subbasin_id = as.integer(parts[, 1]),
    vegetation_id = as.integer(parts[, 2]),
    number_cells = as.numeric(parts[, 3])
  )

  result <- result[
    order(
      result$subbasin_id,
      result$vegetation_id
    ),
  ]

  message(
    "Vegetation IDs inside subbasins: ",
    paste(
      sort(unique(result$vegetation_id)),
      collapse = ", "
    )
  )

  print(result)

  invisible(result)
}


prepare_snow_header <- function(path) {
  writeLines(
    c(
      "eha_id\taspect\trel_alt",
      "1\t2\t1",
      "0\t4\t3"
    ),
    path
  )
}


repair_tc_contains_svc <- function(
    tc_file = "tc.dat",
    svc_definition_file = "soil_vegetation_components.dat",
    reclass_file = "reclass_svc.txt",
    output_file = "tc_contains_svc.dat",
    tolerance = 0.02
) {
  tc <- read_table_checked(
    tc_file,
    header = TRUE,
    sep = ""
  )

  svc_def <- read_table_checked(
    svc_definition_file,
    header = TRUE,
    sep = ""
  )

  required_svc_fields <- c("pid", "veg_id")
  missing_svc_fields <- setdiff(
    required_svc_fields,
    names(svc_def)
  )

  if (length(missing_svc_fields) > 0) {
    stop(
      "Missing field(s) in ",
      svc_definition_file,
      ": ",
      paste(missing_svc_fields, collapse = ", ")
    )
  }

  tc_id_col <- if ("TC" %in% names(tc)) {
    "TC"
  } else if ("tc_id" %in% names(tc)) {
    "tc_id"
  } else {
    names(tc)[1]
  }

  svc_columns <- grep(
    "^svc_c[0-9]+$",
    names(tc),
    value = TRUE
  )

  if (length(svc_columns) == 0) {
    stop(
      "No columns named svc_c1, svc_c2, ... were found in tc.dat.\n",
      "Available columns:\n",
      paste(names(tc), collapse = ", ")
    )
  }

  modified_svc_ids <- as.integer(
    sub("^svc_c", "", svc_columns)
  )

  final_svc_ids <- modified_svc_ids

  if (file.exists(reclass_file)) {
    reclass <- read_table_checked(
      reclass_file,
      header = TRUE,
      sep = ""
    )

    if (
      all(
        c("new_id", "original_id") %in%
          names(reclass)
      )
    ) {
      id_map <- setNames(
        reclass$original_id,
        as.character(reclass$new_id)
      )

      mapped <- as.integer(
        id_map[as.character(modified_svc_ids)]
      )

      if (
        !anyNA(mapped) &&
        all(mapped %in% svc_def$pid)
      ) {
        final_svc_ids <- mapped
      }
    }
  }

  unknown_svc <- setdiff(
    final_svc_ids,
    svc_def$pid
  )

  if (length(unknown_svc) > 0) {
    stop(
      "SVC IDs derived from tc.dat are absent from ",
      svc_definition_file,
      ": ",
      paste(unknown_svc, collapse = ", ")
    )
  }

  fraction_matrix <- as.matrix(
    tc[, svc_columns, drop = FALSE]
  )

  storage.mode(fraction_matrix) <- "numeric"

  result_list <- vector(
    mode = "list",
    length = nrow(tc)
  )

  for (i in seq_len(nrow(tc))) {
    values <- fraction_matrix[i, ]

    keep <- is.finite(values) & values > 0

    if (!any(keep)) {
      stop(
        "Terrain component ",
        tc[[tc_id_col]][i],
        " has no positive SVC fraction in tc.dat."
      )
    }

    result_list[[i]] <- data.frame(
      tc_id = as.integer(tc[[tc_id_col]][i]),
      svc_id = final_svc_ids[keep],
      fraction = as.numeric(values[keep])
    )
  }

  result <- do.call(
    rbind,
    result_list
  )

  # Combine duplicate TC-SVC pairs, if any.
  result <- aggregate(
    fraction ~ tc_id + svc_id,
    data = result,
    FUN = sum
  )

  fraction_sums <- aggregate(
    fraction ~ tc_id,
    data = result,
    FUN = sum
  )

  bad_tc <- fraction_sums[
    abs(fraction_sums$fraction - 1) > tolerance,
  ]

  if (nrow(bad_tc) > 0) {
    stop(
      "SVC fractions in tc.dat are not close to 1 for terrain component(s): ",
      paste(bad_tc$tc_id, collapse = ", ")
    )
  }

  # Normalize small rounding differences to exactly 1.
  sum_map <- setNames(
    fraction_sums$fraction,
    fraction_sums$tc_id
  )

  result$fraction <- result$fraction /
    sum_map[as.character(result$tc_id)]

  result <- result[
    order(result$tc_id, result$svc_id),
  ]

  write.table(
    result,
    output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  linked <- merge(
    result,
    svc_def[, c("pid", "veg_id")],
    by.x = "svc_id",
    by.y = "pid",
    all.x = TRUE
  )

  represented_veg <- sort(
    unique(linked$veg_id)
  )

  expected_veg <- sort(
    unique(svc_def$veg_id)
  )

  missing_veg <- setdiff(
    expected_veg,
    represented_veg
  )

  message(
    "Corrected TC-SVC vegetation IDs: ",
    paste(represented_veg, collapse = ", ")
  )

  if (length(missing_veg) > 0) {
    stop(
      "Vegetation IDs are still missing after rebuilding ",
      output_file,
      ": ",
      paste(missing_veg, collapse = ", ")
    )
  }

  message(
    output_file,
    " rebuilt successfully with ",
    nrow(result),
    " TC-SVC links."
  )

  invisible(result)
}


prepare_lu_database_file <- function(
    lu_file = "lu.dat",
    lu_parameters_file = lupar_ofile,
    output_file = "lu_db.dat"
) {
  luout <- read_table_checked(
    lu_file,
    header = TRUE,
    sep = ""
  )

  lupar <- read_table_checked(
    lu_parameters_file,
    header = TRUE,
    sep = ""
  )

  required_luout <- c(
    "LU-ID",
    "x_length",
    "aspect_p1_c1",
    "aspect_p1_c2",
    "rel_alt_p1"
  )

  missing <- setdiff(
    required_luout,
    names(luout)
  )

  if (length(missing) > 0) {
    stop(
      "Missing required column(s) in lu.dat: ",
      paste(missing, collapse = ", ")
    )
  }

  if (!"pid" %in% names(lupar)) {
    stop(
      "Column 'pid' is missing from ",
      lu_parameters_file
    )
  }

  idx <- match(
    as.integer(lupar$pid),
    as.integer(luout[["LU-ID"]])
  )

  if (anyNA(idx)) {
    stop(
      "LU IDs in ",
      lu_parameters_file,
      " missing from lu.dat: ",
      paste(lupar$pid[is.na(idx)], collapse = ", ")
    )
  }

  if (anyDuplicated(lupar$pid)) {
    stop("Duplicate LU IDs detected in ", lu_parameters_file)
  }

  lupar$slopelength <- luout[["x_length"]][idx]

  aspect_sin <- luout[["aspect_p1_c1"]][idx]
  aspect_cos <- luout[["aspect_p1_c2"]][idx]

  lupar$aspect <- round(
    (
      atan2(aspect_sin, aspect_cos) *
        180 / pi
    ) %% 360,
    2
  )

  lupar$relative_altitude <- round(
    luout[["rel_alt_p1"]][idx],
    2
  )

  lupar$soil_depth <- lu_soil_depth
  lupar$allu_depth <- lu_allu_depth
  lupar$riverbed_depth <- lu_riverbed_depth
  lupar$kf_bedrock <- lu_kf_bedrock
  lupar$gw_dist <- lu_gw_dist
  lupar$frgw_delay <- lu_frgw_delay

  if ("description" %in% names(lupar)) {
    missing_description <- is.na(lupar$description) |
      trimws(lupar$description) == ""

    lupar$description[missing_description] <- paste0(
      "LU_",
      lupar$pid[missing_description]
    )
  }

  if ("sdr_lu" %in% names(lupar)) {
    lupar$sdr_lu[
      is.na(lupar$sdr_lu)
    ] <- lu_sdr
  }

  if (anyNA(lupar$aspect)) {
    stop("NA aspect values detected in lu_db.dat preparation.")
  }

  if (anyNA(lupar$relative_altitude)) {
    stop(
      "NA relative_altitude values detected in ",
      "lu_db.dat preparation."
    )
  }

  if (
    any(
      lupar$aspect < 0 |
      lupar$aspect >= 360
    )
  ) {
    stop("Aspect values must be in the range [0, 360).")
  }

  if (nrow(lupar) != nrow(luout)) {
    stop(
      "LU count mismatch: ",
      nrow(lupar),
      " rows in ",
      lu_parameters_file,
      " and ",
      nrow(luout),
      " rows in lu.dat."
    )
  }

  write.table(
    lupar,
    output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )

  message(
    output_file,
    " created with ",
    nrow(lupar),
    " LUs and snow parameters."
  )

  invisible(lupar)
}


database_status <- function(stage) {
  con <- DBI::dbConnect(
    odbc::odbc(),
    dsn = dsn,
    uid = user,
    pwd = password
  )

  on.exit(DBI::dbDisconnect(con), add = TRUE)

  lu_ids <- DBI::dbGetQuery(
    con,
    "SELECT pid FROM landscape_units ORDER BY pid"
  )$pid

  veg_ids <- DBI::dbGetQuery(
    con,
    "SELECT pid FROM vegetation ORDER BY pid"
  )$pid

  active_veg <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT DISTINCT s.veg_id",
      "FROM soil_veg_components s",
      "INNER JOIN r_tc_contains_svc r",
      "ON r.svc_id = s.pid",
      "ORDER BY s.veg_id"
    )
  )$veg_id

  message("\nDATABASE STATUS: ", stage)
  message("Active LUs: ", length(lu_ids))
  message(
    "Vegetation table IDs: ",
    paste(veg_ids, collapse = ", ")
  )
  message(
    "Vegetation linked to TCs: ",
    paste(active_veg, collapse = ", ")
  )

  invisible(
    list(
      lu_ids = lu_ids,
      vegetation_ids = veg_ids,
      active_vegetation_ids = active_veg
    )
  )
}


extract_first_integer_from_data_lines <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }

  lines <- trimws(
    readLines(path, warn = FALSE)
  )

  data_lines <- lines[
    grepl("^[0-9]+([[:space:]]|$)", lines)
  ]

  if (length(data_lines) == 0) {
    return(integer(0))
  }

  as.integer(
    vapply(
      strsplit(data_lines, "[[:space:]]+"),
      FUN = `[`,
      FUN.VALUE = character(1),
      1
    )
  )
}


extract_hymo_lu_ids <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }

  lines <- trimws(
    readLines(path, warn = FALSE)
  )

  data_lines <- lines[
    grepl("^[0-9]+([[:space:]]|$)", lines)
  ]

  lu_ids <- integer(0)

  for (line in data_lines) {
    values <- strsplit(
      line,
      "[[:space:]]+"
    )[[1]]

    if (length(values) < 4) {
      next
    }

    number_lu <- suppressWarnings(
      as.integer(values[3])
    )

    if (
      is.na(number_lu) ||
      number_lu < 1 ||
      length(values) < 3 + number_lu
    ) {
      stop(
        "Invalid hymo.dat data line: ",
        line
      )
    }

    current_ids <- as.integer(
      values[
        4:(3 + number_lu)
      ]
    )

    lu_ids <- c(
      lu_ids,
      current_ids
    )
  }

  sort(unique(lu_ids))
}


ensure_lu2_from_clean_database <- function() {
  hillslope_dir <- file.path(
    dat_dir,
    "WASA_input",
    "Hillslope"
  )

  dir.create(
    hillslope_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  lu2_file <- file.path(
    hillslope_dir,
    "lu2.dat"
  )

  con <- DBI::dbConnect(
    odbc::odbc(),
    dsn = dsn,
    uid = user,
    pwd = password
  )

  on.exit(DBI::dbDisconnect(con), add = TRUE)

  lu2_db <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT pid, aspect, relative_altitude",
      "FROM landscape_units",
      "ORDER BY pid"
    )
  )

  if (
    nrow(lu2_db) == 0 ||
    anyNA(lu2_db$aspect) ||
    anyNA(lu2_db$relative_altitude)
  ) {
    stop(
      "The cleaned landscape_units table does not contain ",
      "complete snow parameters."
    )
  }

  rewrite <- !file.exists(lu2_file)

  if (!rewrite) {
    current_ids <- extract_first_integer_from_data_lines(
      lu2_file
    )

    rewrite <- !setequal(
      current_ids,
      lu2_db$pid
    )
  }

  if (rewrite) {
    writeLines(
      c(
        "Specification of aspect and relative elevation of landscape units",
        paste(
          "LU-ID[id]",
          "aspect[deg]",
          "mean_altitude_over_subbas_mean[m]",
          sep = "\t"
        )
      ),
      lu2_file
    )

    lu2_out <- data.frame(
      pid = as.integer(lu2_db$pid),
      aspect = round(lu2_db$aspect),
      relative_altitude = round(
        lu2_db$relative_altitude
      )
    )

    write.table(
      lu2_out,
      lu2_file,
      append = TRUE,
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )

    message(
      "lu2.dat regenerated from the cleaned database with ",
      nrow(lu2_out),
      " active LUs."
    )
  }

  invisible(lu2_db)
}



fix_soter_header <- function() {
  soter_file <- file.path(
    dat_dir,
    "WASA_input",
    "Hillslope",
    "soter.dat"
  )

  assert_files_exist(soter_file)

  lines <- readLines(
    soter_file,
    warn = FALSE
  )

  if (length(lines) < 3) {
    stop("soter.dat is incomplete.")
  }

  data_idx <- which(
    grepl(
      "^[[:space:]]*[0-9]+([[:space:]]|$)",
      lines
    )
  )

  if (length(data_idx) == 0) {
    stop("No LU data rows were found in soter.dat.")
  }

  split_fields <- function(x) {
    strsplit(
      trimws(x),
      "[[:space:]]+"
    )[[1]]
  }

  header_fields <- split_fields(lines[2])
  first_data_fields <- split_fields(lines[data_idx[1]])

  if (
    length(first_data_fields) ==
      length(header_fields) + 1
  ) {
    lines[2] <- paste(
      c(
        header_fields,
        "sdr_lu[-]"
      ),
      collapse = "\t"
    )

    writeLines(
      lines,
      soter_file
    )

    message(
      "Added missing sdr_lu[-] label to the soter.dat header."
    )

  } else if (
    length(first_data_fields) !=
      length(header_fields)
  ) {
    stop(
      "Unexpected soter.dat structure: header has ",
      length(header_fields),
      " fields but the first data row has ",
      length(first_data_fields),
      "."
    )
  }

  invisible(soter_file)
}


apply_musle_c_to_vegetation_file <- function(
    vegetation_file = file.path(
      dat_dir,
      "vegetation.txt"
    )
) {
  if (!enable_sediment) {
    message(
      "Hydrology/snow mode: no MUSLE-C values are required."
    )
    return(invisible(NULL))
  }

  veg <- read_table_checked(
    vegetation_file,
    header = TRUE,
    sep = "\t"
  )

  required_columns <- paste0(
    "c_musle_c",
    1:4
  )

  missing_columns <- setdiff(
    required_columns,
    names(veg)
  )

  if (length(missing_columns) > 0) {
    stop(
      "Sediment mode requires these vegetation columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  vegetation_ids <- as.character(
    veg$pid
  )

  missing_mapping <- setdiff(
    vegetation_ids,
    names(musle_c_by_vegetation)
  )

  if (length(missing_mapping) > 0) {
    stop(
      "MUSLE-C mapping is missing vegetation ID(s): ",
      paste(missing_mapping, collapse = ", ")
    )
  }

  mapped_values <- as.numeric(
    musle_c_by_vegetation[vegetation_ids]
  )

  if (
    anyNA(mapped_values) ||
    any(!is.finite(mapped_values))
  ) {
    bad_ids <- vegetation_ids[
      is.na(mapped_values) |
        !is.finite(mapped_values)
    ]

    stop(
      "Sediment mode is enabled, but valid MUSLE-C values ",
      "were not supplied for vegetation ID(s): ",
      paste(bad_ids, collapse = ", ")
    )
  }

  if (
    any(
      mapped_values < 0 |
        mapped_values > 1
    )
  ) {
    stop(
      "MUSLE-C values must be between 0 and 1."
    )
  }

  for (column in required_columns) {
    veg[[column]] <- mapped_values
  }

  write.table(
    veg,
    vegetation_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )

  message(
    "MUSLE-C values written to vegetation.txt for ",
    length(vegetation_ids),
    " vegetation classes."
  )

  invisible(veg)
}


handle_rainy_season_file <- function() {
  output_file <- file.path(
    dat_dir,
    "WASA_input",
    "Hillslope",
    "rainy_season.dat"
  )

  if (!is.null(rainy_season_source)) {
    source_file <- path.expand(
      rainy_season_source
    )

    assert_files_exist(source_file)

    copied <- file.copy(
      source_file,
      output_file,
      overwrite = TRUE
    )

    if (!copied) {
      stop(
        "Could not copy rainy-season file from: ",
        source_file
      )
    }

    message(
      "Copied prepared rainy_season.dat from: ",
      source_file
    )
  }

  assert_files_exist(output_file)

  lines <- trimws(
    readLines(
      output_file,
      warn = FALSE
    )
  )

  data_lines <- lines[
    grepl(
      "^[0-9]+[[:space:]]+[0-9]+[[:space:]]+[0-9]+",
      lines
    )
  ]

  if (length(data_lines) == 0) {
    message_text <- paste0(
      "rainy_season.dat contains no data records. ",
      "Seasonal vegetation interpolation is not fully parameterised."
    )

    if (require_rainy_season) {
      stop(message_text)
    } else {
      warning(
        message_text,
        call. = FALSE
      )
    }

    return(
      invisible(
        list(
          path = output_file,
          records = 0L
        )
      )
    )
  }

  values <- strsplit(
    data_lines,
    "[[:space:]]+"
  )

  field_counts <- lengths(values)

  if (any(field_counts != 7)) {
    stop(
      "Every rainy_season.dat data row must contain exactly ",
      "7 fields: Subbasin, Veg_id, Year, DOY1, DOY2, DOY3, DOY4."
    )
  }

  rainy <- as.data.frame(
    do.call(
      rbind,
      values
    ),
    stringsAsFactors = FALSE
  )

  names(rainy) <- c(
    "subbasin",
    "veg_id",
    "year",
    "doy1",
    "doy2",
    "doy3",
    "doy4"
  )

  rainy[] <- lapply(
    rainy,
    as.integer
  )

  if (
    anyNA(rainy) ||
    any(
      rainy[, c(
        "doy1",
        "doy2",
        "doy3",
        "doy4"
      )] < 1
    ) ||
    any(
      rainy[, c(
        "doy1",
        "doy2",
        "doy3",
        "doy4"
      )] > 366
    )
  ) {
    stop(
      "rainy_season.dat contains missing or invalid day-of-year values."
    )
  }

  if (
    any(
      !(
        rainy$doy1 <= rainy$doy2 &
        rainy$doy2 <= rainy$doy3 &
        rainy$doy3 <= rainy$doy4
      )
    )
  ) {
    stop(
      "rainy_season.dat requires DOY1 <= DOY2 <= DOY3 <= DOY4."
    )
  }

  message(
    "rainy_season.dat validated with ",
    nrow(rainy),
    " records."
  )

  invisible(
    list(
      path = output_file,
      records = nrow(rainy),
      data = rainy
    )
  )
}


validate_sediment_outputs <- function() {
  svc_file_out <- file.path(
    dat_dir,
    "WASA_input",
    "Hillslope",
    "svc.dat"
  )

  assert_files_exist(svc_file_out)

  # svc.dat is tab-delimited and its MUSLE-K header contains spaces.
  # Read it explicitly as a tab-delimited table.
  svc_out <- read.delim(
    svc_file_out,
    header = TRUE,
    skip = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )

  # Column names include units, so identify the required variables by prefix.
  musle_c_col <- grep(
    "^musle_c",
    names(svc_out),
    value = TRUE
  )[1]

  coarse_col <- grep(
    "^coarse_fraction",
    names(svc_out),
    value = TRUE
  )[1]

  if (
    is.na(musle_c_col) ||
    is.na(coarse_col)
  ) {
    stop(
      "Could not identify MUSLE-C or coarse_fraction in exported svc.dat."
    )
  }

  invalid_c <- (
    !is.finite(svc_out[[musle_c_col]]) |
    svc_out[[musle_c_col]] < 0 |
    svc_out[[musle_c_col]] > 1 |
    svc_out[[musle_c_col]] == -9999
  )

  invalid_coarse <- (
    !is.finite(svc_out[[coarse_col]]) |
    svc_out[[coarse_col]] < 0 |
    svc_out[[coarse_col]] > 100 |
    svc_out[[coarse_col]] == -9999
  )

  if (enable_sediment) {
    if (any(invalid_c)) {
      stop(
        "Sediment mode: exported svc.dat contains invalid MUSLE-C values."
      )
    }

    if (any(invalid_coarse)) {
      stop(
        "Sediment mode: exported svc.dat contains invalid ",
        "coarse_fraction values. Correct the soil/horizon parameterisation ",
        "before sediment simulation."
      )
    }

    message(
      "Sediment parameters passed final validation."
    )

  } else {
    if (any(invalid_c)) {
      message(
        "Hydrology/snow mode: MUSLE-C remains undefined; ",
        "this is acceptable only while sediment simulation is disabled."
      )
    }

    if (any(invalid_coarse)) {
      message(
        "Hydrology/snow mode: coarse_fraction remains undefined; ",
        "this is acceptable only while sediment simulation is disabled."
      )
    }
  }

  invisible(svc_out)
}


validate_final_outputs <- function() {
  hillslope_dir <- file.path(
    dat_dir,
    "WASA_input",
    "Hillslope"
  )

  soter_file <- file.path(
    hillslope_dir,
    "soter.dat"
  )

  lu2_file <- file.path(
    hillslope_dir,
    "lu2.dat"
  )

  hymo_file <- file.path(
    hillslope_dir,
    "hymo.dat"
  )

  vegetation_file <- file.path(
    hillslope_dir,
    "vegetation.dat"
  )

  assert_files_exist(
    c(
      soter_file,
      lu2_file,
      hymo_file,
      vegetation_file
    )
  )

  soter_ids <- extract_first_integer_from_data_lines(
    soter_file
  )

  lu2_ids <- extract_first_integer_from_data_lines(
    lu2_file
  )

  hymo_ids <- extract_hymo_lu_ids(
    hymo_file
  )

  if (!setequal(soter_ids, lu2_ids)) {
    stop(
      "Final LU mismatch.\n",
      "Only in soter.dat: ",
      paste(
        setdiff(soter_ids, lu2_ids),
        collapse = ", "
      ),
      "\nOnly in lu2.dat: ",
      paste(
        setdiff(lu2_ids, soter_ids),
        collapse = ", "
      )
    )
  }

  con <- DBI::dbConnect(
    odbc::odbc(),
    dsn = dsn,
    uid = user,
    pwd = password
  )

  on.exit(DBI::dbDisconnect(con), add = TRUE)

  db_lu_ids <- DBI::dbGetQuery(
    con,
    "SELECT pid FROM landscape_units ORDER BY pid"
  )$pid

  db_veg_ids <- DBI::dbGetQuery(
    con,
    "SELECT pid FROM vegetation ORDER BY pid"
  )$pid

  if (!setequal(soter_ids, hymo_ids)) {
    stop(
      "Final LU mismatch between soter.dat and hymo.dat.\n",
      "Only in soter.dat: ",
      paste(
        setdiff(soter_ids, hymo_ids),
        collapse = ", "
      ),
      "\nOnly in hymo.dat: ",
      paste(
        setdiff(hymo_ids, soter_ids),
        collapse = ", "
      )
    )
  }

  if (!setequal(soter_ids, db_lu_ids)) {
    stop(
      "Exported LU IDs do not match the cleaned ",
      "landscape_units database table."
    )
  }

  vegetation_out <- read_table_checked(
    vegetation_file,
    header = TRUE,
    sep = "",
    skip = 1
  )

  exported_veg_ids <- as.integer(
    vegetation_out[[1]]
  )

  if (!setequal(exported_veg_ids, db_veg_ids)) {
    stop(
      "Exported vegetation IDs do not match the ",
      "cleaned vegetation database table.\n",
      "Database IDs: ",
      paste(db_veg_ids, collapse = ", "),
      "\nExported IDs: ",
      paste(exported_veg_ids, collapse = ", ")
    )
  }

  message("\nFINAL VALIDATION PASSED")
  message(
    "LUs in hymo.dat, soter.dat and lu2.dat: ",
    length(soter_ids)
  )
  message(
    "Vegetation IDs exported: ",
    paste(sort(exported_veg_ids), collapse = ", ")
  )

  invisible(
    list(
      lu_ids = sort(soter_ids),
      vegetation_ids = sort(exported_veg_ids)
    )
  )
}


# =============================================================================
# 3. INITIALISE GRASS
# =============================================================================

initGRASS(
  gisBase = gisBase,
  home = getwd(),
  location = grass_location,
  mapset = grass_mapset,
  gisDbase = gisDbase,
  override = TRUE
)

lumpR:::test_grass()

drain_p <- read_VECT(
  vname = vname,
  layer = 1
)

print(drain_p)

drain_p_sp <- as(
  drain_p,
  "Spatial"
)

proj4string(drain_p_sp) <- CRS(
  getLocationProj()
)

print(drain_p_sp)


# =============================================================================
# 4. SUBBASIN DELINEATION
# =============================================================================

execGRASS(
  "g.region",
  raster = dem,
  flags = "a"
)

?calc_subbas # read the documentation!
calc_subbas(
  dem = dem,
  drain_points = drain_p_sp,
  river = river,
  basin_out = subbas,
  stream = stream_pref,
  points_processed = drainp_processed,
  outlet = outlet_row,
  thresh_stream = thresh_stream,
  thresh_sub = thresh_sub,
  snap_dist = snap_dist,
  rm_spurious = rm_spurious,
  keep_temp = keep_temp,
  export_shp = "sub_devoll.shp",
  overwrite = overwrite,
  silent = silent
)

# check the resulting map <subbas> in GRASS

## PREPROCESSING AND HILLSLOPE DEVIATION ####
?lump_grass_prep # read the documentation!
lump_grass_prep(
  # INPUT #
  mask = subbas,
  dem = dem,
  lcov = lcov,
  soil = soil,
  watermask = watermask,
  imperviousmask = imperviousmask,
  # OUTPUT #
  eha=eha,
  flowdir = flowdir,
  flowacc = flowacc,
  stream = stream,
  stream_horton = stream_horton,
  elevriv = elevriv,
  distriv = distriv,
  svc = svc,
  dir_out = getwd(),
  svc_ofile = svc_file,
  eha_thres = eha_thres,
  sizefilter = sizefilter,
  growrad = growrad,
  keep_temp = keep_temp,
  overwrite = overwrite,
  silent = silent,
  addon_path = addon_path
)

vegetation_overlap <- check_vegetation_overlap()

svc_definition <- read_table_checked(
  svc_file,
  header = TRUE,
  sep = ""
)

raster_veg_ids <- sort(
  unique(vegetation_overlap$vegetation_id)
)

svc_veg_ids <- sort(
  unique(svc_definition$veg_id)
)

if (!setequal(raster_veg_ids, svc_veg_ids)) {
  stop(
    "Vegetation IDs differ between the GRASS overlap ",
    "and soil_vegetation_components.dat.\n",
    "Raster IDs: ",
    paste(raster_veg_ids, collapse = ", "),
    "\nSVC IDs: ",
    paste(svc_veg_ids, collapse = ", ")
  )
}


# =============================================================================
# 6. SNOW PREPROCESSING
# =============================================================================

prepare_snow_input(
  dem = dem,
  subbas = subbas,
  eha = eha,
  flowdir = flowdir,
  eha_1d_file = snow_eha_file,
  keep_temp = keep_temp,
  overwrite = overwrite,
  silent = silent
)

assert_files_exist(snow_eha_file)

snow_tab <- read_table_checked(
  snow_eha_file,
  header = TRUE,
  sep = "\t"
)

required_snow_fields <- c(
  "eha_id",
  "aspect.sin",
  "aspect.cos",
  "rel_alt"
)

missing_snow_fields <- setdiff(
  required_snow_fields,
  names(snow_tab)
)

if (length(missing_snow_fields) > 0) {
  stop(
    "Missing field(s) in snow_eha.txt: ",
    paste(missing_snow_fields, collapse = ", ")
  )
}

prepare_snow_header(
  snow_eha_head_file
)


# =============================================================================
# 7. CALCULATE MEAN CATENAS
# =============================================================================

area2catena(
  mask = subbas,
  flowacc = flowacc,
  eha = eha,
  distriv = distriv,
  elevriv = elevriv,
  supp_quant = supp_quant,
  supp_qual = supp_qual,
  dir_out = getwd(),
  eha_2d_file = catena_out,
  eha_2d_head_file = catena_head_out,
  ridge_thresh = 1,
  min_cell_in_slope = min_cell_in_slope,
  min_catena_length = min_catena_length,
  max_riv_dist = max_riv_dist,
  plot_catena = plot_catena,
  grass_files = grass_files,
  ncores = ncores,
  eha_subset = NULL,
  overwrite = overwrite,
  silent = silent,
  allow_debug = TRUE
)


# =============================================================================
# 8. CONFIGURE CLASSIFICATION HEADER
# =============================================================================

header_path <- file.path(
  getwd(),
  catena_head_out
)

header_dat <- readLines(
  header_path,
  warn = FALSE
)

if (length(header_dat) < 9) {
  stop(
    catena_head_out,
    " has fewer than 9 lines."
  )
}

no_classes[1] <- -abs(no_classes[1])

header_dat[8] <- paste(
  no_classes,
  "\t",
  sep = "",
  collapse = ""
)

header_dat[9] <- paste(
  c(
    no_TCs,
    rep(
      0,
      length(no_classes) - 1
    )
  ),
  "\t",
  sep = "",
  collapse = ""
)

writeLines(
  header_dat,
  header_path
)

res_info <- execGRASS(
  "r.info",
  map = dem,
  flags = "g",
  intern = TRUE
)

res_values <- as.numeric(
  gsub(
    "[a-z]*=",
    "",
    grep(
      "nsres|ewres",
      res_info,
      value = TRUE
    )
  )
)

res <- mean(res_values)

if (
  !is.finite(res) ||
  res <= 0
) {
  stop("Could not determine a valid GRASS raster resolution.")
}


# =============================================================================
# 9. CLASSIFY EHAs INTO LUs AND TCs
# =============================================================================

prof_class(
  eha_2d_file = catena_out,
  eha_2d_head_file = catena_head_out,
  svc_column = "svc",
  eha_1d_file = snow_eha_file,
  eha_1d_head_file = snow_eha_head_file,
  dir_out = getwd(),
  luoutfile = "lu.dat",
  tcoutfile = "tc.dat",
  lucontainstcoutfile = "lucontainstc.dat",
  tccontainssvcoutfile = "tc_contains_svc.dat",
  terraincomponentsoutfile = "terraincomponents.dat",
  recl_lu = "reclass_lu.txt",
  saved_clusters = NULL,
  seed = 1312,
  resolution = res,
  classify_type = " ",
  max_com_length = 50,
  com_length = NULL,
  make_plots = plot_profclass,
  eha_subset = NULL,
  overwrite = overwrite,
  silent = silent
)

# IMPORTANT:
# Rebuild the TC-SVC relation from tc.dat. This replaces the incorrect
# relation created by the affected prof_class() workflow.
tc_svc_fixed <- repair_tc_contains_svc(
  tc_file = "tc.dat",
  svc_definition_file = svc_file,
  reclass_file = "reclass_svc.txt",
  output_file = "tc_contains_svc.dat"
)


# =============================================================================
# 10. POST-PROCESS CLASSIFICATION
# =============================================================================

lump_grass_post(
  mask = subbas,
  dem = dem,
  recl_lu = "reclass_lu.txt",
  lu = lu,
  subbasin = subbas,
  eha = eha,
  flowacc = flowacc,
  flowdir = flowdir,
  stream_horton = stream_horton,
  soil_depth = soil_depth,
  sdr = NULL,
  dir_out = getwd(),
  sub_ofile = sub_ofile,
  lu_ofile = lu_ofile,
  lupar_ofile = lupar_ofile,
  fill_holes = TRUE,
  groundwater = 0,
  keep_temp = keep_temp,
  overwrite = overwrite,
  silent = silent
)


# =============================================================================
# 11. RESERVOIR PARAMETERISATION
# =============================================================================

valid_reservoir_modes <- c(
  "none",
  "reuse_existing",
  "generate"
)

if (!reservoir_mode %in% valid_reservoir_modes) {
  stop(
    "Invalid reservoir_mode: ",
    reservoir_mode,
    ". Use one of: ",
    paste(valid_reservoir_modes, collapse = ", ")
  )
}

required_reservoir_files <- c(
  "reservoir.txt",
  "reservoirs_small_classes.dat",
  "r_subbas_contains_reservoirs_small.dat"
)

if (identical(reservoir_mode, "none")) {

  message(
    "No reservoir parameter files will be imported. ",
    "Reservoir database tables will remain empty."
  )

} else if (identical(reservoir_mode, "reuse_existing")) {

  assert_files_exist(
    file.path(
      dat_dir,
      required_reservoir_files
    )
  )

  message(
    "Reusing existing reservoir parameter files; ",
    "no reservoir GRASS processing will be run."
  )

} else if (identical(reservoir_mode, "generate")) {

  strategic_vector_info <- try(
    execGRASS(
      "g.findfile",
      element = "vector",
      file = strategic_reservoir_vector,
      intern = TRUE
    ),
    silent = TRUE
  )

  if (
    inherits(strategic_vector_info, "try-error") ||
    !any(
      grepl(
        "^name=.+",
        strategic_vector_info
      )
    )
  ) {
    stop(
      "Cannot generate strategic reservoirs because GRASS vector '",
      strategic_reservoir_vector,
      "' was not found.\n",
      "Set reservoir_mode <- \"reuse_existing\" to keep the existing ",
      "reservoir parameter files."
    )
  }

  reservoir_outlet(
    flowacc = flowacc,
    dem = dem,
    res_vct = strategic_reservoir_vector,
    outlets_vect = strategic_outlet_vector,
    keep_temp = TRUE,
    overwrite = TRUE
  )

  reservoir_strategic(
    res_vect = strategic_outlet_vector,
    res_file = file.path(
      dat_dir,
      "2bprepared",
      "from_reservoir_inventory",
      "reservoir_pars.csv"
    ),
    reservoir_file = "reservoir.txt",
    dir_out = getwd(),
    overwrite = TRUE,
    subbasin = subbas
  )

  small_vector_info <- try(
    execGRASS(
      "g.findfile",
      element = "vector",
      file = small_reservoir_vector,
      intern = TRUE
    ),
    silent = TRUE
  )

  if (
    inherits(small_vector_info, "try-error") ||
    !any(
      grepl(
        "^name=.+",
        small_vector_info
      )
    )
  ) {
    stop(
      "Cannot generate small-reservoir parameters because GRASS vector '",
      small_reservoir_vector,
      "' was not found.\n",
      "Set reservoir_mode <- \"reuse_existing\" to keep the existing ",
      "reservoir parameter files."
    )
  }

  reservoir_lumped(
    res_vect = small_reservoir_vector,
    subbas = subbas,
    res_vect_class = small_reservoir_class_vector,
    dir_out = getwd(),
    overwrite = TRUE
  )

  assert_files_exist(
    file.path(
      dat_dir,
      required_reservoir_files
    )
  )
}


# =============================================================================
# 12. COPY PREPARED SOIL AND VEGETATION TABLES
# =============================================================================

source_parameter_files <- c(
  file.path(veg_path, "vegetation.txt"),
  file.path(soil_path, "soil.dat"),
  file.path(soil_path, "horizons.dat"),
  file.path(soil_path, "particle_classes.dat"),
  file.path(
    soil_path,
    "r_soil_contains_particles.dat"
  )
)

assert_files_exist(source_parameter_files)

copy_map <- c(
  "vegetation.txt",
  "soil.dat",
  "horizons.dat",
  "particle_classes.dat",
  "r_soil_contains_particles.dat"
)

copy_ok <- file.copy(
  from = source_parameter_files,
  to = file.path(dat_dir, copy_map),
  overwrite = TRUE
)

if (!all(copy_ok)) {
  stop(
    "Failed to copy prepared parameter file(s): ",
    paste(
      source_parameter_files[!copy_ok],
      collapse = ", "
    )
  )
}


# In sediment mode, populate the four seasonal MUSLE-C columns before db_fill().
# In hydrology/snow mode this function makes no changes.
apply_musle_c_to_vegetation_file()


# =============================================================================
# 13. CREATE lu_db.dat ONCE, INCLUDING SNOW PARAMETERS
# =============================================================================

lu_database <- prepare_lu_database_file(
  lu_file = "lu.dat",
  lu_parameters_file = lupar_ofile,
  output_file = "lu_db.dat"
)


# =============================================================================
# 14. CONFIGURE, RESET, CREATE, AND FILL DATABASE
# =============================================================================

configure_odbc_dsn()
reset_database_tables()

db_create(
  dbname = dbname
)

db_update(
  dbname
)

# Core tables are always imported.
db_tables <- c(
  "subbasins",
  "r_subbas_contains_lu",
  "landscape_units",
  "r_lu_contains_tc",
  "terrain_components",
  "r_tc_contains_svc",
  "soils",
  "horizons",
  "soil_veg_components",
  "particle_classes",
  "r_soil_contains_particles",
  "vegetation"
)

db_files <- c(
  "sub_stats.txt",
  "lu_stats.txt",
  "lu_db.dat",
  "lucontainstc.dat",
  "terraincomponents.dat",
  "tc_contains_svc.dat",
  "soil.dat",
  "horizons.dat",
  "soil_vegetation_components.dat",
  "particle_classes.dat",
  "r_soil_contains_particles.dat",
  "vegetation.txt"
)

# Reservoir tables are imported only when reservoir data are available.
if (!identical(reservoir_mode, "none")) {
  db_tables <- c(
    db_tables,
    "reservoirs_strategic",
    "reservoirs_small_classes",
    "r_subbas_contains_reservoirs_small"
  )

  db_files <- c(
    db_files,
    "reservoir.txt",
    "reservoirs_small_classes.dat",
    "r_subbas_contains_reservoirs_small.dat"
  )
}

if (length(db_tables) != length(db_files)) {
  stop(
    "Internal error: db_tables and db_files have different lengths."
  )
}

assert_files_exist(
  file.path(dat_dir, db_files)
)

db_fill(
  dbname = dbname,
  tables = db_tables,
  dat_files = db_files,
  dat_dir = dat_dir,
  overwrite = TRUE,
  verbose = TRUE
)

status_after_fill <- database_status(
  "after db_fill"
)

expected_veg_ids <- sort(
  unique(svc_definition$veg_id)
)

if (
  !setequal(
    status_after_fill$active_vegetation_ids,
    expected_veg_ids
  )
) {
  stop(
    "Vegetation IDs were lost during db_fill.\n",
    "Expected: ",
    paste(expected_veg_ids, collapse = ", "),
    "\nActive after db_fill: ",
    paste(
      status_after_fill$active_vegetation_ids,
      collapse = ", "
    )
  )
}


# =============================================================================
# 15. DATABASE CLEANING
# =============================================================================

db_check(
  dbname,
  check = "check_fix_fractions",
  fix = TRUE,
  verbose = TRUE
)

if (apply_small_area_filter) {
  db_check(
    dbname,
    check = "filter_small_areas",
    option = list(
      area_thresh = small_area_threshold
    ),
    fix = TRUE,
    verbose = TRUE
  )
}

db_check(
  dbname,
  check = "tc_slope",
  option = list(
    treat_slope = c(3, 0.01, 0.1)
  ),
  fix = TRUE,
  verbose = TRUE
)

if (apply_remove_water_svc) {
  db_check(
    dbname,
    check = "remove_water_svc",
    fix = TRUE,
    verbose = TRUE
  )
}

if (apply_compute_rocky_frac) {
  db_check(
    dbname,
    check = "compute_rocky_frac",
    fix = TRUE,
    verbose = TRUE
  )
}

if (apply_remove_impervious_svc) {
  db_check(
    dbname,
    check = "remove_impervious_svc",
    fix = TRUE,
    verbose = TRUE
  )
}

db_check(
  dbname,
  check = "delete_obsolete",
  fix = TRUE,
  verbose = TRUE
)

db_check(
  dbname,
  check = "completeness",
  fix = TRUE,
  verbose = TRUE
)

db_check(
  dbname,
  check = "subbasin_order",
  fix = TRUE,
  verbose = TRUE,
  option = list(
    overwrite = TRUE
  )
)

status_after_cleaning <- database_status(
  "after database cleaning"
)

missing_after_cleaning <- setdiff(
  expected_veg_ids,
  status_after_cleaning$active_vegetation_ids
)

if (length(missing_after_cleaning) > 0) {
  stop(
    "Vegetation classes were removed during database cleaning: ",
    paste(missing_after_cleaning, collapse = ", "),
    "\nReduce small_area_threshold or disable ",
    "apply_small_area_filter."
  )
}


# =============================================================================
# 16. PREPARE MUSLE PARAMETERS
# =============================================================================

db_prepare_musle(
  dbname,
  compute_K = TRUE,
  verbose = TRUE
)

musle_copy_fields <- c(
  "Manning-n",
  "MUSLE-K"
)

if (enable_sediment) {
  musle_copy_fields <- c(
    "MUSLE-C",
    musle_copy_fields
  )
}

db_prepare_musle(
  dbname,
  compute_K = FALSE,
  setP = 1,
  copy_from_other_tables = musle_copy_fields,
  verbose = TRUE
)


# =============================================================================
# 17. EXPORT WASA-SED INPUT FILES ONCE
# =============================================================================

wasa_output_dir <- file.path(
  dat_dir,
  "WASA_input"
)

dir.create(
  wasa_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Remove a stale manually generated lu2.dat before final export.
stale_lu2 <- file.path(
  wasa_output_dir,
  "Hillslope",
  "lu2.dat"
)

if (file.exists(stale_lu2)) {
  file.remove(stale_lu2)
}

wasa_files <- c(
  "info.dat",
  "River/routing.dat",
  "River/response.dat",
  "Hillslope/hymo.dat",
  "Hillslope/soter.dat",
  "Hillslope/terrain.dat",
  "Hillslope/soil_vegetation.dat",
  "Hillslope/soil.dat",
  "Hillslope/lu2.dat",
  "Hillslope/vegetation.dat",
  "Hillslope/svc_in_tc.dat",
  "do.dat",
  "maxdim.dat",
  "part_class.dat",
  "Hillslope/soil_particles.dat",
  "Hillslope/rainy_season.dat",
  "Hillslope/x_seasons.dat",
  "Hillslope/svc.dat",
  "Reservoir/reservoir.dat",
  "Reservoir/lake.dat",
  "Reservoir/lake_number.dat",
  "Reservoir/lake_maxvol.dat"
)

db_wasa_input(
  dbname = dbname,
  dest_dir = wasa_output_dir,
  files = wasa_files,
  overwrite = TRUE,
  verbose = TRUE
)

# Normally db_wasa_input() creates lu2.dat correctly. This conditional safeguard
# rewrites it only if it is absent or does not match the cleaned LU table.
ensure_lu2_from_clean_database()

# db_wasa_input() currently writes one more data field than header label in
# soter.dat when sdr_lu is present. Correct only the missing header label.
fix_soter_header()

# Copy/validate an optional prepared rainy-season file.
rainy_season_status <- handle_rainy_season_file()

# Validate sediment variables according to enable_sediment.
sediment_status <- validate_sediment_outputs()


# =============================================================================
# 18. FINAL CONSISTENCY CHECKS
# =============================================================================

final_status <- validate_final_outputs()

message(
  "\nWorkflow completed successfully.\n",
  "WASA input directory: ",
  wasa_output_dir
)
