# Tutorial workflow for using lumpR for generating geodata input files for the hydrological model WASA-SED
# Jan 2026

# We recommend running this script line by line, so it is easier to understand and fix any errors that may occur.
# For help, please use the inbuilt R-help of the functions (e.g. ?calc_subbas) and the
# the wiki at https://github.com/tpilz/lumpR/wiki/


# PREREQUISITES

# Software:
# - R
# - GRASS GIS, ver. 8, including extensions r.stream.distance and r.stream.order
# - data base (e.g. sqlite, MariaDB, LibreOffice Base or other DBs with ODBC-driver)
# Data:
# - prepared GRASS location + mapset with DEM, soil and vegetation data (see example/data_import.sh for instructions)


# install package and required dependencies ####
#(these need to be run only once)
if (!require("devtools"))
  install.packages("devtools") 
library(devtools)
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE) #tell git_install() to ignore warnings. Otherwise, it gets stuck at each warning
if (!require("lumpR"))
  #install_github("tPilz/lumpR") #install main branch
  install_github("tpilz/lumpR", ref="wo_rgeos_maptools") #install specific branch 
if (!require("rgrass"))
  install.packages("rgrass")


# Load required packages ####
library(lumpR)
library(rgrass)


# SETTINGS ####
## GRASS-GIS-SETTINGS ####

# switch to specified working directory (this is usually the home directory of this very script)
if (!require("rstudioapi", quietly = TRUE)) install.packages("rstudioapi")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#setwd("/home/myname/somefolder") #manual specification. Use "/" instead of "\" in Windows


# initialisation of GRASS session
# Set the path to your GRASS installation and the GRASS location
addon_path=NULL  # normally not required: path to your locally installed GRASS add-ons, e.g. "/home/myuser/.grass8/addons/bin". Must only be given if necessary, see ?lump_grass_prep
gisBase = "/usr/lib/grass83"  # path to GRASS installation, e.g. "/usr/lib/grass83" (use / instead of / under windows, e.g. "d:/programme/GRASS7.0" )
#gisDbase="/home/francke/grassdata"  # path to 'grassdata' directory containing the location specified above and all corresp. data. 
gisDbase="/media/francke/daten/till/uni/grass-db/"  # path to 'grassdata' directory containing the location specified above and all corresp. data. 
#If unknown, run "g.gisenv -n" in GRASS and locate variable GISDBASE
location="isabena_tutorial" # GRASS location
mapset="PERMANENT"    # corresponding GRASS mapset

initGRASS(gisBase=gisBase, 
          home=getwd(), # The directory in which to create the .gisrc file
          location=location,
          mapset=mapset,
          gisDbase=gisDbase,  
          override=TRUE)

#test if GRASS connections working
lumpR:::test_grass()


## INPUT-settings ####
# inputs marked MANDATORY have to be given, the rest can be 'NULL' if not available

# define subbasin outlet points (choose EITHER option A or B, not both!)
# Import the outlet points from a vector layer 

vname = "outlet_points" #specify name of point GRASS vector containing subbasin coordinates
drain_p = read_VECT(vname = vname, layer=1)  #load outlet points from GRASS
print(drain_p)
outlet_row = 1 #specify the row of outlet point of entire watershed


dem = "dem" # DEM raster - MANDATORY
lcov = "vegetation" # land / vegetation cover raster map in GRASS location (e.g. derived from CORINE) - MANDATORY
soil = "soils" # soil raster map in GRASS location - MANDATORY (The package [SoilDataPrep](https://github.com/TillF/SoilDataPrep/tree/master) is designed to help you with the generation of soil data, especially in conjunction with lumpR) 


soil_depth = NULL # soil depth raster map
watermask = NULL # water mask raster map in GRASS location (1=water, 0=no water)
imperviousmask = NULL # impervious surface areas raster map in GRASS location (1=impervious, 0=permeable)
river = NULL # river vector map

# soil and vegetation properties
# path to prepared vegetation parameter table 'vegetation.dat' in WASA format - MANDATORY
#'vegetation.dat' can be e.g. be generated from CORINE data and the mapping table provided in 2bprepared/vegetation_from_corine
veg_path = "2bprepared/vegetation_from_corine" 

# path to prepared soil parameter tables 'horizons.dat, 'particle_classes.dat', 'r_soil_contains_particles.dat', and 'soil.dat' in WASA format - MANDATORY
#these files can be used from prior WASA-parameterizations or generated with https://github.com/TillF/SoilDataPrep
soil_path = "2bprepared/from_SoilDataPrep"


## OUTPUT-settings ####
# outputs marked MANDATORY have to be given, the rest can be NULL
# some outputs are inputs for other functions (e.g. flow accumulation)
# These setting describe the names of the respective maps generated in GRASS. Normally, they do not need to be changed.

subbas = "subbasin" # subbasin raster map generated - MANDATORY

# prefix of calculated stream segments raster and vector maps
stream_pref = "stream_accum"

# prefix of drainage point vector files (given points snapped to river and internally calculated points, respectively) - MANDATORY
drainp_processed = "drain_points"

# elementary hillslope areas raster map - MANDATORY
eha = "eha"

# flow direction raster map - MANDATORY
flowdir = "flowdir"

# flow accumulation raster map - MANDATORY
flowacc = "flowacc"

# stream segments raster map (based on eha calculation; much finer than 'stream_pref' to delineate hillslopes) - MANDATORY
stream = "stream"

# Horton stream order raster map (based on 'stream' above) - MANDATORY
stream_horton = "stream_horton"

# elevation relative to next river cell raster map - MANDATORY
elevriv = "elevriv"

# distance to next river cell raster map - MANDATORY
distriv = "distriv"

# soil vegetation components raster map - MANDATORY
svc = "svc"

# landscape units raster map - MANDATORY
lu = "lu"

# name of file containing svc parameters - MANDATORY
svc_file = "soil_vegetation_components.dat"

# Name of file containing catena information derived from EHAs as output of area2catena() and input to prof_class() - MANDATORY
eha_2d_file = "eha_2d.txt"
# Name of output header file containing meta-information as input for prof_class - MANDATORY
eha_2d_head_file = "eha_2d_head.txt"

#optional file containing additional properties at the EHA-scale. Can be generated e.g. with prepare_snow_input()
eha_1d_file="eha_1d.txt" 
#corresponding header file
eha_1d_head_file="eha_1d_head.txt" 


# Name of subbasin statistics file containing subbasin parameters
sub_ofile = "sub_stats.txt"

# Name of file containing subbasins and the corresponding LUs with their fraction of area in the subbasin
lu_ofile = "lu_stats.txt"

# Name of file containing LUs and related parameters
lupar_ofile = "lu_pars.txt"

# Name for the vector reservoir map created in GRASS location containing information on reservoir size classes in the attribute table
# If NULL it will not be created
res_vect_class = "res_vect_class"



## lumpR- PARAMETERS ####
# STRONGLY CASE-STUDY SPECIFIC! Try out what suits your needs / data

### MISCELLANEOUS PARAMETERS ####
# parameters influencing some outputs but not directly discretisation complexity

# parameters to calc_subbas. For details, see ?calc_subbas
  # Threshhold for derivation of stream from flow accumulation map (in cells). Needs to be set only if river is not set - MANDATORY
  thresh_stream = 8000
  
  # maximum distance for snapping of drain_points to stream cells in the units of your GRASS location (usually meters) - MANDATORY
  snap_dist = 500
  
  # threshold for small spurious subbasins created within calc_subbas() to be removed. 
  rm_spurious = 0.01  

# parameters to lump_grass_prep.  For details, see ?lump_grass_prep  
  # minimum size of EHAs (in map units, usually m²) to be preserved, smaller EHAs (artefacts) are removed; parameter for GRASS function r.reclass.area - MANDATORY
  sizefilter = 30 
  
  # growing radius (in raster cells) to remove artefacts in EHA data; parameter for GRASS function r.grow, see see ?lump_grass_prep - MANDATORY
  growrad = 50
  
# parameters to area2catena  For details, see ?area2catena    
  # minimum number of cells a hillslope area must have, all smaller ones are skipped
  min_cell_in_slope = 30
  
  # minimum number of sampling points (cells) a catena should have. If there are less, the catena is not saved
  min_catena_length = 3

  # maximum distance to river [in cells]: if the closest cell of an EHA is farther than max_riv_dist, the EHA is skipped, otherwise all distances within the EHA are redurced by the distance of the closest cell to river
  max_riv_dist = 15

### LANDSCAPE DISCRETISATION PARAMETERS ####

# Parameter for GRASS function r.watershed defining the minimum size of an exterior watershed basin in number of grid cells. If NULL only the given drainage points are used for subbasin delineation
thresh_sub = NULL

# parameter for GRASS function r.watershed. This is a crucial parameter affecting the size of delineated hillslopes - MANDATORY
eha_thres = 300

# named vector with names of GRASS raster maps containing *quantitative* attributes for LU deviation and
#corresponding number of classes, that should be created from the quantitative attributes
supp_quant_no_classes = c(
  "aspect" = 4, 
  "rel_alt" =2
  ) 

# named vector with GRASS file names of *qualitative* attributes for LU deviation 
# and corresponding number of classes, that should be created from the qualitative attributes
supp_qual_no_classes = c(
  "svc"  =1, # svc has to be defined to generate SVC parameters needed for WASA! Map is generated by lump_grass_prep()
  "soils"=3
  ) 

# number of LUs to be produced per attribute in the order:
no_classes = c(-3, 3, 10, supp_quant_no_classes, supp_qual_no_classes, 1)
names(no_classes) = c("shape", "extent", "weighting", names(supp_quant_no_classes), names(supp_qual_no_classes), "slope_width")
#shape: normalized profile shape of catena
#extent: vertical and horizontal extent of catena
#weighting: weighting of vertical vs. horizontal extent

# number of TCs that are created by algorithm
no_TCs = 3

### RUN TIME PARAMETERS #####

# shall temporary files be kept in the GRASS location, e.g. for debugging or further analyses?
keep_temp = FALSE

# Shall existing outputs of previous calls of functions be deleted? If FALSE the function returns an error if output already exists
overwrite = TRUE

# Shall the function be silent (also suppressing warnings of internally used GRASS functions, etc.)?
silent = FALSE

# produce plots (scatter, mean catena, etc.) for each area / class (written into sub-directory 'plots_area2catena')
plot_catena = TRUE

# produce plots of classification of catenas to landscape units and terrain components
plot_profclass = TRUE

# produce GRASS reclassification files for qualitative raster data
grass_files = TRUE

# number of cores that should be used for parallel computation (where possible)
ncores = 4


## DATABASE-settings ####
# You will need to have set up an ODBC database.  (However, the database is only needed in the very last steps of the processing chain, so you can do this later).
# For details, check https://github.com/tpilz/lumpR/wiki/04-Databases-and-ODBC
# Define the database name

dbname = "tutorial_sqlite.db"

# Construct the ODBC connection string
str_odbc = c(paste0("[", dbname, "]"),
              "Description = lumpR tutorial DB",
              "Driver = SQLite3",
              "Server = localhost",
              "Port = 3306",
              paste0("Database = ", getwd(), "/", dbname),  # Adjust database name as needed
              "User = ",          # Adjust username as needed
              "Password = ")  # Adjust password as needed

# Write the ODBC connection string to .odbc.ini file (required for Linux only once)
#write(str_odbc, file="~/.odbc.ini", ncolumns=1, append=TRUE, sep="\n")


# CALCULATIONS ####
# usually, no adjustments are needed below this line
# run line-by-line to understand what is going on!

drain_p_sp = as(drain_p, "Spatial") # convert SpatVector to SpatialPoints
proj4string(drain_p_sp) = CRS(getLocationProj())
print(drain_p_sp)

##  SUBBASIN DELINEATION ####
# calculate subbasins; one subbasin for each drainage point 

## execGRASS("r.mask", flags = "r") # Set the MASK to null

# Set the GRASS region to your DEM, force alignment
execGRASS("g.region", raster = dem, flags = "a")

# Clean DEM: replace non-finite values (NaN, Inf) with NULL
execGRASS("r.mapcalc",
           expression = sprintf("dem_clean = if(isnull(%s), null(), %s)", dem, dem),
           flags = "overwrite")


?calc_subbas # read the documentation!

calc_subbas(
  # INPUT #
  dem=dem,
  drain_points=drain_p_sp, 
  river=river,
  # OUTPUT #
  basin_out=subbas,
  stream=stream_pref,
  points_processed=drainp_processed,
  # PARAMETERS #
  outlet=outlet_row,
  thresh_stream=thresh_stream,
  thresh_sub=thresh_sub,
  snap_dist=snap_dist,
  rm_spurious=rm_spurious,
  keep_temp=keep_temp,
  export_shp="sub_isabena.shp", #optional, can be used as input to SoilDataPrep
  overwrite=overwrite,
  silent=silent
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
  # PARAMETERS #
  eha_thres = eha_thres,
  sizefilter = sizefilter,
  growrad = growrad,
  keep_temp=keep_temp,
  overwrite=overwrite,
  silent=silent,
  addon_path=addon_path
)

# check the resulting maps <eha>, <stream> in GRASS


## PREPARATION OF SNOW INPUT FILES ####
prepare_snow_input(
  dem=dem,
  subbas=subbas,
  eha=eha,
  flowdir=flowdir,
  eha_1d_file=eha_1d_file,
  keep_temp=T,
  overwrite=T,
  silent=F
) 

# check the file eha_1d.txt, which should contain aspect and elevation information for each EHA

## CALCULATE MEAN CATENA FOR HILLSLOPES ####
?area2catena # read the documentation!
area2catena(
  # INPUT #
  mask=subbas,
  flowacc=flowacc,
  eha=eha,
  distriv=distriv,
  elevriv=elevriv,
  names(supp_quant_no_classes)=names(supp_quant_no_classes),
  names(supp_qual_no_classes)=names(supp_qual_no_classes),
  # OUTPUT #
  dir_out=getwd(),
  eha_2d_file=eha_2d_file,
  eha_2d_head_file=eha_2d_head_file,
  # PARAMETERS #
  ridge_thresh=1,
  min_cell_in_slope=min_cell_in_slope,
  min_catena_length=min_catena_length,
  max_riv_dist=max_riv_dist,
  plot_catena=plot_catena,
  grass_files=grass_files,
  ncores=ncores,
  eha_subset=NULL,
  overwrite=overwrite,
  silent=silent,
  allow_debug = TRUE
)
# Check plots in ./plots_area2catena



## CATENA CLASSIFICATION INTO LANDSCAPE UNITS AND TERRAIN COMPONENTS ####
# change header file according to desired number of classes to create in classification (will be used by prof_class())
modify_eha_head_files <- function(eha_2d_head_file, eha_1d_file=NULL, no_classes) {
  if (!identical(names(no_classes[1:3]), c("shape", "extent", "weighting"))) {
    stop("Argument no_classes needs to be a integer vector with its first three elements named 'shape', 'extent',  'weighting'.")
  }
  
  #read 2D-file
  header_2d_comments = readLines(eha_2d_head_file, n=5) #read header comments only
  
  header_2d_dat = read.table(eha_2d_head_file, sep="\t", skip = 5, fill=TRUE, header=TRUE) #we need to set fill=TRUE as there may be empty fields trailing. However, this shifts the first column to the column names
  if(!all(diff(as.numeric(row.names(header_2d_dat)))==1)) #are the row names not sequentially numbered? Then, fix table structure.
    header_2d_dat[,1:ncol(header_2d_dat)] = cbind(as.numeric(row.names(header_2d_dat)), as.matrix(header_2d_dat)[,-ncol(header_2d_dat)]) #shift columns names back to first column
  attr_2d = names(header_2d_dat)
  attr_2d[1:3] = c("shape", "extent", "weighting") #rename, as they are named differently in the outpu of area2catena
  
  #read 1D-file
  if (!is.null(eha_1d_file))
  {
    data_1d_dat = read.table(eha_1d_file, sep="\t", skip = 0, fill=TRUE, header=TRUE, nrow=1) #read data file to infer its structure
    attr_1d = names(data_1d_dat)[-1] #ignore first column as this contains the EHA-id
    attr_1d = gsub(pattern = "\\..*", repl="", x=attr_1d) #remove parts behind the .  - these are the "parent" attributes, e.g. "aspect.sin" and "aspect.cos" belong to parent "aspect"
    column_counts = table(attr_1d) #count columns used for each parent attribute
    attr_1d = names(column_counts) # remove duplicates
  } else
    attr_1d = NULL
  
  #check for missing attributes
  unspec_attributes = setdiff(c(attr_2d, attr_1d), names(no_classes)) #find attributes in header file not specified in no_classes
  if (length(unspec_attributes)>0) 
    stop(paste0("The following attributes are contained in '", eha_2d_head_file, "' and/or '", eha_1d_file, "' but not contained in 'no_classes':", paste(unspec_attributes, collapse=", "),
                "."))
  
  missing_cols = setdiff(names(no_classes), c(attr_2d, attr_1d)) #find columns not yet contained in the file
  if (length(missing_cols)>0) 
    stop(paste0("The following attributes are contained in 'no_classes', but neither in '", eha_2d_head_file, "' nor '", eha_1d_file, "' :", paste(missing_cols, collapse=", "),
                "."))
  
  #modify 2d-header file
  #header_2d_dat[, missing_cols] = NA #add missing columns
  header_2d_dat[2, 1:3] = no_classes[1:3] #set number of LU-classes for "shape", "extent" and "weighting"
  header_2d_dat[2, 1] = abs(header_2d_dat[1,1]) * -1 #force "successive" mode
  header_2d_dat[2, attr_2d[-(1:3)]] = no_classes[attr_2d[-(1:3)]] #set number of classes for remaining 2D-attributes
  #header_2d_dat[3,] = ... #todo: add writing TC-weighting factors
  
  writeLines(text=header_2d_comments, con=eha_2d_head_file)
  options(warn=-1) #suppress warning about row names)
  write.table(x=header_2d_dat, file=eha_2d_head_file, append=TRUE,  row.names=FALSE, quote=FALSE, sep="\t")
  options(warn=1) #restore warnings
  message(paste0(eha_2d_head_file, " rewritten."))
  
  #write 1d-header file
  if (!is.null(eha_1d_file))
  {
    header_1d_dat = array(NA, dim=c(2, length(attr_1d)+1), dimnames=list(c("columns", "no_classes"), c("eha_id", attr_1d)))
    
    header_1d_dat["columns", 1] = 1 #set number of columns used eha_id
    header_1d_dat["columns", names(column_counts)] = column_counts #set number of columns used for each attribute
    header_1d_dat["no_classes",  names(column_counts)] = no_classes[ names(column_counts)] #set number of classes for 1D-attributes
    eha_1d_head_file = gsub(pattern="\\.([^\\.]*)", repl="_head.\\1", x=eha_1d_file)
    write.table(x=header_1d_dat, file=eha_1d_head_file, row.names=FALSE, quote=FALSE, sep="\t")
    
    message(paste0(eha_1d_head_file, " rewritten."))
  }
}

modify_eha_head_files(eha_2d_head_file, eha_1d_file=eha_1d_file, no_classes) 

# get resolution of GRASS mapset (mean between x and y resolution)
resolution = execGRASS("r.info", map=dem, flags=c("g"), intern=TRUE)
resolution = sum(as.numeric(gsub("[a-z]*=", "", grep("nsres|ewres", resolution, value = TRUE)))) / 2

?prof_class # read the documentation!
prof_class(
  # INPUT #
  eha_2d_file=eha_2d_file,
  eha_2d_head_file=eha_2d_head_file,
  svc_column="svc",
  eha_1d_file=eha_1d_file,
  eha_1d_head_file=eha_1d_head_file,
  # OUTPUT #
  dir_out=getwd(),
  luoutfile="lu.dat",
  tcoutfile="tc.dat",
  lucontainstcoutfile="lucontainstc.dat",
  tccontainssvcoutfile="tc_contains_svc.dat",
  terraincomponentsoutfile="terraincomponents.dat",
  recl_lu="reclass_lu.txt",
  saved_clusters=NULL,
  # PARAMETERS #
  seed=1312,
  resolution=resolution,
  classify_type=' ',
  max_com_length=50,
  com_length=NULL,
  make_plots=plot_profclass,
  eha_subset=NULL,
  #eha_subset=1:100,
  overwrite=overwrite,
  silent=silent,
  plot_silhouette = FALSE
)
# Check plots in ./plots_prof_class



## POST PROCESSING OF CLASSIFICATION ####
?lump_grass_post # read the documentation!
lump_grass_post(
  # INPUT #
  mask = subbas,
  dem = dem,
  recl_lu = "reclass_lu.txt",
  lu = lu,
  subbasin = subbas,
  eha = eha,
  flowacc = flowacc,
  flowdir = flowdir,
  stream_horton = stream_horton,
  soil_depth = NULL,
  sdr=NULL,
  # OUTPUT #
  dir_out = getwd(),
  sub_ofile = sub_ofile,
  lu_ofile = lu_ofile,
  lupar_ofile = lupar_ofile,
  # PARAMETER #
  fill_holes=T,
  groundwater=0,
  keep_temp = keep_temp,
  overwrite = overwrite,
  silent = silent
)


#generate reservoir parameterisation
## RESERVOIR PARAMETERISATION ####
  #strategic/large reservoirs:
  x = reservoir_outlet(flowacc = flowacc, dem = dem, res_vct = "reservoir_vec", outlets_vect = "res_outlets", keep_temp = TRUE, overwrite = TRUE) 
  #generate reservoir parameter file for later import into db
  reservoir_strategic(res_vect = "res_outlets", res_file="./2bprepared/from_reservoir_inventory/reservoir_pars.csv", reservoir_file = "reservoir.txt", dir_out = getwd(), overwrite = TRUE, subbasin = subbas)
  
  #small / distributed reservoirs
  reservoir_lumped(res_vect = "res_small", subbas = subbas, res_vect_class = "res_small_2", dir_out =  "./", overwrite = TRUE)

  
    
# DATABASE ####
# rainy_season not included within this example

# create database
?db_create
odbcDataSources()
db_create(dbname=dbname)
  
# update database (if necessary)
?db_update
db_update(dbname)

# create information file for filling landscape_units into database (SIMPLEST POSSIBLE PARAMETER VALUES USED HEREIN)
luout = read.table("lu.dat", header=TRUE)
#add derived snow parameters
luout$aspect = atan(luout$aspect_p1_c2 / luout$aspect_p1_c1) / 2 / pi * 360 #convert aspect.sin and aspect.cos to aspect [deg]
names(luout)[names(luout)=="x_length"] = "slopelength" #rename column
names(luout)[names(luout)=="rel_alt_p1"] = "relative_altitude" #rename column
lupar = read.table(lupar_ofile, header=TRUE)
lupar$slopelength = NULL
lupar = merge(lupar, luout[,c("LU.ID", "slopelength", "aspect", "relative_altitude")], by.x="pid", by.y="LU.ID") 
lupar$soil_depth = -1 # groundwater option I.1.1 (WASA documentation) 
lupar$allu_depth = -1 # groundwater option I.1.1 (WASA documentation) 
lupar$riverbed_depth = 2000 # riverbed in any case below soil (no information whether this is reasonable or not)
lupar$kf_bedrock = -9999
lupar$gw_dist = -9999
lupar$frgw_delay = -9999
write.table(lupar, "lu_db.dat", quote = F, row.names = F, sep="\t")

# copy soil and vegetation parameter files into output_dir
file.copy(paste(veg_path, "vegetation.txt", sep="/"), "vegetation.txt", overwrite=T)
file.copy(paste(soil_path, "soil.dat", sep="/"), "soil.dat", overwrite=T)
file.copy(paste(soil_path, "horizons.dat", sep="/"), "horizons.dat", overwrite=T)
file.copy(paste(soil_path, "particle_classes.dat", sep="/"), "particle_classes.dat", overwrite=T)
file.copy(paste(soil_path, "r_soil_contains_particles.dat", sep="/"), "r_soil_contains_particles.dat", overwrite=T)

# lumpR output and manually prepared information (e.g. soil parameters) to database
?db_fill
db_fill(dbname=dbname,
        tables = c("subbasins", "r_subbas_contains_lu", 
                   "landscape_units", "r_lu_contains_tc", "terrain_components", "r_tc_contains_svc",
                   "soils", "horizons", "soil_veg_components",
                   "particle_classes", "r_soil_contains_particles","vegetation", "reservoirs_strategic",
                   "reservoirs_small_classes", "r_subbas_contains_reservoirs_small"),
        dat_files=c("sub_stats.txt", "lu_stats.txt", 
                    "lu_db.dat", "lucontainstc.dat", "terraincomponents.dat", "tc_contains_svc.dat",
                    "soil.dat", "horizons.dat", "soil_vegetation_components.dat",
                    "particle_classes.dat", "r_soil_contains_particles.dat", "vegetation.txt", "reservoir.txt",
                    "reservoirs_small_classes.dat", "r_subbas_contains_reservoirs_small.dat"), 
        dat_dir=getwd(),
        overwrite=TRUE, verbose=TRUE)


# Please process these cleaning actions step-by-step according to your needs.
?db_check

db_check(dbname, 
         check=c("check_fix_fractions"), 
         fix=T,
         verbose=T)

db_check(dbname, 
         check=c("filter_small_areas"), 
         option=list(area_thresh=0.01),
         fix=T,
         verbose=T)

db_check(dbname, 
         check=c("tc_slope"), 
         option=list(treat_slope=c(3,0.01,0.1)),
         fix=T,
         verbose=T)

# db_check(dbname, 
#          check=c("special_areas"), 
#          option=list(special_area = data.frame(reference_tbl=c("vegetation", "vegetation", "soils"), ref_id=c(3,4,10), special_id=c(1,1,2))),
#          fix=F,
#          verbose=F)

db_check(dbname, 
         check=c("remove_water_svc"),
         fix=T,
         verbose=T)

db_check(dbname, 
         check=c("compute_rocky_frac"),
         fix=T,
         verbose=T)

db_check(dbname, 
         check=c("remove_impervious_svc"),
         fix=T,
         verbose=T)

# db_check(dbname, 
#          check=c("proxy_frgw_delay"), 
#          option=list(total_mean_delay=50),
#          fix=T,
#          verbose=T)

db_check(dbname,
         check=c("delete_obsolete"),
         fix=T,
         verbose=T)

db_check(dbname,
         check=c("completeness"),
         fix=T,
         verbose=T)

db_check(dbname,
         check=c("subbasin_order"),
         fix=T,
         verbose=T,
         option=list(overwrite=TRUE))



# prepare sediment modelling by computing MUSLE factors
db_prepare_musle(dbname, compute_K = TRUE, verbose = TRUE) #compute K-factor
db_prepare_musle(dbname, compute_K = FALSE, setP = 1,
                 copy_from_other_tables = c("MUSLE-C", "Manning-n", "MUSLE-K"), verbose =
                   TRUE) #copy factors to table 'soil_veg_components'

# Check data in database

# Generate input files for WASA ####
?db_wasa_input
db_wasa_input(dbname = dbname,
              dest_dir = paste(getwd(), "WASA_input", sep="/"),
              files=c("info.dat", "River/routing.dat", "River/response.dat", "Hillslope/hymo.dat",
                      "Hillslope/soter.dat", "Hillslope/terrain.dat", "Hillslope/soil_vegetation.dat",
                      "Hillslope/soil.dat", "Hillslope/vegetation.dat", "Hillslope/svc_in_tc.dat",
                      "do.dat", "maxdim.dat", "part_class.dat", "Hillslope/soil_particles.dat", 
                      "Hillslope/rainy_season.dat", "Hillslope/x_seasons.dat", "Hillslope/svc.dat",
                      "Reservoir/reservoir.dat","Reservoir/lake.dat", "Reservoir/lake_number.dat",
                      "Reservoir/lake_maxvol.dat"),
              overwrite = overwrite, verbose=T)

# Next steps: 
# - manually adjust model input files to your needs ...
# - prepare additional WASA-SED input data (meteo data, rainy_season etc.) See WASA-SED documentation (e.g. https://tillf.github.io/WASA-SED/)


