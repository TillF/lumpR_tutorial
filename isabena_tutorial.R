# Tutorial workflow for using lumpR for generating geodata input files for the hydrological model WASA-SED
# Dec 2025

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
lcov = "vegetation" # land / vegetation cover raster map in GRASS location - MANDATORY
soil = "soils" # soil raster map in GRASS location - MANDATORY (The package [SoilDataPrep](https://github.com/TillF/SoilDataPrep/tree/master) is designed to help you with the generation of soil data, especially in conjunction with lumpR) 


soil_depth = NULL # soil depth raster map
watermask = NULL # water mask raster map in GRASS location (1=water, 0=no water)
imperviousmask = NULL # impervious surface areas raster map in GRASS location (1=impervious, 0=permeable)
river = NULL # river vector map

# soil and vegetation properties
# path to prepared vegetation parameter table 'vegetation.dat' in WASA format - MANDATORY
#'vegetation.dat' can be used from prior WASA-parameterizations or generated with Alban's scripts
veg_path = "2bprepared/from_AlbansScripts" 

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

# Name of output file containing mean catena information as input for prof_class - MANDATORY
catena_out = "rstats.txt"

# Name of output header file containing meta-information as input for prof_class - MANDATORY
catena_head_out = "rstats_head.txt"

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

# vector with names of GRASS raster maps containing *quantitative* attributes for LU deviation (adjust no_classes accordingly)
supp_quant = c()
supp_quant_no_classes = c() #corresponding number of classes, that should be created from the quantitative attributes

# vector with GRASS file names of *qualitative* attributes for LU deviation (adjust no_classes accordingly)
supp_qual = c("svc", "soils") # svc has to be defined to generate SVC parameters needed for WASA! Map is generated by lump_grass_prep()
supp_qual_no_classes = c(1,3) #corresponding number of classes, that should be created from the qualitative attributes

# number of LUs to be produced per attribute in the order:
# c( <shape>, <extent>, <weighting vertical vs. horizontal extent>, <quant. suppl. information>, <qual. suppl. information>, <slope width> )
no_classes = c(-3, 3, 10, supp_quant_no_classes, supp_qual_no_classes, 1)

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

# Write the ODBC connection string to .odbc.ini file (required for Linux only)
write(str_odbc, file="~/.odbc.ini", ncolumns=1, append=TRUE, sep="\n")


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

## CALCULATE MEAN CATENA FOR HILLSLOPES ####
?area2catena # read the documentation!
area2catena(
  # INPUT #
  mask=subbas,
  flowacc=flowacc,
  eha=eha,
  distriv=distriv,
  elevriv=elevriv,
  supp_quant=supp_quant,
  supp_qual=supp_qual,
  # OUTPUT #
  dir_out=getwd(),
  catena_out=catena_out,
  catena_head_out=catena_head_out,
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
# change header file according to desired number of classes to create in classification (input for prof_class)
header_dat = readLines(paste(getwd(), catena_head_out, sep="/"))
no_classes[1] = abs(no_classes[1]) * -1 #force "succesive" mode
header_dat[8] = paste(no_classes, "\t", sep="", collapse="")
header_dat[9] = paste(c(no_TCs, rep(0, length(no_classes)-1)), "\t", sep="", collapse="")
writeLines(header_dat,paste(getwd(), catena_head_out, sep="/"))



# get resolution of GRASS mapset (mean between x and y resolution)
res = execGRASS("r.info", map=dem, flags=c("g"), intern=TRUE)
res = sum(as.numeric(gsub("[a-z]*=", "", grep("nsres|ewres", res, value = T)))) / 2

?prof_class # read the documentation!
prof_class(
  # INPUT #
  catena_file=catena_out,
  catena_head_file=catena_head_out,
  svc_column="svc",
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
  resolution=res,
  classify_type=' ',
  max_com_length=50,
  com_length=NULL,
  make_plots=plot_profclass,
  eha_subset=NULL,
  overwrite=overwrite,
  silent=silent
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
luout = read.table("lu.dat", header=T)
lupar = read.table(lupar_ofile, header=T)
lupar$slopelength = luout$x_length
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
        overwrite=T, verbose=T)


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


