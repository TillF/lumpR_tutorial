# command collection for data input for preparing GRASS for use with LUMP, dataset "tutorial"
# These are just examples! adjust to your needs!
# Most steps can also be performed via the GUI.

#change to working directory
cd /media/francke/daten/till/uni/lehre/wasa/ver5/1_lump/scripts/


#manually create a new location using GUI (Settings -> GRASS working env -> Create New Location) or using the command below:
   r.in.gdal input=./2bimported/dem.tif output=dem -e location=isabena_tutorial 	#import geoTif into GRASS and create new location 

#switch to new location
   g.gisenv set="LOCATION_NAME=isabena_tutorial"

#import overall and subbasin outlet points
	v.in.ogr input=./2bimported/outlet_points.shp output=outlet_points -r

#import soil map (The package [SoilDataPrep](https://github.com/TillF/SoilDataPrep/tree/master) is designed to help you with the generation of soil data, especially in conjunction with lumpR) 
	r.in.gdal input=./2bimported/soils2.tif output=soils -o

#vegetation / landuse
	v.in.ogr input=./2bimported/landuse.shp output=vegetation_vec min_area=0.0001 type=boundary snap=-1 -r
	v.to.rast input=vegetation_vec layer=1 output=vegetation use=attr attribute_column=wasa_id #convert vegetation map to raster

#import large / strategic reservoirs (ploygon with column(s) 'res_id') 
	v.in.ogr input=./2bimported/reservoir.shp output=reservoir_vec type=boundary snap=-1 --overwrite
	#v.to.rast input=reservoir_vec layer=1 output=reservoir use=val val=1

#import small / distributed reservoirs
	v.in.ogr input=./2bimported/reservoirs_small.shp output=res_small type=boundary snap=-1 --overwrite
