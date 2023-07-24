#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 10:12:11 2018

@author: antonis
"""

#from pyhdf.SD import SD,SDC
#from pyresample import kd_tree, geometry,bilinear
#from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset

import datetime as dt
import os
import fnmatch
#import numpy.ma as ma
import numpy as np
#import itertools
import pandas as pd
import calendar
import pyresample as pr
import time
import scipy #pyresample for some reason requires it - at some point it was redacted


scriptstart=dt.datetime.now()

#f = open("/mnt/c/Users/Alexiou_NOA/Desktop/n2_2019.log", "w")

print (scriptstart)
#f.write(str(scriptstart) + '\n')
#f.flush()

#Define number of processors
nprocs=1

#Define here the year
year='2013'

#Define here the name of the satellite
#satellite='MODIS-AQUA' #change to Modis Aqua?
satellite='MODIS-TERRA'

#Define the least significant digits and compression level for the variables which will be stored in the netcdf files
lstsigdig=3
complevel=9 #This varies from 1 to 9 (http://unidata.github.io/netcdf4-python/)

#Define here the radius_of_influence for the regridding procedure
radius_of_influence=25e3 #In meters!!!!

#Define the number of x and y grid points
xgrdpoints=3600
ygrdpoints=1800

#xgrdpoints=3600
#ygrdpoints=1800

#Here defines the grid resolution which will be denoted in the output files
if (xgrdpoints==3600) & (ygrdpoints==1800):
    
    gridres='grid_resolution_0.1'
    lat_0=0.05
    lon_0=0.05
    
elif (xgrdpoints==360) & (ygrdpoints==180):
    
    gridres='grid_resolution_1.0'
    lat_0=0.0
    lon_0=0.0

#Which is the current directory
#cwd=os.getcwd()

#Define here where is the path with the MODIS level 2 data
#modislevel2path='/media/antonis/My Book/MERRA-2/MODIS-SWATHS-WITH-DOD-USING-MERRA2-DUST-RATIO/'+satellite+'/'+year+'/'
modislevel2path='/mnt/nas-1/MIDAS-Processing_Data/MODIS-SWATHS-WITH-DOD-USING-MERRA2-DUST-RATIO/'+satellite+'/'+year+'/'

#Create the folder where the outputs will be stored
#ncfilesfolder='/media/antonis/My Book/MERRA-2/REPROJECTED-MODIS-SWATHS-SYNERGY-WITH-MERRA2/'+satellite+'/'+year+'/'+gridres.upper()+'/'
#ncfilesfolder='/mnt/nas/datasets/MIDAS-Processing_Data/REPROJECTED-MODIS-SWATHS-SYNERGY-WITH-MERRA2/'+satellite+'/'+year+'/'+gridres.upper()+'/'
ncfilesfolder='/mnt/nas-1/MIDAS-Processing_Data/REPROJECTED-MODIS-SWATHS-SYNERGY-WITH-MERRA2/'+satellite+'/'+year+'/'+gridres.upper()+'/'

if not os.path.exists(ncfilesfolder):
    os.makedirs(ncfilesfolder)
    
###############################################################################################################
#Test the pyresample module
#modfilepath='/media/antonis/My Book/MERRA-2/MODIS-SWATHS-WITH-DOD-USING-MERRA2-DUST-RATIO/MODIS-AQUA/2015/'
#modfile='MYD04_L2.A2015032.1135.006.2015033173903-MODIS-SWATH-DOD.nc'

def read_MODIS_level2_swath(inpath,infile):

    modiswath=Dataset(os.path.join(inpath,infile),'r')

    #Extract the coordinates from the swath file
    modlat=modiswath['mod04']['Geolocation_Fields']['Latitude'][:]
    modlon=modiswath['mod04']['Geolocation_Fields']['Longitude'][:]

    modisAOD=modiswath['mod04']['Data_Fields']['AOD_550_Dark_Target_Deep_Blue_Combined'][:]
    modisTOTALDOD=modiswath['mod04']['Data_Fields']['MODIS-TOTAL-Dust-optical-depth-at-550nm'][:]
    
    modisFINEDOD=modiswath['mod04']['Data_Fields']['MODIS-FINE-Dust-optical-depth-at-550nm'][:]
    modisCOARSEDOD=modiswath['mod04']['Data_Fields']['MODIS-COARSE-Dust-optical-depth-at-550nm'][:]
    
    modisALG=modiswath['mod04']['Data_Fields']['AOD_550_Dark_Target_Deep_Blue_Combined_Algorithm_Flag'][:]
    
    modisLANDSEAFLAG=modiswath['mod04']['Data_Fields']['Land_sea_Flag'][:]
    
    modisANGOCEAN1=modiswath['mod04']['Data_Fields']['Angstrom_Exponent_1_Ocean'][:]
    modisANGOCEAN2=modiswath['mod04']['Data_Fields']['Angstrom_Exponent_2_Ocean'][:]
    modisANGLAND=modiswath['mod04']['Data_Fields']['Deep_Blue_Angstrom_Exponent_Land'][:]

    modisSOLZENANGL=modiswath['mod04']['Data_Fields']['Solar_Zenith'][:]
    modisSENZENANGL=modiswath['mod04']['Data_Fields']['Sensor_Zenith'][:]    
    
    
    modiswath.close()
    
    """
    0: Latitude
    1: Longitude
    2: MODIS combined AOD
    3: MODIS TOTAL dust AOD (DOD)
    4: MODIS FINE dust AOD (DOD)
    5: MODIS COARSE dust AOD (AOD)
    6: Which algorithm is applied (DT or DB)
    7: Land-Sea Flag
    8: MODIS Angstrom 1 over ocean (550-860nm)
    9: MODIS Angstrom 2 over ocean (860-2130nm)
    10: MODIS Angstrom over land (Deep Blue product)
    11: Solar zenith angle
    12: Sensor zenith angle
    
    """
    
    return modlat, modlon, modisAOD, modisTOTALDOD, modisFINEDOD, modisCOARSEDOD, modisALG, modisLANDSEAFLAG, modisANGOCEAN1, modisANGOCEAN2, modisANGLAND, modisSOLZENANGL, modisSENZENANGL   
###############################################################################################################
###############################################################################################################
#Extract from netcdf files the dates for each file
def find_MODIS_swath_overpass_time(inpath,infile):

    y=np.int(infile[10:14])
    jd=np.int(infile[14:17])
    hour=np.int(infile[18:20]) 
    minutes=np.int(infile[20:22]) 
    
    month=1
    
    while jd - calendar.monthrange(y,month)[1] > 0 and month <= 12:
        jd = jd - calendar.monthrange(y,month)[1]
        month += 1
    
    modisoverpasstime = dt.datetime(y,month,jd,hour,minutes).strftime("%d:%m:%Y %H:%M:%S")
    
    overpasstimemodis=pd.to_datetime(modisoverpasstime,format="%d:%m:%Y %H:%M:%S")
    
    return overpasstimemodis
###############################################################################################################
###############################################################################################################
#Here create the part of the script where the netcdf files from MODIS swaths will be created

#Read the MODIS Level 2 swaths
for modsubdirs,moddirs,modfiles in sorted(os.walk(modislevel2path)):
    for modfile in sorted(fnmatch.filter(modfiles,'*.nc')):
        #f.write(modfile + '\n')
        #f.flush()
        
        #Create the netcdf file inclduing lat,lot and AOD (filtered based on the function read_MODIS_level2_swath)
        modswathncfile=modfile[:26]+'-MODIS-SWATH-DOD'+'-'+gridres+'.nc'#
        if os.path.exists(os.path.join(ncfilesfolder,modswathncfile)):
            #f.write(modswathncfile + ' exists\n')
            #f.flush()
            continue
            
        try:
           
            modiswath=read_MODIS_level2_swath(modsubdirs,modfile)

            #Working with pyresample in order to regrid the data
            #Definition of the Global domain
            area_def = pr.geometry.AreaDefinition('Global_domain', 'Global_domain', 'Global_domain',
                                           {'lat_0': lat_0,'lon_0': lon_0, 'proj': 'eqc'},
                                           xgrdpoints, ygrdpoints,
                                          (-20037508.34, -10018754.17, 20037508.34, 10018754.17))
            
            #Extract here the lonlats from the defined global domain
            area_def_coords=area_def.get_lonlats(nprocs=nprocs)
            
            globlons=np.round(area_def_coords[0],2)
            globlats=np.round(area_def_coords[1],2)
            
            #Read the coordinates from the swath data    
            swath_def = pr.geometry.SwathDefinition(lons=modiswath[1], lats=modiswath[0])
            
            #Do the reprojection to the global map
            #Nearest for the combined AOD    
            
            resultnearmodisAOD = pr.kd_tree.resample_nearest(swath_def, modiswath[2],area_def,
                                                     fill_value=None, nprocs=nprocs,radius_of_influence=radius_of_influence, epsilon=0.)
            resultnearmodisDOD = pr.kd_tree.resample_nearest(swath_def, modiswath[3],area_def,
                                                     fill_value=None, nprocs=nprocs,radius_of_influence=radius_of_influence, epsilon=0.)
            resultnearmodisDODfine = pr.kd_tree.resample_nearest(swath_def, modiswath[4],area_def,
                                                     fill_value=None, nprocs=nprocs,radius_of_influence=radius_of_influence, epsilon=0.)
            resultnearmodisDODcoarse = pr.kd_tree.resample_nearest(swath_def, modiswath[5],area_def,
                                                     fill_value=None, nprocs=nprocs,radius_of_influence=radius_of_influence, epsilon=0.)
            resultnearmodisALG = pr.kd_tree.resample_nearest(swath_def, modiswath[6],area_def,
                                                     fill_value=None, nprocs=nprocs,radius_of_influence=radius_of_influence, epsilon=0.)
            resultnearmodisLANDSEAFLAG = pr.kd_tree.resample_nearest(swath_def, modiswath[7],area_def,
                                                     fill_value=None, nprocs=nprocs,radius_of_influence=radius_of_influence, epsilon=0.)
            resultnearmodisANGOCEAN1 = pr.kd_tree.resample_nearest(swath_def, modiswath[8],area_def,
                                                     fill_value=None, nprocs=nprocs,radius_of_influence=radius_of_influence, epsilon=0.)
            resultnearmodisANGOCEAN2 = pr.kd_tree.resample_nearest(swath_def, modiswath[9],area_def,
                                                     fill_value=None, nprocs=nprocs,radius_of_influence=radius_of_influence, epsilon=0.)
            resultnearmodisANGLAND = pr.kd_tree.resample_nearest(swath_def, modiswath[10],area_def,
                                                     fill_value=None, nprocs=nprocs,radius_of_influence=radius_of_influence, epsilon=0.)
            resultnearmodisSOLZENANGL = pr.kd_tree.resample_nearest(swath_def, modiswath[11],area_def,
                                                     fill_value=None, nprocs=nprocs,radius_of_influence=radius_of_influence, epsilon=0.)
            resultnearmodisSENZENANGL = pr.kd_tree.resample_nearest(swath_def, modiswath[12],area_def,
                                                     fill_value=None, nprocs=nprocs,radius_of_influence=radius_of_influence, epsilon=0.)
            
            ncmodis=Dataset(os.path.join(ncfilesfolder,modswathncfile),'w',format='NETCDF4')

            #Create Global Attributes
            ncmodis.description=satellite+'-'+'FILTERED-RAW-SWATHS'
            ncmodis.history='Created:  '+time.ctime(time.time())
            ncmodis.source='netCDF4 python module tutorial'
            ncmodis.Conventions = 'CF-1.6'

            #Thresholds which have been used
            #ncmodis.AODocTHRE = str(aodoceanthre)
            #ncmodis.ANGlandTHRE = str(anglandthre)
            #ncmodis.ANGocTHRE = str(angoceanthre)
            #ncmodis.CFaodTHRE = str(aodcfthre)
            
            #Create groups
            maingrp=ncmodis.createGroup('mod04')
            datafieldsgrp=ncmodis.createGroup('/mod04/Data_Fields')
            geolocgrp=ncmodis.createGroup('/mod04/Geolocation_Fields')
            
            
            latitude=ncmodis.createDimension('Latitude',globlats.shape[0])
            longitude=ncmodis.createDimension('Longitude',globlons.shape[1])
            
            #Create the variables which will be stored in the netcdf file
            #Latitude
            Latitude=ncmodis.createVariable('/mod04/Geolocation_Fields/Latitude',np.float32,('Latitude','Longitude'),
                                           fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            Latitude[:]=globlats
            
            Latitude.valid_range=(-90.,90.)
            Latitude.long_name='Geodetic Latitude'
            Latitude.units='degrees_north'
            
            #Longitude
            Longitude=ncmodis.createVariable('/mod04/Geolocation_Fields/Longitude',np.float32,('Latitude','Longitude'),
                                           fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            Longitude[:]=globlons
            
            Longitude.valid_range=(-180.,180.)
            Longitude.long_name='Geodetic Longitude'
            Longitude.units='degrees_east'
            
            #Total AOD
            modaod=ncmodis.createVariable('mod04/Data_Fields/'+'AOD_550_Dark_Target_Deep_Blue_Combined',np.float32,('Latitude','Longitude'),
                                             fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            modaod[:]=resultnearmodisAOD
            
            modaod.long_name="Combined Dark Target, Deep Blue AOT at 0.55 micron for land and ocean."
            modaod.units='None'
            modaod.coordinates='Latitude Longitude'
            
            
            #Total DOD
            moddod=ncmodis.createVariable('mod04/Data_Fields/'+'MODIS-TOTAL-Dust-optical-depth-at-550nm',np.float32,('Latitude','Longitude'),
                                             fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            moddod[:]=resultnearmodisDOD
            
            moddod.long_name="Total Dust Optical Depth at 550nm through MERRA-2 dust-to-total AOD ratio"
            moddod.units='None'
            modaod.coordinates='Latitude Longitude'

            #Fine DOD
            modFINEdod=ncmodis.createVariable('mod04/Data_Fields/'+'MODIS-FINE-Dust-optical-depth-at-550nm',np.float32,('Latitude','Longitude'),
                                             fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            modFINEdod[:]=resultnearmodisDODfine
            
            modFINEdod.long_name="Fine Dust Optical Depth at 550nm through MERRA-2 dust-to-total AOD ratio"
            modFINEdod.units='None'
            modFINEdod.coordinates='Latitude Longitude'
            
            #Coarse DOD
            modCOARSEdod=ncmodis.createVariable('mod04/Data_Fields/'+'MODIS-COARSE-Dust-optical-depth-at-550nm',np.float32,('Latitude','Longitude'),
                                             fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            modCOARSEdod[:]=resultnearmodisDODcoarse
            
            modCOARSEdod.long_name="Coarse Dust Optical Depth at 550nm through MERRA-2 dust-to-total AOD ratio"
            modCOARSEdod.units='None'
            modCOARSEdod.coordinates='Latitude Longitude'
            
            #MODIS Algorithm applied
            modALGORITHM=ncmodis.createVariable('mod04/Data_Fields/'+'AOD_550_Dark_Target_Deep_Blue_Combined_Algorithm_Flag',np.float32,('Latitude','Longitude'),
                                             fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            modALGORITHM[:]=resultnearmodisALG
            
            modALGORITHM.long_name="Combined Dark Target, Deep Blue AOT at 0.55 micron Algorithm Flag (0=Dark Target, 1=Deep Blue, 2=Mixed)"
            modALGORITHM.units='None'
            modALGORITHM.coordinates='Latitude Longitude'
            
            #MODIS Land-Sea mask
            modLANDSEAFLAG=ncmodis.createVariable('mod04/Data_Fields/'+'Land_sea_Flag',np.float32,('Latitude','Longitude'),
                                             fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            modLANDSEAFLAG[:]=resultnearmodisLANDSEAFLAG
            
            modLANDSEAFLAG.long_name="Land_sea_Flag(based on MOD03 Landsea mask 0 = Ocean, 1 = Land and Ephemeral water 2 =Coastal)"
            modLANDSEAFLAG.units='None'
            modLANDSEAFLAG.coordinates='Latitude Longitude'
            
            #Angstrom Ocean 1
            modANGOCEAN1=ncmodis.createVariable('mod04/Data_Fields/'+'Angstrom_Exponent_1_Ocean',np.float32,('Latitude','Longitude'),
                                             fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            modANGOCEAN1[:]=resultnearmodisANGOCEAN1
            
            modANGOCEAN1.long_name="Calculated Angstrom Exponent for 0.55 vs 0.86 micron  for Average Solution"
            modANGOCEAN1.units='None'
            modANGOCEAN1.coordinates='Latitude Longitude'

            #Angstrom Ocean 2
            modANGOCEAN2=ncmodis.createVariable('mod04/Data_Fields/'+'Angstrom_Exponent_2_Ocean',np.float32,('Latitude','Longitude'),
                                             fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            modANGOCEAN2[:]=resultnearmodisANGOCEAN2
            
            modANGOCEAN2.long_name="Calculated Angstrom Exponent for 0.86 vs 2.13 micron for Average Solution"
            modANGOCEAN2.units='None'
            modANGOCEAN2.coordinates='Latitude Longitude'
            
            
            #Angstrom Land (Deep BLue)
            modANGLAND=ncmodis.createVariable('mod04/Data_Fields/'+'Deep_Blue_Angstrom_Exponent_Land',np.float32,('Latitude','Longitude'),
                                             fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            modANGLAND[:]=resultnearmodisANGLAND
            
            modANGLAND.long_name="Deep Blue Angstrom Exponent for land with all quality data (Quality flag=1,2,3)"
            modANGLAND.units='None'
            modANGLAND.coordinates='Latitude Longitude'
            
            #Solar zenith angle
            modSOLZENANGL=ncmodis.createVariable('mod04/Data_Fields/'+'Solar_Zenith',np.float32,('Latitude','Longitude'),
                                             fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            modSOLZENANGL[:]=resultnearmodisSOLZENANGL
            
            modSOLZENANGL.long_name="Solar Zenith Angle, Cell to Sun"
            modSOLZENANGL.units='Degrees'
            modSOLZENANGL.valid_range=(0.,180.)
            modSOLZENANGL.coordinates='Latitude Longitude'
            
            #Sensor zenith angle
            modSENZENANGL=ncmodis.createVariable('mod04/Data_Fields/'+'Sensor_Zenith',np.float32,('Latitude','Longitude'),
                                             fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            modSENZENANGL[:]=resultnearmodisSENZENANGL
            
            modSENZENANGL.long_name="Sensor Zenith Angle, Cell to Sensor"
            modSENZENANGL.units='Degrees'
            modSENZENANGL.valid_range=(0.,180.)
            modSENZENANGL.coordinates='Latitude Longitude'
            
            
            ncmodis.close()
                    
        except:
            
            print ("I cannot read: ",modfile)
            #f.write(modfile + " error/n")
            #f.flush()



#m=Basemap(projection='cyl',llcrnrlat=15,urcrnrlat=45,
#          llcrnrlon=5,urcrnrlon=40,resolution='l')
#
#m.drawcoastlines(linewidth=0.25)#,zorder=2)
#m.drawcountries(linewidth=0.25)#,zorder=2)
#m.fillcontinents(color='Gainsboro',lake_color='Gainsboro',zorder=0)
#m.drawparallels(np.linspace(15,45,7),linewidth=0.25,labels=[1,1,1,1],fmt='%.0f',fontweight='bold',fontsize=8,zorder=3)
#m.drawmeridians(np.linspace(5,40,8),linewidth=0.25,labels=[1,1,1,1],fmt='%.0f',fontweight='bold',fontsize=8,zorder=3)
#m.drawmapboundary(fill_color='Gainsboro',linewidth=0.5)
#
##m.contourf(globlons,globlats,resultnear)
#m.pcolormesh(globlons,globlats,resultnear)


scriptend=dt.datetime.now()
#f.close()

print ("")
print ("")
print ("")
print ("The script ended at:", scriptend)
print ("")
print ("")
print ("")
print ("The execution time was:", scriptend-scriptstart)
print ("")
print ("")
print ("")
