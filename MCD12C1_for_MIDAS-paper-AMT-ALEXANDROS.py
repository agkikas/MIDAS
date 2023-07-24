# %%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 11:56:36 2020

@author: antonis
"""

from pyhdf.SD import SD,SDC
from netCDF4 import Dataset
#from mpl_toolkits.basemap import Basemap
#from matplotlib.patches import Polygon
#from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
#from math import ceil
from scipy import spatial


#import pandas as pd
import datetime as dt
import os
#import fnmatch
#import glob2
import numpy as np
import numpy.ma as ma
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#import cmocean
import time
#import matplotlib.dates as dates
#import matplotlib.ticker
#import calendar

scriptstart=dt.datetime.now()

print ("The script started at:", scriptstart)

# %%
###################################################################################################################################
#Define here what do you need for the path on the raw files
satcol='C061'
satellite='MODIS-AQUA'
satres='GRID_RESOLUTION_0.1'
year = '2020'

if satres=='GRID_RESOLUTION_0.1':
    
    nrow=1800
    ncol=3600
    
    latstep=0.1
    lonpixels=3600
    
elif satres=='GRID_RESOLUTION_1.0':
    
    nrow=180
    ncol=360
    
    latstep=1.0
    lonpixels=360
   
#################################################################################################    
#TODO
#Check always the yearfirst and yearlast!!!!!!
if satellite=="MODIS-AQUA":
    
    yearfirst=2003
    yearlast=2020
    
    years=np.arange(yearfirst,yearlast+1,1)
    
elif satellite=="MODIS-TERRA":

    yearfirst=2001
    yearlast=2020
    
    years=np.arange(yearfirst,yearlast+1,1)    
#%%
################################################################################################# 
#Define here where is the path with the raw files
#inpath='/media/antonis/My Passport/MERRA-2/DAILY-MODIS-LEVEL2-0.1RESOLUTION/'+satcol+'/'+satellite+'/'
inpath='/mnt/nas/datasets/MIDAS-Processing_Data/DAILY-MODIS-FILES-THROUGH-SYNERGY-OF-MODIS-MERRA2/'+satellite+'/'+year+'/'+satres

#Create the output folder
resfolder='/mnt/nas/datasets/MCD12C1/MCD12C1_DATA_FOR_MIDAS/'+year+ '/' #result

if not os.path.exists(resfolder):
    os.makedirs(resfolder)

#List all MODIS files according to the filepattern
#modisfiles = sorted(glob2.glob(inpath+'/**/'+'*'+satres+'*'+'.nc'))

#startdate=dt.datetime.strptime(modisfiles[0][-11:-3],'%Y%m%d')
#enddate=dt.datetime.strptime(modisfiles[-1][-11:-3],'%Y%m%d')

#Create the timeseries
#studyperiod=pd.date_range(startdate,enddate,freq='D')

# winterdays=studyperiod[(studyperiod.month==1) | (studyperiod.month==2) | (studyperiod.month==12)]
# springdays=studyperiod[(studyperiod.month==3) | (studyperiod.month==4) | (studyperiod.month==5)]
# summerdays=studyperiod[(studyperiod.month==6) | (studyperiod.month==7) | (studyperiod.month==8)]
# autumndays=studyperiod[(studyperiod.month==9) | (studyperiod.month==10) | (studyperiod.month==11)]

# alldays=len(studyperiod)
# djfdays=len(winterdays)
# mamdays=len(springdays)
# jjadays=len(summerdays)
# sondays=len(autumndays)

#Define the least significant digits and compression level for the variables which will be stored in the netcdf files
lstsigdig=10
complevel=9 #This varies from 1 to 9 (http://unidata.github.io/netcdf4-python/)

#%%
# # =============================================================================
# # Read the climatological file with dods as well as the coordinates
# # =============================================================================
modisclimpath='/mnt/nas/datasets/MCD12C1/MODIS-CLIMATOLOGICAL-RESULTS-DOD/MAPS/MODIS-AQUA/GRID_RESOLUTION_0.1' #geographical annual
modisclimfile='Geographical_annual_and_seasonal_distributions_AOD_and_DOD_MODIS-AQUA_GRID_RESOLUTION_0.1-01_Jan_2003-31_Dec_2017.nc'

modclim=Dataset(os.path.join(modisclimpath,modisclimfile),'r')

# dodann=modclim["Annual_mean_DOD@550nm"][:] #Annual DOD
# doddjf=modclim["DJF_mean_DOD@550nm"][:] #DJF DOD
# dodmam=modclim["MAM_mean_DOD@550nm"][:] #MAM DOD
# dodjja=modclim["JJA_mean_DOD@550nm"][:] #JJA DOD
# dodson=modclim["SON_mean_DOD@550nm"][:] #SON DOD

modlat=modclim["Latitude"][:] #Latitude
modlon=modclim["Longitude"][:] #Longitude

modclim.close()

#%%
# =============================================================================
# This function is used in order to get IGBP land cover classification
# =============================================================================
def extract_indices_from_MCD12C1 (mcdfilepath,mcdfilename):
    
    #Read the MCD12C1 file for the corresponding year
    modltfile=SD(os.path.join(mcdfilepath,mcdfilename),SDC.READ)
    
    #Create the coordinates for the MCD12C1
    modltlat=np.linspace(90.,-89.95,3600)
    modltlon=np.linspace(-180.,179.95,7200)
    
    #Read the Land Cover Type 1 (Annual IGBP classification)
    sdsobj=modltfile.select("Majority_Land_Cover_Type_1")
    
    fill_value=sdsobj.attributes()['_FillValue']
    
    modlt=sdsobj.get()
    modltfin=ma.masked_equal(modlt,fill_value)
    
    #Mesh the LAND type coordinates
    xx,yy=np.meshgrid(modltlon,modltlat)
    
    #Mesh the MIDAS coordinates
    xxmid,yymid=np.meshgrid(modlon[0,:],modlat[:,0])
    
    #Land type coordinates ravel()
    modltcoord=list(zip(xx.ravel(),yy.ravel()))
    
    #MIDAS coordinates ravel()
    midcoord=list(zip(xxmid.ravel(),yymid.ravel()))
    
    #Find the indices with the minimum distances of midcoord to modltcoord
    tree = spatial.cKDTree(modltcoord)
    
    mindist, minid = tree.query(midcoord)
    
    distids=modltfin.ravel()[minid]
    
    dist2darr=np.reshape(distids,(modlat[:,0].shape[0],modlon[0,:].shape[0]))
    
    modltfile.end()
    
    return dist2darr
#%%
# =============================================================================
# Extract the annual MCD12C1 files for each year
# =============================================================================
#Find the year of the study period
#years=np.sort(np.unique(studyperiod.year))

#Define here where the MCD12C1 data are stored
modltpath='/mnt/nas/datasets/MCD12C1/' +year+ '/' #downloaded

#%%

#Work with each annual file
#for year in years: #!original code
try:#!added
    
    modltfilename='MCD12C1.A'+str(year)+'001.006.2021362215328.hdf' #!check for other runs
    
    print (modltfilename)

    landcover=extract_indices_from_MCD12C1(modltpath,modltfilename)

    # =============================================================================
    # Create the output netcdf file 
    # =============================================================================
    outncfile='MIDAS-'+modltfilename[:-18]+'.nc4'
    
    ncfile=Dataset(os.path.join(resfolder,outncfile),'w',format='NETCDF4')
    
    #Create Global Attributes
    ncfile.description='MCD12C1_FOR_MIDAS'
    ncfile.history='Created:  '+time.ctime(time.time())
    ncfile.source='netCDF4 python module tutorial'
    ncfile.Conventions = 'CF-1.6'
    
    ncfilelat=ncfile.createDimension('Latitude',landcover.shape[0])
    ncfilelon=ncfile.createDimension('Longitude',landcover.shape[1])
    
    #Latitude
    Latitude=ncfile.createVariable('Latitude',np.float32,('Latitude','Longitude'),
                                    fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
        
    Latitude[:]=modlat
    
    Latitude.valid_range=(-90.,90.)
    Latitude.long_name='Geodetic Latitude'
    Latitude.units='degrees_north'
        
    #Longitude
    Longitude=ncfile.createVariable('Longitude',np.float32,('Latitude','Longitude'),
                                       fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel) 
    
    Longitude[:]=modlon
            
    Longitude.valid_range=(-180.,180.)
    Longitude.long_name='Geodetic Longitude'
    Longitude.units='degrees_east'
    
    #Land cover (IGBP)
    Countsann=ncfile.createVariable('IGBP_Land_Cover',np.float32,('Latitude','Longitude'),
                                       fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
    
    Countsann[:]=landcover
    Countsann.long_name='IGBP_most_likely_classification_scheme'
    Countsann.units='None'
    
    #Close the netcdf file
    ncfile.close()
except:
    print('Error in' + year)

###################################################################################################################################

scriptend=dt.datetime.now()

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
# %%
