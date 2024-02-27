#%%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 10:02:06 2020

@author: antonis
"""
#%%
from pyhdf.SD import SD,SDC
from netCDF4 import Dataset
#from mpl_toolkits.basemap import Basemap
#from matplotlib.patches import Polygon
#from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
#from math import ceil
from scipy import spatial

import pandas as pd
import datetime as dt
import os
#import fnmatch
import glob2
import numpy as np
import numpy.ma as ma
import matplotlib as mpl
import matplotlib.pyplot as plt
#import cmocean
import time
#import matplotlib.dates as dates
import matplotlib.ticker
#import calendar
import netCDF4 as ncf

scriptstart=dt.datetime.now()

print ("The script started at:", scriptstart)

###################################################################################################################################
#Define here what do you need for the path on the raw files
satcol='C061'
satellite='MODIS-TERRA'
satres='GRID_RESOLUTION_0.1'
satyear='2022'

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

#%%

#Define here where is the path with the raw files
#inpath='/media/antonis/My Passport/MERRA-2/DAILY-MODIS-LEVEL2-0.1RESOLUTION/'+satcol+'/'+satellite+'/'
inpath='/mnt/nas-1/MIDAS-Processing_Data/DAILY-MODIS-FILES-THROUGH-SYNERGY-OF-MODIS-MERRA2/'+satellite+'/'+satyear+'/'+satres

#Define here the folder where the results will be stored
# resfolder='/media/antonis/My Passport/MERRA-2/'+satcol+'/MEASUREMENT_SAMPLING_TOTAL_UNCERTAINTY/'+satellite+'/'+satres+'/'

# if not os.path.exists(resfolder):
#     os.makedirs(resfolder)
    
#Define here the folder where the DOD uncertainty netcdf files will be stored
#doduncfolder='/media/antonis/My Passport/MERRA-2/'+satcol+'/DOD_UNCERTAINTY_NETCDF_DAILY_FILES/'
doduncfolder='/mnt/nas-1/MCD12C1/DOD_UNCERTAINTY_NETCDF_DAILY_FILES/'

if not os.path.exists(doduncfolder):
    os.makedirs(doduncfolder)

#List all MODIS files according to the filepattern
modisfiles = sorted(glob2.glob(inpath+'/**/'+'*'+satres+'*'+'.nc'))

startdate=dt.datetime.strptime(modisfiles[0][-11:-3],'%Y%m%d')
enddate=dt.datetime.strptime(modisfiles[-1][-11:-3],'%Y%m%d')

#Create the timeseries
studyperiod=pd.date_range(startdate,enddate,freq='D')

winterdays=studyperiod[(studyperiod.month==1) | (studyperiod.month==2) | (studyperiod.month==12)]
springdays=studyperiod[(studyperiod.month==3) | (studyperiod.month==4) | (studyperiod.month==5)]
summerdays=studyperiod[(studyperiod.month==6) | (studyperiod.month==7) | (studyperiod.month==8)]
autumndays=studyperiod[(studyperiod.month==9) | (studyperiod.month==10) | (studyperiod.month==11)]

alldays=len(studyperiod)
djfdays=len(winterdays)
mamdays=len(springdays)
jjadays=len(summerdays)
sondays=len(autumndays)

#Define the least significant digits and compression level for the variables which will be stored in the netcdf files
lstsigdig=10
complevel=9 #This varies from 1 to 9 (http://unidata.github.io/netcdf4-python/)

# =============================================================================
# Map features
# =============================================================================
#Define the borders of your study domain
southlimit=-90.
northlimit=90.
westlimit=-180.
eastlimit=180.

#Define map features
gridres=0.5
gridlw=0.25
landcolor='Gainsboro'
lake_color='Gainsboro'
boundarycolor='Gainsboro'

gridfontsize=8.

#Define the steps for the parallels and the meridians
steppar=30.
stepmer=45.

temin=0.
temax=0.1
teticks=11
tecbarticks=np.linspace(temin,temax,teticks)

tecmap='gist_heat_r'

reltemin=0.
reltemax=100.
relteticks=11
reltecbarticks=np.linspace(reltemin,reltemax,relteticks)

reltecmap='plasma'

meanuncmin=0.
meanuncmax=1.
meanuncticks=11
meanunccbarticks=np.linspace(meanuncmin,meanuncmax,meanuncticks)

meanunccmap='gist_heat_r'

#Define parallels and meridians for the plots
parallels=np.arange(southlimit,northlimit+steppar,steppar)
meridians=np.arange(westlimit,eastlimit+stepmer,stepmer)
# =============================================================================
# Calculate the surface area of each grid cell
# =============================================================================
#Define here the south and north limit of the domain 
nlim=90
slim=-90

r=6371 #In kilometers!!!

lats=np.arange(nlim,slim-latstep,-1.*latstep)

#Create the array with surface areas for each grid cell
mersurf=[]

for nlat,slat in zip(lats[:-1],lats[1:]):
    
    A=2*np.pi*r**2*np.abs(np.sin(np.deg2rad(nlat))-np.sin(np.deg2rad(slat)))
    
    mersurfband=A/lonpixels
    
    mersurf.append(mersurfband)

gridsurfacetmp=np.array(mersurf)    
    
gridsurfacetmp2=np.dstack([gridsurfacetmp]*lonpixels)

gridsurface=gridsurfacetmp2.squeeze()

# =============================================================================
# Read the climatological file with dods as well as the coordinates
# =============================================================================
modisclimpath='/mnt/nas-1/MCD12C1/MODIS-CLIMATOLOGICAL-RESULTS-DOD/MAPS/MODIS-AQUA/GRID_RESOLUTION_0.1'
modisclimfile='Geographical_annual_and_seasonal_distributions_AOD_and_DOD_MODIS-AQUA_GRID_RESOLUTION_0.1-01_Jan_2003-31_Dec_2017.nc'

modclim=Dataset(os.path.join(modisclimpath,modisclimfile),'r')

dodann=modclim["Annual_mean_DOD@550nm"][:] #Annual DOD
doddjf=modclim["DJF_mean_DOD@550nm"][:] #DJF DOD
dodmam=modclim["MAM_mean_DOD@550nm"][:] #MAM DOD
dodjja=modclim["JJA_mean_DOD@550nm"][:] #JJA DOD
dodson=modclim["SON_mean_DOD@550nm"][:] #SON DOD

modlat=modclim["Latitude"][:] #Latitude
modlon=modclim["Longitude"][:] #Longitude

modclim.close()

# =============================================================================
# Set here the equations for AOD and MDF uncertainties
# =============================================================================
#Dark target land
dtlandslope=0.15
dtlandintercept=0.05

#Dark target ocean
dtoceanslope=0.1
dtoceanintercept=0.04

#Deep blue land (Collection 6)
#dblandslope=0.56
#dblandintercept=0.086

#Deep blue land (Collection 6.1)
dblandslopearid=0.61
dblandinterceptarid=0.12

dblandslopeveg=0.67
dblandinterceptveg=0.079

##Provide here the coefficients for the merged DT and DB AOD (Sayer et al., 2014)
##Don't forget to use the square in the formula !!!!!
#mxcoef1=0.084
#mxcoef2=0.0409
#mxcoef3=0.043

#Here is the error for the MDF
#Define here if you are using the CALIOP profiles of the DOD thresholds
#merthre='CALIOP-PROFILES'
merthre='DOD-THRESHOLD'

if merthre=='CALIOP-PROFILES':

    #Provide here the coefficients from the polynomial for MERRA-2 dust fraction
    #You have to check the plots stored in /media/antonis/My Passport/MERRA-2/UNCERTAINTY-DUST-FRACTION-CALIOP-MERRA2
    #Check if you are working with DAY, NIGHT or DAY+NIGHT
    #Check also which number of CALIOP profiles is used as threshold!!!
    fraccoef1=5.254
    fraccoef2=-9.777
    fraccoef3=4.410
    fraccoef4=0.176
    fraccoef5=0.047
        
elif merthre=='DOD-THRESHOLD': 

    #Provide here the coefficients from the polynomial for MERRA-2 dust fraction
    #You have to check the plots stored in /media/antonis/My Passport/MERRA-2/UNCERTAINTY-DUST-FRACTION-CALIOP-MERRA2
    #Check if you are working with DAY, NIGHT or DAY+NIGHT
    #Check also which is the DOD threshold
    fraccoef1=2.282
    fraccoef2=-6.222
    fraccoef3=4.700
    fraccoef4=-0.969
    fraccoef5=0.199
    
# =============================================================================
# Define here the border of the domains
# =============================================================================
#Give here the borders for each region
# westborders=[-180,-180,-180,-45,-90,-11,-13,-17,-13,15,36,36,50,68,76,90,110,135,-180,-124,-73,12,122]
# northborders=[90,90,0,30,30,47,7,30,13,19,37,27,50,30,42,45,47,60,60,43,-35,-19,-18]
# eastborders=[180,180,180,-15,-45,39,9,9,9,19,49,60,65,76,90,110,130,180,-135,-103,-57,25,150]
# southborders=[-90,0,-90,0,10,29,0,13,7,15,27,12,37,24,36,37,30,30,30,32,-56,-30,-30]

# regionnames=['Global','Northern Hemisphere','Southern Hemisphere',
#              'East Tropical Atlantic','West Tropical Atlantic','Mediterranean',
#              'Gulf of Guinea','West Sahara','Sub-Sahel',
#              'Bodele','North Middle East','South Middle East',
#              'Central Asia','Thar','Taklamakan',
#              'Gobi','East Asia','West North Pacific',
#              'East North Pacific','Southwest US','Patagonia',
#              'Kalahari','Australia']

# regacronyms=['GLB','NHE','SHE',
#              'ETA','WTA','MED',
#              'GOG','WSA','SSA',
#              'BOD','NME','SME',
#              'CAS','THA','TAK',
#              'GOB','EAS','WNP',
#              'ENP','SUS','PAT',
#              'KAL','AUS']

## =============================================================================
## This function is used in order to derive the square of the DOD uncertainty
## =============================================================================
##def dod_uncertainty_square(infilepath,infile):
#def dod_modis_outputs(infilepath,infile):
#
#    ncmodis=Dataset(os.path.join(infilepath,infile),'r')
#    
#    modlat=ncmodis['Latitude'][:] #Latitude
#    modlon=ncmodis['Longitude'][:] #Longitude
#    
#    modlsf=ncmodis['Lsf'][:] #Land-Ocean flag
#    
#    modaod=ncmodis['Aod_550_dark_target_deep_blue_combined'][:] #Total AOD
#    moddod=ncmodis['Modis-total-dust-optical-depth-at-550nm'][:] #Total DOD
#    
#    moddfrac=np.abs(moddod/modaod) #DOD-to-AOD ratio
#    
#    modsenzen=np.deg2rad(ncmodis['Sensor_zenith'][:]) #Sensor zenith angle (converted to radians)
#    modsolzen=np.deg2rad(ncmodis['Sza5'][:]) #Solar zentih angle (converted to radians)
#    
#    modamf=(1./np.cos(modsenzen))+(1./np.cos(modsolzen)) #Calculation of the air mass factor
#    
#    modalgflag=ncmodis['Aod_550_dark_target_deep_blue_combined_algorithm_flag'][:] #MODIS Algorithm flag
#    
#    ncmodis.close()
#    
#    #Create nan arrays for uncertainty (each algorithm)
#    dtlandun=np.zeros((modlat.shape[0],modlon.shape[1]))*np.nan
#    dtoceanun=np.zeros((modlat.shape[0],modlon.shape[1]))*np.nan
#    dblandun=np.zeros((modlat.shape[0],modlon.shape[1]))*np.nan
#    dtdblandun=np.zeros((modlat.shape[0],modlon.shape[1]))*np.nan
#    
#    #Calculate the uncertainties for each algorithm
#    dtlandun=np.where((modlsf == 1) & (modalgflag == 0),(dtlandslope*modaod+dtlandintercept)*moddfrac+modaod*(fraccoef1*(moddfrac**4)+fraccoef2*(moddfrac**3)+fraccoef3*(moddfrac**2)+fraccoef4*moddfrac+fraccoef5),dtlandun)
#    dtoceanun=np.where((modlsf == 0) & (modalgflag == 0),(dtoceanslope*modaod+dtoceanintercept)*moddfrac+modaod*(fraccoef1*(moddfrac**4)+fraccoef2*(moddfrac**3)+fraccoef3*(moddfrac**2)+fraccoef4*moddfrac+fraccoef5),dtoceanun)
#    dblandun=np.where((modlsf == 1) & (modalgflag == 1),((dblandslope*modaod+dblandintercept)/modamf)*moddfrac+modaod*(fraccoef1*(moddfrac**4)+fraccoef2*(moddfrac**3)+fraccoef3*(moddfrac**2)+fraccoef4*moddfrac+fraccoef5),dblandun)
#    #C006
#    #dtdblandun=np.where((modlsf == 1) & (modalgflag == 2),np.sqrt((0.0225+(0.3136/modamf**2))*modaod**2 + (0.015+(0.09632/modamf**2))*modaod + 0.0025 + 0.007396/modamf**2)*moddfrac+modaod*(fraccoef1*(moddfrac**4)+fraccoef2*(moddfrac**3)+fraccoef3*(moddfrac**2)+fraccoef4*moddfrac+fraccoef5),dtdblandun)
#    
#    #C061
#    dtdblandun=np.where((modlsf == 1) & (modalgflag == 2),(np.sqrt((0.0225+(0.3136/modamf**2))*modaod**2 + (0.015+(0.09632/modamf**2))*modaod + 0.0025 + 0.007396/modamf**2)/2.)*moddfrac+modaod*(fraccoef1*(moddfrac**4)+fraccoef2*(moddfrac**3)+fraccoef3*(moddfrac**2)+fraccoef4*moddfrac+fraccoef5),dtdblandun)
#        
#    #Merge the four arrays corresponding to each MODIS algorithm
#    dodunctmp=np.nansum(np.dstack((dtlandun,dtoceanun,dblandun,dtdblandun)),2)
#    
#    #Replace zeros with nans in dodunctmp array
#    #dodunctmp1=np.where(modaod.mask==True,np.nan,dodunctmp)
#    dodunctmp1=np.where(modaod.mask==True,0.,dodunctmp)
#        
#    #Final dod uncertainty
#    dodunctmp2=ma.masked_invalid(dodunctmp1)
#    
##    #Reverse the latitudes
##    dodunc=dodunctmp2[::-1,:]
#    
#    #Calculate the square of the dodunc which will be the output of the function
#    #sqrdodunc=dodunc**2
#    sqrdodunc=dodunctmp2**2
#    
#    '''
#    Returns the: (1) square of the DOD uncertainty
#    (2) MODIS DOD
#    (3) DOD uncertainty
#    
#    '''
#    return sqrdodunc,moddod,dodunctmp2

####TEST
#infilepath='/media/antonis/My Passport/MERRA-2/DAILY-MODIS-LEVEL2-0.1RESOLUTION/C061/MODIS-AQUA/2017/GRID_RESOLUTION_0.1/'
#infile='MODIS-AQUA_AOD-and-DOD-GRID_RESOLUTION_0.1-20171231.nc'    
#%%
# =============================================================================
# This function is used in order to derive the square of the DOD uncertainty
# =============================================================================
def dod_modis_outputs(infilepath,infile):
        
    #Find the year of the file
    yearfile=infile[-11:-7]
    
    #Read the mcd12c1file
    #mcdfilepath='/media/antonis/My Passport/MERRA-2/C061/MCD12C1_DATA_FOR_MIDAS/MIDAS-MCD12C1.A'+yearfile+'001.006.nc4'
    mcdfilepath='/mnt/nas-1/MCD12C1/MCD12C1_DATA_FOR_MIDAS/'+satyear+'/'+'MIDAS-MCD12C1.A'+yearfile+'001.006.nc4'
       
    mcdfile=Dataset(mcdfilepath,'r')
    
    mcdvar=mcdfile['IGBP_Land_Cover'][:]
    
    mcdfile.close()
    
    #Work with the MODIS file
    ncmodis=Dataset(os.path.join(infilepath,infile),'r')
    
    modlat=ncmodis['Latitude'][:] #Latitude
    modlon=ncmodis['Longitude'][:] #Longitude
    
    modlsf=ncmodis['Lsf'][:] #Land-Ocean flag
    
    modaod=ncmodis['Aod_550_dark_target_deep_blue_combined'][:] #Total AOD
    moddod=ncmodis['Modis-total-dust-optical-depth-at-550nm'][:] #Total DOD
    
    moddfrac=np.abs(moddod/modaod) #DOD-to-AOD ratio
    
    modsenzen=np.deg2rad(ncmodis['Sensor_zenith'][:]) #Sensor zenith angle (converted to radians)
    modsolzen=np.deg2rad(ncmodis['Sza5'][:]) #Solar zentih angle (converted to radians)
    
    modamf=(1./np.cos(modsenzen))+(1./np.cos(modsolzen)) #Calculation of the air mass factor
    
    modalgflag=ncmodis['Aod_550_dark_target_deep_blue_combined_algorithm_flag'][:] #MODIS Algorithm flag
    
    ncmodis.close()
    
    #Create nan arrays for uncertainty (each algorithm)
    dtlandun=np.zeros((modlat.shape[0],modlon.shape[1]))*np.nan
    dtoceanun=np.zeros((modlat.shape[0],modlon.shape[1]))*np.nan
    dblandunarid=np.zeros((modlat.shape[0],modlon.shape[1]))*np.nan
    dblandunveg=np.zeros((modlat.shape[0],modlon.shape[1]))*np.nan
    dtdblandunarid=np.zeros((modlat.shape[0],modlon.shape[1]))*np.nan
    dtdblandunveg=np.zeros((modlat.shape[0],modlon.shape[1]))*np.nan
    
    #Calculate the uncertainties for each algorithm
    dtlandun=np.where((modlsf == 1) & (modalgflag == 0),(dtlandslope*modaod+dtlandintercept)*moddfrac+modaod*(fraccoef1*(moddfrac**4)+fraccoef2*(moddfrac**3)+fraccoef3*(moddfrac**2)+fraccoef4*moddfrac+fraccoef5),dtlandun)
    dtoceanun=np.where((modlsf == 0) & (modalgflag == 0),(dtoceanslope*modaod+dtoceanintercept)*moddfrac+modaod*(fraccoef1*(moddfrac**4)+fraccoef2*(moddfrac**3)+fraccoef3*(moddfrac**2)+fraccoef4*moddfrac+fraccoef5),dtoceanun)
    
    #DB for arid regions
    dblandunarid=np.where((modlsf == 1) & (modalgflag == 1) & (mcdvar == 16),((dblandslopearid*modaod+dblandinterceptarid)/modamf)*moddfrac+modaod*(fraccoef1*(moddfrac**4)+fraccoef2*(moddfrac**3)+fraccoef3*(moddfrac**2)+fraccoef4*moddfrac+fraccoef5),dblandunarid)
    
    #DB for vegetated areas
    dblandunveg=np.where((modlsf == 1) & (modalgflag == 1) & (mcdvar != 16) & (mcdvar != 0),((dblandslopeveg*modaod+dblandinterceptveg)/modamf)*moddfrac+modaod*(fraccoef1*(moddfrac**4)+fraccoef2*(moddfrac**3)+fraccoef3*(moddfrac**2)+fraccoef4*moddfrac+fraccoef5),dblandunveg)
    
    #C061 arid pixels
    dtdblandunarid=np.where((modlsf == 1) & (modalgflag == 2) & (mcdvar == 16),(np.sqrt((0.0225+(0.3721/modamf**2))*modaod**2 + (0.015+(0.1464/modamf**2))*modaod + 0.0025 + 0.0144/modamf**2)/2.)*moddfrac+modaod*(fraccoef1*(moddfrac**4)+fraccoef2*(moddfrac**3)+fraccoef3*(moddfrac**2)+fraccoef4*moddfrac+fraccoef5),dtdblandunarid)
    
    #C061 vegetated pixels
    dtdblandunveg=np.where((modlsf == 1) & (modalgflag == 2) & (mcdvar != 16) & (mcdvar != 0),(np.sqrt((0.0225+(0.4489/modamf**2))*modaod**2 + (0.015+(0.10586/modamf**2))*modaod + 0.0025 + 0.006241/modamf**2)/2.)*moddfrac+modaod*(fraccoef1*(moddfrac**4)+fraccoef2*(moddfrac**3)+fraccoef3*(moddfrac**2)+fraccoef4*moddfrac+fraccoef5),dtdblandunveg)
    
    #Merge the six arrays corresponding to each MODIS algorithm
    dodunctmp=np.nansum(np.dstack((dtlandun,dtoceanun,dblandunarid,dblandunveg,dtdblandunarid,dtdblandunveg)),2)
    
    #Replace zeros with nans in dodunctmp array
    dodunctmp1=np.where(modaod.mask==True,0.,dodunctmp)
    
    dodunctmp11=np.where(modaod.mask==True,np.nan,dodunctmp)
    
    #Final dod uncertainty
    dodunctmp2=ma.masked_invalid(dodunctmp1)
    
    dodunctmp22=ma.masked_invalid(dodunctmp11)
    
    #Calculate the square of the dodunc which will be the output of the function
    sqrdodunc=dodunctmp2**2
    
    '''
    Returns the: (1) square of the DOD uncertainty
    (2) MODIS DOD
    (3) DOD uncertainty
    (4) DOD uncertainty with masked values
    
    '''
    return sqrdodunc,moddod,dodunctmp2,dodunctmp22

#%%
# =============================================================================
# Read the daily MODIS netcdf file
# =============================================================================
for modisfile in modisfiles:
#for modisfile in modisfiles[:365]:
#for modisfile in modisfiles[:31]:
    
    print (modisfile)
    
    #Here is the modis file date
    modfiledate=dt.datetime.strptime(modisfile[-11:-3],'%Y%m%d')
    
    #Here is the daily path and file
    dpath=modisfile[:-55]
    dfile=modisfile[-55:]
    
    #The DOD uncertainty square for each file    
    #doduncsqr=dod_uncertainty_square(dpath,dfile)
    moddodout=dod_modis_outputs(dpath,dfile)
           
# =============================================================================
#     #Create the DOD uncertainty netcdfile
# =============================================================================
    ncdoduncfolder=doduncfolder+satcol+'/'+satellite+'/'+str(modfiledate.year)+'/'+satres+'/'
    
    if not os.path.exists(ncdoduncfolder):
        os.makedirs(ncdoduncfolder)
    
    outdoduncncfile=satellite+'-DOD-UNCERTAINTY-'+satres+'-'+modisfile[-11:-3]+'.nc'
    
    ncdoduncfile=Dataset(os.path.join(ncdoduncfolder,outdoduncncfile),'w',format='NETCDF4')
    
    #Create Global Attributes
    ncdoduncfile.description='MIDAS_DOD_UNCERTAINTY_WITH_MASKED_ELEMENTS'
    ncdoduncfile.history='Created:  '+time.ctime(time.time())
    ncdoduncfile.source='netCDF4 python module tutorial'
    ncdoduncfile.Conventions = 'CF-1.6'
    
    ncdoduncfilelat=ncdoduncfile.createDimension('Latitude',moddodout[0].shape[0])
    ncdoduncfilelon=ncdoduncfile.createDimension('Longitude',moddodout[0].shape[1])
    ttime=ncdoduncfile.createDimension('Time',None)
    
    #Latitude
    Latitude=ncdoduncfile.createVariable('Latitude',np.float32,('Latitude','Longitude'),
                                    fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
        
    Latitude[:]=modlat
    
    Latitude.valid_range=(-90.,90.)
    Latitude.long_name='Geodetic Latitude'
    Latitude.units='degrees_north'
    
    #Longitude
    Longitude=ncdoduncfile.createVariable('Longitude',np.float32,('Latitude','Longitude'),
                                       fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel) 
    
    Longitude[:]=modlon
            
    Longitude.valid_range=(-180.,180.)
    Longitude.long_name='Geodetic Longitude'
    Longitude.units='degrees_east'
    
    #MIDAS DOD uncertainty
    Middodunc=ncdoduncfile.createVariable('MIDAS_DOD_UNCERTAINTY',np.float32,('Latitude','Longitude'),
                                       fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
    
    Middodunc[:]=moddodout[3]
            
    Middodunc.long_name='MIDAS DOD uncertainty with masked values'
    Middodunc.units='None'
    
    
    #Time
    modtime=ncdoduncfile.createVariable('Time',np.float32,('Time',),zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
           
    modtime[:]=ncf.date2num(modfiledate,units='days since '+modfiledate.strftime('%Y-%m-%d %H:%M:%S'))
    
    #modtime[:]=yeardates[row].strftime('%Y-%m-%d %H:%M:%S')
    
    modtime.long_name='Time'
    modtime.units='days since '+modfiledate.strftime('%Y-%m-%d %H:%M:%S')
    modtime.calendar='julian'
       
    ncdoduncfile.close()
  
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
