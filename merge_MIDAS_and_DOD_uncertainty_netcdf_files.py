#%%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 12:49:33 2020

@author: antonis
"""

from nco import Nco #import
import glob2
#import fnmatch
import os

#%%
#https://github.com/nco/pynco
nco = Nco()

#Define which satellite are you using
satellite='MODIS-AQUA'

#Define here which is the MODIS Collection
modiscollection='C061'

year = '2019'
gridres='GRID_RESOLUTION_0.1'

#Give here the pattern in order to select the appropriate files (i.e. grid resolution)
filepattern='*GRID_RESOLUTION_0.1*'

#Define here the path where the data are stored
modpath='/mnt/nas-1/MIDAS-Processing_Data/DAILY-MODIS-FILES-THROUGH-SYNERGY-OF-MODIS-MERRA2/'+satellite+'/'+year+'/GRID_RESOLUTION_0.1/' #add year,3rd script

#List all MODIS-MIDAS files according to the filepattern
modisfiles = sorted(glob2.glob(modpath+'/**/'+filepattern+'.nc'))

#%%
for modisfile in modisfiles:
    
    print (modisfile)
    
    if satellite == 'MODIS-AQUA':
        
        ncfile1path=modisfile[:119]
        ncfile1=modisfile[119:]
        
    elif satellite == 'MODIS-TERRA':
        
        ncfile1path=modisfile[:122]
        ncfile1=modisfiles[123:]   
        
    #Create the path and the filename of the uncertainty netcdf file
    ncfile2path='/mnt/nas-1/MCD12C1/DOD_UNCERTAINTY_NETCDF_DAILY_FILES/'+modiscollection+'/'+satellite+'/'+ncfile1[-11:-7]+'/'+filepattern[1:-1]#+year, 5th script
                   
    ncfile2= satellite+'-DOD-UNCERTAINTY-'+filepattern[1:-1]+'-'+ncfile1[-11:-3]+'.nc'

    #Merge the two files!!!
    nco.ncks(input=(os.path.join(ncfile2path,ncfile2)),output=(os.path.join(ncfile1path,ncfile1)),options=["-A","--netcdf4"])

# %%
