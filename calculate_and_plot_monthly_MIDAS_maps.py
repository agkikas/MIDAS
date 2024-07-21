#%%
#Import libraries
from pathlib import Path
from nco import Nco
from nco.custom import Rename
from netCDF4 import Dataset
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import os
import numpy as np
import glob
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs  # Projections list
import cartopy.feature as cfeature
import matplotlib as mpl

# %%
#Define here info for the paths
satellite = 'MODIS-AQUA'

#Define here the year
year = '2016'

#Grid resolution
gridres = 'GRID_RESOLUTION_0.1'

#Define the path where the MIDAS files are stored
rootpath = '/home/agkikas/MIDAS-SS/'+satellite+'/'+year+'/'+gridres+'/'

#Define the path where the regional netcdfs will be stored
regncfolder= '/home/agkikas/MIDAS-SS/REGIONAL-NETCDF-FILES/'+satellite+'/'+year+'/'+gridres+'/'

if not os.path.exists(regncfolder):
    os.makedirs(regncfolder)

#%%
#Define here the study domain
studyregion = 'SAHARA_NORTH_ATLANTIC'
northlat = 50
southlat = 0
westlon  = -50
eastlon  = 50

#Colormap
cmap = plt.get_cmap('Oranges')
cbarover = 'red'
cbarextend = 'max'

#Colorbar
cbarlbsize=40
cbarlbpad=50
cbarlbweight='bold'

cbartickssize=26
cbarticksfont='bold'
xticksdist=2.0

#Colorbar text
cbartext = 'DOD@550nm'

#Title features
titletextfw='bold'
titletextfs=40
titletexty=1.05

#Info for the variables range, ticks!
vmin=0
vmax=1
vspaces=11
vspacesticks=5

bounds = np.linspace(vmin, vmax, vspaces)
norm   = mpl.colors.BoundaryNorm(bounds, cmap.N)
ticks  = np.linspace(vmin, vmax, vspaces)

#Features for the timeseries plots!!
#Set figx and figy
figx=20
figy=10

#Steps for meridians and parallels
merstep=5
parstep=5

#Grid lines on the map
grlinethick = 0.5
grlinestyle = '--'
grlinecolor = 'black'
grlinealpha = 1
grlblfsize = 20

#Define the path where the results will be stored
outputfolder='/home/agkikas/MIDAS-SS/MIDAS-MAPS/'+studyregion+'/'+satellite+'/'+year+'/'+gridres+'/'

if not os.path.exists(outputfolder):
    os.makedirs(outputfolder)

# %%
#List all files
midasfiles=sorted(glob.glob(rootpath + '/**/*.nc', recursive=True))

# %%
#Loop in midasfiles
for subdir, dirs, files in os.walk(rootpath):
    for midasfile in sorted(files):
        
        #print (midasfile)

        #Extract date from the midas filename
        filedate = midasfile[-11:-3]  
        midasdate=dt.datetime.strptime(filedate,'%Y%m%d').date()
        
        try:
        
            #Read the netcdf file! Change the coordinates dimension names
            ncmidas = Dataset(os.path.join(rootpath,midasfile),'r+')

            ncmidas.renameDimension(u'Longitude',u'Lon')
            ncmidas.renameDimension(u'Latitude',u'Lat')

            ncmidas.close()
            
        except:
            
            print ('You do not need to change coordinates names in '+midasfile+' !!!')

        #Plot the daily MIDAS file!!!
        ncmidasglob = xr.open_dataset(os.path.join(rootpath,midasfile))

        #Extract the regional domain
        ncmidasreg = ncmidasglob.where((ncmidasglob.Latitude>=southlat) & (ncmidasglob.Latitude<=northlat) & \
                                    (ncmidasglob.Longitude>=westlon) & (ncmidasglob.Longitude<=eastlon),drop=True)
        
        #Save the regional netcdf files
        ncmidasregfilename = studyregion+'_'+midasfile
        
        ncmidasreg.to_netcdf(os.path.join(regncfolder,ncmidasregfilename),mode='w',format='NETCDF4',engine='netcdf4')
        
# %%
#List all regional MIDAS files
regmidasfiles=sorted(glob.glob(regncfolder + '/**/*.nc', recursive=True))

# %%
#Read all files at once!!!
combncs=xr.open_mfdataset(regmidasfiles,parallel=True,coords='minimal',data_vars='minimal',compat='override',concat_dim='Time',combine='nested', name)

# %%
#Select only the MIDAS DOD
midasdod = combncs['Modis-total-dust-optical-depth-at-550nm']

#%%
#Calculate the monthly means
regdommmean = midasdod.resample(Time='1MS').mean(dim='Time')


# %%
