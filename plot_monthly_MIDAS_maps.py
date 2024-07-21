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

#Define the path where the regional netcdfs will be stored
rootpath= '/home/agkikas/MIDAS-SS/REGIONAL-INTRANNUAL-NETCDF-FILES/'+satellite+'/'+year+'/'+gridres+'/'
    
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
merstep=10
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
        
        print (midasfile)
        
#Plot the monthly MIDAS file!!!
regdommmean = xr.open_dataset(os.path.join(rootpath,midasfile))

# %%
#Create an array with the indices of the months of the year
monidxs=np.arange(1,13)
monnames=['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']

for monidx,monname in zip(monidxs,monnames):
    
    print (monname,monidx)

    #Select each month for the regional domain (MIDAS DOD)
    monregdommmean=regdommmean.sel(Time=regdommmean.indexes['Time'].to_datetimeindex().month==monidx)

    ###################################################################################
    #Plot the monthly maps of the MIDAS DOD
    #Map projection
    mapproj=ccrs.PlateCarree()

    fig = plt.figure(figsize=(figx,figy))

    ax = fig.add_subplot(1,1,1,projection=mapproj,facecolor='darkgrey')

    cbar_ax = fig.add_axes([0, 0, 0.1, 0.1])

    ax.set_extent([westlon, eastlon, southlat, northlat], crs=ccrs.PlateCarree())

    posn = ax.get_position()

    cf=ax.pcolormesh(monregdommmean.Longitude.values,monregdommmean.Latitude.values,monregdommmean['Modis-total-dust-optical-depth-at-550nm'].values.squeeze(),
                            transform=ccrs.PlateCarree(),cmap=cmap,norm=norm)

    #cf.cmap.set_under(cbarunder)
    cf.cmap.set_over(cbarover)

    # colorbar and labels
    cb = plt.colorbar(cf,cmap=cmap,spacing='proportional',norm=norm,ticks=bounds,
                    boundaries=bounds,cax=cbar_ax,extend=cbarextend)

    cb.ax.set_ylabel(cbartext,fontsize=cbarlbsize,labelpad=cbarlbpad,fontweight=cbarlbweight,rotation=270) 
    cb.ax.set_position([posn.x0 + posn.width + 0.05, posn.y0,0.04, posn.height])

    titletext = satellite.replace('MODIS','MIDAS')+'  '+monname+' '+year

    ax.set_title(titletext,fontweight=titletextfw,fontsize=titletextfs,y=titletexty)

    for t in cb.ax.get_yticklabels():
                
        t.set_fontsize(cbartickssize)
        t.set_fontweight(cbarticksfont)
        t.set_horizontalalignment('right')   
        t.set_x(xticksdist)

    gl = ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(), lw=0.75,ls='--',color="white",
                    y_inline=False, xlocs=range(westlon,eastlon,merstep), ylocs=range(southlat,northlat,parstep))

    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'color': 'black', 'weight': 'bold', 'size': grlblfsize}
    gl.ylabel_style = {'color': 'black', 'weight': 'bold', 'size': grlblfsize}

    ax.coastlines(linewidth=0.5,color='k')

    resol = '10m'  # use data at this scale
    bodr = cfeature.NaturalEarthFeature(category='cultural', 
        name='admin_0_boundary_lines_land', scale=resol, facecolor='none')
    # land = cfeature.NaturalEarthFeature('physical', 'land', \
    #     scale=resol, edgecolor='k', facecolor=cfeature.COLORS['land'])
    # ocean = cfeature.NaturalEarthFeature('physical', 'ocean', \
    #     scale=resol, edgecolor='none', facecolor=cfeature.COLORS['water'])
    # lakes = cfeature.NaturalEarthFeature('physical', 'lakes', \
    #     scale=resol, edgecolor='b', facecolor=cfeature.COLORS['water'])

    #ax.add_feature(land, facecolor='grey')
    #ax.add_feature(ocean,facecolor='Gainsboro')
    #ax.add_feature(lakes)
    ax.add_feature(bodr, linestyle='--', edgecolor='k', alpha=1)

    outputfilename = satellite.replace('MODIS','MIDAS')+'_'+studyregion+'_DOD550_'+monname+'_'+str(monregdommmean.indexes['Time'].to_datetimeindex().year[0])+'.png'

    plt.savefig(os.path.join(outputfolder,outputfilename),dpi=600,bbox_inches='tight')

    del monregdommmean

    plt.show()

    plt.close()

    fig.clear()



# %%
