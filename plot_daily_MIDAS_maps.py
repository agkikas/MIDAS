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

nco = Nco()

# %%
#Define here info for the paths
satellite = 'MODIS-AQUA'

#Define the path where the MIDAS files are stored
rootpath = '/home/agkikas/MIDAS-PLOTS/MIDAS-FILES/'+satellite+'/'

#%%
#Define here the study domain
studyregion = 'NORTH_ATLANTIC'
northlat = 50
southlat = -10
westlon  = -90
eastlon  = 0

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

# %%
#Define the path where the results will be stored
outputfolder='/home/agkikas/MIDAS-PLOTS/MIDAS-MAPS/'+satellite+'/'+studyregion+'/'

if not os.path.exists(outputfolder):
    os.makedirs(outputfolder)

# %%
#List all files
midasfiles=sorted(glob.glob(rootpath + '/**/*.nc', recursive=True))

#################################################################################################
#################################################################################################
#################################################################################################

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

        #Map projection
        mapproj=ccrs.PlateCarree()

        fig = plt.figure(figsize=(figx,figy))

        ax = fig.add_subplot(1,1,1,projection=mapproj,facecolor='darkgrey')

        cbar_ax = fig.add_axes([0, 0, 0.1, 0.1])

        ax.set_extent([westlon, eastlon, southlat, northlat], crs=ccrs.PlateCarree())

        posn = ax.get_position()

        cf=ax.pcolormesh(ncmidasreg.Longitude,ncmidasreg.Latitude,ncmidasreg['Modis-total-dust-optical-depth-at-550nm'],
                        transform=ccrs.PlateCarree(),cmap=cmap,norm=norm)

        #cf.cmap.set_under(cbarunder)
        cf.cmap.set_over(cbarover)

        # colorbar and labels
        cb = plt.colorbar(cf,cmap=cmap,spacing='proportional',norm=norm,ticks=bounds,
                        boundaries=bounds,cax=cbar_ax,extend=cbarextend)

        cb.ax.set_ylabel(cbartext,fontsize=cbarlbsize,labelpad=cbarlbpad,fontweight=cbarlbweight,rotation=270) 
        cb.ax.set_position([posn.x0 + posn.width + 0.05, posn.y0,0.04, posn.height])

        titletext = satellite.replace('MODIS','MIDAS')+'  '+dt.datetime.strftime(midasdate,'%d %b %Y')

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

        outputfilename = satellite.replace('MODIS','MIDAS')+'_DOD550_'+dt.datetime.strftime(midasdate,'%d%m%Y')+'.png'

        plt.savefig(os.path.join(outputfolder,outputfilename),dpi=600,bbox_inches='tight')

        plt.show()

        plt.close()

        fig.clear()

# %%
