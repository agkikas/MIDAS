#%%
# #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  7 09:02:17 2018

@author: antonis
"""

from pyhdf.SD import SD,SDC
from netCDF4 import Dataset
from scipy.spatial import cKDTree

import datetime as dt
import os
import numpy.ma as ma
import numpy as np
import calendar
import pandas as pd
import time
import fnmatch
import itertools
import time

####################################################################################################################
scriptstart=dt.datetime.now()
#f = open("/mnt/c/Users/Alexiou_NOA/Desktop/n1_2019.log", "w")

print (scriptstart)
#f.write(str(scriptstart) + '\n')
#f.flush()

#%%

#Define here the year
year='2022'  #fix for run!!! (start with 2020 - 2001)

#Define here strings which are required in order to read the MODIS swaths files
shortsatname='myd04-col61'#myd for Aqua -  mod for Terra
longsatname='MYD04_L2' #MYD for Aqua -  MOD for Terra
satellite='MODIS-AQUA' #TERRA/AQUA
pathpart='MODIS-Aqua-Collection6.1-L2'

# shortsatname='mod04-col61'#myd for Aqua -  mod for Terra
# longsatname='MOD04_L2' #MYD for Aqua -  MOD for Terra
# satellite='MODIS-TERRA' #TERRA/AQUA
# pathpart='MODIS-Terra-Collection6.1-L2'

#Define here which version of MERRA2 data you will use depending on the year
if year<='2010':
    
    MERRAVERSION=str(300)
    
elif year>'2010':
    
    MERRAVERSION=str(400) 
    MERRAVERSION2=str(401) #Sept 2020 is 401!

#Threshold for the cloud fraction
aodcfthre=0.8

#Define here the minutes interval that merra provides outputs
merratimestep=60

#Define the least significant digits and compression level for the variables which will be stored in the netcdf files
lstsigdig=10
complevel=9 #This varies from 1 to 9 (http://unidata.github.io/netcdf4-python/)

#Give the path where the MODIS Level2 data are stored 
#modisfilespath='/mnt/nas2/MODIS-Aqua-Collection6.1-L2/'+longsatname+'/'+year+'/'   #Aqua/Terra
modisfilespath='/mnt/nas-2/'+pathpart+'/'+longsatname+'/'+year+'/'   #Aqua/Terra
#? ORIGINAL CODE: /media/antonis/My Passport/MERRA-2/obs/nasa/'+shortsatname+'/original_files/

#Give the path where the MERRA2 outputs (M2T1NXAER) are stored
#merrafilespath='/media/antonis/My Book/MERRA-2/'+year+'/'

#Give the path where the MERRA2 outputs (M2T1NXADG) are stored
#merrafilespath2='/media/antonis/My Book/MERRA-2/'+year+'/M2T1NXADG/'

#Create the folder where the MODIS swaths with DODs will be stored
modswathsdodfolder='/mnt/nas-1/MIDAS-Processing_Data/MODIS-SWATHS-WITH-DOD-USING-MERRA2-DUST-RATIO/'+satellite+'/'+year+'/' 
#? ORIGINAL CODE: /media/antonis/My Passport/MERRA-2/MODIS-SWATHS-WITH-DOD-USING-MERRA2-DUST-RATIO/

if not os.path.exists(modswathsdodfolder):
    os.makedirs(modswathsdodfolder)

#%%

#Function for reading and processing MODIS Level 2 swaths
#The filters of buddy checks, cloud fraction, angstrom and aod above ocean are applied within the function
def read_MODIS_level2_swath(inpath,infile,aodcfthre):
    
    #Here create a list with the 2D variables that you want to extract from the MODIS Level 2 HDF file
    sdslist=['Latitude','Longitude','Aerosol_Cloud_Fraction_Land','Aerosol_Cloud_Fraction_Ocean',
             'AOD_550_Dark_Target_Deep_Blue_Combined','AOD_550_Dark_Target_Deep_Blue_Combined_Algorithm_Flag',
             'AOD_550_Dark_Target_Deep_Blue_Combined_QA_Flag','Sensor_Azimuth','Sensor_Zenith','Solar_Azimuth','Solar_Zenith','Land_sea_Flag',
             'Angstrom_Exponent_1_Ocean','Angstrom_Exponent_2_Ocean','Deep_Blue_Angstrom_Exponent_Land'] #,
             #'Effective_Optical_Depth_0p55um_Ocean','Deep_Blue_Aerosol_Optical_Depth_550_Land_Best_Estimate']
    
    #Read the hdf file
    hdffile=SD(os.path.join(inpath,infile),SDC.READ)
    
    #Create an empty list where all 2D variables will be stored
    vrbllist=[]
    
    #Read each SDS and covert to a physical variable
    for sds in sdslist:
    
        sds_obj=hdffile.select(sds) #Read the sdsvariable
    
        add_offset=sds_obj.attributes()['add_offset']
        scale_factor=sds_obj.attributes()['scale_factor']
        fill_value=sds_obj.attributes()['_FillValue']
    
        #Get the physical variable for each SDS (scientific data set)
        vrbl=sds_obj.get()
        vrbl=ma.masked_equal(vrbl,fill_value)
        vrbl_fin=(vrbl-add_offset)*scale_factor
        
        vrbllist.append(vrbl_fin)
        
    #Convert the list in a masked array
    modisvrbls=ma.array(vrbllist)
    
    #extract the masked arrays for cloud fraction land and ocean as well as the combined AOD
    modiscfland=np.array(modisvrbls[2,:,:])
    modiscfocean=np.array(modisvrbls[3,:,:])
    modisaodextracted=np.array(modisvrbls[4,:,:])
        
    #Extract the algorithm flag as well as the QA flag
    modisalgflag=np.array(modisvrbls[5,:,:])
    modisQAflag=np.array(modisvrbls[6,:,:])
    
    #Extract the solar and viewing angles
    modissensorazimuth=np.array(modisvrbls[7,:,:])
    modissensorzenith=np.array(modisvrbls[8,:,:])
    modissolarazimuth=np.array(modisvrbls[9,:,:])
    modissolarzenith=np.array(modisvrbls[10,:,:])
    
    #Extract the land sea flag
    modislandseaflag=np.array(modisvrbls[11,:,:])
    
    #Extract the angstrom exponent above ocean (2 exponents) and land (deep blue product)
    modisangocean1=np.array(modisvrbls[12,:,:])
    modisangocean2=np.array(modisvrbls[13,:,:])
    modisangland=np.array(modisvrbls[14,:,:])
    
    #Extract the DT-Ocean (QA higher/equal than 1) and DB-Land AODs
    #modisaodDTocean=np.array(modisvrbls[15,:,:])
    #modisaodDBland=np.array(modisvrbls[16,:,:])

###If you want to add 2d variables (listed in the sdslist) continue from here!!!!

####################################################################################################################    
#    #Extract the DT-Land AOD (QA=3)
#    sds_obj_DTland=hdffile.select('Corrected_Optical_Depth_Land') #Read the sdsvariable
#    
#    add_offset_DTland=sds_obj_DTland.attributes()['add_offset']
#    scale_factor_DTland=sds_obj_DTland.attributes()['scale_factor']
#    fill_value_DTland=sds_obj_DTland.attributes()['_FillValue']
#    
#    #Get the physical variable for each SDS (scientific data set)
#    vrbl_DTland=sds_obj_DTland.get()
#    vrbl_DTland=ma.masked_equal(vrbl_DTland,fill_value_DTland)
#    vrbl_fin_DTland=(vrbl_DTland-add_offset_DTland)*scale_factor_DTland    
#
#    modisaodDTland=vrbl_fin_DTland[1,:,:]
####################################################################################################################
        
    #Mask AODs based on cloud fraction (higher than 0.8)
    modisAODtmp=ma.masked_where((modiscfland>aodcfthre) | (modiscfocean>aodcfthre),modisaodextracted)
    
    #Create the final modisaod by excluding the -9999 values
    modisAOD=ma.masked_equal(modisAODtmp,-9999.)
    
    #Apply the buddy check filter (removal of AOD retrievals where all the surrounding pixels are not available)
    for row,col in itertools.product(np.arange(1,modisAOD.shape[0]-1,1),np.arange(1,modisAOD.shape[1]-1,1)):
        
        buddyarray=np.array([modisAOD[row-1,col-1],modisAOD[row-1,col],modisAOD[row-1,col+1],
                           modisAOD[row,col-1],modisAOD[row,col],modisAOD[row,col+1],
                           modisAOD[row+1,col-1],modisAOD[row+1,col],modisAOD[row+1,col+1]])
        
        buddyarray=ma.masked_invalid(buddyarray)
        
        if ((buddyarray[0:4].count()==0) & (buddyarray[5:].count()==0)):
                        
                        modisAOD[row,col]=np.nan
    
    #Mask nans
    modisAOD=ma.masked_invalid(modisAOD)
    
    #Mask the other arrays depending on the mask of MODIS AOD
    modisalgflag=ma.masked_where(modisAOD.mask==True,modisalgflag)
    modisQAflag=ma.masked_where(modisAOD.mask==True,modisQAflag)
    modissensorazimuth=ma.masked_where(modisAOD.mask==True,modissensorazimuth)
    modissensorzenith=ma.masked_where(modisAOD.mask==True,modissensorzenith)
    modissolarazimuth=ma.masked_where(modisAOD.mask==True,modissolarazimuth)
    modissolarzenith=ma.masked_where(modisAOD.mask==True, modissolarzenith)
    modislandseaflag=ma.masked_where(modisAOD.mask==True,modislandseaflag)
    modisangocean1=ma.masked_where(modisAOD.mask==True,modisangocean1)
    modisangocean2=ma.masked_where(modisAOD.mask==True,modisangocean2)
    modisangland=ma.masked_where(modisAOD.mask==True,modisangland)
    #modisaodDTocean=ma.masked_where(modisAOD.mask==True,modisaodDTocean)
    #modisaodDBland=ma.masked_where(modisAOD.mask==True,modisaodDBland)
    #modisaodDTland=ma.masked_where(modisAOD.mask==True,modisaodDTland)    
        
    hdffile.end()
    
    """
    0: Latitude
    1: Longitude
    2: Combined AOD after masking cloud contamination considering buddy checks
    3: MODIS Algorithm Flag
    4: MODIS QA Flag
    5: Sensor_Azimuth
    6: Sensor_Zenith
    7: Solar_Azimuth 
    8: Solar_Zenith
    9: Land Sea Flag
    10: MODIS Angstrom 1 over ocean (550-860nm)
    11: MODIS Angstrom 2 over ocean (860-2130nm)
    12: MODIS Angstrom over land (Deep Blue product)
#########################################################################################################  
    13: Effective_Optical_Depth_0p55um_Ocean
    14: Deep_Blue_Aerosol_Optical_Depth_550_Land_Best_Estimate
    15: Corrected_Optical_Depth_Land at 550nm
    """
    return modisvrbls[0],modisvrbls[1],modisAOD,modisalgflag,modisQAflag,modissensorazimuth,modissensorzenith,modissolarazimuth,modissolarzenith,modislandseaflag,modisangocean1,modisangocean2,modisangland 
    #,modisaodDTocean,modisaodDBland,modisaodDTland

####################################################################################################################
#%%

def extract_attributes_from_MODIS_swath_file(inpath,infile):
    
    sdsvrbls=['Latitude','Longitude','AOD_550_Dark_Target_Deep_Blue_Combined',
              'AOD_550_Dark_Target_Deep_Blue_Combined_Algorithm_Flag','AOD_550_Dark_Target_Deep_Blue_Combined_QA_Flag',
              'Sensor_Azimuth','Sensor_Zenith','Solar_Azimuth','Solar_Zenith','Land_sea_Flag',
              'Angstrom_Exponent_1_Ocean','Angstrom_Exponent_2_Ocean','Deep_Blue_Angstrom_Exponent_Land'] #,
              #'Effective_Optical_Depth_0p55um_Ocean','Deep_Blue_Aerosol_Optical_Depth_550_Land_Best_Estimate']
    
    hdffile=SD(os.path.join(inpath,infile),SDC.READ)
    
    #Here are the lists for variables and their attributes
    vrblslist=[]
    
    vrblshapelist=[]
    
    geolocpointerlist=[]
    longnamelist=[]
    unitslist=[]
    vrangelist=[]
    parametertypelist=[]
    cellalongswathlist=[]
    cellacrossswathlist=[]
    fillvaluelist=[]
    
    for sds in sdsvrbls:   
        
        sds_obj=hdffile.select(sds) #Read the sdsvariable
    
        add_offset=sds_obj.attributes()['add_offset']
        scale_factor=sds_obj.attributes()['scale_factor']
        fill_value=sds_obj.attributes()['_FillValue']  
        
        vrbl=sds_obj.get()
        vrbl=ma.masked_equal(vrbl,fill_value)
        vrbl_fin=(vrbl-add_offset)*scale_factor
    
        vrblslist.append(vrbl_fin)
        
        #Here it finds the shape of each variable
        vrblshape=vrbl.shape
        
        vrblshapelist.append(vrblshape)
        
        #Extract the attributes for each variable
        geolocpointer=sds_obj.attributes()['Geolocation_Pointer']
        longname=sds_obj.attributes()['long_name']
        units=sds_obj.attributes()['units']
        vrange=sds_obj.attributes()['valid_range']        
        parametertyp=sds_obj.attributes()['Parameter_Type']
        cellalongswath=sds_obj.attributes()['Cell_Along_Swath_Sampling']
        cellacrosswath=sds_obj.attributes()['Cell_Across_Swath_Sampling']
        
        geolocpointerlist.append(geolocpointer)
        longnamelist.append(longname)
        unitslist.append(units)
        vrangelist.append(vrange)
        parametertypelist.append(parametertyp)
        cellalongswathlist.append(cellalongswath)
        cellacrossswathlist.append(cellacrosswath)
        fillvaluelist.append(fill_value)    
    hdffile.end()
    
    return vrblshapelist,geolocpointerlist,longnamelist,unitslist,vrangelist,parametertypelist,cellalongswathlist,cellacrossswathlist,fillvaluelist,sdsvrbls
####################################################################################################################
####################################################################################################################
#%%
#This function is used to find the overpass time of MODIS
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
####################################################################################################################
####################################################################################################################
#Work with a test MODIS Level 2 file
#inpath='/media/antonis/My Book/MERRA-2/obs/nasa/myd04-col6/original_files/MYD04_L2/2015/032'
#infile='MYD04_L2.A2015032.1135.006.2015033173903.hdf'

#%%
#Loop through the MODIS swaths files
for subdir,dirs,files in sorted(os.walk(modisfilespath)):
    for filename in sorted(fnmatch.filter(files,'*.hdf')):
    #for file in [x for x in sorted(fnmatch.filter(files,'*.hdf')) if 'A2015032' not in x]:
        #f.write(filename + '\n')
        print (filename)
        #f.flush()
        
        outncfile=filename[:-4]+'-MODIS-SWATH-DOD'+'.nc'
            
        if os.path.exists(os.path.join(modswathsdodfolder,outncfile)):
            #f.write(outncfile + ' exists\n')
            #f.flush()
            continue
        
        try:       

            print (filename)

            #Read the MODIS L2 swath file
            modswath=read_MODIS_level2_swath(subdir,filename,aodcfthre)
            
            #Find the overpass time of MODIS
            modoverpass=find_MODIS_swath_overpass_time(subdir,filename)
            
            #Create the path to point the right MERRA2 (M2T1NXAER) netcdf file  #fix for run, go where Merra Data aer are, make year and month folders!!
            
            merra4modisfile='/mnt/nas-2/Merra/Merra_Archived/'+str(modoverpass.year)+'/'+"{0:0=2d}".format(modoverpass.month)+'/'+'MERRA2_'+MERRAVERSION+'.tavg1_2d_aer_Nx.'+str(modoverpass.year)+"{0:0=2d}".format(modoverpass.month)+"{0:0=2d}".format(modoverpass.day)+'.nc4'

            if os.path.exists (merra4modisfile) == False:

                merra4modisfile='/mnt/nas-2/Merra/Merra_Archived/'+str(modoverpass.year)+'/'+"{0:0=2d}".format(modoverpass.month)+'/'+'MERRA2_'+MERRAVERSION2+'.tavg1_2d_aer_Nx.'+str(modoverpass.year)+"{0:0=2d}".format(modoverpass.month)+"{0:0=2d}".format(modoverpass.day)+'.nc4'
            #? ORIGINAL CODE: /media/antonis/My Passport/MERRA-2/merra_v2/original_files/1hourly/ 
            
            #Create the path to point the right MERRA2 (M2T1NXADG) netcdf file
            
            #merra4modisfile2='/media/antonis/My Book/MERRA-2/'+str(modoverpass.year)+'/M2T1NXADG/'+"{0:0=2d}".format(modoverpass.month)+'/'+'MERRA2_'+MERRAVERSION+'.tavg1_2d_adg_Nx.'+str(modoverpass.year)+"{0:0=2d}".format(modoverpass.month)+"{0:0=2d}".format(modoverpass.day)+'.nc4'

            merra4modisfile2='/mnt/nas-2/Merra/Merra_Archived/'+str(modoverpass.year)+'/'+"{0:0=2d}".format(modoverpass.month)+'/'+'MERRA2_'+MERRAVERSION+'.tavg1_2d_adg_Nx.'+str(modoverpass.year)+"{0:0=2d}".format(modoverpass.month)+"{0:0=2d}".format(modoverpass.day)+'.nc4'
            
            if os.path.exists(merra4modisfile2)==False:

                merra4modisfile2='/mnt/nas-2/Merra/Merra_Archived/'+str(modoverpass.year)+'/'+"{0:0=2d}".format(modoverpass.month)+'/'+'MERRA2_'+MERRAVERSION2+'.tavg1_2d_adg_Nx.'+str(modoverpass.year)+"{0:0=2d}".format(modoverpass.month)+"{0:0=2d}".format(modoverpass.day)+'.nc4'

            #? ORIGINAL CODE: /media/antonis/My Passport/MERRA-2/merra_v2/original_files/1hourly/
            #f.write(merra4modisfile + ' name\n')
            #Read the MERRA-2 annual netcdf file containing daily averages
            ncmerra=Dataset(merra4modisfile,'r')

            #Extract time from the netcdf file
            ncmerraTIME=ncmerra['time']
            ncmerraTIMEunits=ncmerra['time'].units
            ncmerraTIMEstart=ncmerraTIMEunits[14:]
            ncmerraTIMEstartFIN=dt.datetime.strptime(ncmerraTIMEstart,'%Y-%m-%d %H:%M:%S') #This is the first time record in the MERRA2 file

            #Create the list with timesteps in the MERRA2 file
            merra_time_list = [ncmerraTIMEstartFIN + dt.timedelta(minutes=merratimestep*x) for x in range(0, len(ncmerraTIME))]
            #merra_time_list_FIN=[x.strftime('%Y-%m-%d %H:%M:%S') for x in merra_time_list]

            #Find the index from the nearest MERRA2 time output to the MODIS overpass time
            merraMODIStimedifslist=[]
            
            for merratime in merra_time_list:
                
                timediff=np.abs((merratime-modoverpass.to_pydatetime()).total_seconds()/60.)
                
                merraMODIStimedifslist.append(timediff)

            merratimeid=np.argmin(merraMODIStimedifslist)
            #Extract the coordinates from the MERRA2 netcdf file
            #Extract coordinates from ncmerra file
            merralat=ncmerra['lat'][:]
            merralon=ncmerra['lon'][:]

            #Extract the dust and the total extinction AOD from the MERRA2 file for the closest timestep to MODIS overpass
            merradustaod=ncmerra['DUEXTTAU'][merratimeid,:,:]
            merratotaod=ncmerra['TOTEXTTAU'][merratimeid,:,:]

            #Calculate the dust contribution to the total AOD
            merradustcontr=merradustaod/merratotaod

            ncmerra.close()

            #Read the MERRA-2 netcdf file (M2T1NXADG)

            ncmerra2=Dataset(merra4modisfile2,'r')
            merraFINEdustaod=ncmerra2['DUEXTTFM'][merratimeid,:,:]
            ncmerra2.close()

            #Calculate the dust fine and dust coarse contribution to the total dust AOD of MERRA2
            merraFINEdustaodcontr=merraFINEdustaod/merradustaod
            merraCOARSEdustaodcontr=(merradustaod-merraFINEdustaod)/merradustaod   
            
            #Create an array with nans where modswathvrblsattrs[2][4]the modisDOD will be stored
            modswathDOD=np.nan*np.ones((modswath[0].shape[0],modswath[0].shape[1]))
            modswathDODfine=np.nan*np.ones((modswath[0].shape[0],modswath[0].shape[1]))
            modswathDODcoarse=np.nan*np.ones((modswath[0].shape[0],modswath[0].shape[1]))
                        
            ###################################################################
            #Input from Thanasis!!! (13 Jan 2022)
            ###################################################################
            # Create spatial index for merra's lat/lon grid
            # lat and lon are in vectors so we first have to create the grid
            XX, YY = np.meshgrid(merralat, merralon)
            merragrid = np.vstack([XX.ravel(), YY.ravel()]).T
            tr = cKDTree(merragrid) # KDTree spatial index, helps with neighbour search
            
            # Query the nearest neighbour for all MODIS pixels
            # Create x*y*2 array of swath points, there the third dimension is lat/lon
            modswath_points = np.dstack([modswath[0].ravel(), modswath[1].ravel()])
            # Query the nearest neighbour for all MODIS pixels using the KDTree. The reshape brings the
            # result back to the original shape of the swath
            modswath_points_indices = tr.query(modswath_points, k=1)[1].reshape(modswath[0].shape)
            
            # If we wanted the closest MERRA lat/lon for each MODIS pixel, we could use the following:
            # merrapoints = merragrid[modswath_points_indices, :]
            
            # The merragrid array was very long (lat*lon in length). The indices in modswath_points_indices point
            # to this array. We must convert them so they can work in the original MERRA arrays
            indices = np.unravel_index(modswath_points_indices, (merralat.shape[0], merralon.shape[0]), order='F')
            
            # Read data
            merra4modisDUSTRATIOO=merradustcontr[indices]
            merra4modisDUSTfineRATIO=merraFINEdustaodcontr[indices]
            merra4modisDUSTcoarseRATIO=merraCOARSEdustaodcontr[indices]
            
            # Do science
            modswathDOD=modswath[2]*merra4modisDUSTRATIOO
            modswathDODfine=modswathDOD*merra4modisDUSTfineRATIO
            modswathDODcoarse=modswathDOD*merra4modisDUSTcoarseRATIO
            ###################################################################
            ###################################################################
            ###################################################################
                       
            # #Loop in the MODIS Level 2 swath in order to calculate the DUST AOD
            # for modlat in np.arange(modswath[0].shape[0]):
            #     for modlon in np.arange(modswath[0].shape[1]):
                    
            #         #Find the nearest MERRA grid point to MODIS Level 2 cell
            #         merralatpos=np.argmin(np.abs(modswath[0][modlat][modlon]-merralat))
                    
            #         merralonpos=np.argmin(np.abs(modswath[1][modlat][modlon]-merralon))
                    
                    
            #         #print merralatpos,modswath[0][modlat][modlon],merralat[merralatpos]
            #         #print merralonpos,modswath[1][modlat][modlon],merralon[merralonpos]
                    
            #         #Extract the DUST RATIO from the MERRA-2 file (select the nearest gridpoint to MODIS Level 2 cell)
            #         merra4modisDUSTRATIO=merradustcontr[merralatpos,merralonpos]
            #         merra4modisDUSTfineRATIO=merraFINEdustaodcontr[merralatpos,merralonpos]
            #         merra4modisDUSTcoarseRATIO=merraCOARSEdustaodcontr[merralatpos,merralonpos]
                    
            #         modswathDOD[modlat][modlon]=modswath[2][modlat][modlon]*merra4modisDUSTRATIO
                  
                    
            #         modswathDODfine[modlat][modlon]=modswathDOD[modlat][modlon]*merra4modisDUSTfineRATIO
                 
                    
            #         modswathDODcoarse[modlat][modlon]=modswathDOD[modlat][modlon]*merra4modisDUSTcoarseRATIO
                    
            #Here is the final MODIS swath with the DOD       
            modswathDODfin=ma.masked_invalid(modswathDOD)
            modswathDODfinefin=ma.masked_invalid(modswathDODfine)
            modswathDODcoarsefin=ma.masked_invalid(modswathDODcoarse)
            
            #ncmerra.close()
            
            #Run the extract_attributes_from_MODIS_swath_file function to get the variables' attributes
            modswathvrblsattrs=extract_attributes_from_MODIS_swath_file(subdir,filename)
            ####################################################################################################################
            ####################################################################################################################
            #Create a netcdf file with the MODIS DOD swath
            #Create here the netcdf file where the outputs will be stored

            ncmodis=Dataset(os.path.join(modswathsdodfolder,outncfile),'w',format='NETCDF4')
            
            #Create Global Attributes
            ncmodis.description=satellite+'-'+'For data assimilation'
            ncmodis.history='Created:  '+time.ctime(time.time())
            ncmodis.source='netCDF4 python module tutorial'
            ncmodis.Conventions = 'CF-1.6'

            #Create groups
            maingrp=ncmodis.createGroup('mod04')
            datafieldsgrp=ncmodis.createGroup('/mod04/Data_Fields')
            geolocgrp=ncmodis.createGroup('/mod04/Geolocation_Fields')
            

            modislatitude=ncmodis.createDimension('Cell_Along_Swath',modswath[0].shape[0])
            modislongitude=ncmodis.createDimension('Cell_Across_Swath',modswath[0].shape[1])

            #Create the variables which will be stored in the netcdf file
            #Latitude
            Latitude=ncmodis.createVariable('/mod04/Geolocation_Fields/Latitude',np.float32,('Cell_Along_Swath','Cell_Across_Swath'),
                                            fill_value=modswathvrblsattrs[8][0],zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            Latitude[:]=modswath[0]
            
            #Set variable attributes
            Latitude.valid_range=modswathvrblsattrs[4][0]
            Latitude.long_name=modswathvrblsattrs[2][0]
            Latitude.units=modswathvrblsattrs[3][0]
            Latitude.Parameter_Type=modswathvrblsattrs[5][0]
            Latitude.Cell_Along_Swath_Sampling=modswathvrblsattrs[6][0]
            Latitude.Cell_Across_Swath_Sampling=modswathvrblsattrs[7][0]
            Latitude.Geolocation_Pointer=modswathvrblsattrs[1][0]
            Latitude._CoordinateAxisType='Lat'

            #Longitude
            Longitude=ncmodis.createVariable('/mod04/Geolocation_Fields/Longitude',np.float32,('Cell_Along_Swath','Cell_Across_Swath'),
                                            fill_value=modswathvrblsattrs[8][1],zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            Longitude[:]=modswath[1]

            #Set variable attributes
            Longitude.valid_range=modswathvrblsattrs[4][1]
            Longitude.long_name=modswathvrblsattrs[2][1]
            Longitude.units=modswathvrblsattrs[3][1]
            Longitude.Parameter_Type=modswathvrblsattrs[5][1]
            Longitude.Cell_Along_Swath_Sampling=modswathvrblsattrs[6][1]
            Longitude.Cell_Across_Swath_Sampling=modswathvrblsattrs[7][1]
            Longitude.Geolocation_Pointer=modswathvrblsattrs[1][1]
            Longitude._CoordinateAxisType='Lon'

            #Total optical depth
            modisdod=ncmodis.createVariable('/mod04/Data_Fields/AOD_550_Dark_Target_Deep_Blue_Combined',np.float32,('Cell_Along_Swath','Cell_Across_Swath'),
                                            fill_value=modswathvrblsattrs[8][2],zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            modisdod[:]=modswath[2]

            #Set variable attributes
            modisdod.valid_range=modswathvrblsattrs[4][2]
            modisdod.long_name=modswathvrblsattrs[2][2]
            modisdod.units=modswathvrblsattrs[3][2]
            modisdod.Parameter_Type=modswathvrblsattrs[5][2]
            modisdod.Cell_Along_Swath_Sampling=modswathvrblsattrs[6][2]
            modisdod.Cell_Across_Swath_Sampling=modswathvrblsattrs[7][2]
            modisdod.Geolocation_Pointer=modswathvrblsattrs[1][2]

            #Dust optical depth
            modisdod=ncmodis.createVariable('/mod04/Data_Fields/MODIS-TOTAL-Dust-optical-depth-at-550nm',np.float32,('Cell_Along_Swath','Cell_Across_Swath'),
                                            fill_value=modswathvrblsattrs[8][2],zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            modisdod[:]=modswathDODfin

            #Set variable attributes
            modisdod.valid_range=modswathvrblsattrs[4][2]
            modisdod.long_name="MODIS_TOTAL_DOD_at_550nm_based_on_MERRA2_dust-to-total_ratio"
            modisdod.units=modswathvrblsattrs[3][2]
            modisdod.Parameter_Type=modswathvrblsattrs[5][2]
            modisdod.Cell_Along_Swath_Sampling=modswathvrblsattrs[6][2]
            modisdod.Cell_Across_Swath_Sampling=modswathvrblsattrs[7][2]
            modisdod.Geolocation_Pointer=modswathvrblsattrs[1][2]

            #AOD_550_Dark_Target_Deep_Blue_Combined_Algorithm_Flag
            modisalgflag=ncmodis.createVariable('/mod04/Data_Fields/AOD_550_Dark_Target_Deep_Blue_Combined_Algorithm_Flag',np.float32,('Cell_Along_Swath','Cell_Across_Swath'),
                                            fill_value=modswathvrblsattrs[8][3],zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            modisalgflag[:]=modswath[3]
            
            #Set variable attributes
            modisalgflag.valid_range=modswathvrblsattrs[4][3]
            modisalgflag.long_name=modswathvrblsattrs[2][3]
            modisalgflag.units=modswathvrblsattrs[3][3]
            modisalgflag.Parameter_Type=modswathvrblsattrs[5][3]
            modisalgflag.Cell_Along_Swath_Sampling=modswathvrblsattrs[6][3]
            modisalgflag.Cell_Across_Swath_Sampling=modswathvrblsattrs[7][3]
            modisalgflag.Geolocation_Pointer=modswathvrblsattrs[1][3]

            #AOD_550_Dark_Target_Deep_Blue_Combined_QA_Flag
            modisQAflag=ncmodis.createVariable('/mod04/Data_Fields/AOD_550_Dark_Target_Deep_Blue_Combined_QA_Flag',np.float32,('Cell_Along_Swath','Cell_Across_Swath'),
                                            fill_value=modswathvrblsattrs[8][4],zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            modisQAflag[:]=modswath[4]

            #Set variable attributes
            modisQAflag.valid_range=modswathvrblsattrs[4][4]
            modisQAflag.long_name=modswathvrblsattrs[2][4]
            modisQAflag.units=modswathvrblsattrs[3][4]
            modisQAflag.Parameter_Type=modswathvrblsattrs[5][4]
            modisQAflag.Cell_Along_Swath_Sampling=modswathvrblsattrs[6][4]
            modisQAflag.Cell_Across_Swath_Sampling=modswathvrblsattrs[7][4]
            modisQAflag.Geolocation_Pointer=modswathvrblsattrs[1][4]

            #Sensor azimuth
            modisSENSAZIM=ncmodis.createVariable('/mod04/Data_Fields/Sensor_Azimuth',np.float32,('Cell_Along_Swath','Cell_Across_Swath'),
                                            fill_value=modswathvrblsattrs[8][5],zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            modisSENSAZIM[:]=modswath[5]

            #Set variable attributes
            modisSENSAZIM.valid_range=modswathvrblsattrs[4][5]
            modisSENSAZIM.long_name=modswathvrblsattrs[2][5]
            modisSENSAZIM.units=modswathvrblsattrs[3][5]
            modisSENSAZIM.Parameter_Type=modswathvrblsattrs[5][5]
            modisSENSAZIM.Cell_Along_Swath_Sampling=modswathvrblsattrs[6][5]
            modisSENSAZIM.Cell_Across_Swath_Sampling=modswathvrblsattrs[7][5]
            modisSENSAZIM.Geolocation_Pointer=modswathvrblsattrs[1][5]

            #Sensor zenith
            modisSENSZEN=ncmodis.createVariable('/mod04/Data_Fields/Sensor_Zenith',np.float32,('Cell_Along_Swath','Cell_Across_Swath'),
                                            fill_value=modswathvrblsattrs[8][6],zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            modisSENSZEN[:]=modswath[6]

            #Set variable attributes
            modisSENSZEN.valid_range=modswathvrblsattrs[4][6]
            modisSENSZEN.long_name=modswathvrblsattrs[2][6]
            modisSENSZEN.units=modswathvrblsattrs[3][6]
            modisSENSZEN.Parameter_Type=modswathvrblsattrs[5][6]
            modisSENSZEN.Cell_Along_Swath_Sampling=modswathvrblsattrs[6][6]
            modisSENSZEN.Cell_Across_Swath_Sampling=modswathvrblsattrs[7][6]
            modisSENSZEN.Geolocation_Pointer=modswathvrblsattrs[1][6]

            #Solar azimuth
            modisSOLAZIM=ncmodis.createVariable('/mod04/Data_Fields/Solar_Azimuth',np.float32,('Cell_Along_Swath','Cell_Across_Swath'),
                                            fill_value=modswathvrblsattrs[8][7],zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            modisSOLAZIM[:]=modswath[7]

            #Set variable attributes
            modisSOLAZIM.valid_range=modswathvrblsattrs[4][7]
            modisSOLAZIM.long_name=modswathvrblsattrs[2][7]
            modisSOLAZIM.units=modswathvrblsattrs[3][7]
            modisSOLAZIM.Parameter_Type=modswathvrblsattrs[5][7]
            modisSOLAZIM.Cell_Along_Swath_Sampling=modswathvrblsattrs[6][7]
            modisSOLAZIM.Cell_Across_Swath_Sampling=modswathvrblsattrs[7][7]
            modisSOLAZIM.Geolocation_Pointer=modswathvrblsattrs[1][7]

            #Solar zenith
            modisSOLZEN=ncmodis.createVariable('/mod04/Data_Fields/Solar_Zenith',np.float32,('Cell_Along_Swath','Cell_Across_Swath'),
                                            fill_value=modswathvrblsattrs[8][8],zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            modisSOLZEN[:]=modswath[8]

            #Set variable attributes
            modisSOLZEN.valid_range=modswathvrblsattrs[4][8]
            modisSOLZEN.long_name=modswathvrblsattrs[2][8]
            modisSOLZEN.units=modswathvrblsattrs[3][8]
            modisSOLZEN.Parameter_Type=modswathvrblsattrs[5][8]
            modisSOLZEN.Cell_Along_Swath_Sampling=modswathvrblsattrs[6][8]
            modisSOLZEN.Cell_Across_Swath_Sampling=modswathvrblsattrs[7][8]
            modisSOLZEN.Geolocation_Pointer=modswathvrblsattrs[1][8]
            
            #Land and Sea flag
            modisLANDSEAFLAG=ncmodis.createVariable('/mod04/Data_Fields/Land_sea_Flag',np.float32,('Cell_Along_Swath','Cell_Across_Swath'),
                                            fill_value=modswathvrblsattrs[8][9],zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            modisLANDSEAFLAG[:]=modswath[9]
            
            #Set variable attributes
            modisLANDSEAFLAG.valid_range=modswathvrblsattrs[4][9]
            modisLANDSEAFLAG.long_name=modswathvrblsattrs[2][9]
            modisLANDSEAFLAG.units=modswathvrblsattrs[3][9]
            modisLANDSEAFLAG.Parameter_Type=modswathvrblsattrs[5][9]
            modisLANDSEAFLAG.Cell_Along_Swath_Sampling=modswathvrblsattrs[6][9]
            modisLANDSEAFLAG.Cell_Across_Swath_Sampling=modswathvrblsattrs[7][9]
            modisLANDSEAFLAG.Geolocation_Pointer=modswathvrblsattrs[1][9]

            #Fine Dust optical depth
            modisFINEdod=ncmodis.createVariable('/mod04/Data_Fields/MODIS-FINE-Dust-optical-depth-at-550nm',np.float32,('Cell_Along_Swath','Cell_Across_Swath'),
                                            fill_value=modswathvrblsattrs[8][2],zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            modisFINEdod[:]=modswathDODfinefin

            #Set variable attributes
            modisFINEdod.valid_range=modswathvrblsattrs[4][2]
            modisFINEdod.long_name="MODIS_FINE_DOD_at_550nm_based_on_MERRA2_M2T1NXAER_and_M2T1NXADG_outputs"
            modisFINEdod.units=modswathvrblsattrs[3][2]
            modisFINEdod.Parameter_Type=modswathvrblsattrs[5][2]
            modisFINEdod.Cell_Along_Swath_Sampling=modswathvrblsattrs[6][2]
            modisFINEdod.Cell_Across_Swath_Sampling=modswathvrblsattrs[7][2]
            modisFINEdod.Geolocation_Pointer=modswathvrblsattrs[1][2]

            #Coarse Dust optical depth
            modisCOARSEdod=ncmodis.createVariable('/mod04/Data_Fields/MODIS-COARSE-Dust-optical-depth-at-550nm',np.float32,('Cell_Along_Swath','Cell_Across_Swath'),
                                            fill_value=modswathvrblsattrs[8][2],zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            modisCOARSEdod[:]=modswathDODcoarsefin

            #Set variable attributes
            modisCOARSEdod.valid_range=modswathvrblsattrs[4][2]
            modisCOARSEdod.long_name="MODIS_COARSE_DOD_at_550nm_based_on_MERRA2_M2T1NXAER_and_M2T1NXADG_outputs"
            modisCOARSEdod.units=modswathvrblsattrs[3][2]
            modisCOARSEdod.Parameter_Type=modswathvrblsattrs[5][2]
            modisCOARSEdod.Cell_Along_Swath_Sampling=modswathvrblsattrs[6][2]
            modisCOARSEdod.Cell_Across_Swath_Sampling=modswathvrblsattrs[7][2]
            modisCOARSEdod.Geolocation_Pointer=modswathvrblsattrs[1][2]

            #Ocean Angstrom 1 
            modisOCEANANG1=ncmodis.createVariable('/mod04/Data_Fields/Angstrom_Exponent_1_Ocean',np.float32,('Cell_Along_Swath','Cell_Across_Swath'),
                                            fill_value=modswathvrblsattrs[8][10],zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            modisOCEANANG1[:]=modswath[10]

            #Set variable attributes
            modisOCEANANG1.valid_range=modswathvrblsattrs[4][10]
            modisOCEANANG1.long_name=modswathvrblsattrs[2][10]
            modisOCEANANG1.units=modswathvrblsattrs[3][10]
            modisOCEANANG1.Parameter_Type=modswathvrblsattrs[5][10]
            modisOCEANANG1.Cell_Along_Swath_Sampling=modswathvrblsattrs[6][10]
            modisOCEANANG1.Cell_Across_Swath_Sampling=modswathvrblsattrs[7][10]
            modisOCEANANG1.Geolocation_Pointer=modswathvrblsattrs[1][10]

            #Ocean Angstrom 2 
            modisOCEANANG2=ncmodis.createVariable('/mod04/Data_Fields/Angstrom_Exponent_2_Ocean',np.float32,('Cell_Along_Swath','Cell_Across_Swath'),
                                            fill_value=modswathvrblsattrs[8][11],zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            modisOCEANANG2[:]=modswath[11]

            #Set variable attributes
            modisOCEANANG2.valid_range=modswathvrblsattrs[4][11]
            modisOCEANANG2.long_name=modswathvrblsattrs[2][11]
            modisOCEANANG2.units=modswathvrblsattrs[3][11]
            modisOCEANANG2.Parameter_Type=modswathvrblsattrs[5][11]
            modisOCEANANG2.Cell_Along_Swath_Sampling=modswathvrblsattrs[6][11]
            modisOCEANANG2.Cell_Across_Swath_Sampling=modswathvrblsattrs[7][11]
            modisOCEANANG2.Geolocation_Pointer=modswathvrblsattrs[1][11]

            #Land Angstrom (Deep Blue product) 
            modisANGLAND=ncmodis.createVariable('/mod04/Data_Fields/Deep_Blue_Angstrom_Exponent_Land',np.float32,('Cell_Along_Swath','Cell_Across_Swath'),
                                            fill_value=modswathvrblsattrs[8][12],zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
            modisANGLAND[:]=modswath[12]

            #Set variable attributes
            modisANGLAND.valid_range=modswathvrblsattrs[4][12]
            modisANGLAND.long_name=modswathvrblsattrs[2][12]
            modisANGLAND.units=modswathvrblsattrs[3][12]
            modisANGLAND.Parameter_Type=modswathvrblsattrs[5][12]
            modisANGLAND.Cell_Along_Swath_Sampling=modswathvrblsattrs[6][12]
            modisANGLAND.Cell_Across_Swath_Sampling=modswathvrblsattrs[7][12]
            modisANGLAND.Geolocation_Pointer=modswathvrblsattrs[1][12]

###############################################################################################################################################################################################
            
#            #Effective_Optical_Depth_0p55um_Ocean 
#            modisAODDTOCEAN=ncmodis.createVariable('/mod04/Data_Fields/Effective_Optical_Depth_0p55um_Ocean',np.float32,('Cell_Along_Swath','Cell_Across_Swath'),
#                                            fill_value=modswathvrblsattrs[8][13],zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
#            modisAODDTOCEAN[:]=modswath[13]
#            
#            #Set variable attributes
#            modisAODDTOCEAN.valid_range=modswathvrblsattrs[4][13]
#            modisAODDTOCEAN.long_name=modswathvrblsattrs[2][13]
#            modisAODDTOCEAN.units=modswathvrblsattrs[3][13]
#            modisAODDTOCEAN.Parameter_Type=modswathvrblsattrs[5][13]
#            modisAODDTOCEAN.Cell_Along_Swath_Sampling=modswathvrblsattrs[6][13]
#            modisAODDTOCEAN.Cell_Across_Swath_Sampling=modswathvrblsattrs[7][13]
#            modisAODDTOCEAN.Geolocation_Pointer=modswathvrblsattrs[1][13]
#            
#            #Deep_Blue_Aerosol_Optical_Depth_550_Land_Best_Estimate 
#            modisAODDBLAND=ncmodis.createVariable('/mod04/Data_Fields/Deep_Blue_Aerosol_Optical_Depth_550_Land_Best_Estimate',np.float32,('Cell_Along_Swath','Cell_Across_Swath'),
#                                            fill_value=modswathvrblsattrs[8][14],zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
#            modisAODDBLAND[:]=modswath[14]
#            
#            #Set variable attributes
#            modisAODDBLAND.valid_range=modswathvrblsattrs[4][14]
#            modisAODDBLAND.long_name=modswathvrblsattrs[2][14]
#            modisAODDBLAND.units=modswathvrblsattrs[3][14]
#            modisAODDBLAND.Parameter_Type=modswathvrblsattrs[5][14]
#            modisAODDBLAND.Cell_Along_Swath_Sampling=modswathvrblsattrs[6][14]
#            modisAODDBLAND.Cell_Across_Swath_Sampling=modswathvrblsattrs[7][14]
#            modisAODDBLAND.Geolocation_Pointer=modswathvrblsattrs[1][14]
#            
#            #Corrected_Optical_Depth_Land at 550nm 
#            modisAODDTLAND=ncmodis.createVariable('/mod04/Data_Fields/Corrected_Optical_Depth_Land',np.float32,('Cell_Along_Swath','Cell_Across_Swath'),
#                                            fill_value=modswathvrblsattrs[8][14],zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
#            modisAODDTLAND[:]=modswath[15]
#            
#            #Set variable attributes
#            modisAODDTLAND.valid_range=modswathvrblsattrs[4][13]
#            modisAODDTLAND.long_name='Retrieved  AOT 0.55   micron'
#            modisAODDTLAND.units='None'
#            modisAODDTLAND.Parameter_Type='Output'
#            modisAODDTLAND.Cell_Along_Swath_Sampling=modswathvrblsattrs[6][14]
#            modisAODDTLAND.Cell_Across_Swath_Sampling=modswathvrblsattrs[7][14]
#            modisAODDTLAND.Geolocation_Pointer=modswathvrblsattrs[1][14]
#            
#            #If you want to add variable continue from 16!!!!!!!!

            ncmodis.close()
            
        except Exception as ex:
            
            print ("An error occured with this file: ", filename, ex)
            #f.write("An error occured with this file: " + filename + '\n')
            #f.flush()
    
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

# %%
