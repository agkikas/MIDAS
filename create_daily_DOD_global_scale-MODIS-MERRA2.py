#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun May  6 06:22:02 2018

@author: antonis
"""

from netCDF4 import Dataset

import datetime as dt
import pandas as pd
import os
import glob as glob
import numpy as np
import calendar
import numpy.ma as ma
import time
import netCDF4 as ncf

scriptstart=dt.datetime.now()

print (scriptstart)

#Define here the year
year='2002'

#Give the satellite name
satellite='MODIS-TERRA'

#Define here which is the gridresolution that you are working
gridres='GRID_RESOLUTION_0.1'

if gridres == 'GRID_RESOLUTION_1.0':
    
    latdim=180
    londim=360
    
elif gridres == 'GRID_RESOLUTION_0.1':
    
    latdim=1800
    londim=3600

#Give here the description of the output netcdf files
ncfiledescription='MODIS-DOD-through-synergy-of-MODIS-MERRA2'

#Define here which is the keyword that you are searching within the file
keyword='MOD04_L2.A'

#Define the least significant digits and compression level for the variables which will be stored in the netcdf files
lstsigdig=3
complevel=9 #This varies from 1 to 9 (http://unidata.github.io/netcdf4-python/)

#Define here where is the path with Level 2 MODIS DOD swaths
modisdodswaths='/mnt/nas-1/MIDAS-Processing_Data/REPROJECTED-MODIS-SWATHS-SYNERGY-WITH-MERRA2/'+satellite+'/'+year+'/'+gridres+'/'
#/media/antonis/My Passport/MERRA-2/REPROJECTED-MODIS-SWATHS-SYNERGY-WITH-MERRA2/

#List all the modis DOD swath files
mdodswaths=glob.glob(os.path.join(modisdodswaths,'*.nc'))

#Create a folder where the final outputs will be stored
#dailyglobalmodisDODfolder='/media/antonis/My Book/MODIS-OMI-CALIOP/DOD550-GLOBAL-SCALE-'+ncfiledescription+'/'+year+'/'+gridres+'/'
#dailyglobalmodisDODfolder='/home/antonis/Desktop/test-files/MODIS-OMI-CALIOP/DOD550-GLOBAL-SCALE-'+ncfiledescription+'/'+AODTHRE+'/'+year+'/'+gridres+'/'
dailyglobalmodisDODfolder='/mnt/nas-1/MIDAS-Processing_Data/DAILY-MODIS-FILES-THROUGH-SYNERGY-OF-MODIS-MERRA2/'+satellite+'/'+year+'/'+gridres+'/'
#/media/antonis/My Passport/MERRA-2/DAILY-MODIS-FILES-THROUGH-SYNERGY-OF-MODIS-MERRA2

if not os.path.exists(dailyglobalmodisDODfolder):
    os.makedirs(dailyglobalmodisDODfolder)

####################################################################################################################################
#This function is used in order to read the MODIS swath files
#def read_MODIS_swath_files(inpath,infile):
def read_MODIS_swath_files(infile):
    #?create method that hits when problematic files are found
    
    #ncfile=Dataset(os.path.join(inpath,infile),'r')
    ncfile=Dataset(infile,'r')

    #Export coordinates from the ncfile
    lat=ncfile['mod04']['Geolocation_Fields']['Latitude'][:]
    
    lon=ncfile['mod04']['Geolocation_Fields']['Longitude'][:]
    
    #Export the total AOD from the ncfile
    modisaod=ncfile['mod04']['Data_Fields']['AOD_550_Dark_Target_Deep_Blue_Combined'][:]
    
    
    #Export the total, fine and coarse dust AOD from the ncfile
    modisdod=ncfile['mod04']['Data_Fields']['MODIS-TOTAL-Dust-optical-depth-at-550nm'][:]
    modisFINEdod=ncfile['mod04']['Data_Fields']['MODIS-FINE-Dust-optical-depth-at-550nm'][:]
    
    modisCOARSEdod=ncfile['mod04']['Data_Fields']['MODIS-COARSE-Dust-optical-depth-at-550nm'][:]
    
    #Export which MODIS algorithm is applied 
    modisALG=ncfile['mod04']['Data_Fields']['AOD_550_Dark_Target_Deep_Blue_Combined_Algorithm_Flag'][:]
    
    #Export the land-sea flag
   
    modisLANDSEAFLAG=ncfile['mod04']['Data_Fields']['Land_sea_Flag'][:]

    #Export the Angrstrom components
    modisANGOCEAN1=ncfile['mod04']['Data_Fields']['Angstrom_Exponent_1_Ocean'][:]
    modisANGOCEAN2=ncfile['mod04']['Data_Fields']['Angstrom_Exponent_2_Ocean'][:]
    modisANGLAND=ncfile['mod04']['Data_Fields']['Deep_Blue_Angstrom_Exponent_Land'][:]
    
    #Extract the solar and sensor zenith angles
    modisSOLZENANGL=ncfile['mod04']['Data_Fields']['Solar_Zenith'][:]
    
    modisSENZENANGL=ncfile['mod04']['Data_Fields']['Sensor_Zenith'][:]
     
    modisvrbls=ma.array([lat,lon,modisaod,modisdod,modisFINEdod,modisCOARSEdod,modisALG,modisLANDSEAFLAG,
                         modisANGOCEAN1,modisANGOCEAN2,modisANGLAND,modisSOLZENANGL,modisSENZENANGL])
    
    ncfile.close()
    
    """
    0: Latitude
    1: Longitude
    2: MODIS AOD
    3: MODIS TOTAL dust AOD
    4: MODIS FINE dust AOD
    5: MODIS COARSE dust AOD
    6: Which algorithm is applied (DT or DB)
    7: Land-Sea Flag
    8: MODIS Angstrom 1 over ocean (550-860nm)
    9: MODIS Angstrom 2 over ocean (860-2130nm)
    10: MODIS Angstrom over land (Deep Blue product)
    11: Solar zenith angle
    12: Sensor zenith angle
    
    """
    return modisvrbls

####################################################################################################################################
####################################################################################################################################
#This functions is used to find date from each MODIS swath contained in the total list
def find_MODIS_swath_overpass_time(infile):

    y=np.int(infile[:4])
    jd=np.int(infile[4:7])
    hour=np.int(infile[8:10]) 
    minutes=np.int(infile[10:]) 
    
    month=1
    
    while jd - calendar.monthrange(y,month)[1] > 0 and month <= 12:
        jd = jd - calendar.monthrange(y,month)[1]
        month += 1
    
    modisoverpasstime = dt.datetime(y,month,jd,hour,minutes).strftime("%d:%m:%Y %H:%M:%S")
    
    overpasstimemodis=pd.to_datetime(modisoverpasstime,format="%d:%m:%Y %H:%M:%S")
    
    return overpasstimemodis

####################################################################################################################################
#Create a pandas timeseries for the year that you are wroking
startyeardate='1/1/'+year
endyeardate='31/12/'+year
freqyeardate='D'

yeardates=pd.date_range(start=startyeardate,end=endyeardate,freq=freqyeardate)
####################################################################################################################################
####################################################################################################################################
#Create a list with the dates of MODIS swaths
mdodswathdateslist=[]

#Loop in the MODIS swaths filenames in order to extract swath dates
for mdodswath in sorted(mdodswaths):
    
    kwid=mdodswath.find(keyword)
    
    #Find year, date and time from the file
    swathdatestring=mdodswath[kwid+len(keyword):kwid+len(keyword)+12]

    swathdate=find_MODIS_swath_overpass_time(swathdatestring)

    #Append each swathdate to the total list with dates        
    mdodswathdateslist.append(swathdate)

#Convert the dates list to an array
mdodswathdates=np.array(mdodswathdateslist)

#Convert the numpy array in a dataframe
dfdates=pd.DataFrame(mdodswathdates)

dfdates.columns=['Swath_Date'] #Give a name to dfdates dataframe column

dfdates['File_Index']=np.arange(len(dfdates)) #Add a column with the file index

dfdates['Swath_Date']=pd.to_datetime(dfdates['Swath_Date']) #Convert to datetime the Swath_Date column

dfdates.set_index(['Swath_Date'],inplace=True) #Set the dates as index

####################################################################################################################################
for row in np.arange(len(yeardates)):

    #Apply a mask in order to find which files do you want     
    maskdates=(dfdates.index.day==yeardates[row].day) & (dfdates.index.month==yeardates[row].month) & (dfdates.index.year==yeardates[row].year)
    
    #Here is a temporary dataframe which will be used in order to select MODIS swaths from the overall list of files
    daydf=dfdates.loc[maskdates]
       
    #Find the minimum and maximum indices in order to pick the MODIS swath files for each day
    idmin=daydf['File_Index'].min()
    idmax=daydf['File_Index'].max()
    
    print (idmin,idmax)
       
    #Check if idmin or idmax are nans
    if (np.isnan(idmin)==False) & (np.isnan(idmax)==False):
        
        #print idmin,idmax
                   
        #Here selects part of the list with all MODIS swaths that will be used in order to make the calculations for each day
        tmpmdodswaths=sorted(mdodswaths)[idmin:idmax+1]

#####################################################################################################################################             
#        #Here create an empty list where all MODIS swaths with AOD, for each day, will be stored
#        modisaoddaylist=[]
#        
#        #Here create an empty list where all MODIS swaths with DOD, for each day, will be stored
#        modisdoddaylist=[]
#                
#        #Loop in the files list for each day
#        for tmpmdodswath in tmpmdodswaths:
#                        
#            modswath=read_MODIS_swath_files(tmpmdodswath)
#            
#            modisaoddaylist.append(modswath[2])
#               
#            modisdoddaylist.append(modswath[3])
#                   
#        #Convert modisdoddaylist in a masked array
#        modisaoddaytmp=ma.array(modisaoddaylist)
#        modisdoddaytmp=ma.array(modisdoddaylist)                
#              
#        print yeardates[row],modisdoddaytmp.shape
#        
#        #Here is the daily MODIS dust optical depth
#        modisaoddayFIN=np.nanmean(modisaoddaytmp,axis=0)
#        
#        #Here is the daily MODIS dust optical depth
#        modisdoddayFIN=np.nanmean(modisdoddaytmp,axis=0)
##################################################################################################################################### 
#        modisaodFIN=np.nan*np.ones((latdim,londim))
#        
#        modisdodFIN=np.nan*np.ones((latdim,londim))
#        modisdodFINEFIN=np.nan*np.ones((latdim,londim))
#        modisdodCOARSEFIN=np.nan*np.ones((latdim,londim))
#        
#        modisdodALGFIN=np.nan*np.ones((latdim,londim))
#        
#        modisdodLANDSEAFLAGFIN=np.nan*np.ones((latdim,londim))
#        
#        modisdodANGOCEAN1FIN=np.nan*np.ones((latdim,londim))
#        modisdodANGOCEAN2FIN=np.nan*np.ones((latdim,londim))
#        modisdodANGLANDFIN=np.nan*np.ones((latdim,londim))
        
#####################################################################################################################################
        modisaodSUM=np.zeros((latdim,londim))
        
        modisdodSUM=np.zeros((latdim,londim))
        modisdodFINESUM=np.zeros((latdim,londim))
        modisdodCOARSESUM=np.zeros((latdim,londim))
        
        modisdodALGSUM=np.zeros((latdim,londim))
        
        modisdodLANDSEAFLAGSUM=np.zeros((latdim,londim))
        
        modisdodANGOCEAN1SUM=np.zeros((latdim,londim))
        modisdodANGOCEAN2SUM=np.zeros((latdim,londim))
        modisdodANGLANDSUM=np.zeros((latdim,londim))
        
        modisdodSOLZENANGLSUM=np.zeros((latdim,londim))
        modisdodSENZENANGLSUM=np.zeros((latdim,londim))
                
        modisaodCOUNTS=np.zeros((latdim,londim))
        
        modisdodCOUNTS=np.zeros((latdim,londim))
        modisdodFINECOUNTS=np.zeros((latdim,londim))
        modisdodCOARSECOUNTS=np.zeros((latdim,londim))
        
        modisdodALGCOUNTS=np.zeros((latdim,londim))
        
        modisdodLANDSEAFLAGCOUNTS=np.zeros((latdim,londim))
        
        modisdodANGOCEAN1COUNTS=np.zeros((latdim,londim))
        modisdodANGOCEAN2COUNTS=np.zeros((latdim,londim))
        modisdodANGLANDCOUNTS=np.zeros((latdim,londim))
        
        modisdodSOLZENANGLCOUNTS=np.zeros((latdim,londim))
        modisdodSENZENANGLCOUNTS=np.zeros((latdim,londim))
                
        #Loop in the files list for each day
        for tmpmdodswath in tmpmdodswaths:
            
            try: 
            
                modswath=read_MODIS_swath_files(tmpmdodswath)
                    
                maskelemAOD=modswath[2].mask
                maskelemDOD=modswath[3].mask  
                maskelemDODfine=modswath[4].mask
                maskelemDODcoarse=modswath[5].mask
                maskelemALG=modswath[6].mask
                maskelemLANDSEAFLAG=modswath[7].mask
                maskelemANGOCEAN1=modswath[8].mask
                maskelemANGOCEAN2=modswath[9].mask
                maskelemANGLAND=modswath[10].mask
                maskelemSOLZENANGL=modswath[11].mask
                maskelemSENZENANGL=modswath[12].mask
                        
                maskelemAODbin=maskelemAOD.astype(int)
                maskelemDODbin=maskelemDOD.astype(int)
                maskelemDODfinebin=maskelemDODfine.astype(int)
                maskelemDODcoarsebin=maskelemDODcoarse.astype(int)
                maskelemALGbin=maskelemALG.astype(int)
                maskelemLANDSEAFLAGbin=maskelemLANDSEAFLAG.astype(int)
                maskelemANGOCEAN1bin=maskelemANGOCEAN1.astype(int)
                maskelemANGOCEAN2bin=maskelemANGOCEAN2.astype(int)
                maskelemANGLANDbin=maskelemANGLAND.astype(int)
                maskelemSOLZENANGLbin=maskelemSOLZENANGL.astype(int)
                maskelemSENZENANGLbin=maskelemSOLZENANGL.astype(int)
            
                nonmaskidxsAOD=np.where(maskelemAODbin==0)
                nonmaskidxsDOD=np.where(maskelemDODbin==0)            
                nonmaskidxsDODfine=np.where(maskelemDODfinebin==0)
                nonmaskidxsDODcoarse=np.where(maskelemDODcoarsebin==0)
                nonmaskidxsALG=np.where(maskelemALGbin==0)
                nonmaskidxsLANDSEAFLAG=np.where(maskelemLANDSEAFLAGbin==0)
                nonmaskidxsANGOCEAN1=np.where(maskelemANGOCEAN1bin==0)
                nonmaskidxsANGOCEAN2=np.where(maskelemANGOCEAN2bin==0)
                nonmaskidxsANGLAND=np.where(maskelemANGLANDbin==0)
                nonmaskidxsSOLZENANGL=np.where(maskelemSOLZENANGLbin==0)
                nonmaskidxsSENZENANGL=np.where(maskelemSENZENANGLbin==0)
            
#####################################################################################################################################                
#            modisaodFIN[nonmaskidxsAOD[0],nonmaskidxsAOD[1]]=modswath[2][nonmaskidxsAOD[0],nonmaskidxsAOD[1]]
#            
#            modisdodFIN[nonmaskidxsDOD[0],nonmaskidxsDOD[1]]=modswath[3][nonmaskidxsDOD[0],nonmaskidxsDOD[1]]
#            modisdodFINEFIN[nonmaskidxsDODfine[0],nonmaskidxsDODfine[1]]=modswath[4][nonmaskidxsDODfine[0],nonmaskidxsDODfine[1]]
#            modisdodCOARSEFIN[nonmaskidxsDODcoarse[0],nonmaskidxsDODcoarse[1]]=modswath[5][nonmaskidxsDODcoarse[0],nonmaskidxsDODcoarse[1]]
#            
#            modisdodALGFIN[nonmaskidxsALG[0],nonmaskidxsALG[1]]=modswath[6][nonmaskidxsALG[0],nonmaskidxsALG[1]]
#            
#            modisdodLANDSEAFLAGFIN[nonmaskidxsLANDSEAFLAG[0],nonmaskidxsLANDSEAFLAG[1]]=modswath[7][nonmaskidxsLANDSEAFLAG[0],nonmaskidxsLANDSEAFLAG[1]]
#            
#            modisdodANGOCEAN1FIN[nonmaskidxsANGOCEAN1[0],nonmaskidxsANGOCEAN1[1]]=modswath[8][nonmaskidxsANGOCEAN1[0],nonmaskidxsANGOCEAN1[1]]
#            modisdodANGOCEAN2FIN[nonmaskidxsANGOCEAN2[0],nonmaskidxsANGOCEAN2[1]]=modswath[9][nonmaskidxsANGOCEAN2[0],nonmaskidxsANGOCEAN2[1]]
#            modisdodANGLANDFIN[nonmaskidxsANGLAND[0],nonmaskidxsANGLAND[1]]=modswath[10][nonmaskidxsANGLAND[0],nonmaskidxsANGLAND[1]]
#####################################################################################################################################
                #Sum the values
                modisaodSUM[nonmaskidxsAOD[0],nonmaskidxsAOD[1]]+=modswath[2][nonmaskidxsAOD[0],nonmaskidxsAOD[1]]
                modisdodSUM[nonmaskidxsDOD[0],nonmaskidxsDOD[1]]+=modswath[3][nonmaskidxsDOD[0],nonmaskidxsDOD[1]]
                modisdodFINESUM[nonmaskidxsDODfine[0],nonmaskidxsDODfine[1]]+=modswath[4][nonmaskidxsDODfine[0],nonmaskidxsDODfine[1]]
                modisdodCOARSESUM[nonmaskidxsDODcoarse[0],nonmaskidxsDODcoarse[1]]+=modswath[5][nonmaskidxsDODcoarse[0],nonmaskidxsDODcoarse[1]]
                modisdodALGSUM[nonmaskidxsALG[0],nonmaskidxsALG[1]]+=modswath[6][nonmaskidxsALG[0],nonmaskidxsALG[1]]
                modisdodLANDSEAFLAGSUM[nonmaskidxsLANDSEAFLAG[0],nonmaskidxsLANDSEAFLAG[1]]+=modswath[7][nonmaskidxsLANDSEAFLAG[0],nonmaskidxsLANDSEAFLAG[1]]
                modisdodANGOCEAN1SUM[nonmaskidxsANGOCEAN1[0],nonmaskidxsANGOCEAN1[1]]+=modswath[8][nonmaskidxsANGOCEAN1[0],nonmaskidxsANGOCEAN1[1]]
                modisdodANGOCEAN2SUM[nonmaskidxsANGOCEAN2[0],nonmaskidxsANGOCEAN2[1]]+=modswath[9][nonmaskidxsANGOCEAN2[0],nonmaskidxsANGOCEAN2[1]]
                modisdodANGLANDSUM[nonmaskidxsANGLAND[0],nonmaskidxsANGLAND[1]]+=modswath[10][nonmaskidxsANGLAND[0],nonmaskidxsANGLAND[1]]
                modisdodSOLZENANGLSUM[nonmaskidxsSOLZENANGL[0],nonmaskidxsSOLZENANGL[1]]+=modswath[11][nonmaskidxsSOLZENANGL[0],nonmaskidxsSOLZENANGL[1]]
                modisdodSENZENANGLSUM[nonmaskidxsSENZENANGL[0],nonmaskidxsSENZENANGL[1]]+=modswath[12][nonmaskidxsSENZENANGL[0],nonmaskidxsSENZENANGL[1]]
                        
                #Sum the counts
                modisaodCOUNTS[nonmaskidxsAOD[0],nonmaskidxsAOD[1]]+=1
                modisdodCOUNTS[nonmaskidxsDOD[0],nonmaskidxsDOD[1]]+=1
                modisdodFINECOUNTS[nonmaskidxsDODfine[0],nonmaskidxsDODfine[1]]+=1
                modisdodCOARSECOUNTS[nonmaskidxsDODcoarse[0],nonmaskidxsDODcoarse[1]]+=1
                modisdodALGCOUNTS[nonmaskidxsALG[0],nonmaskidxsALG[1]]+=1
                modisdodLANDSEAFLAGCOUNTS[nonmaskidxsLANDSEAFLAG[0],nonmaskidxsLANDSEAFLAG[1]]+=1
                modisdodANGOCEAN1COUNTS[nonmaskidxsANGOCEAN1[0],nonmaskidxsANGOCEAN1[1]]+=1
                modisdodANGOCEAN2COUNTS[nonmaskidxsANGOCEAN2[0],nonmaskidxsANGOCEAN2[1]]+=1
                modisdodANGLANDCOUNTS[nonmaskidxsANGLAND[0],nonmaskidxsANGLAND[1]]+=1            
                modisdodSOLZENANGLCOUNTS[nonmaskidxsSOLZENANGL[0],nonmaskidxsSOLZENANGL[1]]+=1
                modisdodSENZENANGLCOUNTS[nonmaskidxsSENZENANGL[0],nonmaskidxsSENZENANGL[1]]+=1
                
            except:
                
                print ('Something is missing in the '+tmpmdodswath+'!!!!')

        #Here is the calculation of the mean value. This is done for the overlapping areas between MODIS swaths       
        modisaodFIN=modisaodSUM/modisaodCOUNTS
        modisdodFIN=modisdodSUM/modisdodCOUNTS
        modisdodFINEFIN=modisdodFINESUM/modisdodFINECOUNTS
        modisdodCOARSEFIN=modisdodCOARSESUM/modisdodCOARSECOUNTS
        modisdodALGFIN=modisdodALGSUM/modisdodALGCOUNTS
        modisdodLANDSEAFLAGFIN=modisdodLANDSEAFLAGSUM/modisdodLANDSEAFLAGCOUNTS
        modisdodANGOCEAN1FIN=modisdodANGOCEAN1SUM/modisdodANGOCEAN1COUNTS
        modisdodANGOCEAN2FIN=modisdodANGOCEAN2SUM/modisdodANGOCEAN2COUNTS
        modisdodANGLANDFIN=modisdodANGLANDSUM/modisdodANGLANDCOUNTS
        modisdodSOLZENANGLFIN=modisdodSOLZENANGLSUM/modisdodSOLZENANGLCOUNTS
        modisdodSENZENANGLFIN=modisdodSENZENANGLSUM/modisdodSENZENANGLCOUNTS
          
#####################################################################################################################################         
        modisaoddayFIN=ma.masked_invalid(modisaodFIN)        
        modisdoddayFIN=ma.masked_invalid(modisdodFIN)
        modisdaydodFINEFIN=ma.masked_invalid(modisdodFINEFIN)
        modisdaydodCOARSEFIN=ma.masked_invalid(modisdodCOARSEFIN)        
        modisdaydodALGFIN=ma.masked_invalid(modisdodALGFIN)        
        modisdaydodLANDSEAFLAGFIN=ma.masked_invalid(modisdodLANDSEAFLAGFIN)
        modisdaydodANGOCEAN1FIN=ma.masked_invalid(modisdodANGOCEAN1FIN)
        modisdaydodANGOCEAN2FIN=ma.masked_invalid(modisdodANGOCEAN2FIN)
        modisdaydodANGLANDFIN=ma.masked_invalid(modisdodANGLANDFIN)
        modisdaydodSOLZENANGLFIN=ma.masked_invalid(modisdodSOLZENANGLFIN)
        modisdaydodSENZENANGLFIN=ma.masked_invalid(modisdodSENZENANGLFIN)
####################################################################################################################################                                 
        #Create the nectdf file
        #outncfile=satellite+'_'+'AOD-and-DOD-'+gridres+'-'+yeardates[row].strftime('%d_%b_%Y')+'.nc'
        outncfile=satellite+'_'+'AOD-and-DOD-'+gridres+'-'+yeardates[row].strftime('%Y%m%d')+'.nc'
        
        ncmodis=Dataset(os.path.join(dailyglobalmodisDODfolder,outncfile),'w',format='NETCDF4')
        
        #Create Global Attributes
        ncmodis.description=ncfiledescription
        ncmodis.history='Created:  '+time.ctime(time.time())
        ncmodis.source='netCDF4 python module tutorial'
        ncmodis.Conventions = 'CF-1.6'
        ncmodis.time='days since '+yeardates[row].strftime('%Y-%m-%d %H:%M:%S')
        
        #Create groups
        #maingrp=ncmodis.createGroup('mod04')
        #datafieldsgrp=ncmodis.createGroup('/mod04/Data_Fields')
        #geolocgrp=ncmodis.createGroup('/mod04/Geolocation_Fields')
        
        latitude=ncmodis.createDimension('Latitude',modswath[0].shape[0])
        longitude=ncmodis.createDimension('Longitude',modswath[1].shape[1])
        ttime=ncmodis.createDimension('Time',None)
               
        #Create the variables which will be stored in the netcdf file
        #Latitude
        #Latitude=ncmodis.createVariable('/mod04/Geolocation_Fields/Latitude',np.float32,('Latitude','Longitude'),
        #                                fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
        
        Latitude=ncmodis.createVariable('Latitude',np.float32,('Latitude','Longitude'),
                                        fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
        
        Latitude[:]=modswath[0][::-1,:]
        
        Latitude.valid_range=(-90.,90.)
        Latitude.long_name='Geodetic Latitude'
        Latitude.units='degrees_north'
        
        #Longitude
        #Longitude=ncmodis.createVariable('/mod04/Geolocation_Fields/Longitude',np.float32,('Latitude','Longitude'),
        #                                 fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
        
        Longitude=ncmodis.createVariable('Longitude',np.float32,('Latitude','Longitude'),
                                         fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
        
        Longitude[:]=modswath[1][::-1,:]
        
        Longitude.valid_range=(-180.,180.)
        Longitude.long_name='Geodetic Longitude'
        Longitude.units='degrees_east'
               
        #Aerosol optical depth
        #modaod=ncmodis.createVariable('/mod04/Data_Fields/'+'AOD_550_Dark_Target_Deep_Blue_Combined',np.float32,('Latitude','Longitude'),
        #                                 fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
        
        modaod=ncmodis.createVariable('Aod_550_dark_target_deep_blue_combined',np.float32,('Latitude','Longitude'),
                                      fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)       
        
        
        modaod[:]=modisaoddayFIN[::-1,:]
        
        modaod.long_name='AOD_550_Dark_Target_Deep_Blue_Combined_(masked_CF>0.8_and_isolated_retrievals)'
        modaod.coordinates="Latitude Longitude"
        modaod.units='None'
        
        #Dust optical depth
        #moddod=ncmodis.createVariable('/mod04/Data_Fields/'+'MODIS-TOTAL-Dust-optical-depth-at-550nm',np.float32,('Latitude','Longitude'),
        #                                 fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
        
        moddod=ncmodis.createVariable('Modis-total-dust-optical-depth-at-550nm',np.float32,('Latitude','Longitude'),
                                      fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
        
        moddod[:]=modisdoddayFIN[::-1,:]
        
        moddod.long_name='MODIS-TOTAL-Dust-optical-depth-at-550nm_through_MERRA2'
        moddod.coordinates="Latitude Longitude"
        moddod.units='None'
        
        #Fine Dust optical depth
        #modFINEdod=ncmodis.createVariable('/mod04/Data_Fields/'+'MODIS-FINE-Dust-optical-depth-at-550nm',np.float32,('Latitude','Longitude'),
        #                                 fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
        
        modFINEdod=ncmodis.createVariable('Modis-fine-dust-optical-depth-at-550nm',np.float32,('Latitude','Longitude'),
                                          fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
                
        modFINEdod[:]=modisdaydodFINEFIN[::-1,:]
        
        modFINEdod.long_name='MODIS-FINE-Dust-optical-depth-at-550nm_through_MERRA2'
        modFINEdod.coordinates="Latitude Longitude"
        modFINEdod.units='None'
        
        #Coarse Dust optical depth
        #modCOARSEdod=ncmodis.createVariable('/mod04/Data_Fields/'+'MODIS-COARSE-Dust-optical-depth-at-550nm',np.float32,('Latitude','Longitude'),
        #                                 fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
        
        modCOARSEdod=ncmodis.createVariable('Modis-coarse-dust-optical-depth-at-550nm',np.float32,('Latitude','Longitude'),
                                            fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
             
        modCOARSEdod[:]=modisdaydodCOARSEFIN[::-1,:]
        
        modCOARSEdod.long_name='MODIS-COARSE-Dust-optical-depth-at-550nm_through_MERRA2'
        modCOARSEdod.coordinates="Latitude Longitude"
        modCOARSEdod.units='None'
        
        #MODIS algorithm
        #modalgorithm=ncmodis.createVariable('/mod04/Data_Fields/'+'AOD_550_Dark_Target_Deep_Blue_Combined_Algorithm_Flag',np.float32,('Latitude','Longitude'),
        #                                 fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
        
        modalgorithm=ncmodis.createVariable('Aod_550_dark_target_deep_blue_combined_algorithm_flag',np.float32,('Latitude','Longitude'),
                                            fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
               
        modalgorithm[:]=modisdaydodALGFIN[::-1,:]
        
        modalgorithm.long_name='Combined Dark Target, Deep Blue AOT at 0.55 micron Algorithm Flag (0=Dark Target, 1=Deep Blue, 2=Mixed)'
        modalgorithm.coordinates="Latitude Longitude"
        modalgorithm.units='None'
        
        #MODIS Land-Sea flag
        #modlandseaflag=ncmodis.createVariable('/mod04/Data_Fields/'+'Land_sea_Flag',np.float32,('Latitude','Longitude'),
        #                                 fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
        
        modlandseaflag=ncmodis.createVariable('Lsf',np.float32,('Latitude','Longitude'),
                                         fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)       
        
        
        modlandseaflag[:]=modisdaydodLANDSEAFLAGFIN[::-1,:]
        
        modlandseaflag.long_name='Land_sea_Flag(based on MOD03 Landsea mask 0 = Ocean, 1 = Land and Ephemeral water 2 =Coastal)'
        modlandseaflag.coordinates="Latitude Longitude"
        modlandseaflag.units='None'
        
        #Angstrom exponent ocean 1
        #modangocean1=ncmodis.createVariable('/mod04/Data_Fields/'+'Angstrom_Exponent_1_Ocean',np.float32,('Latitude','Longitude'),
        #                                 fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)

        modangocean1=ncmodis.createVariable('Angstrom_exponent_1_ocean',np.float32,('Latitude','Longitude'),
                                            fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
                  
        modangocean1[:]=modisdaydodANGOCEAN1FIN[::-1,:]
        
        modangocean1.long_name='Calculated Angstrom Exponent for 0.55 vs 0.86 micron  for Average Solution'
        modangocean1.coordinates="Latitude Longitude"
        modangocean1.units='None'
        
        #Angstrom exponent ocean 2
        #modangocean2=ncmodis.createVariable('/mod04/Data_Fields/'+'Angstrom_Exponent_2_Ocean',np.float32,('Latitude','Longitude'),
        #                                 fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
        
        modangocean2=ncmodis.createVariable('Angstrom_exponent_2_ocean',np.float32,('Latitude','Longitude'),
                                            fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
        
        modangocean2[:]=modisdaydodANGOCEAN2FIN[::-1,:]
        
        modangocean2.long_name='Calculated Angstrom Exponent for 0.86 vs 2.13 micron for Average Solution'
        modangocean2.coordinates="Latitude Longitude"
        modangocean2.units='None'
        
        #Angstrom exponent land
        #modangland=ncmodis.createVariable('/mod04/Data_Fields/'+'Deep_Blue_Angstrom_Exponent_Land',np.float32,('Latitude','Longitude'),
        #                                 fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
        
        modangland=ncmodis.createVariable('Deep_blue_angstrom_exponent_land',np.float32,('Latitude','Longitude'),
                                          fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
        
        modangland[:]=modisdaydodANGLANDFIN[::-1,:]
        
        modangland.long_name='Deep Blue Angstrom Exponent (412-470nm)'
        modangland.coordinates="Latitude Longitude"
        modangland.units='None'
        
        #Solar zenith angle
        #modsolzenangl=ncmodis.createVariable('/mod04/Data_Fields/'+'Solar_Zenith',np.float32,('Latitude','Longitude'),
        #                                 fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
        
        modsolzenangl=ncmodis.createVariable('Sza5',np.float32,('Latitude','Longitude'),
                                             fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
        
        modsolzenangl[:]=modisdaydodSOLZENANGLFIN[::-1,:]
        
        modsolzenangl.long_name='Solar Zenith Angle, Cell to Sun'
        modsolzenangl.coordinates="Latitude Longitude"
        modsolzenangl.units='Degrees'
        
        #Sensor zenith angle
        #modsenzenangl=ncmodis.createVariable('/mod04/Data_Fields/'+'Sensor_Zenith',np.float32,('Latitude','Longitude'),
        #                                 fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
        
        modsenzenangl=ncmodis.createVariable('Sensor_zenith',np.float32,('Latitude','Longitude'),
                                             fill_value=np.nan,zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
        
        modsenzenangl[:]=modisdaydodSENZENANGLFIN[::-1,:]
        
        modsenzenangl.long_name='Sensor Zenith Angle, Cell to Sensor'
        modsenzenangl.coordinates="Latitude Longitude"
        modsenzenangl.units='Degrees'
        
        #Time
        modtime=ncmodis.createVariable('Time',np.float32,('Time',),zlib=True,least_significant_digit=lstsigdig,complevel=complevel)
               
        modtime[:]=ncf.date2num(yeardates[row],units='days since '+yeardates[row].strftime('%Y-%m-%d %H:%M:%S'))
        
        #modtime[:]=yeardates[row].strftime('%Y-%m-%d %H:%M:%S')
        
        modtime.long_name='Time'
        modtime.units='days since '+yeardates[row].strftime('%Y-%m-%d %H:%M:%S')
        modtime.calendar='julian'

        ncmodis.close()

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
