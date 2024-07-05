import numpy as np
import firegrib as fg
import os 
import sys
import matplotlib.pyplot as plt
import multiprocessing as mp
import datetime as dt
import numpy.ma as ma
import pickle
from statistics import linear_regression

def readGFASField(species,hour,hourX,spuriousFile=None, region='global', landcover=None,Rerun=True):
    diff = ((hourX-hour).days + 1)
    path = '/xnilu_wrk/users/dgho/gfas/archive/mirrorMARS/0001/'
    if os.path.isfile(species+'field' + ".pickle") and (not Rerun):
        with open(species+'field' + ".pickle",'rb') as fp:
            FinalField = pickle.load(fp)
        f = fg.Field(os.path.join(path,hour.strftime('%Y%m%d'),hour.strftime('%H%M'),"%d.grb.gz"%fg.PARAMID[species]))
    else:
        FinalField = np.zeros((1800,3600))
        while hour <= hourX:
            print(hour)
            f = fg.Field(os.path.join(path,hour.strftime('%Y%m%d'),hour.strftime('%H%M'),"%d.grb.gz"%fg.PARAMID[species]))
            FinalField += f.val
            hour += dt.timedelta(days=1)
        with open(species+'field' + ".pickle",'wb') as fp:
            pickle.dump(FinalField,fp)

    if spuriousFile != None:
        spu = fg.Field(spuriousFile)
        FinalField[spu.val>0] = 0
    
    f.val = (FinalField)/diff
    f.coarser(10)
    maximum = np.max(FinalField)
    f.plot(outFile="TotalFieldPlots/GFASNEW_%s.png"%species,minValue=1e-30,logScale=True,countries='black',cmap='turbo',title="GFAS-NEIVA 2003-2020 %s"%species)


#Read GFED field and sum fields from hour to hourX
def ReadTotalGFEDField(species,hour,hourX, region='global', landcover=None,Rerun=True):

    diff = ((hourX-hour).days + 1)
    path = 'GFED_DM/archive2/'
    
    pickleName = species + 'GFED' + 'Totalfield'+'.pickle'
    if os.path.isfile(pickleName) and (not Rerun):
        with open(pickleName,'rb') as fp:
            FinalField = pickle.load(fp)
        f = fg.Field(os.path.join(path,hour.strftime('%Y%m%d'),"%s.grb"%species))
    else:
        FinalField = np.zeros((720,1440))
        while hour <= hourX:
            print(hour)
            f = fg.Field(os.path.join(path,hour.strftime('%Y%m%d'),"%s.grb"%species))
            FinalField += f.val
            hour += dt.timedelta(days=1)
        with open(pickleName,'wb') as fp:
            pickle.dump(FinalField,fp)
    f.val = (FinalField)/diff
    f.coarser(4)
    maximum = np.max(FinalField)
    f.plot(outFile="TotalFieldPlots/GFED_%s.png"%species,minValue=1e-30,logScale=True,countries='black',cmap='turbo',title="GFED5 2003-2020 %s"%species)



#Read COP field and sum fields from hour to hourX
def ReadTotalCOPField(species,hour,hourX, region='global', landcover=None,Rerun=True):
    diff = ((hourX-hour).days + 1)
    pickleName = species + 'COP' + 'Totalfield'+'.pickle'
    if os.path.isfile(pickleName) and (not Rerun):
        with open(pickleName,'rb') as fp:
            FinalField = pickle.load(fp)
        f = fg.Field("fire_biomes_with_peat_percentcover_2018_0p1deg_v1p0_compatibalised2.grb2")
    else:
        path = 'copernicus'
        FinalField = np.zeros((1800,3600))
        with open(os.path.join(path,"%s.grb"%species),'rb') as fp:
            while hour <= hourX:
                print(hour)
                f = fg.Field(fp)
                FinalField += f.val
                hour += dt.timedelta(days=1)
        with open(pickleName,'wb') as fp:
            pickle.dump(FinalField,fp)
    f.val = (FinalField)/diff
    maximum = np.max(FinalField)
    f.coarser(10)
    f.plot(outFile="TotalFieldPlots/COP_%s.png"%species,minValue=1e-30,logScale=True,countries='black',cmap='turbo',title="GFASv1.2 2003-2020 %s"%species)





names = ['C3H6O','NH3','BC','CO2','CO','C2H6S','C2H6','CH2O','C5H8','CH4','NOx','N2O','NMHC','OC','C3H8','C3H6','SO2','DM']

hour = dt.datetime(2003,1,2,hour=12)
hourX = dt.datetime(2020,12,31,hour=12)
args = []
for species in names:
    args.append([species,hour,hourX])

for species in names:
    ReadTotalCOPField(species,hour,hourX,Rerun=False)
    ReadTotalGFEDField(species,hour,hourX,Rerun=False)
    readGFASField(species,hour,hourX,Rerun=False)