import numpy as np
import firegrib as fg
import os 
import sys
import matplotlib.pyplot as plt
import multiprocessing as mp
import datetime as dt
import numpy.ma as ma
import tracemalloc
from copy import deepcopy
import tracemalloc
import pickle
from statistics import linear_regression

landcovers = {'tropical_bare': 1970, 'temperate_bare': 2970, 'boreal_bare': 3970,
             'tropical_sparse': 1100, 'temperate_sparse': 2100, 'boreal_sparse': 3100,
             'high_woodland_savanna': 1080, 'low_woodland_savanna': 1090,
             'tropical_grassland': 1060, 'temperate_grassland': 2060, 'boreal_grassland': 3060,
             'temperate_shrubland': 2070, 'boreal_shrubland': 3070,
             'tropical_forest': 1020, 'temperate_broadleaf_evergreen': 2020, 'boreal_broadleaf_evergreen': 3020,
             'temperate_broadleaf_deciduous': 2030, 'boreal_broadleaf_deciduous': 3030,
             'temperate_needleleaf_evergreen': 2040, 'boreal_needleleaf_evergreen': 3040,
             'temperate_needleleaf_deciduous': 2050, 'boreal_needleleaf_deciduous': 3050,
             'tropical_cropland': 1010, 'temperate_cropland': 2010, 'boreal_cropland':3010}

dlc = fg.Field('fire_biomes_with_peat_percentcover_2018_0p1deg_v1p0_compatibalised.grb2')
dlcmasks = {}
peatmask = ((dlc.val % 10) == 9)
mixed_peatmask = ((dlc.val % 10) > 1)
for land in landcovers:
    dlcmasks[land] = np.logical_and((dlc.val // 10) == (landcovers[land] // 10), np.logical_not(mixed_peatmask))
dlcmasks['peat'] = peatmask
landcovers['peat'] = 'R'

def fetch_GFAS_frp(hour):
    path = '/xnilu_wrk/users/dgho/gfas/archive/mirrorMARS/0001/'

    spuriousName = "GFAS_MCD14ML_VNP14ML" + '_' +  hour.strftime('%Y') + '.grb2'       
    spu = fg.Field(os.path.join('SpuriousMap',spuriousName))
    filename = str(fg.PARAMID['FRP'])+".grb.gz"

    f = fg.Field(os.path.join(path,hour.strftime('%Y%m%d'),hour.strftime('%H%M'),filename))
    f.val[spu.val>0] = 0
    return f

def fetch_GFED_DM(hour):
    path = '/xnilu_wrk/users/dgho/gfas/satDataAna/GFED_DM/archive2/'

    filename = 'DM.grb'
    f = fg.Field(os.path.join(path,hour.strftime('%Y%m%d'),filename))

    return f


#This is preprocessing for scattering all points in a landcover type. it removes points with zero DM and FRP for each and cover mask and compresses the data.
def remove_zeros(hour,hourX):
    numberNonZero = 0
    savedFRP = {}
    savedDM = {}
    for land in landcovers:
        savedFRP[land] = np.array([0])
        savedDM[land] = np.array([0])
    while hour <= hourX:
        FRP = fetch_GFAS_frp(hour).val
        DM = fetch_GFED_DM(hour).val
        print(hour)
        

        for land in landcovers:
            mask = np.logical_or(np.logical_and(FRP < 1e-7,DM < 1e-13),np.logical_not(dlcmasks[land]))
            FRP_m = ma.masked_array(FRP,mask)
            DM_m = ma.masked_array(DM,mask)
            savedFRP[land] = np.append(savedFRP[land],FRP_m.compressed())
            savedDM[land] = np.append(savedDM[land],DM_m.compressed())
 
        
        hour += dt.timedelta(days=1)
    hour -= dt.timedelta(days=1)
    for land in landcovers:
        path = "SavedPairs/%s"%(land)
        pairs = np.array([savedFRP[land],savedDM[land]])
        if not os.path.exists(path):
            os.makedirs(path)
        with open(os.path.join(path,hour.strftime("%Y.npy")),"wb") as fp:
            np.save(fp, pairs)

def Calculate_timestep(hour,deltat):
    savedFRP = {}
    compressedDLC = {}
    print(hour)

    iterator = dt.datetime(hour.year, hour.month, hour.day, hour=12)

    tempFRP = np.zeros((1800,3600))
    tempDM = np.zeros((720,1440))
    while iterator < hour + dt.timedelta(days=deltat):
        print(iterator)
        DM = fetch_GFED_DM(iterator)
        FRP = fetch_GFAS_frp(iterator)
        tempDM += DM.val
        tempFRP += FRP.val
        iterator += dt.timedelta(days=1)
    
    DM.val = tempDM/deltat
    FRP.val = tempFRP/deltat
    FRP2 = fg.Field('fire_biomes_with_peat_percentcover_2018_0p1deg_v1p0_compatibalised.grb2')
    FRP2.val = FRP.val
    FRP2.coarser(5)
    DM.coarser(2)
    mask = np.logical_or(FRP2.val > 0, DM.val > 0)
    FRPlc = fg.Field('fire_biomes_with_peat_percentcover_2018_0p1deg_v1p0_compatibalised.grb2')
    DLCmaskField = fg.Field('fire_biomes_with_peat_percentcover_2018_0p1deg_v1p0_compatibalised.grb2')
    for land in landcovers:
        DLCmaskField.val = dlcmasks[land]
        DLCmaskField.lon = FRP.lon
        DLCmaskField.lat = FRP.lat
        DLCmaskField.area = FRP.area
        DLCmaskField.coarser(5,algo='sum')
        DLCmaskField.val = DLCmaskField.val > 20
        FRPlc.val = FRP.val*dlcmasks[land]
        FRPlc.lon = FRP.lon
        FRPlc.lat = FRP.lat
        FRPlc.area = FRP.area
        FRPlc.coarser(5)
        FRPcompressed = FRPlc.val[mask]
        compressedDLC[land] = DLCmaskField.val[mask]

        savedFRP[land] = FRPcompressed
        print(land)
    
    savedDM = DM.val[mask]

    return [savedFRP,savedDM,compressedDLC]


def SetupRegression(hour,hourX,deltat,savefile="DM_Regression_Equation_System.pickle",PARALLELL=False):
    """
    Removes (hourX-hour)mod(deltat) days from the end of the timeseries to ensure uniform time intervals
    """

    diff = (hourX-hour).days
    timesteps = [(hour + dt.timedelta(days=(i*deltat)),deltat) for i in range(diff // deltat)]
    print(len(timesteps))
    if PARALLELL:
        pool = mp.Pool()
        results = pool.starmap(Calculate_timestep,timesteps)
        pool.close()
        pool.join()
        
        print("done")
    else:
        results = []
        for i in range(len(timesteps)):
            a = Calculate_timestep(timesteps[i][0],deltat)
            print("watsup")
            results.append(a)
    
    savedFRP = {}
    savedDM = np.array([0])
    compressedDLC = {}
    for land in landcovers:
        savedFRP[land] = np.array([0])
        compressedDLC[land] = np.array([True])
    for i in range(len(results)):
        savedDM = np.append(savedDM,results[i][1])
        for land in landcovers:
            savedFRP[land] = np.append(savedFRP[land],results[i][0][land])
            compressedDLC[land] = np.append(compressedDLC[land],results[i][2][land])

    with open(savefile,"wb") as fp:
        pickle.dump([savedFRP,savedDM,compressedDLC],fp,protocol=pickle.HIGHEST_PROTOCOL)


    return savedFRP,savedDM,compressedDLC

def Cf_Regression_Plot(hour,hourX,deltat,savefile="DM_Regression_Equation_System.pickle"):
    if os.path.exists(savefile):
        with open(savefile,'rb') as fp:
            FRPl, DM, compressedDLC = pickle.load(fp)


    else:
        FRPl,DM,compressedDLC = SetupRegression(hour,hourX,deltat,savefile)
    
    landcover_types = []
    FRP = []
    i = 0
    for land in landcovers:
        landcover_types.append(land)
        FRP.append(FRPl[land])
    FRP = np.array(FRP)
    FRPs = np.sort(FRP)
    firstNonZero = (FRPs != 0).argmax(axis=-1)
    #Remove highest 5% for each landcover
    max_values = []
    maskHighest = np.ones(len(FRP[1]),dtype=bool)
    for i in range(len(FRP)):
 #       print(len(FRPs[1]), firstNonZero[i] + int((len(FRPs[1])-firstNonZero[i])*0.95))
        max_values.append(FRPs[i][firstNonZero[i] + int((len(FRPs[1])-firstNonZero[i])*0.95)])
        maskHighest = np.logical_and(maskHighest,FRP[i] <max_values[i])
    
    #Remove low DM values
    maskHighest = np.logical_and(maskHighest, DM > 0.5e-10)
    print(len(FRP[0]))
    FRPreg = FRP[:,maskHighest]
    print(len(FRP[0]))
    DMreg = DM[maskHighest]

    #Regression
    FRPreg = np.transpose(FRPreg)
 #   print(FRP.shape,DM.shape)
    Cf = np.linalg.lstsq(FRPreg,DMreg)
#    print(np.sum(np.dot(FRP,Cf[0])),np.sum(DM))
    print(1 - Cf[1]/np.sum((DM-np.sum(DM)/len(DM))**2))
  # print(Cf[0])


#Plotting and printing
    plt.plot([2,3],[3,4])
    plt.show()
    plt.clf()
    logaritmic = False
    print("Cf DM_Regression DM_GFED landcover")
    for land in landcovers:
        mask2 = np.logical_and(maskHighest,compressedDLC[land])
        FRPplot = FRP[landcover_types.index(land)][mask2]
        DMplot = DM[mask2]
        print("%.3e"%Cf[0][landcover_types.index(land)], "%.3e"%np.sum(FRPplot*Cf[0][landcover_types.index(land)]), "%.3e"%np.sum(DMplot),land)
        if len(FRPplot) < 2:
            continue
        try:
            Cf2, b = linear_regression(FRPplot,DMplot,proportional=True)
        except:
            Cf2 = 1
        res = np.sum((Cf[0][landcover_types.index(land)]*FRPplot -DMplot)**2)
        restot = np.sum((DMplot-np.sum(DMplot)/len(DMplot))**2)
        R2 = 1- res/restot
        if logaritmic:
            mask = np.logical_not(np.logical_or(FRPplot < 1e-15, DMplot < 1e-15))
            FRPplot = np.log10(FRPplot[mask])
            DMplot = np.log10(DMplot[mask])
            plt.plot(FRPplot,FRPplot + np.log10(Cf[0][landcover_types.index(land)]),label="Total fit",color="black")

            plt.plot(FRPplot,FRPplot + np.log10(Cf2),":",color = "red",label="Plotted points fit")
        else:
            plt.plot(FRPplot,FRPplot*Cf2,":",color = "red",label="Plotted points fit")
            plt.plot(FRPplot,FRPplot*Cf[0][landcover_types.index(land)],label="Total fit",color="black")
        plt.hist2d(FRPplot,DMplot,bins=100,norm="log")
        if logaritmic:
            plt.ylim([np.log10(0.5e-10),-6])
            plt.xlim([-5,-1])
        plt.title(r"%s, $R^2 = $%f"%(land,R2))
        plt.colorbar()
        plt.xlabel(r"FRP [$W/m^2s$]")
        plt.ylabel(r"DM [$kg/m^2s$]")
        plt.legend()
        dir = str(deltat) + "plot"
        if not os.path.exists(dir):
            os.makedirs(dir)
        plt.savefig(os.path.join(dir,'%s_scattered.png'%land))
        plt.clf()
        
            

hour = dt.datetime(2003,1,2,hour = 12)
hourX = dt.datetime(2020,12,31,hour = 12)


# times = [(hour,dt.datetime(2003,12,31,hour=12))]
# for i in range(2004,2021):
#     times.append((dt.datetime(i,1,1,hour=12),dt.datetime(i,12,31,hour=12)))
timedelta = 10
savefile = "DM_Regression_Equation_System_dt=%d.pickle"%timedelta
oldfile = "DM_Regression_Equation_System.pickle"

# Run SetupRegression to produce the pickle file with the plot data
# Run Cf_Regression_Plot to plot the data from the pickle file. If no pickle file exists it will run SetupRegression with PARALLELL=False, which can be very slow.

#SetupRegression(hour,hourX,10,savefile=savefile,PARALLELL=True)
Cf_Regression_Plot(hour,hourX,10,savefile=savefile)




