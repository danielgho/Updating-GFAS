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

#Creating DLC masks
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

def fix_latlon(field1,field2):
    field1.lat = field2.lat
    field1.lon = field2.lon
    field1.area = field2.area


def daily_contribution(hour):
    print(hour,flush=True)
    #Define Variables
    DM = fetch_GFED_DM(hour)
    FRP = fetch_GFAS_frp(hour)
    contributing_DM = {}
    contributing_FRP = {}
    Cells = {}

    #Initialize Field objects
    coarse_FRP = fg.Field('fire_biomes_with_peat_percentcover_2018_0p1deg_v1p0_compatibalised.grb2')
    coarse_DLC_FRP = fg.Field('fire_biomes_with_peat_percentcover_2018_0p1deg_v1p0_compatibalised.grb2')
    coarse_DLC_mask = fg.Field('fire_biomes_with_peat_percentcover_2018_0p1deg_v1p0_compatibalised.grb2')

    DM.coarser(2)

    coarse_FRP.val = FRP.val
    coarse_FRP.coarser(5)
    
    for land in landcovers:
        fix_latlon(coarse_DLC_FRP,FRP)
        fix_latlon(coarse_DLC_mask,FRP)

        coarse_DLC_mask.val = dlcmasks[land]
        coarse_DLC_mask.coarser(5,algo='sum')
        coarse_DLC_mask.val = coarse_DLC_mask.val > 0

        coarse_DLC_FRP.val = dlcmasks[land]*FRP.val
        coarse_DLC_FRP.coarser(5)

        contributing_cells = coarse_DLC_FRP.val > 0.9*coarse_FRP.val
        contributing_cells = np.logical_and(contributing_cells, DM.val > 0)
        contributing_DM[land] = contributing_cells*DM.val
        contributing_FRP[land] = contributing_cells*coarse_FRP.val #Could change to DLC_FRP instead of FRP
        Cells[land] = contributing_cells
    
    return contributing_FRP, contributing_DM, Cells

def process_timestep(hours):
    total_DM = {}
    total_FRP = {}
    total_contribution = {}

    for land in landcovers:
        total_DM[land] = np.zeros((360,720))
        total_FRP[land] = np.zeros((360,720))
        total_contribution[land] = np.zeros((360,720))


    for hour in hours:
        results = daily_contribution(hour)
        for land in landcovers:
            total_FRP[land] += results[0][land]
            total_DM[land] += results[1][land]
            total_contribution[land] += results[2][land]
    
    return total_FRP, total_DM, total_contribution

def process_total_field(hours_list, Parallell = False,savefile="Temporal_average_scattering.pickle"):
    total_DM = {}
    total_FRP = {}
    total_contribution = {}

    for land in landcovers:
        total_DM[land] = np.zeros((360,720))
        total_FRP[land] = np.zeros((360,720))
        total_contribution[land] = np.zeros((360,720))
    
    if Parallell:
        pool = mp.Pool()
        results = pool.map(process_timestep,list(hours_list))
    else:
        results = []
        for hours in hours_list:
            results.append(process_timestep(hours))
    
    for i in range(len(results)):
        for land in landcovers:
            total_FRP[land] += results[i][0][land]
            total_DM[land] += results[i][1][land]
            total_contribution[land] += results[i][2][land]
            

    


#Averaging FRP and DM over number of days included, and removing highest 5% DM and FRP values.
    for land in landcovers:
        total_FRP[land] = np.divide(total_FRP[land], total_contribution[land], out=np.zeros_like(total_FRP[land]), where=total_contribution[land]!=0)
        total_DM[land] = np.divide(total_DM[land], total_contribution[land], out=np.zeros_like(total_DM[land]), where=total_contribution[land]!=0)
        

        FRP_sorted = np.sort(total_FRP[land].flatten())
        DM_sorted = np.sort(total_DM[land].flatten())

        firstNonZero = (FRP_sorted != 0).argmax()
        max_frp = FRP_sorted[firstNonZero + int((len(FRP_sorted)-firstNonZero)*0.95)]

        firstNonZero = (DM_sorted != 0).argmax()
        max_dm = DM_sorted[firstNonZero + int((len(FRP_sorted)-firstNonZero)*0.95)]

        lowMask = np.logical_and(total_FRP[land] <= max_frp, total_DM[land] <= max_dm)
        total_FRP[land] = total_FRP[land][lowMask]
        total_DM[land] = total_DM[land][lowMask]

    
    with open(savefile,"wb") as fp:
        pickle.dump([total_FRP,total_DM], fp)
    
    return total_FRP, total_DM

def regression_scatter(hours_list,savefile="Temporal_average_scattering.pickle"):
    if os.path.exists(savefile):
        with open(savefile,"rb") as fp:
            results = pickle.load(fp)
            FRP = results[0]
            DM = results[1]

    else:
        FRP,DM = process_total_field(FRP,DM)
    
    for land in landcovers:
        mask = np.logical_not(np.logical_or(np.logical_and(FRP[land] == 0, DM[land] == 0),np.logical_or(np.isnan(FRP[land]), np.isnan(DM[land]))))
        FRPl = FRP[land][mask]
        DMl = DM[land][mask]
        try:
            cf,m = linear_regression(FRPl,DMl,proportional=True)
        except:
            print("hello")
            continue
        cf2 = np.sum(DMl)/np.sum(FRPl)
        print("%.4e"%cf,"%.4e"%cf2,land)

        logarithmic = False
        if logarithmic:
            mask = np.logical_not(np.logical_or(FRPl == 0, DMl == 0))
            FRPl = np.log(FRPl[mask])
            DMl = np.log(DMl[mask])
            plt.plot(FRPl,FRPl + np.log(cf),":", label="least square",color="black")
            plt.plot(FRPl,FRPl + np.log(cf2),label="average",color="black")
        else:
            plt.plot(FRPl,cf*FRPl,":", label="least square",color="black")
            plt.plot(FRPl,cf2*FRPl,label="average",color="red")
        #plt.scatter(FRPl,DMl)
        plt.hist2d(FRPl,DMl,bins=50,norm="log")
        if logarithmic:
            plt.ylim([-24,-13])
            plt.xlim([-8,-2])
        plt.colorbar()
        plt.legend()
        plt.title("%s"%land)
        plt.savefig("TemporalScatterPlotsNoZero/%s.png"%land)
        plt.clf()



hour = dt.datetime(2003,1,2,hour=12)
hourX = dt.datetime(2020,12,31, hour=12)
t = [hour + dt.timedelta(days=i) for i in range((hourX-hour).days + 1)]

#file = "test.pickle"
file = "NoZeroes.pickle"

hour_list = np.array_split(t,29)
hour_list = list(hour_list)

#process_total_field(hour_list,savefile=file,Parallell=True)
regression_scatter(hour_list,file)



