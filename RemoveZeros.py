import numpy as np
import firegrib as fg
import os 
import sys
import matplotlib.pyplot as plt
import multiprocessing as mp
import datetime as dt
import numpy.ma as ma
from copy import deepcopy
import tracemalloc
import pickle
from statistics import linear_regression



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
