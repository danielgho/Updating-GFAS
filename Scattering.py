import matplotlib.pyplot as plt
import numpy as np
import os
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



def scatter(land):
    path = "SavedPairs/%s"%(land)
    FRP = np.zeros(1)
    DM = np.zeros(1)
    for years in range(2003,2021):
     #   print(years)
        with open(os.path.join(path,"%d.npy"%years),"rb") as fp:
            new = np.load(fp)
        FRP = np.append(FRP,new[0])
        DM = np.append(DM,new[1])

    DM = DM[1::]
    FRP = FRP[1::]

    slope, intercept = linear_regression(FRP, DM, proportional=True)

    plt.plot(FRP,slope*FRP,color='black')
 #   print(slope, land, len(FRP))
 #   b = np.sum(FRP*DM)/np.sum(FRP*FRP)
    print(np.sum(FRP*slope),np.sum(DM), len(DM),land, slope)
    plt.hist2d(FRP,DM,bins=100,norm="log")
    plt.title("%s"%land)
    plt.colorbar()
    plt.xlabel(r"FRP [$W/m^2s$]")
    plt.ylabel(r"DM [$kg/m^2s$]")
    plt.savefig('scatterAll/%s_scattered.png'%land)
    plt.clf()
    return land, slope


with open("CFout.dat","w") as fp:
    for land in landcovers:
        l,s = scatter(land)
        fp.write("%s %f\n"%(l,s))
        

