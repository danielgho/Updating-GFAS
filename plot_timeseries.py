import matplotlib.pyplot as plt
import numpy as np
import firegrib as fg
import datetime as dt
import os
import pickle
from statistics import linear_regression
import sys
import multiprocessing as mp

param = "DM"
Region = 'global'
#landcovers = [boreal forest]
landcovers = None


def landCoverMask(landcovers):
    if not isinstance(landcovers, list):
        landcovers = [landcovers]
    dlc = fg.Field('fire_biomes_with_peat_percentcover_2018_0p1deg_v1p0_compatibalised.grb2')
#    g = fg.Field('fire_biomes_with_peat_percentcover_2018_0p1deg_v1p0_compatibalised.grb2')
    land2emi = {}
    with open('CAMS2NEIVA.dat', 'r') as fp:
        for line in fp:
            temporaryList = []
            for value in line.split()[1::]:
                temporaryList.append(value)
            land2emi[line.split()[0]] = temporaryList[:]

    RelativePeat = False
    masks = {}
    if "peat" in land2emi.keys():
        if land2emi["peat"][0] == "R":
            del land2emi["peat"]
            RelativePeat = True
            peatmask = (dlc.val % 10 == 9)     
            for land in landcovers:
                if land == 'peat':
                    continue
                mask = np.full(dlc.val.shape, False)

                for id in land2emi[land]:
                    #Checking 3 leftmost digits in the land cover value. (Rightmost digit represents peat and is ignored)         
                    mask = np.logical_or(mask, (np.logical_and(dlc.val // 10 == int(id) // 10, np.logical_not(peatmask))))              
                masks[land] = mask
            masks["peat"] = peatmask

    if not RelativePeat:
        for land in landcovers:
            mask = np.full(dlc.val.shape, False)
            
            for id in land2emi[land]:         
                mask = np.logical_or(mask, dlc.val == int(id))
            masks[land] = mask


    mask = np.full(dlc.val.shape, False)
    for land in landcovers:
        mask = np.logical_or(mask, masks[land])
    return mask


## Old way of reading in copernicus data:

# def readTimeSeries(param,hour,hourX,lsmFile=None, spuriousFile=None, region='global', landcover=None):
#     print('reading time series from files...')
#     b = []
#     fields = []
    

#     while hour <= hourX:
#         directory = param
#         print(hour)
#         with open(os.path.join(directory,"%d.grib"%hour.day)) as fp:
#             try:
#                 f = fg.field(fp)
#             except:
#                 break
#             fields.append(f)
#             hour += dt.timedelta(days=1)
    
#     for f in fields:
#         if lsmFile:
#             lsm = fg.field(lsmFile)
#             f.val *= lsm.val
#         if spuriousFile:
#             spu = fg.field(spuriousFile)
#             f.val[spu.val>0] = 0
#         b.append(f.budget(region=region))
        
#     return b

#Read Copernicus daily time series
def readCopernicusTimeSeries(filename,hour,hourX,spuriousFile=None, region='global', landcover=None):
    path = '/xnilu_wrk/users/dgho/gfas/satDataAna/copernicus/archive'
    if spuriousFile != None:
        spu = fg.field(spuriousFile)
    if landcover != None:
        dlcmask = landCoverMask(landcover)
    b = []
    t = []
    hourtemp = hour
    while hourtemp <= hourX:
        t.append(hourtemp)
        hourtemp += dt.timedelta(days=1)

    while hour <= hourX:
        print(hour)
        f = fg.Field(os.path.join(path,hour.strftime('%Y%m%d'),filename))
        hour += dt.timedelta(days=1)
        if landcover:
            f.val = f.val*dlcmask
        if spuriousFile != None:
            spu = fg.field(spuriousFile)
            f.val[spu.val>0] = 0

        b.append(f.budget(region=region))


    return t, b



#Read GFAS daily time series
def readArchiveTimeSeries(filename,hour,hourX,spuriousFile=None, region='global', landcover=None):
    path = '/xnilu_wrk/users/dgho/gfas/archive/mirrorMARS/0001/'
    if spuriousFile != None:
        spu = fg.field(spuriousFile)
    if landcover != None:
        dlcmask = landCoverMask(landcover)
    b = []
    t = []
    hourtemp = hour
    while hourtemp <= hourX:
        t.append(hourtemp)
        hourtemp += dt.timedelta(days=1)

    while hour <= hourX:
        print(hour)
        f = fg.Field(os.path.join(path,hour.strftime('%Y%m%d'),hour.strftime('%H%M'),filename))
        hour += dt.timedelta(days=1)
        if landcover:
            f.val = f.val*dlcmask
        if spuriousFile != None:
            spu = fg.field(spuriousFile)
            f.val[spu.val>0] = 0

        b.append(f.budget(region=region))


    return t, b
#Read GFED daily time series from hour to hourX
def readGFEDTimeSeries(filename,hour,hourX, region='global', landcover=None):
    if landcover:
        dlcmask = landCoverMask(landcover)
    t = []
    b = []
    hourtemp = hour
    while hourtemp <= hourX:
        t.append(hourtemp)
        hourtemp += dt.timedelta(days=1)

    path = '/xnilu_wrk/users/dgho/gfas/satDataAna/GFED_DM/archive2/'
    while hour <= hourX:
        print(hour)
        f = fg.Field(os.path.join(path,hour.strftime('%Y%m%d'),filename))

        if landcover:
            f.val = f.val*dlcmask
        hour += dt.timedelta(days=1)
    
        b.append(f.budget(region=region))

    return t,b


#Read GFED field and sum fields from hour to hourX
def ReadTotalGFEDField(filename,hour,hourX, region='global', landcover=None):
    
    pickleName = filename + 'GFED' + 'field'+'.pickle'
    if os.path.isfile(pickleName):
        with open(pickleName,'rb') as fp:
            FinalField = pickle.load(fp)
    else:
        path = '/xnilu_wrk/users/dgho/gfas/satDataAna/GFED_DM/archive/'

        FinalField = np.zeros((1800,3600))
        while hour <= hourX:
            print(hour)
            f = fg.Field(os.path.join(path,hour.strftime('%Y%m%d'),filename))
            FinalField += f.val
            hour += dt.timedelta(days=1)
        with open(pickleName,'wb') as fp:
            pickle.dump(FinalField,fp)
    
    if landcover != None:
        FinalField = FinalField*landCoverMask(landcover)

    return FinalField

#Read GFAS field and sum fields from hour to hourX
def readArchiveField(filename,hour,hourX,spuriousFile=None, region='global', landcover=None):
    path = '/xnilu_wrk/users/dgho/gfas/archive/mirrorMARS/0001/'
    if os.path.isfile(filename+'field' + ".pickle"):
        with open(filename+'field' + ".pickle",'rb') as fp:
            FinalField = pickle.load(fp)
    else:
        FinalField = np.zeros((1800,3600))
        while hour <= hourX:
            print(hour)
            f = fg.Field(os.path.join(path,hour.strftime('%Y%m%d'),hour.strftime('%H%M'),filename))
            FinalField += f.val
            hour += dt.timedelta(days=1)
        with open(filename+'field' + ".pickle",'wb') as fp:
            pickle.dump(FinalField,fp)

    if landcover != None:
        FinalField = FinalField*landCoverMask(landcover)
    if spuriousFile != None:
        spu = fg.Field(spuriousFile)
        FinalField[spu.val>0] = 0
    

    return FinalField



#Old plotting code:
def plotTimeSeries(t,y, param=param, title='', lab='GFAS'):
    lss = ('solid', 'dashdot', 'dashed', 'dotted', '', '', '', '')
    mks = ('', '', '', '', '.', '*', 'x', '+', 'X' )
    if param == 'FRP':
        units = 'MW'
        scale = 1e9
    else:
        units = 'kg/s'
        scale = 1e3
    plt.ylabel(f'global {param} [{units}]')
    plt.plot(t, y,label=lab)

#Scatter daily budgets:
def scatterPlot(species,landcover,hour,hourX,spuriousFile='GFAS_MCD14ML_VNP14ML_2020.grb2',region='global'):

    if False:#os.path.exists("scatter.pickle"):
        with open("scatter.pickle",'rb') as fp:
            t,b1,b2 = pickle.load(fp)
    
    else:
        t,b1 = readArchiveTimeSeries(str(fg.PARAMID[species[0]]) + ".grb.gz",hour,hourX,spuriousFile=spuriousFile, region=region,landcover=landcover)
        t,b2 = readGFEDTimeSeries(species[1] + '.grb',hour,hourX, region=region, landcover=landcover)
#        with open("scatter.pickle",'wb') as fp:
#            pickle.dump([t,b1,b2],fp)
    

    
    b1 = np.array(b1)
    b2 = np.array(b2)
    
    if landcover == None:
        landcover = 'All'
    if species[0] == 'FRP':
        units = '[GW/s]'
        b1 /= 1e9
    else:
        units = '[Mg/s]'
        b1 /= 1000
    b2 /= 1000
    max_val1 = np.max(b1)
    max_val2 = np.max(b2)
    check = max_val1 > max_val2
    max_val = check*max_val1 + (not check)*max_val2
    max_val = 1.1*max_val
    plt.xlim(0,max_val)
    plt.ylim(0,max_val)
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    slope, intercept = linear_regression(b1, b2, proportional=True)
    plt.scatter(b1,b2)
    b1[0] = 0
    plt.plot(b1,slope*b1,color='black')
    plt.xlabel('GFAS %s %s'%(species[0],units))
    plt.ylabel('GFED %s [Mg/s]'%species[1])
    plt.title('GFED %s vs GFAS %s, %s'%(species[1],species[0], landcover))
    plt.savefig('GFED %s vs GFAS %s,%s.png'%(species[1],species[0],landcover))
    
    plt.clf()
    return slope
    

def monthly_series(b,hour,hourX):
    b3 = []
    t3 = []
    i = 0
    temp1 = 0
    while hour <= hourX:
        print(b[i])
        temp1 += b[i]
        if (hour + dt.timedelta(days=1)).month != hour.month:
            if hour.year == 2003 and hour.month == 1:
                timediff = hour.day - 1
            else:
                timediff = hour.day
            b3.append(temp1/timediff)
            t3.append(dt.datetime(hour.year,hour.month,15))
            temp1 = 0
            print(hour)
        print(hour)
        hour += dt.timedelta(days=1)
        i += 1
    return t3,b3

def plot_GFAS_GFED_COP(species,hour,hourX,pickleName="GFEDGFASCOP",reRead=False,monthly=True,spuriousFile='GFAS_MCD14ML_VNP14ML_2003_2021_merged.grb2'):
    name = "%s%s.pickle"%(pickleName,species)
    if os.path.exists(name) and (not reRead):
        with open(name,"rb") as fp:
            t,GFAS,GFED,COP = pickle.load(fp)
    else:
        t,GFAS = readArchiveTimeSeries("%d.grb.gz"%fg.PARAMID[species],hour,hourX,spuriousFile=spuriousFile)
        t,GFED = readGFEDTimeSeries("%s.grb"%species,hour,hourX)
        t,COP = readCopernicusTimeSeries("%s.grb"%species,hour,hourX,spuriousFile=spuriousFile)
        if monthly:
            t,GFAS = monthly_series(GFAS,hour,hourX)
            t,GFED = monthly_series(GFED,hour,hourX)
            t,COP = monthly_series(COP,hour,hourX)

        with open(name,"wb") as fp:
            pickle.dump([t,GFAS,GFED,COP],fp)
    
    maximum = max(np.max(GFAS),np.max(GFED),np.max(COP))

    plt.plot(t,GFAS,label="GFAS New EF",color='blue')
    plt.plot(t,GFED,label="GFED5",color='red')
    plt.plot(t,COP,label="GFAS Copernicus",color='yellow')


    plt.ylim(0,maximum*1.1)
    plt.xticks(rotation=45, ha='right')
    plt.legend()
    plt.title('Global budget of %s [kg/s]'%species)
    plt.ylabel("%s"%species)
    plt.savefig('TimeSeriesPlots/%s%s.png'%(pickleName,species),dpi=300, bbox_inches = "tight")
    plt.show()
    plt.clf()   

    

# Species = ['C3H6O','NH3','BC','CO2','CO','C2H6S','C2H6','CH2O','C5H8','CH4','NOx','N2O']#['DM','NMHC','OC','C3H8','C3H6','SO2']
# hour = dt.datetime(2003,1,2,hour=12,minute=0)
# hourX = dt.datetime(2020,12,31,hour=12,minute=0)
# pickleName = "FullRun"
# args = []
# for species in Species:
#     args.append([species,hour,hourX,pickleName])


# pool = mp.Pool(2)
# pool.starmap(plot_GFAS_GFED_COP,args)

hour = dt.datetime(2020,1,1,hour=12,minute=0)
hourX = dt.datetime(2020,3,1,hour=12,minute=0)
name = str(fg.PARAMID['NOx']) + ".grb.gz"
b3 = readTimeSeries('NOx',hour,hourX)
t,b = readArchiveTimeSeries(name,hour,hourX)
t,b2 = readCopernicusTimeSeries('NOx.grb',hour,hourX)


plt.plot(t,b,label="GFAS New EF",color='blue')
plt.plot(t,b2,label="new COP",color='red')
plt.plot(t,b3,label="GFAS Copernicus",color='yellow')
plt.legend()
plt.savefig("test.png")







# hour = dt.datetime(2003,1,2,hour=12,minute=0)
# hourX = dt.datetime(2020,12,31,hour=12,minute=0)

# t2,b2 = readArchiveTimeSeries('210100.grb.gz',hour,hourX,spuriousFile='GFAS_MCD14ML_VNP14ML_2003_2021_merged.grb2')
# t1,b1 = readGFEDTimeSeries('DM.grb',hour,hourX)
# #t1,b1 = readTimeSeries(hour,hourX,spuriousFile='GFAS_MCD14ML_VNP14ML_2003_2021_merged.grb2')


# # #monthly values:
# temp1 = 0
# temp2 = 0
# b3 = []
# b4 = []
# t3 = []
# hour = dt.datetime(2003,1,2,hour=0,minute=0)
# hourX = dt.datetime(2020,12,31,hour=0,minute=0)
# i = 0
# while hour <= hourX:
#     temp1 += b1[i]
#     temp2 += b2[i]
#     if (hour + dt.timedelta(days=1)).month != hour.month:
#         if hour.year == 2003 and hour.month == 1:
#             timediff = hour.day - 1
#         else:
#             timediff = hour.day
#         b3.append(temp1/timediff)
#         b4.append(temp2/timediff)
#         t3.append(dt.datetime(hour.year,hour.month,15))
#         temp1 = 0
#         temp2 = 0
#         print(hour)
#     print(hour)
#     hour += dt.timedelta(days=1)
#     i += 1

# b1 = b3
# b2 = b4
# t2 = t3

# b3 = np.array(b3)
# b4 = np.array(b4)
# print('Total Emission GFED:', np.sum(b3))
# print('Total Emission GFAS:', np.sum(b4))

# plotTimeSeries(t2,b3,lab='GFED. Avg: %.1f'%(np.sum(b3)/len(b3)))
# plotTimeSeries(t2,b4,lab='GFAS. Avg: %.1f'%(np.sum(b4)/len(b4)))

# plt.ylim(0,max(np.max(b3),np.max(b4))*1.1)
# plt.xticks(rotation=45, ha='right')
# plt.legend()
# plt.title('Global budget of DM')
# plt.savefig('2003-2020Comparison.png',dpi=300, bbox_inches = "tight")
# plt.show()
# plt.clf()   



# speciess = [['DM','DM'],['FRP','DM']]
# landcovers = ['savanna', 'boreal_forest', 'crop_residue', 'peat', 'temperate_forest', 'tropical_forest',None]
# args = []
# for i in range(len(speciess)):
#     for j in range(len(landcovers)):
#         args.append((speciess[i], landcovers[j],hour,hourX))
# with open("CF.dat",'w') as fp:
#     for species in speciess:
#         for landcover in landcovers:
#             m = scatterPlot(species,landcover,hour,hourX)
#             fp.write(str(m) + " ")
#             if landcover == None:
#                 landcover = "All"
#             fp.write(landcover)
#             fp.write("\n")





# pool = mp.Pool()

# a = pool.starmap(scatterPlot, args)

# with open("CF.dat",'w') as fp:
#     for i in range(len(args)):
#         if args[i][0][0] == 'FRP':
#             print(a[i],args[i],flush=True)
#             fp.write(str(a[i]) + " " + args[i][1])
















#retrieve(e, c)

#t1,b1 = readGFEDTimeSeries('DM.grb',hour,hourX)
# t1,b1 = readTimeSeries(hour,hourX,spuriousFile='GFAS_MCD14ML_VNP14ML_2003_2021_merged.grb2')
# t2,b2 = readArchiveTimeSeries('210100.grb.gz',hour,hourX,spuriousFile='GFAS_MCD14ML_VNP14ML_2003_2021_merged.grb2')

# b2 = np.array(b2)
# b1 = np.array(b1)
# print('Total Emission GFAS:', np.sum(b2))
# print('Total Emission Copernicus:', np.sum(b1))

# plotTimeSeries(t2,b1,lab='COP')
# plotTimeSeries(t2,b2)

# plt.legend()
# plt.title('Global budget of DM')
# plt.savefig('%s_Cop_comparison.png'%param)
# plt.show()
# plt.clf()    

# hour = dt.datetime(2020,1,1,hour=12)
# hourX = dt.datetime(2020,12,31,hour=12)
# f1 = ReadTotalGFEDField('DM.grb',hour,hourX,landcover=landcovers)
# f2 = readArchiveField('210100.grb.gz',hour,hourX,spuriousFile=None, region='global', landcover=landcovers)
# f1 *= 1000
# f2 *= 1000


# F = (f2-f1)

# PlottingField = fg.Field('GFAS_MCD14ML_VNP14ML_2003_2021_merged.grb2')
# PlottingField.val = F
# check = np.max(F) > np.abs(np.min(F))
# print(check*np.max(F))
# F[0,0] = -(check*np.max(F)) - (np.min(F)*(not check))
# PlottingField.plot(outFile=('global' +'differences'+ param +".png"),countries='black',cmap='seismic')

# PlottingField.val = f1
# PlottingField.plot(outFile=('global' +'GFED'+ param +".png"),countries='black')

# PlottingField.val = f2
# PlottingField.plot(outFile=('Same_scale_as_GFED_'+'global' +'GFAS'+ param +".png"),countries='black',maxValue=0.075,cmap='gnuplot2_r')






