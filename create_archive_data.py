import firegrib as fg
import numpy as np
import datetime as dt
import os
import utilities as ut


datadir = 'Area2015'
archivePath = '/xnilu_wrk/users/dgho/gfas/archive/mirrorMARS/0001'
f = fg.Field(os.path.join(datadir,'1.grib'))

YYYYMMDD = f.getKey('date')
hour = dt.datetime(YYYYMMDD // 10000, (YYYYMMDD // 100) % 100, YYYYMMDD % 100)

for i in range(1,366):
    print(i)
    f = fg.Field(os.path.join(datadir,'%d.grib'%i))
    f.setKey('time', 1200)
    ut.archive(f,archivePath)


    
