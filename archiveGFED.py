import os
import datetime as dt
import firegrib as fg


speciess = ['DM']
for species in speciess:
    for i in range(1,13):
        filename = species + "%02d"%i + '.grb'
        path = os.path.join("GFED_DM/2015",filename)
        with open(path,'rb') as fp:
            while True:
                try:
                    f = fg.Field(fp)
                except:
                    print(path)
                    break
                f.val = f.val/(1000*24*3600) #converting from g/day to kg/s
                param = fg.PARAMID[species]
                time = f.getKey('date')
                outdir = os.path.join('GFED_DM','archive',str(time))
                if not os.path.exists(outdir):
                    os.makedirs(outdir)
                f.streamout(os.path.join(outdir,species + ".grb"))
                print(time)

