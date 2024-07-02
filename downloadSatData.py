import cdsapi
import datetime as dt
import os
import multiprocessing as mp
import sys

c = cdsapi.Client()

#Getting species from commandline in case I need to batch many at once.
s = sys.argv[1]

# species = ['wildfire_flux_of_acetone', 'wildfire_flux_of_ammonia', 'wildfire_flux_of_black_carbon',
#             'wildfire_flux_of_carbon_dioxide', 'wildfire_flux_of_carbon_monoxide', 'wildfire_flux_of_dimethyl_sulfide',
#             'wildfire_flux_of_ethane', 'wildfire_flux_of_formaldehyde', 'wildfire_flux_of_isoprene',
#             'wildfire_flux_of_methane', 'wildfire_flux_of_nitrogen_oxides', 'wildfire_flux_of_nitrous_oxide',
#             'wildfire_flux_of_non_methane_hydrocarbons', 'wildfire_flux_of_organic_carbon', 'wildfire_flux_of_propane',
#             'wildfire_flux_of_propene', 'wildfire_flux_of_sulphur_dioxide']


dicts = []
temp = {}

hour = dt.datetime(2003,1,2)
hourX = dt.datetime(2020,12,31)


date = hour.strftime("%Y-%m-%d/") + hourX.strftime("%Y-%m-%d")

c.retrieve(
    'cams-global-fire-emissions-gfas',
    {
        'date': date,
        'format': 'grib',
        'variable': s,
    },
    "copernicus/%s.grib"%(s))
