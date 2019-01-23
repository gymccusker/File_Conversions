##--------------------------------------------------------------------------
##
##			Script to read in WRF output files, extract necessary data,
##			then save into new NetCDF file (reduces file size for archiving)
##					-- GYoung
##
##--------------------------------------------------------------------------

from netCDF4 import Dataset
import numpy as np
import time 
from datetime import datetime, timedelta 
from netCDF4 import num2date, date2num 
import matplotlib.pyplot as plt

##--------------------------------------------------------------------------
##--------------------------------------------------------------------------
##---------------				IN
##--------------------------------------------------------------------------
##--------------------------------------------------------------------------

###################################
###################################
## LEM OUTPUT
###################################
###################################

# run1 = 120
# info1 = 'C86_Ocean'

# run1 = 135
# info1 = 'D10_Ocean'

run1 = 143
info1 = 'ACC_Ocean'

# run7='120';		YES
# info7='C86';

# run8='135';		YES
# info8='D10';

# run9='143';		YES
# info9='ACC';

## /gws/nopw/j04/ncas_weather/gyoung/ACCACIA/ModellingScripts/mod_obs_comp_noNice.m


###################################
# Define time dump separation
###################################

hours = np.array([7,11,15,19,23,27,31,35,39,43,47,51,55,59,63,67,71,75,79,83,87,91,95,99])

# os.chdir("../LEM/r1")

L_vap = 2.5e6    # J/kg
L_sub = 2.836e6  # J/kg
cp = 1004.6      # J/kg.K

###################################
# Load in data
###################################

nc1 = {}			# define nc1 as a dictionary
strg1 = "%2.f" % run1 

filedir = '/gws/nopw/j04/ncas_weather/gyoung/ACCACIA/LEM/r'
rundir = "".join([filedir,strg1,'/'])

for i in range(0, len(hours)):
	strg2 = "%02d" % hours[i] # string of hour index
	a1 = ''.join([rundir,'RUN0',strg1,'_00',strg2,'.nc']) # string of filename
	strgi = "%1.f" % (i+1) # string of hour number
	nc1[strgi] = Dataset(a1,'r')

###################################
## Read in NetCDF variables to usable variables
###################################

## Define data
icenum1 = np.zeros([24,104])
largeice1 = np.zeros([24,104])
liqmass1 = np.zeros([24,104])
watvap1 = np.zeros([24,104])
tempK1 = np.zeros([24,104])
pres1 = np.zeros([24,104])
evs1 = np.zeros([24,104])
qvs1 = np.zeros([24,104])
rh1 = np.zeros([24,104])
incloud1 = np.zeros([24,104])
rh1 = np.zeros([24,104,130,130])
for i in range(0, len(hours)):
	strgi = "%1.f" % (i+1) # string of hour number
	icenum1[i,:] = (nc1[strgi]['QBAR07'][1:]+nc1[strgi]['QBAR08'][1:]+nc1[strgi]['QBAR09'][1:])/1000
	largeice1[i,:] = (nc1[strgi]['ALL_Ni100'][1:]+nc1[strgi]['ALL_Nis100'][1:])/1000
	liqmass1[i,:] = nc1[strgi]['QBAR02'][1:]*1000
	watvap1[i,:] = nc1[strgi]['QBAR01'][1:]*1000
	tempK1[i,:] = nc1[strgi]['ALL_TEMP'][1:]
	pres1[i,:] = nc1[strgi]['PREFN'][1:]
	evs1[i,:] = (0.611*np.exp(17.27*(tempK1[i,:]-273.15)/((tempK1[i,:]-273.15)+237.3)))*1000
	qvs1[i,:] = (0.622*evs1[i,:])/(pres1[i,:]-evs1[i,:])
	rh1[i,:] = ((watvap1[i,:]/1000)/qvs1[i,:])*100
	# incloud1[i,:] = (rh1[i,:]>=100).nonzero()
	ice1[i,:,:,:] = (nc1[strgi]['Q07'][1:]+nc1[strgi]['Q08'][1:]+nc1[strgi]['Q09'][1:])/1000
timesec1 = (nc1['24']['TIMES'][:])/3600
Z1 = nc1['24']['ZN'][1:]
# X1 = nc1['24']['XN'][:]
# Y1 = nc1['24']['YN'][:]
times1=np.arange(1,25)

##--------------------------------------------------------------------------
##--------------------------------------------------------------------------
##---------------				OUT
##--------------------------------------------------------------------------
##--------------------------------------------------------------------------

###################################
## Open File
###################################
outfile = "".join([info1,'.nc'])
dataset =  Dataset(outfile, 'w', format ='NETCDF4_CLASSIC') 

print dataset.file_format 

###################################
## Global Attributes
###################################
# info2 = ' Cooper et al., 1986 parametrization used for primary ice nucleation. Model assumes oceanic surface.'
# info2 = ' Approximation of the DeMott et al., 2010 parametrization used for primary ice nucleation (see Young et al., 2017 (ACP) for details). Model assumes oceanic surface.'
info2 = ' Parametrization derived from ACCACIA ice number concentration measurements (measured using a 2-Dimensional Stereo particle imaging probe) used for primary ice nucleation (see Young et al., 2017 (ACP) for details). Model assumes oceanic surface.'
desc = info1 + ' simulation from Young et al., 2017 (ACP). x/y grid size = 130x130 grid points (120m grid size) with 104 vertical levels (20m resolution up to 1500m, then 50m resolution between 1500m and 3000m). Domain size = 16km x 16km, centred on 75.0N, 24.5E. Model initialised with radiosonde data (sonde number 5) from ACCACIA flight B762 (23-MAR-2013).' + info2
dataset.description = desc
dataset.history = 'Created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S')
dataset.source = 'UK Met Office Large Eddy Model (LEM), version 2.4, coupled with the Morrison et al., 2005 (JAS) microphysics scheme (ported from the Weather Research and Forecasting model).' 
dataset.references = 'First published in Young et al., 2017 (ACP): Microphysical sensitivity of coupled springtime Arctic stratocumulus to modelled primary ice over the ice pack, marginal ice, and ocean. (doi:10.5194/acp-17-4209-2017)'
dataset.project = 'Aerosol-Cloud Coupling and Climate Interactions in the Arctic (ACCACIA), funded by the UK Natural Environment Research Council (Grant no. NE/I028696/1).'
dataset.comment = 'Other LEM variables from this simulation are archived locally at the University of Manchester. Contact Gillian Young (G.Young1@leeds.ac.uk) for details.'
dataset.institution = 'University of Manchester.'

###################################
## Switch off automatic filling 
###################################
dataset.set_fill_off()

###################################
## Data dimensions
###################################
time = dataset.createDimension('time', np.size(icenum1,0))
level = dataset.createDimension('level', np.size(icenum1,1)) 

###################################
## Dimensions variables
###################################
#### Times
times = dataset.createVariable('time', np.float32, ('time',),fill_value='-9999') 
times.comment = 'hourly data dumps'
times.units = ['hours since 09:00:00 on 23-MAR-2013']
times[:] = times1[:]

#### Levels
Z = dataset.createVariable('Z', np.float32, ('level',),fill_value='-9999') 
Z.long_name = 'Altitude'
Z.comment = 'Levels spaced by 20m up to 1500m, then spaced by 50m between 1500m and 3000m'
Z.units = 'm'
Z[:] = Z1[:]

###################################
###################################
## Create 2-d variables
###################################
###################################
#### Nisg
nisg = dataset.createVariable('nisg', np.float32, ('time','level',),fill_value='-9999')
nisg.long_name = 'total ice number concentration'
nisg.comment = 'Sum of ice, snow, and graupel particles'
nisg.units = 'L-1'
nisg[:] = icenum1[:]

#### Nisg100
nisg100 = dataset.createVariable('nisg100', np.float32, ('time','level',),fill_value='-9999')
nisg100.long_name = 'total ice number concentration greater than 100micron'
nisg100.comment = 'Sum of ice, snow, and graupel particles of sizes greater than 100micron'
nisg100.units = 'L-1'
nisg100[:] = largeice1[:]

#### Qliq
qliq = dataset.createVariable('qliq', np.float32, ('time','level',),fill_value='-9999')
qliq.long_name = 'cloud liquid mass mixing ratio'
qliq.comment = 'Only cloud liquid field included; rain category is separate.'
qliq.units = 'g kg-1'
qliq[:] = liqmass1[:]

#### Qvap
qvap = dataset.createVariable('qvap', np.float32, ('time','level',),fill_value='-9999')
qvap.long_name = 'water vapour mixing ratio'
qvap.comment = ''
qvap.units = 'g kg-1'
qvap[:] = watvap1[:]

#### Temperature
temp = dataset.createVariable('temp', np.float32, ('time','level',),fill_value='-9999')
temp.long_name = 'temperature'
temp.comment = ''
temp.units = 'K'
temp[:] = tempK1[:]


###################################
## Write out file
###################################
# dataset.close()



###################################
## TESTING
###################################
# nc1 = Dataset('C86_Ocean.nc','r')
# nc2 = Dataset('D10_Ocean.nc','r')
# nc3 = Dataset('ACC_Ocean.nc','r')

# fig = plt.figure(figsize=(4,6))

# plt.subplot(211)
# plt.plot(nc1.variables['nisg100'][6,:],nc1.variables['Z'][:],'m')
# plt.plot(nc2.variables['nisg100'][6,:],nc2.variables['Z'][:],'g')
# plt.plot(nc3.variables['nisg100'][6,:],nc3.variables['Z'][:],'b')
# plt.xlim([0,1.5])
# plt.ylim([0,2000])
# plt.grid('on')

# plt.subplot(212)
# plt.plot(nc1.variables['qliq'][6,:],nc1.variables['Z'][:],'m')
# plt.plot(nc2.variables['qliq'][6,:],nc2.variables['Z'][:],'g')
# plt.plot(nc3.variables['qliq'][6,:],nc3.variables['Z'][:],'b')
# plt.xlim([0,0.7])
# plt.ylim([0,2000])
# plt.grid('on')
# plt.show()