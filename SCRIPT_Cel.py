import matplotlib
matplotlib.use('Agg')
import pylab as plt
import numpy as np
#from matplotlib import pyplot as plt
import yt
import trident
from trident import LightRay
import illustris_python as il
from astropy import units as u
from astropy import constants as c
from yt.utilities.cosmology import Cosmology
#from matplotlib import pyplot as plt
import sys
import h5py
from constants import *

def getField(snapFile,partType,fieldName,toDouble=False):
        if toDouble==False:
                return snapFile[partType][fieldName][:]
        if toDouble==True:
                return snapFile[partType][fieldName][:].astype('f8')


def gadgetDens2SI(dens):
        return dens*1E10*Msol_si/hubble/(kpc_si/hubble)**3

snapshots = [99,91,84,78,72,67,59,50,40,33,25,21,17]#,13,11,8,6,4,3,2]
chunks = np.arange(0,9,1)

#simString = "L205n2500TNG" #large volume
#nFiles = 600
simString = 'L205n625TNG'
snapshotPath="/n/holylfs/LABS/hernquist_lab/IllustrisTNG/Runs/"+simString+"/output"
######################################
#initialise arrays that hold all data#
######################################

#redshift array
allsnaps_redshifts=[]

#baryon arrays
#PartType0 (gas)
allsnaps_tot_PT0_arr = []
allsnaps_hal_PT0_arr = []
allsnaps_fil_PT0_arr = []
allsnaps_voi_PT0_arr = []
#PartType4 (stars)
allsnaps_tot_PT4_arr = []
allsnaps_hal_PT4_arr = []
allsnaps_fil_PT4_arr = []
allsnaps_voi_PT4_arr = []

#check arrays
allsnaps_dens_crit_array = []
allsnaps_m_p_array = []

################################
#loop over snapshots and chunks#
################################

for i in range(len(snapshots)):
    
    ######################################################
    #initialise arrays that hold individual snapshot data#
    ######################################################

    #redshift array
    redshifts=[]

    #baryon arrays
    #PartType0 (gas)
    tot_PT0_arr = []
    hal_PT0_arr = []
    fil_PT0_arr = []
    voi_PT0_arr = []
    #PartType4 (stars)
    tot_PT4_arr = []
    hal_PT4_arr = []
    fil_PT4_arr = []
    voi_PT4_arr = []

    #check arrays
    dens_crit_array = []
    m_p_array = []
    
    for j in range(len(chunks)):
    
        ###########################
        #select snapshot and chunk#
        ###########################

        snapshot = snapshots[i]
        chunk = chunks[j]
        print('processsing snapshot {0:02d}, chunk {1}'.format(snapshot,chunk))

        ######################
        #load snapshot header#
        ######################

        print('loading snapshot header')
        basepath = '/n/holylfs/LABS/hernquist_lab/IllustrisTNG/Runs/'+simString+'/output/'#'/virgo/simulations/IllustrisTNG/TNG100-3/output/'
        header=il.groupcat.loadHeader(basepath,snapshot)
        print('header loaded.')
        

        #simString = "L205n2500TNG" #large volume
        #nFiles = 600
        #snapshotPath="/n/holylfs/LABS/hernquist_lab/IllustrisTNG/Runs/"+simString+"/output"

        fname = snapshotPath+'/snapdir_'+str(snapshot).zfill(3)+'/'+'snap_'+str(snapshot).zfill(3)+'.'+str(j)+'.hdf5'
        print fname

        snapFile = h5py.File(fname,'r')       

        ########################################################
        #extract redshift from header, append to redshift array#
        ########################################################

        redshifts.append(header['Redshift'])

        #####################
        #load Illustris data#
        #####################

        #data='/n/holylfs/LABS/hernquist_lab/IllustrisTNG/Runs/L205n625TNG/output/snapdir_0{0:02d}/snap_0{0:02d}.{1}.hdf5'.format(snapshot,chunk)
        #print('loading dataset: {0}'.format(data))
        #ds=yt.load(data)
        #print('data loaded.')

        ####################################
        #define proton mass in solar masses#
        ####################################

        #print('calculating proton mass:')
        #m_p = ds.quan((c.m_p).to(u.solMass).value,'Msun')
        #print ('proton mass: {0}'.format(m_p))
        #m_p_array.append(m_p)

        #ds.index

        ############################################
        #calculate critical density of the Universe#
        ############################################

        print('Calculating critical density of the Universe')
        #co=Cosmology()
        #grav=ds.quan(6.6743e-11,'m**3/(kg*s**2)')
        #H=co.hubble_parameter(0).in_units('km/s/Mpc')
        #print('grav={0}'.format(grav))
        #print('H0 = {0}'.format(H))
        #dens_crit = ((3 * H**2)/(8*np.pi* grav))#.in_units('h**2/(code_length**3/(Msun*s**2))')
        #print ('critical density is: {0}'.format(dens_crit))
        #print('critical density in kg/m^3: {0}'.format(dens_crit.in_units('kg/m**3')))
        #dens_crit_array.append(dens_crit)
        dens_crit = 8.61975454419e-27 #Celeste:extracted from "constant.py"
        hubble = 0.677
        #ad = ds.all_data()

        ##########################################################
        #convert dark matter density to units of critical density#
        ##########################################################

        #PartType0 (gas)
        #dark_units_PT0 = ad['PartType0','SubfindDMDensity']*ds.quan(1e10,"(Msun/h)/((code_length)**3)")/dens_crit
        gasMass = getField(snapFile,'PartType0','Masses',toDouble=True)
        #dark_units_PT0 = getField(snapFile,'PartType0','SubfindDMDensity',toDouble=True)*1E10/hubble/dens_crit
        DMDens0 = getField(snapFile,'PartType0','SubfindDMDensity',toDouble=True)
        DMDens0 = gadgetDens2SI(DMDens0)/dens_crit
        #PartType4 (stars)
        #dark_units_PT4 = ad['PartType4','SubfindDMDensity']*ds.quan(1e10,"(Msun/h)/((code_length)**3)")/dens_crit
        stellarMass = getField(snapFile,'PartType4','Masses',toDouble=True)
        #dark_units_PT4 = getField(snapFile,'PartType4','SubfindDMDensity',toDouble=True)*1E10/hubble/dens_crit
        DMDens4 = getField(snapFile,'PartType4','SubfindDMDensity',toDouble=True)
        DMDens4 = gadgetDens2SI(DMDens4)/dens_crit 
        ###############################################################
        #create Large-Scale Structure (LSS) masks for each matter type#
        ###############################################################

        mask0 = np.logical_and(DMDens0 >= 1E-11, DMDens0 < 1)
        mask4 = np.logical_and(DMDens4 >= 1E-11, DMDens4 < 1)
        print "Check something outside", len(DMDens4[mask4]), len(DMDens0[mask0])

        #PartType0 (gas)
        voi_mask_PT0 = DMDens0 < 0.1
        fil_mask_PT0 = np.logical_and(DMDens0 >= 0.1, DMDens0 < 57)#CELESTE:CORRECTED
        hal_mask_PT0 = DMDens0 >= 57 
        #PartType4 (stars)
        voi_mask_PT4 = DMDens4 < 0.1
        fil_mask_PT4 = np.logical_and(DMDens4 >= 0.1, DMDens4 < 57)#CELESTE:CORRECTED
        hal_mask_PT4 = DMDens4 >= 57

        ##########################################################
        #calculate mass in each LSS for each baryonic matter type#
        ##########################################################

        #PartType0 (gas)
        tot_PT0 = np.sum(gasMass) #np.sum(ad['PartType0','Masses'])
        hal_PT0 = np.sum(gasMass[hal_mask_PT0])#np.sum(ad['PartType0','Masses'][hal_mask_PT0])
        fil_PT0 = np.sum(gasMass[fil_mask_PT0])#np.sum(ad['PartType0','Masses'][fil_mask_PT0])
        voi_PT0 = np.sum(gasMass[hal_mask_PT0])#np.sum(ad['PartType0','Masses'][voi_mask_PT0])
        #PartType4 (stars)
        tot_PT4 = np.sum(stellarMass)#np.sum(ad['PartType4','Masses'])
        hal_PT4 = np.sum(stellarMass[hal_mask_PT4])#np.sum(ad['PartType4','Masses'][hal_mask_PT4])
        fil_PT4 = np.sum(stellarMass[fil_mask_PT4])#np.sum(ad['PartType4','Masses'][fil_mask_PT4])
        voi_PT4 = np.sum(stellarMass[voi_mask_PT4])#np.sum(ad['PartType4','Masses'][voi_mask_PT4])

        ##########################
        #append results to arrays#
        ##########################

        #PartType0 (gas)
        tot_PT0_arr.append(tot_PT0)
        hal_PT0_arr.append(hal_PT0)
        fil_PT0_arr.append(fil_PT0) 
        voi_PT0_arr.append(voi_PT0) 
        #PartType4 (gas)
        tot_PT4_arr.append(tot_PT4)
        hal_PT4_arr.append(hal_PT4)
        fil_PT4_arr.append(fil_PT4) 
        voi_PT4_arr.append(voi_PT4)
        
    ######################################################    
    #append all results for a snapshot to the final array#
    ######################################################
    
    #PartType0 (gas)
    allsnaps_tot_PT0_arr.append(tot_PT0_arr)
    allsnaps_hal_PT0_arr.append(hal_PT0_arr)
    allsnaps_fil_PT0_arr.append(fil_PT0_arr) 
    allsnaps_voi_PT0_arr.append(voi_PT0_arr) 
    #PartType4 (gas)
    allsnaps_tot_PT4_arr.append(tot_PT4_arr)
    allsnaps_hal_PT4_arr.append(hal_PT4_arr)
    allsnaps_fil_PT4_arr.append(fil_PT4_arr) 
    allsnaps_voi_PT4_arr.append(voi_PT4_arr)    
    
    allsnaps_m_p_array.append(m_p_array)
    allsnaps_dens_crit_array.append(dens_crit_array)
    allsnaps_redshifts.append(redshifts)

    
    print('filled baryon arrays')


#########################
#convert to numpy arrays#
#########################

#tot_PT0_arr=np.array(tot_PT0_arr)
#hal_PT0_arr=np.array(hal_PT0_arr)
#fil_PT0_arr=np.array(fil_PT0_arr)
#voi_PT0_arr=np.array(voi_PT0_arr)

#tot_PT4_arr=np.array(tot_PT4_arr)
#hal_PT4_arr=np.array(hal_PT4_arr)
#fil_PT4_arr=np.array(fil_PT4_arr)
#voi_PT4_arr=np.array(voi_PT4_arr)

allsnaps_tot_PT0_arr = np.array(allsnaps_tot_PT0_arr)
allsnaps_hal_PT0_arr = np.array(allsnaps_hal_PT0_arr)
allsnaps_fil_PT0_arr = np.array(allsnaps_fil_PT0_arr)
allsnaps_voi_PT0_arr = np.array(allsnaps_voi_PT0_arr)

allsnaps_tot_PT4_arr = np.array(allsnaps_tot_PT4_arr)
allsnaps_hal_PT4_arr = np.array(allsnaps_hal_PT4_arr)
allsnaps_fil_PT4_arr = np.array(allsnaps_fil_PT4_arr)
allsnaps_voi_PT4_arr = np.array(allsnaps_voi_PT4_arr)

allsnaps_redshifts = np.array(allsnaps_redshifts)
allsnaps_m_p_array = np.array(allsnaps_m_p_array)
allsnaps_dens_crit_array = np.array(allsnaps_dens_crit_array)

print("critical_density", allsnaps_dens_crit_array)
###############
#check lengths#
###############

print(len(snapshots))
print(allsnaps_redshifts.shape)
print(allsnaps_m_p_array.shape)
print(allsnaps_dens_crit_array.shape)
print(allsnaps_hal_PT0_arr.shape)
print(allsnaps_fil_PT0_arr.shape)
print(allsnaps_voi_PT0_arr.shape)
print(allsnaps_hal_PT4_arr.shape)
print(allsnaps_fil_PT4_arr.shape)
print(allsnaps_voi_PT4_arr.shape)


print(allsnaps_tot_PT0_arr.sum(axis=1).shape) #summing chunks for each snapshot
print(allsnaps_redshifts.mean(axis=1)) #recovering the redshift for each snapshot


##################################################
#for each of the 13 snapshots, sum the masses in #
#each of the 7 simulation subchunks to obtain    #
#the total mass for the entire redshift          #
##################################################

tot_PT0_arr_summed = allsnaps_tot_PT0_arr.sum(axis=1)
hal_PT0_arr_summed = allsnaps_hal_PT0_arr.sum(axis=1)
fil_PT0_arr_summed = allsnaps_fil_PT0_arr.sum(axis=1)
voi_PT0_arr_summed = allsnaps_voi_PT0_arr.sum(axis=1)

tot_PT4_arr_summed = allsnaps_tot_PT4_arr.sum(axis=1)
hal_PT4_arr_summed = allsnaps_hal_PT4_arr.sum(axis=1)
fil_PT4_arr_summed = allsnaps_fil_PT4_arr.sum(axis=1)
voi_PT4_arr_summed = allsnaps_voi_PT4_arr.sum(axis=1)

print(allsnaps_hal_PT0_arr.sum(axis=1).shape)


####################################################
#for each of the 13 snapshots, sum the masses in   #
#PT0 and PT4 to obtain the total baryonic mass for #
#each redshift                                     #
####################################################

tot_PT14_arr_summed = tot_PT0_arr_summed + tot_PT4_arr_summed
hal_PT14_arr_summed = hal_PT0_arr_summed + hal_PT4_arr_summed
fil_PT14_arr_summed = fil_PT0_arr_summed + fil_PT4_arr_summed
voi_PT14_arr_summed = voi_PT0_arr_summed + voi_PT4_arr_summed

print(hal_PT14_arr_summed.shape)


#######################
#create mass fractions#
#######################

#hal_frac = (hal_PT0_arr + hal_PT4_arr)/(hal_PT0_arr + hal_PT4_arr + fil_PT0_arr + fil_PT4_arr + voi_PT0_arr + voi_PT4_arr)
#fil_frac = (fil_PT0_arr + fil_PT4_arr)/(hal_PT0_arr + hal_PT4_arr + fil_PT0_arr + fil_PT4_arr + voi_PT0_arr + voi_PT4_arr)
#voi_frac = (voi_PT0_arr + voi_PT4_arr)/(hal_PT0_arr + hal_PT4_arr + fil_PT0_arr + fil_PT4_arr + voi_PT0_arr + voi_PT4_arr)


hal_frac = (hal_PT14_arr_summed)/(hal_PT14_arr_summed + fil_PT14_arr_summed + voi_PT14_arr_summed)
fil_frac = (fil_PT14_arr_summed)/(hal_PT14_arr_summed + fil_PT14_arr_summed + voi_PT14_arr_summed)
voi_frac = (voi_PT14_arr_summed)/(hal_PT14_arr_summed + fil_PT14_arr_summed + voi_PT14_arr_summed)

print(hal_frac.shape)

#####################
#plot mass fractions#
#####################

fig=plt.figure(figsize=(12,9))

ax1=fig.add_subplot(111)
ax1.set_xlabel('z',fontsize=15)
ax1.set_ylabel('mf',fontsize=15)
ax1.set_yscale('log')
ax1.set_ylim([6E-3,1])
ax1.set_xlim([0.,6.])

ax1.plot((allsnaps_redshifts.mean(axis=1))[0:13],hal_frac[0:13],label='halo')
ax1.plot((allsnaps_redshifts.mean(axis=1))[0:13],fil_frac[0:13],label='filament')
ax1.plot((allsnaps_redshifts.mean(axis=1))[0:13],voi_frac[0:13],label='void')


ax1.legend(fontsize=15)
plt.savefig('evolving_mass_fraction.png')

