import numpy as np
import scipy.special
from netCDF4 import Dataset

#--- This program was created to replicate how the effective radius
#--- of different hydrometeor species is calculated in the Morrison 
#--- Microphysical Parameterization

#--- Assumed Constants for the Morrison-Hail Scheme
rhow = 997    #--- Density of Liquid Water
rhoh = 900.   #--- Density of Hail (Assumed for this simulation)
rhos = 100.   #--- Density of Snow 
rhoi = 500.   #--- Density of Cloud Ice
ndcsnt = 250. #--- The assumed number of cloud droplets per cm^-3 (namelist parameter in CM1)

def cal_lamda(alpha,cx,ntx,qx):

   """
   Calculate the slope parameter of the Particle Size Distribution

   alpha, cx, ntx, qx, rhoa are floats

   Required arguments:
      alpha:   The shape parameter (can be float or array) 
      cx:      (pi/6) * hydrometeor density (kg m^-3) [float/array]
      ntx:     The number concentration (stones / kg^-1) [array]
      qx:      The mixing ratio  (kg/kg)

   Returns:
       lamda:  The slope parameter (m^-1)
   """

   gammaOne = scipy.special.gamma(1.+alpha)
   gammaFour = scipy.special.gamma(4.+alpha)
   lamda = ((gammaFour/gammaOne)*cx*ntx/qx)**(1./3.)
   lamda = np.where(qx>0.,lamda,0.)

   return lamda


def cal_eff_rad_morrison(qx,ntx,rhox,rhoa,hydro_type):
   """
   Calculate the effective radius of a hydrometeor speciec
   
   Required Arguments:
       qx: Hydrometeor mixing ratio (kg/kg) [nhydro,nz,ny,nx]
       ntx: Hydrometeor Number Concentration (kg^-1) [nhydro,nz,ny,nx]
       rhox: Hydrometeor density (kg m^-3)
       rhoa: Air Density (kg m^-3) [nz,ny,nx]
       hydro_type: A string for the hydrometeor type you desire
                   (cloud,rain,snow,hail,ice,snowice)

   Returns: Particle Effective Radius (Microns)
   """
   micron_conv = 1.e6  #--- Convert from m to microns
   QSMALL = 1.e-14     #--- Smallest allowed mass concentration

   #--- Sanity Check, make sure requested hydrometeor type exists
   hydro_type = hydro_type.lower() 
   if hydro_type in ['cloud','rain','snow','ice','hail','snowice']:
      print('The hydro type is ...',hydro_type)
   else:
      print('The desired hydrometeor type does not work...')


   #--- Most Hydrometeor Types assime alpha == 0 
   #--- except cloud water
   if 'cloud' in hydro_type:
      pgam = 0.0005714*(ntx[hydro_type]/1.E6*rhoa)+0.2714
      alpha = 1./(pgam**2)-1.
      alpha[alpha<2] = 2.
      alpha[alpha>10] = 10.
   else:
      alpha = 0.

   gammaThree = scipy.special.gamma(3.+alpha)
   gammaFour  = scipy.special.gamma(4.+alpha)

   #--- Creating a combined snow and ice category
   #--- Commonly ingested by radiation parameterizations
   if 'snowice' in hydro_type:
      #--- Snow
      cxs = (np.pi/6.)*rhox['snow']
      lamdas = cal_lamda(alpha,cxs,ntx['snow'],qx['snow'])
      lammins = 1./10.E-6
      lammaxs = 1./2000.E-6
      lamdas = np.where(lamdas>lammaxs,lammaxs,lamdas)
      lamdas = np.where(lamdas<lammins,lammins,lamdas)

      #--- Cloud Ice
      cxi = (np.pi/6.)*rhox['ice']
      lamdai = cal_lamda(alpha,cxi,ntx['ice'],qx['ice'])
      lammini = 1./(2.*125.E-6+100.E-6)
      lammaxi = 1./1.E-6
      lamdai = np.where(lamdai>lammaxi,lammaxi,lamdai)
      lamdai = np.where(lamdai<lammini,lammini,lamdai)

      #--- Create an effective radius for all situations where snow and cloud ice are present
      eff_rad = np.ones((qx['snow'].shape))*25.
      eff_rad = np.where((qx['snow']>QSMALL) & (qx['ice']>QSMALL),
                          micron_conv*(gammaFour/(2.*gammaThree))*((1./(lamdas**4.))+(1./(lamdai**4.)))/((1./(lamdai**3.))+(1./(lamdas**3.))),
                          eff_rad)
      eff_rad = np.where((qx['snow']>QSMALL) & (qx['ice']<QSMALL),micron_conv*gammaFour/(2.*lamdas*gammaThree),eff_rad)
      eff_rad = np.where((qx['ice']>QSMALL) & (qx['snow']<QSMALL),micron_conv*gammaFour/(2.*lamdai*gammaThree),eff_rad)

   else: 
      cx = (np.pi/6.)*rhox[hydro_type]
      lamda = cal_lamda(alpha,cx,ntx[hydro_type],qx[hydro_type])
      if 'cloud' in hydro_type:
         lammin = (alpha+1.)/60.E-6
         lammax = (alpha+1.)/1.E-6
      elif 'snow' in hydro_type:
         lammin = 1./10.E-6
         lammax = 1./2000.E-6
      elif 'hail' in hydro_type:
         lammin = 1./2000.E-6
         lammax = 1./20.E-6
      elif 'rain' in hydro_type:
         lammin = 1./2800.E-6
         lammax = 1./20.E-6
      elif 'ice' in hydro_type:
         lammin = 1./(2.*125.E-6+100.E-6)   
         lammax = 1./1.E-6
      lamda = np.where(lamda>lammax,lammax,lamda)
      lamda = np.where(lamda<lammin,lammin,lamda)
      eff_rad = np.where(qx[hydro_type]>QSMALL,micron_conv*gammaFour/(2.*lamda*gammaThree),25.)
   
   #--- Maximum effective radius for certain hydrometeor species in Morrison
   if hydro_type in ['snow','ice','snowice']:
      eff_rad[eff_rad>400.] = 400.

   return eff_rad


#---- Beginning Executable Code
#--- First Read In the Hydrometeor Variables and Define Values
filepath = '/scratch/jonathan.labriola/Large_Domain/QLCS/0.5km_Simulations/mod_snd/ens001_fric_rad_dry_nocorr_100km_10kmwave_jackson/cm1out_000024.nc'
dumpfile = Dataset(filepath)
rhoa = dumpfile.variables['rho'][0,:,:,:] 

#--- Create Dictionaries to Store Hydrometeor information
qx = {}   # kg kg^-1
ntx = {}  # kg^-1
rhox = {} # kg m^-3
#--- Rain
qx['rain'] = dumpfile.variables['qr'][0,:,:,:]
ntx['rain'] = dumpfile.variables['ncr'][0,:,:,:]
rhox['rain'] = rhow
#--- Snow
qx['snow']  = dumpfile.variables['qs'][0,:,:,:]
ntx['snow'] = dumpfile.variables['ncs'][0,:,:,:]
rhox['snow'] = rhos
#--- ice 
qx['ice'] = dumpfile.variables['qi'][0,:,:,:]
ntx['ice'] = dumpfile.variables['nci'][0,:,:,:]
rhox['ice'] = rhoi
#--- Hail
qx['hail'] = dumpfile.variables['qg'][0,:,:,:]
ntx['hail'] = dumpfile.variables['ncg'][0,:,:,:]
rhox['hail'] = rhoh
#---Cloud 
qx['cloud'] = dumpfile.variables['qc'][0,:,:,:]
ntx['cloud'] = (ndcsnt * 1E6)/rhoa 
rhox['cloud'] = rhow

#--- Calculate Effective Radii For Each Desired Species
effc = cal_eff_rad_morrison(qx,ntx,rhox,rhoa,'cloud')
print('Cloud Mean Effective Radius = ',np.mean(effc))
effr = cal_eff_rad_morrison(qx,ntx,rhox,rhoa,'rain')
print('Rain Mean Effective Radius = ',np.mean(effr))
effh = cal_eff_rad_morrison(qx,ntx,rhox,rhoa,'hail')
print('Hail Mean Effective Radius = ',np.mean(effh))
effs = cal_eff_rad_morrison(qx,ntx,rhox,rhoa,'snow')
print('Snow Mean Effective Radius = ',np.mean(effs))
effi = cal_eff_rad_morrison(qx,ntx,rhox,rhoa,'ice')
print('Ice Mean Effective Radius = ',np.mean(effi))
effsi = cal_eff_rad_morrison(qx,ntx,rhox,rhoa,'snowice')
print('Snow Ice Mean Effective Radius = ',np.mean(effsi))


