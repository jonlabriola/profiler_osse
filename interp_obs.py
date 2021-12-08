import post.io
import argparse
import os
import numpy as np
from netCDF4 import Dataset
#from scipy import interpolate
from scipy.interpolate import RectBivariateSpline

obvarname = 'shs'
lev = 0
hgt = 1.0 
spec_hgt = True


#--- Opening the forecast file to interpolate Observations to
filepath = '/scratch/jonathan.labriola/osse/500m_NR/MORR_WoFS_blend/ens_all/ens001/mem001/cm1out_000002.nc'
dumpfile = Dataset(filepath,"r",fortmat="NETCDF4")
xh_coarse = post.io.cm1.readvar(dumpfile,'xh')
yh_coarse = post.io.cm1.readvar(dumpfile,'yh')
dumpfile.close()

for mindex,mem in enumerate(range(1,2)):
   for time in range(1,24):
      filepath = '/scratch/jonathan.labriola/osse/500m_NR/MORR_WoFS_blend/ens_all/obs%03d/cm1out_%06d.nc'%(mem,time)
      dumpfile = Dataset(filepath,"r",fortmat="NETCDF4")
      if obvarname in ['shs','sus','rain']:
         var_fine = post.io.cm1.readvar(dumpfile,obvarname)
      else:
         if spec_hgt: 
            var_fine = post.io.cm1.readvar(dumpfile,obvarname,interp_hgt=hgt)
         else:
            var_fine = post.io.cm1.readvar(dumpfile,obvarname,lev=lev)
      xh_fine =  post.io.cm1.readvar(dumpfile,'xh')
      yh_fine =  post.io.cm1.readvar(dumpfile,'yh')
      dumpfile.close()



      #--- Interpolate
      f = RectBivariateSpline(yh_fine, xh_fine, var_fine)
      var_coarse = f(yh_coarse,xh_coarse)
      outpath = '/scratch/jonathan.labriola/osse/500m_NR/MORR_WoFS_blend/ens_all/obs%03d/%s_3km_%06d.npz'%(mem,obvarname,time)
      print('Outputting %s'%outpath)
      np.savez(outpath, xh=xh_coarse, yh=yh_coarse,var=var_coarse)
