import os
from netCDF4 import Dataset
import readvar,verification,ns_colormap,jl_modules
import argparse
import matplotlib,verification
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import numpy as np
import plotting_lib
matplotlib.use('Agg')
#----------------------------------------------------#
#            PARSE USER COMMAND LINE INPUT           #
#----------------------------------------------------#

parser = argparse.ArgumentParser()
parser.add_argument("file_loc",type=str,help='microphysical option 1,2,3,4')
parser.add_argument("var",type=str,help='The var to be calculated uh, dbz')
parser.add_argument("--grdbas",type=str,help='The gridbse file to be read in')
parser.add_argument("--name",type=str,default='',help='Picture titles')
parser.add_argument("--threshold",type=float,default=45,help='The threshold to verify (Default is 45 dBZ)')
parser.add_argument("--kernel",type=float,default=2.5,help='Kernel Neighborhood Radius (in gridpoints)')
parser.add_argument("--tstart",type=int,default=1,help='The desired verification method, Reliability, AUC, Performance')
parser.add_argument("--tstop",type=int,default=60,help='The desired verification method, Reliability, AUC, Performance')
parser.add_argument("--nen",type=int,default=10,help='The desired verification method, Reliability, AUC, Performance')
parser.add_argument("--lev",type=int,default=1,help='The plotted level')
parser.add_argument("--interp_hgt",type=int,default=-1,help='The plotted level')
parser.add_argument("--xmin",default=-1,type=int,help =' Minimum X value')
parser.add_argument("--xmax",default=-1,type=int,help =' Maximum X value')
parser.add_argument("--ymin",default=-1,type=int,help =' Minimum Y value')
parser.add_argument("--ymax",default=-1,type=int,help =' Maximum Y value')
parser.add_argument("--local",action='store_true',help='Save file locally')
parser.add_argument("--hov",action='store_true',help='Hoevmoller')
parser.add_argument("--elreno",action='store_true',help='El Reno Case')
parser.add_argument("--linex",default='-1,-1',type=str,help='Use if you want to make a line')
parser.add_argument("--liney",default='-1,-1',type=str,help='Use if you want to make a line')
parser.add_argument("--obs",default='N/A',type=str,help='The path towards the observation file')
parser.add_argument("--refl_overlay",default=0.,type=float,help='Reflectivity Overlay threshold')
parser.add_argument("--rst",action='store_true',help='Restart')
parser.add_argument("--paintball",action='store_true',help='Paintball')
parser.add_argument("--pmmean",action='store_true',help='Probability Matched Mean')
parser.add_argument("--percentile",action='store_true',help='Probability Matched Mean')
parser.add_argument("--mean",action='store_true',help='Mean')
parser.add_argument("--max",action='store_true',help='Max')
parser.add_argument("--reliability",action='store_true',help='Reliability')
parser.add_argument("--variance",action='store_true',help='Plot Variance')
parser.add_argument("--plot_title",action='store_true',help='Add a time stamp to title')
parser.add_argument('--small_ens',action='store_true',help='Small Ens')
parser.add_argument("--subtract",action='store_true',help='Subtract Value')
parser.add_argument("--hist",action='store_true',help='Plot observation')
parser.add_argument("--save",action='store_true',help='save file')
arguments = parser.parse_args()

#Get Grid Dimensions of Output
#path = '%s/ENF001/cm1out_000001.nc'%arguments.file_loc
#dumpfile = Dataset(path,"r",fortmat="NETCDF4")
#t_dim #--- Number of Time Steps (For HOV Plots)

dumpfile = Dataset(arguments.grdbas,"r",fortmat="NETCDF4") 
grdbas = readvar.readGrdbas(dumpfile,rst=arguments.rst,arguments=arguments)
grdbas.nt = (arguments.tstop - arguments.tstart) + 1

#nx = readvar.open_file(dumpfile,'nx')
#ny = readvar.open_file(dumpfile,'ny')
#nz = readvar.open_file(dumpfile,'nz')
#nt = (arguments.tstop - arguments.tstart) + 1
#x_h = readvar.open_file(dumpfile,'xh')
#y_h = readvar.open_file(dumpfile,'yh')
#z_h = readvar.open_file(dumpfile,'zh')
#dx = (x_h[1] - x_h[0])*1000.


if arguments.xmin > -1:
   grdbas.nx = arguments.xmax - arguments.xmin
   grdbas.ny = arguments.ymax - arguments.ymin
   xmax = arguments.xmax
   xmin = arguments.xmin
   ymax = arguments.ymax
   ymin = arguments.ymin
else:
   xmax = grdbas.nx
   xmin = 0
   ymax = grdbas.ny
   ymin = 0
 
if ymax - ymin <= 1:
   ny = grdbas.nz
   nx = grdbas.nx
   vert_plot = True
elif xmax - xmin <= 1:
   ny = grdbas.nz
   nx = grdbas.ny
   vert_plot = True
else:
   ny = grdbas.ny
   nx = grdbas.nx
   vert_plot = False  
print('nz = %s, ny = %s, nz = %s'%(str(grdbas.nx),str(grdbas.ny),str(grdbas.nz)))





#--- Read in Files if Verifying HOV Diagrams
if arguments.hov:
   var = np.zeros((arguments.nen,grdbas.ny,grdbas.nt))
   for mindex,member in enumerate(range(1,arguments.nen+1)):
      for tindex,time in enumerate(range(arguments.tstart,arguments.tstop+1)):
         try:
            if arguments.rst:
               path = arguments.file_loc+'/LES_%02d/cm1rst_%06d.nc'%(member,time)
            else:
              path = arguments.file_loc+'/LES_%02d/cm1out_%06d.nc'%(member,time)
            dumpfile = Dataset(path,"r",fortmat="NETCDF4")
         except:
            if arguments.rst:
               path = arguments.file_loc+'/mem%03d/cm1rst_%06d.nc'%(member,time)
            else:
               try:
                  path = arguments.file_loc+'/obs%03d/cm1out_%06d.nc'%(member,time)
                  dumpfile = Dataset(path,"r",fortmat="NETCDF4")
               except:
                  try:
                     path = arguments.file_loc+'/mem%03d/cm1out_%06d.nc'%(member,time)
                     dumpfile = Dataset(path,"r",fortmat="NETCDF4")
                  except:
                     path = arguments.file_loc+'/mem%03d_periodic/cm1out_%06d.nc'%(member,time)
                     dumpfile = Dataset(path,"r",fortmat="NETCDF4")
          
         var_tmp = np.amax(readvar.open_file(dumpfile,arguments.var,xmax=xmax,xmin=xmin,ymax=ymax,ymin=ymin),axis=1) 
         var[mindex,:,tindex] = var_tmp 
         dumpfile.close()

#--- Read in Files if Verifying Full Domain
else:
   if arguments.small_ens:
      members =  np.array([2,3,4,5,6,7,8,9,10,11,13,14,16,17,18,19,22,24,25,26,27,28,29,30,31,33])
      arguments.nen = 26
   else:
      members = np.arange(1,arguments.nen+1)
 
   #if arguments.ymax == arguments.ymin:
   #   var = np.zeros((arguments.nen,grdbas.nz,grdbas.nx))
   #elif arguments.xmax == arguments.xmin:
   #   var = np.zeros((arguments.nen,grdbas.nz,grdbas.ny))
   #else:
   var = np.zeros((arguments.nen,ny,nx))
   varTwo = np.zeros(var.shape)
   if vert_plot:  # -- JDL new
      if arguments.elreno:
         var_1km = np.zeros((arguments.nen,600,250))
      else:
         var_1km = np.zeros((arguments.nen,800,200))
      var_2km = np.zeros(var_1km.shape)
      var_5km = np.zeros(var_1km.shape)
   for mindex,member in enumerate(members):
      for tindex in range(arguments.tstart,arguments.tstop+1):
         #--- Different Observation Directories
         if os.path.isdir(arguments.file_loc+'/LES_%02d/'%member):
            base_path = arguments.file_loc+'/LES_%02d/'%member
         elif os.path.isdir(arguments.file_loc+'/dir_model%03d/'%member):
            base_path = arguments.file_loc+'/dir_model%03d/'%member
         elif os.path.isdir(arguments.file_loc+'/obs%03d/'%member):
            base_path = arguments.file_loc+'/obs%03d/'%member
         elif os.path.isdir(arguments.file_loc+'/mem%03d/'%member):
            base_path = arguments.file_loc+'/mem%03d/'%member
         elif os.path.isdir(arguments.file_loc+'/ens%03d/'%member):
            base_path = arguments.file_loc+'/ens%03d/'%member
         elif os.path.isdir(arguments.file_loc+'/mem%03d_periodic/'%member):
            base_path = arguments.file_loc+'/mem%03d_periodic/'%member
         else:
            print('No Directories')

         if arguments.rst:
            path = '%s/cm1out_rst_%06d.nc'%(base_path,tindex)
         else:
            path = '%s/cm1out_%06d.nc'%(base_path,tindex)
         dumpfile = Dataset(path,"r",fortmat="NETCDF4")

         print('Opening ... ',path)
         if vert_plot :
            var_tmp = readvar.open_file(dumpfile,arguments.var,xmax=xmax,xmin=xmin,ymax=ymax,ymin=ymin)
            if arguments.refl_overlay > 0:
               var_tmp_two = readvar.open_file(dumpfile,'dbz',xmax=xmax,xmin=xmin,ymax=ymax,ymin=ymin)
            if arguments.elreno:
               temp = readvar.open_file(dumpfile,arguments.var,interp_hgt=1)
               var_1km_tmp = readvar.open_file(dumpfile,arguments.var,xmin=750,xmax=1000,ymin=0,ymax=600,interp_hgt=1)
               var_2km_tmp =  readvar.open_file(dumpfile,arguments.var,xmin=750,xmax=1000,ymin=0,ymax=600,interp_hgt=2)
               var_5km_tmp =  readvar.open_file(dumpfile,arguments.var,xmin=750,xmax=1000,ymin=0,ymax=600,interp_hgt=0.5)
            else:
               var_1km_tmp = readvar.open_file(dumpfile,arguments.var,xmin=600,xmax=800,ymin=0,ymax=800,interp_hgt=1)
               var_2km_tmp =  readvar.open_file(dumpfile,arguments.var,xmin=600,xmax=800,ymin=0,ymax=800,interp_hgt=2)
               var_5km_tmp =  readvar.open_file(dumpfile,arguments.var,xmin=600,xmax=800,ymin=0,ymax=800,interp_hgt=0.5)
         else:
            if arguments.interp_hgt > 0:
               var_tmp = readvar.open_file(dumpfile,arguments.var,xmax=xmax,xmin=xmin,ymax=ymax,ymin=ymin,interp_hgt=arguments.interp_hgt)
            else:
               var_tmp = readvar.open_file(dumpfile,arguments.var,xmax=xmax,xmin=xmin,ymax=ymax,ymin=ymin,lev=arguments.lev)
            if arguments.refl_overlay > 0:
               var_tmp_two = readvar.open_file(dumpfile,'dbz',xmax=xmax,xmin=xmin,ymax=ymax,ymin=ymin,lev=-1)



         if arguments.variance:
            var[mindex] = var_tmp
         else:
            var_tmp = jl_modules.smoothObs(var_tmp)
            var[mindex] = np.where(var[mindex] > var_tmp, var[mindex],var_tmp)

         if vert_plot:
            var_1km[mindex] = var_1km_tmp
            var_2km[mindex] = var_2km_tmp
            var_5km[mindex] = var_5km_tmp
         if arguments.refl_overlay > 0:
            varTwo[mindex ] = np.where(varTwo[mindex] > var_tmp_two, varTwo[mindex],var_tmp_two)
         dumpfile.close()

      #--- Create the Subtraction
      if arguments.subtract:
         path = '%s/cm1out_rst_000001.nc'%base_path
         dumpfile = Dataset(path,"r",fortmat="NETCDF4")
         if arguments.interp_hgt > 0:
            var_tmp_sub = readvar.open_file(dumpfile,arguments.var,xmax=xmax,xmin=xmin,ymax=ymax,ymin=ymin,interp_hgt=arguments.interp_hgt)
         else:
            var_tmp_sub = readvar.open_file(dumpfile,arguments.var,xmax=xmax,xmin=xmin,ymax=ymax,ymin=ymin,lev=arguments.lev)
         var[mindex] += (-1*var_tmp_sub)

varTwo=np.mean(varTwo,axis=0)
varTwo_threshold = arguments.refl_overlay

if arguments.mean or arguments.pmmean or arguments.max > 0. or arguments.percentile:
   if arguments.mean:
      var = np.mean(var,axis=0)
   elif arguments.pmmean:
      var = verification.pmMean(var) 
   elif arguments.percentile:
      var = np.percentile(var,90,axis=0)
   elif arguments.max:
      var = np.amax(var,axis=0)
elif arguments.variance:
   var = np.var(var,axis=0) 
   print('Max variance == ',np.amax(var))
   print('Mean variance == ',np.mean(var))
   if vert_plot:
      var_1km = np.var(var_1km,axis=0)
      var_2km = np.var(var_2km,axis=0)
      var_5km = np.var(var_5km,axis=0)
else:
   #--- NEP
   if arguments.hov:
      binary_ens = np.zeros((arguments.nen,grdbas.ny,grdbas.nt))
   else:  
      binary_ens = np.zeros((arguments.nen,ny,nx))
   for ens in range(0,arguments.nen):
      if arguments.var == 'T': # JDL
         binary_ens[ens] = verification.localThresh(-1.*var[ens]+1E-4,arguments.kernel,-1.*arguments.threshold)
      else:
         binary_ens[ens] = verification.localThresh(var[ens]+1E-4,arguments.kernel,arguments.threshold)
   binary_ens[binary_ens>0.0] = 1.0
   print('The max of binary ens is ',np.amax(binary_ens))
   if arguments.paintball:
      for index in range(0,arguments.nen):
         binary_ens[index] = np.where(binary_ens[index]==1,index+1,np.nan)
      var = binary_ens
   else:
      NMEP = np.sum(binary_ens,axis=0)/float(arguments.nen)
      var = gaussian_filter(NMEP,sigma = arguments.kernel)


figure = plt.figure()

#--- Time To Get the Observation
#--- Observations can be on a different grid for plotting purposes
#--- Observations must still be on the same grid
if arguments.obs != 'N/A' or arguments.reliability:
   if arguments.tstart == arguments.tstop:
      try:
         #--- Read from saved file
         path = '%s/'
      except:
         #--- Read from model output file
         path = '%s/cm1out_%06d.nc'%(arguments.obs,arguments.tstart)
         print('Opening ... %s'%path)
         dumpfile = Dataset(path,"r",fortmat="NETCDF4")
         grdbas_obs = readvar.readGrdbas(dumpfile)
         xh_obs = readvar.open_file(dumpfile,'xh')
         yh_obs = readvar.open_file(dumpfile,'yh')
         dx_obs =  (xh_obs[1] - xh_obs[0])
         dx_conversion = grdbas.dx / dx_obs  #--- What is the difference
         xmin_obs = int(xmin * dx_conversion)
         xmax_obs = int(xmax * dx_conversion)
         ymin_obs = int(ymin * dx_conversion)
         ymax_obs = int(ymax * dx_conversion) 
         ny_obs = ymax_obs - ymin_obs
         nx_obs = xmax_obs - xmin_obs


      if arguments.var == 'shs':
         #obs_threshold = 300
         ob_var = 'sws'
         obs_threshold = 33.45
      else:
         ob_var = 'sws'
         obs_threshold = 33.45

      if arguments.interp_hgt > 0:
         varTwo = readvar.open_file(dumpfile,ob_var,interp_hgt=arguments.interp_hgt)
      else:
         varTwo = readvar.open_file(dumpfile,ob_var,lev=arguments.lev)

      if arguments.subtract:
         sub_path = '%s/cm1out_000001.nc'%(arguments.obs)
         print('Opening ... %s'%sub_path)
         sub_dumpfile = Dataset(sub_path,"r",fortmat="NETCDF4")
         if arguments.interp_hgt > 0:
            sub_var = readvar.open_file(sub_dumpfile,arguments.var,interp_hgt=arguments.interp_hgt)
         else:
            sub_var = readvar.open_file(sub_dumpfile,arguments.var,lev=arguments.lev)
         varTwo += (-1*sub_var)

      print('JDL the new threshold is ... ',obs_threshold)
      varTwo_threshold = obs_threshold #arguments.threshold
      kernel_obs = arguments.kernel * int(dx_conversion)
      varTwo = jl_modules.smoothObs(varTwo)
      varTwo = verification.localThresh(varTwo,arguments.kernel*dx_conversion,varTwo_threshold) 
      varTwo[varTwo>0] = varTwo_threshold+1.
   else:
      varTwo = np.zeros(var.shape)
      pass #--- write eventually to get many obs times
   xx_obs,yy_obs = np.meshgrid(xh_obs[xmin_obs:xmax_obs],yh_obs[ymin_obs:ymax_obs]) #,xh_obs[xmin_obs:xmax_obs])
   plt.contour(xx_obs,yy_obs,varTwo[ymin_obs:ymax_obs,xmin_obs:xmax_obs],[0.,arguments.threshold],colors='k')
#--- Titles
if arguments.reliability:
    sample_climo,no_skill,obs_frequency,bin_centers,bin_climo,reliability_bins,bin_count = verification.createReliability(var,varTwo,arguments.threshold)
    BSS = verification.BSS(var,varTwo,arguments.threshold)
    #plotting_lib.plotReliability(sample_climo,no_skill,obs_frequency,bin_centers,bin_climo,linecolor = colors[index],label='%sM (%0.3f)'%(res,BSS))



varname = arguments.var
if arguments.mean:
   output_dir = '/work/jonathan.labriola/python_scripts/plots/ensemble/mean/' 
   title = 'Mean_%s_%s_%02d_%02d_lev%03d'%(arguments.name,arguments.var,arguments.tstart,arguments.tstop,arguments.lev)
elif arguments.max:
   output_dir = '/work/jonathan.labriola/python_scripts/plots/ensemble/max/'
   title = 'Max_%s_%s_%02d_%02d_lev%03d'%(arguments.name,arguments.var,arguments.tstart,arguments.tstop,arguments.lev)
elif arguments.pmmean:
   output_dir = '/work/jonathan.labriola/python_scripts/plots/ensemble/pmmean/'
   title = 'pmmean_%s_%s_%02d_%02d_lev%03d'%(arguments.name,arguments.var,arguments.tstart,arguments.tstop,arguments.lev)
elif arguments.percentile:
   output_dir = '/work/jonathan.labriola/python_scripts/plots/ensemble/percentile/'
   title = 'percentile_%s_%s_%02d_%02d_lev%03d'%(arguments.name,arguments.var,arguments.tstart,arguments.tstop,arguments.lev)
elif arguments.variance:
   output_dir = '/work/jonathan.labriola/python_scripts/plots/ensemble/variance/'
   title = 'Variance_%s_%s_%02d_%02d_lev%03d'%(arguments.name,arguments.var,arguments.tstart,arguments.tstop,arguments.lev)
   varname = arguments.var+'_var'
elif arguments.paintball:
   output_dir = '/work/jonathan.labriola/python_scripts/plots/ensemble/paintball/'
   title = 'paintball_%s_%s_%02d_%02d_lev%03d'%(arguments.name,arguments.var,arguments.tstart,arguments.tstop,arguments.lev)
elif arguments.reliability:
   output_dir = '/work/jonathan.labriola/python_scripts/plots/ensemble/reliability/'
   varname = 'reliability'
   title = 'reliability_%s_%s_%02d_%02d_lev%03d'%(arguments.name,arguments.var,arguments.tstart,arguments.tstop,arguments.lev) 
else: #--- NEP
   output_dir = '/work/jonathan.labriola/python_scripts/plots/ensemble/NEP/'
   varname = 'NEP'
   if arguments.hov:
      title = 'NEP_%s_%s_%02d_%02d_thresh%02d_HOV_lev%03d'%(arguments.name,arguments.var,arguments.tstart,arguments.tstop,arguments.threshold,arguments.lev)
   else:
      title = 'NEP_%s_%s_%02d_%02d_thresh%02d_lev%03d'%(arguments.name,arguments.var,arguments.tstart,arguments.tstop,arguments.threshold,arguments.lev)


if arguments.local:
   output_dir = './'
if vert_plot:
   title=title+'_vert'

if arguments.linex:
   lines = arguments.linex.split(',')
   if int(lines[0]) >= 0:
      linex = [grdbas.xh[int(lines[0])],grdbas.xh[int(lines[1])]]
      lines = arguments.liney.split(',')
      liney = [grdbas.yh[int(lines[0])],grdbas.yh[int(lines[1])]] 
   else:
      linex = [np.nan,np.nan]
      liney = [np.nan,np.nan]
   

outpath = '%s/%s.png'%(output_dir,title)
#print(np.array([0.,arguments.refl_overlay]))
print('Second Threshold = ',varTwo_threshold)
print('varTwo shape = ',varTwo.shape)
print('var shape = ',var.shape)

if arguments.save:
   if arguments.variance:
      if vert_plot and ymax - ymin <= 1 :
         #--- JDL New Addition
         mean_variance = np.zeros((3))
         mean_variance[0] = np.mean(var_5km)
         mean_variance[1] = np.mean(var_1km)
         mean_variance[2] = np.mean(var_2km)
         np.savez('variance_%s_%02d_%02d_%s_vert.npz'%(arguments.var,arguments.tstart,arguments.tstop,arguments.name), field=var, xh=grdbas.xh[xmin:xmax], yh=grdbas.zh, mean=mean_variance)
      elif vert_plot and xmax -xmin <= 1 :
         mean_variance = np.zeros((2))
         mean_variance[0] = np.mean(var_5km)
         mean_variance[1] = np.mean(var_1km)
         mean_variance[2] = np.mean(var_2km)
         np.savez('variance_%s_%02d_%02d_%s_vert.npz'%(arguments.var,arguments.tstart,arguments.tstop,arguments.name), field=var, xh=grdbas.yh[ymin:ymax], yh=grdbas.zh, mean=mean_variance)
      else:
         np.savez('variance_%s_%02d_%02d_%s.npz'%(arguments.var,arguments.tstart,arguments.tstop,arguments.name), field=var, xh=grdbas.xh, yh=grdbas.yh)
   elif arguments.paintball:
      np.savez('paintball_%02d_%02d_%s.npz'%(arguments.tstart,arguments.tstop,arguments.name),field=var,xh=grdbas.xh[arguments.xmin:arguments.xmax],yh=grdbas.yh[arguments.ymin:arguments.ymax])
   else:
      np.savez('NEP_%02d_%02d_%s.npz'%(arguments.tstart,arguments.tstop,arguments.name), field=var, xh=grdbas.xh, yh=grdbas.yh)

   #fn = '%s.nc'%title
   #ds = Dataset(fn, 'w', format='NETCDF4')
   #yh = ds.createDimension('yh', ymax-ymin)
   #xh = ds.createDimension('xh', xmax-xmin)
   #yh = ds.createVariable('yh', 'f4', ('yh',))
   #xh = ds.createVariable('xh', 'f4', ('xh',))
   #probs = ds.createVariable('probs', 'f4', ('yh', 'xh',))
   #probs.units = 'Unknown'
   #yh[:] = grdbas.yh
   #xh[:] = grdbas.xh
   #probs[:] = var
plt.tight_layout()

if arguments.hist and varname =='NEP':
   outpath = '%s/%s_hist.png'%(output_dir,title)
   bins = np.arange(0.05,1.05,0.05)
   weights = np.zeros((len(var.flatten())))
   weights[:] = grdbas.dx**2
   plt.hist(var.flatten(),bins = bins, weights = weights,histtype='step')
   plt.xlim([.05,1.0])
   plt.ylim([0.,1000.])
   plt.savefig(outpath)
elif arguments.reliability:
   plotting_lib.plotReliability(sample_climo,no_skill,obs_frequency,bin_centers,bin_climo,linecolor = 'r',label='BSS = (%0.3f)'%(BSS))   
   plt.legend()
   plt.savefig(outpath)
else:
   outpath = '%s/%s.png'%(output_dir,title)
   if arguments.paintball:
      plotting_lib.plotVar(var,'paintball',grdbas=grdbas,outpath=outpath,contour_output=varTwo,contours= np.array([varTwo_threshold]),linex=linex,liney=liney)
   else:
      print('JDL the varname is = ',varname)
      #plotting_lib.plotVar(var,varname,grdbas=grdbas,outpath=outpath,contour_output=varTwo,contours= np.array([varTwo_threshold]),linex=linex,liney=liney)
      if arguments.obs != 'N/A':
         plotting_lib.plotVar(var,varname,grdbas=grdbas,outpath=outpath,contour_output=varTwo,contours= np.array([varTwo_threshold]),linex=linex,liney=liney,grdbas_obs=grdbas_obs)
      else:
         plotting_lib.plotVar(var,varname,grdbas=grdbas,outpath=outpath,contour_output=varTwo,contours= np.array([varTwo_threshold]),linex=linex,liney=liney)
      
      #plotting_lib.plotVar(var,varname,grdbas=grdbas,outpath=outpath,linex=linex,liney=liney)


#plotting_lib.plotVar(var,varname,grdbas=grdbas,outpath=outpath,contour_output=varTwo,contours= np.array([varTwo_threshold]))

#if arguments.plot_title: #--- Read times
#   hdftime = (arguments.tstart - 1) * 300.
#   plt.title('%s s'%str(int(hdftime)))

#print('Saving... %s/%s'%(output_dir,title))
#plt.savefig('%s/%s'%(output_dir,title))
