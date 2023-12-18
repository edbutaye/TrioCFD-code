import sys
import numpy as np
from numpy import *
from pylab import *


fileseq = ["FTD_TCL_parallel_BCtheta_couplingThermal_pb1_Diffusion_chaleur.face" \
           , "FTD_TCL_parallel_BCtheta_couplingThermal_pb1_Diffusion_chaleur_vapeur.face" \
           , "FTD_TCL_parallel_BCtheta_couplingThermal_pb2_Diffusion_chaleur.face"]

filepar = ["FTD_TCL_parallel_BCtheta_couplingThermal_PAR_pb1_Diffusion_chaleur.face" \
           , "FTD_TCL_parallel_BCtheta_couplingThermal_PAR_pb1_Diffusion_chaleur_vapeur.face" \
           , "FTD_TCL_parallel_BCtheta_couplingThermal_PAR_pb2_Diffusion_chaleur.face"]
name = ["liq_","vap_","sol_"]


if (len(sys.argv) != 2):
   raise Exception("One argument is compulsory to set, seq or par")

# print(f"post-procesing  {sys.argv[1]}")
linput = str(sys.argv[1])
seq = (linput=='seq')
  
if seq:
   myfile = fileseq
   ldir = 'SEQ/'
   lopdir = "RUN_SEQ/"
   npc = 1
   print("post-procesing seq")
else:
   myfile = filepar
   ldir = 'PAR/'
   npc = 6
   lopdir = "RUN_PAR/"
   print("post-procesing par")
   
for ncase, lfile in enumerate(myfile):
   f = open(lopdir+lfile)
   lines = f.readlines()
   nb_lines = len(lines)

   cnt =0
   for ibl, line in enumerate(lines):
      if 'temps' in line:
         cnt+=1
   nt=int(cnt/npc)
   nx=int((nb_lines-cnt)/nt)
   mat = zeros((nt,nx,5))
   it = -1
   ix = 0
   t = []
   for line in lines:
      if 'temps' in line:
         tt = float(line.split('temps')[1].replace(":", '').strip())
         if tt not in t:
            t.append(tt)
            it +=1
            ix = 0
               #print "Reading time t=", tt
      else:
         li = line.split(' ')
         x = li[4]
         y= li[6]
         s= li[8]
         phi=li[10] # flux_par_surface(W/m2)
         pui=li[12] # flux(W)
         mat[it,ix,0] = x
         mat[it,ix,1] = y
         mat[it,ix,2] = s
         mat[it,ix,3] = phi 
         mat[it,ix,4] = pui
         ix +=1
         pass
      pass
   print("Number of timesteps was correct: ", nt == len(t), " (nt=%d)"%nt)


   xx = mat[0,:,0]
   savetxt(ldir+name[ncase]+'xx.txt', xx)
   savetxt(ldir+name[ncase]+'mphit.txt', mat[:,:,3])
   savetxt(ldir+name[ncase]+'xphi.txt', np.c_[xx,mat[-1,:,3]])

   

   cum=mat*0
   cumsum(mat, axis=1, out=cum)
         #plot(xx,cum[-1,:,4])

   itdemi=int(nt/2)
   print("saving medium step : t_medium=",t[itdemi])
   savetxt(ldir+name[ncase]+'xphi_half.txt', np.c_[xx,mat[itdemi,:,3],cum[itdemi,:,4]])
   print("saving last step : t_end=",t[nt-1])
   savetxt(ldir+name[ncase]+'xphi_last.txt', np.c_[xx,mat[-1,:,3],cum[-1,:,4]])

   
liq_ = np.loadtxt(ldir+'liq_xphi_half.txt')
vap_ = np.loadtxt(ldir+'vap_xphi_half.txt')
sol_ = np.loadtxt(ldir+'sol_xphi_half.txt')

# Extract the second column (all rows, second column)
sum_flux = liq_[:, 1] + vap_[:, 1] +sol_[:, 1]
# check the 
conserve = np.all(sum_flux == 0)
if conserve:
   print(" :) heat flux CONSERVES at the fluid-solid connexion")
else:
   print(" :< ATTENTION: heat flux does NOT converve at fluid-solid connexion")
