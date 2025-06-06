#! /usr/bin/python
# -*- coding:utf8 -*-

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import sys
import time
import multiprocessing
import gc

fontsize=20
plt.rc("font",size=fontsize)
plt.rc('font', family='serif')
plt.rc('text', usetex=True)
#matplotlib.use('TkAgg')

color=['#ff0000','#ff6600','#00ff00','#006600','#00ffff','#0000ff','#cc66ff','k']
fmt='os>*^hv'

clock=time.time()

##################
### PARAMETERS ###
##################

beta1=0.75
beta2=0.75
epsilon=0.9
rho0=5.
mag0=0.
gamma=0.5
LX=1024
init=0
rhomax=12
t0=5000
tmax=10000
NCPU=4
multi=True
movie=False

for arg in sys.argv[1:]:
	if "-beta1=" in arg:
		beta1=float(arg[7:])
	elif "-beta2=" in arg:
		beta2=float(arg[7:])
	elif "-epsilon=" in arg:
		epsilon=float(arg[9:])
	elif "-rho0=" in arg:
		rho0=float(arg[6:])
	elif "-mag0=" in arg:
		mag0=float(arg[6:])
	elif "-gamma=" in arg:
		gamma=float(arg[7:])
	elif "-LX=" in arg:
		LX=int(arg[4:])
	elif "-rhomax=" in arg:
		rhomax=float(arg[8:])
	elif "-t0=" in arg:
		t0=int(arg[4:])
	elif "-tmax=" in arg:
		tmax=int(arg[6:])
	elif "-init=" in arg:
		init=int(arg[6:])
	elif "-NCPU=" in arg:
		NCPU=int(arg[6:])
	elif "-movie" in arg:
		movie=True
	else:
		print("Bad Argument: ",arg)
		sys.exit(1)
		
if NCPU==1:
	multi=False
elif NCPU>1:
	multi=True
elif NCPU<1:
	print("Bad value of NCPU: ",NCPU)
	sys.exit(1)

DT=5
dpi=240

#############################
### CREATION OF SNAPSHOTS ###
#############################

def numberSnaps():
	j=0	
	SNAPS=[]
	while os.path.isfile('snapshots/figure_TSAIM_species_profile_beta1=%.8g_beta2=%.8g_epsilon=%.8g_rho0=%.8g_mag0=%.8g_gamma=%.8g_LX=%d_init=%d_%d.png'%(beta1,beta2,epsilon,rho0,mag0,gamma,LX,init,j)) and t0+j*DT<=tmax:
		SNAPS.append(j)
		j+=1
	return len(SNAPS)

def Snapshot(i):
	if not os.path.isfile('snapshots/figure_TSAIM_species_profile_beta1=%.8g_beta2=%.8g_epsilon=%.8g_rho0=%.8g_mag0=%.8g_gamma=%.8g_LX=%d_init=%d_%d.png'%(beta1,beta2,epsilon,rho0,mag0,gamma,LX,init,i)):
		t=t0+i*DT
		data=np.loadtxt('data_TSAIM_species_dynamics1d/TSAIM_species_profile_beta1=%.8g_beta2=%.8g_epsilon=%.8g_rho0=%.8g_mag0=%.8g_gamma=%.8g_LX=%d_init=%d_t=%.8g.txt'%(beta1,beta2,epsilon,rho0,mag0,gamma,LX,init,t))
		X=data[:,0]
		RHOA=data[:,1]
		RHOB=data[:,2]
		RHO=RHOA+RHOB
		
		fig=plt.figure(figsize=(7,6))
		gs=matplotlib.gridspec.GridSpec(1,1,width_ratios=[1],height_ratios=[1],left=0.12,right=0.96,bottom=0.12,top=0.94,hspace=0.1,wspace=0.1)

		ax=plt.subplot(gs[0,0])
		
		plt.plot(X,RHOA,'-r')
		plt.fill_between(X,0,RHOA,color='r',alpha=0.2)
		plt.plot(X,RHOB,'-b')
		plt.fill_between(X,0,RHOB,color='b',alpha=0.2)
		
		plt.xlim([0,LX])
		plt.ylim([0,rhomax])
		plt.xlabel('$x$')
		plt.ylabel('$\\rho_s(x)$')
		
		plt.xticks([0,0.25*LX,0.5*LX,0.75*LX,LX])
		ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('$%.8g$'))
		ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('$%.8g$'))
		
		plt.text(0,1.03*rhomax,'$t=%d$'%(t),ha='left',va='center',fontsize=17)
		plt.text(LX,1.03*rhomax,'$\\beta_1=%.8g$, $\\beta_2=%.8g$, $\\varepsilon=%.8g$, $\\rho_0=%.8g$, $m_0=%.8g$'%(beta1,beta2,epsilon,rho0,mag0),ha='right',va='center',fontsize=17)
		
		plt.savefig('snapshots/figure_TSAIM_species_profile_beta1=%.8g_beta2=%.8g_epsilon=%.8g_rho0=%.8g_mag0=%.8g_gamma=%.8g_LX=%d_init=%d_%d.png'%(beta1,beta2,epsilon,rho0,mag0,gamma,LX,init,i),dpi=dpi)
		plt.close()
		
		print('-snap=%d/%d -t=%d -rmin=%.4g -rmax=%.4g -tcpu=%d'%(i+1,Nsnap,t,RHO.min(),RHO.max(),time.time()-clock))
		del fig,data


os.system('mkdir -p snapshots/')

i=0

ARG=[]
while os.path.isfile('data_TSAIM_species_dynamics1d/TSAIM_species_profile_beta1=%.8g_beta2=%.8g_epsilon=%.8g_rho0=%.8g_mag0=%.8g_gamma=%.8g_LX=%d_init=%d_t=%.8g.txt'%(beta1,beta2,epsilon,rho0,mag0,gamma,LX,init,t0+i*DT)) and t0+i*DT<=tmax:
	ARG.append(i)
	i+=1
	
Nsnap=len(ARG)
Nsnapdone=numberSnaps()
	
print('%d Snapshots (%d already created)'%(Nsnap,Nsnapdone))

if multi:
	pool=multiprocessing.Pool(NCPU)
	results=pool.imap_unordered(Snapshot,ARG)
	pool.close()
	pool.join()
else:
	for i in ARG:
		Snapshot(i)
		
Nsnapdone=numberSnaps()

#############################
### CREATION OF THE MOVIE ###
#############################

if movie:
	os.system('mkdir -p movies')	
	os.system('ffmpeg -v quiet -stats -y -r 25/1 -i snapshots/figure_TSAIM_species_profile_beta1=%.8g_beta2=%.8g_epsilon=%.8g_rho0=%.8g_mag0=%.8g_gamma=%.8g_LX=%d_init=%d_%%01d.png -c:v h264 -r 25 -crf 30 -s %dx%d movies/movie_TSAIM_species_hydro_profile_beta1=%.8g_beta2=%.8g_epsilon=%.8g_rho0=%.8g_mag0=%.8g_gamma=%.8g_LX=%d_init=%d.mp4'%(beta1,beta2,epsilon,rho0,mag0,gamma,LX,init,7*dpi,6*dpi,beta1,beta2,epsilon,rho0,mag0,gamma,LX,init))
	if Nsnap>(tmax-t0)/DT and Nsnap==Nsnapdone:
		os.system('rm snapshots/figure_TSAIM_species_profile_beta1=%.8g_beta2=%.8g_epsilon=%.8g_rho0=%.8g_mag0=%.8g_gamma=%.8g_LX=%d_init=%d_*.png'%(beta1,beta2,epsilon,rho0,mag0,gamma,LX,init))

print('OK - time=%d sec'%(time.time()-clock))
