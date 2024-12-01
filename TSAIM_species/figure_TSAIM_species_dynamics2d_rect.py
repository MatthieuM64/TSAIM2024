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
beta2=1.5
D=1.
v=1.8
theta=0.
rho0=10.
mag0=0.
gamma=0.5
LX=1024
LY=128
init=2
ran=0
rhomax=16.
t0=0.
tmax=600000
NCPU=4
multi=True
movie=False

for arg in sys.argv[1:]:
	if "-beta1=" in arg:
		beta1=float(arg[7:])
	elif "-beta2=" in arg:
		beta2=float(arg[7:])
	elif "-D=" in arg:
		D=float(arg[3:])
	elif "-v=" in arg:
		v=float(arg[3:])
	elif "-theta=" in arg:
		theta=float(arg[7:])
	elif "-rho0=" in arg:
		rho0=float(arg[6:])
	elif "-mag0=" in arg:
		mag0=float(arg[6:])
	elif "-gamma=" in arg:
		gamma=float(arg[7:])
	elif "-LX=" in arg:
		LX=int(arg[4:])
	elif "-LY=" in arg:
		LY=int(arg[4:])
	elif "-gamma=" in arg:
		gamma=float(arg[7:])
	elif "-rhomax=" in arg:
		rhomax=float(arg[8:])
	elif "-t0=" in arg:
		t0=int(arg[4:])
	elif "-tmax=" in arg:
		tmax=int(arg[6:])
	elif "-init=" in arg:
		init=int(arg[6:])
	elif "-ran=" in arg:
		ran=int(arg[5:])
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
	
dpi=240
DT=150
DeltaT=1./(4*D+theta*v+np.exp(2*beta1)+gamma*np.exp(2*beta2))

#############################
### CREATION OF SNAPSHOTS ###
#############################

def Snapshot(i):
	if not os.path.isfile('snapshots/figure_TSAIM_species_density_beta1=%.8g_beta2=%.8g_D=%.8g_v=%.8g_theta=%.8g_rho0=%.8g_mag0=%.8g_gamma=%.8g_LX=%d_LY=%d_init=%d_ran=%d_%d.png'%(beta1,beta2,D,v,theta,rho0,mag0,gamma,LX,LY,init,ran,i)):
		t=t0+i*DT
		RHO=np.loadtxt('data_TSAIM_species_dynamics2d/TSAIM_species_density_beta1=%.8g_beta2=%.8g_D=%.8g_v=%.8g_theta=%.8g_rho0=%.8g_mag0=%.8g_gamma=%.8g_LX=%d_LY=%d_init=%d_ran=%d_t=%d.txt'%(beta1,beta2,D,v,theta,rho0,mag0,gamma,LX,LY,init,ran,t))
		
		x=np.linspace(0,LX,LX)
		y=np.linspace(0,LY,LY)

		fig=plt.figure(figsize=(8,1.5))
		gs=matplotlib.gridspec.GridSpec(1,1,width_ratios=[1],height_ratios=[1],left=0.07,right=1.09, bottom=0.225,top=0.775,hspace=0.1,wspace=0.1)

		ax=plt.subplot(gs[0,0])
		
		cmap=plt.get_cmap('bwr')
		plt.pcolormesh(x,y,RHO,vmin=-rhomax,vmax=rhomax,rasterized=True,cmap=cmap)
		cb=plt.colorbar(ticks=[-rhomax,0,rhomax],aspect=5)
		cb.ax.set_yticklabels(['$%.4g$'%rhomax,'$0$','$%.4g$'%rhomax])
		cb.solids.set_rasterized(True)
		
		#plt.axis('equal')
		plt.xlim([0,LX])
		plt.ylim([0,LY])
		plt.xticks([0,0.25*LX,0.5*LX,0.75*LX,LX])
		plt.yticks([0,0.5*LY,LY])
		ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('$%.8g$'))
		ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('$%.8g$'))
		
		plt.text(0,1.2*LY,'$t=%d$'%(t*DeltaT),ha='left',va='center')
		plt.text(LX,1.2*LY,'$\\beta_1=%.8g$, $\\beta_2=%.8g$, $D=%.8g$, $v=%.8g$, $\\rho_0=%.8g$'%(beta1,beta2,D,v,rho0),ha='right',va='center')
		
		plt.savefig('snapshots/figure_TSAIM_species_density_beta1=%.8g_beta2=%.8g_D=%.8g_v=%.8g_theta=%.8g_rho0=%.8g_mag0=%.8g_gamma=%.8g_LX=%d_LY=%d_init=%d_ran=%d_%d.png'%(beta1,beta2,D,v,theta,rho0,mag0,gamma,LX,LY,init,ran,i),dpi=dpi)
		plt.close()
		
		print('-snap=%d/%d -t=%d -rhomin=%.4g -rhomax=%.4g -tcpu=%d'%(i+1,Nsnap,t,RHO.min(),RHO.max(),time.time()-clock))
		del fig,RHO


os.system('mkdir -p snapshots/')

i=0
ARG=[]
while os.path.isfile('data_TSAIM_species_dynamics2d/TSAIM_species_density_beta1=%.8g_beta2=%.8g_D=%.8g_v=%.8g_theta=%.8g_rho0=%.8g_mag0=%.8g_gamma=%.8g_LX=%d_LY=%d_init=%d_ran=%d_t=%d.txt'%(beta1,beta2,D,v,theta,rho0,mag0,gamma,LX,LY,init,ran,t0+i*DT)) and t0+i*DT<=tmax:
	ARG.append(i)
	i+=1
	
Nsnap=len(ARG)
	
print('%d Snapshots'%len(ARG))
if multi:
	pool=multiprocessing.Pool(NCPU)
	results=pool.imap_unordered(Snapshot,ARG[::-1])
	pool.close()
	pool.join()
else:
	for i in ARG:
		Snapshot(i)

#############################
### CREATION OF THE MOVIE ###
#############################

if movie:
	os.system('mkdir -p movies')	
	os.system('ffmpeg -v quiet -stats -y -r 25/1 -i snapshots/figure_TSAIM_species_density_beta1=%.8g_beta2=%.8g_D=%.8g_v=%.8g_theta=%.8g_rho0=%.8g_mag0=%.8g_gamma=%.8g_LX=%d_LY=%d_init=%d_ran=%d_%%01d.png -c:v h264 -r 25 -crf 30 -s %dx%d movies/movie_TSAIM_species_density_beta1=%.8g_beta2=%.8g_D=%.8g_v=%.8g_theta=%.8g_rho0=%.8g_mag0=%.8g_gamma=%.8g_LX=%d_LY=%d_init=%d_ran=%d.mp4'%(beta1,beta2,D,v,theta,rho0,mag0,gamma,LX,LY,init,ran,8*dpi,1.5*dpi,beta1,beta2,D,v,theta,rho0,mag0,gamma,LX,LY,init,ran))
	if Nsnap>(tmax-t0)/DT:
		os.system('rm snapshots/figure_TSAIM_species_density_beta1=%.8g_beta2=%.8g_D=%.8g_v=%.8g_theta=%.8g_rho0=%.8g_mag0=%.8g_gamma=%.8g_LX=%d_LY=%d_init=%d_ran=%d_*.png'%(beta1,beta2,D,v,theta,rho0,mag0,gamma,LX,LY,init,ran))

print('OK - time=%d sec'%(time.time()-clock))
