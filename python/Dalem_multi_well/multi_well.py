#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
import time
#######  First define the functions to calculate 
	return model_dd
	s_model=np.empty_like(obs_time)
	max_obs_time=max(obs_time)
	tvec=make_tvec(params[0],params[1],min(r),max_obs_time)
	Q=makeQ_of_t(Q_data,tvec)
	beginidx=int(0)
		model_dd=fun_pump(params,Q,tvec,r[i],obs_time[beginidx:(beginidx+datalength[i])])
		s_model[beginidx:(beginidx+datalength[i])]=model_dd
		beginidx=beginidx+datalength[i]
		
	ax.cla()
	ax.plot(obs_time,obs_dd,'o',markersize=5, fillstyle = 'none')
	plt.legend(["model", "data"], loc ="upper left")
	plt.pause(0.1)
	weights=1.0/obs_dd
	resid=weights*(s_model - obs_dd)
def makeQ_of_t(Q_data,tvec):
	return Q
def make_tvec(S,T,r,tmax):
        tvec[0]=1e-20

        return tvec	
#########  Main Program - user inputs required in this section ###############

#  Read in observed drawdown in columns of time and drawdown
r=[30,60,90,120]    # Radii for all wells to be read in
nwells=len(r)
datalength=np.empty(nwells,dtype=int)

#  Read in observed drawdown in columns of time and drawdown and consolidate
c0=np.loadtxt("Dalem_r30_dd.txt")
datalength[0]=len(c0[:,0])
obs_data=c0
c1=np.loadtxt("Dalem_r60_dd.txt")
datalength[1]=len(c1[:,0])
obs_data=np.append(obs_data,c1,axis=0)
c2=np.loadtxt("Dalem_r90_dd.txt")
c3=np.loadtxt("Dalem_r120_dd.txt")
datalength[3]=len(c3[:,0])
obs_data=np.append(obs_data,c3,axis=0)

obs_dd=obs_data[:,1]
print(datalength)
# Define starting value of parameters, and lower and upper bounds
# For Walton's Leaky well fn, this is [S,T,c]

# Get a figure going for visualizing the evolving model fit ...
fig, ax=plt.subplots(figsize=(4,5))

# Do the nonlinear least squares 
final_params=res_lsq.x
#  Get the final modeled drawdown and plot against data
#  Make some half-way decent plots.  
tplotmin=0.9*min(obs_time)
tplotmax=1.1*max(obs_time)
splotmin=0.9*min(obs_dd)
splotmax=1.1*max(obs_dd)

fig2, (ax1, ax2) = plt.subplots(2)

beginidx=int(0)
	model_dd=fun_pump(final_params,Q,tvec,r[i],obs_time[beginidx:(beginidx+datalength[i])])
	ax1.plot(obs_time[beginidx:(beginidx+datalength[i])],model_dd,'r-',label="model")

ax1.plot(obs_time,obs_dd,'o',markersize=5, fillstyle = 'none', label="data")
ax2.loglog(obs_time,obs_dd,'o',markersize=5, fillstyle = 'none', label="data")

ax1.set_xlim(tplotmin,tplotmax)
ax1.set_ylim(splotmin,splotmax)
ax1.set_xlabel("Time (d)")
ax1.set_ylabel("Drawdown (ft)")
ax1.legend()
ax2.set_xlim(tplotmin,tplotmax)
ax2.set_ylim(splotmin,splotmax)
ax2.set_xlabel("Time (d)")
ax2.set_ylabel("Drawdown (ft)")
ax2.legend(loc ="upper left")
plt.show()