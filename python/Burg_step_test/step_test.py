#!/usr/bin/env python3
import numpy as np
import time
#######  First define the functions to calculate 
def resid(params,obs_dd,Q_data,r,obs_time):
	ax.cla()
	ax.plot(obs_time,obs_dd,'o',markersize=5, fillstyle = 'none')
	plt.legend(["model", "data"], loc ="upper right")
	plt.pause(0.01)
	weights=np.ones_like(obs_dd)
#	weights=np.abs(1.0/obs_dd)
#	weights[np.isinf(weights)]=-999
#	weights[weights==-999]=max(weights)
	resid=weights*(model_dd - obs_dd)
def makeQ_of_t(Q_data,tvec):
def make_tvec(S,T,r,tmax):
# This may need to be adjusted based on Green's function:

#########  Main Program - user inputs required in this section ###############

#  Read in observed drawdown in columns of time and drawdown
obs_data= np.loadtxt("Burg_dd.txt")
# Read in the pumping rate data (columns of start time for Q_i, then Q_i)
# Define starting value of parameters, and lower and upper bounds of vector [S,T,C,p]

##############################################################################
# Get a figure going for visualizing the evolving model fit ...
fig, ax=plt.subplots(figsize=(4,5))

# Do the nonlinear least squares 
final_params=res_lsq.x
#  Get the final modeled drawdown and plot against data
print('Execution time (h:m:s):', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
tplotmin=0.9*min(obs_time)
tplotmax=1.1*max(obs_time)
splotmin=0.9*min(obs_dd)
splotmax=1.1*max(obs_dd)

fig2, ax2 = plt.subplots(figsize=(5,5))

ax2.loglog(tvec,final_dd,'r-',label="model")
ax2.loglog(obs_time,obs_dd,'o',markersize=5, fillstyle = 'none', label="data")
ax2.set_xlim(tplotmin,tplotmax)
ax2.set_ylim(splotmin,splotmax)
ax2.set_xlabel("Time (d)")
ax2.set_ylabel("Drawdown (ft)")
ax2.legend(loc ="upper left")

fig3, ax3 = plt.subplots(figsize=(5,5))
ax3.plot(tvec,final_dd,'r-',label="model")
ax3.plot(obs_time,obs_dd,'o',markersize=5, fillstyle = 'none', label="data")
ax3.set_xlim(tplotmin,tplotmax)
ax3.set_ylim(splotmin,splotmax)
ax3.set_xlabel("Time (d)")
ax3.set_ylabel("Drawdown (ft)")
ax3.legend()

plt.show()