#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
import time
#######  First define the functions to calculate def fun_pump(params,Q,tvec,r,obs_time):	S=params[0]	print('S = ' +str(S))	T=params[1]	print('T = ' +str(T))	c=params[2]	print('c = ' +str(c))#	p=params[3]#	print('p = ' +str(p))	L=np.sqrt(T*c)	dt=tvec[2]-tvec[1]	green=(np.exp(-(r*r*S/4./T/tvec) - (T*tvec/L/L/S) ) )/tvec	green[0]=0.0	lwant=len(tvec)	ltot=len(green)+len(Q)	print('number of convolution points = '+str(ltot)+'  working on it ...')	z=np.zeros(ltot-lwant)	s1=np.fft.irfft( np.fft.rfft(np.append(green,z))*np.fft.rfft(np.append(Q,z)) )	#s1=np.convolve(Q,green)	s1=dt*s1[0:len(tvec)]	model_dd=s1/(4*np.pi*T)   #+(C*Q**p)	model_dd=np.interp(obs_time,tvec,model_dd)
	return model_dddef resid(params,obs_dd,Q_data,r,obs_time,datalength,nwells):
	s_model=np.empty_like(obs_time)
	max_obs_time=max(obs_time)
	tvec=make_tvec(params[0],params[1],min(r),max_obs_time)
	Q=makeQ_of_t(Q_data,tvec)
	beginidx=int(0)	for i in range(0,nwells):
		model_dd=fun_pump(params,Q,tvec,r[i],obs_time[beginidx:(beginidx+datalength[i])])
		s_model[beginidx:(beginidx+datalength[i])]=model_dd
		beginidx=beginidx+datalength[i]
		
	ax.cla()	ax.plot(obs_time,s_model,'.')
	ax.plot(obs_time,obs_dd,'o',markersize=5, fillstyle = 'none')
	plt.legend(["model", "data"], loc ="upper left")
	plt.pause(0.1)
	weights=1.0/obs_dd
	resid=weights*(s_model - obs_dd)	print('Classical RMSE = '+str(np.sqrt(np.mean((s_model - obs_dd)**2)))+'\n')		return resid
def makeQ_of_t(Q_data,tvec):	Q=np.zeros_like(tvec)	for k in range(0,len(Q_data[:,0])):		Qnow=Q_data[k,1]		tnow=Q_data[k,0]		Q[tvec>tnow]=Qnow
	return Q
def make_tvec(S,T,r,tmax):        extrat=1.2        tmin=0.0 # This may need to be adjusted based on Green's function        tpeak=r*r*S/4.0/T        dt=tpeak/4.0 # Divisor is points between zero and peak of Green's fn. Check fig 3 for a decent peak!        ntpoints = max(5000,1+int(np.ceil((extrat*tmax-tmin)/dt)))        tvec=np.linspace(tmin,extrat*tmax,ntpoints)
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
c2=np.loadtxt("Dalem_r90_dd.txt")datalength[2]=len(c2[:,0])obs_data=np.append(obs_data,c2,axis=0)
c3=np.loadtxt("Dalem_r120_dd.txt")
datalength[3]=len(c3[:,0])
obs_data=np.append(obs_data,c3,axis=0)

obs_dd=obs_data[:,1]obs_time=obs_data[:,0]print(obs_time)
print(datalength)# Read in the pumping rate data (columns of start time for Q_i, then Q_i)Q_data = np.loadtxt("Dalem_Q.txt")print(Q_data[:,1])
# Define starting value of parameters, and lower and upper bounds
# For Walton's Leaky well fn, this is [S,T,c]ST=np.array( [0.05, 400, 400])lower=np.array([1e-5, 1.e-6, 1e-6])upper=np.array([0.5, 10000., 5000.]) 
##############################################################################
# Get a figure going for visualizing the evolving model fit ...
fig, ax=plt.subplots(figsize=(4,5))

# Do the nonlinear least squares res_lsq=least_squares(resid,ST,bounds=(lower, upper),args=(obs_dd,Q_data,r,obs_time,datalength,nwells))
final_params=res_lsq.xprint('----- Final Parameters -----')
#  Get the final modeled drawdown and plot against datatvec=make_tvec(final_params[0],final_params[1],min(r),max(obs_time))Q=makeQ_of_t(Q_data,tvec)
#  Make some half-way decent plots.  
tplotmin=0.9*min(obs_time)
tplotmax=1.1*max(obs_time)
splotmin=0.9*min(obs_dd)
splotmax=1.1*max(obs_dd)

fig2, (ax1, ax2) = plt.subplots(2)

beginidx=int(0)for i in range(0,nwells):
	model_dd=fun_pump(final_params,Q,tvec,r[i],obs_time[beginidx:(beginidx+datalength[i])])
	ax1.plot(obs_time[beginidx:(beginidx+datalength[i])],model_dd,'r-',label="model")	ax2.loglog(obs_time[beginidx:(beginidx+datalength[i])],model_dd,'r-',label="model")	beginidx=beginidx+datalength[i]

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
