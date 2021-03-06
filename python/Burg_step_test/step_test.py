#!/usr/bin/env python3
import numpy as npimport matplotlib.pyplot as pltfrom scipy.optimize import least_squares
import time
#######  First define the functions to calculate def fun_pump(params,Q,tvec,r,obs_time):	S=np.exp(params[0])	print('S = ' +str(S))	T=params[1]	print('T = ' +str(T))	C=np.exp(params[2])	print('C = ' +str(C))	p=params[3]	print('p = ' +str(p))	dt=tvec[2]-tvec[1]	green=(np.exp(-(r*r*S/4./T/tvec)))/tvec	green[0]=0.0	lwant=len(tvec)	ltot=len(green)+len(Q)	print('number of convolution points = '+str(ltot)+'  working on it ...')	z=np.zeros(ltot-lwant)	s1=np.fft.irfft( np.fft.rfft(np.append(green,z))*np.fft.rfft(np.append(Q,z)) )	#s1=np.convolve(Q,green)	s1=dt*s1[0:len(tvec)]	model_dd=s1/(4*np.pi*T)+(C*Q**p)	model_dd=np.interp(obs_time,tvec,model_dd)	return model_dd
def resid(params,obs_dd,Q_data,r,obs_time):	max_obs_time=max(obs_time)	tvec=make_tvec(np.exp(params[0]),params[1],r,max_obs_time)	Q=makeQ_of_t(Q_data,tvec)	model_dd=fun_pump(params,Q,tvec,r,obs_time)
	ax.cla()	ax.plot(obs_time,model_dd)
	ax.plot(obs_time,obs_dd,'o',markersize=5, fillstyle = 'none')
	plt.legend(["model", "data"], loc ="upper right")
	plt.pause(0.01)
	weights=np.ones_like(obs_dd)
#	weights=np.abs(1.0/obs_dd)
#	weights[np.isinf(weights)]=-999
#	weights[weights==-999]=max(weights)
	resid=weights*(model_dd - obs_dd)	print('Classical RMSE = '+str(np.sqrt(np.mean((model_dd - obs_dd)**2)))+'\n')		return resid
def makeQ_of_t(Q_data,tvec):	Q=np.zeros_like(tvec)	for k in range(0,len(Q_data[:,0])):		Qnow=Q_data[k,1]		tnow=Q_data[k,0]		Q[tvec>tnow]=Qnow	return Q
def make_tvec(S,T,r,tmax):        extrat=1.2        tmin=0.0 
# This may need to be adjusted based on Green's function:        tpeak=r*r*S/4.0/T        dt=tpeak/4.0 # Divisor is points between zero and peak of Green's fn. Check fig 3 for a decent peak!        ntpoints = 1+int(np.ceil((extrat*tmax-tmin)/dt))        tvec=np.linspace(tmin,extrat*tmax,ntpoints)        tvec[0]=1e-20        return tvec	
start_time = time.process_time()
#########  Main Program - user inputs required in this section ###############

#  Read in observed drawdown in columns of time and drawdown
obs_data= np.loadtxt("Burg_dd.txt")print(obs_data)obs_time=obs_data[:,0]obs_dd=obs_data[:,1]r=0.1					# Well radius
# Read in the pumping rate data (columns of start time for Q_i, then Q_i)Q_data = np.loadtxt("Burg_Q.txt")print(Q_data[:,1])
# Define starting value of parameters, and lower and upper bounds of vector [S,T,C,p]ST=np.array( [np.log(0.2), 4500, np.log(1e-7), 2.0])lower=np.array([np.log(1.e-3), 1.e-6, np.log(1.e-25), 1.])upper=np.array([np.log(0.5), 5000., np.log(1.), 5.]) 

##############################################################################
# Get a figure going for visualizing the evolving model fit ...
fig, ax=plt.subplots(figsize=(4,5))

# Do the nonlinear least squares #res_lsq=least_squares(resid,ST,tr_solver='lsmr',loss='soft_l1',args=(obs_dd,Q,tvec,r,obs_time))res_lsq=least_squares(resid,ST,bounds=(lower, upper),args=(obs_dd,Q_data,r,obs_time))
final_params=res_lsq.xprint('Final Parameters'+'\n')
#  Get the final modeled drawdown and plot against datatvec=make_tvec(np.exp(final_params[0]),final_params[1],r,max(obs_time))Q=makeQ_of_t(Q_data,tvec)final_dd=fun_pump(final_params,Q,tvec,r,tvec)elapsed_time = start_time - time.process_time()
print('Execution time (h:m:s):', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))#  Make some decent plots.  
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
