clear all
% Initial guess of T and S:
tic
% Some log-transformed parameters * Also check the function that calculates drawdown!
%ST=[log(0.2) 40 log(1e-7) 2]
%lb=[log(1e-5) .1 log(1e-20) 1];
%ub=[log(0.25) 1e5 log(1) 5]

% Linear parameters * Also check the function that calculates drawdown!
% In order: [S T C p]
ST=[0.01 20 2e-10 4];
lb=[1e-5 1 1e-20 1];
ub=[0.25 5000 1 5]
% optimization options
%options = optimset('PlotFcns',@optimplotfval,'Algorithm','levenberg-marquardt');
options = optimset('PlotFcns',@optimplotfval);

[params,resnorm,resid,flag,output,lambda,jacobian] = lsqnonlin(@Q_t_Theis_obj,ST,lb,ub,options)
%[params,R,J,CovB,MSE,errorinfo]=nlinfit()
% Just for grins, get the parameter estimate covariance matrix and the 
% correlation coeffs. for paramters
covp = inv(jacobian'*jacobian)*mean(resid.*resid)
sd   = sqrt(diag(covp))
corr = full(covp./(sd*sd'))
toc