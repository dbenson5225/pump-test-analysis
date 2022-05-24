clear all
% Initial guess of T and S:

% Some log-transformed parameters * Also check the function that calculates drawdown!
%ST=[log(0.001) 20 log(2e-10) 4.0]
%lb=[log(1e-4) .1 log(1e-20) 1];
%ub=[log(0.5) 1e6 log(1) 5]

% Linear parameters * Also check the function that calculates drawdown!
% In order: [S T C p]
%ST=[0.05 4000 1e-7 2.0]
%lb=[1e-5 2000 1e-20 1];
%ub=[0.5 500000 1 5]

%  For observation wells only (usinf Walton's Leaky Aquifer W(u,r/B)
% Remember that r/B is a function of other parameters except a new unique parmater "c":
% In order: [S T c]
ST=[0.05 400 500]
lb=[1e-5 1 1];
ub=[0.5 5000 5000]

% optimization options
%options = optimset('PlotFcns',@optimplotfval,'Algorithm','levenberg-marquardt');
options = optimset('PlotFcns',@optimplotfval);

[params,resnorm,resid,flag,output,lambda,jacobian] = lsqnonlin(@Q_t_Theis_obj_multi,ST,lb,ub,options)
%[params,R,J,CovB,MSE,errorinfo]=nlinfit()
% Just for grins, get the parameter estimate covariance matrix and the 
% correlation coeffs. for paramters
covp = inv(jacobian'*jacobian)*mean(resid.*resid)
sd   = sqrt(diag(covp))
corr = full(covp./(sd*sd'))