function resid = Q_t_Theis_obj(ST)
%S=exp(ST(1,1))  
S=ST(1,1)  
T=ST(1,2)
%C=exp(ST(1,3))
C=ST(1,3)
p=ST(1,4)  
r=0.33;   % radius from pumping well (of r of pumping well casing for step test)
makeplots=true;  % makeplots=false;  % False is faster but cannot see progress.

Green =@(r,S,T,tvec) exp(-(r*r*S/4/T./tvec))./tvec;

welldata = importdata('ET1_dd.txt'); % Assumes columns of (t,s)
Q = importdata('ET1_Q.txt');         % Assumes columns of (tstart_i,Q_i)

%Q(:,1)=Q(:,1)/24/60;                 % convert minutes to days if needed
%welldata(:,1)=welldata(:,1)/24/60;   % convert minutes to days if needed

tpeak=r*r*S/4/T;  %Note that small r^2*S and/or large T makes the convolution slow!
tmin=0; 
tmax=max(welldata(:,1));
dt=tpeak/4; % Divisor is points between zero and peak of Green's fn. Check fig 3 for a decent peak!
ntpoints = 1+ceil((tmax-tmin)/dt);
tvec=linspace(tmin,tmax,ntpoints);
dt=tvec(2)-tvec(1);

time=welldata(:,1);
obs1=welldata(:,2);

% Construct the Q vector, at same times as Green's, i.e., use same tvec:
Qvec=zeros(size(tvec)); 
for k=1:size(Q,1)
    tstart=Q(k,1); 
    Qnow=Q(k,2);
    Qvec(tvec>=tstart) = Qnow;
end
%Qvec=interp1(Q(:,1),Q(:,2),tvec);  % Or a linear interpolation of Q's.

Green1=Green(r,S,T,tvec); Green1(1)=0;

ltot=length(Green1)+length(Qvec)-1
s1=ifft(fft([Green1 zeros(1,ltot-length(Green1))]).*fft([Qvec zeros(1,ltot-length(Qvec))]));
s1=(dt/(4*pi*T))*s1(1:ntpoints);  % This is linear Theis response
%s1=dt*conv(Qvec,Green1,'full')/(4*pi*T); s1=s1(1:ntpoints);  % Slower method
s1=s1+C*Qvec.^p;                  % Add in the nonlinear well losses
s1atdata=interp1(tvec,s1,time);   % Interplolate the model to the observation times
Qatdata=interp1(tvec,Qvec,time);  % For use in weights (if you want 1/Q weights)

weights=ones(size(obs1));    % Just use this for no weighting
weights=abs(1./obs1);             % Assumes VAR(measurements) propto measurements 
%weights=abs(1./Qatdata);        % Assumes VAR(measurements) propto Q
%  Get rid of bad values (from zero Q or s)
weights(isinf(weights))=-9999; weights(weights==-9999)=max(weights);

resid=weights.*(s1atdata-obs1);
classicalRMSE= sqrt(mean((s1atdata-obs1).^2))  % Display for a quick check

%  Make some plots to see progress of parameter estimation (Really slow!!)
% These axis limits will give full log base 10 limits
if(makeplots)
minplottime=10^(floor(log10(min(time))));
maxplottime=10^(ceil(log10(max(time))));
minplots=1;
%minplots=10^(floor(log10(min(obs1))));
maxplots=10^(ceil(log10(max(obs1))));

s1(s1<0)=1e-20;

figure(1)
semilogy(time,obs1,'bo',tvec,s1,'b-')
axis([0 1.1*tmax minplots maxplots])
xlabel('Time (d)'); ylabel('s (m)')
legend('data','model','Location','NE')

figure(2)
loglog(time,obs1,'bo',tvec,s1,'b-')
axis([minplottime maxplottime minplots maxplots])
xlabel('Time (d)'); ylabel('s (m)')
legend('data','model','Location','NE')

figure(3)
semilogx(time,obs1,'bo',tvec,s1,'b-')
axis([minplottime maxplottime min(obs1) 1.1*max(obs1)])
xlabel('Time (d)'); ylabel('s (m)')
legend('data','model','Location','NE')

figure(4)
loglog(tvec(2:end),Green1(2:end)/(4*pi*T),'-o',tvec,Qvec,'-d')
axis([tmin tmax 1e-6 1e6])
xlabel('Time (d)'); ylabel('d/m^2')
legend('Greens function','Q(t)');
title('Final plot should resolve peak nicely')
end
