function resid = Q_t_Theis_obj_multi(ST)
%S=exp(ST(1,1))  
S=ST(1,1)  
T=ST(1,2)
L=sqrt(T*ST(1,3))

r=[30 60 90 120];   % radius from pumping well (of r of pumping well casing for step test)
%Green =@(r,S,T,tvec) exp(-(r*r*S/4/T./tvec))./tvec;
Green =@(r,S,T,L,tvec) exp(-(r*r*S/4/T./tvec)-(T*tvec/L/L/S))./tvec;

welldata1 = importdata('Dalem_r30_dd.txt');  % Assumes columns of (t,s)
welldata2 = importdata('Dalem_r60_dd.txt');
welldata3 = importdata('Dalem_r90_dd.txt');
welldata4 = importdata('Dalem_r120_dd.txt');
Q = importdata('Dalem_Q.txt');         % Assumes columns of (tstart_i,Q_i)
%Q(:,1)=Q(:,1)/24/60;                 % convert minutes to days if needed
%welldata1(:,1)=welldata1(:,1)/24/60;   % convert minutes to days if needed
%welldata2(:,1)=welldata2(:,1)/24/60;   % convert minutes to days if needed
%welldata3(:,1)=welldata3(:,1)/24/60;   % convert minutes to days if needed

tpeak=min(r)*min(r)*S/4/T;  %Note that small r^2*S and/or large T makes the convolution slow!
tmin=0; 
tmax=max([max(welldata1(:,1))  max(welldata2(:,1))  max(welldata3(:,1)) max(welldata4(:,1))]);
dt=tpeak/6; % Divisor is points between zero and peak of Green's fn. Check fig 3 for a decent peak!
ntpoints = 1+ceil((tmax-tmin)/dt);
tvec=linspace(tmin,tmax,ntpoints);
dt=tvec(2)-tvec(1);

time1=welldata1(:,1); time2=welldata2(:,1); time3=welldata3(:,1); time4=welldata4(:,1);
obs1 =welldata1(:,2); obs2 =welldata2(:,2); obs3 =welldata3(:,2); obs4 =welldata4(:,2);
allobs=[obs1' obs2' obs3' obs4']';

% Construct the Q vector, at same times as Green's, i.e., use same tvec:
Qvec=zeros(size(tvec)); 
for k=1:size(Q,1)
    tstart=Q(k,1); 
    Qnow=Q(k,2);
    Qvec(tvec>tstart) = Qnow;
end
%Qvec=interp1(Q(:,1),Q(:,2),tvec);  % Or a linear interpolation of Q's.

Green1=Green(r(1),S,T,L,tvec); Green1(1)=0;
Green2=Green(r(2),S,T,L,tvec); Green2(1)=0;
Green3=Green(r(3),S,T,L,tvec); Green3(1)=0;
Green4=Green(r(4),S,T,L,tvec); Green4(1)=0;

ltot=length(Green1)+length(Qvec)-1   % total length of fft vectors
s1=ifft(fft([Green1 zeros(1,ltot-length(Green1))]).*fft([Qvec zeros(1,ltot-length(Qvec))]));
s2=ifft(fft([Green2 zeros(1,ltot-length(Green2))]).*fft([Qvec zeros(1,ltot-length(Qvec))]));
s3=ifft(fft([Green3 zeros(1,ltot-length(Green3))]).*fft([Qvec zeros(1,ltot-length(Qvec))]));
s4=ifft(fft([Green4 zeros(1,ltot-length(Green3))]).*fft([Qvec zeros(1,ltot-length(Qvec))]));

s1=(dt/(4*pi*T))*s1(1:ntpoints);  % This is linear Theis response
s2=(dt/(4*pi*T))*s2(1:ntpoints);  % This is linear Theis response
s3=(dt/(4*pi*T))*s3(1:ntpoints);  % This is linear Theis response
s4=(dt/(4*pi*T))*s4(1:ntpoints); 
%  Make log plots not throw an error:
s1(s1<=0)=1e-20; s2(s2<=0)=1e-20; s3(s3<=0)=1e-20; s4(s4<=0)=1e-20;

s1atdata=interp1(tvec,s1,time1);   % Interplolate the model to the observation times
s2atdata=interp1(tvec,s2,time2);   % Interplolate the model to the observation times
s3atdata=interp1(tvec,s3,time3);   % Interplolate the model to the observation times
s4atdata=interp1(tvec,s4,time4);   % Interplolate the model to the observation times
smodel=[s1atdata' s2atdata' s3atdata' s4atdata']';

%Qatdata=interp1(tvec,Qvec,time);

weights=ones(size(smodel));    % Just use this for no weighting
%weights=1./allobs;             % Assumes VAR(measurements) propto measurements 
%weights=(1./Qatdata);        % Assumes VAR(measurements) propto Q
%  Get rid of bad values (from zero Q or s)
%weights(isinf(weights))=-9999; weights(weights==-9999)=max(weights);

resid=weights.*(smodel-allobs);
classicalRMSE= sqrt(mean((smodel-allobs).^2))

%  Make some plots to see progress of parameter estimation
% These axis limits will give full log base 10 limits
minobstime=min([min(time1) min(time2) min(time3)]) ;
minplottime=10^(floor(log10(minobstime)));
maxplottime=10^(ceil(log10(tmax)));
minplots=.01;
%minplots=10^(floor(log10(min(obs1))));
maxplots=10^(ceil(log10(max(obs1))));

figure(1)
semilogy(time1,obs1,'bo',tvec,s1,'b-');
hold on
plot(time2,obs2,'gd',tvec,s2,'g-');
plot(time3,obs3,'rsq',tvec,s3,'r-');
plot(time4,obs4,'ko',tvec,s4,'k-');
axis([0 1.1*tmax minplots maxplots])
xlabel('Time (d)'); ylabel('s (m)')
legend('data','model','Location','NW')
drawnow; hold off;

figure(2)
loglog(time1,obs1,'bo',tvec,s1,'b-');
hold on;
plot(time2,obs2,'gd',tvec,s2,'g-');
plot(time3,obs3,'rsq',tvec,s3,'r-');
plot(time4,obs4,'ko',tvec,s4,'k-');
axis([minplottime maxplottime minplots maxplots])
xlabel('Time (d)'); ylabel('s (m)')
legend('data','model','Location','NW')
drawnow; hold off

figure(3)
semilogx(time1,obs1,'bo',tvec,s1,'b-');
hold on;
plot(time2,obs2,'gd',tvec,s2,'g-');
plot(time3,obs3,'rsq',tvec,s3,'r-');
plot(time4,obs4,'ko',tvec,s4,'k-');

axis([minplottime maxplottime 0 1.1*max(obs1)])
xlabel('Time (d)'); ylabel('s (m)')
legend('data','model','Location','NW');
drawnow; hold off;

figure(4)
loglog(tvec,Green1/(4*pi*T),'-o',tvec,Qvec,'-d')
axis([tmin tmax 1e-6 1e6])
xlabel('Time (d)'); ylabel('d/m^2')
legend('Greens function','Q(t)');
title('Final plot should resolve peak nicely')




