%========================================================================%
%                           IN THE NAME OF GOD                           %
%                           PROJECT OF  MATLAB                           %
%                         BY: MOHAMMAD JAVAD ADEL                        %
%                               9621010042                               %
%                             DATE: 98/11/7                              %
%========================================================================%
clc        % Clear Command Window
clear      % Remove items from workspace, freeing up system memory
close all  % closes all figures
%========================================================================%
Fs=20000;               %frequency of sampling
Ts=1/Fs;                %time of sampling
T=0.05;                 
t=0:Ts:T;               %time variable
Am=1;                   %amplitude of time signal
Fm=100;                 %frequency of time signal
X_t=Am*sin(2*pi*Fm.*t); %time signal
figure;                 %opening a window for ploting
subplot(2,2,[1 2]);     %ploting in one figure
plot(t,X_t);            %plot time signal
grid on ;               
xlim([0,T]);
ylim([-Am,Am]);
ylabel('x(t)');
title('X(t) Graph');
%fourier transform of x(t)=> X(f) and ploting:
q=2^(floor(log2(length(X_t))+5));                
X_f=fft(X_t,q);
f=-Fs/2:Fs/q:+Fs/2-Fs/q;
X_f=(Ts*fftshift(X_f)).*exp(-1j*2*pi*(-T)*f);
subplot(2,2,[3 4]);
plot(f,abs(X_f));
grid on ;
xlim([-2*Fm,2*Fm]);
ylim([0,T/2]);
ylabel('X(f)');
title('X(f) Graph');
%========================================================================%
u=0.5;
Ac=10;
Fc1=2000;
figure;
Xc_t=Ac*(1+u*X_t).*cos(2*pi*Fc1.*t);
An=1;
n_t=An*rand(size(Xc_t));
Xc_t=n_t+Xc_t;
subplot(2,2,[1 2]);
plot(t,Xc_t);
grid on ;
xlim([0,T]);
ylim([-(Ac+u*Am*Ac),(Ac+u*Am*Ac)]);
xlabel('time (second)');
ylabel('Xc(t)');
title('Xc(t) Graph with noise');
%fourier transform of Xc(t)=> Xc(f) and ploting:
q=2^(floor(log2(length(Xc_t))+5));                
Xc_f=fft(Xc_t,q);
f=-Fs/2:Fs/q:+Fs/2-Fs/q;
Xc_f=(Ts*fftshift(Xc_f)).*exp(-1j*2*pi*(-T)*f);
subplot(2,2,[3 4]);
plot(f,abs(Xc_f));
grid on ;
xlim([-2*Fc1,2*Fc1]);
ylim([0,(u*Am*Ac)/20]);
ylabel('Xc(f)');
title('Xc(f) Graph');
%========================================================================%
Alo=4;
Fc2=Fc1;
Xlo_t=Alo*cos(2*pi*Fc2.*t);
Xlo_t=Xc_t.*Xlo_t;
figure;
subplot(2,2,[1 2]);
plot(t,Xlo_t);
grid on ;
xlim([0,T]);
ylim([0,(Ac+u*Am*Ac)*Alo]);
ylabel('Xlo(t)');
title('Xlo(t) Graph');
%fourier transform of Xlo(t)=> Xlo(f) and ploting:
q=2^(floor(log2(length(Xlo_t))+5));                
Xlo_f=fft(Xlo_t,q);
f=-Fs/2:Fs/q:+Fs/2-Fs/q;
Xlo_f=(Ts*fftshift(Xlo_f)).*exp(-1j*2*pi*(-0.05)*f);
subplot(2,2,[3 4]);
plot(f,abs(Xlo_f));
grid on ;
xlim([-Fs/2,Fs/2]);
ylim([0,1]);
ylabel('Xlo(t)');
title('Xlo(f) Graph');
%========================================================================%
fL=Fm;
n=110;
figure;
g0=0.8;
b = fir1(n,fL/(Fs/2),'low'); % Lowpass filter design
z_t=filter(b,1,Xlo_t); % filtering
A=mean(Xlo_t);
z_t=(z_t - A);
m=2/(u*Ac*Alo*g0);
z_t=(z_t)*m;
subplot(2,2,[1 2]);
plot(t,z_t);
grid on ;
xlim([0,T]);
ylim([-Am,Am]);
xlabel('time (second)');
ylabel('z(t)');
title('z(t) Graph');
%fourier transform of z(t)=> Z(f) and ploting:
q=2^(floor(log2(length(z_t))+5));                
z_f=fft(z_t,q);
f=-Fs/2:Fs/q:+Fs/2-Fs/q;
z_f=(Ts*fftshift(z_f)).*exp(-1j*2*pi*(-T)*f);
subplot(2,2,[3 4]);
plot(f,abs(z_f));
grid on ;
xlim([-2*Fm,2*Fm]);
ylim([0,0.03]);
ylabel('z(f)');
title('z(f) Graph');
%========================================================================%
A_t=abs(Xc_t);
figure;
subplot(2,2,[1 2]);
fL=300;
n=40;
b = fir1(n,fL/(Fs/2),'low'); % Lowpass filter design
z_t=filter(b,1,A_t); % filtering
A=mean(A_t);
z_t=(z_t - A)/(u*Ac*0.65);
%subplot(3,2,[3 4]);
plot(t,z_t);
grid on ;
xlim([0,T]);
ylim([-Am,Am]);
ylabel('Z(t)');
title('Z(t) Graph');
%fourier transform of Z(t)=> Z(f) and ploting:
q=2^(floor(log2(length(z_t))+5));                
z_f=fft(z_t,q);
f=-Fs/2:Fs/q:+Fs/2-Fs/q;
z_f=(Ts*fftshift(z_f)).*exp(-1j*2*pi*(-T)*f);
subplot(2,2,[3 4]);
plot(f,abs(z_f));
grid on ;
xlim([-2*Fm,2*Fm]);
ylim([0,T/2]);
ylabel('Z(f)');
title('Z(f) Graph');
%========================================================================%
