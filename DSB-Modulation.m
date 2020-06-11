%========================================================================%
%                           IN THE NAME OF GOD                           %
%                           PROJECT OF  MATLAB                           %
%                         BY: MOHAMMAD JAVAD ADEL                        %
%                               9621010042                               %
%                             DATE: 98/10/14                             %
%========================================================================%
clc      % Clear Command Window
clear      % Remove items from workspace, freeing up system memory
close all  % closes all figures
%========================================================================%
%difine time signal x(t):
Fq=10000;                       %sampling frequency
T=0.01;
t=-T:1/Fq:2*T ;                 % time variable
X_t=0.5*(1-cos((2*pi*t)/T)) ;
Pulse=rectpuls(t-T/2,T);
X_t=X_t.*Pulse;
%========================================================================%
%plot time signal x(t):
figure;
subplot(3,3,[1,2,3]);
plot(t,X_t);
grid on ;
xlim([-T,2*T]);
ylim([-0.1,1]);
xlabel('time (second)')
ylabel('x(t)')
title('X(t) Graph')
%========================================================================%
%fourier transform of x(t)=X(f):
n=2^(floor(log2(length(X_t))+5)); 
X_f=fft(X_t,n);
f=-Fq/2:Fq/n:+Fq/2-Fq/n;
X_f=((1/Fq)*fftshift(X_f)).*exp(-1j*2*pi*(-T)*f);
%========================================================================%
%plot Domain and phase of X(f):
subplot(3,3,[4,5,6]);
plot(f,abs(X_f),'r');
grid on;
xlim([-Fq/2,Fq/2]);
ylim([-T/2,T]);
xlabel('f(Hz)');
ylabel('|X(f)|');
hold on;
subplot(3,3,[7,8,9]);
plot(f,angle(X_f),'r');
grid on;
xlim([-1000,1000]);
ylim([-2*pi,2*pi]);
xlabel('f(Hz)');
ylabel('<X(f)');
hold on;
%========================================================================%
%plot Equivalent of X(f):
X_ff=(0.5*T*(sinc(f*T)./(1-(f.^2)*(T.^2)))).*exp(-1j*pi*f*T);
subplot(3,3,[4,5,6]);
plot(f,abs(X_ff),'b');
grid on;
xlim([-Fq/2,Fq/2]);
ylim([-T/2,T]);
xlabel('f(Hz)');
ylabel('|X(f)|');
subplot(3,3,[7,8,9]);
plot(f,angle(X_ff),'b');
grid on;
xlim([-1000,1000]);
ylim([-2*pi,2*pi]);
xlabel('f(Hz)');
ylabel('<X(f)');
%========================================================================%
%difine Xc(t) for DSB Modolation:
fc=1000;
Ac=10;
Xc_t=Ac*X_t.*cos(2*pi*fc*t);
figure;
subplot(3,3,[1,2,3]);
plot(t,Xc_t);
grid on ;
xlim([-T,2*T]);
ylim([-10,10]);
xlabel('time (second)');
ylabel('x(t)');
title('X(t) Graph');
%========================================================================%
%fourier transform of DCB Modolation:
n=2^(floor(log2(length(Xc_t))+5));
f=-Fq/2:Fq/n:+Fq/2-Fq/n;
Xc_f=fft(Xc_t,n);
Xc_f=(1/Fq)*fftshift(Xc_f).*exp(-1j*2*pi*(-T)*f);
%========================================================================%
%plot Xc(f):
subplot(3,3,[4,5,6]);
plot(f,abs(Xc_f),'r');
grid on;
xlim([-Fq/2,Fq/2]);
ylim([-0.03,0.03]);
xlabel('f(Hz)');
ylabel('|Xc(f)|');
subplot(3,3,[7,8,9]);
plot(f,angle(Xc_f),'r');
grid on;
xlim([-1000,1000]);
ylim([-2*pi,2*pi]);
xlabel('f(Hz)');
ylabel('Xc(f)');
%========================================================================%
%                                 END                                    %