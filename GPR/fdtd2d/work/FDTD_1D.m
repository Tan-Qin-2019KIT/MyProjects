%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Scott Hudson, WSU Tri-Cities
%1D electromagnetic finite-difference time-domain (FDTD) program.
%Assumes Ey and Hz field components propagating in the x direction.
%Fields, permittivity, permeability, and conductivity
%are functions of x. Try changing the value of "profile".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;

L = 5; %domain length in meters
N = 505; %# spatial samples in domain
Niter = 2.^10; %# of iterations to perform
fs = 300e6; %source frequency in Hz
ds = L/N; %spatial step in meters
dt = ds/300e6; %"magic time step"
eps0 = 8.854e-12; %permittivity of free space
mu0 = pi*4e-7; %permeability of free space
x = linspace(0,L,N); %x coordinate of spatial samples
showWKB=0; %if =1 then show WKB appromination at end

%scale factors for E and H
ae = ones(N,1)*dt/(ds*eps0);
am = ones(N,1)*dt/(ds*mu0);
as = ones(N,1);
epsr = ones(N,1);
mur= ones(N,1);
sigma = zeros(N,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Here we specify the epsilon, sigma, and mu profiles. I've
%predefined some interesting examples. Try profile = 1,2,3,4,5,6 in sequence.
%You can define epsr(i), mur(i) (relative permittivity and permeability)
%and sigma(i) (conductivity) to be anything you want.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
profile = 1;
for i=1:N
   epsr(i) = 1;
   mur(i) = 1;
   w1 = 0.5;
   w2 = 1.5;
   if (profile==1) %dielectric window
      if (abs(x(i)-L/2)<0.5) epsr(i)=4; end
   end
   if (profile==2)%dielectric window with smooth transition
   	if (abs(x(i)-L/2)<1.5) epsr(i)=1+3*(1+cos(pi*(abs(x(i)-L/2)-w1)/(w2-w1)))/2; end
      if (abs(x(i)-L/2)<0.5) epsr(i)=4; end
   end
   if (profile==3) %dielectric discontinuity
      if (x(i)>L/2) epsr(i) = 9; end
   end
  	if (profile==4) %dielectric disontinuity with 1/4-wave matching layer
      if (x(i)>(L/2-0.1443)) epsr(i) = 3; end 
      if (x(i)>L/2) epsr(i) = 9; end
   end
   if (profile==5) %conducting half space
      if (x(i)>L/2) sigma(i) = 0.005; end
   end
   if (profile==6) %sinusoidal dielectric
      epsr(i) = 1+sin(2*pi*x(i)/L)^2;
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ae = ae./epsr;
am = am./mur;
ae = ae./(1+dt*(sigma./epsr)/(2*eps0));
as = (1-dt*(sigma./epsr)/(2*eps0))./(1+dt*(sigma./epsr)/(2*eps0));

%plot the permittivity, permeability, and conductivity profiles
figure(1)
subplot(3,1,1);
plot(x,epsr);
grid on;
axis([3*ds L min(epsr)*0.9 max(epsr)*1.1]);
title('relative permittivity');
subplot(3,1,2);
plot(x,mur);
grid on;
axis([3*ds L min(mur)*0.9 max(mur)*1.1]);
title('relative permeabiliity');
subplot(3,1,3);
plot(x,sigma);
grid on;
axis([3*ds L min(sigma)*0.9-0.001 max(sigma)*1.1+0.001]);
title('conductivity');

%initialize fields to zero
Hz = zeros(N,1);
Ey = zeros(N,1);
figure(2);
set(gcf,'doublebuffer','on'); %set double buffering on for smoother graphics
plot(Ey);
grid on;

t=(0:Niter)*dt;
% w=blackharrispulse(fs,t);
[t w]=func_ricker(fs,Niter,5e-9,dt);
for iter=1:Niter
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %The next 10 or so lines of code are where we actually integrate Maxwell's
   %equations. All the rest of the program is basically bookkeeping and plotting.
   %"smooth turn on" sinusoidal source 
%    Ey(3) = Ey(3)+2*(1-exp(-((iter-1)/50)^2))*sin(2*pi*fs*dt*iter); 
   Ey(3) = Ey(3)+w(iter);
   Hz(1) = Hz(2); %absorbing boundary conditions for left-propagating waves
   for i=2:N-1 %update H field
      Hz(i) = Hz(i)-am(i)*(Ey(i+1)-Ey(i));
   end
   Ey(N) = Ey(N-1); %absorbing boundary conditions for right propagating waves
   for i=2:N-1 %update E field
      Ey(i) = as(i)*Ey(i)-ae(i)*(Hz(i)-Hz(i-1));
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if mod(iter,5)==0
       figure(2)
       hold off
       plot(x,Ey,'b');
       axis([3*ds L -2 2]);
       grid on;
       title('E (blue) and 377*H (red)');
       hold on
       plot(x,377*Hz,'r');
       xlabel('x (m)');
       pause(0.01);
   end
end

%WKB prediction for the fields
if (showWKB==1)
	phase = cumsum((epsr).^0.5)*ds;
	beta0 = 2*pi*fs/(300e6);
   theory = sin(2*pi*fs*(Niter+4)*dt-beta0*phase)./(epsr.^0.25);
   input('press enter to show WKB theory');
   plot(x,theory,'b.');
   theory = sin(2*pi*fs*(Niter+4)*dt-beta0*phase).*(epsr.^0.25);
   plot(x,theory,'r.');
   title('E (blue), 377*H (red), WKB theory (points)');
end


