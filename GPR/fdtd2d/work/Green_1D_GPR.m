clear all; close all;
epsilon_0=8.8541878176e-12; mu_0=4*pi*1e-7;
epsilon_r=9;sigma=0.01;
x2=1;

fc=450e6; number=2.^9; f_zoom=1;
t0=5e-9; dt=0.05e-9;

[ t,w_t ] = func_ricker( fc,number,t0,dt );

df=1/dt;
f=(0:number-1)*df/number;
f_interp=f/f_zoom;
W_f=fft(w_t,number);  % calculate the frequency spectrum
W_f_interp = interp1(f,W_f,f_interp,'spline');

figure(1)
subplot(221);
plot(t*1e9,w_t,'.-');grid on;
title('Ricker-wavelet');
xlabel('t£¨ns£©');
ylabel('A-t');
subplot(222)
plot(f_interp*1e-6,abs(W_f_interp),'.-');grid on;
title('Ricker-frequency');
xlabel('f£¨MHz£©');
ylabel('A-f');

omega=2*pi*f;
Z_0=sqrt(mu_0/epsilon_0); Z_1=Z_0/sqrt(epsilon_r);
c_0=1/sqrt(epsilon_0*mu_0); c_1=c_0/sqrt(epsilon_r);
G_f_0=Z_0./(2*pi*(1-epsilon_r)*x2.^2)*exp(-j*omega*abs(x2)/c_0);
G_f_1=-sqrt(epsilon_r)*Z_0./(2*pi*(1-epsilon_r)*x2.^2)*...
    exp(-j*omega*abs(x2)/c_1-sigma*Z_1./2*abs(x2));
E_f_0=G_f_0.*W_f; E_f_1=G_f_1.*W_f;
E_f=E_f_0+E_f_1;
e_t=ifft(E_f,number);
e_t_0=ifft(E_f_0,number);
e_t_1=ifft(E_f_1,number);
E_f_interp = interp1(f,E_f,f_interp,'spline');
E_f_0_interp = interp1(f,E_f_0,f_interp,'spline');
E_f_1_interp = interp1(f,E_f_1,f_interp,'spline');

w_t_c0=[zeros(1,floor(abs(x2)/c_0/dt)) w_t(1:end-floor(abs(x2)/c_0/dt))];
w_t_c1=[zeros(1,floor(abs(x2)/c_1/dt)) w_t(1:end-floor(abs(x2)/c_1/dt))];
e_t_c0=Z_0./(2*pi*(1-epsilon_r)*x2.^2)*w_t_c0;
e_t_c1=-sqrt(epsilon_r)*Z_0./(2*pi*(1-epsilon_r)*x2.^2)*...
    exp(-sigma*Z_1./2*abs(x2))*w_t_c1;
e_t_c01=e_t_c0+e_t_c1;

subplot(223)
plot(t*1e9,real(e_t_0),'r*');hold on;
plot(t*1e9,real(e_t_1),'go');hold on;
plot(t*1e9,real(e_t),'.');hold on;
plot(t*1e9,real(e_t_c01),'-');grid on;
title('Synthetic Record');
xlabel('t£¨ns£©');
ylabel('A-t');
legend('air wave','ground wave','synthetic wave f','synthetic wave t');legend('boxoff');
subplot(224);
plot(f_interp*1e-6,abs(E_f_0_interp),'r*');hold on;
plot(f_interp*1e-6,abs(E_f_1_interp),'go');hold on;
plot(f_interp*1e-6,abs(E_f_interp),'.');grid on;
title('Synthetic Record-frequency');
xlabel('f£¨MHz£©');
ylabel('A-f');

%     v=1./(sqrt(mu_r*mu_0.*epsilon_r*epsilon_0/2.*(sqrt(1+(sigma./...
%         (epsilon_r.*epsilon_0*omega)).^2)+1)));
%     k=omega./v;
%     k12=sqrt(k(1)^2+k(2)^2);
%     zeta_0=j*omega*mu_0;
%     eta=sigma+j*omega*epsilon_r*epsilon_0;
%     Gamma=sqrt(k(1)^2+k(2)^2+eta*zeta_0);
%     rho=sqrt(x1^2+x2^2);
%     Phi=atan(x2/x1);
