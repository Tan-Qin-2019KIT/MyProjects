%% this code is used for the exercise 1 of seismic modelling
%   "Modelling of seismic wavesby using Green¡¯s function"
%   it is made by Tan Qin in GPI
clear all;close all;

%% define the basic parameters
v=1500;
xs=1250; ys=xs; zs=xs;
xr=100:100:2500; yr=1250; zr=1250;
fc=50; t0=600e-3; dt=1e-3; tmax=1.2;
dx=50; x=0:dx:2500; y=x; z=x; t=0:dt:tmax;
[x_2D y_2D]=meshgrid(x,y);
[x_3D y_3D z_3D]=meshgrid(x,y,z);
r_1D=abs(x-xs); % receiver-source offset
r_2D=sqrt((x_2D-xs).^2+(y_2D-ys).^2);
r_3D=sqrt((x_3D-xs).^2+(y_3D-ys).^2+(z_3D-zs).^2);

%% symfunction which might be used in  the calculation of analytical solution
% syms x_sym t_sym f_sym
% func_r_1D(x_sym)=abs(x_sym-xs);
% func_r_2D(x_sym)=sqrt((x_sym-xs).^2+(y_2D-ys).^2);
% func_r_3D(x_sym)=sqrt((x_sym-xs).^2+(y_3D-ys).^2+(z_3D-zs).^2);
% func_s(t_sym)=(1-2*(pi*fc*(t_sym-t0)).^2).*exp(-(pi*fc*(t_sym-t0)).^2);
% func_S(f_sym)=fourier(func_s(t_sym),t_sym,f_sym);
% func_Green_1D(x_sym,t_sym)=v/2*heaviside(t_sym-func_r_1D(x_sym)/v);
% func_GREEN_1D(x_sym,f_sym)=fourier(func_Green_1D(x_sym,t_sym),t_sym,f_sym);
% func_p_1D(x_sym,t_sym)=ifourier(func_GREEN_1D(x_sym,f_sym) * func_S(f_sym),f_sym,t_sym);
% func_Green_2D(x_sym,t_sym)=1/(2*pi)*heaviside(t_sym-func_r_2D/v)./sqrt(t_sym.^2-(func_r_2D/v).^2);
% func_GREEN_2D(x_sym,f_sym)=fourier(func_Green_2D(x_sym,t_sym),t_sym,f_sym);
% func_p_2D(x_sym,t_sym)=ifourier(func_GREEN_2D(x_sym,f_sym) * func_S(f_sym),f_sym,t_sym);
% func_Green_3D(x_sym,t_sym)=1./(4*pi*func_r_3D).*dirac(t_sym-func_r_3D(x_sym)/v);
% func_GREEN_3D(x_sym,f_sym)=fourier(func_Green_3D(x_sym,t_sym),t_sym,f_sym);
% func_p_3D(x_sym,t_sym)=ifourier(func_GREEN_3D(x_sym,f_sym) * func_S(f_sym),f_sym,t_sym);

%% define a ricker wavelet and plot it
s=(1-2*(pi*fc*(t-t0)).^2).*exp(-(pi*fc*(t-t0)).^2);
% s=double(subs(func_s,t_sym,t));
figure(1)
plot(t./dt,s,'.-');
set(gca,'YLim',[-2 2]);
xlabel('t in s')
ylabel('Amplitude')
title('Ricker wavelet')
pause(0.001)

%% update the waveform field by time
t_step=1;
for i=1:t_step:length(t)
    % Green function calculation, pleasse note that we should check the
    % data wheather there are infinite or NAN value.
    Green_1D=v/2*heaviside(t(i)-r_1D/v);
    Green_2D=1/(2*pi)*heaviside(t(i)-r_2D/v)./sqrt(t(i).^2-(r_2D/v).^2);
    delta_t=t(i)-r_3D/v;
    [index_zero]=find(abs(delta_t)<=50*dt);
    delta_t(index_zero)=0;
    Green_3D=1./(4*pi*r_3D).*dirac(delta_t);
    %     Green_2D=double(subs(func_Green_2D,{x_sym,t_sym},{x_2D,t(i)}));
    %     Green_1D=double(subs(func_Green_1D,{x_sym,t_sym},{x,t(i)}));
    %     Green_3D=double(subs(func_Green_3D,{x_sym,t_sym},{x_3D,t(i)}));
    [index_inf]=isinf(Green_1D);
    Green_1D(index_inf)=1;
    [index_nan]=isnan(Green_1D);
    Green_1D(index_nan)=0;
    [index_inf]=isinf(Green_2D);
    Green_2D(index_inf)=1;
    [index_nan]=isnan(Green_2D);
    Green_2D(index_nan)=0;
    [index_inf]=isinf(Green_3D);
    Green_3D(index_inf)=1;
    [index_nan]=isnan(Green_3D);
    Green_3D(index_nan)=0;
    
    % pressure field is the convolution of Green function and the source
    p_1D=conv(Green_1D,s,'same');
    p_2D=conv2(Green_2D,s,'same');
    p_3D=convn(Green_3D,s,'same');
    %     p_1D=double(subs(func_p_1D,{x_sym,t_sym},{x,t(i)}));
    %     p_2D=double(subs(func_p_2D,{x_sym,t_sym},{x_2D,t(i)}));
    %     p_3D=double(subs(func_p_3D,{x_sym,t_sym},{x_3D,t(i)}));
    [index_inf]=isinf(p_1D);
    p_1D(index_inf)=1;
    [index_nan]=isnan(p_1D);
    p_1D(index_nan)=0;
    [index_inf]=isinf(p_2D);
    p_2D(index_inf)=1;
    [index_nan]=isnan(p_2D);
    p_2D(index_nan)=0;
    [index_inf]=isinf(p_3D);
    p_3D(index_inf)=1;
    [index_nan]=isnan(p_3D);
    p_3D(index_nan)=0;
    
    % record the signal in receiver's position
    t_trace((i-1)/t_step+1)=t(i);
    trace_1D((i-1)/t_step+1,1:length(xr))=p_1D(xr./dx);
    trace_2D((i-1)/t_step+1,1:length(xr))=p_2D(xr./dx,yr/dx);
    trace_3D((i-1)/t_step+1,1:length(xr))=p_3D(xr./dx,yr/dx,zr/dx);
    
    % imaging region
    figure(2)
    subplot(2,3,1)
    plot(x, Green_1D,'.')
    set(gca,'XLim',[0 2500]);
    xlabel('x in m')
    ylabel('Green function')
    title('Green function 1D')
    caxis([-max(abs(Green_1D)),max(abs(Green_1D))]);
    
    subplot(2,3,4)
    plot(x,p_1D./max(abs(p_1D)),'.-');
    set(gca,'XLim',[0 2500]);
    set(gca,'YLim',[-2 2]);
    xlabel('x in m')
    ylabel('Amplitude')
    title('P-wave 1D')
    
    subplot(2,3,2)
    imagesc(x,y,Green_2D);
    axis([0 2500 0 2500]);
    xlabel('x in m')
    ylabel('y in m')
    title('Green function 2D')
    axis equal
    colormap(jet);colorbar;
    caxis([-max(max(abs(Green_2D))),max(max(abs(Green_2D)))]);
     
    subplot(2,3,5)
    imagesc(x,y,p_2D./max(max(abs(p_2D))));
    axis([0 2500 0 2500]);
    xlabel('x in m')
    ylabel('y in m')
    title('P-wave 2D')
    axis equal
    colormap(jet);colorbar;caxis([-2,2]);
    
    subplot(2,3,3)
    slice(x_3D,y_3D,z_3D,Green_3D,1250,1250,1250);
    shading interp
    xlabel('x in m')
    ylabel('y in m')
    zlabel('z in m')
    title('Green function 3D')
    colormap(jet);colorbar;
    caxis([-max(max(max(abs(Green_3D)))),max(max(max(abs(Green_3D))))]);
    
    subplot(2,3,6)
    slice(x_3D,y_3D,z_3D,p_3D./max(max(max(abs(p_3D)))),1250,1250,1250);
    shading interp
    xlabel('x in m')
    ylabel('y in m')
    zlabel('z in m')
    title('P-wave 3D')
    colormap(jet);colorbar;caxis([-2,2]);
    pause(0.001)
end

figure(3)
subplot(2,3,1)
for(i=1:length(xr))
    % pick the arrival time and calculate the P-wave velocity
    [trace_1D_max index_max]=max(abs(trace_1D(:,i)));
    v_pick_1D(i)=abs(xr(i)-xs)./(t_trace(index_max)-0);
    plot((trace_1D(:,i)./2./max(abs(trace_1D(:,i)))+i)*(xr(2)-xr(1)),t_trace,'k'); hold on
    plot((trace_1D(index_max,i)./2./max(abs(trace_1D(:,i)))+i)*(xr(2)-xr(1)),t_trace(index_max),'ro');
    hold on
end
set(gca,'XLim',[0 2500+(xr(2)-xr(1))]);
set(gca,'ydir','reverse');
xlabel('x in m')
ylabel('t in s')
title('recorded signals 1D')
subplot(2,3,4)
plot(xr,v_pick_1D,'ko');
set(gca,'XLim',[0 2500+(xr(2)-xr(1))]);
set(gca,'ydir','normal');
xlabel('x in m')
ylabel('v in m/s')
title('velocity calculated by 1D')

subplot(2,3,2)
for(i=1:length(xr))
    [trace_2D_max index_max]=max(abs(trace_2D(:,i)));
    v_pick_2D(i)=sqrt((xr(i)-xs).^2+(yr-ys).^2)./(t_trace(index_max)-0);
    plot((trace_2D(:,i)./2./max(abs(trace_2D(:,i)))+i)*(xr(2)-xr(1)),t_trace,'k'); hold on
    plot((trace_2D(index_max,i)./2./max(abs(trace_2D(:,i)))+i)*(xr(2)-xr(1)),t_trace(index_max),'ro');
    hold on
end
set(gca,'XLim',[0 2500]);
set(gca,'ydir','reverse');
xlabel('x in m')
ylabel('t in s')
title('recorded signals 2D')
subplot(2,3,5)
plot(xr,v_pick_2D,'ko');
set(gca,'XLim',[0 2500+(xr(2)-xr(1))]);
set(gca,'ydir','normal');
xlabel('x in m')
ylabel('v in m/s')
title('velocity calculated by 2D')

subplot(2,3,3)
for(i=1:length(xr))
    [trace_3D_max index_max]=max(abs(trace_3D(:,i)));
    v_pick_3D(i)=sqrt((xr(i)-xs).^2+(yr-ys).^2+(zr-zs).^2)./(t_trace(index_max)-0);
    plot((trace_3D(:,i)./2./max(abs(trace_3D(:,i)))+i)*(xr(2)-xr(1)),t_trace,'k'); hold on
    plot((trace_3D(index_max,i)./2./max(abs(trace_3D(:,i)))+i)*(xr(2)-xr(1)),t_trace(index_max),'ro');
    hold on
end
set(gca,'XLim',[0 2500]);
set(gca,'ydir','reverse');
xlabel('x in m')
ylabel('t in s')
title('recorded signals 3D')
subplot(2,3,6)
plot(xr,v_pick_3D,'ko');
set(gca,'XLim',[0 2500+(xr(2)-xr(1))]);
set(gca,'ydir','normal');
xlabel('x in m')
ylabel('v in m/s')
title('velocity calculated by 3D')