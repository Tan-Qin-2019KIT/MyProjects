%% this is the test code of exercise 1 - Visualization of Tomograms - 10pts
% it is compiled by Tan Qin in GPI KIT
clear all;close all;clc

%% to load the data first
load('layer6.dat');
X=layer6(:,1)';
Y=layer6(:,2)';
delta_vel=layer6(:,3)';

% to define the spatial sampling interval

delta_v_min=min(delta_vel);
delta_v_max=max(delta_vel);

% X_min=min(X);
% X_max=max(X);
% Y_min=min(Y);
% Y_max=max(Y);
% x=X_min-50:5:X_max+50;
% y=Y_min-50:5:Y_max+50;

x=-550:5:700;
y=-600:5:650;
[x_grid y_grid]=meshgrid(x,y);

%% interpolation of using two methods

% method 1: function griddata
% ZZ=griddata(X,Y,delta_vel,x_grid,y_grid,'linear');
% ZZ=griddata(X,Y,delta_vel,x_grid,y_grid,'nearest');
% ZZ=griddata(X,Y,delta_vel,x_grid,y_grid,'natural');
% ZZ=griddata(X,Y,delta_vel,x_grid,y_grid,'cubic');
% ZZ=griddata(X,Y,delta_vel,x_grid,y_grid,'v4');

% method 2: function spline2d_int
% t=0.01;
t=0.3;
% t=0.6;
% t=0.99;
delta_vel_grid = spline2d_int (x_grid, y_grid, X, Y, delta_vel, t );

%% image ragion
figure(1)
subplot(221)
mesh(x_grid,y_grid,delta_vel_grid); hold on
plot3(X,Y,delta_vel,'+');hold off
title('velocity perturbation')
xlabel('X in km')
ylabel('Y in km')
zlabel('Z in km/s')
set(gcs,'YDir','normal')

subplot(222)
imagesc(x,y,delta_vel_grid);hold on
plot(X,Y,'+');hold off
title('velocity perturbation (km/s)')
xlabel('X in km')
ylabel('Y in km')
set(gcs,'YDir','normal')
colormap(jet);colorbar;
caxis([0.5*delta_v_min 0.5*delta_v_max]);

subplot(223)
imagesc(x,y,delta_vel_grid);hold on
plot(X,Y,'+');hold off
title('velocity perturbation (km/s)')
xlabel('X in km')
ylabel('Y in km')
set(gcs,'YDir','normal')
colormap(jet);colorbar;
caxis([1.5*delta_v_min 1.5*delta_v_max]);

subplot(224)
imagesc(x,y,delta_vel_grid);hold on
plot(X,Y,'+');hold off
title('velocity perturbation (km/s)')
xlabel('X in km')
ylabel('Y in km')
set(gcs,'YDir','normal')
colormap(jet);colorbar;
caxis([3.0*delta_v_min 3.0*delta_v_max]);

%% the end of this code
