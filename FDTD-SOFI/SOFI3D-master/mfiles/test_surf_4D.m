clc;clear; close all;

rpol=25/2000;

nthe=51;

npsi=51;

the=linspace(0,pi,nthe);

psi=linspace(0,2*pi,npsi);

% Matriculated;

[psi,the]=meshgrid(psi,the);

rpol=rpol*ones(size(the));

% Convert spherical coordinate to Cartesian Coordinate;

x_mm=rpol.*cos(the)*1000;                     % z in transformational relation;

y_mm=rpol.*sin(the).*cos(psi)*1000;           % x in transformational relation;

z_mm=rpol.*sin(the).*sin(psi)*1000;           % y in transformational relation;

tic 

fval=2*x_mm.^2+y_mm.^2-z_mm.^2;    %% this is a example !!;

toc

%% matrix ( x y z fval ) used to figure fval is color;

%% figure

figure('Renderer','zbuffer','Color',[1 1 1]);

surf(x_mm,y_mm,z_mm,fval);

shading interp;light;lighting gouraud;

colorbar

axis equal;

xlabel('x(mm)');ylabel('y(mm)');zlabel('z(mm)');