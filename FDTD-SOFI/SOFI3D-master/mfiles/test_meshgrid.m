close all;
clear all;
clc;
[x,y,z] = meshgrid(-2:.2:2, -2:.2:4, -2:.2:6);
v = x .* exp(-x.^2 - y.^2 - z.^2);
slice(x,y,z,v,[-1.2 .8 2],2,[-2 -.2])
        xlabel('x in m');
        ylabel('y in m');
        zlabel('Depth z in m');