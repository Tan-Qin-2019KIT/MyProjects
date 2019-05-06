clear all;close all;
syms f x y
f=x^2+y^2+2*x*y;
symvar(f) %该函数返回的是符号函数中的自变量
g=matlabFunction(f);
z=g(1,1)