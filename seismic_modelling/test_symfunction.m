clear all;close all;
syms f x y
f=x^2+y^2+2*x*y;
symvar(f) %�ú������ص��Ƿ��ź����е��Ա���
g=matlabFunction(f);
z=g(1,1)