clear all; close all;
format long;
x=(0:0.1:20)';
y_0=besselj(0,x);
y_1=besselj(1,x);
y_2=besselj(2,x);
figure(1)
h=plot(x,y_0,'',x,y_1,'',x,y_2,''); grid on;
axis([0,20,-1,1]);
xlabel('Variable X');
ylabel('Variable J');
title('Bessel function');
legend('Zeroth order','1st order','2nd order');