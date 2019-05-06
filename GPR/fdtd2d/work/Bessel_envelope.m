clear ,clc;
format long
x=(0:0.01:100)';
y_0=besselj(0,x);
% y_1=besselj(1,x); %一阶，这里只画了0阶
% y_2=besselj(2,x); %二阶


plot(x,y_0);grid on;
axis([0,100,-1,1]);
title('0阶贝塞尔函数曲线图');
xlabel('Variable X');
ylabel('Variable Y');

%画包络线
hold on;
[up,down] = envelope(x,y_0,'spline');
plot(x, up, 'r');
plot(x, down, 'r');
