clear ,clc;
format long
x=(0:0.01:100)';
y_0=besselj(0,x);
% y_1=besselj(1,x); %һ�ף�����ֻ����0��
% y_2=besselj(2,x); %����


plot(x,y_0);grid on;
axis([0,100,-1,1]);
title('0�ױ�������������ͼ');
xlabel('Variable X');
ylabel('Variable Y');

%��������
hold on;
[up,down] = envelope(x,y_0,'spline');
plot(x, up, 'r');
plot(x, down, 'r');
