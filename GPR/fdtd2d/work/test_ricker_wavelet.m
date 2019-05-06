clear all; close all;
fm=100;%��Ƶ
dt=0.0005;%ʱ����������
t0=0.05;
number=201;%��������
t=(1:number)*dt;
a=(1-2*(pi*fm*(t-t0)).^2).*exp(-(pi*fm*(t-t0)).^2);
subplot(2,1,1);
plot(t,a);
title('Ricker-�׿��Ӳ�');
xlabel('ʱ��t��ms��');
ylabel('��ֵA');

df=1/dt;
f=(1:(number+1)/2)*df/number;
source_freq=abs(fft(a));  % calculate the frequency spectrum
Y=source_freq(1:(number+1)/2);
% Y=abs(fft(a));%fourier�任��ȡ�����
subplot(2,1,2)
plot(f,Y);
title('Ricker�Ӳ��������');
xlabel('Ƶ��f��hz��');
ylabel('�����');

%ע��ʱ����������Ϊ0.001s����������Ϊ100�㣬�ܵ�ʱ�䳤��Ϊ0.1s����Ƶ����������Ϊ10hz��