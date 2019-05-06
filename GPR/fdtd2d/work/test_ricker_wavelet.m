clear all; close all;
fm=100;%主频
dt=0.0005;%时间域采样间隔
t0=0.05;
number=201;%采样点数
t=(1:number)*dt;
a=(1-2*(pi*fm*(t-t0)).^2).*exp(-(pi*fm*(t-t0)).^2);
subplot(2,1,1);
plot(t,a);
title('Ricker-雷克子波');
xlabel('时间t（ms）');
ylabel('幅值A');

df=1/dt;
f=(1:(number+1)/2)*df/number;
source_freq=abs(fft(a));  % calculate the frequency spectrum
Y=source_freq(1:(number+1)/2);
% Y=abs(fft(a));%fourier变换，取振幅谱
subplot(2,1,2)
plot(f,Y);
title('Ricker子波的振幅谱');
xlabel('频率f（hz）');
ylabel('振幅谱');

%注：时间域采样间隔为0.001s，采样点数为100点，总的时间长度为0.1s，则频率域采样间隔为10hz。