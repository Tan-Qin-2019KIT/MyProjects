clear all; close all;
Fs = 8; % 采样频率
Ts = 1 / Fs; % 采样时间间隔
L = 32; % length of signal
t = (0 : (L - 1)) * Ts; % discrete time
x = 2 + 3 * cos(2 * pi * 1 * t - 30 * pi / 180); % original signal

N = 2 ^ nextpow2(L); % 采样点数
X = fft(x, N); 
f = Fs / N * (0 : (N - 1)); % 频率
x_hat = ifft(X, N);
subplot(311)
stem(t, x);
title('original signal'),xlabel('time'),ylabel('amplitude')

subplot(312)
X(1) = X(1) / 2;
stem(f, abs(X / N * 2));
title('amplitude spectrum'),xlabel('frequency'),ylabel('amplitude')

subplot(313)
stem(t, x_hat);
title('recovered signal'),xlabel('time'),ylabel('amplitude')
