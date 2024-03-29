%Exercicio 1

close all;
clear all;
clc

%1.1
wmax = 65 * pi;

w0 = 5 * pi;
f0 = w0 / (2 * pi);
T0 = 1/f0;

ws = 2 * wmax + 2 * w0;
fs = ws / 2 * pi;
Ts = 1/ fs;

t = linspace(0, T0, 500);

xt = -1 + 3 * sin(30 * pi * t) + 4 * cos(20 * pi * t + (pi / 4)) .* sin(45 * pi * t);

%1.2
N = T0 / Ts;
disp(N);
n = 0 : N-1;

xn = -1 + 3 * sin(30 * pi * n * Ts) + 4 * cos(20 * pi * n * Ts + (pi / 4)) .* sin(45 * pi * n * Ts);

%1.3
disp('1.3');
figure(1), plot(t, xt, 'r', n * Ts, xn, 'o'), title('xt e xn');

pause();

%1.4
disp('1.4');
x = fftshift(fft(xn));

if mod(N,2) == 0
    W = linspace(-ws/2, ws/2 - ws/N, N);
else
    W = linspace(-ws/2 + ws/(2*N), ws/2 - ws/(2*N), N);
end

figure(2), subplot(212), stem(W, abs(x)), title('DFT em m�dulo');

subplot(211), plot(W, angle(x), '.'), title('DFT em �ngulo');

pause();

%1.5
disp('1.5');
cm = x / N;

figure(3), subplot(212), stem(W, abs(cm)), title('cm em m�dulo');

subplot(211), plot(W, angle(cm), '.'), title('cm em �ngulo');

pause();

%1.6
disp('1.6');
mid = floor(N/2);
Cm = abs(cm(mid+1:N))*2;
Cm(1) = (Cm(1))/2;
Teta = angle(cm(mid+1:N));

figure(4), subplot(212), stem(0:mid-1,Cm), title('Cm');

subplot(211), plot(0:mid-1,Teta,'.'), title('Tetam');

pause();

%1.7
disp('1.7');
xtR = zeros(size(t));

for m = 0:mid-1
    xtR = xtR + Cm(m+1) * cos(m*w0*t + Teta(m+1));
end

figure(5);
plot(t, xt, t, xtR, 'ro');
legend('Sinal original', 'Sinal reconstru�do');
title('Sinal reconstru�do a partir dos coeficientes');