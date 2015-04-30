%Exercicio 1

close all;
clear all;
clc

%1.1
wmax = 65 * pi;

w0 = pi;
f0 = w0 / 2 * pi;
T0 = 1/f0;

ws = 2 * wmax + 6 * w0;
fs = ws / 2 * pi;
Ts = 1/ fs;

t = linspace(0, T0, 500);

xt = -1 + 3 * sin(30 * pi * t) + 4 * cos(20 * pi * t + (pi / 4)) .* sin(45 * pi * t);

%1.2
N = Ts / T0;
n = 1 : N-1;

xn = -1 + 3 * sin(30 * pi * n * Ts) + 4 * cos(20 * pi * n * Ts + (pi / 4)) .* sin(45 * pi * n * Ts);

%1.3
figure(1), plot(t, xt, 'r', n * Ts, xn, 'o'), title('xt e xn');
