%Exercicio 2

close all;
clear all;
clc

%2.1
disp('2.1');

[y, Fs] = audioread('saxriff.wav');
%sound(y, Fs);

ws= 2 * pi * Fs;
f = linspace( -ws/2, ws/2, length(y));

x = fftshift(fft(y));

% Plot single-sided amplitude spectrum.
plot(f, abs(x)) 
title('Espetro em função da frequência')

pause();

%2.2
disp('2.2');

maxAmp = max(abs(x));
indF = find(abs(x) == maxAmp);
freqMax = f(indF(2))/ (2 * pi);

disp('Amplitude Máxima: ');
disp(freqMax);

%2.3
disp('2.3');


