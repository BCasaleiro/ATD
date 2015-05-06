%Exercicio 2

close all;
clear all;
clc

%2.1
disp('2.1');
nbits = 16;

[y, Fs] = audioread('saxriff.wav');
sound(y, Fs);

pause();

ws= 2 * pi * Fs;
N = 2^nbits;

if mod(N,2) == 0
    f = linspace(-ws/2, ws/2 - ws/N, N);
else
    f = linspace(-ws/2 + ws/(2*N), ws/2 - ws/(2*N), N);
end
    
x = fftshift(fft(y, N));

figure(1), plot(f, abs(x)), title('Espetro em função da frequência');

pause();

%2.2
disp('2.2');

maxAmp = max(abs(x));
indF = find(abs(x) == maxAmp);
freqMax = f(indF(2))/ (2 * pi);

disp('Amplitude Máxima: ');
disp(freqMax);

pause();

%2.3
disp('2.3');
freq_min = 8500;
freq_max = 9500;

ind_pos = find(f >= 2*pi*freq_min & f <= 2*pi*freq_max);
ind_neg = find(f >= -2*pi*freq_max & f <= -2*pi*freq_min);

x_maxAmp = 0.1 * maxAmp * rand(size(ind_pos));
x_angle =  2*pi*rand(size(ind_pos))-pi;

x_ruido_pos = x_maxAmp .* (cos(x_angle) + 1i * sin(x_angle)) ;
x_ruido_neg = conj(x_ruido_pos);

x_ruido = x;
x_ruido(ind_pos) = x(ind_pos)' + x_ruido_pos;
x_ruido(ind_neg) = x(ind_neg)' + x_ruido_neg(end:-1:1) ;

figure(2), plot(f, abs(x_ruido)), title('Espectro do sinal com ruído entre 8500 e 9500'); 

pause();

%2.4
disp('2.4');
x_ruido_aux = real(ifft(ifftshift(x_ruido)));
sound(x_ruido_aux, Fs);

pause();

%2.5
disp('2.5');
fc = 8000;
wn = 2*fc/Fs;
Nf = 6;

[b, a] = butter(Nf, wn);
disp('Numerador:');
disp(b);
disp('Denominador:');
disp(a);

x_filt = filter(b, a, x_ruido_aux);
x_filtrado = fftshift(fft(x_filt));
polos = roots(a);
zeros = roots(b);

disp('Pólos:');
disp(polos);
disp('Zeros:');
disp(zeros);

figure(3), zplane(b, a);

pause();

figure(4), plot(f, abs(x_filtrado)), title('Espectro do sinal com ruído filtrado');

pause();

sound(real(x_filt), Fs);

%2.6
disp('2.6.1');
freq_min = 2000;
freq_max = 3000;

ind_pos = find(f >= 2*pi*freq_min & f <= 2*pi*freq_max);
ind_neg = find(f >= -2*pi*freq_max & f <= -2*pi*freq_min);

x_maxAmp = 0.1 * maxAmp * rand(size(ind_pos));
x_angle =  2*pi*rand(size(ind_pos))-pi;

x_ruido_pos = x_maxAmp .* (cos(x_angle) + 1i * sin(x_angle)) ;
x_ruido_neg = conj(x_ruido_pos);

x_ruido = x;
x_ruido(ind_pos) = x(ind_pos)' + x_ruido_pos;
x_ruido(ind_neg) = x(ind_neg)' + x_ruido_neg(end:-1:1) ;

figure(5), plot(f, abs(x_ruido)), title('Espectro do sinal com ruído entre 2000 e 3000'); 

pause();

disp('2.6.2');

x_ruido_aux = real(ifft(ifftshift(x_ruido)));
sound(x_ruido_aux, Fs);

pause();

disp('2.6.3');

fc_low = 2000;
fc_high = 3000;
wn = [2*fc_low/Fs, 2*fc_high/Fs];
Nf = 6;

[b, a] = butter(Nf, wn, 'stop');
disp('Numerador:');
disp(b);
disp('Denominador:');
disp(a);

x_filt = filter(b, a, x_ruido_aux);
x_filtrado = fftshift(fft(x_filt));
polos = roots(a);
zeros = roots(b);

disp('Pólos:');
disp(polos);
disp('Zeros:');
disp(zeros);

figure(6), zplane(b, a);

pause();

figure(7), plot(f, abs(x_filtrado)), title('Espectro do sinal com ruído filtrado');

pause();

sound(real(x_filt), Fs);