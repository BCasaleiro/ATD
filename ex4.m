%Exercicio 4

close all;
clear all;
clc

%4.1
disp('4.1');
file = 'saxriff.wav';

info = audioinfo(file);
nbits = info.BitsPerSample;
 
[y, fs] = audioread(file);
sound(y, fs);

pause();

ws= 2 * pi * fs;

size_y = length(y);
x = fftshift(fft(y));

if mod(size_y, 2) == 0
    f = -ws/2 : ws/size_y : (ws/2) - (ws/size_y);
else
    f = -ws/2 + ((ws/2)/size_y): ws/size_y: (ws/2) - ((ws/2)/size_y);
end

figure(1), plot(f, abs(x)), title('Magnitude do espectro do sinal');
 
pause();

%4.2
disp('4.2');
time = 46.44;
step = 5.8;

window = round(time * fs/1000);
vstep = round(step * fs/1000);
N = window;
 
if mod(size_y,2) == 0
    f = linspace(-fs/2, fs/2 - fs/window, window);
else
    f = linspace(-fs/2 + fs/(2*window), fs/2 - fs/(2*window), window);
end

matrixsize = 1 : N - vstep : size_y - N;
%Janelas
freqs = zeros(size(matrixsize));
amps = zeros(size(matrixsize));


xaprox = zeros(size(y));

Ts = 1/fs;
t = 0:Ts:size_y*Ts-Ts;

j=1;
for i=1 : N-vstep : size_y-N
    janelaX = fftshift(fft(y(i : i+N-1) .* hamming(N), N));
    Xabsmax = max(abs(janelaX));
    ind = find(abs(janelaX) == Xabsmax);
    
    amps(j) = Xabsmax;
    freqs(j) = f(ind(end));   
    
    A = 10*abs(amps(j))/N;
    xaprox(i : i+N-1) = A * sin(2 * pi * freqs(j) * t(i : i+N-1));
    
    j=j+1;   
end

freqs(freqs <= 100) = 0;

figure(2), plot(1 : time-step : size_y/fs*1000-time, freqs, '.'), title('Frequências fundamentais do sistema');

pause();

%4.3
disp('4.3');

option = menu('Escolha uma opção:', 'Mediana de 5', 'Mediana de 7', 'Mediana de 9');

if(option == 1)
    med = 5;
elseif(option == 2)
    med = 7;
else
    med = 9;
end

x2 = zeros(size_y, 1);

for i=1 : size_y-med
	x2(i) = median(y(i : i+med-1));
end

j=1;
for i=1 : N-vstep : size_y-N
    janelaX = fftshift(fft(x2(i : i+N-1).*hamming(N), N));
    Xabsmax = max(abs(janelaX));
    ind = find(abs(janelaX) == Xabsmax);
    
    amps(j) = Xabsmax;
    freqs(j) = f(ind(end));   
    
    A = 10*abs(amps(j))/N;
    xaprox(i : i+N-1) = A*sin(2*pi*freqs(j)*t(i : i+N-1));
    
    j=j+1;   
end

freqs(freqs <= 100) = 0;

figure(3), plot(1 : time-step : size_y/fs*1000-time, freqs, '.'), title('Frequências fundamentais do sinal com filtro da mediana escolhida');

pause();

%4.4
disp('4.4');
audiowrite(sprintf('som_med_%d.wav', med), xaprox, fs);
somR = audioread(sprintf('som_med_%d.wav', med));
sound(somR, fs);

pause();

figure(4), subplot(211), plot(x2), title('saxriff.wav com filtro da mediana escolhida');

subplot(212), plot(somR), title('saxriff.wav reconstruído com filtro da mediana escolhida');
