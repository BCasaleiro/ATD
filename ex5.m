%5.1 at? 5.4
option = menu('Escolha um som:', 'escala.wav', 'piano.wav', 'flauta.wav');

if(option == 1)
    som = 'escala.wav';
elseif(option == 2)
    som = 'piano.wav';
elseif(option == 3)
    som = 'flauta.wav';
end

[x, fs] = audioread(som);
info = audioinfo(som);
nbits = info.BitsPerSample;

if size(x, 2) == 2
    temp = zeros(length(x), 1);
    for i = 1:length(x)
        temp(i) = (x(i, 1) + x(i, 2)) / 2;
    end
    x = temp;
end

sound(x, fs);
X = fftshift(fft(x));
sizeofX = length(x);

if mod(sizeofX, 2) == 0
    w = -fs*pi : 2*pi*fs/sizeofX : fs*pi-2*pi*fs/sizeofX;
else
    w = -fs*pi+fs*pi/sizeofX : 2*pi*fs/sizeofX : fs*pi - fs*pi/sizeofX;
end

subplot(2, 1, 1);
plot(x);
title(sprintf('sinal: %s', som));
subplot(2, 1, 2);
plot(w, abs(X), '.');
title('Magnitude');

disp('Press any key to continue');
pause();

frequencia_notas = [262 227 294 311 330 349 370 392 415 440 466 494];
nome_notas = {'Do   ';'Do#  ';'Re   ';'Re#  ';'Mi   ';'Fa   ';'Fa#  ';'Sol  ';'Sol# ';'La   ';'La#  ';'Si   '};

janela = round(100/1000*fs);
sobreposicao = round(12.5/1000*fs);

N = janela;

if (mod(sizeofX, 2) == 0)
    f = linspace(-fs/2, fs/2 - fs/N, N);
else
    f = linspace(-fs/2 + fs/(2*N), fs/2 - fs/(2*N), N);
end


matrix_sizeofX = 1: N-sobreposicao: sizeofX-N;

freqs = zeros(size(matrix_sizeofX));
amps = zeros(size(matrix_sizeofX));

j=1;
for i = 1: N-sobreposicao: sizeofX-N
    janela_x = fftshift(fft(x(i : i+N-1).* hamming(N), N));
    x_abs_max = max(abs(janela_x));
    ind = find(abs(janela_x) == x_abs_max);
    
    amps(j) = x_abs_max;
    freqs(j) = f(ind(floor(end/2) + 1));
    
    j = j + 1;
end

figure(2);
plot(1: janela*1000/fs-sobreposicao*1000/fs : sizeofX/fs*1000 - janela*1000/fs, freqs, '.');
title('frequencias fundamentais');
xlabel('ms');
ylabel('Hz');
disp('Press any key to continue');
pause();

fprintf('Resolucao em Frequencia do som %s = %.2f \n', som, fs/N);
disp('Press any key to continue');
pause();

fprintf('Notas do som %s\n', som);
for i=1 : length(freqs)
    freq = freqs(i);
    if freq ~= 0
        while(freq < frequencia_notas(1))
            freq = freq*2;
        end
        while(freq > frequencia_notas(end))
            freq = freq/2;
        end
        for j=1: length(frequencia_notas)
            if(freq < frequencia_notas(j))
                break;
            end
        end
        if(j ~= 1)
            if(abs(freq-frequencia_notas(j-1)) <= abs(freq - frequencia_notas(j)))
                j = j-1;
            end
        end
        fprintf('%s \n', nome_notas{j});
    end
end
        
        
        
        
        
        







