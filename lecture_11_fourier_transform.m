%% fourier transform of signal

Fs  = 100; %sampling frequency, Hz
Ss = 1/Fs; %sampling period (rate), Samples per second
time_base = -0.5:Ss:0.5; %time base -T/2 to T/2

% original signal
x = sin(2*pi*3*time_base) + .33 * sin(2*pi*9*time_base) ...
    + .2 * sin(2.*pi*15*time_base); 

% What happens when you add some phase shift to the signal?
x = sin(2*pi*3*time_base) + .33 * sin(2*pi*9*time_base+(2*pi/7)) ...
    + .2 * sin(2.*pi*15*time_base+(2*pi/9)); 

figure;
plot(time_base,x);

% incresing size of frequency basis
% does this provide more frequency information?
NFFT = length(x); % n = M
% NFFT = 2^nextpow2(length(x)); % using powers of 2 allows fast calculation
% NFFT = 2^10;

freq_base = Fs/2*linspace(0,1,NFFT/2+1); % single sided spectrum
X = fft(x, NFFT);

figure;
stem(freq_base,abs(X(1:NFFT/2+1)),'k');
xlabel('Frequency (Hz)')
ylabel('Power')

x1 = ifft(X);
figure;
plot(time_base,x1);

% figure;
% plot(-atan(real(X)./imag(X)));

%% trade off between time and frequency domains


x = zeros(1,250);
n = [0 5 10 25 125];% 0 = delta function, larger numbers = larger step/box functions
Fs = 500; %sampling frequency

for i=1:5 
     
    % generate delta/step functions of different sizes
    if i == 5
        x(length(x)/2-n(i)+1:length(x)/2+n(i)) = 1;
    else
        x(length(x)/2-n(i):length(x)/2+n(i)) = 1;
    end;
    
    L = length(x); % lenght of signal
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y, fastest calculation
    
    f = Fs*linspace(-1,1,NFFT); % setup frequncy base
    X = fft(x,NFFT)/L; % calculate fourier transform
    X = fftshift(X); % center spectrum
   
     % plot function & power spectrum.
    figure;
    
    ax1 = subplot(2,1,1);
    stem(x,'k');
    title('Time Domain')
    xlabel('Time')
    ylabel('Amplitude')
    
    ax2 = subplot(2,1,2);
    plot(f,abs(X),'r')
    title('Frequency Domain')
    xlabel('Frequency (Hz)')
    ylabel('Power')
end;