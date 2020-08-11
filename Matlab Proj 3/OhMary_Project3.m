%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Project 3
%   Author: Mary Oh
%   Date: 2016/11/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
%close all;
%clear all;
%clc;

ASUID = 4; % Input the last digit of your ASU ID number.

%% Question 1

% Variable and signal initialization
t = 0:0.002:30; % time vector
mean_noise = 0; % mean of noise distribution
sigma_noise = ((ASUID+3)*3)/4; % standard deviation of noise distribution.
noise = mean_noise + sigma_noise*randn(size(t)); % Noise vector

s1 = sin(2*pi*12*t) + sin(2*pi*25*t) + sin(2*pi*60*t) + sin(2*pi*110*t) + noise; %signal

figure (1);
subplot(2,1,1)
plot(t,s1) % Plot Noisy Signal
xlabel('Time (sec)','fontsize',12)
ylabel('Amplitude x(t)','fontsize',12)
title('Noisy signal 1')

subplot(2,1,2)
plot(-250:1/30:250,abs(fftshift(fft(s1)))) % Plots the frequency profile of the signal s1
xlabel('Frequency (Hz)','fontsize',12)
ylabel('Amplitude X(j\omega)','fontsize',12)
title('FT of Noisy signal 1')
xlim([0 250]);

save('signal1','s1') % Save s1 signal in current directory
save('t','t'); % Save time vector in current directory

%% 1a

% Create 4 filters with sampling rate (Fs) 500Hz and the cut-off frequencies:
%   10Hz, 25Hz, 70Hz, 120Hz
% For the other parameters, use the parameters shown in the Tutorial 4
% lecture.

% Export SOS Matrix and Scale Values to the workspace.
% Use the following code to save the variables to .mat files.
% Change the names of the .mat filename for each filter

save('coeffs_10Hz.mat', 'SOS') % Change SOS to the name of your SOS Matrix variable
save('G_10Hz.mat', 'G') % Change G to the name of your Scale Values variable

save('coeffs_25Hz.mat', 'SOS1')
save('G_25Hz.mat', 'G1')

save('coeffs_70Hz.mat', 'SOS2')
save('G_70Hz.mat', 'G2')

save('coeffs_120Hz.mat', 'SOS3')
save('G_120Hz.mat', 'G3')
%% 1b
load signal1.mat % Load the saved noisy signal
load t.mat % Load the saved t
load coeffs_10Hz.mat; % Load to the filter coefficients. Change name according to what you saved the filer as.
load G_10Hz.mat; % Load the filter gain. Change name of the G variable according to what you saved it as.

Fs = 1/0.002; % sampling frequency
f = linspace(-Fs/2,Fs/2,length(s1)); % Vector of frequencies to plot signals

[num, denom] = sos2tf(SOS, G); % Change "SOS" and "G" to the names you chose while exporting the filter
y = filter(num, denom, s1); % applies filter to noisy signal

figure (2); % Time-Domain Figure
subplot(2,1,1)
plot(t, s1)
xlabel('Time (ms)')
ylabel('x(t)','fontsize',12)
ylim([min(s1)-1 max(s1)+1])
title('Original Signal - Mary Oh 1208315416')

subplot(2,1,2)
plot(t, y, 'r') % Plot filtered signal
ylim([min(s1)-1 max(s1)+1])
xlabel('Time (ms)')
ylabel('y(t)','fontsize',12)
title('Filtered signal. Fc = 10 Hz - Mary Oh 1208315416') % Change title according to the Fc used for each signal

figure (3); % Frequency-Domain Figure
subplot(2,1,1)
signal_FT = abs(fftshift(fft(s1))); % Compute spectrum of noisy signal
plot(f, signal_FT)
xlim([0 Fs/2]);
ylim([min(signal_FT)-100 max(signal_FT)+100]);
xlabel('Frequency (Hz)','fontsize',12)
ylabel('Amplitude X(w)','fontsize',12)
title('FT of Original Noisy signal 1 - Mary Oh 1208315416')

subplot(2,1,2)
filt_signal_FT = abs(fftshift(fft(y))); % Compute spectrum of filtered signal
plot(f, filt_signal_FT,'r')
xlim([0 Fs/2]);
ylim([min(signal_FT)-100 max(signal_FT)+100]);
xlabel('Frequency (Hz)','fontsize',12)
ylabel('Amplitude Y(w)','fontsize',12)
title('FT of Filtered signal 1, Fc = 10Hz - Mary Oh 1208315416') % Change title according to the Fc used for each signal


load coeffs_25Hz.mat; % Load to the filter coefficients. Change name according to what you saved the filer as.
load G_25Hz.mat; % Load the filter gain. Change name of the G variable according to what you saved it as.

[num1, denom1] = sos2tf(SOS1, G1); % Change "SOS" and "G" to the names you chose while exporting the filter
y1 = filter(num1, denom1, s1); % applies filter to noisy signal

figure (4); % Time-Domain Figure
subplot(2,1,1)
plot(t, s1)
xlabel('Time (ms)')
ylabel('x(t)','fontsize',12)
ylim([min(s1)-1 max(s1)+1])
title('Original Signal - Mary Oh 1208315416')

subplot(2,1,2)
plot(t, y1, 'r') % Plot filtered signal
ylim([min(s1)-1 max(s1)+1])
xlabel('Time (ms)')
ylabel('y(t)','fontsize',12)
title('Filtered signal. Fc = 25 Hz - Mary Oh 1208315416') % Change title according to the Fc used for each signal

figure (5); % Frequency-Domain Figure
subplot(2,1,1)
signal_FT = abs(fftshift(fft(s1))); % Compute spectrum of noisy signal
plot(f, signal_FT)
xlim([0 Fs/2]);
ylim([min(signal_FT)-100 max(signal_FT)+100]);
xlabel('Frequency (Hz)','fontsize',12)
ylabel('Amplitude X(w)','fontsize',12)
title('FT of Original Noisy signal 1 - Mary Oh 1208315416')

subplot(2,1,2)
filt_y1_FT = abs(fftshift(fft(y1))); % Compute spectrum of filtered signal
plot(f, filt_y1_FT,'r')
xlim([0 Fs/2]);
ylim([min(signal_FT)-100 max(signal_FT)+100]);
xlabel('Frequency (Hz)','fontsize',12)
ylabel('Amplitude Y(w)','fontsize',12)
title('FT of Filtered signal 1, Fc = 25Hz - Mary Oh 1208315416') % Change title according to the Fc used for each signal


load coeffs_70Hz.mat; % Load to the filter coefficients. Change name according to what you saved the filer as.
load G_70Hz.mat; % Load the filter gain. Change name of the G variable according to what you saved it as.

[num2, denom2] = sos2tf(SOS2, G2); % Change "SOS" and "G" to the names you chose while exporting the filter
y2 = filter(num2, denom2, s1); % applies filter to noisy signal

figure (6); % Time-Domain Figure
subplot(2,1,1)
plot(t, s1)
xlabel('Time (ms)')
ylabel('x(t)','fontsize',12)
ylim([min(s1)-1 max(s1)+1])
title('Original Signal - Mary Oh 1208315416')

subplot(2,1,2)
plot(t, y2, 'r') % Plot filtered signal
ylim([min(s1)-1 max(s1)+1])
xlabel('Time (ms)')
ylabel('y(t)','fontsize',12)
title('Filtered signal. Fc = 70Hz - Mary Oh 1208315416') % Change title according to the Fc used for each signal

figure (7); % Frequency-Domain Figure
subplot(2,1,1)
signal_FT = abs(fftshift(fft(s1))); % Compute spectrum of noisy signal
plot(f, signal_FT)
xlim([0 Fs/2]);
ylim([min(signal_FT)-100 max(signal_FT)+100]);
xlabel('Frequency (Hz)','fontsize',12)
ylabel('Amplitude X(w)','fontsize',12)
title('FT of Original Noisy signal 1 - Mary Oh 1208315416')

subplot(2,1,2)
filt_y2_FT = abs(fftshift(fft(y2))); % Compute spectrum of filtered signal
plot(f, filt_y2_FT,'r')
xlim([0 Fs/2]);
ylim([min(signal_FT)-100 max(signal_FT)+100]);
xlabel('Frequency (Hz)','fontsize',12)
ylabel('Amplitude Y(w)','fontsize',12)
title('FT of Filtered signal 1, Fc = 70Hz - Mary Oh 1208315416') % Change title according to the Fc used for each signal


load coeffs_120Hz.mat; % Load to the filter coefficients. Change name according to what you saved the filer as.
load G_120Hz.mat; % Load the filter gain. Change name of the G variable according to what you saved it as.

[num3, denom3] = sos2tf(SOS3, G3); % Change "SOS" and "G" to the names you chose while exporting the filter
y3 = filter(num3, denom3, s1); % applies filter to noisy signal

figure (8); % Time-Domain Figure
subplot(2,1,1)
plot(t, s1)
xlabel('Time (ms)')
ylabel('x(t)','fontsize',12)
ylim([min(s1)-1 max(s1)+1])
title('Original Signal - Mary Oh 1208315416')

subplot(2,1,2)
plot(t, y3, 'r') % Plot filtered signal
ylim([min(s1)-1 max(s1)+1])
xlabel('Time (ms)')
ylabel('y(t)','fontsize',12)
title('Filtered signal. Fc = 120Hz - Mary Oh 1208315416') % Change title according to the Fc used for each signal

figure (9); % Frequency-Domain Figure
subplot(2,1,1)
signal_FT = abs(fftshift(fft(s1))); % Compute spectrum of noisy signal
plot(f, signal_FT)
xlim([0 Fs/2]);
ylim([min(signal_FT)-100 max(signal_FT)+100]);
xlabel('Frequency (Hz)','fontsize',12)
ylabel('Amplitude X(w)','fontsize',12)
title('FT of Original Noisy signal 1 - Mary Oh 1208315416')

subplot(2,1,2)
filt_y3_FT = abs(fftshift(fft(y3))); % Compute spectrum of filtered signal
plot(f, filt_y3_FT,'r')
xlim([0 Fs/2]);
ylim([min(signal_FT)-100 max(signal_FT)+100]);
xlabel('Frequency (Hz)','fontsize',12)
ylabel('Amplitude Y(w)','fontsize',12)
title('FT of Filtered signal 1, Fc = 120Hz - Mary Oh 1208315416') % Change title according to the Fc used for each signal

%% Question 2
load emg.mat % Loads EMG data vector and its time vector

figure (10)
subplot(2,1,1)
plot(t,emg1)
xlim([0 t(end)])
xlabel('Time (sec)')
ylabel('Voltage (mV)')
title('EMG1 Signal','fontsize',13)

subplot(2,1,2)
plot(t,emg2)
xlim([0 t(end)])
xlabel('Time (sec)')
ylabel('Voltage (mV)')
title('EMG2 Signal','fontsize',13)

%% 2a

Fs = 1/0.002; % sampling frequency
f1 = linspace(-Fs/2,Fs/2,length(emg1));
f2 = linspace(-Fs/2,Fs/2,length(emg2));

emg1_FT = abs(fftshift(fft(emg1)));
emg2_FT = abs(fftshift(fft(emg2)));

figure (11)
subplot(2,1,1)
plot(t, emg1)
xlim([0 t(end)])
xlabel('Time (sec)')
ylabel('Voltage (mV)')
title('EMG1 Signal')

subplot(2,1,2)
plot(f1, emg1_FT)
xlim([0 550])
xlabel('Frequency (Hz)')
ylabel('Amplitude X(w)')
title('FT of EMG1 Signal')

figure (12)
subplot(2,1,1)
plot(t, emg2)
xlim([0 t(end)])
xlabel('Time (sec)')
ylabel('Voltage (mV)')
title('EMG2 Signal')

subplot(2,1,2)
plot(f2, emg2_FT)
xlim([0 550])
xlabel('Frequency (Hz)')
ylabel('Amplitude X(w)')
title('FT of EMG2 Signal')

%% 2b

MNF1 = sum(abs(f1).*emg1_FT)/sum(f1.*t) % Change "ft_emg" to your variable name
% that contains fourier transform of emg1 and emg2
MNF2 = sum(abs(f2).*emg2_FT)/sum(f2.*t)

%% 2c

save('coeffs_emg.mat', 'SOS_emg')
save('G_emg.mat', 'G_emg')

load coeffs_emg.mat;
load G_emg.mat;

[num_emg, denom_emg] = sos2tf(SOS_emg, G_emg); % Change "SOS" and "G" to the names you chose while exporting the filter
y_emg1 = filter(num_emg, denom_emg, emg1);
y_emg2 = filter(num_emg, denom_emg, emg2);

figure (13); % Time-Domain Figure
subplot(2,1,1)
plot(t, emg1)
xlabel('Time (sec)')
ylabel('Voltage (mV)','fontsize',12)
ylim([min(emg1)-1 max(emg1)+1])
title('Original EMG1 Signal - Mary Oh 1208315416')

subplot(2,1,2)
plot(f1, y_emg1, 'r') % Plot filtered signal
ylim([min(emg1)-1 max(emg1)+1])
xlabel('Time (sec)')
ylabel('Voltage (mV)','fontsize',12)
title('Filtered EMG1 signal - Mary Oh 1208315416') % Change title according to the Fc used for each signal

figure (14); % Frequency-Domain Figure
subplot(2,1,1)
plot(f1, emg1_FT)
xlim([0 Fs/2]);
ylim([min(emg1_FT)-100 max(emg1_FT)+100]);
xlabel('Frequency (Hz)','fontsize',12)
ylabel('Amplitude X(w)','fontsize',12)
title('FT of Original Noisy signal 1 - Mary Oh 1208315416')

subplot(2,1,2)
filt_emg1_FT = abs(fftshift(fft(y_emg1))); % Compute spectrum of filtered signal
plot(f1, filt_emg1_FT,'r')
xlim([0 Fs/2]);
ylim([min(emg1_FT)-100 max(emg1_FT)+100]);
xlabel('Frequency (Hz)','fontsize',12)
ylabel('Amplitude Y(w)','fontsize',12)
title('FT of Filtered EMG1 signal - Mary Oh 1208315416')

figure (15); % Time-Domain Figure
subplot(2,1,1)
plot(t, emg2)
xlabel('Time (sec)')
ylabel('Voltage (mV)','fontsize',12)
ylim([min(emg2)-1 max(emg2)+1])
title('Original EMG2 Signal - Mary Oh 1208315416')

subplot(2,1,2)
plot(f2, y_emg2, 'r') % Plot filtered signal
ylim([min(emg2)-1 max(emg2)+1])
xlabel('Time (sec)')
ylabel('Voltage (mV)','fontsize',12)
title('Filtered EMG2 signal - Mary Oh 1208315416') % Change title according to the Fc used for each signal

figure (16); % Frequency-Domain Figure
subplot(2,1,1)
plot(f2, emg2_FT)
xlim([0 Fs/2]);
ylim([min(emg2_FT)-100 max(emg2_FT)+100]);
xlabel('Frequency (Hz)','fontsize',12)
ylabel('Amplitude X(w)','fontsize',12)
title('FT of Original Noisy EMG2 signal - Mary Oh 1208315416')

subplot(2,1,2)
filt_emg2_FT = abs(fftshift(fft(y_emg2))); % Compute spectrum of filtered signal
plot(f2, filt_emg2_FT,'r')
xlim([0 Fs/2]);
ylim([min(emg2_FT)-100 max(emg2_FT)+100]);
xlabel('Frequency (Hz)','fontsize',12)
ylabel('Amplitude Y(w)','fontsize',12)
title('FT of Filtered EMG2 signal - Mary Oh 1208315416')

MNF_FT1 = sum(abs(f1).*filt_emg1_FT)/sum(f1.*t)
MNF_FT2 = sum(abs(f2).*filt_emg2_FT)/sum(f2.*t)

%% 2d

save('coeffs_258Hz.mat', 'SOS_d1')
save('G_258Hz.mat', 'G_d1')

load coeffs_258Hz.mat;
load G_258Hz.mat;

[num_d1, denom_d1] = sos2tf(SOS_d1, G_d1); % Change "SOS" and "G" to the names you chose while exporting the filter
y_d11 = filter(num_d1, denom_d1, emg1);
y_d12 = filter(num_d1, denom_d1, emg2);

figure (17); % Time-Domain Figure
subplot(2,1,1)
plot(t, emg1)
xlabel('Time (sec)')
ylabel('Voltage (mV)','fontsize',12)
ylim([min(emg1)-1 max(emg1)+1])
title('Original EMG1 Signal')

subplot(2,1,2)
plot(f1, y_d11, 'r') % Plot filtered signal
ylim([min(emg1)-1 max(emg1)+1])
xlabel('Time (sec)')
ylabel('Voltage (mV)','fontsize',12)
title('Filtered EMG1 signal - Fhigh = ((ASUID*2)+250)Hz') % Change title according to the Fc used for each signal

figure (18); % Frequency-Domain Figure
subplot(2,1,1)
plot(f1, emg1_FT)
xlim([0 Fs/2]);
ylim([min(emg1_FT)-100 max(emg1_FT)+100]);
xlabel('Frequency (Hz)','fontsize',12)
ylabel('Amplitude X(w)','fontsize',12)
title('FT of Original Noisy EMG1 signal')

subplot(2,1,2)
filt_d11_FT = abs(fftshift(fft(y_d11))); % Compute spectrum of filtered signal
plot(f1, filt_d11_FT,'r')
xlim([0 Fs/2]);
ylim([min(emg1_FT)-100 max(emg1_FT)+100]);
xlabel('Frequency (Hz)','fontsize',12)
ylabel('Amplitude Y(w)','fontsize',12)
title('FT of Filtered EMG1 signal - Fhigh = ((ASUID*2)+250)Hz')

figure (19); % Time-Domain Figure
subplot(2,1,1)
plot(t, emg2)
xlabel('Time (sec)')
ylabel('Voltage (mV)','fontsize',12)
ylim([min(emg2)-1 max(emg2)+1])
title('Original EMG2 Signal')

subplot(2,1,2)
plot(f2, y_d12, 'r') % Plot filtered signal
ylim([min(emg2)-1 max(emg2)+1])
xlabel('Time (sec)')
ylabel('Voltage (mV)','fontsize',12)
title('Filtered EMG2 signal - Fhigh = ((ASUID*2)+250)Hz') % Change title according to the Fc used for each signal

figure (20); % Frequency-Domain Figure
subplot(2,1,1)
plot(f2, emg2_FT)
xlim([0 Fs/2]);
ylim([min(emg2_FT)-100 max(emg2_FT)+100]);
xlabel('Frequency (Hz)','fontsize',12)
ylabel('Amplitude X(w)','fontsize',12)
title('FT of Original Noisy EMG2 signal')

subplot(2,1,2)
filt_d12_FT = abs(fftshift(fft(y_d12))); % Compute spectrum of filtered signal
plot(f2, filt_d12_FT,'r')
xlim([0 Fs/2]);
ylim([min(emg2_FT)-100 max(emg2_FT)+100]);
xlabel('Frequency (Hz)','fontsize',12)
ylabel('Amplitude Y(w)','fontsize',12)
title('FT of Filtered EMG2 signal - Fhigh = ((ASUID*2)+250)Hz')

MNF_d11 = sum(abs(f1).*filt_d11_FT)/sum(f1.*t)
MNF_d12 = sum(abs(f2).*filt_d12_FT)/sum(f2.*t)



save('coeffs_308Hz.mat', 'SOS_d2')
save('G_308Hz.mat', 'G_d2')

load coeffs_308Hz.mat;
load G_308Hz.mat;

[num_d2, denom_d2] = sos2tf(SOS_d2, G_d2); % Change "SOS" and "G" to the names you chose while exporting the filter
y_d21 = filter(num_d2, denom_d2, emg1);
y_d22 = filter(num_d2, denom_d2, emg2);

figure (21); % Time-Domain Figure
subplot(2,1,1)
plot(t, emg1)
xlabel('Time (sec)')
ylabel('Voltage (mV)','fontsize',12)
ylim([min(emg1)-1 max(emg1)+1])
title('Original EMG1 Signal')

subplot(2,1,2)
plot(f1, y_d21, 'r') % Plot filtered signal
ylim([min(emg1)-1 max(emg1)+1])
xlabel('Time (sec)')
ylabel('Voltage (mV)','fontsize',12)
title('Filtered EMG1 signal - Fhigh = ((ASUID*2)+300)') % Change title according to the Fc used for each signal

figure (22); % Frequency-Domain Figure
subplot(2,1,1)
plot(f1, emg1_FT)
xlim([0 Fs/2]);
ylim([min(emg1_FT)-100 max(emg1_FT)+100]);
xlabel('Frequency (Hz)','fontsize',12)
ylabel('Amplitude X(w)','fontsize',12)
title('FT of Original Noisy EMG1 signal')

subplot(2,1,2)
filt_d21_FT = abs(fftshift(fft(y_d21))); % Compute spectrum of filtered signal
plot(f1, filt_d21_FT,'r')
xlim([0 Fs/2]);
ylim([min(emg1_FT)-100 max(emg1_FT)+100]);
xlabel('Frequency (Hz)','fontsize',12)
ylabel('Amplitude Y(w)','fontsize',12)
title('FT of Filtered EMG1 signal - Fhigh = ((ASUID*2)+300)')

figure (23); % Time-Domain Figure
subplot(2,1,1)
plot(t, emg2)
xlabel('Time (sec)')
ylabel('Voltage (mV)','fontsize',12)
ylim([min(emg2)-1 max(emg2)+1])
title('Original EMG2 Signal')

subplot(2,1,2)
plot(f2, y_d22, 'r') % Plot filtered signal
ylim([min(emg2)-1 max(emg2)+1])
xlabel('Time (sec)')
ylabel('Voltage (mV)','fontsize',12)
title('Filtered EMG2 signal - Fhigh = ((ASUID*2)+300)') % Change title according to the Fc used for each signal

figure (24); % Frequency-Domain Figure
subplot(2,1,1)
plot(f2, emg2_FT)
xlim([0 Fs/2]);
ylim([min(emg2_FT)-100 max(emg2_FT)+100]);
xlabel('Frequency (Hz)','fontsize',12)
ylabel('Amplitude X(w)','fontsize',12)
title('FT of Original Noisy EMG2 signal')

subplot(2,1,2)
filt_d22_FT = abs(fftshift(fft(y_d22))); % Compute spectrum of filtered signal
plot(f2, filt_d22_FT,'r')
xlim([0 Fs/2]);
ylim([min(emg2_FT)-100 max(emg2_FT)+100]);
xlabel('Frequency (Hz)','fontsize',12)
ylabel('Amplitude Y(w)','fontsize',12)
title('FT of Filtered EMG2 signal - Fhigh = ((ASUID*2)+300)')

MNF_d21 = sum(abs(f1).*filt_d21_FT)/sum(f1.*t)
MNF_d22 = sum(abs(f2).*filt_d22_FT)/sum(f2.*t)
% EOF