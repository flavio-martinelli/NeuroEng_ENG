%% Working with the ENG rat data recorded from the experiments with cuff electrodes

clear all
close all
clc

%% Data Import 

d=importdata('Brush.txt');


%% Raw data and labels
signal=d(:,1);
labels=d(:,2);
labels=round(labels); % Because there could be some noise in the analogical signal. 
Fs = 20000;

figure(1)
plot(zscore(signal));
hold on;
title('ENG original signal and stimuli application labels')
xlim([1 , length(signal)]);
ylabel('ENG');
xlabel('Samples')

yyaxis right
plot(zscore(labels),'r');
ylim([-5 5])
ylabel('Labels: -1 = rest, +1 = activity');


%% Examination of the spectral context of the signal-PSD

[rows_act,cols_act,values_act] = find(labels>0);
[rows_rest,cols_rest,values_rest] = find(labels==0);

figure(2)
signalOfInterest=signal(rows_act);
notOfInterest=signal(rows_rest);
h = spectrum.welch; % creates the Welch spectrum estimator
SOIf=psd(h,signalOfInterest,'Fs',Fs); % calculates and plot the one sided PSD
NOIf=psd(h,notOfInterest,'Fs',Fs);
plot(SOIf); % Plot the one-sided PSD. 
hold on;
plot(NOIf);
title('Spectrum of SOI, and NOI for the stimuli application')

%% Filtering 

%Determine the Low pass cutoff frequency for your filter
Fl=800;
%Determine the High pass cutoff frequency for your filter
Fh=2200;


[Signal_filtered]=filtra(signal',Fs,Fl,Fh);

figure(3)
soi_f=Signal_filtered(rows_act);
noi_f=Signal_filtered(rows_rest);
h1 = spectrum.welch; % creates the Welch spectrum estimator
SOI_f=psd(h,soi_f,'Fs',Fs); % calculates and plot the one sided PSD
NOI_f=psd(h,noi_f,'Fs',Fs);
plot(SOI_f); % Plot the one-sided PSD. 
hold on;
plot(NOI_f);
title('Spectrum of SOI, and NOI for VF')

%% Feature extraction 

% Determine the step size

step=0.1 * Fs; %  step size for 100 ms window

% Compute the MAV and Zero-Crossing features for each time window                      
% Hint: MAV and ZeroCross can be functions in separate files         


L2=length(signal);

for i=1:step:(L2-step)
    % Compute MAV and ZCross features and store them in a 2D
    % "features" vector
    MAVI((i-1)/step +1) = MAV(signal(i:i+step));
    ZCrossi((i-1)/step +1) = ZCross(signal(i:i+step));
    
end

% Resizing the label vector to the size of the features vector
for i=1:step:(L2-step)  
    z=ceil(i/step);
    labels_resized(z)=labels(i);
end

% plot the zscored features along with the corresponding labels



%% Signal-to-Noise Estimation

% Extract the features for the non zero labels

% Compute Activation 

% Extract the features for the labels = zero

% Compute the noise, SNR and SNR in dB
