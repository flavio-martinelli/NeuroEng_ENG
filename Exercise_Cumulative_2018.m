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
plot(zscore(labels),'r');
title('ENG original signal and stimuli application labels')

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
Fl=;
%Determine the High pass cutoff frequency for your filter
Fh=;


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
step=; 

% Compute the MAV and Zero-Crossing features for each time window                      
% Hint: MAV and ZeroCross can be functions in separate files         


L2=length(signal);
for i=1:step:(L2-step)
    % Compute MAV and ZCross features and store them in a 2D
    % "features" vector
end

% Resizing the label vector to the size of the features vector
for i=1:step:(L2-step);  
    z=ceil(i/step);
    labels_resized(z)=labels(i);
end

% plot the zscored features along with the corresponding labels



%% Signal-to-Noise Estimation

% Extract the features for the non zero labels

% Compute Activation 

% Extract the features for the labels = zero

% Compute the noise, SNR and SNR in dB
