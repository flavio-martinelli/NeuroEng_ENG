%% Working with the ENG rat data recorded from the experiments with cuff electrodes
clear all
close all
clc

%% Data Import
filenames = {'Brush.txt','FootFlexion.txt','Touch.txt'};
names = {'VF filament','Flexion','Touch'};

f_map = containers.Map(filenames,names);

for k = 1:length(filenames)
   
    d=importdata(filenames{k});
     
    %% Raw data and labels
    signal=d(:,1);
    labels=d(:,2);
    labels=round(labels); % Because there could be some noise in the analogical signal.
    Fs = 20000;
    
    figure()
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
    
    figure()
    signalOfInterest = signal(rows_act);
    notOfInterest    = signal(rows_rest);
    h = spectrum.welch; % creates the Welch spectrum estimator
    SOIf=psd(h,signalOfInterest,'Fs',Fs); % calculates and plot the one sided PSD
    NOIf=psd(h,notOfInterest,'Fs',Fs);
    plot(SOIf) % Plot the one-sided PSD.
    hold on
    plot(NOIf)
    title('Spectrum of SOI and NOI for the stimuli application')
    legend('SOI (upper curve)','NOI (lower curve)')

    ciao = SOIf.copy();
    
    figure()
    plot(SOIf.getdata() - NOIf.getdata());
    xlim([0 130])
    xticks(0:13:130)
    xticklabels(xticks/13)
    ylabel('Power/frequency (dB/Hz)')
    xlabel('Frequency (kHz)')
    title('Difference spectrum of SOI and NOI for the stimuli application');
    grid on
    
    %% Filtering
    
    %Determine the Low pass cutoff frequency for your filter
    Fl=800;
    %Determine the High pass cutoff frequency for your filter
    Fh=2200;
    
    
    [Signal_filtered]=filtra(signal',Fs,Fl,Fh);
    
    figure()
    soi_f=Signal_filtered(rows_act);
    noi_f=Signal_filtered(rows_rest);
    h1 = spectrum.welch; % creates the Welch spectrum estimator
    SOI_f=psd(h,soi_f,'Fs',Fs); % calculates and plot the one sided PSD
    NOI_f=psd(h,noi_f,'Fs',Fs);
    plot(SOI_f); % Plot the one-sided PSD.
    hold on;
    plot(NOI_f);
    title('Spectrum of SOI and NOI after filtering')
    legend('SOI (upper curve)','NOI (lower curve)')

    %% Feature extraction
    
    % Determine the step size
    step=Fs*0.1;
    
    % Compute the MAV and Zero-Crossing features for each time window
    % Hint: MAV and ZeroCross can be functions in separate files
    
    
    L2=length(signal);
    
    for i=1:step:(L2-step)
        % Compute MAV and ZCross features and store them in a 2D
        % "features" vector
        MAVI((i-1)/step +1) = MAV(signal(i:i+step));
        Wavelengthi((i-1)/step +1) = Wavelength(signal(i:i+step));
        ZCrossi((i-1)/step + 1) =  ZCross(signal(i:i+step));
    end
    
    % Resizing the label vector to the size of the features vector
    for i=1:step:(L2-step)
        z=ceil(i/step);
        labels_resized(z)=labels(i);
    end
    
    % plot the zscored features along with the corresponding labels
    
    zMAVI =  zscore(MAVI);
    zWavelengthi = zscore(Wavelengthi);
    zZCrossi  = zscore(ZCrossi);
    
    figure
    plot(zMAVI)
    hold on
    plot(zWavelengthi)
    plot(zZCrossi)

    stairs(labels_resized,'LineWidth',1.2)
    grid on
    grid minor
    xlabel('ROW')
    ylabel('Normalized feature value')
    legend('MAV','Wavelength','ZCross','Labels')
    title(['Normalized features for ',f_map(filenames{k})])
    
    %% Signal-to-Noise Estimation
    
    % Extract the features for the non zero labels
    
    MAVI_SOI = MAVI(labels_resized  ~= 0);
    Wavelengthi_SOI = Wavelengthi(labels_resized ~= 0);
    ZCrossi_SOI = ZCrossi(labels_resized ~= 0);
    
    % Compute Activation
    activation = signal(labels ~=0);
    % Extract the features for the labels = zero
    
    MAVI_NOI = MAVI(labels_resized  == 0);
    Wavelengthi_NOI = Wavelengthi(labels_resized == 0);
    ZCrossi_NOI = ZCrossi(labels_resized == 0);
    
    % Compute the noise, SNR and SNR in dB
    noise = signal(labels == 0);
    
    SNR_signal = mean(activation)/mean(noise);
    SNR_MAVI = mean(MAVI_SOI)/mean(MAVI_NOI);
    SNR_Wavelengthi = mean(Wavelengthi_SOI)/mean(Wavelengthi_NOI);
    SNR_ZCrossi = mean(ZCrossi_SOI)/mean(ZCrossi_NOI);
    
    SNR_signal_dB(k) = 20 * log10 (SNR_signal);
    SNR_MAVI_dB(k)= 20 * log10(SNR_MAVI);
    SNR_Wavelengthi_dB(k) = 20 * log10(SNR_Wavelengthi);
    SNR_ZCrossi_dB(k) = 20*log10(SNR_ZCrossi);
    
end
