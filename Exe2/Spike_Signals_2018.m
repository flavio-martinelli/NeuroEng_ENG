clear all
close all
neuralData=load ('Proprio_1.txt');
Fs=24414.06;

% Compute the time vector from the frequency value
time = (1:length(neuralData))/Fs;


%% Pre-processing of the signal 
figure()
plot(time,neuralData')
title('Original Signal'); 
xlabel('Time (s)')
ylabel('Potential (V)')

%low pass 
Fb=5000;
%High Pass
Fh=500;

neuralData = filtra(neuralData,Fs,Fb,Fh);

figure()
plot(time,neuralData')
title('Filtered Signal'); 
xlabel('Time (s)')
ylabel('Potential (V)')


%% Spike Extraction 

%%% Vary the threshold from 2 to 6  of the spike extracting algorithm
%%% (spike_extract.m)  and observe how this affects the detection of spike
%%% events and the computation time 
threshold = 2:0.5:6;
timeWindow = 0.002;

for i = 1 : length(threshold)
    init = tic;
    [filteredSpikes{i}, spikesIndex{i}]=spike_extract(neuralData,threshold(i),Fs,timeWindow);
    endt = toc(init);
    
    comptime(i) = endt;
    
    n_spikes(i) = length(filteredSpikes{i});
    time2 =( 1:size(filteredSpikes{i},2))*1000/Fs; %in ms
    
    figure()
    plot(time2,filteredSpikes{i}')
    xlim([0 2])
    xlabel('Time (ms)')
    ylabel('Potential')
    title(['Spikes with threshold ',num2str(threshold(i)),' standard deviations'])
    grid on
    grid minor
    
end

% Plot the waveform of the extracted spikes (filteredSpikes) for each
% threshold and choose the optimal threshold value

figure
plot(threshold,comptime)
xlabel('Threshold value (n. of std)')
ylabel('Computation time (s)')
grid on
grid minor


figure
plot(n_spikes,comptime)
xlabel('Number of spikes detected')
ylabel('Computation time (s)')
grid on
grid minor
%% Principal component analysis

% Extract the filteredSpikes using the optimal threshold (found previously)

threshold = 4;
[filteredSpikes, spikesIndex]=spike_extract(neuralData,threshold,Fs,timeWindow);

% Compute the principal components of the filteredSpikes (size n= number of spikes)
% Hint: use Matlab's function pca

[coeff,score] = pca(filteredSpikes);

% Store the representation of the filteredSpikes on the 2 first PC -  in a nx2 matrix 

pc2 = score(:,1:2);


%% Spike sorting

% Cluster the spikes into 2 groups using the kmeans algorithm (already
% implemented in matlab) and store the cluster indices in a vector "idx".
[idx,C] = kmeans(pc2,2);
% Plot the two clusters (two distinct colours) of data and the centroids of each cluster 
figure
scatter(pc2(idx==1,1),pc2(idx==1,2),'ro')
hold on
scatter(pc2(idx ==2,1),pc2(idx==2,2),'bo')
scatter(C(1,1),C(1,2),100,'k*')
scatter(C(2,1),C(2,2),100,'g*')
title('Distribution of the data along the two PCs')
xlabel('PC1')
ylabel('PC2')
grid on
grid minor
legend('Cluster 1','Cluster 2')

% Plot the  waveforms contained in the two clusters in two distinct colours 

figure
plot(time2,filteredSpikes(idx==1,:),'b')
hold on
plot(time2,filteredSpikes(idx==2,:),'r')
xlim([0,2])
xlabel('Time (ms)')
ylabel('Potential')
title('Spikes sorted along clusters')
grid on
grid minor
legend('Cluster 1','Cluster 2')



%% Calculation of the Firing Rate 
% Compute the firing rate using a time window of 200 ms (and the optimal
% threshold found previously)
threshold = 4;
shift = round(Fs*0.2);

[spikes,spikeidx] = spike_extract(neuralData,threshold,Fs,timeWindow);

n_windows = ceil(length(neuralData)/shift);
fires = zeros(1,n_windows);

for i = 1:length(spikeidx)
     window_of_interest  = ceil(spikeidx(i)/shift);
     fires(window_of_interest) = fires(window_of_interest) + 1; 
end

firing_rate = fires / 0.2;
% Plot the original neural signal overlapped with the detected spikes
spikes_signal = zeros(1,length(neuralData));
% for i = 1:length(spikeidx)
%     spikes_signal(spikeidx(i):spikeidx(i)+48)=spikes(i,:);
% end

figure
plot(time,neuralData)
hold on
plot((spikeidx)/Fs, neuralData(spikeidx),'*')
plot(time,spikes_signal)
xlabel('Time (s)')
ylabel('Potential (mV)')
legend('Recording','Spike')

% Plot the original neural signal overlapped with the calculated firing
% rate
x_axis =  (1:n_windows)*0.2;

figure
yyaxis left
plot(time,neuralData)
hold on
yyaxis right
stairs(x_axis,firing_rate,'LineWidth',1.2)
xlabel('Time (s)')
ylabel('Potential (mV)')
ylim([-15 15])