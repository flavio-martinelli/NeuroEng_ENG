function [spikes, spikesIndex]=spike_extract(neuralData,threshold,Fs,timeWindow)

% threshold: number of std as a threshold, 
% Fs: sampling frequency in Hz, 
% timeWindow: Time window for spike observation in seconds 


%% Thresholding
artifactThreshold = 40;
neuralData=zscore(neuralData);
% Artifact removal
neuralData(find(abs(neuralData)>artifactThreshold))=0;
neuralData=zscore(neuralData);

th=threshold*std(neuralData); % the threshold is given as the number of the standard deviations. 
m=mean(neuralData);
n=round(Fs*timeWindow); % the length of the observation window for the spike

% Spike detection
spikes=zeros(1,n);
k=0;
i=1;
while(i<=(length(neuralData)-round(n/2)))
    if(neuralData(i)<=m-th )
        k = k+1;
        spikes(k,:)=neuralData(round((i-n/2)):round((i+n/2-1)));
        spikesIndex(k)=i;
        i=i+round(n); % with this we avoid that two spikes overlap.    
    else
        i=i+1;
    end
end

% Spikes with excessive amplitude are removed
a=1;
for(i=1:k)
    if(isempty(find(abs(spikes(i,:))>m+2*th)))
        filteredSpikes(a,:)=spikes(i,:);
        a=a+1;
    end
end
k=a-1;


end