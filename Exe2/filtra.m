function [signal]=filtra(signal,Fs,Fb,Fh)
%The function implements:high pass, low pass// signal, Fs= sampling Hz, Fb=freq low pass, Fh=freq high pass

Fc=2*Fb/Fs;
Fch=2*Fh/Fs;
order=3;
[nChannel N]=size(signal);

%% HIGH-PASS
d=fdesign.highpass('N,Fc',order,Fch);
h=design(d,'butter');

for(i=1:nChannel)
signal(i,:)=filtfilt(h.sosMatrix,h.ScaleValues,signal(i,:));

end

%% LOW-PASS
d=fdesign.lowpass('N,Fc',order,Fc);
h=design(d,'butter');

for(i=1:nChannel)
signal(i,:)=filtfilt(h.sosMatrix,h.ScaleValues,signal(i,:));
end

end