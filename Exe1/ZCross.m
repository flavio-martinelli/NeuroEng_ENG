function [ZCrossi] = ZCross(signal)
    ZCrossi = 0;
    for i = 1:length(signal)-1
        if signal(i+1)*signal(i)<=0
            ZCrossi = ZCrossi+1;
        end
    end
    ZCrossi = ZCrossi/length(signal);
end