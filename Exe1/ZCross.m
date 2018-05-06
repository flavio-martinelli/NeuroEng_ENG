function [ZCROSS] = ZCross(signal)
    
    ZCROSS = norm(diff(signal),1)/length(signal);
    
end