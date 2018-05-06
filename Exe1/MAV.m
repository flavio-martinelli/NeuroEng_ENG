function [Mav] = MAV(signal)
    
    Mav = norm(signal,1)/length(signal);
    
end