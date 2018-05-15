function [Wavelengthi] = Wavelength(signal)
    
    Wavelengthi = (diff(signal),1)/length(signal);
    
end