function [fx, F] = power_spectrum(data,fs)

    %-- Function which computes the power spectrum of ultrasound signal
    %-- Function prototype: [fx,F] = power_spectrum(data,fs)

    %-- Authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %--          Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)

    %-- $Date: 2016/03/01 $          
    
    Nfft = size(data,1);
    F = fftshift(fft(data, Nfft));
    F = abs(squeeze(mean(mean(mean(abs(F(:,:,:,:)),4),3),2)));
    F = F/max(F);
    fx = linspace(-fs/2,fs/2, Nfft);
    
end
