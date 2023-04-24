function [filtered_p,h,w] = low_pass(p,Fs,F)

    %-- Function which implements a low-pass filtering with Kaiser FIR filter
    %-- Function prototype: filtered_p = low_pass(p,Fs,[upper_freq_on upper_freq_off])
    
    %-- Authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %--          Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)

    %-- $Date: 2016/03/01 $  
    
    %-- filter specification
    A = [1 0];                                          %-- band type: 0='stop', 1='pass'
    dev = [1e-3 1e-3];                                  %-- ripple/attenuation spec
    [M,Wn,beta,typ] = kaiserord(F,A,dev,Fs);            %-- window parameters
    b = fir1(M,Wn,typ,kaiser(M+1,beta),'noscale');      %-- filter design
    [h,w] = freqz(b);
    
    %-- filtering
    filt_delay = round((length(b)-1)/2);
    filtered_p = filter(b,1,[p; zeros(filt_delay,size(p,2),size(p,3),size(p,4))],[],1);

    %-- correcting the delay
    filtered_p = filtered_p((filt_delay+1):end,:,:);
    
end
