function [fc, bw] = estimate_frequency(data,Fs)

    %-- Function 
    
    %-- Authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %--          Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)

    %-- $Date: 2016/03/01 $  

    [fx] = tools.power_spectrum(data,Fs);
    fpw = filter(ones(1,26)./26,1,pw);fpw=[fpw(13:end); zeros(12,1)];
    [dc] = max(fpw.*(fx>0).'); fc=fx(ic);
    bw_up = min(fx((fx>fc)&(fpw<dc/2).'));        %-- -6dB upper limit
    bw_do = max(fx((fx<fc)&(fpw<dc/2).'));        %-- -6dB down limit
    fc = (bw_up+bw_do)/2;                         %-- center frequency
    bw = 2*(bw_up-fc);                            %-- bandwidth
    
end
