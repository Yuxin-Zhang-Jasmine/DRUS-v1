function [data, decimated_time_vector] = demodulate(time_vector,input_data,modulation_frequency,downsample_factor)

    %-- Function which performs the demodulation of RF ultrasound signals
    %-- Function prototype: [data, decimated_time_vector] = demodulate(time_vector,input_data,modulation_frequency,downsample_factor)
    
    %-- Authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %--          Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)

    %-- $Date: 2016/03/01 $  
    

    %-- complex envelope 
    data = hilbert(input_data(:,:)); % Find pre-envelope
    
    %-- demodulate
    demodVect = exp(-1i*2*pi*modulation_frequency*time_vector);
    
    %-- reshape
    siz = size(input_data);
    data = reshape(bsxfun(@times, data, demodVect), siz);

    %-- downsampling
    assert(downsample_factor==round(downsample_factor),'The downsample factor must be an integer');
    deciVect = 1:downsample_factor:siz(1);              %-- decimating vector
    decimated_time_vector = time_vector(deciVect);      %-- resampled time vector
    data = data(deciVect,:,:);                          %-- resampled data

end








