function env = envelope(p)

    %-- Function which computes the envelope of a modulated signal
    %-- Function prototype: env = envelope(p)
    
    %-- Authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %--          Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)

    %-- $Date: 2016/03/01 $  

    %-- complex demodulation
    env = reshape(abs(hilbert(p(:,:))),size(p));

end
