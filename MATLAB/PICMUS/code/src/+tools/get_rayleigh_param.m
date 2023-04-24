function sig = get_rayleigh_param(sample)

    %-- Function which estimates the parameter of the Rayleigh distribution
    %-- Function prototype: sig = get_rayleigh_param(sample)
    
    %-- Authors: Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)
    %--          Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)

    %-- $Date: 2016/03/01 $  

    sample = sample(:);
    M = sum( power(sample,2) );
    N = length(sample);
    sig = M / (2*N);

end