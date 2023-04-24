function apod = apodization(distance,aperture,window)

    %-- Function which assigns different apodization to a set of pixels and elements
    %-- Function prototype: apodization = apodization(distance,aperture,window)
    
    %-- Authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %--          Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)

    %-- $Date: 2016/03/01 $  

    switch(window)
        case 'none' 
            apod = ones(size(distance)); 
        case 'boxcar' 
            apod = double(distance<=aperture/2); 
        case 'hanning'
            apod = double(distance<=aperture/2).*(0.5 + 0.5*cos(2*pi*distance./aperture)); 
        case 'hamming'
            apod = double(distance<=aperture/2).*(0.53836 + 0.46164*cos(2*pi*distance./aperture)); 
        case 'tukey25'
            roll=0.25;
            apod =(distance<(aperture/2*(1-roll))) + (distance>(aperture/2*(1-roll))).*(distance<(aperture/2)).*0.5.*(1+cos(2*pi/roll*(distance./aperture-roll/2-1/2)));                               
        case 'tukey50'
            roll=0.5;
            apod=(distance<(aperture/2*(1-roll))) + (distance>(aperture/2*(1-roll))).*(distance<(aperture/2)).*0.5.*(1+cos(2*pi/roll*(distance./aperture-roll/2-1/2)));                               
        case 'tukey75'
            roll=0.75;
            apod=(distance<(aperture/2*(1-roll))) + (distance>(aperture/2*(1-roll))).*(distance<(aperture/2)).*0.5.*(1+cos(2*pi/roll*(distance./aperture-roll/2-1/2)));                               
        otherwise
            error('Unknown window type. Known types are: boxcar, hamming, hanning, tukey25, tukey50, tukey75.');
    end
    
end
