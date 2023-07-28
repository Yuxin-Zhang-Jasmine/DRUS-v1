function [FWHMA,FWHML, CNR, gCNR, SNR, KS] = evaluation(path_phantom,image,flagDisplay)
    %-- Function used for the evaluation in ultrasound imaging
    
    %-- Authors: yuxin Zhang

    %-- $Date: 2023/June/07 $  

    %-- Initialization
    FWHMA = 0; FWHML = 0; CNR = 0; gCNR = 0; SNR = 0; KS = -100;

    %-- Read input data
    pht = us_phantom();
    pht.read_file(path_phantom);

    if pht.N_scatterers < 30
        close all
        %-- Perform testing for FWHM
        FWHM = us_resolution();
        FWHM.pht = pht;
        FWHM.scan = image.scan;
        FWHM.image = image;
        FWHM.flagDisplay = flagDisplay;
        FWHM.evaluate();
        FWHM = FWHM.score;
        FWHMA = mean(FWHM(:,:,1));
        FWHML = mean(FWHM(:,:,2));
    end

    if pht.occlusionDiameter ~= 0
        close all
        %-- Perform testing for CNR
        CNR = us_contrast();
        CNR.pht = pht;
        CNR.scan = image.scan;
        CNR.image = image;
        CNR.flagDisplay = flagDisplay;
        CNR.evaluate();
%         %--keep all of the values
%         [M, N] = size(CNR.score);
%         CNRres = cell(M,N);
%         for i = 1: M*N
%             CNRres{i} = CNR.score(i);
%         end
%         CNR = CNRres;
        %--keep only the average 
        CNR = mean(CNR.score(:));

        close all
        %-- Perform testing for gCNR
        gCNR = us_gcnr();
        gCNR.pht = pht;
        gCNR.scan = image.scan;
        gCNR.image = image;
        gCNR.flagDisplay = flagDisplay;
        gCNR.evaluate();
%         %--keep all of the values
%         [M, N] = size(gCNR.score);
%         gCNRres = cell(M,N);
%         for i = 1: M*N
%             gCNRres{i} = gCNR.score(i);
%         end
%         gCNR = gCNRres;
        %--keep only the average        
        gCNR = mean(gCNR.score(:));   
    end
        
    %-- Perform testing for speckle quality
    if pht.RoiPsfTimeX ~=0
        close all
        KS = us_speckle_quality();
        KS.pht = pht;
        KS.scan = image.scan;
        KS.image = image;
        KS.flagDisplay = flagDisplay;
        KS.evaluate(); 
%         %--keep all of the values
%         [M, N] = size(KS.score);
%         KSres = cell(M,N);
%         for i = 1: M*N
%             KSres{i} = KS.score(i);
%         end
%         KS = KSres;
        %--keep only the average
        KS = mean(KS.score(:));

        close all    
        SNR = us_speckle_SNR();
        SNR.pht = pht;
        SNR.scan = image.scan;
        SNR.image = image;
        SNR.flagDisplay = flagDisplay;
        SNR.evaluate(); 
%         %--keep all of the values
%         [M, N] = size(SNR.score);
%         SNRres = cell(M,N);
%         for i = 1: M*N
%             SNRres{i} = SNR.score(i);
%         end
%         SNR = SNRres;  
        %--keep only the average
        SNR = mean(SNR.score(:));
    end


    
end



    