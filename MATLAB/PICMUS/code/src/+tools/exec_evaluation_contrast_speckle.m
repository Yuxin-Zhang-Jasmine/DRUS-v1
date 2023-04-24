function exec_evaluation_contrast_speckle(path_scan,path_phantom,path_img,flag_simu,flag_display,path_output)
  
    %-- Function used for the evaluation of contrast and speckle quality in ultrasound imaging
    
    %-- Authors: Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)
    %--          Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)

    %-- $Date: 2016/03/01 $  


    %-- Parameters
    category = [1 11 75];


    %-- Convert input argument received as string
    flag_simu = uint8(str2num(flag_simu));   		%-- convert string back to int
    flag_display = uint8(str2num(flag_display));   	%-- convert string back to int

    
    %-- Read input data
    scan = linear_scan();
    scan.read_file(path_scan);
    pht = us_phantom();
    pht.read_file(path_phantom);
    image = us_image();
    image.read_file(path_img);


    %-- Perform testing for contrast
    testing_contrast = us_contrast();
    testing_contrast.pht = pht;
    testing_contrast.scan = scan;
    testing_contrast.image = image;
    testing_contrast.flagDisplay = flag_display;
    testing_contrast.evaluate();
    
    
    %-- Perform testing for speckle quality
    testing_speckle = us_speckle_quality();
    testing_speckle.pht = pht;
    testing_speckle.scan = scan;
    testing_speckle.image = image;
    testing_speckle.flagDisplay = flag_display;
    testing_speckle.evaluate();    
    
    
    %-- Check categories
    pw_list = testing_contrast.getNumberPlaneWavesList();
    indexCategory = [0 0 0];
    selected_pw = [];    
    
    
    %-- Extract which category among the 3 main ones
    its = 0;
    for c=1:length(category)
        idx = find(pw_list==category(c));
        if (~isempty(idx))
            indexCategory(c) = idx;
            its = its+1;
            selected_pw(its) = category(c);
        end
    end    
    
    
    %-- Extract bonus category
    pw_bonus = setdiff(pw_list,selected_pw);
    if ( (length(pw_bonus)>1) || isempty(pw_bonus) )
        pw_bonus = 0;
    end
        
    
    %-- Get output result and store it in a file
    fid = fopen(path_output,'w');    
    
    %-- Store main category results
    for c=1:length(category)
        %-- Write result
        fprintf(fid,'#C%d \n\n',category(c));
        if ( indexCategory(c)~=0 )
            idx = indexCategory(c);
            score_contrast = squeeze(testing_contrast.score(idx,:,:));
            writeScoreContrast(fid,score_contrast,flag_simu);
            score_speckle = squeeze(testing_speckle.score(idx,:));
            writeScoreSpeckle(fid,score_speckle);
        end            
        fprintf(fid,'\n#END\n\n');            
    end    
    
    %-- Store bonus category results        
    if (pw_bonus == 0)            
        fprintf(fid,'#BONUS-NONE \n\n');
    else
        fprintf(fid,'#BONUS-C%d \n\n',pw_bonus);
        idx = find(pw_list==pw_bonus);
        score = squeeze(testing_contrast.score(idx(1),:,:));
        writeScoreContrast(fid,score,flag_simu);
        score_speckle = squeeze(testing_speckle.score(idx(1),:));
        writeScoreSpeckle(fid,score_speckle);        
    end
    fprintf(fid,'\n#END\n\n');

    fclose(fid);        
    
end
    
    
    
%----------------------------------------------------------------
%----------------------------------------------------------------
%----------------------------------------------------------------
%-- Auxillary functions

%-- Write contrast results 
function writeScoreContrast(fid,score,flag_simu)

    if (flag_simu==1)
        writeSimulationScoreContrast(fid,score);
    else
        writeExperimentScoreContrast(fid,score);
    end
    
end


%-- Write contrast results corresponding to numerical phantom
function writeSimulationScoreContrast(fid,score)
        
    %-- Mean axial and lateral values
    fprintf(fid,'Mean contrast scores (dB): %.2f \n\n',mean(score));
    
%     fprintf(fid,'Mean contrast scores (dB) \n');
%     fprintf(fid,'%.2f \n\n',mean(score));  

    %-- Left column
    fprintf(fid,'Left column targets \n');
    for k=4:6
        fprintf(fid,'%.2f\n',score(k));  
    end
    fprintf(fid,'\n');

    %-- Middle column
    fprintf(fid,'Middle column targets \n');
    for k=1:3
        fprintf(fid,'%.2f\n',score(k));  
    end
    fprintf(fid,'\n');

    %-- Right column
    fprintf(fid,'Right column targets \n');
    for k=7:9
        fprintf(fid,'%.2f\n',score(k));  
    end
    fprintf(fid,'\n');
           
end


%-- Write contrast results corresponding to ex-vivo phantom
function writeExperimentScoreContrast(fid,score)

    %-- Mean axial and lateral values
    fprintf(fid,'Mean contrast scores (dB): %.2f \n\n',mean(score)); 
    
%     fprintf(fid,'Mean contrast scores (dB) \n');
%     fprintf(fid,'%.2f\n\n',mean(score));  

    %-- Middle column
    fprintf(fid,'Middle column targets \n');
    for k=1:2
        fprintf(fid,'%.2f\n',score(k));  
    end
    fprintf(fid,'\n');  

end


%-- Write speckle results 
function writeScoreSpeckle(fid,score)

    %-- Speckle score display (0: test failed || 1: test passed)
    fprintf(fid,'PENALITY: ');
    if ~isempty(find(score==0,1))
        fprintf(fid,'-40\n');
    else
        fprintf(fid,'0\n');
    end
    fprintf(fid,'Speckle quality (from Kolmogorov-Smirnov test) \n');    
    
end

