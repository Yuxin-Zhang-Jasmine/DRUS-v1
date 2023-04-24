function exec_evaluation_resolution_distorsion(path_scan,path_phantom,path_img,flag_simu,flag_display,path_output)
  
    %-- Function used for the evaluation of resolution in ultrasound imaging
    
    %-- Authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %--          Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)

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


    %-- Perform testing for resolution
    testing_resolution = us_resolution();
    testing_resolution.pht = pht;
    testing_resolution.scan = scan;
    testing_resolution.image = image;
    testing_resolution.flagDisplay = flag_display;
    testing_resolution.evaluate();
    
    
    %-- Perform testing for distorsion
    if (flag_simu==1)
        testing_distorsion = us_distorsion();
        testing_distorsion.pht = pht;
        testing_distorsion.scan = scan;
        testing_distorsion.image = image;
        testing_distorsion.flagDisplay = flag_display;
        testing_distorsion.evaluate();
    end
    
    
    %-- Check categories
    pw_list = testing_resolution.getNumberPlaneWavesList();
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
            score_resolution = squeeze(testing_resolution.score(idx,:,:));
            writeScoreResolution(fid,score_resolution,flag_simu);
            if (flag_simu==1)
                score_distorsion = squeeze(testing_distorsion.score(idx,:));
            else
                score_distorsion = nan;
            end
            writeScoreDistorsion(fid,score_distorsion,flag_simu);
        end
        fprintf(fid,'\n#END\n\n');
    end

    %-- Store bonus category results        
    if (pw_bonus == 0)            
        fprintf(fid,'#BONUS-NONE \n\n');
    else
        fprintf(fid,'#BONUS-C%d \n\n',pw_bonus);
        idx = find(pw_list==pw_bonus);
        score_resolution = squeeze(testing_resolution.score(idx(1),:,:));
        writeScoreResolution(fid,score_resolution,flag_simu);
        if (flag_simu==1)
            score_distorsion = squeeze(testing_distorsion.score(idx,:));
        else
            score_distorsion = nan;
        end
        writeScoreDistorsion(fid,score_distorsion,flag_simu);        
    end
    fprintf(fid,'\n#END\n\n');
    
    fclose(fid);       
   
end



%----------------------------------------------------------------
%----------------------------------------------------------------
%----------------------------------------------------------------
%-- Auxillary functions

function writeScoreResolution(fid,score,flag_simu)

    if (flag_simu==1)
        writeSimulationScoreResolution(fid,score);
    else
        writeExperimentScoreResolution(fid,score);
    end
    
end


%-- Write results corresponding to numerical phantom
function writeSimulationScoreResolution(fid,score)

    %-- Mean axial and lateral values
    fprintf(fid,'Mean resolution scores (axial & lateral): %.2f \t %.2f\n\n',mean(score(:,1)),mean(score(:,2)));

    %-- Vertical target scores
    fprintf(fid,'Vertical targets (axial & lateral) \n');
    for k=1:8
        fprintf(fid,'%.2f \t %.2f\n',score(k,1),score(k,2));  
    end
    fprintf(fid,'\n');

    %-- Horizontal 2-cm target scores
    fprintf(fid,'Horizontal targets at 2 cm (axial & lateral) \n');
    for k=9:11
        fprintf(fid,'%.2f \t %.2f\n',score(k,1),score(k,2));  
    end    
    k = 3;
    fprintf(fid,'%.2f \t %.2f\n',score(k,1),score(k,2)); 
    for k=12:14
        fprintf(fid,'%.2f \t %.2f\n',score(k,1),score(k,2));  
    end
    fprintf(fid,'\n');

    %-- Horizontal 4-cm target scores
    fprintf(fid,'Horizontal targets at 4 cm (axial & lateral) \n');
    for k=15:17
        fprintf(fid,'%.2f \t %.2f\n',score(k,1),score(k,2));  
    end    
    k = 7;
    fprintf(fid,'%.2f \t %.2f\n',score(k,1),score(k,2)); 
    for k=18:20
        fprintf(fid,'%.2f \t %.2f\n',score(k,1),score(k,2));  
    end
    fprintf(fid,'\n');

end


%-- Write results corresponding to ex-vivo phantom
function writeExperimentScoreResolution(fid,score)

    %-- Mean axial and lateral values
    fprintf(fid,'Mean resolution scores (axial & lateral): %.2f \t %.2f\n\n',mean(score(:,1)),mean(score(:,2)));
    
    %-- Vertical target scores
    fprintf(fid,'Vertical targets (axial & lateral) \n');
    for k=1:3
        fprintf(fid,'%.2f \t %.2f\n',score(k,1),score(k,2));  
    end
    fprintf(fid,'\n');

    %-- Horizontal 4-cm target scores
    fprintf(fid,'Horizontal targets near 4 cm (axial & lateral) \n');
    k=4;
    fprintf(fid,'%.2f \t %.2f\n',score(k,1),score(k,2));
    k=3;
    fprintf(fid,'%.2f \t %.2f\n',score(k,1),score(k,2));
    k=5;
    fprintf(fid,'%.2f \t %.2f\n\n',score(k,1),score(k,2));

end


%-- Write distorsion results 
function writeScoreDistorsion(fid,score,flag_simu)


    fprintf(fid,'PENALITY: ');
    if (flag_simu==0)
        fprintf(fid,'NONE\n');
    else
        if ~isempty(find(score==0,1))
            fprintf(fid,'-40\n');
        else
            fprintf(fid,'0\n');            
        end
    end
    fprintf(fid,'Distorsion \n');

end


