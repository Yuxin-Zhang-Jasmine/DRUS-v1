classdef us_distorsion < handle
    
    %-- Class defining a testing procedure to assess distorsion artefacts
    %-- eventually produced by beamforming techniques in ultrasound imaging
    
    %-- Authors: Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)
    %--          Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)

    %-- $Date: 2016/03/01 $      
    
    
    %------------------------------------------------------
    %------------------------------------------------------
    %-- Class attributes
    
    properties (SetAccess = public)        
        %-- administration
        name                         %-- String containing the name of the us_field_simulation
        author                       %-- String containing the name of the author(s) of the beamformed image
        affiliation                  %-- String containing the affiliation of the author(s) of the beamformed image
        creation_date                %-- String containing the date the reconstruction was created

        %-- Common attributes
        pht
        scan
        image
        score
        flagDisplay
        dynamic_range
        z_correction
        
    end
    
    %-- data
    properties (SetAccess = private)
        padInside
        padROI 
        indiceToCheck
    end
    
    
    %------------------------------------------------------
    %------------------------------------------------------
    %-- Class methods
    
    %-- Constructor
    methods (Access = public)

        function h = us_distorsion(input_name)
            
            if exist('input_name') 
                h.name= input_name; 
            else
                h.name = ' ';
            end            
            h.creation_date = sprintf('%d/%02d/%d %02d:%02d:%02.2f',clock);
            h.author = ' ';
            h.affiliation = ' ';
            h.dynamic_range = 60; 
            h.padInside = 2.9570e-04;   %-- this value corresponds to lambda (c/f0)
            h.padROI = 1.8e-3;
            h.flagDisplay = 0;
            h.indiceToCheck = [1 5 8 9 14 15 20];
            h.z_correction = 0.2e-3;	%-- due to lens curvature and 0.5 dB / [MHz . cm] of the medium
        end
        
    end
    
    %-- set methods, input format check
    methods  
        
        %-- name
        function set.name(h,input)
            assert(isstr(input), 'Wrong format of the beamformed data name. It should be a string.');
            h.name = input;
        end
        %-- author
        function set.author(h,input)
            assert(isstr(input), 'Wrong format of the author. It should be a string.');
            h.author = input;
        end
        %-- affiliation
        function set.affiliation(h,input)
            assert(isstr(input), 'Wrong format of the affiliation. It should be a string.');
            h.affiliation = input;
        end
        %-- creation_date
        function set.creation_date(h,input_date)
            assert(isstr(input_date), 'Wrong format of the creation date. It should be a string.');
            h.creation_date = input_date;
        end
        %-- pht
        function set.pht(h,input)
            assert(isa(input,'us_phantom'), 'Wrong format of the phantom. It should be a US_PHANTOM class.');
            h.pht = input;
        end
        %-- scan
        function set.scan(h,input)
            assert(isa(input,'linear_scan'), 'Wrong format of the scan. It should be a LINEAR_SCAN class.');
            h.scan = input;
        end
        %-- image
        function set.image(h,input)
            assert(isa(input,'us_image'), 'Wrong format of the image. It should be a US_IMAGE class.');
            h.image = input;
        end
        %-- flag for display
        function set.flagDisplay(h,input)
            assert(numel(input)==1&&isnumeric(input), 'Wrong format of the flag for display. It should be a numeric scalar');
            h.flagDisplay = input;
        end
        %-- flag for display
        function set.dynamic_range(h,input)
            assert(numel(input)==1&&isnumeric(input), 'Wrong format of the dynamic range. It should be a numeric scalar');
            h.dynamic_range = input;
        end        

    end
        
    %-- Main method call to generate phantom
    methods (Access = public)
        
        function evaluate(h)
  
            %-- Define parameters / variables
            nb_frames = length(h.image.number_plane_waves);
            frame_list = 1:nb_frames;
            h.score = zeros(nb_frames,length(h.indiceToCheck));
            maskInside = zeros(size(h.scan.x_matrix));
            maskROI = zeros(size(h.scan.x_matrix));
            
            %-- Apply update on h.pht.sca -> z in order to take into
            %-- account the elevation focus + attenuation effects 
            %-- (in our case 0.2 mm effect)
            h.pht.sca(:,3) = h.pht.sca(:,3) + h.z_correction;
            
            for k=1:size(h.pht.sca,1)
                
                %-- Compute mask inside
                x = h.pht.sca(k,1);
                z = h.pht.sca(k,3);
                mask = k * ( (h.scan.x_matrix > (x-h.padInside)) & ...
                (h.scan.x_matrix < (x+h.padInside)) & ...
                (h.scan.z_matrix > (z-h.padInside)) & ...
                (h.scan.z_matrix < (z+h.padInside)) );
                maskInside = maskInside + mask;
                
                %-- Compute mask outside
                mask = k * ( (h.scan.x_matrix > (x-h.padROI)) & ...
                (h.scan.x_matrix < (x+h.padROI)) & ...
                (h.scan.z_matrix > (z-h.padROI)) & ...
                (h.scan.z_matrix < (z+h.padROI)) );
                maskROI = maskROI + mask;
                
            end

            %-- Setting axis limits (mm)
            x_lim = [min(h.scan.x_matrix(:)) max(h.scan.x_matrix(:))]*1e3; 
            z_lim = [min(h.scan.z_matrix(:)) max(h.scan.z_matrix(:))]*1e3;             
            
            %-- Setting dynamic range for visualization
            vrange = [-h.dynamic_range 0];            
            
            %-- Loop over frames
            for f=frame_list            
            
                %-- Compute dB values
                env = h.image.data(:,:,f);
                bmode = 20*log10(env./max(env(:)));

                %-- Ploting image reconstruction
                if (h.flagDisplay==1)
                    hid = figure; subplot(1,2,1); set(gca,'fontsize',16);
                    imagesc((h.scan.x_axis)*1e3,(h.scan.z_axis)*1e3,bmode); 
                    shading flat; colormap gray; caxis(vrange); colorbar; hold on;
                    axis equal manual; xlabel('x [mm]'); ylabel('z [mm]'); axis([x_lim z_lim]);
                    set(gca,'YDir','reverse');
                    set(gca,'fontsize',16);                           
                    title(sprintf('%s\n %d plane waves',char(h.image.name),h.image.number_plane_waves(f)));
                    set(hid,'position',[239 425 1142 413]);
                end

                %-- Perform distortion measurements
                ind = 0;
                for k=1:size(h.pht.sca,1)

                    if find(k==h.indiceToCheck)

                        ind = ind+1;
                        %-- Concentrate on region of interest
                        imTest = bmode;
                        imTest(maskROI~=k) = min(bmode(:));
                        maskTestROI = maskROI;
                        maskTestROI(maskROI~=k) = 0;
                        maskTestIn = maskInside;
                        maskTestIn(maskInside~=k) = 0;

                        %-- Perform test
                        [idz,idx] = find(imTest==max(imTest(:)));    
                        if ( maskTestIn(idz,idx) == k )
                            h.score(f,ind) = 1;
                        else
                            h.score(f,ind) = 0;
                        end    

                        %-- Display intermediate testing
                        if (h.flagDisplay==1)                

                            figure(hid); subplot(1,2,1); 
                            hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,maskTestROI,[1 1],'r-');  
                            hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,maskTestIn,[1 1],'g-');
                            [idzz,idxx] = find(maskTestROI==k);
                            x_lim_test = [h.scan.x_axis(min(idxx)) h.scan.x_axis(max(idxx))]*1e3;
                            z_lim_test = [h.scan.z_axis(min(idzz)) h.scan.z_axis(max(idzz))]*1e3;
                            figure(hid); subplot(1,2,2); imagesc(h.scan.x_axis*1e3,h.scan.z_axis*1e3,imTest); shading flat; colormap gray; caxis(vrange); colorbar;
                            axis equal manual; xlabel('x [mm]'); ylabel('z [mm]'); set(gca,'YDir','reverse'); set(gca,'fontsize',16); axis([x_lim_test z_lim_test]);
                            hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,maskTestIn,[1 1],'g-');
                            hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,maskTestROI,[1 1],'r-');
                            hold on; plot(h.scan.x_axis(idx)*1e3,h.scan.z_axis(idz)*1e3,'ob','linewidth',2);
                            if (h.score(f,ind) == 1)
                                title(['Succeeded | Pt = [',num2str(round(h.scan.x_axis(idx)*1e4)/10),' , ',num2str(round(h.scan.z_axis(idz)*1e4)/10),'] mm']);
                            else
                                title(['Failed | Pt = [',num2str(round(h.scan.x_axis(idx)*1e4)/10),' , ',num2str(round(h.scan.z_axis(idz)*1e4)/10),'] mm']);
                            end    
                            pause(1);
                        end
                    end 
                end                
            end                
        end        
    end    
end
