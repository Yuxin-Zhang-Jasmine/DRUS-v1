classdef us_resolution < handle
    
    %-- Class defining a testing procedure to assess resolution performance of 
    %-- beamforming techniques in ultrasound imaging
        
    %-- Authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %--          Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)

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
        
    end
    
    %-- data
    properties (SetAccess = private)
        padROI
    end
    
    
    %------------------------------------------------------
    %------------------------------------------------------
    %-- Class methods
    
    %-- Constructor
    methods (Access = public)
        
        function h = us_resolution(input_name)
            
            if exist('input_name') 
                h.name= input_name; 
            else
                h.name = ' ';
            end            
            h.creation_date=sprintf('%d/%02d/%d %02d:%02d:%02.2f',clock);
            h.author = ' ';
            h.affiliation = ' ';
            h.dynamic_range = 60;
            h.padROI = 1.8e-3;
            h.flagDisplay = 0;
            h.score = [];
            
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
            h.score = zeros(nb_frames,size(h.pht.sca,1),2);
            maskROI = zeros(size(h.scan.x_matrix));
            for k=1:size(h.pht.sca,1)                
                %-- Compute mask inside
                x = h.pht.sca(k,1);
                z = h.pht.sca(k,3);                
                %-- Compute mask ROI
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
                    close all;
                    hid = figure(1); subplot(1,3,1); set(gca,'fontsize',16);
                    imagesc((h.scan.x_axis)*1e3,(h.scan.z_axis)*1e3,bmode); 
                    shading flat; colormap gray; caxis(vrange); colorbar; hold on;
                    axis equal manual; xlabel('x [mm]'); ylabel('z [mm]'); axis([x_lim z_lim]);
                    set(gca,'YDir','reverse');
                    set(gca,'fontsize',16);                
                    title(sprintf('%s\n %d plane waves',char(h.image.name),h.image.number_plane_waves(f)));
                    set(hid,'position',[124 175 1142 413]);
                    hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,maskROI,[1 1],'b-');
                end

                %-- Perform  resolution measurements
                for k=1:size(h.pht.sca,1)

                    %-- Concentrate on region of interest
                    patchImg = bmode;
                    patchImg(maskROI~=k) = min(bmode(:));
                    patchMask = maskROI;
                    patchMask(maskROI~=k) = 0;

                    %-- Extract region of interest
                    [idzz,idxx] = find(patchMask==k);
                    x_lim_patch = [h.scan.x_axis(min(idxx)) h.scan.x_axis(max(idxx))]*1e3;
                    z_lim_patch = [h.scan.z_axis(min(idzz)) h.scan.z_axis(max(idzz))]*1e3;
                    x_patch = h.scan.x_axis(min(idxx):1:max(idxx))*1e3;
                    z_patch = h.scan.z_axis(min(idzz):1:max(idzz))*1e3;

                    %-- Extract maximum point coordinates
                    [idz,idx] = find(patchImg==max(patchImg(:)));
                    signalLateral = patchImg(idz,min(idxx):max(idxx));
                    signalAxial = patchImg(min(idzz):max(idzz),idx);

                    %-- Display intermediate testing
                    if (h.flagDisplay==1)
                        pause(1);
                        %-- Center display on current patch
                        figure(1); subplot(1,3,1); 
                        hold on; plot([x_lim_patch(1) x_lim_patch(end)],[h.scan.z_axis(idz) h.scan.z_axis(idz)]*1e3,'-b','linewidth',1);
                        hold on; plot([h.scan.x_axis(idx) h.scan.x_axis(idx)]*1e3,[z_lim_patch(1) z_lim_patch(end)],'-r','linewidth',1);
                        axis([x_lim_patch z_lim_patch]);
                    end    

                    %-- Compute resolutions
                    [res_axial] = h.Compute_6dB_Resolution(z_patch,signalAxial,2,'-r');                    
                    [res_lateral] = h.Compute_6dB_Resolution(x_patch,signalLateral,3,'-b');

                    %-- Store score
                    h.score(f,k,1) = res_axial;
                    h.score(f,k,2) = res_lateral;

                end
                
            end
        end       
    end
    
    methods (Access = private)
        
        %-- Compute resolution from -6dB definition
        function [res] = Compute_6dB_Resolution(h,x_axis,y_signal,num,color)   
            
            %-- Perform interpolation
            coeff = 10;
            nb_sample = length(x_axis);
            nb_interp = nb_sample * coeff;
            x_interp = linspace(x_axis(1),x_axis(end),nb_interp);
            y_interp = interp1(x_axis,y_signal,x_interp);
            
            ind = find(y_interp >= (max(y_interp)-6) );
            idx1 = min(ind);
            idx2 = max(ind);
            res = x_interp(idx2) - x_interp(idx1);            
            
            %-- Display profil
            if (h.flagDisplay==1) 
                figure(1); subplot(1,3,num);
                plot(x_interp,y_interp,color,'linewidth',2);    
                hold on; plot([x_interp(idx1) x_interp(idx1)],[-100 0],'-k','linewidth',1);
                hold on; plot([x_interp(idx2) x_interp(idx2)],[-100 0],'-k','linewidth',1);            
                hold off; ylabel('Amp [dB]');
                if (num==2)
                    title(sprintf('Axial res = %02.2f [mm]',res));
                    xlabel('z [mm]');
                else
                     xlabel('x [mm]');
                    title(sprintf('Lateral res = %02.2f [mm]',res));
                end
            end
            
        end
        
    end
    
    
    %------------------------------------------------------
    %-- Get methods

    %-- Methods call to get attribute
    methods (Access = public)
       
        function data = getNumberPlaneWavesList(h)            
            data = h.image.number_plane_waves;
        end
        
    end
    
end
