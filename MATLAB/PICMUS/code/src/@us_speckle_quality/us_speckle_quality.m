classdef us_speckle_quality < handle
    
    %-- Class defining a testing procedure to assess contrast performance of 
    %-- beamforming techniques in ultrasound imaging
    
    %-- Authors: Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)
    %--          Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %--          Adrien Besson (adrien.besson@epfl.ch)

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
    
    
    %------------------------------------------------------
    %------------------------------------------------------
    %-- Class methods
    
    %-- Constructor
    methods (Access = public)
        
        function h = us_speckle_quality(input_name)
            if exist('input_name') 
                h.name= input_name; 
            else
                h.name = ' ';
            end            
            h.creation_date=sprintf('%d/%02d/%d %02d:%02d:%02.2f',clock);
            h.author = ' ';
            h.affiliation = ' ';
            h.dynamic_range = 60;
            h.score = 0;
            h.flagDisplay = 0;
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
            h.score = zeros(nb_frames,length(h.pht.RoiCenterX));

            %-- Setting axis limits (mm)
            x_lim = [min(h.scan.x_matrix(:)) max(h.scan.x_matrix(:))]*1e3; 
            z_lim = [min(h.scan.z_matrix(:)) max(h.scan.z_matrix(:))]*1e3; 
            
            %-- Setting dynamic range for visualization
            vrange = [-h.dynamic_range 0];            
            
            %-- padding variable
            padROIx = h.pht.RoiPsfTimeX * h.pht.lateralResolution;
            padROIz = h.pht.RoiPsfTimeZ * h.pht.axialResolution;            
            
            %-- Loop over frames
            for f=frame_list            
            
                %-- Compute dB values
                env = h.image.data(:,:,f);
                bmode = 20*log10(env./max(env(:)));                
                %yuxin edit
                [bmode,~] = envelope(bmode,1,'peak');
                
                %-- Ploting image reconstruction
                if (h.flagDisplay==1)
                    figure(1); set(gca,'fontsize',16);
                    imagesc((h.scan.x_axis)*1e3,(h.scan.z_axis)*1e3,bmode); 
                    shading flat; colormap gray; caxis(vrange); colorbar; hold on;
                    axis equal manual; xlabel('x [mm]'); ylabel('z [mm]'); axis([x_lim z_lim]);
                    set(gca,'YDir','reverse');
                    set(gca,'fontsize',16);
                    title(sprintf('%s\n %d plane waves',char(h.image.name),h.image.number_plane_waves(f)));
                    pause(0.1);
                end

                %-- Main loop
                for k=1:length(h.pht.RoiCenterX)

                    %-- Compute mask inside
                    x = h.pht.RoiCenterX(k);
                    z = h.pht.RoiCenterZ(k);
                    
                    %-- Compute mask ROI
                    maskROI = k * ( (h.scan.x_matrix > (x-padROIx(k))) & ...
                                 (h.scan.x_matrix < (x+padROIx(k))) & ...
                                 (h.scan.z_matrix > (z-padROIz(k))) & ...
                                 (h.scan.z_matrix < (z+padROIz(k))) );

                    %-- Extract corresponding enveloppe ROI
                    [idzz,idxx] = find(maskROI==k);
                    envRoi = env(min(idzz):max(idzz),min(idxx):max(idxx));
                    
                    %-- Downsample the block to ensure statistical independency
                    sample = envRoi(1:5:end,1:5:end);
                    sample = sample(:);         
                    
                    %-- Empirical evaluation of the variance on the selected block
                    var_block = tools.get_rayleigh_param(sample);
                    
                    %-- Application of the KS test against the Rayleigh pdf with 5% confidence interval
                    testres = 1-kstest(sample, 'CDF', [sample, raylcdf(sample, sqrt(var_block))], 'alpha', 0.05);
                    
                    %-- Ploting Roi contour along with a code color
                    %-- green: test passed || red: test failed                    
                    if (h.flagDisplay==1)
                        figure(1);
                        if (testres==0)
                            hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,maskROI,[1 1],'r-','linewidth',2);
                        else
                            hold on; contour(h.scan.x_axis*1e3,h.scan.z_axis*1e3,maskROI,[1 1],'g-','linewidth',2);
                        end
                        pause(0.1);
                    end
                    h.score(f,k) = testres;

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
