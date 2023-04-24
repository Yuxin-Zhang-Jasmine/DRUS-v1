classdef us_phantom < handle
    
    %-- Class defining a standard way to generate numerical phantom in
    %-- ultrasound
    
    %-- Authors: Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)
    %--          Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    
    
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
        N_scatterers
        amp
        sca
        xLimits
        zLimits
        bckDensity
                
        %-- Occlusion attributes
        occlusionCenterX
        occlusionCenterZ
        occlusionDiameter
        
        %-- Region to check speckle attributes
        RoiCenterX
        RoiCenterZ
        RoiPsfTimeX
        RoiPsfTimeZ
        
        %-- Point scatterers
        xPts
        zPts
                
        %-- system resolution (needed to fix the number of scatterers)
        axialResolution
        lateralResolution        
        
    end
    
    %-- data
    properties (SetAccess = private)
        mode
    end
    
    
    %------------------------------------------------------
    %------------------------------------------------------
    %-- Class methods
    
    %-- Constructor
    methods (Access = public)
        
        function h = us_phantom(input_name)
            %-- Constructor of the US_PHANTOM class.
            %-- Syntax:
            %-- US_PHANTOM(name) 
            %-- name: Name of the reconstruction
            
            if exist('input_name') 
                h.name= input_name; 
            else
                h.name = ' ';
            end            
            h.creation_date = sprintf('%d/%02d/%d %02d:%02d:%02.2f',clock);
            h.author = ' ';
            h.affiliation = ' ';
            h.xLimits = 0;
            h.zLimits = 0;
            h.bckDensity = 20;
            h.occlusionCenterX = 0;
            h.occlusionCenterZ = 0;
            h.occlusionDiameter = 0;
            h.RoiCenterX = 0;
            h.RoiCenterZ = 0;
            h.RoiPsfTimeX = 0;
            h.RoiPsfTimeZ = 0;
            
            h.mode = 0;
            h.xPts = 0;
            h.zPts = 20e-3;   
            h.axialResolution = 0.45e-3;
            h.lateralResolution = 0.45e-3;
            
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
        %-- xPts
        function set.xPts(h,input_xPts)
            assert(isnumeric(input_xPts), 'Wrong format of the x-axis point scatterers xPts. It should be a numeric vector in (m)');
            h.xPts = input_xPts;
        end
        %-- zPts
        function set.zPts(h,input_zPts)
            assert(isnumeric(input_zPts), 'Wrong format of the x-axis point scatterers zPts. It should be a numeric vector in (m)');
            h.zPts = input_zPts;
        end
        %-- xLimits
        function set.xLimits(h,input_xLimits)
            assert(isnumeric(input_xLimits), 'Wrong format of the x-axis limits xLimits. It should be a numeric vector in (m)');
            h.xLimits = input_xLimits;
        end 
        %-- zLimits
        function set.zLimits(h,input_zLimits)
            assert(isnumeric(input_zLimits), 'Wrong format of the z-axis limits xLimits. It should be a numeric vector in (m)');
            h.zLimits = input_zLimits;
        end
        %-- background density
        function set.bckDensity(h,input_bckDensity)
            assert(numel(input_bckDensity)==1&&isnumeric(input_bckDensity), 'Wrong format of the background density bckDensity. It should be a numeric scalar');
            h.bckDensity = input_bckDensity;
        end
        %-- occlusion center X
        function set.occlusionCenterX(h,input)
            assert(isnumeric(input), 'Wrong format of the occlusion center occlusionCenterX. It should be a numeric vector in (m)');
            h.occlusionCenterX = input;
        end
        %-- occlusion center Z
        function set.occlusionCenterZ(h,input)
            assert(isnumeric(input), 'Wrong format of the occlusion center occlusionCenterZ. It should be a numeric vector in (m)');
            h.occlusionCenterZ = input;
        end        
        %-- occlusion diameter
        function set.occlusionDiameter(h,input)
            assert(isnumeric(input), 'Wrong format of the occlusion diameter occlusion. It should be a numeric vector in (m)');
            h.occlusionDiameter = input;
        end
        %-- RoiCenterX
        function set.RoiCenterX(h,input)
            assert(isnumeric(input), 'Wrong format of the ROI center RoiCenterX. It should be a numeric vector in (m)');
            h.RoiCenterX = input;
        end     
        %-- RoiCenterZ
        function set.RoiCenterZ(h,input)
            assert(isnumeric(input), 'Wrong format of the ROI center RoiCenterZ. It should be a numeric vector in (m)');
            h.RoiCenterZ = input;
        end     
        %-- RoiPsfTimeX
        function set.RoiPsfTimeX(h,input)
            assert(isnumeric(input), 'Wrong format of the ROI psf times RoiPsfTimeX. It should be a numeric vector in (m)');
            h.RoiPsfTimeX = input;
        end        
        %-- RoiPsfTimeZ
        function set.RoiPsfTimeZ(h,input)
            assert(isnumeric(input), 'Wrong format of the ROI psf times RoiPsfTimeZ. It should be a numeric vector in (m)');
            h.RoiPsfTimeZ = input;
        end              
        %-- axial resolution of the system
        function set.axialResolution(h,input)
            assert(numel(input)==1&&isnumeric(input), 'Wrong format of the axial resolution value. It should be a numeric scalar in (m)');
            h.axialResolution = input;
        end
        %-- lateral resolution of the system
        function set.lateralResolution(h,input)
            assert(numel(input)==1&&isnumeric(input), 'Wrong format of the lateral resolution value. It should be a numeric scalar in (m)');
            h.lateralResolution = input;
        end        
        
    end
        
    %-- Main method call to generate phantom
    methods (Access = public)
       
        function setPtsScatterersMode(h)            
            h.mode = 0;
        end  
        
        function setOcclusionMode(h)            
            h.mode = 1;
        end
                
        function generate(h)
            switch (h.mode)
                case 0
                    h.generatePtsScatterersPht();                    
                case 1
                    h.generateOcclusionPht();
            end
        end
        
        
        function read_file(h,filename)
            
            %-- Reads all the information from a mat or hdf5 file
            %-- Syntax:
            %-- read_file(file_name)
            %-- file_name: Name of the mat or hdf5 file
            
            [pathstr, name, ext] = fileparts(filename); 
            switch ext
                case '.mat'
                    h.read_file_mat(filename);
                case '.hdf5'
                    h.read_file_hdf5(filename);
                otherwise
                    error('Unknown signal format!');
            end
            
        end
        
        
        function write_file(h,filename)
            
            %-- Write all the information into a mat or hdf5 file
            %-- Syntax:
            %-- write_file(file_name)
            %-- file_name: Name of the mat or hdf5 file
            
            [pathstr, name, ext] = fileparts(filename); 
            switch ext
                case '.mat'
                    h.write_file_mat(filename);
                case '.hdf5'
                    h.write_file_hdf5(filename);
                otherwise
                    error('Unknown signal format!');
            end
            
        end        
        
        %-- HDF5 Ultrasound File Format
        function write_file_hdf5(h,filename)
            
            %-- write HDF5 version in the root group
            attr_details.Name = 'version';
            attr_details.AttachedTo = '/';
            attr_details.AttachType = 'group';
            hdf5write(filename, attr_details, 'v.0.0.39');
            
            %-- We create the /US metagroup in case it is not there
            try
                h5info(filename,'/US')
            catch 
                fid = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
                gid = H5G.create(fid,'US','H5P_DEFAULT','H5P_DEFAULT','H5P_DEFAULT');
                H5G.close(gid);
                H5F.close(fid);
            end
            
            %-- We create a unique us_dataset group
            group_name='US_DATASET0000';
            fid = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
            gid = H5G.open(fid,'/US');
            s_gid = H5G.create(gid,group_name,'H5P_DEFAULT','H5P_DEFAULT','H5P_DEFAULT');
            H5G.close(s_gid);
            H5G.close(gid);
            H5F.close(fid);
            location = ['/US/' group_name];
            
            %-- Attributes 
            %-- Dataset type
            file = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
            filetype = H5T.enum_create ('H5T_NATIVE_INT');
            H5T.enum_insert (filetype, 'US', 0); 
            H5T.enum_insert (filetype, 'SR', 1); 
            gid = H5G.open(file,location);
            space = H5S.create_simple (1,1,[]);

            attr = H5A.create (gid, 'type', filetype, space, 'H5P_DEFAULT');
            H5A.write(attr, filetype, uint32(0));  % <--- US
            H5A.close(attr);
            H5G.close(gid);    
            H5S.close(space);
            H5T.close(filetype);
            H5F.close(file);
            
            %-- us_dataset subtype
            file = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
            filetype = H5T.enum_create ('H5T_NATIVE_INT');
            H5T.enum_insert(filetype, 'STA', 0); 
            H5T.enum_insert(filetype, 'CPW', 1); 
            H5T.enum_insert(filetype, 'VS', 2); 
            H5T.enum_insert(filetype, 'BS', 3); 
            gid = H5G.open(file,location);
            space = H5S.create_simple (1,1,[]);
            attr = H5A.create(gid, 'subtype', filetype, space, 'H5P_DEFAULT');            
            H5A.write(attr, filetype, uint32(1)); % <---- TYPE CPWC
                        
            %-- Signal format 
            file = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
            filetype = H5T.enum_create ('H5T_NATIVE_INT');
            H5T.enum_insert (filetype, 'RF', 0); 
            H5T.enum_insert (filetype, 'IQ', 1); 
            gid = H5G.open(file,location);
            space = H5S.create_simple (1,1,[]);
            attr = H5A.create (gid, 'signal_format', filetype, space, 'H5P_DEFAULT');
            H5A.write(attr, filetype, uint32(1));  % <--- IQ
            H5A.close(attr);
            H5G.close(gid);    
            H5S.close(space);
            H5T.close(filetype);
            H5F.close(file);

            %-- add name
            attr_details.Name = 'name';
            attr_details.AttachedTo = location;
            attr_details.AttachType = 'group';
            hdf5write(filename, attr_details, h.name, 'WriteMode', 'append');

            %-- add author
            attr_details.Name = 'author';
            attr_details.AttachedTo = location;
            attr_details.AttachType = 'group';
            hdf5write(filename, attr_details, h.author, 'WriteMode', 'append');            
            
            %-- add affiliation
            attr_details.Name = 'affiliation';
            attr_details.AttachedTo = location;
            attr_details.AttachType = 'group';
            hdf5write(filename, attr_details, h.affiliation, 'WriteMode', 'append');            
                        
            %-- add date
            attr_details.Name = 'creation_date';
            attr_details.AttachedTo = location;
            attr_details.AttachType = 'group';
            hdf5write(filename, attr_details, h.creation_date, 'WriteMode', 'append');

            %-- Data
            
            %-- Common attributes
            
            dset_details.Location = location;
            dset_details.Name = 'nb_scatterers';
            hdf5write(filename, dset_details, single(h.N_scatterers), 'WriteMode', 'append');           
            %--
            dset_details.Location = location;
            dset_details.Name = 'scatterers_amplitude';
            hdf5write(filename, dset_details, single(h.amp), 'WriteMode', 'append');            
            %--
            dset_details.Location = location;
            dset_details.Name = 'scatterers_positions';
            hdf5write(filename, dset_details, single(h.sca), 'WriteMode', 'append');               
            %--
            dset_details.Location = location;
            dset_details.Name = 'phantom_xLimits';
            hdf5write(filename, dset_details, single(h.xLimits), 'WriteMode', 'append');               
            %--
            dset_details.Location = location;
            dset_details.Name = 'phantom_zLimits';
            hdf5write(filename, dset_details, single(h.zLimits), 'WriteMode', 'append');               
            %--
            dset_details.Location = location;
            dset_details.Name = 'phantom_bckDensity';
            hdf5write(filename, dset_details, single(h.bckDensity), 'WriteMode', 'append');  
            
            %-- Occlusion attributes
            dset_details.Location = location;
            dset_details.Name = 'phantom_occlusionCenterX';
            hdf5write(filename, dset_details, single(h.occlusionCenterX), 'WriteMode', 'append');
            %--
            dset_details.Location = location;
            dset_details.Name = 'phantom_occlusionCenterZ';
            hdf5write(filename, dset_details, single(h.occlusionCenterZ), 'WriteMode', 'append');
            %--            
            dset_details.Location = location;
            dset_details.Name = 'phantom_occlusionDiameter';
            hdf5write(filename, dset_details, single(h.occlusionDiameter), 'WriteMode', 'append');
            
            %-- Roi attributes
            dset_details.Location = location;
            dset_details.Name = 'phantom_RoiCenterX';
            hdf5write(filename, dset_details, single(h.RoiCenterX), 'WriteMode', 'append');
            %--
            dset_details.Location = location;
            dset_details.Name = 'phantom_RoiCenterZ';
            hdf5write(filename, dset_details, single(h.RoiCenterZ), 'WriteMode', 'append');
            %--
            dset_details.Location = location;
            dset_details.Name = 'phantom_RoiPsfTimeX';
            hdf5write(filename, dset_details, single(h.RoiPsfTimeX), 'WriteMode', 'append');
            %--
            dset_details.Location = location;
            dset_details.Name = 'phantom_RoiPsfTimeZ';
            hdf5write(filename, dset_details, single(h.RoiPsfTimeZ), 'WriteMode', 'append');
            
            %-- Point scatterers
            dset_details.Location = location;
            dset_details.Name = 'phantom_xPts';
            hdf5write(filename, dset_details, single(h.xPts), 'WriteMode', 'append');
            %--
            dset_details.Location = location;
            dset_details.Name = 'phantom_zPts';
            hdf5write(filename, dset_details, single(h.zPts), 'WriteMode', 'append');       
            
            %-- System resolution
            %-- Axial
            dset_details.Location = location;
            dset_details.Name = 'phantom_axialResolution';
            hdf5write(filename, dset_details, single(h.axialResolution), 'WriteMode', 'append');               
            %-- Lateral
            dset_details.Location = location;
            dset_details.Name = 'phantom_lateralResolution';
            hdf5write(filename, dset_details, single(h.lateralResolution), 'WriteMode', 'append');            
                       
        end    
    
   
        function read_file_hdf5(h,filename)
            
            %-- read US metagroup
            info=h5info(filename,'/US');

            %-- read the groups in the metagroup
            for n=1:length(info.Groups)
                location=info.Groups(n).Name;
                dstype=h5readatt(filename,location,'type');
                if strcmp(dstype,'US')
                    subtype=h5readatt(filename,location,'subtype');
                    if strcmp(subtype{1},'CPW')


                        %-- Attributes
                        %-- read name
                        a=h5readatt(filename,location,'name');
                        h.name = a{1}(1:end-1);
                        
                        %-- read author
                        a = h5readatt(filename,location,'author');
                        h.author = a{1}(1:end-1);
                        
                        %-- read affiliation
                        a = h5readatt(filename,location,'affiliation'); 
                        h.affiliation = a{1}(1:end-1);
                        
                        %-- read date
                        a = h5readatt(filename,location,'creation_date'); 
                        h.creation_date = a{1}(1:end-1);
                        

                        %-- Data
                        %-- read N_scatterers
                        h.N_scatterers = h5read(filename,[location '/nb_scatterers']);
 
                        %-- read scatterers_amplitude
                        h.amp = h5read(filename,[location '/scatterers_amplitude']); 
                        h.amp = double(h.amp);
                        
                        %-- read scatterers_amplitude
                        h.sca = h5read(filename,[location '/scatterers_positions']);
                        h.sca = double(h.sca);
                        if ( (size(h.sca,1)==3) && (size(h.sca,2)==1) )
                            h.sca = h.sca';
                        end
                        
                        %-- read scatterers_amplitude
                        h.xLimits = h5read(filename,[location '/phantom_xLimits']);                        
                        
                        %-- read scatterers_amplitude
                        h.zLimits = h5read(filename,[location '/phantom_zLimits']);                        

                        %-- read scatterers_amplitude
                        h.bckDensity = h5read(filename,[location '/phantom_bckDensity']);
                        
                        %-- read phantom_occlusionCenterX
                        h.occlusionCenterX = h5read(filename,[location '/phantom_occlusionCenterX']);
                        
                        %-- read phantom_occlusionCenterZ
                        h.occlusionCenterZ = h5read(filename,[location '/phantom_occlusionCenterZ']);                        
                        
                        %-- read phantom_occlusionDiameter
                        h.occlusionDiameter = h5read(filename,[location '/phantom_occlusionDiameter']);

                        %-- read phantom_RoiCenterX
                        h.RoiCenterX = h5read(filename,[location '/phantom_RoiCenterX']);
                        
                        %-- read phantom_RoiCenterZ
                        h.RoiCenterZ = h5read(filename,[location '/phantom_RoiCenterZ']);
                        
                        %-- read phantom_RoiPsfTimeX
                        h.RoiPsfTimeX = h5read(filename,[location '/phantom_RoiPsfTimeX']);
                        
                        %-- read phantom_RoiPsfTimeZ
                        h.RoiPsfTimeZ = h5read(filename,[location '/phantom_RoiPsfTimeZ']);
                        
                        %-- read phantom_xPts
                        h.xPts = h5read(filename,[location '/phantom_xPts']);
                        
                        %-- read phantom_zPts
                        h.zPts = h5read(filename,[location '/phantom_zPts']);               
                         
                        %-- read axial resolution 
                        h.axialResolution = h5read(filename,[location '/phantom_axialResolution']);
                        
                        %-- read lateral resolution
                        h.lateralResolution = h5read(filename,[location '/phantom_lateralResolution']);
                        
                    end
                end
            end
        end
        
        
        function read_file_mat(h,filename)
           
            %-- Load mat file
            load(filename);
            
            %-- Attributes
            %-- read name
            h.name = PARAM.name;

            %-- read author
            h.author = PARAM.author;

            %-- read affiliation
            h.affiliation = PARAM.affiliation;

            %-- read date
            h.creation_date = PARAM.creation_date;

            %-- Data
            %-- read N_scatterers
            h.N_scatterers = PARAM.N_scatterers;

            %-- read scatterers_amplitude
            h.amp = PARAM.amp;

            %-- read scatterers_amplitude
            h.sca = PARAM.sca;
            if ( (size(h.sca,1)==3) && (size(h.sca,2)==1) )
                h.sca = h.sca';
            end

            %-- read scatterers_amplitude
            h.xLimits = PARAM.xLimits;

            %-- read scatterers_amplitude
            h.zLimits = PARAM.zLimits;

            %-- read scatterers_amplitude
            h.bckDensity = PARAM.bckDensity;

            %-- read phantom_occlusionCenterX
            h.occlusionCenterX = PARAM.occlusionCenterX;

            %-- read phantom_occlusionCenterZ
            h.occlusionCenterZ = PARAM.occlusionCenterZ;

            %-- read phantom_occlusionDiameter
            h.occlusionDiameter = PARAM.occlusionDiameter;

            %-- read phantom_RoiCenterX
            h.RoiCenterX = PARAM.RoiCenterX;

            %-- read phantom_RoiCenterZ
            h.RoiCenterZ = PARAM.RoiCenterZ;

            %-- read phantom_RoiPsfTimeX
            h.RoiPsfTimeX = PARAM.RoiPsfTimeX;

            %-- read phantom_RoiPsfTimeZ
            h.RoiPsfTimeZ = PARAM.RoiPsfTimeZ;

            %-- read phantom_xPts
            h.xPts = PARAM.xPts;

            %-- read phantom_zPts
            h.zPts = PARAM.zPts;

            %-- read axial resolution 
            h.axialResolution = PARAM.axialResolution;

            %-- read lateral resolution
            h.lateralResolution = PARAM.lateralResolution;
                        
        end
        
        
        function write_file_mat(h,filename)        
        
            %-- Attributes
            %-- add name
            PARAM.name = h.name;

            %-- add author
            PARAM.author = h.author;

            %-- add affiliation
            PARAM.affiliation = h.affiliation;

            %-- add date
            PARAM.creation_date = h.creation_date;

            %-- Data
            %-- add N_scatterers
            PARAM.N_scatterers = h.N_scatterers;

            %-- add scatterers_amplitude
            PARAM.amp = single(h.amp);

            %-- add scatterers_amplitude
            PARAM.sca = single(h.sca);

            %-- add x limit
            PARAM.xLimits = single(h.xLimits);

            %-- add z limit
            PARAM.zLimits = single(h.zLimits);

            %-- add scatterers_amplitude
            PARAM.bckDensity = single(h.bckDensity);

            %-- add phantom_occlusionCenterX
            PARAM.occlusionCenterX = single(h.occlusionCenterX);

            %-- add phantom_occlusionCenterZ
            PARAM.occlusionCenterZ = single(h.occlusionCenterZ);

            %-- add phantom_occlusionDiameter
            PARAM.occlusionDiameter = single(h.occlusionDiameter);

            %-- add phantom_RoiCenterX
            PARAM.RoiCenterX = single(h.RoiCenterX);

            %-- add phantom_RoiCenterZ
            PARAM.RoiCenterZ = single(h.RoiCenterZ);

            %-- add phantom_RoiPsfTimeX
            PARAM.RoiPsfTimeX = single(h.RoiPsfTimeX);

            %-- add phantom_RoiPsfTimeZ
            PARAM.RoiPsfTimeZ = single(h.RoiPsfTimeZ);

            %-- read phantom_xPts
            PARAM.xPts = single(h.xPts);

            %-- read phantom_zPts
            PARAM.zPts = single(h.zPts);

            %-- read axial resolution 
            PARAM.axialResolution = single(h.axialResolution);

            %-- read lateral resolution
            PARAM.lateralResolution = single(h.lateralResolution);
                    
            %-- Write mat file
            save(filename,'PARAM');
            
        end

        
    end

    
    methods (Access = private)
        
        function generatePtsScatterersPht(h) 
            
            if ( size(h.xPts,1)==1 )
                xxp = h.xPts';
                zzp = h.zPts';
            else
                xxp = h.xPts;
                zzp = h.zPts;
            end
            h.N_scatterers = length(zzp(:));                    %-- total number of scatterers
            h.sca = [xxp(:) zeros(h.N_scatterers,1) zzp(:)];    %-- list with the scatterers coordinates [m]
            h.amp = ones(h.N_scatterers,1);                
            
        end
        
        function generateOcclusionPht(h)
            
            psf_length = h.axialResolution;
            psf_width = h.lateralResolution;
            width = h.xLimits(2) - h.xLimits(1);
            depth = h.zLimits(2) - h.zLimits(1);
            h.N_scatterers = round( h.bckDensity * depth * width / psf_length / psf_width );
            x = width * (rand(h.N_scatterers,1)-0.5);
            y = zeros(size(x));
            z = depth * (rand(h.N_scatterers,1))+h.zLimits(1);
            amp = randn(h.N_scatterers,1);
            
            for k=1:length(h.occlusionDiameter)           
                if (h.occlusionDiameter(k)>0)
                    r = h.occlusionDiameter(k) / 2;          
                    xc = h.occlusionCenterX(k);
                    zc = h.occlusionCenterZ(k);  
                    outside = ( ((x-xc).^2 + (z-zc).^2) >= r^2);
                    x = x(outside);
                    y = y(outside);
                    z = z(outside);
                    amp = amp(outside);
                end
            end
            h.sca = [x y z];
            h.amp = amp;
            
        end
        
    end
    
end
