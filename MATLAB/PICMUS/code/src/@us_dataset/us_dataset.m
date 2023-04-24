classdef us_dataset < handle
   
    %-- Class defining a standard way to store ultrasound rawdata information
    
    %-- Authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %--          Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)

    %-- $Date: 2016/03/01 $    

    
    %-- public properties
    properties (SetAccess = public)
        name                    %-- String containing the name of the us_dataset
        creation_date           %-- String containing the date the us_dataset class was created
        probe_geometry          %-- Matrix containing the position of the acoustic center of each element in the probe [x, y, z] (m)
        data                    %-- Matrix containing the us_dataset in the order [time_samples, channels, firings, frames]
        c0                      %-- reference speed of sound (m/s)
        initial_time            %-- Time interval between the transmit event and the first sample acquired (s)
        sampling_frequency      %-- Sampling frequency (Hz)
        modulation_frequency    %-- Modulation frequency (Hz)
        PRF                     %-- Pulse repetition frequency, time interval between consecutive transmit events (Hz)
        angles                  %-- vector containing the angles in the plane-wave sequence (rad)
    end
   
    %-- internal properties
    properties  (SetAccess = protected)   
        samples                 %-- number of samples in the fast time direction
        channels                %-- number of channels in the transducer
        firings                 %-- number of firings in the sequence (i.e. number of plane waves in CPWC)
        frames                  %-- number of frames in the us_dataset        
        version                 %-- version of this us_dataset class, needed for back compatibility
    end
    
    %-- constructors
    methods (Access = public)
        
        function h = us_dataset(input_name)
            
            %-- Constructor of the CPWC_US_DATASET class.
            %-- Syntax:
            %-- CPW(name) 
            %-- name: Name of the us_dataset
        
            if exist('input_name') h.name = input_name; end
            h.creation_date = sprintf('%d/%02d/%d %02d:%02d:%02.2f',clock);
            h.version = 'v1.96';
            h.PRF = 100;
            
        end
        
    end

    %-- set methods, input format check
    methods  
    
        %-- name
        function set.name(h,input_name)
            assert(isstr(input_name), 'Wrong format of the us_dataset name. It should be a string.');
            h.name = input_name;
        end
        
        %-- creation_date
        function set.creation_date(h,input_date)
            assert(isstr(input_date), 'Wrong format of the us_dataset creation date. It should be a string.');
            h.creation_date = input_date;
        end
        %-- probe geometry
        function set.probe_geometry(h,input_geometry)
            assert(ndims(input_geometry)==2&&size(input_geometry,2)==3&&isnumeric(input_geometry), 'Wrong format of the probe geometry. It should be a numeric three column matrix [x, y, z] in (m)');
            h.probe_geometry = input_geometry;
        end
        %-- c0
        function set.c0(h,input_c0)
            assert(numel(input_c0)==1&&isnumeric(input_c0), 'Wrong format of the reference speed of sound c0. It should be a numeric scalar in (m/s)');
            h.c0 = input_c0;
        end
        %-- angles
        function set.angles(h,input_angles)
            assert(size(input_angles,1)==numel(input_angles)&&isnumeric(input_angles), 'Wrong format of the angle vector. It should be a numeric column vector in (rad)');
            h.angles = input_angles;
        end
        %-- initial_time
        function set.initial_time(h,input_t0)
            assert(numel(input_t0)==1&&isnumeric(input_t0), 'Wrong format of the initial time. It should be a numeric scalar in (s)');
            h.initial_time = input_t0;
        end
        %-- sampling_frequency
        function set.sampling_frequency(h,input_Fs)
            assert(numel(input_Fs)==1&&isnumeric(input_Fs), 'Wrong format of the sampling frequency. It should be a numeric scalar in (Hz)');
            h.sampling_frequency=input_Fs;
        end
        %-- PRF
        function set.PRF(h,input_PRF)
            assert(numel(input_PRF)==1&&isnumeric(input_PRF), 'Wrong format of the pulse repetition frequency. It should be a numeric scalar in (Hz)');
            h.PRF = input_PRF;
        end                
        %-- modulation_frequency
        function set.modulation_frequency(h,input_Fd)
            assert(numel(input_Fd)==1&&isnumeric(input_Fd), 'Wrong format of the modulation frequency. It should be a numeric scalar in (Hz)');
            h.modulation_frequency = input_Fd;
        end
        %-- data
        function set.data(h,input_data)            
            assert(isnumeric(input_data)&&ndims(input_data)<5, 'Wrong format of the data. It should be a numeric array of dimensions [samples, channels, firings, frames].');
            %-- dimensions
            h.samples = size(input_data,1);
            h.channels = size(input_data,2);
            h.firings = size(input_data,3);
            h.frames = size(input_data,4);            
            %-- data
            h.data = input_data;            
            %-- some initial checks
            if(isreal(h.data)&(h.modulation_frequency~=0)) 
                warning('The inserted data is real (RF format), but the modulation frequency is not zero.'); 
            end
            if(~isreal(h.data)&(h.modulation_frequency==0)) 
                warning('The inserted data is complex (IQ format), but the modulation frequency is zero.'); 
            end
            if(h.channels~=size(h.probe_geometry,1)) 
                warning(sprintf('The number of channels in the inserted data (%d) does not match the size of the probe geometry matrix (%d).',h.channels, size(h.probe_geometry,1))); 
            end
            if(h.firings~=size(h.angles,1)) 
                warning(sprintf('The number of firings in the inserted data (%d) does not match the size of the angle vector (%d).',h.firings, size(h.angles,1))); 
            end            
        end
        
    end
      
    %-- HDF5 Ultrasound File Format
    methods (Access = public)    
    
        function write_file_hdf5(h,filename)
            
            %-- Writes all the information in the us_dataset to a HUFF (HDF5 Ultrasound File Format) file
            %-- Syntax:
            %-- write_file_hdf5(file_name,group_name)
            %-- file_name: Name of the hdf5 file
            %-- group_name: Name of the group destination
            
            %-- write HUFF version in the root group
            attr_details.Name = 'version';
            attr_details.AttachedTo = '/';
            attr_details.AttachType = 'group';
            hdf5write(filename, attr_details, 'v.0.0.40');
            
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
            group_name = 'US_DATASET0000';
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
            filetype = H5T.enum_create('H5T_NATIVE_INT');
            H5T.enum_insert(filetype,'US', 0); 
            H5T.enum_insert(filetype,'SR', 1); 
            gid = H5G.open(file,location);
            space = H5S.create_simple(1,1,[]);
            attr = H5A.create (gid,'type',filetype,space,'H5P_DEFAULT');
            H5A.write(attr, filetype, uint32(0));  % <--- US
            H5A.close(attr);
            H5G.close(gid);    
            H5S.close(space);
            H5T.close(filetype);
            H5F.close(file);
            
            %-- us_dataset subtype
            file = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
            filetype = H5T.enum_create('H5T_NATIVE_INT');
            H5T.enum_insert (filetype,'STA',0); 
            H5T.enum_insert (filetype,'CPW',1); 
            H5T.enum_insert (filetype,'VS',2); 
            H5T.enum_insert (filetype,'BS',3); 
            gid = H5G.open(file,location);
            space = H5S.create_simple (1,1,[]);
            attr = H5A.create (gid,'subtype',filetype,space,'H5P_DEFAULT');            
            H5A.write(attr,filetype,uint32(1)); % <---- TYPE CPWC
                        
            %-- Signal format 
            file = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
            filetype = H5T.enum_create('H5T_NATIVE_INT');
            H5T.enum_insert(filetype,'RF',0); 
            H5T.enum_insert(filetype,'IQ',1); 
            gid = H5G.open(file,location);
            space = H5S.create_simple(1,1,[]);

            attr = H5A.create(gid,'signal_format',filetype,space,'H5P_DEFAULT');
            if(h.modulation_frequency>0)
                H5A.write(attr,filetype,uint32(1));  % <--- IQ
            else
                H5A.write(attr,filetype,uint32(0));  % <--- RF
            end
            H5A.close(attr);
            H5G.close(gid);    
            H5S.close(space);
            H5T.close(filetype);
            H5F.close(file);

            %-- add name
            attr_details.Name = 'name';
            attr_details.AttachedTo = location;
            attr_details.AttachType = 'group';
            hdf5write(filename,attr_details,h.name,'WriteMode','append');

            %-- add version of the us_dataset
            attr_details.Name = 'version';
            attr_details.AttachedTo = location;
            attr_details.AttachType = 'group';
            hdf5write(filename,attr_details,h.version,'WriteMode','append');
            
            %-- add date
            attr_details.Name = 'creation_date';
            attr_details.AttachedTo = location;
            attr_details.AttachType = 'group';
            hdf5write(filename,attr_details,h.creation_date,'WriteMode','append');

            %-- Data
            
            %-- speed of sound
            dset_details.Location = location;
            dset_details.Name = 'sound_speed';
            hdf5write(filename, dset_details, single(h.c0), 'WriteMode', 'append');
            
            %-- add initial_time
            dset_details.Location = location;
            dset_details.Name = 'initial_time';
            hdf5write(filename, dset_details, single(h.initial_time), 'WriteMode', 'append');

            %-- add sampling_frequency
            dset_details.Location = location;
            dset_details.Name = 'sampling_frequency';
            hdf5write(filename, dset_details, single(h.sampling_frequency), 'WriteMode', 'append');

            %-- add PRF
            dset_details.Location = location;
            dset_details.Name = 'PRF';
            hdf5write(filename, dset_details, single(h.PRF), 'WriteMode', 'append');
            
            %-- add modulation frequency
            dset_details.Location = location;
            dset_details.Name = 'modulation_frequency';
            if(h.modulation_frequency>0)
                hdf5write(filename, dset_details, single(h.modulation_frequency), 'WriteMode', 'append');
            else
                hdf5write(filename, dset_details, single(0), 'WriteMode', 'append');
            end

            %-- add probe geometry
            dset_details.Location = location;
            dset_details.Name = 'probe_geometry';
            hdf5write(filename, dset_details, single(h.probe_geometry), 'WriteMode', 'append');
            
            %-- add angles
            dset_details.Location = location;
            dset_details.Name = 'angles';
            hdf5write(filename, dset_details, single(h.angles), 'WriteMode', 'append');
            
            %-- add data
            dset_details.Location = [location '/data'];
            dset_details.Name = 'real';
            hdf5write(filename, dset_details, single(real(h.data)), 'WriteMode', 'append');
            dset_details.Name = 'imag';
            hdf5write(filename, dset_details, single(imag(h.data)), 'WriteMode', 'append');
            
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
                
        
        function read_file_hdf5(h,filename)

            %-- Reads all the information from a HUFF (HDF5 Ultrasound File Format) file
            %-- Syntax:
            %-- read_file_hdf5(file_name)
            %-- file_name: Name of the hdf5 file
            
            %-- read US metagroup
            info = h5info(filename,'/US');

            %-- read the groups in the metagroup
            for n=1:length(info.Groups)
                location=info.Groups(n).Name;
                dstype=h5readatt(filename,location,'type');
                if strcmp(dstype,'US')
                    subtype=h5readatt(filename,location,'subtype');
                    if strcmp(subtype{1},'CPW')

                        %-- subtype
                        dataset_subtype=h5readatt(filename,location,'subtype');
                        assert(strcmp(dataset_subtype,'CPW'),'Only CPWC us_dataset are supported!');
                        
                        %-- read signal format 
                        signal_format=h5readatt(filename,location,'signal_format');
                        
                        %-- read modulation frequency
                        h.modulation_frequency=h5read(filename,[location '/modulation_frequency']);

                        %-- check format
                        switch(signal_format{1})
                            case 'RF'
                                assert(h.modulation_frequency==0,'RF dataset cannot have a modulation frequency');
                            case 'IQ'
                                assert(h.modulation_frequency>0,'IQ dataset cannot have a null modulation frequency');
                            otherwise
                                error('Unknown signal format!');
                        end

                        %-- Attributes
                        %-- read name
                        a = h5readatt(filename,location,'name'); h.name=a{1}(1:end-1);

                        %-- read date
                        a = h5readatt(filename,location,'creation_date'); h.creation_date=a{1}(1:end-1);
                        
                        %-- read version
                        a = h5readatt(filename,location,'version'); h.version=a{1}(1:end-1);

                        %-- Data
                        %-- read speed of sound
                        h.c0 = h5read(filename,[location '/sound_speed']);

                        %-- read initial_time
                        h.initial_time = h5read(filename,[location '/initial_time']);

                        %-- read sampling_frequency
                        h.sampling_frequency = h5read(filename,[location '/sampling_frequency']);

                        %-- read sampling_frequency
                        h.PRF = h5read(filename,[location '/PRF']);

                        %-- read transducer geometry
                        h.probe_geometry = h5read(filename,[location '/probe_geometry']);

                        %-- read angles
                        h.angles = h5read(filename,[location '/angles']);
                        
                        %-- read data
                        real_part = h5read(filename,[location '/data/real']);
                        imag_part = h5read(filename,[location '/data/imag']);
                        h.data = real_part+1i*imag_part;
                    end
                end
            end
        end
        
        function read_file_mat(h,filename)

            %-- Reads all the information from a mat file
            %-- Syntax:
            %-- read_file_mat(file_name)
            %-- file_name: Name of the mat file

            %-- Load saved structured
            load(filename);
            
            %-- read signal format 
            signal_format = PARAM.signal_format;

            %-- read modulation frequency
            h.modulation_frequency = PARAM.modulation_frequency;

            %-- check format
            switch(signal_format)
                case 'RF'
                    assert(h.modulation_frequency==0,'RF dataset cannot have a modulation frequency');
                case 'IQ'
                    assert(h.modulation_frequency>0,'IQ dataset cannot have a null modulation frequency');
                otherwise
                    error('Unknown signal format!');
            end

            %-- Attributes
            %-- read name
            h.name = PARAM.name;

            %-- read date
            h.creation_date = PARAM.creation_date;

            %-- read version
            h.version = PARAM.version;

            %-- Data
            %-- read speed of sound
            h.c0 = PARAM.c0;

            %-- read initial_time
            h.initial_time = PARAM.initial_time;

            %-- read sampling_frequency
            h.sampling_frequency = PARAM.sampling_frequency;

            %-- read sampling_frequency
            h.PRF = PARAM.PRF;

            %-- read transducer geometry
            h.probe_geometry = PARAM.probe_geometry;

            %-- read angles
            h.angles = PARAM.angles;

            %-- read data
            real_part = PARAM.real;
            imag_part = PARAM.imag;
            h.data = real_part+1i*imag_part;
            
        end 
        
        
        function write_file_mat(h,filename)

            %-- Writes all the information in the us_dataset to a mat file
            %-- Syntax:
            %-- write_file_mat(file_name,group_name)
            %-- file_name: Name of the mat file
            %-- group_name: Name of the group destination

            %-- Signal format 
            if(h.modulation_frequency>0)
                PARAM.signal_format = 'IQ';
            else
                PARAM.signal_format = 'RF';
            end

            %-- add name
            PARAM.name = h.name;

            %-- add version of the us_dataset
            PARAM.version = h.version;

            %-- add date
            PARAM.creation_date = h.creation_date;


            %-- Data
            
            %-- speed of sound
            PARAM.c0 = single(h.c0);
            
            % add initial_time
            PARAM.initial_time = single(h.initial_time);
            
            %-- add sampling_frequency
            PARAM.sampling_frequency = single(h.sampling_frequency);
            
            %-- add PRF
            PARAM.PRF = single(h.PRF);

            %-- add modulation frequency
            if(h.modulation_frequency>0)
                PARAM.modulation_frequency = single(h.modulation_frequency);
            else
                PARAM.modulation_frequency = single(0);
            end

            %-- add probe geometry
            PARAM.probe_geometry = single(h.probe_geometry);

            %-- add angles
            PARAM.angles = single(h.angles);

            %-- add data
            PARAM.real = single(real(h.data));
            PARAM.imag = single(imag(h.data));

            save(filename,'PARAM');
            
        end         
                
    end
    
end
