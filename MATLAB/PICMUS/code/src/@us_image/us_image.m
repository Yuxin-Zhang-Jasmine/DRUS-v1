classdef us_image < handle

    %-- Class containing the beamformed data
    
    %-- Authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %--          Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)    

    properties  (SetAccess = public)
        
        %-- administration
        name                        %-- String containing the name of the beamformed image
        author                      %-- String containing the name of the author(s) of the beamformed image
        affiliation                 %-- String containing the affiliation of the author(s) of the beamformed image
        algorithm                   %-- String containing the name of the algorithm used to create the beamformed image
        creation_date               %-- String containing the date the reconstruction was created
        
        %-- input data
        scan                        %-- SCAN class defining the scan area 
        
        %-- output data
        number_plane_waves          %-- Vector containing number of plane waves used in each reconstructed frame
        data                        %-- Matrix containing the envelope of the reconstructed signal 
        
        %- information
        transmit_f_number           %-- Scalar of the F-number used on transmit
        transmit_apodization_window %-- String describing the transmit apodization window
        receive_f_number            %-- Scalar of the F-number used on receive
        receive_apodization_window  %-- String describing the receive apodization window 
        
    end

    %-- version
    properties  (SetAccess = private)
        version                     %-- version of this us_image class, needed for back compatibility
    end
    
    %-- constructor
    methods (Access = public)
        
        function h = us_image(name)
            
            %-- Constructor of the US_IMAGE class.
            %-- Syntax:
            %-- US_IMAGE(name) 
            %-- name: Name of the reconstruction

            h.creation_date = sprintf('%d/%02d/%d %02d:%02d:%02.2f',clock);
            h.version = 'v0.10';
            
            if exist('name') 
                h.name=name;
            else
                h.name='';
            end
            author = '';
            affiliation = '';
            algorithm = '';
            creation_date = '';
            number_plane_waves = NaN;
            transmit_f_number = NaN;
            transmit_apodization_window = '';
            receive_f_number = NaN;
            receive_apodization_window = '';
            
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
        
        %-- algorithm
        function set.algorithm(h,input)
            assert(isstr(input), 'Wrong format of the algorithm name. It should be a string.');
            h.algorithm = input;
        end  
        
        %-- creation_date
        function set.creation_date(h,input_date)
            assert(isstr(input_date), 'Wrong format of the creation date. It should be a string.');
            h.creation_date = input_date;
        end
        
        %-- set scan
        function set.scan(h,input)
            assert(isa(input,'linear_scan'), 'Wrong format of the scan. It should be a LINEAR_SCAN class.');
            h.scan = input;
        end
        
        %-- data
        function set.data(h,input_data)
            assert(isnumeric(input_data)&&ndims(input_data)>1&&ndims(input_data)<4, 'Wrong format of the data. It should be a numeric array of dimensions [z_axis, x_axis, frames].');
            %-- scan
            if(~isempty(h.scan))
                if(size(input_data,1)~=length(h.scan.z_axis)) error('The data size in the z-axis does not match the size of the scan area');end
                if(size(input_data,2)~=length(h.scan.x_axis)) error('The data size in the x-axis does not match the size of the scan area');end
            end            
            %-- number of plane waves used
            assert(~any(isnan(h.number_plane_waves)),'The number of plane waves must be set before assigning the data');
            if(size(input_data,3)~=length(h.number_plane_waves)) error('The number of images size(data,3) should match the number of elements in the vector number_plane_waves');end            
            %-- data
            h.data = input_data;
            %-- data check
            assert(isreal(h.data)&~any(h.data(:)<0),'Data provided must be the envelope of the beamformed signal, i.e. the magnitude of the beamformed IQ signal or the magnitude of the Hilbert transform of a RF beamformed signal. Do not apply compression to the signal envelope.'); 
        end
        
    end
    
    
    %-- presentation methods
    methods (Access = public)
        
        function im = show(h,dynamic_range,frame_list)
            
            %-- Plots the envelope of the beamformed data and returns a copy of the image
            %-- Syntax:
            %-- image = show(dynamic_range)
            %-- image: Matrix with the output signal
            %-- dynamic_range: Desired dynamic range of the displayed images in dB
            %-- frame_list: List of frames to be displayed
       
            if ~exist('dynamic_range') dynamic_range=60; end
            
            %-- setting axis limits (mm)
            x_lim = [min(h.scan.x_matrix(:)) max(h.scan.x_matrix(:))]*1e3; 
            z_lim = [min(h.scan.z_matrix(:)) max(h.scan.z_matrix(:))]*1e3; 
            
            %-- ploting image reconstruction
            figure; set(gca,'fontsize',16);
            if ~exist('frame_list') || isempty(frame_list) 
                frame_list = 1:size(h.data,3); 
            end
            for f=frame_list
                
                    %-- compute dB values
                    env = h.data(:,:,f);
                    im = 20*log10(env./max(env(:)));
                    vrange = [-dynamic_range 0];

                    %-- display image
                    imagesc((h.scan.x_axis)*1e3,(h.scan.z_axis)*1e3,im); 
                    shading flat; colormap gray; caxis(vrange); colorbar; hold on;
                    axis equal manual;
                    xlabel('x [mm]');
                    ylabel('z [mm]');
                    set(gca,'YDir','reverse');
                    set(gca,'fontsize',16);
                    axis([x_lim z_lim]);
                    title(sprintf('%s\n %d plane waves',char(h.name),h.number_plane_waves(f)));
                    drawnow; hold off;
                    pause(0.5);
                    
            end
            
        end
        
    end
    
    %-- HDF5 Ultrasound File Format
    methods (Access = public)

        
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
        
        
        function write_file_hdf5(h,filename)
            
            %-- Writes all the information in the us_dataset to a HUFF (HDF5 Ultrasound File Format) file
            %-- Syntax:
            %-- write_file(file_name,group_name)
            %-- file_name: Name of the hdf5 file
            %-- group_name: Name of the group destination
            
            %-- write HDF5 version in the root group
            attr_details.Name = 'version';
            attr_details.AttachedTo = '/';
            attr_details.AttachType = 'group';
            hdf5write(filename, attr_details, 'v.0.0.41');
            
            %-- create the /US metagroup in case it is not there
            try
                h5info(filename,'/US')
            catch 
                fid = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
                gid = H5G.create(fid,'US','H5P_DEFAULT','H5P_DEFAULT','H5P_DEFAULT');
                H5G.close(gid);
                H5F.close(fid);
            end
            
            %-- create a unique us_dataset group
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
            filetype = H5T.enum_create ('H5T_NATIVE_INT');
                H5T.enum_insert (filetype, 'US', 0); 
                H5T.enum_insert (filetype, 'SR', 1); 
            gid = H5G.open(file,location);
            space = H5S.create_simple (1,1,[]);

            attr = H5A.create (gid, 'type', filetype, space, 'H5P_DEFAULT');
            H5A.write (attr, filetype, uint32(1));  % <--- SR
            H5A.close (attr);
            H5G.close(gid);    
            H5S.close (space);
            H5T.close (filetype);
            H5F.close (file);
            
            %-- Signal format 
            file = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
            filetype = H5T.enum_create ('H5T_NATIVE_INT');
                H5T.enum_insert (filetype, 'RF', 0); 
                H5T.enum_insert (filetype, 'IQ', 1); 
                H5T.enum_insert (filetype, 'ENV', 2);                 
            gid = H5G.open(file,location);
            space = H5S.create_simple (1,1,[]);

            attr = H5A.create (gid, 'signal_format', filetype, space, 'H5P_DEFAULT');
            H5A.write (attr, filetype, uint32(2));  % <--- ENVELOPE
            H5A.close (attr);
            H5G.close(gid);    
            H5S.close (space);
            H5T.close (filetype);
            H5F.close (file);

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

            %-- add algorithm
            attr_details.Name = 'algorithm';
            attr_details.AttachedTo = location;
            attr_details.AttachType = 'group';
            hdf5write(filename, attr_details, h.algorithm, 'WriteMode', 'append');
            
            %-- add version of the us_dataset
            attr_details.Name = 'version';
            attr_details.AttachedTo = location;
            attr_details.AttachType = 'group';
            hdf5write(filename, attr_details, h.version, 'WriteMode', 'append');
            
            %-- add date
            attr_details.Name = 'creation_date';
            attr_details.AttachedTo = location;
            attr_details.AttachType = 'group';
            hdf5write(filename, attr_details, h.creation_date, 'WriteMode', 'append');

            %-- Data

            %-- add linear scan
            dset_details.Location = [location '/scan'];
            dset_details.Name = 'x_axis';
            hdf5write(filename, dset_details, single(h.scan.x_axis), 'WriteMode', 'append');
            
            dset_details.Location = [location '/scan'];
            dset_details.Name = 'z_axis';
            hdf5write(filename, dset_details, single(h.scan.z_axis), 'WriteMode', 'append');

            %-- F number
            dset_details.Location = location;
            dset_details.Name = 'transmit_f_number';
            hdf5write(filename, dset_details, single(h.transmit_f_number), 'WriteMode', 'append');
            dset_details.Location = location;
            dset_details.Name = 'receive_f_number';
            hdf5write(filename, dset_details, single(h.receive_f_number), 'WriteMode', 'append');
            
            %-- apodization windows
            dset_details.Location = location;
            dset_details.Name = 'transmit_apodization_window';
            hdf5write(filename, dset_details, char(h.transmit_apodization_window), 'WriteMode', 'append');
            dset_details.Location = location;
            dset_details.Name = 'receive_apodization_window';
            hdf5write(filename, dset_details, char(h.receive_apodization_window), 'WriteMode', 'append');
            
            %-- add number of plane waves
            dset_details.Location = location;
            dset_details.Name = 'number_plane_waves';
            hdf5write(filename, dset_details, single(h.number_plane_waves), 'WriteMode', 'append');
                        
            %-- add data
            dset_details.Location = [location '/data'];
            dset_details.Name = 'real';
            hdf5write(filename, dset_details, single(real(h.data)), 'WriteMode', 'append');
            dset_details.Name = 'imag';
            hdf5write(filename, dset_details, single(imag(h.data)), 'WriteMode', 'append');
            
        end
        
        function read_file_hdf5(h,filename)
            
            %-- Reads all the information from a HUFF (HDF5 Ultrasound File Format) file
            %-- Syntax:
            %-- read_file(file_name,location)
            %-- file_name                    Name of the hdf5 file
            
            %-- check the version here
            version_cell = h5readatt(filename,'/','version');
            version = version_cell{1}(1:8);
            
            switch(version)
                case 'v.0.0.40'
                    warning('The version of the hdf5 is obsolete. Please generate the image again.');
                    
                    %-- read US metagroup
                    info = h5info(filename,'/US');

                    %-- read the groups in the metagroup
                    for n=1:length(info.Groups)
                        location = info.Groups(n).Name;
                        dstype = h5readatt(filename,location,'type');
                        if strcmp(dstype,'SR')
                            %-- read signal format 
                            signal_format = h5readatt(filename,location,'signal_format');
                            switch(signal_format{1})
                                case 'RF'
                                    error('RF format not available!');
                                case 'IQ'
                                    error('RF format not available!');
                                case 'ENV'
                                    ;
                                otherwise
                                    error('Unknown signal format!');
                            end

                            %-- Attributes
                            %-- read name
                            a = h5readatt(filename,location,'name'); h.name=a{1}(1:end-1);

                            %-- read author
                            a = h5readatt(filename,location,'author'); h.author=a{1}(1:end-1);

                            %-- read affiliation
                            a = h5readatt(filename,location,'affiliation'); h.affiliation=a{1}(1:end-1);

                            %-- read algorithm
                            a = h5readatt(filename,location,'algorithm'); h.algorithm=a{1}(1:end-1);

                            %-- read date
                            a = h5readatt(filename,location,'creation_date'); h.creation_date=a{1}(1:end-1);

                            %-- read version
                            a = h5readatt(filename,location,'version'); h.version=a{1}(1:end-1);

                            %-- Data

                            %-- read scan
                            x_axis = h5read(filename,[location '/scan/x_axis']);
                            z_axis = h5read(filename,[location '/scan/z_axis']);
                            h.scan = linear_scan(x_axis,z_axis);

                            %-- F-numbers
                            h.transmit_f_number = h5read(filename,[location '/transmit_f_number']);
                            h.receive_f_number = h5read(filename,[location '/receive_f_number']);

                            %-- Apodization window
                            h.transmit_apodization_window = h5read(filename,[location '/transmit_apodization_window']);
                            h.receive_apodization_window = h5read(filename,[location '/receive_apodization_window']);

                            %-- read data
                            real_part = h5read(filename,[location '/data/real']);
                            imag_part = h5read(filename,[location '/data/imag']);
                            h.number_plane_waves = zeros(1,size(real_part,3));
                            h.data = real_part+1i*imag_part;
                            
                        end
                        
                    end
                    
                case 'v.0.0.41'
                    %-- read US metagroup
                    info = h5info(filename,'/US');

                    %-- read the groups in the metagroup
                    for n=1:length(info.Groups)
                        location = info.Groups(n).Name;
                        dstype = h5readatt(filename,location,'type');
                        if strcmp(dstype,'SR')
                            %-- read signal format 
                            signal_format = h5readatt(filename,location,'signal_format');
                            switch(signal_format{1})
                                case 'RF'
                                    error('RF format not available!');
                                case 'IQ'
                                    error('RF format not available!');
                                case 'ENV'
                                    ;
                                otherwise
                                    error('Unknown signal format!');
                            end

                            %-- Attributes
                            %-- read name
                            a = h5readatt(filename,location,'name'); h.name=a{1}(1:end-1);

                            %-- read author
                            a = h5readatt(filename,location,'author'); h.author=a{1}(1:end-1);

                            %-- read affiliation
                            a = h5readatt(filename,location,'affiliation'); h.affiliation=a{1}(1:end-1);

                            %-- read algorithm
                            a = h5readatt(filename,location,'algorithm'); h.algorithm=a{1}(1:end-1);

                            %-- read date
                            a = h5readatt(filename,location,'creation_date'); h.creation_date=a{1}(1:end-1);

                            %-- read version
                            a = h5readatt(filename,location,'version'); h.version=a{1}(1:end-1);

                            %-- Data

                            %-- read scan
                            x_axis = h5read(filename,[location '/scan/x_axis']);
                            z_axis = h5read(filename,[location '/scan/z_axis']);
                            h.scan = linear_scan(x_axis,z_axis);

                            %-- F-numbers
                            h.transmit_f_number = h5read(filename,[location '/transmit_f_number']);
                            h.receive_f_number = h5read(filename,[location '/receive_f_number']);

                            %-- Apodization window
                            h.transmit_apodization_window = h5read(filename,[location '/transmit_apodization_window']);
                            h.receive_apodization_window = h5read(filename,[location '/receive_apodization_window']);

                            %-- read number of plane waves
                            h.number_plane_waves = h5read(filename,[location '/number_plane_waves']);

                            %-- read data
                            real_part = h5read(filename,[location '/data/real']);
                            imag_part = h5read(filename,[location '/data/imag']);
                            h.data = real_part+1i*imag_part;
                        end
                    end
                otherwise
                    error(sprintf('HUFF version %s not supported',version));
            end  
            
        end
        

        function read_file_mat(h,filename)
           
            %-- Read mat file
            load(filename);
            
            %-- Attributes
            %-- read name
            h.name = PARAM.name;

            %-- read author
            h.author = PARAM.author;

            %-- read affiliation
            h.affiliation = PARAM.affiliation;

            %-- read algorithm
            h.algorithm = PARAM.algorithm;

            %-- read date
            h.creation_date = PARAM.creation_date;

            %-- read version
            h.version = PARAM.version;

            %-- Data

            %-- read scan
            x_axis = PARAM.x_axis;
            z_axis = PARAM.z_axis;
            h.scan = linear_scan(x_axis,z_axis);

            %-- F-numbers
            h.transmit_f_number = PARAM.transmit_f_number;
            h.receive_f_number = PARAM.receive_f_number;

            %-- Apodization window
            h.transmit_apodization_window = PARAM.transmit_apodization_window;
            h.receive_apodization_window = PARAM.receive_apodization_window;

            %-- read number of plane waves
            h.number_plane_waves = PARAM.number_plane_waves;

            %-- read data
            real_part = PARAM.real_part;
            imag_part = PARAM.imag_part;
            h.data = real_part+1i*imag_part;

        end
        
        
        function write_file_mat(h,filename)
            
            %-- Attributes
            %-- add name
            PARAM.name = h.name;

            %-- add author
            PARAM.author = h.author;

            %-- add affiliation
            PARAM.affiliation = h.affiliation;

            %-- add algorithm
            PARAM.algorithm = h.algorithm;

            %-- add date
            PARAM.creation_date = h.creation_date;

            %-- add version
            PARAM.version = h.version;


            %-- Data

            %-- add scan
            PARAM.x_axis = single(h.scan.x_axis);
            PARAM.z_axis = single(h.scan.z_axis);


            %-- F-numbers
            PARAM.transmit_f_number = single(h.transmit_f_number);
            PARAM.receive_f_number = single(h.receive_f_number);

            %-- Apodization window
            PARAM.transmit_apodization_window = char(h.transmit_apodization_window);
            PARAM.receive_apodization_window = char(h.receive_apodization_window);

            %-- read number of plane waves
            PARAM.number_plane_waves = single(h.number_plane_waves);

            %-- read data
            PARAM.real_part = single(real(h.data));
            PARAM.imag_part = single(imag(h.data));
            
            %-- write mat file
            save(filename,'PARAM');
            
        end
        
    end
    
end

