classdef linear_scan < handle

    %-- Class defining a linear scan area
    
    %-- Authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %--          Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)

    %-- $Date: 2016/03/01 $      

    properties  (SetAccess = public)
        x_axis      %-- Vector defining the x coordinates of each row of pixels
        z_axis      %-- Vector defining the z coordinates of each column of pixels
    end
    
    properties  (SetAccess = private)
        x_matrix    %-- Matrix containing the x coordinate of each pixel in the matrix
        z_matrix    %-- Matrix containing the z coordinate of each pixel in the matrix
        x           %-- Vector containing the x coordinate of each pixel in the matrix
        z           %-- Vector containing the z coordinate of each pixel in the matrix
        dx          %-- Spatial step in x-axis
        dz          %-- Spatial step in z-axis
        Nx          %-- Number of samples in x-axis
        Nz          %-- Number of samples in z-axis
        pixels      %-- total number of pixels in the matrix
    end
    
    
    methods (Access = public)
    
        %-- Constructor
        function h = linear_scan(input_x,input_z)
            
            %-- Constructor of linear_scan class
            %-- Syntax:
            %-- h = linear_scan(x_vector,z_vector)
            %-- x_vector: Vector defining the x coordinates of each row of pixels
            %-- z_vector: Vector defining the z coordinates of each column of pixels
            
            if nargin>0
                h.x_axis=input_x;
            end
            if nargin>1
                h.z_axis=input_z;
            end
            
        end
        
        %-- lateral_distance
        function xd = lateral_distance(h,x0,z0,steer_angle)
            
            %-- Calculates the lateral distance from the center of
            %-- the apodization window for a specific scanning mode
            %-- Syntax:
            %-- h = lateral_distance(element_position,steering_angle)
            %-- x0: Vector containing the x coordinates of the probe elements (either real or virtual) [m]
            %-- z0 : Vector containing the x coordinates of the probe elements (either real or virtual) [m]
            %-- steering_angle: Steerin angle [rad]
            xd = abs(x0-h.x+h.z*tan(steer_angle));
        end
    
        %-- Spatial step in the beam direction
        function dz = depth_step(h)
            
            %-- Calculates the spatial step for a given scan area
            %-- Syntax:
            %-- dz = depth_step()
            %-- dz: Spatial step in the beam direction
            dz = mean(diff(h.z_axis)); %-- spatial step in the beam direction
        end
        
    end
    
    %-- Set methods
    methods
        
        function h = set.x_axis(h,input_vector)
            assert(size(input_vector,1)>size(input_vector,2), 'The x vector must be a column vector!')
            h.x_axis = input_vector;
            h.dx = h.x_axis(2)-h.x_axis(1);
            h.Nx = numel(h.x_axis);
            [h.x_matrix,h.z_matrix] = meshgrid(h.x_axis,h.z_axis); 
            h.x = h.x_matrix(:);
            h.z = h.z_matrix(:);
            h.pixels = length(h.x);
        end
        
        function h = set.z_axis(h,input_vector)
            assert(size(input_vector,1)>size(input_vector,2), 'The z vector must be a column vector!')
            h.z_axis = input_vector;
            h.dz = h.z_axis(2)-h.z_axis(1);
            h.Nz = numel(h.z_axis);
            [h.x_matrix,h.z_matrix] = meshgrid(h.x_axis,h.z_axis); 
            h.x = h.x_matrix(:);
            h.z = h.z_matrix(:);
            h.pixels = length(h.x);
        end
        
    end
    
    %-- HDF5 file management
    methods (Access = public)
        
        function write_file_hdf5(h,filename)
            
            %-- write HUFF version in the root group
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
            H5T.enum_insert(filetype, 'US', 0); 
            H5T.enum_insert(filetype, 'SR', 1); 
            gid = H5G.open(file,location);
            space = H5S.create_simple(1,1,[]);

            attr = H5A.create (gid,'type',filetype,space,'H5P_DEFAULT');
            H5A.write (attr,filetype, uint32(0));  % <--- US
            H5A.close (attr);
            H5G.close(gid);    
            H5S.close(space);
            H5T.close(filetype);
            H5F.close(file);
            
            %-- us_dataset subtype
            file = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
            filetype = H5T.enum_create('H5T_NATIVE_INT');
            H5T.enum_insert (filetype,'STA', 0);
            H5T.enum_insert (filetype,'CPW', 1);
            H5T.enum_insert (filetype,'VS', 2);
            H5T.enum_insert (filetype,'BS', 3);
            gid = H5G.open(file,location);
            space = H5S.create_simple (1,1,[]);
            attr = H5A.create (gid,'subtype',filetype,space,'H5P_DEFAULT');            
            H5A.write (attr,filetype,uint32(1)); % <---- TYPE CPWC
                        
            %-- Signal format 
            file = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
            filetype = H5T.enum_create('H5T_NATIVE_INT');
            H5T.enum_insert (filetype,'RF',0); 
            H5T.enum_insert (filetype,'IQ',1); 
            gid = H5G.open(file,location);
            space = H5S.create_simple(1,1,[]);

            attr = H5A.create (gid, 'signal_format', filetype, space, 'H5P_DEFAULT');
            H5A.write(attr, filetype, uint32(1));  % <--- IQ
            H5A.close(attr);
            H5G.close(gid);    
            H5S.close(space);
            H5T.close(filetype);
            H5F.close(file);

            %-- Common attributes
      
            %-- add x-axis
            dset_details.Location = location;
            dset_details.Name = 'x_axis';
            hdf5write(filename, dset_details, h.x_axis, 'WriteMode', 'append');

            %-- add z-axis
            dset_details.Location = location;
            dset_details.Name = 'z_axis';
            hdf5write(filename, dset_details, h.z_axis, 'WriteMode', 'append');
            
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
           
            %-- read US metagroup
            info = h5info(filename,'/US');

            %-- read the groups in the metagroup
            for n=1:length(info.Groups)                
                location = info.Groups(n).Name;
                dstype = h5readatt(filename,location,'type');
                if strcmp(dstype,'US')
                    subtype = h5readatt(filename,location,'subtype');
                    if strcmp(subtype{1},'CPW')                        
                        %-- subtype
                        dataset_subtype = h5readatt(filename,location,'subtype');
                        assert(strcmp(dataset_subtype,'CPW'),'Only CPWC us_dataset are supported!');                        
                        %-- read signal format 
                        signal_format = h5readatt(filename,location,'signal_format');
                        switch(signal_format{1})
                            case 'RF'
                                error('RF format not available!');
                            case 'IQ'
                                ;
                            otherwise
                                error('Unknown signal format!');
                        end
                        %-- Data
                        %-- read x_axis
                        h.x_axis = h5read(filename,[location '/x_axis']);
                        %-- read z_axis
                        h.z_axis = h5read(filename,[location '/z_axis']);                                                
                    end                    
                end                
            end
            
        end

        function read_file_mat(h,filename)

            %-- Load mat file
            load(filename);
            
            %-- Data

            %-- read x_axis
            h.x_axis = PARAM.x_axis;

            %-- read z_axis
            h.z_axis = PARAM.z_axis;
        
        end
        
        function write_file_mat(h,filename)
           
            %-- Common attributes
      
            %-- add x-axis            
            PARAM.x_axis = h.x_axis;

            %-- add z-axis
            PARAM.z_axis = h.z_axis;

            %-- write mat file
            save(filename,'PARAM');
            
        end
        
    end
    
end

