function [dataset, dasSaveName] = datasetReader(acquisition_type,phantom_type,data_type)

%Parameters
% acquisition_type       %-- 1 = simulation || 2 = experiments || 3 = vivo
% phantom_type           %-- 1 = resolution & distorsion || 2 = contrast & speckle quality || 3 = cross || 4 = long 
% data_type              %-- 1 = IQ || 2 = RF


clc;
addpath(genpath('../PICMUS/code/src'));

%-- Parsing parameter choices
switch acquisition_type    
    case 1
        acquisition = 'simulation';
        acqui = 'simu';
    case 2
        acquisition = 'experiments';
        acqui = 'expe';
    case 3
        acquisition = 'in_vivo';
        acqui = 'expe';       
end

switch phantom_type    
    case 1
        phantom = 'resolution_distorsion';
        phan = 'reso';
    case 2
        phantom = 'contrast_speckle';
        phan = 'cont';
    case 3
        phantom = 'carotid_cross';
        phan = 'cross';
    case 4
        phantom = 'carotid_long';
        phan = 'long';
end

switch data_type    
    case 1
        data = 'iq';
    case 2
        data = 'rf';      
end


%-- Create path to load corresponding files
path_dataset = ['../PICMUS/database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_dataset_',data,'.hdf5'];


%-- Read the corresponding dataset and the region where to reconstruct the image
dataset = us_dataset();
dataset.read_file(path_dataset);
dasSaveName = [acqui '_' phan];
end