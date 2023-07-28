
% This script gives the values in Tab.1 in the paper submitted to DGM4MICCAI
% -- DAS1 vs. DAS11 vs. DAS75 vs. DRUS vs. WDRUS
% -- on the picmus datasets SR, SC, ER, EC

close all
clear
clc
addpath(genpath('../PICMUS'))

% parameters
flag_display = 0;
phantoms = [1,2,3,4];     % 1-SR || 2-SC || 3-ER  || 4-EC
methods = {'DAS1', 'DAS11', 'DAS75', 'BH50', 'CBH50', 'BH50_fineTune', 'CBH50_fineTune'}; %'DAS1', 'DAS11', 'DAS75', 
pathPicmusResults = '02_picmus/'; 

dynamicRange = 60;
scan = linear_scan(linspace(-0.018,0.018,256).', linspace(0.01,0.036+0.01,256).');
numPhantoms  = length(phantoms);
numMethods = length(methods);

% Initialization
ImageSet = cell(numPhantoms,numMethods);
FWHMAs = zeros(numPhantoms,numMethods);
FWHMLs = zeros(numPhantoms,numMethods);
CNRs   = zeros(numPhantoms,numMethods);
gCNRs  = zeros(numPhantoms,numMethods);
SNRs   = zeros(numPhantoms,numMethods);
KSs    = zeros(numPhantoms,numMethods);

for phantomIdx = phantoms   

    path_phantom = [pathPicmusResults 'phantoms/picmus_phantom_' num2str(phantomIdx) '.hdf5'];


    % load the DAS images (DAS1, 11, 75)
    methodIdx = 1;
    for i = [1,3,4] 
        ImageSet{phantomIdx, methodIdx} = us_image();
        ImageSet{phantomIdx, methodIdx}.read_file([pathPicmusResults 'DAS/',num2str(phantomIdx),'.hdf5']);
        ImageSet{phantomIdx, methodIdx}.number_plane_waves = ImageSet{phantomIdx, methodIdx}.number_plane_waves(i);
        ImageSet{phantomIdx, methodIdx}.data = ImageSet{phantomIdx, methodIdx}.data(:,:,i);
        methodIdx = methodIdx+1;
    end
    
    % load the restored images
    for i = methodIdx : numMethods
        ImageSet{phantomIdx, i} = us_image(methods{i});
        ImageSet{phantomIdx, i}.scan = scan;
        ImageSet{phantomIdx, i}.number_plane_waves=1;
        if ismember(phantomIdx, [2,4])
            ImageSet{phantomIdx, i}.postenv = 1;
        end
        if strcmp(methods{i}(1),'B')
            currentFolder = 'BH';
        elseif strcmp(methods{i}(1), 'C')
            currentFolder = 'CBH';
        end
        load([pathPicmusResults currentFolder '/Results/' methods{i} '/' num2str(phantomIdx) '_-1.mat']);
        x1 = x(1,:,:); 
        x2 = x(2,:,:); 
        x3 = x(3,:,:); 
        x = (x1+x2+x3) ./ 3;
        x = squeeze(x);
        ImageSet{phantomIdx, i}.data = reshape(abs(x),scan.Nz, scan.Nx);  
    end
    
    
    % Evaluation
    for i = 1: numMethods
        [FWHMAs(phantomIdx,i), FWHMLs(phantomIdx,i), CNRs(phantomIdx,i), gCNRs(phantomIdx,i), SNRs(phantomIdx,i), KSs(phantomIdx,i)] = evaluation(path_phantom, ImageSet{phantomIdx, i}, flag_display);
    end

end

clearvars -except FWHMAs FWHMLs CNRs gCNRs SNRs KSs

