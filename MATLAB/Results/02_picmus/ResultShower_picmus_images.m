% This script plots the Fig.3 in the paper submitted to DGM4MICCAI
% -- DAS1 vs. DAS11 vs. DAS75 vs. DRUS vs. WDRUS
% -- on the picmus datasets SR, SC, ER, EC

close all
clear
clc
addpath(genpath('../PICMUS/code/src'));


% parameters
numPhan = 4;
model = {'DAS1','DAS11', 'DAS75', 'BH50','CBH50', 'BH50_fineTune','CBH50_fineTune'};
numModel = length(model);
scan = linear_scan(linspace(-0.018,0.018,256).', linspace(0.01,0.036+0.01,256).');


% initialization
X = cell(numPhan,numModel);
x_lim = [min(scan.x_matrix(:)) max(scan.x_matrix(:))]*1e3; 
z_lim = [min(scan.z_matrix(:)) max(scan.z_matrix(:))]*1e3;


% load and quantitize each restored image
for i = 1:numModel
    for j = 1:numPhan
        if strcmp(model{i},'DAS1')
            img = us_image;
            img.read_file_hdf5(['02_picmus/DAS/', num2str(j), '.hdf5'])
            x = img.data(:,:,1);
        elseif strcmp(model{i},'DAS11')
            img = us_image;
            img.read_file_hdf5(['02_picmus/DAS/', num2str(j), '.hdf5'])
            x = img.data(:,:,3);  
        elseif strcmp(model{i},'DAS75')
            img = us_image;
            img.read_file_hdf5(['02_picmus/DAS/', num2str(j), '.hdf5'])
            x = img.data(:,:,4);
        elseif (strcmp(model{i},'BH50_fineTune') || strcmp(model{i},'BH50'))
            phanIdx = num2str(j); 
            pathus = ['02_picmus/BH/Results/', model{i} ,'/']; 
            load([pathus phanIdx '_-1.mat'])
            x1 = x(1,:,:); 
            x2 = x(2,:,:); 
            x3 = x(3,:,:); 
            x = (x1+x2+x3) ./ 3;
            x = squeeze(x);
        elseif (strcmp(model{i},'CBH50_fineTune') || strcmp(model{i},'CBH50'))
            phanIdx = num2str(j); 
            pathus = ['02_picmus/CBH/Results/', model{i} ,'/']; 
            load([pathus phanIdx '_-1.mat'])
            x1 = x(1,:,:); 
            x2 = x(2,:,:); 
            x3 = x(3,:,:); 
            x = (x1+x2+x3) ./ 3;
            x = squeeze(x);
        end
        X{j,i} = x;
    end
end


% plot the reconstructed images
figure('Position', [10 10 numModel*200 numPhan*190]); 
h = tiledlayout(numPhan,numModel,'TileSpacing','none','Padding','none');
for j = 1: numPhan
    for i = 1:  numModel
        nexttile
        if ((~any(i == [1,2,3])) && (~any(j == [1,3])))
            Image_realScale_envelope(X{j,i}, '', scan)
            axis off; colorbar off;
        else
            Image_realScale(X{j,i}, '', scan)
            axis off; colorbar off;
        end
    end
end
h.InnerPosition=[0.035 0.00 0.89 0.925];


% add annotations
name = {'Before fine-tuning','After fine-tuning'};
for i = 1: length(name)
    annstr = (name{i}); % annotation text
    annpos = [(i-1)*0.255+0.45 0.96 0.2,0.05]; % annotation position in figure coordinates
    ha = annotation('textbox',annpos,'string',annstr,'Interpreter','latex');
    ha.HorizontalAlignment = 'center';
    ha.BackgroundColor = 'none'; % make the box opaque with some color
    ha.EdgeColor = 'none';
    ha.FontSize = 25;
end

name = {'DAS1','DAS11', 'DAS75', 'DRUS','WDRUS', 'DRUS','WDRUS'};
%name = {'DAS 1PW','DAS 11PWs', 'DAS 75PWs', 'DRUS 1PW','WDRUS 1PW', 'DRUS 1PW','WDRUS 1PW'};
for i = 1: numModel
    annstr = (name{i}); % annotation text
    annpos = [(i-1)*0.126+0.028 0.92 0.15,0.05]; % annotation position in figure coordinates
    ha = annotation('textbox',annpos,'string',annstr,'Interpreter','latex');
    ha.HorizontalAlignment = 'center';
    ha.BackgroundColor = 'none'; % make the box opaque with some color
    ha.EdgeColor = 'none';
    ha.FontSize = 25;
end
colorbar('Location','eastoutside','Position',[1.02 0.045 0.02 0.99],'FontSize',20)


name = {'EC', 'ER','SC', 'SR'};
for i = 1: length(name)
    annstr = (name{i}); % annotation text
    annpos = [ 0.0001, (i-1)*0.24+0.1 0.03,0.05]; % annotation position in figure coordinates
    ha = annotation('textbox',annpos,'string',annstr);
    ha.HorizontalAlignment = 'center';
    ha.BackgroundColor = 'none'; % make the box opaque with some color
    ha.EdgeColor = 'none';
    ha.FontSize = 25;
end
