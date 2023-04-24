%% picmus single image
clc
step = '50'; model = 'CBH'; phanIdx = 4;

folder = [model,step];
path = ['picmus/', model ,'/Results/' folder '/']; 
phanIdx = num2str(phanIdx); 
load([path phanIdx '_-1.mat'])  
x1 = x(1,:,:); 
x2 = x(2,:,:); 
x3 = x(3,:,:); 
x = (x1+x2+x3) ./ 3;
x = squeeze(x);

addpath(genpath('../PICMUS/code/src'));
scan = linear_scan(linspace(-0.018,0.018,256).', linspace(0.01,0.036+0.01,256).');
figure,
subplot(2,1,1);Image_realScale_envelope(x,'with envelope', scan); 
subplot(2,1,2);Image_realScale(x,'no envelope', scan)


%% picmus summary
numPhan = 4;
model = {'DAS1','DAS11', 'DAS75', 'BH30','BH50', 'CBH30','CBH50'};
numModel = length(model);
addpath(genpath('../PICMUS/code/src'));
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
            img.read_file_hdf5(['picmus/Evaluation/database/DAS/', num2str(j), '.hdf5'])
            x = img.data(:,:,1);
        elseif strcmp(model{i},'DAS11')
            img = us_image;
            img.read_file_hdf5(['picmus/Evaluation/database/DAS/', num2str(j), '.hdf5'])
            x = img.data(:,:,3);  
        elseif strcmp(model{i},'DAS75')
            img = us_image;
            img.read_file_hdf5(['picmus/Evaluation/database/DAS/', num2str(j), '.hdf5'])
            x = img.data(:,:,4);
        elseif (strcmp(model{i},'BH30') || strcmp(model{i},'BH50'))
            phanIdx = num2str(j); 
            pathus = ['picmus/BH/Results/', model{i} ,'/']; 
            load([pathus phanIdx '_-1.mat'])
            x1 = x(1,:,:); 
            x2 = x(2,:,:); 
            x3 = x(3,:,:); 
            x = (x1+x2+x3) ./ 3;
            x = squeeze(x);
        elseif (strcmp(model{i},'CBH30') || strcmp(model{i},'CBH50'))
            phanIdx = num2str(j); 
            pathus = ['picmus/CBH/Results/', model{i} ,'/']; 
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

figure('Position', [10 10 numModel*250 numPhan*300]); 
h = tiledlayout(numPhan,numModel,'TileSpacing','tight','Padding','compact');
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
h.InnerPosition=[0.02 0.025 0.9091 0.9075];

name = {'DAS 1PW','DAS 11PWs', 'DAS 75PWs', 'DRUS 1PW $it$=30','DRUS 1PW $it$=50', 'WDRUS 1PW $it$=30','WDRUS 1PW $it$=50'};
for i = 1: numModel
    annstr = (name{i}); % annotation text
    annpos = [(i-1)*0.133+0.025 0.94 0.111,0.05]; % annotation position in figure coordinates
    ha = annotation('textbox',annpos,'string',annstr,'Interpreter','latex');
    ha.HorizontalAlignment = 'center';
    ha.BackgroundColor = 'none'; % make the box opaque with some color
    ha.EdgeColor = 'none';
    ha.FontSize = 20;
end
colorbar('Location','eastoutside','Position',[0.97 0.06 0.008 0.9],'FontSize',20)



