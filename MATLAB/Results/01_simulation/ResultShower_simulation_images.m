% This script plots Fig.2 (a-b) in the paper submitted to DGM4MICCAI
% -- Baseline vs. DRUS vs. WDRUS vs. Ground Truth
% -- on the sythetic datasets fetus & SynVitro


close all
clear
clc
addpath(genpath('../PICMUS'));


% ---------------------- set parameters ----------------------
phantomIdx = 3;                    % 1-kidney || 2-fetus || 3-simuComplex
gamma = {0.3, 0.7, 1, 1.5, 2, 2.5, 3, 3.5};
NoiseLevels = [1,3,6];             % for selecting a subset of gamma levels
numNoiseLevels = length(NoiseLevels);

model = {'DAS','Hty_50it','CHty_50it'};
lm = length(model);
dynamicRange = 60;
J1 = customcolormap([0, 0.4, 0.5, 0.6, 1], [1 1 1; 1 0 0; 0 0 0;0 0 1; 1 1 1]);
J2 = customcolormap([0, 1], [1 1 1; 0 0 0]);



% ------------------------ load image data -----------------------

% load the ground Truth x_sign_linear
pathtrue = ['01_simulation/SimulationResults/',num2str(phantomIdx),'/']; %addpath(path)
load([pathtrue 'o_orig_' num2str(phantomIdx) '.mat'])
x_sign_linear = reshape(single(x_orig),[256,256]);


% x_sign_dB_abs
x_sign_dB_abs = 20*log10(abs(x_sign_linear)./max(abs(x_sign_linear(:))))+dynamicRange;
x_sign_dB_abs(x_sign_dB_abs<0) = 0;
x_sign_dB_abs = x_sign_dB_abs ./ dynamicRange;


% initialization
X_dB_abs = cell(length(gamma),lm);

% load each restored image
for i = 1:lm
    for idx = 1: numNoiseLevels
        j = NoiseLevels(idx);
        pathus = [pathtrue,'us_',model{i} '/']; 
        load([pathus num2str(j) '_-1.mat'])
        if strcmp(model{i},'DAS')
            x = o_Hty(1:65536,:);
            x = reshape(x,[256,256]);
        else
            x1 = x(1,:,:); 
            x2 = x(2,:,:); 
            x3 = x(3,:,:); 
            x = (x1+x2+x3) ./ 3;
            x = squeeze(x);
        end
        x = x./max(abs(x(:)));

        % --x_sign_dB_abs
        temp = 20*log10(abs(x)) ./ dynamicRange + 1;
        temp(temp<0) = 0;
        X_dB_abs{j,i} = temp;
             
    end
end

% ---------------------- plot the images ------------------------

figure('Position', [10 10 (numNoiseLevels+1)*220+100 lm*220]); 
h = tiledlayout(lm,numNoiseLevels+1,'TileSpacing','none','Padding','none');
nexttile
if phantomIdx == 3
imagesc(1:256,1:256,x_sign_dB_abs*(60)-60);colormap(J1); caxis([-60,0]); axis equal manual; axis off;
for i = 1: lm
    for idx = 1: numNoiseLevels
        j = NoiseLevels(idx);
        nexttile((i-1)*(numNoiseLevels+1)+idx+1)
        imagesc(1:256,1:256,X_dB_abs{j,i}*(60)-60); colormap(J2); caxis([-60,0]);
        axis equal manual;axis off;
    end
end
elseif ismember(phantomIdx, [1,2])
imagesc(1:256,1:256,x_sign_dB_abs);colormap(J1); caxis([0,1]); axis equal manual; axis off;
for i = 1: lm
    for idx = 1: numNoiseLevels
        j = NoiseLevels(idx);
        nexttile((i-1)*(numNoiseLevels+1)+idx+1)
        imagesc(1:256,1:256,X_dB_abs{j,i}); colormap(J2); caxis([0,1]);
        axis equal manual;axis off;
    end
end
end

colorbar('Location','eastoutside','Position',[0.13 0.1 0.02 0.52],'FontSize',20)



% ---------------------- add annotations ----------------------

h.InnerPosition=[0.005 0.01 0.84 0.908];
name = {'Ground Truth','\gamma='+string(gamma(NoiseLevels(1))), '\gamma='+string(gamma(NoiseLevels(2))), '\gamma='+string(gamma(NoiseLevels(3)))};
for i = 1: 4
    annstr = (name{i}); % annotation text
    annpos = [(i-1)*0.21+0.022 0.94 0.2,0.05]; % annotation position in figure coordinates
    ha = annotation('textbox',annpos,'string',annstr);
    ha.HorizontalAlignment = 'center';
    ha.BackgroundColor = 'none'; % make the box opaque with some color
    ha.EdgeColor = 'none';
    ha.FontSize = 20;
end

name = {'WDRUS', 'DRUS','Baseline'};
for i = 1: 3
    annstr = (name{i}); % annotation text
    annpos = [ 0.8, (i-1)*0.35+0.1 0.2,0.05]; % annotation position in figure coordinates
    ha = annotation('textbox',annpos,'string',annstr);
    ha.HorizontalAlignment = 'center';
    ha.BackgroundColor = 'none'; % make the box opaque with some color
    ha.EdgeColor = 'none';
    ha.FontSize = 20;
end