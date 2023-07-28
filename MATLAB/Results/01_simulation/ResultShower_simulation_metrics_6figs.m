% This script plots Fig.2 (c) in the paper submitted to DGM4MICCAI
% -- Baseline vs. DRUS vs. WDRUS vs. Ground Truth
% -- on the sythetic datasets fetus & SynVitro

close all
clear
clc
addpath(genpath('../PICMUS'));

% -------------------- for the phantom fetus-------------------------

% parameters
phantomIdx = 2;       
gamma = [0.3, 0.7, 1, 1.5, 2, 2.5];%, 3, 3.5];
NoiseLevels = [1,2,3,4,5,6];%,7,8]; 
numNoiseLevels = length(NoiseLevels);
model = {'DAS','Hty_50it','CHty_50it'};
lm = length(model);
name = {'DAS','DRUS ','WDRUS'};
dynamicRange = 60;


% load the ground Truth x_sign_linear
pathtrue = ['01_simulation/SimulationResults/',num2str(phantomIdx),'/']; %addpath(path)
load([pathtrue 'o_orig_' num2str(phantomIdx) '.mat'])
gt_sign_linear = reshape(single(x_orig),[256,256]);

% --gt_dB_abs
gt_dB_abs = 20*log10(abs(gt_sign_linear)./max(abs(gt_sign_linear(:))))+dynamicRange;
gt_dB_abs(gt_dB_abs<0) = 0;
gt_dB_abs = gt_dB_abs ./ dynamicRange;


% initialization
SSIM_dB_abs = zeros(numNoiseLevels,lm);
PSNR_dB_abs = zeros(numNoiseLevels,lm);


% load and quantitize each restored image
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
        
        % --gt_dB_abs
        temp = 20*log10(abs(x)) ./ dynamicRange + 1;
        temp(temp<0) = 0;
        SSIM_dB_abs(j,i) = ssim(temp,gt_dB_abs);
        PSNR_dB_abs(j,i) = psnr(temp,gt_dB_abs);
             
    end
end


% -------------------- for the phantom SynVitro-------------------------
phantomIdx = 3; 

flag_display = 0;
path_phantom = '01_simulation/SimulationResults/3/phantom_3.hdf5';
scan = linear_scan(linspace(-0.018,0.018,256).', linspace(0.01,0.036+0.01,256).');
image = us_image(); image.scan = scan; image.number_plane_waves=1;
dynamicRange = 60;


% load the ground Truth x_sign_linear
pathtrue = ['01_simulation/SimulationResults/',num2str(phantomIdx),'/']; %addpath(path)
load([pathtrue 'o_orig_' num2str(phantomIdx) '.mat'])
gt_sign_linear = reshape(single(x_orig),[256,256]);


% metrics of the Groud Truth
image.data = reshape(abs(gt_sign_linear), scan.Nz, scan.Nx);
[gtFWHMA, gtFWHML, gtCNR, gtgCNR, ~, ~] = evaluation(path_phantom, image, flag_display);
        

% initialization
FWHMA = zeros(numNoiseLevels,lm+1); FWHMA(:,end) = gtFWHMA;
FWHML = zeros(numNoiseLevels,lm+1); FWHML(:,end) = gtFWHML;
CNR = zeros(numNoiseLevels,lm+1);   CNR(:,end) = gtCNR;
gCNR = zeros(numNoiseLevels,lm+1);  gCNR(:,end) = gtgCNR;
%SNR = zeros(numNoiseLevels,lm+1);   SNR(:,end) = gtSNR;


% load and quantitize each restored image
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
        image.data = reshape(abs(x),scan.Nz, scan.Nx);
        
        % --evaluation
        [FWHMA(j,i), FWHML(j,i), CNR(j,i), gCNR(j,i), ~, ~] = evaluation(path_phantom, image, flag_display);
             
    end
end


% ------------------- plot the figures of scores --------------------

fontSize = 18; linewidth = 3;
figure('Position',  [100, 100, 2000, 350]); 
t = tiledlayout(1,6, 'Padding', 'compact', 'TileSpacing', 'compact', 'Position',[0.035,0.33,0.95,0.55]); 

lengedName = {'Baseline','DRUS ','WDRUS','Ground Truth'};
marVec = {"o", '+','*','^'};

nexttile    
for i  = 1:  length(lengedName)
plot(gamma,FWHMA(NoiseLevels,i),'Marker',marVec{i},'LineWidth',linewidth, 'MarkerSize',10); 
hold on
end
hold off
title('FWHM_{axial} [mm]'); xlabel('\gamma');
xlim([min(gamma),max(gamma)]);%xticks(gamma);
set(gca,'fontsize',fontSize);

nexttile
for i  = 1:  length(lengedName)
plot(gamma,FWHML(NoiseLevels,i),'Marker',marVec{i},'LineWidth',linewidth, 'MarkerSize',10); 
hold on
end
hold off
title('FWHM_{lateral} [mm]'); xlabel('\gamma');
xlim([min(gamma),max(gamma)]);%xticks(gamma);
set(gca,'fontsize',fontSize);

nexttile
for i  = 1:  length(lengedName)
plot(gamma,gCNR(NoiseLevels,i),'Marker',marVec{i},'LineWidth',linewidth, 'MarkerSize',10); 
hold on
end
hold off
title('gCNR'); xlabel('\gamma');
xlim([min(gamma),max(gamma)]);%xticks(gamma);
set(gca,'fontsize',fontSize);

nexttile
for i  = 1:  length(lengedName)
plot(gamma,CNR(NoiseLevels,i),'Marker',marVec{i},'LineWidth',linewidth, 'MarkerSize',10); 
hold on
end
hold off
title('CNR [dB]'); xlabel('\gamma');
xlim([min(gamma),max(gamma)]);%xticks(gamma);
set(gca,'fontsize',fontSize);
legend(lengedName,'Location','southoutside','Orientation','horizontal','FontSize',fontSize+3);
legend('boxoff');

nexttile    
for i  = 1:  lm
plot(gamma,SSIM_dB_abs(NoiseLevels,i),'Marker',marVec{i},'LineWidth',linewidth, 'MarkerSize',10); 
hold on
end
hold off
title('SSIM(fetus)'); xlabel('\gamma');
xlim([min(gamma),max(gamma)]);%xticks(gamma);
set(gca,'fontsize',fontSize);

nexttile
for i  = 1:  lm
plot(gamma,PSNR_dB_abs(NoiseLevels,i),'Marker',marVec{i},'LineWidth',linewidth, 'MarkerSize',10); 
hold on
end
hold off
title('PSNR(fetus)'); xlabel('\gamma');
xlim([min(gamma),max(gamma)]); %xticks(gamma);
set(gca,'fontsize',fontSize);

