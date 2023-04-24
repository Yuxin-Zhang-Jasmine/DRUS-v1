

% compare CHtH and HtH on the simulated datasets
close all
clear
clc

% parameters
folderName = 'medical_Bivertical'; 
phantomIdx = 2; 
gamma = [0.3, 0.7, 1, 1.5, 2, 2.5];%, 3, 3.5];
NoiseLevels = [1,2,3,4,5,6];%,7,8]; 
numNoiseLevels = length(NoiseLevels);
model = {'DAS','CHty','Hty_20it'};
lm = length(model);
name = {'DAS','WDRUS ','DRUS'};

dynamicRange = 60;
J1 = customcolormap([0, 0.4, 0.5, 0.6, 1], [1 1 1; 1 0 0; 0 0 0;0 0 1; 1 1 1]);
J2 = customcolormap([0, 1], [1 1 1; 0 0 0]);

% load the ground Truth x_sign_linear
pathtrue = ['simulation/SimulationResults/' folderName '/',num2str(phantomIdx),'/']; %addpath(path)
load([pathtrue 'o_orig_' num2str(phantomIdx) '.mat'])
x_sign_linear = reshape(single(x_orig),[256,256]);

% --x_sign_dB
x_signPos = x_sign_linear;
x_signPos(x_signPos<0) = 0;
x_signPos = 20*log10(abs(x_signPos)./max(abs(x_sign_linear(:))))+dynamicRange;
x_signPos(x_signPos<0) = 0;
x_signPos = x_signPos ./ dynamicRange;

x_signNeg = x_sign_linear;
x_signNeg(x_signNeg>=0) = 0;
x_signNeg = 20*log10(abs(x_signNeg)./max(abs(x_sign_linear(:))))+dynamicRange;
x_signNeg(x_signNeg<0) = 0;
x_signNeg = x_signNeg ./ dynamicRange;

x_sign_dB = -x_signNeg + x_signPos;

% --x_sign_dB_abs
x_sign_dB_abs = 20*log10(abs(x_sign_linear)./max(abs(x_sign_linear(:))))+dynamicRange;
x_sign_dB_abs(x_sign_dB_abs<0) = 0;
x_sign_dB_abs = x_sign_dB_abs ./ dynamicRange;


% initialization
X_linear = cell(numNoiseLevels,lm);
MAE_linear = zeros(numNoiseLevels,lm);
MSE_linear = zeros(numNoiseLevels,lm);

X_dB = cell(numNoiseLevels,lm);
MAE_dB = zeros(numNoiseLevels,lm);
MSE_dB = zeros(numNoiseLevels,lm);

X_dB_abs = cell(numNoiseLevels,lm);
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
        
        % --x_sign_linear
        X_linear{j,i} = x;
        MAE_linear(j,i) = sum(abs(x(:)-x_sign_linear(:)))/numel(x);
        MSE_linear(j,i) = sum((x(:)-x_sign_linear(:)).^2)/numel(x);
        
        % --x_sign_dB
        x_signPos = x;
        x_signPos(x_signPos<0) = 0;
        x_signPos = 20*log10(abs(x_signPos)./max(abs(x(:))))+dynamicRange;
        x_signPos(x_signPos<0) = 0;
        x_signPos = x_signPos ./ dynamicRange;
        
        x_signNeg = x;
        x_signNeg(x_signNeg>=0) = 0;
        x_signNeg = 20*log10(abs(x_signNeg)./max(abs(x(:))))+dynamicRange;
        x_signNeg(x_signNeg<0) = 0;
        x_signNeg = x_signNeg ./ dynamicRange;
        
        temp = -x_signNeg + x_signPos;
        X_dB{j,i} = temp;
        MAE_dB(j,i) = sum(abs(temp(:)-x_sign_dB(:)))/numel(x);
        MSE_dB(j,i) = sum((temp(:)-x_sign_dB(:)).^2)/numel(x);


        % --x_sign_dB_abs
        temp = 20*log10(abs(x)) ./ dynamicRange + 1;
        temp(temp<0) = 0;
        X_dB_abs{j,i} = temp;
        SSIM_dB_abs(j,i) = ssim(temp,x_sign_dB_abs);
        PSNR_dB_abs(j,i) = psnr(temp,x_sign_dB_abs);
             
    end
end

fontSize = 20;
figure('Position',  [100, 100, 2000, 500]); 
t = tiledlayout(1,4, 'Padding', 'compact', 'TileSpacing', 'compact'); 
nexttile    
plot(gamma,SSIM_dB_abs(NoiseLevels,:),'LineWidth',1.5); title('SSIM(fetus)'); 
xlabel('\gamma');xticks(gamma);xlim([min(gamma),max(gamma)]);
legend(name,'Location','northeast');
set(gca,'fontsize',fontSize);
nexttile
plot(gamma,PSNR_dB_abs(NoiseLevels,:),'LineWidth',1.5); title('PSNR(fetus)'); xlabel('\gamma');
xticks(gamma);xlim([min(gamma),max(gamma)]);
%legend(name,'Location','northeast'); 
set(gca,'fontsize',fontSize);


% parameters
folderName = 'medical_Bivertical'; 
phantomIdx = 1; 
gamma = [0.3, 0.7, 1, 1.5, 2, 2.5];%, 3, 3.5];
NoiseLevels = [1,2,3,4,5,6];%,7,8]; 
numNoiseLevels = length(NoiseLevels);
model = {'DAS','CHty','Hty_20it'};
lm = length(model);


dynamicRange = 60;
J1 = customcolormap([0, 0.4, 0.5, 0.6, 1], [1 1 1; 1 0 0; 0 0 0;0 0 1; 1 1 1]);
J2 = customcolormap([0, 1], [1 1 1; 0 0 0]);

% load the ground Truth x_sign_linear
pathtrue = ['simulation/SimulationResults/' folderName '/',num2str(phantomIdx),'/']; %addpath(path)
load([pathtrue 'o_orig_' num2str(phantomIdx) '.mat'])
x_sign_linear = reshape(single(x_orig),[256,256]);

% --x_sign_dB
x_signPos = x_sign_linear;
x_signPos(x_signPos<0) = 0;
x_signPos = 20*log10(abs(x_signPos)./max(abs(x_sign_linear(:))))+dynamicRange;
x_signPos(x_signPos<0) = 0;
x_signPos = x_signPos ./ dynamicRange;

x_signNeg = x_sign_linear;
x_signNeg(x_signNeg>=0) = 0;
x_signNeg = 20*log10(abs(x_signNeg)./max(abs(x_sign_linear(:))))+dynamicRange;
x_signNeg(x_signNeg<0) = 0;
x_signNeg = x_signNeg ./ dynamicRange;

x_sign_dB = -x_signNeg + x_signPos;

% --x_sign_dB_abs
x_sign_dB_abs = 20*log10(abs(x_sign_linear)./max(abs(x_sign_linear(:))))+dynamicRange;
x_sign_dB_abs(x_sign_dB_abs<0) = 0;
x_sign_dB_abs = x_sign_dB_abs ./ dynamicRange;


% initialization
X_linear = cell(numNoiseLevels,lm);
MAE_linear = zeros(numNoiseLevels,lm);
MSE_linear = zeros(numNoiseLevels,lm);

X_dB = cell(numNoiseLevels,lm);
MAE_dB = zeros(numNoiseLevels,lm);
MSE_dB = zeros(numNoiseLevels,lm);

X_dB_abs = cell(numNoiseLevels,lm);
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
        
        % --x_sign_linear
        X_linear{j,i} = x;
        MAE_linear(j,i) = sum(abs(x(:)-x_sign_linear(:)))/numel(x);
        MSE_linear(j,i) = sum((x(:)-x_sign_linear(:)).^2)/numel(x);
        
        % --x_sign_dB
        x_signPos = x;
        x_signPos(x_signPos<0) = 0;
        x_signPos = 20*log10(abs(x_signPos)./max(abs(x(:))))+dynamicRange;
        x_signPos(x_signPos<0) = 0;
        x_signPos = x_signPos ./ dynamicRange;
        
        x_signNeg = x;
        x_signNeg(x_signNeg>=0) = 0;
        x_signNeg = 20*log10(abs(x_signNeg)./max(abs(x(:))))+dynamicRange;
        x_signNeg(x_signNeg<0) = 0;
        x_signNeg = x_signNeg ./ dynamicRange;
        
        temp = -x_signNeg + x_signPos;
        X_dB{j,i} = temp;
        MAE_dB(j,i) = sum(abs(temp(:)-x_sign_dB(:)))/numel(x);
        MSE_dB(j,i) = sum((temp(:)-x_sign_dB(:)).^2)/numel(x);


        % --x_sign_dB_abs
        temp = 20*log10(abs(x)) ./ dynamicRange + 1;
        temp(temp<0) = 0;
        X_dB_abs{j,i} = temp;
        SSIM_dB_abs(j,i) = ssim(temp,x_sign_dB_abs);
        PSNR_dB_abs(j,i) = psnr(temp,x_sign_dB_abs);
             
    end
end

nexttile    
plot(gamma,SSIM_dB_abs(NoiseLevels,:),'LineWidth',1.5); title('SSIM(kidney)'); xlabel('\gamma');
xlim([min(gamma),max(gamma)]);xticks(gamma);
%legend(name,'Location','northeast');
set(gca,'fontsize',fontSize);
nexttile
plot(gamma,PSNR_dB_abs(NoiseLevels,:),'LineWidth',1.5); title('PSNR(kidney)'); xlabel('\gamma');
xticks(gamma);xlim([min(gamma),max(gamma)]);
%legend(name,'Location','northeast');
set(gca,'fontsize',fontSize);