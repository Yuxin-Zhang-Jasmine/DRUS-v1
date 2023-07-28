
close all
clc
addpath(genpath('../PICMUS/code/src'));

% ---------------
% Params
%----------------
RxApodName =  'hanning'; %'tukey25';%
rx_f_number = 0.5;
[dataset,~] = datasetReader(1,1,2);      

% K = dataset.samples;
Ne = size(dataset.probe_geometry,1);
fs = dataset.sampling_frequency;
c = dataset.c0;

% -------------------------
% Waveform PICMUS
%--------------------------
ros = 10; dth = 1/fs/ros;

fe = 5.208e6; BWRe = 0.67; ae = 4*log(2)/(pi*fe*BWRe).^2; Ncycles = 2.5; 
te = (-Ncycles/fe/2: dth :Ncycles/fe/2)';
he = exp(-(te).^2/ae).*cos(2*pi*fe*(te));

f0 = 5.1333e6; BWR0 = 0.65; a0 = 4*log(2)/(pi*f0*BWR0).^2; 
t0 = (-2/f0: dth :2/f0)';
h0 = exp(-(t0).^2/a0).*cos(2*pi*f0*(t0));

one_way_ir = conv(h0, he);
h = conv(one_way_ir, h0);
h = h ./ max(h(:));

Nth = length(h); h = h(fix(Nth/4): ceil(Nth/4*3)); Nth = length(h);
lagt = (Nth/2) *dth ;

% figure, 
% subplot(311), plot(te,he) ;hold on; plot(t0, h0); hold off
% subplot(312), plot(abs(fft(he))); hold on; plot(abs(fft(h0))); hold off
% subplot(313), plot(te,he) ;hold on; plot(t0, h0); plot((0:Nth-1)*dth - lagt, h); grid on;

Kh = fix(Nth/(ros+1)); 
h = reshape(h(1:Kh*(ros+1)), ros+1, Kh)';
h = flip(h,2);


% ----------------------
% Geometry and TimeDelay
%-----------------------
scan = linear_scan(linspace(-0.018,0.018,256).', linspace(0.01,0.036+0.01,256).');
scan_left = linear_scan(linspace(-0.018,scan.x_axis(128),128).', linspace(0.01,0.036+0.01,256).');

% time_vector = (dataset.initial_time + (0:(K-1)) ./ fs).';
% z_axis = time_vector .* c ./ 2;
% id1 = find((z_axis-0.01) < 0.005); id1 = id1(1);
% id2 = find((z_axis<0.025  & z_axis>0.0245); id2 = id2(end);
% scan0 = linear_scan(scan.x_axis,z_axis(id1:id2));

xe = dataset.probe_geometry(:,1).';         % 1 x Ne
x = scan_left.x ;                               % pixels x 1
z = scan_left.z ;                               % pixels x 1

idxi = linspace(1, Ne, 4*Ne);               % 1 x 4Ne
xTi = interp1(xe, idxi, 'spline');          % 1 x 4Ne
dTX = min(sqrt((xTi-x).^2 + z.^2), [], 2);  % pixels x 1
dRX = sqrt((xe-x).^2 + z.^2);               % pixels x Ne

tau = (dTX+dRX)/c - lagt;                   % pixels x Ne  [unit: s]
delay = double(tau * fs + 1);               % pixels x Ne

K = ceil(max(delay(:))) + Kh ;
% ---------------------
% Forward Model (left)
%----------------------
rx_aperture = scan_left.z/rx_f_number;
rx_aperture_distance = abs(scan_left.x*ones(1,dataset.channels)-ones(scan_left.pixels,1)*dataset.probe_geometry(:,1).');
receive_apodization = tools.apodization(rx_aperture_distance,rx_aperture*ones(1,dataset.channels),RxApodName);
receive_apodization = repmat(receive_apodization(:).', [Kh,1]);

tic
di = fix(delay);
dd = delay - di;

r = di+(0:Ne-1)*K;
i = (r(:)+(0:Kh-1))';
j = repmat((1:scan_left.pixels),[Kh,Ne]);
s = h(:, round(dd(:)'*ros)+1);
s = s .* receive_apodization;

H_left = sparse(i(:),j(:),double(s(:)),K*Ne,scan_left.pixels);
toc

% -------------------------
%% Model Matrix constructiuon
%---------------------------

temp = [1:256*128];
temp = reshape(temp, 256,128);
temp = temp(:,end:-1:1);
pb = temp(:);


temp = [1:K*Ne];
temp = reshape(temp, K,Ne);
temp = temp(:,end:-1:1);
pa = temp(:);

H_right = H_left(:,pb);
H_right = H_right(pa,:);

H = [H_left, H_right];
Ht = H';


%% Caculate the whiten matrix C
load('../SVD/01_simulation/svd/V.mat')
load('../SVD/01_simulation/svd/lbd.mat')
C = V' ./ lbd;

% -------------------------------------------------------------
%% DDRM Input Builder (Simulation of 3-channel synthetic data) 
%--------------------------------------------------------------
% Param
close all
clearvars -except C  H  Ht scan K Ne
% Ground Truth
ImgIdx = 3;

path = ['01_simulation/SimulationResults/' num2str(ImgIdx) '/'];
try 
    load([path 'o_orig_' num2str(ImgIdx) '.mat'])
catch
    warning('o_orig does not exist. Load orig and create o_orig...');

    load([path 'orig_' num2str(ImgIdx) '.mat'])
    x_orig = squeeze(x_orig);
    
    x_orig = x_orig(:);
    x_orig = exp(x_orig*2.5);
    x_orig = x_orig-min(min(x_orig));
    x_orig = x_orig/max(max(x_orig));
    x_orig = x_orig.*randn(65536,1);
    x_orig = x_orig/max(max(x_orig));
    save([path 'o_orig_' num2str(ImgIdx) '.mat'], 'x_orig');
    %return
end
o1 = double(x_orig(:));
o2 = double(x_orig(:));
o3 = double(x_orig(:));
o = [o1, o2, o3];

%% --Create the channel data
count = 1;
for std_noisey = [0.3, 0.7, 1, 1.5, 2, 2.5, 3, 3.5] %gamma levels  

    saveNamePre = ['simulation', num2str(std_noisey)];
    
    noisey = std_noisey*randn(K*Ne,3); 
    y = [H*o1, H*o2, H*o3] + noisey;
    
    % --Beamforming
    o_Hty= [ Ht*y(:,1), Ht*y(:,2), Ht*y(:,3)];
    o_Hty = single(o_Hty);
    
    % --Calculate the noise standard derivation with Ht 
    noisex = [Ht*noisey(:,1); Ht*noisey(:,2); Ht*noisey(:,3)];
    pd = fitdist(noisex,'Normal')
    coe_std = pd.sigma / std_noisey            
    
    % --postChecking
%     figure,
%     subplot(2,3,1); Imagenet(mean(o,2),'Ground truth', scan)
%     subplot(2,3,4); Imagenet(mean(o_Hty,2),['Ht*y \gamma=',num2str(std_noisey)], scan)
%     subplot(2,3,2); Image_realScale(mean(o,2),'Ground truth', scan)
%     subplot(2,3,5); Image_realScale(mean(o_Hty,2),['Ht*y \gamma=',num2str(std_noisey)], scan)
%     subplot(2,3,3); plot(mean(o,2));title('mean(o)');
%     subplot(2,3,6); plot(mean(o_Hty,2));title(['mean(Ht*y) \gamma=',num2str(std_noisey)]);
     
    
    % --Whiten the observation
    o_CHty = C * o_Hty;
    o_CHty = o_CHty(:);
    %figure; plot(o_CHty); title('CHty')
    
    % --save yd
    o_Hty = o_Hty(:);
    save([path, 'yd/', saveNamePre, '_Hty_',  num2str(ImgIdx),'.mat'], 'o_Hty');
    save([path, 'us_DAS/', num2str(count),'_-1.mat'], 'o_Hty');
    count = count + 1;
    save([path, 'yd/', saveNamePre, '_CHty_', num2str(ImgIdx), '.mat'], 'o_CHty');
    'saved'
    %pause(1);
    
end