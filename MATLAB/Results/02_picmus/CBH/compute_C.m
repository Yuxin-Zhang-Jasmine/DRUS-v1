
clearvars 
close all
clc
addpath(genpath('../PICMUS/code/src'));

% ---------------
%% Params
%----------------
RxApodName =  'tukey25';  
rx_f_number = 1.4;       
[dataset,~] = datasetReader(1,1,2);      

Ne = size(dataset.probe_geometry,1);
fs = dataset.sampling_frequency;
c = dataset.c0;

% -------------------------
%% Waveform PICMUS
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
%% Geometry and TimeDelay
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

k = floor(min(delay(:)));
K = ceil(max(delay(:))) + Kh ;

% ---------------------
%% Forward Model (left)
%----------------------
tic
di = fix(delay);
dd = delay - di;

r = di+(0:Ne-1)*K;
i = (r(:)+(0:Kh-1))';
j = repmat((1:scan_left.pixels),[Kh,Ne]);
s = h(:, round(dd(:)'*ros)+1);

H_left = sparse(i(:),j(:),double(s(:)),K*Ne,scan_left.pixels);
toc

% ---------------------
%% Backward Model (left)
%----------------------
tic
rx_aperture = scan_left.z/rx_f_number;
rx_aperture_distance = abs(scan_left.x*ones(1,dataset.channels)-ones(scan_left.pixels,1)*dataset.probe_geometry(:,1).');
receive_apodization = tools.apodization(rx_aperture_distance,rx_aperture*ones(1,dataset.channels),RxApodName);
receive_apodization = repmat(receive_apodization(:).', [Kh,1]);
bs = conj(s) .* receive_apodization;

HT_apod_left = sparse(j(:),i(:),double(bs(:)),scan_left.pixels,K*Ne);
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

HT_apod_right = HT_apod_left(pb,:);
HT_apod_right = HT_apod_right(:,pa);

H = [H_left, H_right];
B = [HT_apod_left; HT_apod_right];

% ----------------------
%% Computing C
%-----------------------
tic;   BBt = B * B.';  toc

clearvars -except BBt
BBt = single(full(BBt));


% 92 mins eig(BBt)
tic;  [V,LBD] = eig(BBt, 'vector');  toc
clearvars BBt
pause(60); 


% Sorting LBD and V in descending order
[LBD,ind] = sort(LBD,'descend'); save LBD.mat LBD -mat;
V = V(:,ind); save V.mat V -mat;

% take the square root of LBD and compute C
lbd=sqrt(LBD); save lbd.mat lbd -mat;
C = V.' ./ lbd;
