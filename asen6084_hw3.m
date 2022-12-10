% asen6084_hw3
close all; clc; clear all;
rng(20220929);
fontSizeVal = 14;
%% Inputs
% 5 objects with photometric SNR
SNR = [ 300 30 ]';
nObjs = numel(SNR);
% Pixel dimensions
mxi_pxl = 500;
mxj_pxl = 500;

m = mxi_pxl * mxj_pxl;
      
% focal distance, aperture diameter, pixel variance, pixel pitch, wavelength 
f_m = 4;
D_m = 0.5;
sig_readnoise = 10;
p_m = 1.5E-6;
lam_m = 500E-9;

%% 1. Create Sequence of Images

% $\sigma_{\theta} - assume 6 sigma theta here assumed here
sixSigTheta_rd = asin(1.22*lam_m*f_m/D_m);

% IFOV
IFOV_rd = 2*atan(p_m/(2*f_m));

% Pixel space theta
sigPsfPixel_pxl = sixSigTheta_rd / (6*IFOV_rd);

% Tracking window params
nSigMult = 10; %<---------ADHOC, prof recommended 6 (good), 10 seems good 12, 15?
               % Note: this will change the repeatable "random" locations
               % because we feed the random number generation different
               % counts
mxiTW_pxl = 2*floor((nSigMult*sigPsfPixel_pxl)/2)+1;
mxjTW_pxl = 2*floor((nSigMult*sigPsfPixel_pxl)/2)+1;
mTW = mxjTW_pxl * mxiTW_pxl;

% Get Number of counts for each object
S = SNR .* sqrt(mTW) .* sig_readnoise; % <--- tracking window mTW


% Poisson ratio, given by professor
lambda_poisson = sig_readnoise^2;

% Create an image with no noise
IM_noBNoise = zeros(mxj_pxl,mxi_pxl);

% Creat an image background with poisson distribution
%  I_{b}
IM_BNoise   = poissrnd(lambda_poisson,mxi_pxl,mxj_pxl);
figure("Name","Image with Background Only");
imshow(mat2gray(IM_BNoise));
set(gca,'FontSize',fontSizeVal) % Creates an axes and sets its FontSize 

% Image with objs and noise initialization
IM_BNoise_Objs = IM_noBNoise;% ---> Do this later: IM_BNoise_Objs = IM_BNoise_Objs + IM_BNoise; % start here

% Setup data retention for each object's own image
IM_Obj_i = repmat(IM_noBNoise,1,1,nObjs);

% RGB is the synthetic image with markers
RGB = IM_noBNoise;

% The random center coordinates [ x_c,  y_c ]
osi_pxl = zeros(nObjs,2);

nFrames = 50;

xStreak = linspace(2,482,nFrames);%linspace(250,250,nFrames);
% yStreak = 4/5.*xStreak;%linspace(2,482,nFrames);
yStreak = 1/3.*xStreak+1;%linspace(2,482,nFrames);

osi_pxl= [xStreak',yStreak'];

ImSeq = repmat(struct('IDATA',zeros([mxi_pxl, mxj_pxl, nFrames ]),...
                      'IDATAwithMarkers',zeros([mxi_pxl, mxj_pxl, nFrames ]),...
                      'BDATA',zeros([mxi_pxl, mxj_pxl, nFrames ])...
                      ),nObjs,1);

% Loop through objects and super position objects
for iObj = 1 : nObjs 
% iObj = 1; % <-- This was for debugging
showNFrames = [10 20 30 40 50 ];
iShow = 1;
for iFrame = 1 : nFrames
   
   
   [tempIM_true, tempRGB, IM_BNoise] = generate_image(SNR(iObj), mxi_pxl, mxj_pxl, f_m, D_m, sig_readnoise, p_m,...
   lam_m, osi_pxl(iFrame,:));
   ImSeq(iObj).IDATA(:,:,iFrame) = tempIM_true;
%    ImSeq(iObj).IDATAwithMarkers(:,:,iFrame) = tempRGB;
   ImSeq(iObj).BDATA(:,:,iFrame) = IM_BNoise;
   IM_BNoise_Objs = tempIM_true;

   % Plot
   if iFrame == showNFrames(iShow)
%    figure("Name",['Image number = ',num2str(iFrame,'%d')]);
% %    imshow(mat2gray(IM_BNoise_Objs));
%    imshowpair(mat2gray(tempIM_true),mat2gray(tempRGB),'montage');
   % h = gca;
   % h.Visible = 'on';
   % It is expected that after each iteration, the image in 
   % the figure drowns out the previous iteration's object
   % because mat2gray normalizes the matrix.
%    hold;
   iShow= iShow + 1;
   end
%    pause(2.5)
end

% Show the syntheic image
% figure("Name","Image with Objects with Background Noise Added");
% IM_BNoise_Objs = IM_BNoise_Objs + IM_BNoise;
% imshow(mat2gray(IM_BNoise_Objs));
% h = gca;
% h.Visible = 'on';

% IM_true = mat2gray(IM_BNoise_Objs); % Older going to use it differently
% here
IM_true = ImSeq(1).IDATA(:,:,end);%IM_BNoise_Objs;
figure("Name","True Synthetic Image with Markers");

OsXY = [osi_pxl(:,2),osi_pxl(:,1)]; % Its a bit funky here but this lines up
RGB = insertMarker(mat2gray(IM_true),OsXY,"plus","Size",12,"Color","magenta");
% hax =  axes('BoxStyle','full');
% ax1 = axes('Position',[0.1 0.1 .6 .6],'Box','on');
% ax2 = axes('Position',[.35 .35 .6 .6],'Box','on');
imshow(RGB)
% imshowpair(IM_true,mat2gray(RGB),'montage');
h = gca;
h.Visible = 'on';
% xticks([0:100:1000]);
% lab = split(num2str([0:100:500,100:100:500]))
% xticklabels(lab);
set(gca,'FontSize',fontSizeVal) % Creates an axes and sets its FontSize 
end


%% 2. Thresholding & Clustering

Bhat = sig_readnoise^2;

Ib = IM_true - Bhat;
figure("Name","Background subtracted Image I - sigma_readnoise");
imshow(mat2gray(Ib));
IB_hat = Ib;
figure("Name","Background Subtracted Image"); imshow(mat2gray(IB_hat))
h = gca;
h.Visible = 'on';
% Plots
figure("Name","True Image next to Background Subtracted Image")
imshowpair(mat2gray(RGB),mat2gray(IB_hat),'montage')
h = gca;
h.Visible = 'on';


%% Match Filter

% Pixel space theta
% sigPsfPixel_pxl

% Tracking window params
nSigMult = 3; %<---------ADHOC, prof recommended 6 (good), 10 seems good 12, 15?
               % Note: this will change the repeatable "random" locations
               % because we feed the random number generation different
               % counts
% mxiTW_pxl = 2*floor((nSigMult*sigPsfPixel_pxl)/2)+1;
% mxjTW_pxl = 2*floor((nSigMult*sigPsfPixel_pxl)/2)+1;
% this is kind of what anthony did
mxiTW_pxl = floor((nSigMult*sigPsfPixel_pxl))+1;
mxjTW_pxl = floor((nSigMult*sigPsfPixel_pxl))+1;
mTW = mxjTW_pxl * mxiTW_pxl;
xrng = 1:mxiTW_pxl;
yrng = 1:mxjTW_pxl;
midpnt = median(xrng);
centMTWrng = xrng - midpnt; % square window

% Convolution function
g = exp(-((centMTWrng/2).^2+(centMTWrng'/2).^2)./(2*sigPsfPixel_pxl^2)); 
% 


% MF Score
% convTic = tic;
MFscore = conv2(Ib, g, 'same');
% convToc = toc(convTic);
% fprintf('Elapsed time for convolution (s): %f\n',convToc)
% 
% 
% % MF Variance
sigmaMF = sqrt(var(MFscore(:)));
% 
% % MF Threshold
nSigMFMult = 4;
T_MF = sigmaMF * nSigMFMult;
% 
% % MF Detections
idxDetMF= MFscore >= T_MF;
DetMFMat = zeros(size(IB_hat));
DetMFMat(idxDetMF) = ones(size(IB_hat(idxDetMF)));

figure("Name","Match Filter Detections Image")
imshow(mat2gray(DetMFMat))
set(gca,'FontSize',fontSizeVal) % Creates an axes and sets its FontSize to 18



[xMFdetc_pxl,yMFdetc_pxl] = find(DetMFMat); %<--- Works


detMF_coord_pxl = [xMFdetc_pxl,yMFdetc_pxl];
disp('Detection coordinates (first few):')
disp(detMF_coord_pxl (1:5,:))

ticDbscan = tic;
[IdxClustMF, C] = dbscan(detMF_coord_pxl ,6*sigPsfPixel_pxl,30);
tocDbscan = toc(ticDbscan);
fprintf('Elapsed time for dbscan (s): %f\n',tocDbscan);
CentroidsMFHat = calc_centroids_from_clusters(yMFdetc_pxl,xMFdetc_pxl,IdxClustMF,DetMFMat);
% 
%% MF Clustering Results Figure
figure("Name","MF Clustering Results")
hold on;
detfig = gscatter(yMFdetc_pxl,xMFdetc_pxl,IdxClustMF);
% detfig(1).Color='k';
% % detfig(2).Color=[0 1 0];
% % detfig(3).Color=[0;139;69]./norm([0;139;69]);%[178;58;238]./norm([178;58;238]);
% % detfig(4).Color=[154;205;50]./norm([154;205;50]);%
% % detfig(5).Color=[102;205;0]./norm([102;205;0]);%
% 
detfig_ax = gca;
detfig_ax.YDir ="reverse";
legend('Object Detection Group',Location='northeast')
xlabel(''); ylabel('')
detfig_ax.FontSize = 14;

hold on;

plot(CentroidsMFHat(:,2),CentroidsMFHat(:,1),'+g','MarkerSize',12,'LineWidth',2,'DisplayName','Est. Centroid');
plot(osi_pxl(:,2),osi_pxl(:,1),'+m','MarkerSize',12,'LineWidth',2,'DisplayName','True Centroid');
axis equal
xlim([0 500])
ylim([0 500])
grid minor;


%% BINARY HYPOTHIS WITH MATCH FILTER
minPntsMat = [[3 3 3]
                [3 10 100]];
nsigmultiterArray=[6 10];
iMinPntsArray = 1;
for iObj = 1:nObjs
   nSigMultIter = nsigmultiterArray(iObj);
   minPntsArray = minPntsMat(iObj,:);
   thisIb = squeeze(ImSeq(iObj).IDATA(:,:,end)-Bhat);
for p_FA = [1e-4 0.01 0.10]

% F0 = cdf('Normal',x,mu,sigPsfPixel_pxl)
% f0 = pdf('Normal',x,mu,sigma)
gamma = icdf('Normal',1-p_FA,0,sigmaMF);

% gamma is now our MF threshold that helps us determine 
% if a given pixel is a detection or not because if our 
% MF score at a pixel is higher than gamma, we can reject H0, right?
% Yes that gives us gamma, or in the MF application that gives us T_MF
% And then you threshold as before, creating a logical matrix based on MF_score >T_MF

% Convolution function
g = exp(-((centMTWrng/2).^2+(centMTWrng'/2).^2)./(2*sigPsfPixel_pxl^2)); 

% MF Score
% convTic = tic;
MFscore = conv2(thisIb, g, 'same');
% convToc = toc(convTic);
% fprintf('Elapsed time for convolution (s): %f\n',convToc)
 
% % MF Variance
sigmaMF = sqrt(var(MFscore(:)));
% 
% % MF Threshold
nSigMFMult = 1;
% T_MF = sigmaMF * nSigMFMult;
T_MF = nSigMFMult *gamma;
% 
% % MF Detections
idxDetMF= MFscore >= T_MF;
DetMFMat = zeros(size(thisIb));
DetMFMat(idxDetMF) = ones(size(thisIb(idxDetMF)));

figure("Name","Match Filter Detections Image")
imshow(mat2gray(DetMFMat))
set(gca,'FontSize',fontSizeVal) % Creates an axes and sets its FontSize to 18



[xMFdetc_pxl,yMFdetc_pxl] = find(DetMFMat); %<--- Works


detMF_coord_pxl = [xMFdetc_pxl,yMFdetc_pxl];
disp('Detection coordinates (first few):')
disp(detMF_coord_pxl (1:5,:))
%%
try
ticDbscan = tic;
[IdxClustMF, C] = dbscan(detMF_coord_pxl ,nSigMultIter*sigPsfPixel_pxl,minPntsArray(iMinPntsArray));
tocDbscan = toc(ticDbscan);
fprintf('Elapsed time for dbscan (s): %f\n',tocDbscan);
CentroidsMFHat = calc_centroids_from_clusters(yMFdetc_pxl,xMFdetc_pxl,IdxClustMF,DetMFMat);

%% MF+BINARY TESTING Clustering Results Figure
figure("Name",['MF+BINARY Clustering Results p_FA = ',num2str(p_FA,'%0.2e')])
hold on;
detfig = gscatter(yMFdetc_pxl,xMFdetc_pxl,IdxClustMF);
% detfig(1).Color='k';
% % detfig(2).Color=[0 1 0];
% % detfig(3).Color=[0;139;69]./norm([0;139;69]);%[178;58;238]./norm([178;58;238]);
% % detfig(4).Color=[154;205;50]./norm([154;205;50]);%
% % detfig(5).Color=[102;205;0]./norm([102;205;0]);%
% 
detfig_ax = gca;
detfig_ax.YDir ="reverse";
legend('Object Detection Group',Location='northeast')
xlabel(''); ylabel('')
detfig_ax.FontSize = 14;

hold on;

plot(CentroidsMFHat(:,2),CentroidsMFHat(:,1),'+g','MarkerSize',12,'LineWidth',2,'DisplayName','Est. Centroid');
plot(osi_pxl(:,2),osi_pxl(:,1),'+m','MarkerSize',12,'LineWidth',2,'DisplayName','True Centroid');
axis equal
xlim([0 500])
ylim([0 500])
grid minor;
catch
   fprintf('dbscan did not work for p_FA = %0.3e\n',p_FA);
end
   iMinPntsArray = iMinPntsArray + 1;
end
end