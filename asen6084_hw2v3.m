% asen6084_hw2v2
close all; clc; clear all;
rng(20220902)
%% Inputs
% 5 objects with photometric SNR
SNR = [ 3 10 30 100 300 ]';
nObjs = numel(SNR);
% Pixel dimensions
mxi_pxl = 500;
mxj_pxl = 500;

m = mxi_pxl * mxj_pxl;
      
% focal distance, aperture diameter, pixel variance, pixel pitch, wavelength 
f_m = 4;
D_m = 0.5;
sig_pxl = 10;
p_m = 1.5E-6;
lam_m = 500E-9;

%% 1. Create Image

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
S = SNR .* sqrt(mTW) .* sig_pxl; % <--- tracking window mTW


% Poisson ratio, given by professor
lambda_poisson = sig_pxl^2;

% Create an image with no noise
IM_noBNoise = zeros(mxj_pxl,mxi_pxl);

% Creat an image background with poisson distribution
%  I_{b}
IM_BNoise   = poissrnd(lambda_poisson,mxi_pxl,mxj_pxl);
imshow(mat2gray(IM_BNoise));

% Image with objs and noise initialization
IM_BNoise_Objs = IM_BNoise; % start here

% Setup data retention for each object's own image
IM_Obj_i = repmat(IM_noBNoise,1,1,nObjs);

% RGB is the synthetic image with markers
RGB = IM_noBNoise;

% The random center coordinates [ x_c,  y_c ]
osi_pxl = zeros(nObjs,2);

% Loop through objects and super position objects
for iObj = 1 : nObjs 
   
   
   % Object origin randomized
   osi_pxl(iObj,:) = [randi(mxi_pxl),randi(mxj_pxl)]; 

   % This object's count
   nCounts = S(iObj);

   % Create the object origin based on the counts i.e. counts= rows
   Osi_pxl = repmat(osi_pxl(iObj,:),nCounts,1);
   
   % Initialize the the objects XY coordinates based on the counts (rows)
   % % REF: a + (b-a)*sig, normal dist between a and b
   XYs1_bar = round(repmat(-sigPsfPixel_pxl,nCounts,2) + ...
      (sigPsfPixel_pxl*2) .* randn(nCounts,2) + Osi_pxl);

   % Loop through the counts and increment 
   for iCount = 1 : nCounts
      % Check if the XYcoord is bound by 500x500 image and positive
      if all([XYs1_bar(iCount,1),XYs1_bar(iCount,2)]>0) && ...
            XYs1_bar(iCount,1)<=mxi_pxl && XYs1_bar(iCount,2)<=mxj_pxl
         % Increment the value of the xy indecies (intensity)
         IM_Obj_i(XYs1_bar(iCount,1),XYs1_bar(iCount,2),iObj) = ...
            IM_Obj_i(XYs1_bar(iCount,1),XYs1_bar(iCount,2),iObj) + 1;
      end
   end

   % Super position images together
   IM_BNoise_Objs = IM_Obj_i(:,:,iObj) + IM_BNoise_Objs;

   % Plot
   %imshow(mat2gray(IM_BNoise_Objs));
   %hold;
   %pause(2.5)
end
% hold off;

% Show the syntheic image
imshow(mat2gray(IM_BNoise_Objs));
h = gca;
h.Visible = 'on';

IM_true = mat2gray(IM_BNoise_Objs);
figure("Name","True Synthetic Image");

OsXY = [osi_pxl(:,2),osi_pxl(:,1)]; % Its a bit funky here but this lines up
RGB = insertMarker(IM_true,OsXY,"plus","Size",12,"Color","magenta");
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

%% 2. Thresholding & Clustering

% Build Image mask
I_m = nan(size(IM_BNoise_Objs));
nmult = 6;
[ixb, iyb] =find( IM_BNoise_Objs <= (nmult * sig_pxl)  );
I_m(ixb,iyb) = IM_BNoise_Objs(ixb,iyb);
Bhat = median(I_m,"all","omitnan"); % This is the same as lambda
% Bhat = median(IM_true,"all","omitnan"); % This is the same as lambda

% Estimate backgound
% IB_hat = imboxfilt(IM_true,11);
kern_func = @(A) median(A,"all");
% kern_func = @(A) norm(A,1);
% kern_func = @(A) norm(A,2);
% kern_func = @(A) norm(A,"inf");

% Here we would start to iterate over the size of the window? 
nSigMult = 60; %10; % 60
% for nSigMult = 1 : 2 : 10 % This searching which multiplier for the
% tracking window would be good, dont forget the end

mxiTW_pxl = 2*floor((nSigMult*sigPsfPixel_pxl)/2)+1;
mxjTW_pxl = 2*floor((nSigMult*sigPsfPixel_pxl)/2)+1;
mTW = mxjTW_pxl * mxiTW_pxl;

% Subtract the background

% Calculate local mean
% IB_hat = nlfilter(IM_true,[mxiTW_pxl mxjTW_pxl],kern_func);
% IB_hat = mat2gray(IM_BNoise_Objs - kern_func(IM_BNoise_Objs));% - 100);
% IB_hat = IM_BNoise_Objs - Bhat;
IB_hat = IM_BNoise_Objs - Bhat;
IB_hat =  IB_hat + min(IB_hat,[],'all');
idxBhatLtZ = find(IB_hat<0);
IB_hat(idxBhatLtZ) = zeros(size(idxBhatLtZ));

figure("Name","Background Subtracted Image"); imshow(mat2gray(IB_hat))
h = gca;
h.Visible = 'on';
% Plots
figure("Name","True Image next to Background Subtracted Image")
imshowpair(mat2gray(RGB),IB_hat,'montage')
h = gca;
h.Visible = 'on';

% IB_hat = mat2gray(IB_hat);
stdPxlS = floor(std(IB_hat(1:mxiTW_pxl,1:mxjTW_pxl),[],"all"));
% stdPxlS = floor(std(IB_hat,[],"all"));
% avgPxlS = mean(IB_hat,"all");
avgPxlS = mean(IB_hat(1:mxiTW_pxl,1:mxjTW_pxl),"all");
% medPxlS = median(IB_hat,"all");
medPxlS = median(IB_hat(1:mxiTW_pxl,1:mxjTW_pxl),"all");
% IBm = im2bw(IB_hat,1-4.89*stdPxlS);
% idxIBm = im2bw(IB_hat,avgPxlS+4.59*stdPxlS);

%--------------
% Thresholding here
%--------------
% idxIBm = im2bw(mat2gray(IB_hat),1);
% idxIBm = im2bw(IB_hat,avgPxlS);
% idxIBm = find(IB_hat>=(stdPxlS));
idxIBm = find(IB_hat>=(avgPxlS));
% idxIBm = find(IB_hat>=(medPxlS));
% idxIBm = find(IB_hat>0);
% idxIBm = find(IB_hat>= 4*sig_pxl)'


% IB_hat = IM_true;
IB_hat(idxIBm) = ones(size(IB_hat(idxIBm)));
IB_hat(~idxIBm) = zeros(size(IB_hat(~idxIBm)));
figure("Name","Image Threshold Image"); imshow(IB_hat)
h = gca;
h.Visible = 'on';

% sig_n = floor(std(double(IB_hat),0,"all"));
% idxThresh = IB_hat >= (sig_n/3);
% IB_hat(idxThresh) = uint8(0); 

% Plots
figure("Name","True Image next to Threshold Image")
imshowpair(mat2gray(RGB),IB_hat,'montage')
h = gca;
h.Visible = 'on';
xticks(0:100:1000);
lab = split(num2str([0:100:500,100:100:500]));
xticklabels(lab);
title('Truth Image vs Threshold')
% figure("Name","Truth Image diff with Threshold")
% imshowpair(mat2gray(RGB),IB_hat,'diff')
% h = gca;
% h.Visible = 'on';
% title('Truth Image diff Threshold')
% figure("Name","Truth Image diff with Threshold")
% imshowpair(mat2gray(RGB),IB_hat,'falsecolor')
% h = gca;
% h.Visible = 'on';
% title('Truth Image falsecolor nlfilter()')

% end % from iterating on nSigMult for the mTW tracking window

%--------------
% Clustering
%--------------
% IB_hat_flat = IB_hat(:);%IB_hat(:);

[xdetc_pxl,ydetc_pxl] = find(IB_hat);

det_coord_pxl = [xdetc_pxl,ydetc_pxl];
disp('Detection coordinates (first few):')
disp(det_coord_pxl(1:5,:))

% [IdxClust, C] = dbscan(det_coord_pxl,6*sigPsfPixel_pxl,80);
[IdxClust, C] = dbscan(det_coord_pxl,6*sigPsfPixel_pxl,3);
CentroidsHat = calc_centroids_from_clusters(ydetc_pxl,xdetc_pxl,IdxClust,IB_hat);
idxTrueMatchClust = [3 4 2 5];
idxTrueMatchClust = [ 4 3 5 2];
for icent = 1 : numel(CentroidsHat(:,1))
   idxT = idxTrueMatchClust(icent);
   thisAbsErrMag = sqrt(sum((osi_pxl(idxT,:) - CentroidsHat(icent,:)).^2));
   fprintf('Object %d & %7.3f & %7.3f & %7.3f & %7.3f & %7.3f \n',...
      idxT,CentroidsHat(icent,:),osi_pxl(idxT,:),thisAbsErrMag)
%    fprintf('Object %d & %7.3f & %7.3f \\\\ \n',icent,CentroidsHat(icent,:));
end
fprintf('Object %d & n/a & n/a & %7.3f & %7.3f & n/a \n',...
      1,osi_pxl(idxT,:))
%%
figure("Name","Clustering Results")
% scatter(xdetc_pxl,ydetc_pxl,'.k','DisplayName','All Detections');
hold on;
detfig = gscatter(ydetc_pxl,xdetc_pxl,IdxClust);
detfig(1).Color='k';
% detfig(2).Color=[0 1 0];
% detfig(3).Color=[0;139;69]./norm([0;139;69]);%[178;58;238]./norm([178;58;238]);
% detfig(4).Color=[154;205;50]./norm([154;205;50]);%
% detfig(5).Color=[102;205;0]./norm([102;205;0]);%

detfig_ax = gca;
detfig_ax.YDir ="reverse";
legend('Outlier','1st obj','2nd obj','3rd obj','4th obj',Location='northeast')
xlabel(''); ylabel('')

hold on;

plot(CentroidsHat(:,2),CentroidsHat(:,1),'+g','MarkerSize',12,'LineWidth',2,'DisplayName','Est. Centroid');
plot(osi_pxl(:,2),osi_pxl(:,1),'+m','MarkerSize',12,'LineWidth',2,'DisplayName','True Centroid');
axis equal
xlim([0 500])
ylim([0 500])
grid minor;
