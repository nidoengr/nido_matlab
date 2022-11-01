function IMDATA = project_data_match_filtering(IMDATA,nPntsReqInClust,plotFlag)
if nargin<2
nPntsReqInClust = 500;
end
nIMs = numel(IMDATA);

% focal distance, aperture diameter, pixel variance, pixel pitch, wavelength 
% warning('FIX ME! Assuming camera parameters. Change to inputs for the function.')

% From the website brochure of the optical system used in
% project
f_m =  4; % 0.600;%
D_m = 0.5; % 0.200;%
% sig_readnoise = 10;
p_m = 1.5E-6;
lam_m = 500E-9;

% $\sigma_{\theta} - assume 6 sigma theta here assumed here
sixSigTheta_rd = asin(1.22*lam_m*f_m/D_m);

% IFOV
IFOV_rd = 2*atan(p_m/(2*f_m));

% Pixel space theta
sigPsfPixel_pxl = sixSigTheta_rd / (6*IFOV_rd);

%% Match Filter

% Tracking window params
% warning('FIX ME! Hardcoded tracking window multiplier for convolution.')
nSigMult = 10; %<---------ADHOC, prof recommended 6 (good), 10 seems good 12, 15?
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
% This seems to be the integrand of one of the first few equations in Lecture 8...

for iIM = 1:nIMs

   % MF Score
   % convTic = tic;
   MFscore = conv2(IMDATA(iIM).Ib, g, 'same');
   IMDATA(iIM).MFscore = MFscore;
   % convToc = toc(convTic);
   % fprintf('Elapsed time for convolution (s): %f\n',convToc)
   %
   %
   % % MF Variance
   sigmaMF = sqrt(var(MFscore(:)));
   IMDATA(iIM).MFscore = sigmaMF;
   %
   % % MF Threshold
   nSigMFMult = 8;
   T_MF = sigmaMF * nSigMFMult;
   %
   % % MF Detections
   idxDetMF = MFscore >= T_MF;
   DetMFMat = zeros(size(IMDATA(iIM).Ib));

   DetMFMat(idxDetMF) = ones(size(IMDATA(iIM).Ib(idxDetMF)));
   IMDATA(iIM).MFDetMat = DetMFMat;

%    if iIM<2
%       figure("Name","Match Filter Detections Image (1)")
%       imshow(mat2gray(DetMFMat))
%    end


   [xMFdetc_pxl,yMFdetc_pxl] = find(DetMFMat); %<--- Works


   detMF_coord_pxl = [xMFdetc_pxl,yMFdetc_pxl];
   IMDATA(iIM).detMF_coord_pxl  = detMF_coord_pxl;
%    disp('Detection coordinates (first few):')
%    disp(detMF_coord_pxl(1:5,:))

   ticDbscan = tic;
   try
      [IdxClustMF, C] = dbscan(detMF_coord_pxl ,6*sigPsfPixel_pxl,nPntsReqInClust);
      tocDbscan = toc(ticDbscan);
%       fprintf('Elapsed time for dbscan (s): %f\n',tocDbscan);
      CentroidsMFHat = calc_centroids_from_clusters(yMFdetc_pxl,xMFdetc_pxl,IdxClustMF,DetMFMat);
%       PhotoMetCetHat = cal_photomet_cent_from_image_and_clust(IMDATA(iIM).Ib,DetMFMat,CentroidsMFHat,...
%          6*sigPsfPixel_pxl); % <-----MAKE INTO A FUNCTION
%          TO CALCULATE THE PHOTOMETRIC CENTROID!
      IMDATA(iIM).IdxClustMF  = IdxClustMF;
      IMDATA(iIM).CentroidsMFHat = CentroidsMFHat;
   catch
      warning(['Image ',num2str(iIM,'%d'),' dbscan not successful']);
      IMDATA(iIM).IdxClustMF  = nan;
      IMDATA(iIM).CentroidsMFHat = nan;
   end
   %
   %% MF Clustering Results Figure
   if plotFlag
   if iIM<5 && ~all(isnan(IMDATA(iIM).IdxClustMF))
      tempFig = figure("Name","MF Clustering Results");
      hold on;
      detfig = gscatter(yMFdetc_pxl,xMFdetc_pxl,IdxClustMF);
      detfig(1).Color='k';
      % % detfig(2).Color=[0 1 0];
      % % detfig(3).Color=[0;139;69]./norm([0;139;69]);%[178;58;238]./norm([178;58;238]);
      % % detfig(4).Color=[154;205;50]./norm([154;205;50]);%
      % % detfig(5).Color=[102;205;0]./norm([102;205;0]);%
      %
      detfig_ax = gca;
      detfig_ax.YDir ="reverse";
      % legend('Outlier','1st obj','2nd obj','3rd obj','4th obj','5th obj',Location='northeast')
      legend('off');
      % xlabel(''); ylabel('')

      hold on;

      plot(CentroidsMFHat(:,2),CentroidsMFHat(:,1),'+g','MarkerSize',12,'LineWidth',2,'DisplayName','Est. Centroid');
      % plot(osi_pxl(:,2),osi_pxl(:,1),'+m','MarkerSize',12,'LineWidth',2,'DisplayName','True Centroid');
      axis equal
      % xlim([0 500])
      % ylim([0 500])
      grid minor;
      drawnow;
      pause(1)
      saveas(tempFig,['C:\Users\joses\OneDrive\Pictures\test_run_asen_6084_proj\im_mf_proc_',num2str(iIM,'%2d')],'png');
   end
   end
end