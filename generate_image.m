function  [IM_true, RGB, IM_BNoise] = generate_image(SNR, mxi_pxl, mxj_pxl, f_m, D_m, sig_readnoise, p_m,...
   lam_m, obj_loc)

if nargin<1
   % Pixel dimensions
   mxi_pxl = 500;
   mxj_pxl = 500;
   % focal distance, aperture diameter, pixel variance, pixel pitch, wavelength 
   SNR = 300;
   f_m = 4;
   D_m = 0.5;
   sig_readnoise = 10;
   p_m = 1.5E-6;
   lam_m = 500E-9;
end
randomize_obj_loc = false;
if isempty(obj_loc) || nargin<9
   randomize_obj_loc = true;
end

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
nObjs = numel(SNR);

% Poisson ratio, given by professor
lambda_poisson = sig_readnoise^2;

% Creat an image background with poisson distribution
%  I_{b}
IM_BNoise   = poissrnd(lambda_poisson,mxi_pxl,mxj_pxl);
% figure("Name","Image with Background Only");
% imshow(mat2gray(IM_BNoise));

% Create an image with no noise
IM_noBNoise = zeros(mxj_pxl,mxi_pxl);

% Image with objs and noise initialization
IM_BNoise_Objs = IM_noBNoise;% ---> Do this later: IM_BNoise_Objs = IM_BNoise_Objs + IM_BNoise; % start here

% Setup data retention for each object's own image
IM_Obj_i = repmat(IM_noBNoise,1,1,nObjs);

% The random center coordinates [ x_c,  y_c ]
osi_pxl = zeros(nObjs,2);

% figure("Name","Image with Objects no Background noise");

% Loop through objects and super position objects
for iObj = 1 : nObjs 
   % This object's count
   nCounts = S(iObj);

   if randomize_obj_loc
      % Object origin randomized
      osi_pxl(iObj,:) = [randi(mxi_pxl),randi(mxj_pxl)];      
   else
      osi_pxl(iObj,:) = obj_loc(iObj,:);
   end


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
%    imshow(mat2gray(IM_BNoise_Objs));
   % It is expected that after each iteration, the image in 
   % the figure drowns out the previous iteration's object
   % because mat2gray normalizes the matrix.
%    hold;
%    pause(2.5)
end
% hold off;
% Add background noise
IM_true = IM_BNoise_Objs +IM_BNoise;

% Show the syntheic image
% figure("Name","Image with Objects with Background Noise Added");
% IM_BNoise_Objs = IM_BNoise_Objs + IM_BNoise;
% imshow(mat2gray(IM_BNoise_Objs));
% h = gca;
% h.Visible = 'on';

% IM_true = mat2gray(IM_BNoise_Objs); % Older going to use it differently
% here
% IM_true = IM_BNoise_Objs;
% figure("Name","True Synthetic Image with Markers");
% 
OsXY = [osi_pxl(:,2),osi_pxl(:,1)]; % Its a bit funky here but this lines up
RGB = insertMarker(mat2gray(IM_true),OsXY,"plus","Size",12,"Color","magenta");
% hax =  axes('BoxStyle','full');
% ax1 = axes('Position',[0.1 0.1 .6 .6],'Box','on');
% ax2 = axes('Position',[.35 .35 .6 .6],'Box','on');
% imshow(RGB)
% imshowpair(IM_true,mat2gray(RGB),'montage');
% h = gca;
% h.Visible = 'on';
% xticks([0:100:1000]);
% lab = split(num2str([0:100:500,100:100:500]))
% xticklabels(lab);

