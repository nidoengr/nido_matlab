function PixelCoordMat = generate_expected_image_counts(D_m,f_m,lam_m,p_m,nPhotons,figOptions)
%
% D_m:   Apeture diameter (m)
% f_m:   focal distance (m)
% lam_m: wavelength (m)
% p_m:   pixel pitch
rng(1);
includePixelNoise = true;
sigPixel_pxl = 5;
if nargin<6
   figOptions.figNoiseFlag = false;
   figOptions.subplotFlag  = false;
%    figOptions.fig_pxl_xy   = figure("Name","AiryDisk");
   figOptions.contourFlag  = true;
%    figOptions.countourFig  = figure("Name","Contour");

else
   figure(figOptions.fig_pxl_xy);
   if figOptions.subplotFlag

      subplot(figOptions.subplot);
   end
end


% Wave number
k= 2*pi*lam_m;

% Focal ratio
N = f_m / D_m; % focal distance to apeture diameter ratio
A = 1/(2*N);

% Airy Disk Angular edge approximation
% Saying three sigma because its approximately most of the photons in that
threeSigTheta_deg = asind(1.22.*(lam_m./(2.*A)));
sigPsfTheta_deg = threeSigTheta_deg ./ 3;

% Instantaneous field of view
IFOV = 2*atan(p_m/(2*N*D_m));

% Centroid uncertainty of the measurement
ST2 = lam_m^2 / ( 16 * pi^2 * (D_m/4)^2 * nPhotons);

% Define the circle in pixel space
sigPsfPixel_pxl = round(( 6.* sigPsfTheta_deg ) / IFOV );

% Here take the 3.*sigPsfPixel as the radius of the airy disk
rAiryDisk_pxl = 6.*sigPsfPixel_pxl;

% loop through and create a circle of rAiryDisk_pxl in 2d plot
t=linspace(0,2*pi,100)';
airyDiskCoordTheo = rAiryDisk_pxl.*[sin(t),cos(t)];
airyDiskCoord_pxl = round(airyDiskCoordTheo);
linWid = 2;
% plot(airyDiskCoordTheo(:,1),airyDiskCoordTheo(:,2),"-.","DisplayName","Theoretical $3 \sigma_{PSF,pixel}$","LineWidth",linWid); hold on;
% stairs(airyDiskCoord_pxl(:,1),airyDiskCoord_pxl(:,2),"DisplayName","Rounded $3 \sigma_{PSF,pixel}$","LineWidth",linWid);


% Determine the number of photons
PixelCoordMat = nan(nPhotons,2);
xbar_c = 0; ybar_c = 0;

xhat_c = round(xbar_c + sigPsfPixel_pxl .* randn([nPhotons,1]));%( 1 + sig1 * randn ) + sig2 * randn;
yhat_c = round(ybar_c + sigPsfPixel_pxl .* randn([nPhotons,1]));% ( 1 + sig1 * randn ) + sig2 * randn;

if includePixelNoise
   xnoise_pxl = round(randn(size(xhat_c))).*sigPixel_pxl;
   ynoise_pxl = round(randn(size(xhat_c))).*sigPixel_pxl;
   xhat_c = xhat_c + xnoise_pxl ;
   yhat_c = yhat_c + ynoise_pxl;
   PixelNoiseCoord = round(rand(size(PixelCoordMat))) .* sigPixel_pxl;
%    plot(PixelNoiseCoord(:,1),PixelNoiseCoord(:,2),"sk","DisplayName","Pixel noise location (uniform)","MarkerSize",1,"MarkerFaceColor","auto")
end

PixelCoordMat(:,:) = [xhat_c, yhat_c];

% plot(xhat_c,yhat_c,"s","DisplayName","Photon arrival $\sim\mathcal{N}(0,\sigma_{PSF,pixel}^2)$","MarkerSize",2.5,"MarkerFaceColor","auto")
% grid minor; axis equal;
% legend("location","best","Interpreter","latex")
% title(['$n_{photons} = ',num2str(nPhotons,'%d'),'$'],"Interpreter","latex");


% if figOptions.figNoiseFlag
% figure(figOptions.figNoise);
% subplot(figOptions.subplot);
% plot(xnoise_pxl); hold on;
% plot(ynoise_pxl);
% legend("Noise X", "Noise Y");
% grid minor
% xlabel('n counts')
% ylabel('Pixel noise (pxl unit)')

% end

% if figOptions.contourFlag
%    figure(figOptions.countourFig);
%    if figOptions.subplotFlag
%       subplot(figOptions.subplot);
%    end
   npro_at_points = ones([nPhotons,1])';
%    for iM = 1 
%    [X,Y] = meshgrid(xhat_c,yhat_c);
%    Z = ones(size(X));
   %    Z = sqrt( X.^2 + Y.^2 );
%       [X,Y,Z] = meshgrid(xhat_c,yhat_c,npro_at_points);
   %    % sqrt(sum(PixelCoordMat.^2,2))

%    contour(X,Y,Z);

% end
% Tracking Window

% m = 1/( 12 * sig3^2 * 2 * sqrt(pi));

end