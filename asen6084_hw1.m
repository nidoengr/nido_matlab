% asen6084_hw1.m
close all;

D_m   = 1;
f_m   = 8;
lam_m = 500E-9;
p_m   = 2E-6;

figOptions.fig_pxl_xy = figure("Name","AiryDisk");
set(gca,'FontSize',12)
figOptions.subplotFlag = true;
% figOptions.countourFig = figure("Name","Contour");
figOptions.figNoiseFlag = true;
figOptions.figNoise = figure("Name","Noise");
subplotSeq = [311 312 313];
i=1;
for nPhotons = [10 100 1000]
   figOptions.subplot = subplotSeq(i); i=i+1;
   figOptions.contourFlag = false;
   generate_expected_image_counts(D_m,f_m,lam_m,p_m,nPhotons,figOptions)
end
%%
figure(figOptions.fig_pxl_xy )
for i = [311 312 313], subplot(i); 
   xlim([-2000 9000]); ylim([-1500 1500]);
   legend("Location","east","Interpreter","latex");
   ax = gca;
   ax.FontSize = 14;
end

% set(0,'units','pixels');
% X1Y1X2Y2=get(0,'ScreenSize');

% axesHandles = findall(figOptions.fig_pxl_xy,'type','axes');
% get(figOptions.fig_pxl_xy,'position')


movegui("east")

figOptions.fig_pxl_xy.Position = [769.8000   49.8000  766.4000  732.8000];

%% save
% saveas(figOptions.fig_pxl_xy,'C:\Users\joses\Dropbox\cu_boulder\cu_2022_02\hw1p3_figure_v2.png')