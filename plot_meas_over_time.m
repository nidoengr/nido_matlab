function plot_meas_over_time(IMDATA)

%%
fig_meas_vs_time = figure("Name","Measurements over time");
% clrtime=cool(40);
clrtime=spring(40);

for i=1:40
%    plot(IMDATA(i).detMF_coord_pxl(:,1),...
%       IMDATA(i).detMF_coord_pxl(:,2),...
%       '+','Color',clrtime(i,:),'MarkerSize',3);
   plot(IMDATA(i).CentroidsMFHat(:,1),...
      IMDATA(i).CentroidsMFHat(:,2),...
      '+','Color',clrtime(i,:),...
      'MarkerSize',6);
   if i<2 || i==40
         hold all;
         figgca = gca;
         figgca.YDir ="reverse";
         grid minor;
         %figgca
         axis equal;
         xlim([0 2048]);
         ylim([0 2048]);
      end

   drawnow();
end
grid minor;
xlim([0 2048]);
ylim([0 2048]);
fig_meas_vs_time.Colormap = clrtime;
clrbr = colorbar;
clrbr.LimitsMode ="manual";
clrbr.TicksMode = "manual";
clrbr.TickLabelsMode = "manual";
clrbr.Limits= [0 1];
clrbr.Ticks = (0:4:40)./40;
clrbrticks = 0:4:40;
ticklabels = cell(numel(clrbr.Ticks),1);
for i=1:numel(clrbr.Ticks)
ticklabels{i} = num2str(clrbrticks(i),'%d');
end
clrbr.TickLabels = ticklabels;

% fig_meas_vs_time.WindowState ="fullscreen";
figgca.FontSize = 20;






