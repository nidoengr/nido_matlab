function plot_confirmed_tracks_over_time(TRACKS)


clrtime=spring(40);
fig_confirmed_tracks_v_time = figure("Name","Confirmed tracks over time");
for i=1:40
   nCt = numel(TRACKS(i).confirmedTracks);
   for j=1:nCt
      hold all;
      xplt = TRACKS(i).confirmedTracks(j).State(1);
      yplt = TRACKS(i).confirmedTracks(j).State(3);
      plot(xplt,yplt,'^','Color',clrtime(i,:));
      if i<2 || i==40
         figgca = gca;
         figgca.YDir ="reverse";
         %grid minor;
         xlim([0 2048]);
         ylim([0 2048]);
         axis equal;
      end
   end

%    if i==1
%       keyboard();
%    else
%       pause(0.1);
%    end
   drawnow;
end

grid minor;
xlim([0 2048]);
ylim([0 2048]);
axis equal;
drawnow;
% grid minor;
% xlim([0 2048]);
% ylim([0 2048]);
fig_confirmed_tracks_v_time.Colormap = clrtime;
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
figgca.FontSize = 20;