function save_all_figures_as_pngs(preFix,suFix)
if nargin<1
   preFix = [];
   suFix = [];
end

allFigsHandle = findobj('Type','figure');

nFigs = numel(allFigsHandle);

for iFig = 1 : nFigs

   thisFigNameRaw = allFigsHandle(iFig).Name;
   thisFigNameProc = ['Fig_',num2str(iFig,'%2d'),'_',...
                      strrep(thisFigNameRaw,' ','_')];
   thisFigNameProc = strrep(thisFigNameProc,'.','p');
   if ~isempty(preFix)
      thisFigNameProc = [preFix,thisFigNameProc];
   end
   if ~isempty(suFix)
      thisFigNameProc = [thisFigNameProc,suFix];
   end
   thisFigNameProc = [thisFigNameProc,'.png'];
   saveas(allFigsHandle(iFig),thisFigNameProc);
end
