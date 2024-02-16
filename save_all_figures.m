function save_all_figures(preFix,suFix,ext)
arguments
   preFix (1,:) char = '';
   suFix  (1,:) char = '';
   ext    (1,:) char = '.fig'; % .png, .jpeg
end

allFigsHandle = findobj('Type','figure');

nFigs = numel(allFigsHandle);

for iFig = 1 : nFigs

   thisFigNameRaw = allFigsHandle(iFig).Name;
   % TODO better handling of the name
   thisFigNameProc = ['Fig_',num2str(iFig,'%2d'),'_',...
                      strrep(thisFigNameRaw,' ','_')];
   thisFigNameProc = strrep(thisFigNameProc,'.','p');

   if ~isempty(preFix)
      thisFigNameProc = sprintf('%s%s',preFix,thisFigNameProc);
   end
   if ~isempty(suFix)
      thisFigNameProc = sprintf('%s%s',thisFigNameProc,suFix);
   end
   thisFigNameProc = sprintf('%s%s',thisFigNameProc,ext];
   % TODO Better handling of the name
   thisFigNameProc = strrep(thisFigNameProc,':','_');
   thisFigNameProc = strrep(thisFigNameProc,' ','_');
   saveas(allFigsHandle(iFig),thisFigNameProc);
end
