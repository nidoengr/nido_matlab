%% This part works - need to finish updloading/downloading truth data (fits images to astrometry.net)
astrometrycalspath = 'C:\Users\joses\Documents\AstrometryDotNet';
nImgs=40;
IMTRUE = repelem(struct('astrometryCorrTable',nan),nImgs,1);
for i=1:nImgs
filename=sprintf('corr_%02d.fits',i);
fullpath=fullfile(astrometrycalspath,filename);
bindata=fitsread(fullpath,'binarytable');
binarray=cell2mat(bindata);
fileinfo=fitsinfo(fullpath);
nfields = fileinfo.BinaryTable.NFields;
tit = cell(nfields,1);

charary = '';
itit = 1;
for j=1:height(fileinfo.BinaryTable.Keywords)
   comparisonstr = sprintf('TTYPE%d',j);
   idxHead = strcmp(fileinfo.BinaryTable.Keywords(:,1),comparisonstr);
   if any(idxHead)
      disp([comparisonstr,' -> ',fileinfo.BinaryTable.Keywords{idxHead,1}]);
      thisstr = sprintf('bindata{%d},',itit);
      tit{itit} = fileinfo.BinaryTable.Keywords{idxHead,2};
      charary = strcat(charary,thisstr);
      itit = itit+1;
   end
end
IMTRUE(i).astrometryCorrTable = eval(['table(',charary(1:end-1),')'])
IMTRUE(i).astrometryCorrTable.Properties.VariableNames = tit;
end


%% PLOT
hold all; 
for i=1:nImgs
%    plot(IMTRUE(i).astrometryCorrTable.field_y,IMTRUE(i).astrometryCorrTable.field_x,'ko','MarkerSize',10);
%    plot(IMTRUE(i).astrometryCorrTable.field_y,IMTRUE(i).astrometryCorrTable.field_x,'r.','MarkerSize',10);
   plot(IMTRUE(i).astrometryCorrTable.index_y,IMTRUE(i).astrometryCorrTable.index_x,'g.','MarkerSize',5);
end
thisgca = gca;
% set(thisgca,'Color',[0;43;43]./255)%
thisgca.Color = [0;43;43]./255;
% thisgca.ColorOrder = ...
%    [0.6350    0.0780    0.1840
%     0.8500    0.3250    0.0980
%     0.9290    0.6940    0.1250
%     0.4940    0.1840    0.5560
%     0.0000    0.4470    0.7410
%     0.4660    0.6740    0.1880
%     0.3010    0.7450    0.9330];
grid minor
% % thisgca.GridColor = [173;255;47]./255;%[0;43;43]./255;'w'; 


% thisgca.GridLineStyle
%%
figure();
clrtime=spring(40);
for i=1:nImgs
   x=10.*cos(IMTRUE(i).astrometryCorrTable.index_ra).*...
      cos(IMTRUE(i).astrometryCorrTable.index_dec);
   y=10.*sin(IMTRUE(i).astrometryCorrTable.index_ra).*...
      cos(IMTRUE(i).astrometryCorrTable.index_dec);
   z=10.*sin(IMTRUE(i).astrometryCorrTable.index_dec);
   plot3(x,y,z,'d','Color',clrtime(i,:),'MarkerFaceColor',clrtime(i,:));
   
   hold all;
end
[X,Y,Z] = sphere();
surf(X,Y,Z);
axis equal
%% OLD NOT USED
% tit = cell(13,1);
% % nchars = nfields * numel('bindata{00},') - 1;
% % charary = strcat(1:nchars)
% charary = '';
% itit = 1;
% for i=1:height(fileinfo.BinaryTable.Keywords)
%    comparisonstr = sprintf('TTYPE%d',i);
%    idxHead = strcmp(fileinfo.BinaryTable.Keywords(:,1),comparisonstr);
%    if any(idxHead)
%       disp([comparisonstr,' -> ',fileinfo.BinaryTable.Keywords{idxHead,1}]);
%       thisstr = sprintf('bindata{%d},',itit);
%       tit{itit} = fileinfo.BinaryTable.Keywords{idxHead,2};
%       charary = strcat(charary,thisstr);
%       itit = itit+1;
%    end
% end
% 
% % IMTRUE.table = eval(['table(',charary(1:end-1),')'])
% % IMTRUE.table.Properties.VariableNames = tit;
% 
% IMTRUE.table = table(bindata{1 },...
%                      bindata{2 },...
%                      bindata{3 },...
%                      bindata{4 },...
%                      bindata{5 },...
%                      bindata{6 },...
%                      bindata{7 },...
%                      bindata{8 },...
%                      bindata{9 },...
%                      bindata{10},...
%                      bindata{11},...
%                      bindata{12},...
%                      bindata{13}...
%                      )
% 
% IMTRUE.table.Properties.VariableNames = tit;