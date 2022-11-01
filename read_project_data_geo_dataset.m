function READDATA = read_project_data_geo_dataset()
%read_project_data_geo_dataset reads the
%   GEO_dataset.zip/MultiTarget fits files
%   and places them in a struct array
%
% Inputs: none;
% Outputs: READATA  an nFitsFiles x 1 column struct array
%                   with field FITSREADDATA
%
%          READATA.FITSREADDATA is the output of the
%                               fitsread function
%

tic();
fitsDir = 'C:\Users\joses\Documents\MultiTarget';

listDir = dir(fitsDir);

% Loop through directory list, note idx 1 and 2 are . and ..
countFits =0;
for i=3:numel(listDir)
   if contains(listDir(i).name(end-3:end),'.fit')
      countFits=countFits+1;
      %fprintf('%s\n',listDir(i).name);
   end
end

nFits = countFits;

fitsFullPathsCell = cell([nFits,1]);

idxCell=1;
for i=3:numel(listDir)
   if contains(listDir(i).name(end-3:end),'.fit')
      fitsFullPathsCell{idxCell} =  fullfile(listDir(i).folder,listDir(i).name);
      idxCell=idxCell+1;
      %fprintf('%s\n',listDir(i).name);
   end
end
% disp(fitsFullPathsCell{1:3,:})

READDATA = repmat(struct('FITSREADDATA',nan,...
                         'FITSINFO',nan),nFits,1);
for i=1:nFits
   READDATA(i).FITSREADDATA = fitsread(fitsFullPathsCell{i,:});
%    READDATA(i).FITSINFO = fitsinfo(fitsFullPathsCell{i,:});
end
fprintf('Dataset read time: %f\n',toc)
% READDATA(1).FITSREADDATA
% 
% imshow(READDATA(1).FITSREADDATA)
figure("Name","First Image")
imshow(mat2gray(READDATA(1).FITSREADDATA))