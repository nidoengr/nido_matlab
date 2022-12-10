function IMDATA = image_proc_and_det(IMDATA)
ticFunc = tic();
%% Remove Backgrounds
sig_readnoise = 10;%sqrt(mean(IMDATA(1).FITSREADDATA,'all'));%
IMDATA = project_data_background_removal(IMDATA, sig_readnoise);


%% MF Thrsholding
% TODO! MAKE these parameters to pass into the functions, rather than hardcoded within the function
nSigMult = 3;
nSigMFMult = 1/2;

% A loop to find how the number of points required in the cluster impact the results
for nPntsInCluster = 3%10%100%450%[10,50:50:600] % see ..\..\cu_2022_02\project\nClusters_phase_space_mapping_on_dbscan_GEO_dataset.txt
                                    % The idea here being
                                    % start tracking smaller
                                    % number of clusters to
                                    % prove the algorithms,
                                    % then increase the
                                    % clusters to track more
                                    % objects...
    IMDATA = project_data_match_filtering(IMDATA,nPntsInCluster,false,sig_readnoise);
    
    %Stats?
    disp('-----------------------------------------------')
    fprintf('run with nReqPntsInCluster = %d\n',nPntsInCluster)
    nImgs = numel(IMDATA);
    elmsArray = nan([nImgs,1]);
    for i=1:nImgs
        % fprintf('nClusters = %d\n',size(IMDATA(i).CentroidsMFHat,1));
        elmsArray(i) = size(IMDATA(i).CentroidsMFHat,1);
    end
    
    if all(isnan(elmsArray)), fprintf('breaking loop at %d pnts in cluster',nPntsInCluster);break;end

    fprintf('average clusters = %5.2f\n', mean(  elmsArray,'all','omitnan'));
    fprintf('median clusters  = %5.2f\n', median(elmsArray,'all','omitnan'));
    fprintf('max clusters     = %5.2f\n', max(   elmsArray));
    fprintf('min clusters     = %5.2f\n', min(   elmsArray));
end
fprintf('Dataset image processing and detections elapsed time: %f\n',toc(ticFunc));