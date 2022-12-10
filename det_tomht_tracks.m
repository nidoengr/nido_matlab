function TRACKS = det_tomht_tracks(IMDATA)
ticFun = tic();

nImgs = numel(IMDATA);


numTracks = 900; % Maximum number of tracks, 20 orig
gate = 10;      % Association gate, 45 originally
vol = 1e9;      % Sensor bin volume
beta = 1e-14;   % Rate of new targets in a unit volume
pd = 0.95;       % Probability of detection
far = 1e-8;     % False alarm rate
maxNumHypothesis = 7; % 10

% For reference:
% tracker = trackerTOMHT( ...
%     'FilterInitializationFcn',@initIMMFilter, ...
%     'MaxNumTracks', numTracks, ...
%     'MaxNumSensors', 1, ...
%     'AssignmentThreshold', [0.2, 1, 1]*gate,...
%     'DetectionProbability', pd, 'FalseAlarmRate', far, ...
%     'Volume', vol, 'Beta', beta, ...
%     'MaxNumHistoryScans', 10,...
%     'MaxNumTrackBranches', 5,...
%     'NScanPruning', 'Hypothesis', ...
%     'OutputRepresentation', 'Tracks');

tracker = trackerTOMHT('FilterInitializationFcn',@initcvkf, ...
    'ConfirmationThreshold',20, ...
    'DeletionThreshold',-7, ...
    'AssignmentThreshold', [0.2, 1, 1]*gate,... this is new 0.2 1 1 was the matlab
    'NScanPruning','Hypothesis',... this is new
    'DetectionProbability', pd, 'FalseAlarmRate', far, ... this is new
    'MaxNumHypotheses',maxNumHypothesis,...
    'MaxNumTracks',numTracks);

timeArray_UTC = NaN([nImgs,1]);
idxDateObs = 8;
TRACKS = repelem(struct(...
   'confirmedTracks',objectTrack,...
   'tentativeTracks',objectTrack,...
   'allTracks',objectTrack...
   ),nImgs,1);
for iImg = 1:nImgs
   %    IMDATA(iImg).CentroidsMFHat
   %    idxDateObs = find(contains(IMDATA(1).FITSINFO.PrimaryData.Keywords(:,1),'DATE-OBS'));
   %    if idxDateObs ~= 8
   %       break;
   %    end
   PREstrDateVal = IMDATA(iImg).FITSINFO.PrimaryData.Keywords(idxDateObs,2);
   strDateVal = strrep(PREstrDateVal{:},'T',' ');

%    timeArray(iImg) = juliandate(datetime(datevec(strDateVal))); 
   timeArray_UTC(iImg) = datenum(datevec(strDateVal)); % epoch seconds matlab time
   timeArray_UTC = (timeArray_UTC - timeArray_UTC(1))*86146;

   nDet = height(IMDATA(iImg).CentroidsMFHat);
   detections = cell([nDet,1]);
   for iDet = 1 : nDet
      detections{iDet} = objectDetection(....
         timeArray_UTC(iImg),IMDATA(iImg).CentroidsMFHat(iDet,:),...
         'MeasurementNoise',diag(IMDATA(iImg).CentroidsUncMFHat(iDet,:).^2),... read noise?
         'SensorIndex',1, ...
         'ObjectClassID', iDet,...
         'ObjectAttributes',{struct('ID',iDet)});
   end
%    [confirmedTracks,tentativeTracks,allTracks] = tracker(detections,timeArray(iImg));
   [confirmedTracks,~,~] = tracker(detections,timeArray_UTC(iImg));
   TRACKS(iImg).confirmedTracks = confirmedTracks;
%    TRACKS(iImg).tentativeTracks = tentativeTracks;
%    TRACKS(iImg).allTracks = allTracks;
end

% TRACKS.time_jd = timeArray;

fprintf('Dataset TOMHT elapsed time: %f\n',toc(ticFun));