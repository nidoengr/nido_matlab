function [XDATA, PDATA, residDATA, idxMeasAvailable] = ckf_process(tArray_s, deltax0, P0, Xref0,...
   muEarth_km3ps2, J2, J3, GsLatLon, GMT0_deg, Y, R, opt, sigvdot_kmps2)
%ckf_process is a function for the sequential estimation algorith
%   aka the classical kalman filter.
%
% Input:
% tArray_s  (1,nTimeElms) the time array of specific elements
% xbar0     (6,1)   An estimate of the state error [x y z vx vy vz]
% P0        (6,6)   The state covariance matrix
% Xref0     (6,1)   The initial reference state
% muEarth_km3ps2 (1,1) The Earth gravitation parameter
% J2        (1,1)
% J3        (1,1)
% GSlatlon  (nStations,2)
% GMT0_det  (1,1)   Initial greenwhich time in degrees
% Y         (nMeas,nTimeElms)  A measurement matrix with
%                 measurements in each row corresponding to the
%                 time column.
% R         (nMeas,nMeas) Measurement covariance matrix
% opt       ( . )   Integration options class
%
%%
if nargin < 13
   sigvdot_kmps2 = 0;
end
%% Constants
I3x3 = eye(3);
I6x6 = eye(6);

%% Preprocess inputs
nMeasPerFrame = numel(Y(:,1));

idxMeasAvailable = ~isnan(Y(1,:)) | ...
   ~isnan(Y(2,:)) | ...
   ~isnan(Y(3,:)) | ...
   ~isnan(Y(4,:)) | ...
   ~isnan(Y(5,:)) | ...
   ~isnan(Y(6,:));

% idxMeasAvailable = false([1,numel(Y(1,:))]);
% for idxY = 1:numel(Y(1,:))
%    if ~all(isnan(Y(:,idxY))) %any(~isnan(Y(:,idxY)))
%       idxMeasAvailable(idxY) = true;
%    end
% end

tMeasAvail_s = tArray_s(idxMeasAvailable);
YMeasAvail = Y(:,idxMeasAvailable);

nTimeElms = numel(tMeasAvail_s);
% nStates = numel(X0ref0);

% nTimeElms = numel(tArray_s);
nStates = numel(deltax0);
nGS = numel(GsLatLon(:,1));

% STM0 = eye(nStates);
% STM0_flat = reshape(STM0,nStates^2,1);

% Use one with an extra state for J2 to use keplerJ2_wPhi_ODE
% function [x y z xdot ydot zdot J2 STM_FLAT] (56,1)
STM0wJ2 = eye(nStates+1);
STM0wJ2_flat = reshape(STM0wJ2,(nStates+1)^2,1);

Z0wJ2STM = [Xref0;J2;STM0wJ2_flat];

%% Initialize data storage
ZDATAwJ2 = nan([numel(Z0wJ2STM), nTimeElms]);
PDATA = nan([nStates,nStates,nTimeElms]);
XDATA = nan([nStates,nTimeElms]);
xDATA = nan([nStates,nTimeElms]);
residDATA = nan([nMeasPerFrame,nTimeElms]);

% Assign 1st states
ZDATAwJ2(:,1) = Z0wJ2STM;
PDATA(:,:,1) = P0;
XDATA(:,1) = Xref0;
xDATA(:,1) = deltax0;

%%

for k = 2 : nTimeElms
   % Assign previous index
   tkm1 = tMeasAvail_s(k-1);
   tk   = tMeasAvail_s(k);
   dT = tk-tkm1;
   deltaxkm1 = xDATA(:,k-1);
   Pkm1 = PDATA(:,:,k-1);
   Zkm1wJ2 = ZDATAwJ2(:,k-1);
   Xkm1 = Zkm1wJ2(1:6,1);
   
   % Propagate reference state (non-linear dynamics)
   % Integrate from tkm1 to tk
   [~, ZkwJ2temp] = ode113(@keplerJ2_wPhi_ODE,[tkm1 tk],Zkm1wJ2,opt);
   ZkwJ2 = ZkwJ2temp(end,:)';
   Xkm = ZkwJ2(1:6,1);
   
   % Pull out STM(tk,t0)
   STMkwJ2_flat = ZkwJ2(8:end,1);
   STMkwJ2 = reshape(STMkwJ2_flat,nStates+1,nStates+1);
   STMk = STMkwJ2([1:6],[1:6]);
   
   %
   Q = sigvdot_kmps2^2 .* [1/3*dT^3*I3x3, 1/2*dT^2*I3x3
      1/2*dT^2*I3x3, dT*I3x3];

   deltaxkm = STMk*deltaxkm1;
   Pkm = STMk * Pkm1 * STMk' + Q; % + Proc Noise
   
   % Process measurements in pairs of rho rhodot
   nMeas = numel(YMeasAvail(:,k));
   
   % Grab the measurements that are not nans
   idxYValid = ~isnan(YMeasAvail(:,k));
   if all(~idxYValid), keyboard, end
   
   % Get indicies of coresponding ground stations
   idxGSValid = idxYValid(1:2:nMeas-1);
   
   % Count number of valid measurements
   nValidMeas = sum(idxYValid);
   
   % Assign this step measurement vector
   yk = YMeasAvail(idxYValid,k);
   GSValidLatLon = GsLatLon(idxGSValid,:);
   
   % Grab the measurement variance
   Rk = R(idxYValid,idxYValid);
   
   if mod(nValidMeas,2)~=0
      error('Measurement at k are not pairs');
   else
      nMeasPairs = nValidMeas/2;
   end
   
   % Initialization predicted measurement vector
   ybark = nan(size(yk)); % Make same size as nValidMeas
   
   % Initialize sensitivity matrix
   Hk = nan(nValidMeas,nStates);
   
   % Predict measurement
   idxyk = [1 2];
   for iMeasPair = 1:nMeasPairs
      
      % Get observer state in ECI
      XkObs_ECI = lat_lon_to_ECI(GSValidLatLon(iMeasPair,1),...
         GSValidLatLon(iMeasPair,2),GMT0_deg, tk);
      
      % Measurement partials with respect to state
      Htilde = Htilde_sc_rho_rhod(Xkm, XkObs_ECI);
      Hk(idxyk,:) = Htilde;
      
      % predicted meas
      ybark(idxyk,1) = compute_range_rangerate(Xkm,...
         GSValidLatLon(iMeasPair,1),GSValidLatLon(iMeasPair,2), GMT0_deg, tk)';
      
      % update indices
      idxyk = idxyk + 2;
      
   end
   yresidk = yk - ybark;
   residDATA(idxYValid,k) = yresidk;
   
   K = Pkm*Hk' / ( Hk * Pkm * Hk' + Rk);
   
   deltaxk = deltaxkm + K*(yresidk - Hk*deltaxkm);
   Pk = (eye(6) - K*Hk)*Pkm*(eye(6) - K*Hk)' + K*Rk*K';
   
   Xk = Xkm + deltaxk;

   % Assing new data storage
   ZDATAwJ2(:,k) = [Xkm;J2;reshape(eye(nStates+1),(nStates+1)^2,1)];
   
   % Assing Pk and xhatk to data storage
   PDATA(:,:,k) = Pk;
   xDATA(:,k) = deltaxk;
   XDATA(:,k) = Xk;
   
end
% XDATA = ZDATAwJ2(1:6,:);