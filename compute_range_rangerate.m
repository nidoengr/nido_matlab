function obs = compute_range_rangerate(x_sc, latitude, longitude,...
   theta0, obsTime, includeElevationFlag)
%
%
% Inputs:
% x_sc        (6,1)   ECI position velocity of SC
% latitude    (1,1)   ground station latitude in degrees
% longitude   (1,1)   ground station longitude in degress
% theta0      (1,1)   the initial earth rotation angle in degrees
%                     at t=0 
% obsTime     (1,1)   time since epoch in seconds
% obs         (1,2)   measurement vector with range range rate in
%                     km and km/s

if nargin<6
   includeElevationFlag = false;
end

% Unpack:
%--------
rSc_ECI_km = x_sc(1:3,1);
vSc_ECI_kmps = x_sc(4:6,1);

% Assumptions:
%-------------
% Earth rotation rate
wdot = 7.2921158553E-5; % Earth Rotation Rate, rad/s

% spherical earth
RE_km = 6378;

% Transformation matrix from lat lon to ECEF
lat_rd = latitude * pi/180;
lon_rd = longitude * pi/180;

clat = cos(lat_rd);
clon = cos(lon_rd);

slat = sin(lat_rd);
slon = sin(lon_rd);

T_LL2ECEF = [ clat*clon, -slat*clon, -slon;...
              clat*slon, -slat*slon,  clon;...
                   slat,       clat,     0];
                
% Transfomation matrix from ECI 2 ECEF
alphaG_rd = wdot*obsTime + (theta0*pi/180);

c = cos( alphaG_rd );
s = sin( alphaG_rd );

T_ECEF2ECI = [c, -s, 0; s, c, 0; 0, 0, 1];

% Ground station location in ECEF
rGs_ECEF_km = T_LL2ECEF * [RE_km;0;0];
% rGs_ECEF_km = RE_km * [clat*clon;clat*slon;slat];
vGs_ECEF_kmps = [0;0;0];

rGs_ECI_km = T_ECEF2ECI * rGs_ECEF_km;
vGs_ECI_kmps = cross( [0;0;wdot], rGs_ECI_km );

% Alternatively
rRel_ECI_km   = rSc_ECI_km - rGs_ECI_km;
vRel_ECI_kmps = vSc_ECI_kmps - vGs_ECI_kmps;

% rho_km = sqrt( sum( rRel_ECI_km.^2 ) );
rho_km = norm( rRel_ECI_km );

rhodot_kmps = dot( rRel_ECI_km, vRel_ECI_kmps )/rho_km;

% Check elevation angle
rGs_UECI = rGs_ECI_km./norm(rGs_ECI_km);
rRel_UECI = rRel_ECI_km./rho_km;

elFromSite_deg = 90 - acosd( dot(rGs_UECI,rRel_UECI));

if includeElevationFlag
   obs = [ rho_km, rhodot_kmps, elFromSite_deg ];
else
   obs = [ rho_km, rhodot_kmps ];
end
