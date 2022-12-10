%ASEN6084_PROJ_SIM.M is a top level script that reads the image data set, applies thresholding and clustering. TODO! implement TO-MHT, photometric centroid function,
clear all; close all; clc;
disp('Final Project Script');
ticSim = tic();


% Location of sensor: CU Boulder Aerospace Building (lat lon alt)=(deg deg m)
sen_LLA_deg_deg_m = [40.010051204761496, -105.2438703830158, 1655];

%% Read Data
IMDATA = read_project_data_geo_dataset();

%% Process Images and Extract Detections
IMDATA = image_proc_and_det(IMDATA);

%%  TO-HOMHT
TRACKS = det_tomht_tracks(IMDATA);
plot_confirmed_tracks_over_time(TRACKS)
plot_meas_over_time(IMDATA)
%%
tocSim = toc(ticSim);
fprintf('Simulation elapsed time (s): %f\n',tocSim);