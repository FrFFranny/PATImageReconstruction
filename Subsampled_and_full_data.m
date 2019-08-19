clc;
clear;
close all;

% This script is used for comparing the performance of image reconstruction
% based on subsampled data and full data via the time reversal operator in
% PAT and the results of adjoint operator are also presented.

%% simulation with mask

% load the test image
p0 = double(rgb2gray(imread('vesselTestImage.png')))/255;
[m,n] = size(p0);

% define a disc mask
radius = 100; % grid points
mask = makeDisc(m, n, m/2, n/2, radius);
p0(mask == 0) = 0;

% visulize the original figure
figure()
imagesc(p0)
title("original figure")


%%%%%%%%%%%%%%% k-wave setting %%%%%%%%%%%%%%%%%%%%

% assign the grid size and create the computational grid
[m,n] = size(p0); % size of image
Nx = m;    % number of grid points in the x direction
Ny = n;    % number of grid points in the y direction
x = 10e-3;                  % total grid size [m]
y = x;                  % total grid size [m]
dx = x / Nx;                % grid point spacing in the x direction [m]
dy = y / Ny;                % grid point spacing in the y direction [m]
setting.kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
setting.sound_speed = 1500;  % [m/s] default
setting.PML_size = 20;
setting.medium_density = 1;

% define equally spaced point sources lying on a circle centred at the
% centre of the detector face
sensor_radius = (radius + 10) * dx;
sensor_angle = 2*pi;

% construct sensor mask
radius_grid_points = round(sensor_radius/setting.kgrid.dx);
binary_sensor_mask = makeCircle(setting.kgrid.Nx,setting.kgrid.Ny,...
    setting.kgrid.Nx/2,setting.kgrid.Ny/2,radius_grid_points,sensor_angle);
setting.sensor_mask = logical(binary_sensor_mask);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% implement the k-wave forward model
sensor_data = kwaveForwardSimulation(p0,setting,[]);


%% visualization the full sensor data
temp.mask = setting.sensor_mask;
reorder = reorderSensorData(setting.kgrid,temp,sensor_data);
figure();
imagesc(reorder)
title("sensor data - full data")

%% subsample

sub_percent = 0.25; % subsample percent
sensor_points = size(sensor_data,1); % obtain the sensor points
sample = fix(sub_percent * sensor_points); % amount of subset data
sample = randperm(sensor_points, sensor_points - sample);
subsampled_data = sensor_data;
subsampled_data(sample,:) = 0;

% visualization the subsampled sensor data
temp.mask = setting.sensor_mask;
reorder_subsample = reorderSensorData(setting.kgrid,temp,subsampled_data);
figure();
imagesc(reorder_subsample)
title([num2str(sub_percent) + " subsampled data"])


%% k-wave reconstruction

% full data reconstruction
recon_fulldata = kwaveTRReconstruction(sensor_data,mask,setting);
figure();
imagesc(recon_fulldata)
title("Time reversal - full data")

%%
% subsampled data reconstruction
recon_subsample = kwaveTRReconstruction(subsampled_data,mask,setting);
figure();
imagesc(recon_subsample)
title(["Time reversal - "+num2str(sub_percent) + " subsampled data"])

%% k-wave adjoint results

% obtain the adjoint result based on the kwave adjoint operator
adj_fulldata = kwaveAdjoint(sensor_data,setting);
adj_fulldata(adj_fulldata<0) = 0;
adj_fulldata(mask == 0) = 0;

adj_subsampleddata = kwaveAdjoint(subsampled_data,setting);
adj_subsampleddata(adj_subsampleddata<0) = 0;
adj_subsampleddata(mask == 0) = 0;

% visualization
figure();
imagesc(adj_fulldata)
title("Adjoint - full data")

figure();
imagesc(adj_subsampleddata)
title(["Adjoint - " + num2str(sub_percent) + " subsampled data"])
