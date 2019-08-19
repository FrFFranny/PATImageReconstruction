clc;
clear;
close all;

% This script is used for reconstructing the PAT images for full or
% subsampled data via FISTA algorithm
%%
% load the original figure
fig = double(rgb2gray(imread('vesselTestImage.png')))/255;
% size of image
[m,n] = size(fig);

% define a disc mask
radius = 100; % grid points
mask = makeDisc(m, n, m/2, n/2, radius);
fig(mask == 0) = 0;

%%
%%%%%%%%%%%%%%% k-wave setting %%%%%%%%%%%%%%%%%%%%

% assign the grid size and create the computational grid
Nx = m;    % number of grid points in the x direction
Ny = n;    % number of grid points in the y direction
x = 10e-3;                  % total grid size [m]
y = x;                  % total grid size [m]
dx = x / Nx;                 % grid point spacing in the x direction [m]
dy = y / Ny;                % grid point spacing in the y direction [m]
setting.kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
setting.sound_speed = 1500;  % [m/s] default
setting.PML_size = 20;
setting.medium_density = 1;

% define equally spaced point sources lying on a circle centred at the
% centre of the detector face
sensor_radius = (radius+10) * dx; % there should be a gap between the circular sensors and the mask;
sensor_angle = 2*pi;

% construct sensor mask
radius_grid_points = round(sensor_radius/setting.kgrid.dx);
binary_sensor_mask = makeCircle(setting.kgrid.Nx,setting.kgrid.Ny,...
    setting.kgrid.Nx/2,setting.kgrid.Ny/2,radius_grid_points,sensor_angle);
setting.sensor_mask = logical(binary_sensor_mask);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% obtain y: the recorded sensor data
y = kwaveForwardSimulation(fig,setting,[]);

% obtain the sensor points and define the subsample of sensor data
sensor_points = size(y,1);
sub_percent = 0.125; % subsample percent (0,1]
sample = fix(sub_percent * sensor_points); % amount of subset data
sample = randperm(sensor_points, sensor_points - sample);

% subsampled y
y(sample,:) = 0;

% obtain Ax based on the kwave forward model
A = @(p0) kwaveForwardSimulation(p0,setting,sample);

% obtain the adjoint operator for kwave model
A_adj = @(sensor_data) kwaveAdjoint(sensor_data,setting);

%% reconstruction via FISTA algorithm

% setting parameters
info.a = 50;
info.tau = 4e-3;
info.radius = '1'; % 1, 2 or Inf
info.Nx = m;
info.Ny = n;
info.maxIter = 100;
info.tol = 1e-4;

% obtain the sigma for Gaussian Filter
width1 = 2; sigma1 = width1 / 2;
width2 = 10; sigma2 =  width2 / 2;

% parameters for FISTA+ADMM algorithm
info.rho = 0.3;
info.tau1 = 4e-3*info.rho;
info.tau2 = 4e-3*info.rho;

% Image reconstruction via FISTA / FISTA + ADMM
[image_recon, normErr] = reconstruction(y,A,A_adj,info,[sigma1,sigma2]);

% apply the mask to the reconstructed image
image_recon(mask == 0) = 0;

%% visualization
figure()
imagesc(fig)
colorbar()
title('original figure')
figure()
imagesc(image_recon)
colorbar()
title(["reconstruction (" + num2str(sub_percent*100) + "% subsampled)"])

