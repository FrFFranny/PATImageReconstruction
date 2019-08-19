function [sensor_data] = kwaveForwardSimulation(p0,setting,sample)
% Use k-wave toolbox to implement the time-domain forward model of acoustic
% wave propagation for acoustically heterogeneous media with power law 
% absorption.

%   p0: the original image
%   setting: the setting for k-wave simulation
%            setting.kgrid: the computational grid for k-wave simulation
%            setting.sound_speed: property of the propagation medium
%            setting.medium_density: property of the propagation medium
%            setting.PML_size: the size of perfectly matched layer 
%            setting.sensor_mask: the logical sensor mask
%   sample: a vector for setting selected rows of output sensor data to be
%           zeros based on the entries of 'sample' vector

% define the computational grid
kgrid = setting.kgrid;

% define the properties of the propagation medium
medium.sound_speed = setting.sound_speed;  % [m/s]
medium.density = 1;
% resize the input image to the desired number of grid points
p0 = resize(p0, [kgrid.Nx, kgrid.Ny]);

% assign to the source structure
source.p0 = p0;

% define sensor structure
sensor.mask = setting.sensor_mask;

% create the time array
kgrid.makeTime(medium.sound_speed);

% assign the parameters
medium.sound_speed = setting.sound_speed;
medium.density = setting.medium_density;

% set the input options
input_args = {'Smooth', false, 'PMLSize',setting.PML_size,'PMLInside', false, 'PlotPML', false,'PlotSim', false};

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% add noise to the recorded sensor data
signal_to_noise_ratio = 40;	% [dB]
sensor_data = addNoise(sensor_data, signal_to_noise_ratio, 'peak');

sensor_data(sample,:) = 0;
end

