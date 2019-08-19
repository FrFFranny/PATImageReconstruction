function [p0_recon] = kwaveTRReconstruction(sensor_data,mask,info)
% Implement the k-wave time-reversal reconstruction

% Input:
%   sensor_data: the recorded sensor data
%   mask: the mask for original image
%   info: the setting for k-wave simulation
%            info.kgrid: the computational grid for k-wave simulation
%            info.sound_speed: property of the propagation medium
%            info.medium_density: property of the propagation medium
%            info.PML_size: the size of perfectly matched layer 
%            info.sensor_mask: the logical sensor mask

% Output: the time reversal reconstructed result

% define the computation grid for the reconstruction to avoid the
% inverse crime
kgrid_recon = info.kgrid;

% define the properties of the propagation medium
medium.sound_speed = info.sound_speed;  % [m/s]

% use the same time array for the reconstruction
kgrid_recon.setTime(info.kgrid.Nt, info.kgrid.dt);

% set the initial pressure
source.p0 = 0;

% assign the time reversal data
sensor.time_reversal_boundary_data = sensor_data;

sensor.mask = info.sensor_mask;

% set the input options
input_args = {'Smooth', false, 'PMLInside', false, 'PlotPML', false,'PlotSim', false};

% % run the time-reversal reconstruction
p0_recon = kspaceFirstOrder2D(kgrid_recon, medium, source, sensor, input_args{:});

% add the mask and the non negative condition
p0_recon(mask == 0) = 0;
p0_recon(p0_recon < 0) = 0;
end

