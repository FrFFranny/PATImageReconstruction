function p0_recon = kwaveAdjoint(f,setting)
% function to calculate the adjoint of forward operator acting on p0 
%
% Inputs:  f: PAT data
%          setting: the setting for k-wave simulation
%            setting.kgrid: the computational grid for k-wave simulation
%            setting.sound_speed: property of the propagation medium
%            setting.medium_density: property of the propagation medium
%            setting.PML_size: the size of perfectly matched layer 
%            setting.sensor_mask: the logical sensor mask
%
% Ouputs:  p0_recon: reconstructed p0
%
% Using cartesian sensor mask to have nice angle ordering. 
% In TR explicitely, in Adjoint through defining time varying sources using cartesian sensor mask. 
 
% define the computation grid for the reconstruction to avoid the
% inverse crime
PML_size = setting.PML_size; % default
kgrid_recon = setting.kgrid;


% attach the original time array
medium.sound_speed = setting.sound_speed;
medium.density = setting.medium_density; % [kg/m^3] default;
[kgrid.t_array, dt] = makeTime(kgrid_recon, medium.sound_speed);
kgrid_recon.t_array = kgrid.t_array;
 
% a special embeding of the adjoint source
source = [];
source.p_mask = setting.sensor_mask;
source.p = [f(:,end:-1:1),zeros(size(f,1),1)] + [zeros(size(f,1),1),f(:,end:-1:1)];
source.p(:,end-1) = source.p(:,end-1) + source.p(:,end);
source.p = source.p(:,1:end-1);
if(isscalar(medium.sound_speed))
    source.p = source.p * (medium.sound_speed * kgrid_recon.dx/(4 * kgrid_recon.dt));
else
    source.p = bsxfun(@times,source.p, medium.sound_speed(source.p_mask)) * kgrid_recon.dx/(4 * kgrid_recon.dt);
end

% rescale from pressure to density
if (isscalar(medium.density))
    source.p = medium.density * source.p;
else
    source.p = bsxfun(@times,source.p, medium.desntiy(source.p_mask));
end
source.p_mode = 'additive';

% we are interested in the pressure at the last time point
sensor = [];
sensor.mask = ones(kgrid_recon.Nx,kgrid_recon.Ny);
sensor.record = {'p','p_final'};

% set the input options
input_args = {'PMLAlpha',2,'PMLSize', PML_size, 'Smooth', false, 'PMLInside', false, 'PlotPML', false,'PlotSim', false};

kWaveResult = kspaceFirstOrder2D(kgrid_recon, medium, source, sensor, input_args{:});

result = gather(kWaveResult.p_final)./(medium.density .* medium.sound_speed.^2);
p0_recon = result;
%p0_recon = smooth(kgrid_recon,result,false,'Blackman');
end