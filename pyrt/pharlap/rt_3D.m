%% Description
%
% Author : Chakraborty, Shibaji
% Name : 3D Ray trace model
%
% Purpose :
%   Estimate the SuperDARN radar rays (for only one beam) due to
%   change in ionospheric conditions during solar flare events using the
%   the spherical Earth 3D NRT raytrace_3d_sp for a fan of rays.
%   Ray trajectories are saved under 'local' directory
%
% Calling sequence :
%   rt_3d
%
% Inputs :
%   Radar code
%   Time index
%   Beams num
%   Directory
%   
% Outputs :
%   .mat files with ray trajectory

%% Close and clear all panes and variables and command line window and load PharLap
close all
startup_3D

load(char(strcat('../', dic, 'bearing_', compose('%02d', bm), '.mat')))
load(char(strcat('../', dic, compose('%02d',UT(4)), '.', compose('%02d',UT(5)), '_', compose('%02d', bm), '.mat')))

elevs = double(elevs);
ray_bears = zeros(size(elevs));
iono_grid_parms = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
      ht_start, ht_inc, num_ht, ];
      
B_num_ht = ceil(num_ht .* ht_inc ./ B_ht_inc);
B_num_lat = ceil(num_lat .* lat_inc ./ B_lat_inc);
B_num_lon = ceil(num_lon .* lon_inc ./ B_lon_inc);
geomag_grid_parms = [B_lat_start, B_lat_inc, B_num_lat, B_lon_start, ...
      B_lon_inc, B_num_lon, B_ht_start, B_ht_inc, B_num_ht];

tic
fprintf('Generating ionospheric and geomag grids... ')
[iono_pf_grid, iono_pf_grid_5, collision_freq, Bx, By, Bz] = ...
    gen_iono_grid_3d(UT, R12, iono_grid_parms, ...
                     geomag_grid_parms, doppler_flag);
toc
fprintf('\n')

% convert plasma frequency grid to  electron density in electrons/cm^3
%iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
%iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;

tic
[ray_data_O, ray_O, ray_state_vec_O] = ...
  raytrace_3d(origin_lat, origin_long, origin_ht, elevs, ray_bears, freqs, ...
              OX_mode, nhops, tol, iono_en_grid, iono_en_grid_5, ...
	          collision_freq, iono_grid_parms, Bx, By, Bz, ...
	          geomag_grid_parms);
	      
NRT_total_time = toc;
fprintf('\n   NRT-only execution time = %f, Total mex execution time = %f\n\n', ...
        [ray_data_O.NRT_elapsed_time], NRT_total_time)
for rayId=1:length(elevs)
  num = length(ray_O(rayId).lat);
  ground_range = zeros(1, num);
  lat = ray_O(rayId).lat;
  lon = ray_O(rayId).lon; 
  ground_range(2:num) = latlon2raz(lat(2:num), lon(2:num), origin_lat, ...
      origin_long,'wgs84')/1000.0;
  ray_O(rayId).ground_range = ground_range;
end

fprintf('\n   Save to file = %s \n\n', char(strcat('../', fname)))
save(char(strcat('../', fname)), 'ray_data_O', 'ray_O', 'ray_state_vec_O');