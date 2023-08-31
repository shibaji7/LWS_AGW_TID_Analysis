%% Description
%
% Author : Chakraborty, Shibaji
% Name : 3D Ray trace model
%
% Purpose :
%   Estimate the SuperDARN radar rays (for only one beam) due to
%   change in ionospheric conditions during solar flare events using the
%   the spherical Earth 3D NRT raytrace_2d_sp for a fan of rays.
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

speed_of_light = 2.99792458e8;
load 
