% metagrating_design.m  
% 1D Metagrating Design - GLOnet solver analog assuming given img and known 
%                         desired tilt angle, wavelength, diffraction angle,
%                         and thickness
% ------------------------------------------------------------------------
% Inputs:
%  
% Outputs:
% 
% Dependencies:
% 
% Notes:
%  * Reference code adapted from Jiaqi Jiang & Jonathan A. Fan, 
%    Simulator-based training of generative models for the inverse design
%    of metasurfaces
%
% Author: Nathan Miller
% 
% Date: 2021-07-13
% ------------------------------------------------------------------------

clear;
clc;
close all;

%% Forward Simulation

% input image training data
img = [0,0,0,1,1,1,1,1,1,1,0,1,1,1,0,1,1,0,0,1,1,1,1,1,0,0,0];

% operations on input image training data
% img = (img / 2.0) + 0.5; % push devices to possess binary refractive index values in the vector x ∈ {−1,1}^N_segments, where −1 represents air and +1 represents silicon.

% minimum feature size
min_feat = 50; % nm

% light parameters
desired_angle = 60; % deg, in air
incident_angle = 0; % deg; direction of incoming light
wavelength = 900; % nm
period = wavelength / (sind(desired_angle));
diff_order = 1; % target diffraction order

% grating thickness - to be optimized later
thickness = 325; % nm; thickness of silicon grating

% define array of layers
layers = 1:100; % layers in grating

% side wall tilt angle parameters
theta_tilt = 20; % deg

if abs(theta_tilt) <= 1
    theta_tilt = 0;
end

t_layer = thickness / numel(layers); % thickness of each layer
d = t_layer * tand(theta_tilt); % distance each layer is shifted by
t_layer = t_layer * ones(1,(size(layers,2) - 2));

% z sampling parameters
z_step = min(wavelength) / 40; % z step for z discretization
num_z_disc = ones(1,(size(layers,2) - 2)); % exclude top and bottom layers (only consider the silicon)
num_z_disc = num_z_disc * (ceil(t_layer(1) / z_step)); % z plane locations of the computed fields for silicon divided into each layer width

% material and fluid properties
n_air = 1; % refractive index of air
n_Si = 3.45; % refractive index of silicon
n_sub = 1.45; % refractive index of SiO_2

k_par_forw = sind(incident_angle) * n_sub; % normalized parallel incident wave vector for forward simulation

parm = res0(-1); % TM polarization; For TE: parm = res0(1); default parameters based on polarization type
parm.res1.champ = 1; % the electromagnetic field is calculated accurately
parm.res3.sens = -1; % the grating is illuminated from the *bottom* (direction of incident plane wave)
parm.res3.trace = 1; % automatic plot (shows ALL the components of the EM fields and grating refractive index distribution)
parm.res3.npts = [0,num_z_disc,0]; % z locations of the computed fields, vector whose length is equal to the number of layers. Default is 10 z = constant planes per layer
parm.res3.champs = 0; % only produce the refractive index plot in automatic plotting

% Fourier harmonics
nn = ceil((12 * period) / min(wavelength)); % defines the set of Fourier harmonics retained for the computation

% erase temporary files
retio([],inf*1i);

% create textures for all layers
textures = cell(1,numel(layers));

% top layer (air)
textures{1} = n_air; % uniform texture

% middle layers (silicon)
nlength = length(img);
dx = period / nlength;
xvec = ([1:nlength] * dx) - (0.5 * period);
nvec = (img * (n_Si - n_air)) + n_air;

for i = (layers(end) - 1):-1:2
textures{i} = {xvec - ((i - 1) * d),nvec}; % tilt to the right requires starting from the top-most layer to the bottom and shifting to the left
end

% bottom layer (SiO_2)
textures{layers(end)} = n_sub; % uniform texture

% calculating the eigenmodes associated to all textures
aa = res1(wavelength,period,textures,nn,k_par_forw,parm); % contains all the information on the eigenmodes of all textures

profile = {[0,t_layer,0],layers}; % contains successive thickness and texture-label information relative to every layer, starting from the top layer and finishing at the bottom layer

% computing the diffracted waves
result1 = res2(aa,profile); % contains all the information on the diffracted fields (waves)
wave_results_forw = result1.inc_bottom_transmitted; % contains all the information concerning the propagative reflected waves for an incident wave from the bottom layer of the grating.

% extract deflection efficiency for the +1 diffraction order
eff_transmitted = wave_results_forw.efficiency(diff_order);

% field calculation
n_points_x = round(period / (d / 2)); % fine enough to know the layer is shifting
x = linspace(-period/2,period/2,n_points_x); % x coordinates (z-coordinates are determined by res3.m)
einc = 1; % defines the y component of the complex amplitude of the incident electric (in TE polarisation) or magnetic field (in TM polarisation) field
[e_forw,~,~] = res3(x,aa,profile,einc,parm); % The variable “e” contains all the electromagnetic field quantities, “z” is the vector containing the z-coordinate of the sampling points, “index” is the complex refractive index of the considered grating

%% Backward (Adjoint) Simulation

k_par_adj = -wave_results_forw.K(diff_order,1); % normalized parallel incident wave vector for adjoint simulation

% calculating the eigenmodes associated to all textures
aa = res1(wavelength,period,textures,nn,k_par_adj,parm);

% computing the diffracted waves
result2 = res2(aa,profile);
wave_results_adj = result2.inc_top_transmitted; % contains all the information concerning the propagative reflected waves for an incident wave from the top layer of the grating.

parm.res1.champ = 1; % the electromagnetic field is calculated accurately
parm.res3.sens = 1; % the grating is illuminated from the top

%% FOM gradient (Figure of Merit)

% for TM polarized light:
[e_adj,~,~] = res3(x,aa,profile,exp(1i*angle(conj(wave_results_forw.H(diff_order,2)))),parm);
EEpar = e_forw(:,:,3) .* e_adj(:,:,3);
EEnorm = e_forw(:,:,2) .* e_adj(:,:,2);
gr0all = real(1i * (EEpar + EEnorm));

% image(gr0all,'CDataMapping','scaled')
% colorbar

% compute average electric field in each refractive index slice (at evenly spaced xvec)
% n_points_per_slice = round(length(x) / nlength)

% mean(gr0all(textures{2:layers(end)}

grall = mean(gr0all,1) .* nvec;
grtot = grall;
grtot(1:round(min_feat / 2)) = 0;
grtot((nlength - round(min_feat / 2)):nlength) = 0;
grtot((img == 1) & (grtot > 0)) = 0;
grtot((img == 0) & (grtot < 0)) = 0;

grtot = grtot/max(abs(grtot));

% gradient output
Gr = grtot * 2.0;

%%
image(gr0all,'CDataMapping','scaled')
colorbar