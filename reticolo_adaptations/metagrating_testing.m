%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE 1D (TE or TM) %
%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
close all;

wavelength=8;
period=10;% same unit as wavelength
n_air=1;% refractive index of the top layer
n_Si=1.5;% refractive index of the bottom layer

angle_theta0=-10;
k_parallel=n_air*sin(angle_theta0*pi/180);

parm=res0(1);% TE polarization. For TM : parm=res0(-1)
parm.res1.champ=1;% the electromagnetic field is calculated accurately

nn=40;% Fourier harmonics run from [-40,40]

% textures for all layers including the top and bottom layers
textures=cell(1,3);
textures{1}= n_air; % uniform texture
textures{2}= n_Si; % uniform texture
textures{3}={[period,period + 0.1],[n_air,n_Si] };

aa= res1(wavelength,period,textures,nn,k_parallel,parm);

profile ={[4.1,5.2,4.1],[1,3,2]};

one_D_TE=res2(aa,profile);
eff=one_D_TE.inc_top_reflected.efficiency{-1};
J=one_D_TE.Jones.inc_top_reflected{-1};% Jones’coefficients
abs(J)^2 % first order efficiency for an illumination from the top layer

% field calculation
x=linspace(-period/2,period/2,51);% x coordinates(z-coordinates are determined by res3.m)
einc=1;
parm.res3.trace=1; % plotting automatically
parm.res3.npts=[50,50,50];
[e,z,index]=res3(x,aa,profile,einc,parm);
figure;pcolor(x,z,real(squeeze(e(:,:,1)))); % user plotting
shading flat;xlabel('x');ylabel('y');axis equal;title('Real(Ey)');

% Loss calculation
textures{3}={[-2.5,2.5],[n_air,.1+5i] };
aa_loss=res1(wavelength,period,textures,nn,k_parallel,parm);
one_D_loss=res2(aa_loss,profile);
parm.res3.npts=[[0,10,0];[1,3,1]];
einc=one_D_loss.inc_top.PlaneWave_E(2);
[e,z,index,wZ,loss_per_layer,loss_of_Z,loss_of_Z_X,X,wX]=res3([-period/2,period/2],aa_loss,profile,einc,parm);
energy_conservation = sum(one_D_loss.inc_top_reflected.efficiency)+sum(one_D_loss.inc_top_transmitted.efficiency)+sum(loss_per_layer)/(.5* period)-1

retio % erase temporary files

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1D Metagrating Design - unit cell example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
close all;

% light parameters
desired_angle = 60; % deg, in air
incident_angle = 0; % deg; direction of incoming light
wavelength = 8;
% wavelength = 900; % nm
period = wavelength / sind(desired_angle);
% period = 10; % same unit as wavelength

thickness = 4; % nm; thickness of silicon grating
z_step = min(wavelength) / 40; % z step for z discritization
num_z_disc = round(thickness / z_step);

n_air = 1; % refractive index of air
n_Si = 3.45; % refractive index of silicon
n_sub = 1.45; % refractive index of SiO_2

k_par_f = sind(incident_angle) * n_sub; % normalized parallel incident wave vector

parm = res0(-1); % TM polarization; For TE: parm = res0(1)
parm.res1.champ = 1;% the electromagnetic field is calculated accurately
parm.res3.sens = -1; % the grating is illuminated from the *bottom* (direction of incident plane wave)
parm.res3.trace = 1; % automatic plot (shows ALL the components of the EM fields and grating refractive index distribution)
parm.res3.npts = [0,num_z_disc, 0];
% parm.res3.npts = [100,100,100];

nn = ceil((12 * period) / min(wavelength)); % defines the set of Fourier harmonics retained for the computation

% % textures for all layers including the top and bottom layers
textures = cell(1,5);
% number of layers 
textures{1} = n_incident_medium; % uniform texture
textures{2} = {[-2.1,2.9],[n_air,n_Si]};
textures{3} = {[-2.3,2.7],[n_air,n_Si]};
textures{4} = {[-2.5,2.5],[n_air,n_Si]};
textures{5} = n_substrate; % uniform texture

% % textures for all layers including the top and bottom layers
% textures = cell(1,3);
% textures{1} = n_air; % uniform texture
% textures{2} = {[-2.5,2.5],[n_air,n_Si]}; % 1st vector contains all the x-values of the discontinuities, 2nd vector contains the refractive indices of the material between the discontinuities
% textures{3} = n_sub; % uniform texture

aa = res1(wavelength,period,textures,nn,k_par_f,parm);

% profile = {[4.1,thickness,1.2],[1,2,3]};
profile = {[0,thickness,0],[1,2,3]};

result1 = res2(aa,profile);

% extract efficiencies for the +1 diffraction order
eff_transmitted = result1.inc_bottom_transmitted.efficiency{1};
eff_reflected = result1.inc_bottom_reflected.efficiency{1};

% field calculation
x = linspace(-period/2,period/2,51);% x coordinates(z-coordinates are determined by res3.m)
einc = 1;
[e,z,index] = res3(x,aa,profile,einc,parm);
figure;pcolor(x,z,real(squeeze(e(:,:,1)))); % user plotting
shading flat;xlabel('x');ylabel('y');axis equal;title('Real(Ey)');

retio % erase temporary files

%% 1D Metagrating Design - GLOnet solver analog
% Nathan Miller
% 7/13/2021

clear;
clc;
close all;

% input image training data
% img = [0,0,0,1,1,1,1,0,0,1,1,0,1,1,1];
img = [0,0,0,1,1,1,1,0,0,1,1,0,1,1,1,0,1,0,1,1,1,0,0,0,0,0,1];

% light parameters
desired_angle = 60; % deg, in air
incident_angle = 0; % deg; direction of incoming light
wavelength = 900; % nm
period = wavelength / (sind(desired_angle));

% grating thickness - to be optimized later
thickness = 325; % nm; thickness of silicon grating
layers = 1:10; % number of layers in grating
z_step = min(wavelength) / 40; % z step for z discretization
num_z_disc = ones(1,(size(layers,2) - 2)); % exclude top and bottom layers (only consider the silicon)
num_z_disc = num_z_disc * (round(thickness / z_step)); % z plane locations of the computed fields for silicon

% material and fluid properties
n_air = 1; % refractive index of air
n_Si = 3.45; % refractive index of silicon
n_sub = 1.45; % refractive index of SiO_2

k_par_f = sind(incident_angle) * n_sub; % normalized parallel incident wave vector

parm = res0(-1); % TM polarization; For TE: parm = res0(1)
parm.res1.champ = 1;% the electromagnetic field is calculated accurately
parm.res3.sens = -1; % the grating is illuminated from the *bottom* (direction of incident plane wave)
parm.res3.trace = 1; % automatic plot (shows ALL the components of the EM fields and grating refractive index distribution)
parm.res3.npts = [0,num_z_disc,0]; % z locations of the computed fields, vector whose length is equal to the number of layers. Default is 10 z = constant planes per layer

% Fourier harmonics
nn = ceil((12 * period) / min(wavelength)); % defines the set of Fourier harmonics retained for the computation

% side wall tilt angle parameters
theta_tilt = 20; % deg
d = 20; % distance each layer is shifted by
t_layer = thickness / numel(layers); % thickness of each layer
t_layer = t_layer * ones(1,(size(layers,2) - 2));
% theta_tilt = atand(d / t_layer);

% textures for all layers
textures = cell(1,numel(layers));

% top layer (air)
textures{1} = n_air; % uniform texture

% middle layers (silicon
nlength = length(img);
dx = period / nlength;
xvec = ([1:nlength] * dx) - (0.5 * period);
nvec = (img * (n_Si - n_air)) + n_air;

% for i = 2:(n_layers(end) - 1)
% textures{i} = {xvec + ((i - 1) * d),nvec}; % for tilt to the left
% end

for i = (layers(end) - 1):-1:2
textures{i} = {xvec - ((i - 1) * d),nvec}; % tilt to the right requires starting from the top-most layer to the bottom and shifting to the left
end

% bottom layer (SiO_2)
textures{layers(end)} = n_sub; % uniform texture

aa = res1(wavelength,period,textures,nn,k_par_f,parm); % contains all the information on the eigenmodes of all textures

profile = {[0,t_layer,0],layers};

result1 = res2(aa,profile);

% extract efficiencies for the +1 diffraction order
eff_transmitted = result1.inc_bottom_transmitted.efficiency{1};
eff_reflected = result1.inc_bottom_reflected.efficiency{1};

% field calculation
x = linspace(-period/2,period/2,51);% x coordinates(z-coordinates are determined by res3.m)
einc = 1;
[e,z,index]=res3(x,aa,profile,einc,parm);
figure;pcolor(x,z,real(squeeze(e(:,:,1)))); % user plotting
shading flat;xlabel('x');ylabel('y');axis equal;title('Real(Ey)');

retio % erase temporary files