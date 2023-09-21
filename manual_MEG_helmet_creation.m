%%%%%%%%%%%%%%%%%%
%code from field trip
%"https://www.fieldtriptoolbox.org/example/compute_forward_simulated_data_and_apply_a_beamformer_scan/"
% modified to create a full 360 degree array of sensors
%using some of Wanjin's code "CC_main.m" and "gensphere.mat"
%%%%%%%%%%%%%%%%%%%
clear
import mne.*

%% create spherical spacing and setup
%meg_geometry=load("meg_geometry.mat");
spacing = 18; 
radius = 0.12; %0.12 m spherical surface
[xsens,ysens,zsens] = gensph(radius,spacing); % spacing^2 # sensors 
ch_types = ones(size(xsens,1),1); % 1 = magnetometers, 0 = gradiometers for S_in code
nchan=size(xsens,1);
d = (21e-3)/2; % square sensor with width d/2
R = [xsens, ysens, zsens]';
RT=transpose(R);

%% sensor orientations
% can randomise normal sensor orientations by up to 0.05 m in each axis to ensure S_in is lin indep
%normr(R(:,k)' + 0.05*normr(2*rand(1,3)-1));
RR = [];
for k = 1:size(R,2) 
    RR(k,:) = normr(R(:,k)'); 
end
for k = 1:size(R,2) % define the other (local) sensor orientations
    nullspace = null(RR(k,:));    
    EX(k,:) = nullspace(:,1); 
    EY(k,:) = nullspace(:,2);
    EZ(k,:) = normc(RR(k,:)'); 
end

%% create grad struct
%R=sensor positions, EZ=normal orientations for grad.coilori
grad = [];
grad.coilpos = RT;
grad.coilori= EZ;
for i=1:nchan 
  grad.label{i} = sprintf('chan%03d', i);
end

%% plot grad to check
figure(1)
title('Spherical Helmet')
ft_plot_sens(grad);
grid on
rotate3d
view(135, 20);

%% specify cfg
% create a spherical volume conductor with 10cm radius
vol.r = 10;
vol.o = [0 0 0];

cfg = [];
cfg.headmodel = vol;
cfg.grad = grad;
%cfg.nchan=306;
cfg.magscale=100;
cfg.dip.pos = [1 0 0];    % cm in x axis
cfg.dip.mom = [0 0 1];   % dopole should not point along x axis, this is along z axis
cfg.relnoise = 0; 
cfg.ntrials = 1; %trials were relevant for beamforming but not for SSS calculations
data = ft_dipolesimulation(cfg);
%data_5cm_v3 = ft_dipolesimulation(cfg);

[x,y,z,sphcord]=gensph_equidist(radius,5);
figure(3)
title('Spherical Helmet Equidistant')
plot3(x,y,z,'o');
rotate3d
view(135, 20);