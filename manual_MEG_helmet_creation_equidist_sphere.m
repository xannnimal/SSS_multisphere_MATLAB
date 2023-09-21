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
[xsens,ysens,zsens]=gensph_equidist(radius,5); % spacing^2 # sensors 
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
title('Spherical Helmet Equidistant')
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

%% construct phi from data for SSS calcs
magscale = 100;
% phi is data.trial{1,1} from collumn 1 to 100
phi= data.trial{1,1}(:,:);

%rescale phi data, every third row of data needs to be multiplied by 100
%corresponding to magnetometer channels
phi_0=[];
for i=(1:nchan)
    if mod(i,3)==0
        for j=(1:250)
            phi_0(i,j)=phi(i,j)*magscale;
        end
    else
        for j=(1:250)
            phi_0(i,j)=phi(i,j);
        end
    end
end

%% calculate SSS basis using "sss_script_for_Xan.m"
Lin = 8; % Truncation order of the internal VSH basis
Lout = 3; % Truncation order of the external VSH basis
vsh_origin = [0;0;0]; % Typical origin in the device coordinate system
dim_in = (Lin+1)^2 - 1; % Dimension of the internal SSS basis, should be 80

[Sin,SNin] = Sin_vsh_vv(vsh_origin,R,EX',EY',EZ',ch_types,Lin);
[Sout,SNout] = Sout_vsh_vv(vsh_origin,R,EX',EY',EZ',ch_types,Lout);

pS=pinv([SNin SNout]);
XN=pS*phi_0;
%reconstrct internal phi using SNin
data_rec=real(SNin*XN(1:dim_in,:));
angle_phi_0= subspace(phi_0,data_rec)*180/pi; 

%% plot data to check
chan_num=3; %3 corresponds to MEG0111
data_time=data.time{1,1};
data_chan_num=data.trial{1,1}(chan_num,:); 

figure(8);
hold on;
plot(data_time(:,1:100), data_chan_num(:,1:100))
plot(data_time(:,1:100),data_rec(chan_num,1:100))
title('MEG0121 Data before and after reconstruction, dipole 1cm x, spherical')
xlabel('time')
ylabel('MEG0121')
%ylim([-8e-12 8e-12])
legend({'Raw Data','Reconstructed'},'location','northwest')
hold off

