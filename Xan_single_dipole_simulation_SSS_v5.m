%%%%%%%%%%%%%%%%%%
%code from field trip
%"https://www.fieldtriptoolbox.org/example/compute_forward_simulated_data_and_apply_a_beamformer_scan/"
% modified to do single dipole simulation 
%%%%%%%%%%%%%%%%%%%
clear
import mne.*
nchan=306;
mags = 3:3:nchan; % Magnetometer channels
magscale = 100; % Numerical scaling factor between magnetometer and gradiometer signals
coordsys = 'device'; % Let's use the device coordinate system (instead of the head coordinate system) in our calculations for now
rawfile = "sample_audvis_raw.fif";
grad_neuromag306 = ft_read_sens('Program Files/MATLAB/fieldtrip-20230613/template/gradiometer/neuromag306.mat','senstype', 'meg');

%% specify grad
% read in "neuromag306.mat" for MEG sensors
%grad = ft_read_sens('Program Files/MATLAB/fieldtrip-20230613/template/gradiometer/neuromag306.mat','senstype', 'meg');
%create grad using sensor positions from fiff file
[R,EX,EY,EZ] = fiff_getpos(rawfile,coordsys);
ch_types = ones(size(EX',1),1); % 1 = magnetometers, 0 = gradiometers for S_in code
%fieldTrip needs 306x3, not 3x306 to work so need to transpose
RT=transpose(R);
EXT=transpose(EX);
EYT=transpose(EY);
EZT=transpose(EZ);
%R=sensor positions, EZ=normal orientations for grad.coilori
grad = [];
grad.coilpos = RT;
grad.coilori= EZT;
grad.senstype = 'meg';
grad.tra= eye(size(RT,1));
%grad.tra= grad_neuromag306.tra;
for i=1:nchan
  grad.label{i} = sprintf('MEG%03d', i);
end

%% plot grad to check helmet/sensor layout
% figure(1)
% ft_plot_sens(grad, 'label', 'yes');
% grid on
% rotate3d
% view(135, 20);

%% specify cfg
% create a spherical volume conductor with 10cm radius
vol.r = 10;
vol.o = [0 0 0];

cfg = [];
cfg.headmodel = vol;
cfg.channel = {'MEG'};
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
% phi is data.trial{1,1} from collumn 1 to 100
phi= data.trial{1,1}(:,:);
% phi=[];
% for i=(1:time_size)
%     phi(:,i)=data_0.trial{1,1}(:,i);
% end

%rescale phi data, every third row of data needs to be multiplied by 100
%corresponding to magnetometer channels
phi_0=[];
for i=(1:nchan)
    if mod(i,3)==0
        for j=(1:250)
            phi_0(i,j)=phi(i,j);
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
%convert simulated data to fiff file
%%%%%%%%% ISSUE: grad.chanunit, some are T/m, not T, and the function
%%%%%%%%% "fieltrip2fiff" can't work if its T/m, issue is currently being
%%%%%%%%% worked on
% fiff_file  = 'ctf-raw.fif';
% fieldtrip2fiff(fiff_file, data)
% rawfile = "ctf-raw.fif";
% [data2,times] = mne_ex_read_raw(rawfile,1,60);
% phi = data2(:,100);


% [Sin,SNin] = Sin_basic(rawfile,vsh_origin,magscale,coordsys,Lin);
% [Sout,SNout] = Sout_basic(rawfile,vsh_origin,magscale,coordsys,Lout); %this originally said "Lin"
%mags = mags_phi(rawfile,vsh_origin,magscale,coordsys,Lin);

[Sin,SNin] = Sin_vsh_vv(vsh_origin,R,EX,EY,EZ,ch_types,Lin);
[Sout,SNout] = Sout_vsh_vv(vsh_origin,R,EX,EY,EZ,ch_types,Lout);
pS=pinv([SNin SNout]);
XN=pS*phi_0;
%reconstrct internal phi using SNin
data_rec=real(SNin*XN(1:dim_in,:));
angle_phi_0= subspace(phi_0,data_rec)*180/pi; 


%% plot data to check
%plot data from single channel
chan_num=6; 
data_time=data.time{1,1};
data_chan_num=data.trial{1,1}(chan_num,:); 
angle=subspace(SNin,data.trial{1,1}(:,5))*180/pi;

figure(2);
hold on;
plot(data_time(:,1:100), data_chan_num(:,1:100))
plot(data_time(:,1:100),data_rec(chan_num,1:100))
title('MEG0121 Data before and after reconstruction, dipole 1cm x')
xlabel('time')
ylabel('MEG0121')
%ylim([-8e-12 8e-12])
legend({'Raw Data','Reconstructed'},'location','northwest')
hold off

% plot data from all channels
figure(3);
plot(data_time, data.trial{1,1})
title('Raw Data dipole 1cm')

figure(4);
plot(data_time, data_rec)
title('Reconstructed Data dipole 1cm')