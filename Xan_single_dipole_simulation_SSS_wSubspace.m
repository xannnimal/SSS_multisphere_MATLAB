%%%%%%%%%%%%%%%%%%
%code from field trip
%"https://www.fieldtriptoolbox.org/example/compute_forward_simulated_data_and_apply_a_beamformer_scan/"
% modified to create a manual config 306 magnetometer helmet in the same config as an MEG
% helmet with no gradiometers. Code also configures a single dipole
% the helmet is rotated around the z axis and the SSS interior basis is
% calculated (using manual config) for each rotation. Then, subspace angles between the
% un-rotated and each rotated Sin basis are calculated and shown in a
% surface plot. 
%%%%%%%%%%%%%%%%%%%
clear
import mne.*
nchan=306;
mags = 3:3:nchan; % Magnetometer channels
magscale = 100; % Numerical scaling factor between magnetometer and gradiometer signals
coordsys = 'device'; % Let's use the device coordinate system (instead of the head coordinate system) in our calculations for now
rawfile = "sample_audvis_raw.fif";
Lin = 8; % Truncation order of the internal VSH basis
Lout = 3; % Truncation order of the external VSH basis
vsh_origin = [0;0;0]; % Typical origin in the device coordinate system
dim_in = (Lin+1)^2 - 1; % Dimension of the internal SSS basis, should be 80

%% specify grad
%create grad using sensor positions from fiff file
[R,EX,EY,EZ] = fiff_getpos(rawfile,coordsys);
ch_types = ones(size(EX',1),1); % 1 = magnetometers, 0 = gradiometers for S_in code
%fieldTrip needs 306x3, not 3x306 to work so need to transpose
RT=transpose(R);
EXT=transpose(EX);
EYT=transpose(EY);
EZT=transpose(EZ);

grad = [];
grad.coilori= EZT;
grad.senstype = 'meg';
grad.tra= eye(size(RT,1));
for i=1:nchan
  grad.label{i} = sprintf('MEG%03d', i);
end

%% specify cfg
% create a spherical volume conductor with 10cm radius
%vol.r = 10;
%vol.o = [0 0 0];

cfg = [];
%cfg.headmodel = vol;
cfg.channel = {'MEG'};
cfg.magscale=100;
cfg.relnoise = 0; 
cfg.ntrials = 1; %trials were relevant for beam
%cfg.dip.pos = [1 0 0];
%cfg.dip.mom = [0 0 1];   % dopole should not point along x axis, this is along z axis

% angle = 150;
% theta = deg2rad(angle);
% c = cos(theta);
% s = sin(theta);
% rot = [c, -s, 0; s, c, 0; 0, 0, 1];
% rot_RT=zeros([size(RT,1) size(RT,2)]);
% for i=1:size(RT,1)
%     rot_rt = rot*RT(i,:)';
%     rot_RT(i,:)=rot_rt';
% end
% grad.coilpos = rot_RT;
% cfg.grad = grad;
% 
% %% simulate data from dipole
% data = ft_dipolesimulation(cfg);
% phi_0= data.trial{1,1}(:,:);
% 
% %% calculate SSS basis using "sss_script_for_Xan.m"
% 
% [Sin_150check,SNin] = Sin_vsh_vv(vsh_origin,rot_RT',EX,EY,EZ,ch_types,Lin);
% filename_Sin=sprintf("Sin_%02d.mat",angle);
% %     %save(filename_Sin,"Sin","-mat")
% %     %[Sout,SNout] = Sout_vsh_vv(vsh_origin,R,EX,EY,EZ,ch_types,Lout);
% %     %pS=pinv([SNin SNout]);
% %     %XN=pS*phi_0;
% %     %reconstrct internal phi using SNin
% %     %data_rec=real(SNin*XN(1:dim_in,:));

%% put this in a for loop
all_Sin=cell(72,1);
for n=0:72 %360/5=72, angle incriments in 5 deg from 0 to 360
    clear("Sin")
    % rotate the coilpos
    angle = n*5;
    theta = deg2rad(angle);
    c = cos(theta);
    s = sin(theta);
    rot = [c, -s, 0; s, c, 0; 0, 0, 1];
    rot_RT=zeros([size(RT,1) size(RT,2)]);
    for i=1:size(RT,1)
        rot_rt = rot*RT(i,:)';
        rot_RT(i,:)=rot_rt';
    end
    grad.coilpos = rot_RT;
    cfg.grad = grad;

%% simulate data from dipole
    % data = ft_dipolesimulation(cfg);
    % phi_0= data.trial{1,1}(:,:);

%% calculate SSS basis using "sss_script_for_Xan.m"
    [Sin,SNin] = Sin_vsh_vv(vsh_origin,rot_RT',EX,EY,EZ,ch_types,Lin);
    all_Sin{n+1}=Sin;
    n
    %filename_Sin=sprintf("Sin_%02d.mat",angle);
    %save(filename_Sin,"Sin","-mat")
    %[Sout,SNout] = Sout_vsh_vv(vsh_origin,R,EX,EY,EZ,ch_types,Lout);
    %pS=pinv([SNin SNout]);
    %XN=pS*phi_0;
    %reconstrct internal phi using SNin
    %data_rec=real(SNin*XN(1:dim_in,:));
end

all_angles=cell(73,1);
for i=1:73
    for j=1:80
        thetas = subspace(all_Sin{1,1},all_Sin{i,1}(:,j));
        all_angles{i}=[all_angles{i},thetas/pi*180];
    end
end

all_angles_mat=[];
for i=1:73
    all_angles_mat(i,:) = all_angles{i};
end

%% plot
all_col=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80];
all_angles=[0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200,205,210,215,220,225,230,235,240,245,250,255,260,265,270,275,280,285,290,295,300,305,310,315,320,325,330,335,340,345,350,355,360];
figure(1);
surf(all_angles,all_col,all_angles_mat')
xticks([0 45 90 135 180 225 270 315 360])
title('Subspace Angles, Rotations Around the Z Axis, 306 Magnetometer Helmet')
xlabel('Rotation Angle')

% chan_num=6; 
% data_time=data.time{1,1};
% data_chan_num=data.trial{1,1}(chan_num,:); 
% angle=subspace(SNin,data.trial{1,1}(:,5))*180/pi;
% 
% figure(2);
% hold on;
% plot(data_time(:,1:100), data_chan_num(:,1:100))
% plot(data_time(:,1:100),data_rec(chan_num,1:100))
% title('MEG0121 Data before and after reconstruction, dipole 1cm x')
% xlabel('time')
% ylabel('MEG0121')
% %ylim([-8e-12 8e-12])
% legend({'Raw Data','Reconstructed'},'location','northwest')
% hold off
% 
% % plot data from all channels
% figure(3);
% plot(data_time, data.trial{1,1})
% title('Raw Data dipole 1cm')
% 
% figure(4);
% plot(data_time, data_rec)
% title('Reconstructed Data dipole 1cm')