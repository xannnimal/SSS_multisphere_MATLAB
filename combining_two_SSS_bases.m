%% combining multiple SSS bases into one
% use projection from tSSS but in the spatial dimenson to find a common
% basis between the two spheres (bulk), then create a basis containing the
% bulk, parts of S1 that aren't included in bulk, and parts of S2 that
% aren't included in bulk
%%%%%%%%%%%%%%%%%
clear
import mne.*

nchan=306;
mags = 3:3:nchan; % Magnetometer channels
magscale = 100; % Numerical scaling factor between magnetometer and gradiometer signals
coordsys = 'device'; % Let's use the device coordinate system (instead of the head coordinate system) in our calculations for now
rawfile = "sample_audvis_raw.fif";
Lin = 8; % Truncation order of the internal VSH basis
Lout = 3; % Truncation order of the external VSH basis
dim_in = (Lin+1)^2 - 1; % Dimension of the internal SSS basis, should be 80
corr_limit = 0.9; %usually 0.9 to 0.98

%% specify grad
%create grad using sensor positions from fiff file
[R,EX,EY,EZ] = fiff_getpos(rawfile,coordsys);
ch_types = ones(size(EX',1),1); % 1 = magnetometers, 0 = gradiometers for S_in code
%fieldTrip needs 306x3, not 3x306 to work so need to transpose
RT=transpose(R);
EXT=transpose(EX);
EYT=transpose(EY);
EZT=transpose(EZ);

% grad = [];
% grad.coilori= EZT;
% grad.coilpos= RT;
% grad.senstype = 'meg';
% grad.tra= eye(size(RT,1));
% for i=1:nchan
%   grad.label{i} = sprintf('MEG%03d', i);
% end
% specify cfg
% create a spherical volume conductor with 10cm radius
% vol.r = 10;
% vol.o = [0 0 0];
% 
% cfg = [];
% cfg.headmodel = vol;
% cfg.channel = {'MEG'};
% cfg.magscale=100;
% cfg.relnoise = 0; 
% %cfg.ntrials = 1; %trials were relevant for beam
% cfg.dip.pos = [1 0 0];
% cfg.dip.mom = [0 0 1]; 
% cfg.grad = grad;

%% simulate data from dipole
% data = ft_dipolesimulation(cfg);
% phi= data.trial{1,1}(:,:);

%% calculate SSS basis based on grad
for i=1:51 %0cm to 5cm is 0mm to 50mm
    l=i-1;
    x2=(l)/10;
    for j=1:51
        y2=(j-1)/10;
        vsh_origin_1 = [0;0;0];% Typical origin in the device coordinate system
        vsh_origin_2 = [x2;y2;0];%cm or m??
        coord_sphere_2x(i,j)=x2;
        coord_sphere_2y(i,j)=y2;

        [Sin_1,SNin_1] = Sin_vsh_vv(vsh_origin_1,R,EX,EY,EZ,ch_types,Lin);
        [Sout_1,SNout_1] = Sout_vsh_vv(vsh_origin_1,R,EX,EY,EZ,ch_types,Lout);
        [Sin_2,SNin_2] = Sin_vsh_vv(vsh_origin_2,R,EX,EY,EZ,ch_types,Lin);
        [Sout_2,SNout_2] = Sout_vsh_vv(vsh_origin_2,R,EX,EY,EZ,ch_types,Lout);
        %pSN_1=pinv([SNin_1 SNout_1]);
        %pSN_2=pinv([SNin_2 SNout_2]);

        Ein_1 = orth(SNin_1); %orth returns an orthonormal basis for range of Bin
        Ein_2 = orth(SNin_2);
        [Q1,ignore] = qr(Ein_1,0); % upper triangular matrix 'ignore' of the same dimension as Ein and a unitary matrix QA so that X = Q*R.
        [Q2,ignore] = qr(Ein_2,0);
        [U,S,V] = svd(Q1'*Q2);
        Up = Q1*U';
        s = diag(S); %cosines, principle angles
        s=s';
        num=0;
        len=length(s);
        for k=1:len
            if s(k)>0.999 %some cosines show as 1, others as 1.000 for some reason
                num=num+1;
            end
            numofones(i,j)=num;
        end  
        j
    end
    i
end

outputpath1 = 'C:\Users\xanmc\OneDrive\Documents\UofW\RESEARCH\Sin origin shift data';
save(fullfile(outputpath1,...
    sprintf('X_coord_Sin2',coord_sphere_2x)),'coord_sphere_2x')
outputpath2 = 'C:\Users\xanmc\OneDrive\Documents\UofW\RESEARCH\Sin origin shift data';
save(fullfile(outputpath2,...
    sprintf('Y_coord_Sin2',coord_sphere_2y)),'coord_sphere_2y')
outputpath3 = 'C:\Users\xanmc\OneDrive\Documents\UofW\RESEARCH\Sin origin shift data';
save(fullfile(outputpath3,...
    sprintf('Num_of_zeros',numofones)),'numofones')

%% code from "tsss_rawdata_ctf_test_22070.m" lines 55 to 77
%X = pSN*data(meg_chs,:);
%data_in = real(SNin*X(1:size(SNin,2),:));
%data_out = real(SNout*X(size(SNin,2)+1:end,:));
%data_res = data(meg_chs,:) - (data_in+data_out);
% X_1=pSN_1*phi; %could be norm or not
% X_2=pSN_2*phi;
% %reconstrct internal phi using SNin
% data_in_1 = real(SNin_1*X_1(1:dim_in,:));
% data_out_1 = real(SNout_1*X_1(size(SNin_1,2)+1:end,:));
% data_in_2 = real(SNin_2*X_2(1:dim_in,:));
% data_out_2 = real(SNout_2*X_2(size(SNin_2,2)+1:end,:));
% data_res = phi - (data_in_1+data_in_2);
%C=intersect(X_1,X_2,'stable','rows');
%C_data=intersect(data_in_1,data_in_2,'stable','rows');
% intersection of Phi_in and Phi_res not Cin and Cout as in paper
% Bin = data_in/norm(data_in);
% Bres = data_res/norm(data_res);
% Ein = orth(Bin'); %orth returns an orthonormal basis for range of Bin
% Eres = orth(Bres');
% [QA,ignore] = qr(Ein,0); % upper triangular matrix 'ignore' of the same dimension as Ein and a unitary matrix QA so that X = Q*R.
% [QB,ignore] = qr(Eres,0);
% [U,S,V] = svd(QA'*QB);
% Up = QA*U';
% Vp = QB*V;
% s = diag(S);
% inter_indices = find(s>corr_limit); 
% max(s), length(inter_indices) %%%
% Eproj = Vp(:,inter_indices);
% P = Eproj*Eproj';
% Bp = ((eye(size(P,1))-P)*data_in')';
% data_phi = Bp;