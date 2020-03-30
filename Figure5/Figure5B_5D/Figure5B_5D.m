function [flag] = Figure5B_5D()
flag = 1;
%% The point of this code is to plot gray overlays
%% THE INPUT PARAMETERS: CMK2
% set #2 all of them
load('Params_kinMOD_g1vskonTFkoff_plane_2.mat');
paramset2 = paramset;

load('Params_kinMOD_g1vskonTFkoff_plane.mat');
paramset1 = paramset;

%% Figure 5B CMK2
% set #2 pCMK2 - dynamic protein timecourse
% combined
param10000 = [paramset2,paramset1];

% PLOT ON THE g1 VS KON*T + KOFF PLANE: CMK2
T = 2.6;
kon = param10000(1,:); koff = param10000(2,:); g1 = param10000(4,:);
figure(11); plot(log10(kon.*T + koff),log10(g1), 'k.','markersize',10); hold on;
xlabel('log10(kon*T + koff)'); ylabel('log10(g1)'); title('pCMK2 fits to dynamic data');
xlim([-2 2]); ylim([-2 1]);

% OVERLAY THE OUTPUT PARAMETER FITS ASSOCIATED WITH THESE 3 INPUT PARAMETER
% FITS: CMK2
% from set #1 and #2
load('ParamFits_CMK_g1vskonTFkoff_P.mat','param10000_P');

% PLOT ON THE g1 VS KON*T + KOFF PLANE: CMK2
kon_P = param10000_P(1,:); koff_P = param10000_P(2,:); g1_P = param10000_P(4,:);
figure(11); plot(log10(kon_P.*T + koff_P),log10(g1_P), 'c.','markersize',10); % set #1 and #2 
xlabel('log10(kon*T + koff)'); ylabel('log10(g1)'); title('pCMK2 fits to dynamic protein data 807');
axis([-2 2 -2 1])

% GET THE LOWER BOUNDS FOR G1, KON, AND KOFF - for pCMK2
g1 = 0.1115;
kon = 0.0010;
koff = 0.0518; 
% CALC THE TS values
% 1/g1 + 1/kon*TF + koff = 1/Ts
TsInv = 1./param10000_P(4,:) + 1./(param10000_P(1,:).*2.6+param10000_P(2,:));
TsInv_mean = mean(TsInv); % 7.66
TsInv_max = max(TsInv); % 9
TsInv_min = min(TsInv); % 6.9
%TsInv_std = std(TsInv); % 0.386
Ts = 1./TsInv;
Ts_mean = mean(Ts); % 0.1309
Ts_std = std(Ts); % 0.0065

%% Figure 5B YPS1
% THE INPUT PARAMETERS: YPS1
% set #2 all of them
load('Params_kinMOD_g1vskonTFkoff_plane_2.mat');
paramset2 = paramset;

load('Params_kinMOD_g1vskonTFkoff_plane.mat');
paramset1 = paramset;

% combined
param10000 = [paramset2,paramset1];

% PLOT ON THE g1 VS KON*T + KOFF PLANE: YPS1
T = 2.6;
kon = param10000(1,:); koff = param10000(2,:); g1 = param10000(4,:);
figure(2); plot(log10(kon.*T + koff),log10(g1), 'k.','markersize',10); hold on;
xlabel('log10(kon*T + koff)'); ylabel('log10(g1)'); title('pYPS1 fits to dynamic data');
xlim([-2 2]); ylim([-2 1]);
% OVERLAY THE OUTPUT PARAMETER FITS ASSOCIATED WITH THESE 3 INPUT PARAMETER
% FITS: YPS1
% set #1
% OVERLAY THE OUTPUT PARAMETER FITS ASSOCIATED WITH THESE 3 INPUT PARAMETER
% FITS: CMK2
% from set #1 and #2
load('ParamFits_YPS_g1vskonTFkoff_P.mat','param10000_P');
% addl sets that all fit
T = 2.6;
kon_YPS = param10000_P(1,:); koff_YPS = param10000_P(2,:); g1_YPS = param10000_P(4,:);
%V_pYPS = 1./g1_YPS + 1./(kon_YPS.*T + koff_YPS);
%V = mean(V_pYPS);

% PLOT ON THE g1 VS KON*T + KOFF PLANE: YPS1
kon_P = param10000_P(1,:); koff_P = param10000_P(2,:); g1_P = param10000_P(4,:);
figure(2); plot(log10(kon_P.*T + koff_P),log10(g1_P), 'm.','markersize',10); % set #1 and #2 
xlabel('log10(kon*T + koff)'); ylabel('log10(g1)'); title('pYPS1 fits to dynamic protein data 2355');

% GET THE LOWER BOUNDS FOR G1, KON, AND KOFF - for pCMK2
g1 = 0.03;
kon = 0.001;
koff = 0.006; 
% CALC THE TS values
% 1/g1 + 1/kon*TF + koff = 1/Ts
TsInv = 1./param10000_P(4,:) + 1./(param10000_P(1,:).*2.6+param10000_P(2,:));
TsInv_mean = mean(TsInv); % 23.8
TsInv_max = max(TsInv); % 18
TsInv_min = min(TsInv); % 32.8


%% Figure 5D CMK2 O-O
% set #2 pCMK2 - Output-Occpuancy plots
% set #2 all of them
load('Params_kinMOD_g1vskonTFkoff_plane_2.mat');
paramset2 = paramset;
load('Params_kinMOD_g1vskonTFkoff_plane.mat');
paramset1 = paramset;

param10000 = [paramset2,paramset1];

% PLOT ON THE g1 VS KON*T + KOFF PLANE: CMK2
T = 2.6;
kon = param10000(1,:); koff = param10000(2,:); g1 = param10000(4,:);
figure(1111); plot(log10(kon.*T + koff),log10(g1), 'k.','markersize',10); hold on;
xlabel('log10(kon*T + koff)'); ylabel('log10(g1)'); title('pCMK2 fits to dynamic data');
xlim([-2 2]); ylim([-2 1]);


% from fit outputs
load('ParamFits_CMK_g1vskonTFkoff_PE.mat','param10000_PE');

% PLOT ON THE g1 VS KON*T + KOFF PLANE: CMK2
kon_PE = param10000_PE(1,:); koff_PE = param10000_PE(2,:); g1_PE = param10000_PE(4,:);
figure(1111); plot(log10(kon_P.*T + koff_P),log10(g1_P), 'c.','markersize',10); % set #1 and #2 
xlabel('log10(kon*T + koff)'); ylabel('log10(g1)'); title('pCMK2 fits to O-O data');
axis([-2 2 -2 1])

%% Figure 5D YPS1 O-O
% THE INPUT PARAMETERS: YPS1
% set #1
load('Params_kinMOD_g1vskonTFkoff_plane_2.mat');
paramset2 = paramset;
load('Params_kinMOD_g1vskonTFkoff_plane.mat');
paramset1 = paramset;

param10000 = [paramset2,paramset1];

% PLOT ON THE g1 VS KON*T + KOFF PLANE: YPS1
T = 2.6;
kon = param10000(1,:); koff = param10000(2,:); g1 = param10000(4,:);
figure(4); plot(log10(kon.*T + koff),log10(g1), 'k.','markersize',10); hold on;
xlabel('log10(kon*T + koff)'); ylabel('log10(g1)'); title('pYPS1 fits to dynamic data 0');
xlim([-2 2]); ylim([-2 1]);

% LOAD FITS
load('ParamFits_YPS_g1vskonTFkoff_PE.mat','param10000_PE');

T = 2.6;
kon_YPS = param10000_PE(1,:); koff_YPS = param10000_PE(2,:); g1_YPS = param10000_PE(4,:);


% PLOT ON THE g1 VS KON*T + KOFF PLANE: YPS1
kon_PE = param10000_PE(1,:); koff_PE = param10000_PE(2,:); g1_PE = param10000_PE(4,:);
figure(4); plot(log10(kon_PE.*T + koff_PE),log10(g1_PE), 'm.','markersize',10); % set #1 and #2 
xlabel('log10(kon*T + koff)'); ylabel('log10(g1)'); title('pYPS1 fits to Output-Occupancy data 300');

% GET THE LOWER BOUNDS FOR G1, KON, AND KOFF - for pCMK2
g1 = param10000_PE(4,:); %0.03;
kon = param10000_PE(1,:); %0.001;
koff = param10000_PE(2,:); %0.006; 
konTFkoff = kon.*2.6+koff;
max_konTFkoff = max(konTFkoff)
min_konTFkoff = min(konTFkoff)

end