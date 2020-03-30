function [flag] = Figure5A_5C_5E_CMK()
flag = 1;
%% Parameter Relationships
load('Params_kinMOD_g1vskonTFkoff_plane_2.mat');
paramset2 = paramset;
load('Params_kinMOD_g1vskonTFkoff_plane.mat');
paramset1 = paramset;

param10000 = [paramset2,paramset1];

%% DON'T CHANGE STUFF -- this is the script to generate PARTS OF FIGURE 5
%clear all
% load trajectories - first set
load('Outputs_g1vskonTFkoff_plane_2.mat','allsGen','allsGen_mrna')
X3 = allsGen(:,1:12);
Z3 = allsGen_mrna(:,13);
clear allsGen; clear allsGen_mrna;
load('Outputs_g1vskonTFkoff_plane.mat','allsGen','allsGen_mrna')
X2 = allsGen(:,1:12);
Z2 = allsGen_mrna(:,13);
clear allsGen; clear allsGen_mrna;

X_all=[X3;X2];
Z_all=[Z3;Z2];

%%
allsGen=[X_all];
allsGen_mrna = [Z_all];

% model end point
val = 60; % 5hrs  
mod2 = {allsGen(:,1),allsGen(:,2),allsGen(:,3),allsGen(:,4),allsGen(:,5)};
mod40 = {allsGen(:,1),allsGen(:,6),allsGen(:,7),allsGen(:,8),allsGen(:,9)};
% model transferfunction
mod_transffn= {allsGen(:,1),allsGen(:,11),allsGen(:,10),allsGen(:,12)};
% model dynamics
mod_dynam = allsGen(:,10); % includes dynamON and dynamOFF
% model mRNA dynamics
mod_mrna_dynam = allsGen_mrna; 
% dynamics ON - (1:120)
% dynamics OFF - (240:360)
mod_dynam{1,1}(1:120);
mod_dynam{1,1}(240:360);

%%
% load data (end point)
load('Data_ypsgypcmk_2mTPvsdurTrajs_adhcrz19Azl_ENDPT','GE2m_cmk','GEdur_cmk','sGE2m_cmk','sGEdur_cmk')
GE2m = GE2m_cmk'; GEdur = GEdur_cmk'; sGE2m = sGE2m_cmk'; sGEdur = sGEdur_cmk';
% define Xaxis
load('Input_idealizedLight.mat','lightInputs')
kfact = 1000;
orig_lightInputs = lightInputs;
m2_xaxis=[orig_lightInputs(1).auc,orig_lightInputs(2).auc./kfact,orig_lightInputs(3).auc./kfact,orig_lightInputs(4).auc./kfact,orig_lightInputs(5).auc./kfact];
single_xaxis=[orig_lightInputs(1).auc,orig_lightInputs(6).auc./kfact,orig_lightInputs(7).auc./kfact,orig_lightInputs(8).auc./kfact,orig_lightInputs(9).auc./kfact];

%% protein Dynamics -- load data
load('Data_ROBOT_meanstd','YPS1CMK2_key','YPS1CMK2_mean','YPS1CMK2_std','GYP7_key','GYP7_mean','GYP7_std');
pYPS1_y = YPS1CMK2_mean(2,10:end); % take out the photobleaching initial part
pYPS1_sy = YPS1CMK2_std(2,10:end);
pYPS1_x = YPS1CMK2_mean(1,1:end-9).*60; pYPS1_x(1) = 1;

pCMK2_y = YPS1CMK2_mean(4,10:end); 
pCMK2_sy = YPS1CMK2_std(4,10:end);
pCMK2_x = YPS1CMK2_mean(1,1:end-9).*60; pCMK2_x(1) = 1;

% normalize from 0-1 irrespective
mm = 14;
norm01_pYPS = (pYPS1_y(3:mm)-min(pYPS1_y(3:mm)))./(max(pYPS1_y(3:mm))-min(pYPS1_y(3:mm)));
norm01_pCMK = (pCMK2_y(3:mm)-min(pCMK2_y(3:mm)))./(max(pCMK2_y(3:mm))-min(pCMK2_y(3:mm)));

% figure(4); errorbar(pYPS1_x(3:mm),norm01_pYPS, pYPS1_sy(3:mm), 'm.','markersize',20); 
% hold on; errorbar(pCMK2_x(3:mm),norm01_pCMK, pCMK2_sy(3:mm), 'b.','markersize',20); legend('pYPS1','pCMK2');
% xlabel('time(min)'); ylabel('norm FITC/SSC'); title('Robot Time course'); axis tight; box off;
% xlim([0 120])

%%
% scale all trajs to experimental values
% best fit line through the data FOR ENDPOINT and DOSE RESPONSE data
% 2min
p = polyfit(m2_xaxis,GE2m,1);
x2 = polyval(p,m2_xaxis);
% figure(10); errorbar(m2_xaxis,GE2m,sGE2m,'.r','markersize',20,'linewidth',2);...
%     hold on; plot(m2_xaxis,x2,'r'); xlabel('cum nuc occpuancy'); ylabel('FITC/SSC');
% 40min
p = polyfit(single_xaxis,GEdur,1);
xdur = polyval(p,single_xaxis);
% errorbar(single_xaxis,GEdur,sGEdur,'.b','markersize',20,'linewidth',2);...
%     plot(single_xaxis,xdur,'b'); axis tight; title('best fit line to data'); box off; axis tight;

% ENDPOINT: scale to max and min of best fit line
a = min(xdur); b = max(xdur);
scaled_2_endpt = [];
scaled_C_endpt = [];
mod_trans = [];
for i = 1:length(mod2{1}) 
    % find error for endpoint
    x = [mod2{1}{i}(end-val),mod2{2}{i}(end-val),mod2{3}{i}(end-val),mod2{4}{i}(end-val),mod2{5}{i}(end-val)];
    y = [mod40{1}{i}(end-val),mod40{2}{i}(end-val),mod40{3}{i}(end-val),mod40{4}{i}(end-val),mod40{5}{i}(end-val)];
    ix = (b-a).*(x-min(y))./(max(y)-min(y))+a;
    iy = (b-a).*(y-min(y))./(max(y)-min(y))+a;
    scaled_2_endpt(i,:) = ix;
    scaled_C_endpt(i,:) = iy; 
end

% plot out to check (just as expected, scaling this way b0 needs to be set,
% also need to exactly match traj & dose response) - could see if you
% actually get any fits (scale to max worse than scale to max and min)
modTemp2 = scaled_2_endpt(1,:);
modTempS = scaled_C_endpt(1,:); 
% figure(17);
% plot(m2_xaxis,modTemp2,'r'); hold on; 
% plot(single_xaxis,modTempS,'b');
% errorbar(m2_xaxis,GE2m,sGE2m,'.r','markersize',20,'linewidth',2); hold on;
% errorbar(single_xaxis,GEdur,sGEdur,'.b','markersize',20,'linewidth',2);
% axis tight; box off
% xlabel('TF AUC'); ylabel('FITC/SSC');

%% scale the model outputs to match the experiment (interpolate also)
% PROTEIN DYNAMICS
% 
mod_dynam = allsGen(:,10); % includes dynamON and dynamOFF
mod_dynamON = []; 
imod_dynamON = [];
imod_dynam_x = [];
imod_dynam_YPS = [];
imod_dynam_CMK = [];
for ii = 1:length(mod_dynam);

    mod_dynamON(ii,:) = mod_dynam{ii,1}(1:120); % dynamON

    % 1. align model with experiment (note that the X axis is the same for
    % pYPS1 and pCMK2

    imod_dynam_x = [pCMK2_x(1:12)];
    %imod_dynamON(ii,:) = interp1([1:120],mod_dynamON(ii,1:120),imod_dynam_x); % make pYPS1 and pCMK2 the same -1
    imod_dynamON(ii,:) = interp1([1:100],mod_dynamON(ii,1:100),imod_dynam_x(1:12)); % make pYPS1 and pCMK2 the same -1
    % plot to make sure interpolation worked
    
    % 3. scale the model to the data (pYPS and pCMK separately) - NOT
    % FOLDCHANGE
    a = min(norm01_pCMK); b = max(norm01_pCMK);
    imod_dynam_CMK(ii,:) = (b-a).*(imod_dynamON(ii,:)-min(imod_dynamON(ii,:)))./(max(imod_dynamON(ii,:))-min(imod_dynamON(ii,:))) + a;
    %figure(1); plot(imod_dynam_x,imod_dynam_YPS(ii,:),'b');

end

%% Sequential Fitting - how to carry through the ordering of variables
% initial set of parameters and trajectories in ORDER1 (b4 any least sq)
% - param10000 
% - scaled_2_endpt (endpt)
% - scaled_C_endpt (endpt)
% = mod_trans (dose response)
% - imod_dynam_YPS (protein)

%% DYNAMIC PROTEIN DATA (ACTUAL NEW FIT INDIVIDUALLY)
% this is for pYPS1, should do for pCMK2
uppYPS = norm01_pCMK(1:12) + pCMK2_sy(3:14);
lowYPS = norm01_pCMK(1:12) - pCMK2_sy(3:14);

numThru = []; %err_transf = [];
for i = 1:length(imod_dynam_CMK)
    ix = imod_dynam_CMK(i,1:12);

    uppx = ix <= uppYPS;
    lowx = ix >= lowYPS;
    
    numThru(i) = sum(uppx+lowx);

end

indxFits_prot = find(numThru>=20); % pYPS1 = 17 % if protein only pCMK2 = 20; 19 otherwise

% plot to overlay with data, see goodness of fit
figure(50); 
%meanPROT = mean(imod_dynam_CMK(indxFits_prot,:));
%stdPROT = std(imod_dynam_CMK(indxFits_prot,:));
%shadedErrorBar(pYPS1_x(1:12)+15,meanPROT(1:12),stdPROT(1:12),'c'); hold on;
plot(pYPS1_x(1:12)+15,imod_dynam_CMK(indxFits_prot,1:12),'c','markersize',20); hold on;
%plot(pYPS1_x(1:12)+15,imod_dynam_CMK(indx_prot(paraNum),1:12),'b','markersize',20); hold on;
errorbar(pCMK2_x(3:14),norm01_pCMK(1:12),pCMK2_sy(3:14),'k.','markersize',20,'linewidth',2);
xlabel('time(min)'); ylabel('norm FITC/SSC'); 
%title(['pCMK2 Protein Dynamics ',num2str(numel(paraNum)),' params']); box off; axis tight; 
title(['pCMK2 Protein Dynamics ',num2str(numel(indxFits_prot)),' params']); box off; axis tight; 


param10000_P = param10000(:,indxFits_prot);
scaled_2_endpt_P = scaled_2_endpt(indxFits_prot,:);
scaled_C_endpt_P = scaled_C_endpt(indxFits_prot,:);
Xall_P = X_all(indxFits_prot,:);
imod_dynam_CMK_P = imod_dynam_CMK(indxFits_prot,1:12);

%% Endpoint
% compute straightness of traj
strTrajScore = [];
straightX = []; straightY = [];
for i = 1:length(scaled_2_endpt_P)
    ix = scaled_2_endpt_P(i,:);
    iy = scaled_C_endpt_P(i,:);
    
    % find residual, how close to straight line
    p = polyfit(m2_xaxis,ix,1); % 2min
    x2 = polyval(p,m2_xaxis);
    straightX(i) = sum((x2-ix).^2);
    p = polyfit(single_xaxis,iy,1); % constant
    y2 = polyval(p,single_xaxis);
    straightY(i) = sum((y2-iy).^2);
    strTrajScore(i) = straightY(i) + straightX(i);
end

[A,indx_strTraj]=sort(strTrajScore);
straightNum = 750; %600%800;%1312;
f_scaled_2_endpt_P = scaled_2_endpt_P(indx_strTraj(1:straightNum),:);
f_scaled_C_endpt_P = scaled_C_endpt_P(indx_strTraj(1:straightNum),:);
f_param10000_P = param10000_P(:,indx_strTraj(1:straightNum));
f_imod_dynam_CMK_P = imod_dynam_CMK_P(indx_strTraj(1:straightNum),1:12);
f_Xall_P = Xall_P(indx_strTraj(1:straightNum),:);
%%
% compute least square error in transfer function + endpt
upp_red = GE2m + sGE2m;
low_red = GE2m - sGE2m;
upp_blue = GEdur + sGEdur;
low_blue = GEdur - sGEdur;

numThru = [];
for i = 1:length(f_scaled_2_endpt_P)
    ix = f_scaled_2_endpt_P(i,:);
    iy = f_scaled_C_endpt_P(i,:);
    
    uppx = ix <= upp_red; lowx = ix >= low_red;
    uppy = iy <= upp_blue; lowy = iy >= low_blue;
    
    numThru(i) = sum(uppx+lowx) + sum(uppy+lowy);
end

indxFits_endpt = find(numThru >= 13); % for pYPS1 = 16; pCMK2 is 13
ff_scaled_2_endpt_P = f_scaled_2_endpt_P(indxFits_endpt,:);
ff_scaled_C_endpt_P = f_scaled_C_endpt_P(indxFits_endpt,:);
ff_param10000_P = f_param10000_P(:,indxFits_endpt);
ff_imod_dynam_CMK_P = f_imod_dynam_CMK_P(indxFits_endpt,1:12);
ff_Xall_P = f_Xall_P(indxFits_endpt,:);


figure(19);
%paraNum = length(indx_endpt);
mm2 = mean(f_scaled_2_endpt_P(indxFits_endpt,:));
sm2 = std(f_scaled_2_endpt_P(indxFits_endpt,:));
mSing = mean(f_scaled_C_endpt_P(indxFits_endpt,:));
sSing = std(f_scaled_C_endpt_P(indxFits_endpt,:));
%shadedErrorBar(m2_xaxis, mm2,sm2,'r'); hold on;
%shadedErrorBar(single_xaxis, mSing,sSing,'b'); 
plot(m2_xaxis, f_scaled_2_endpt_P(indxFits_endpt,:),'r'); hold on;
plot(single_xaxis, f_scaled_C_endpt_P(indxFits_endpt,:),'b'); 
%plot(m2_xaxis, scaled_2_endpt_P(indx_endpt(1:paraNum),:),'r'); hold on;
%plot(single_xaxis, scaled_C_endpt_P(indx_endpt(1:paraNum),:),'b'); 
errorbar(m2_xaxis,GE2m,sGE2m,'.k','markersize',20,'linewidth',2); hold on;
errorbar(single_xaxis,GEdur,sGEdur,'.k','markersize',20,'linewidth',2);
axis tight; box off; axis on;
title(['pYPS1 endpt ', num2str(numel(indxFits_endpt)),' params']); xlabel('nuclear occupancy'); ylabel('FITC/SSC');
%title(['pCMK2 endpt ', num2str(numel(1:paraNum)),' params']); xlabel('nuclear occupancy'); ylabel('FITC/SSC');


param10000_PE = f_param10000_P(:,indxFits_endpt);
scaled_2_endpt_PE = f_scaled_2_endpt_P(indxFits_endpt,:);
scaled_C_endpt_PE = f_scaled_C_endpt_P(indxFits_endpt,:);
imod_dynam_CMK_PE = f_imod_dynam_CMK_P(indxFits_endpt,1:12);
Xall_PE = f_Xall_P(indxFits_endpt,:);


%% DYNAMIC - VALIDATE ON 2min pulsed dynamic data
% PLOT the VALIDATION - 2min6min with 2hours!!!

% plotting
i_f = 4; %4;
i_e = 31; %31;
j_f = 12;
j_e = 23;
a = YPS1CMK2_mean(5,i_f:end-i_e);
aa = YPS1CMK2_mean(5,i_f:end-i_e)./a(1);
b = YPS1CMK2_mean(4,j_f:end-j_e);
bb = YPS1CMK2_mean(4,j_f:end-j_e)./b(1);
cmk_2mDATA = (aa - min(aa))./(max(bb) - min(bb));
cmk_consDATA = (bb - min(bb))./(max(bb) - min(bb));

cmk_2m = [];
cmk_cons = [];
for i = 1:length(Xall_PE(:,5))
    a = Xall_PE{i,5}(1:100); % 120
    b = Xall_PE{i,9}(61:160); % 180
    uppB = max(b);
    lowB = min(b);
    cmk_2m(i,:) = (a - lowB)./(uppB - lowB);
    cmk_cons(i,:) = (b - lowB)./(uppB - lowB);
end

figure(1); 
mean2m = mean(cmk_2m); std2m = std(cmk_2m);
meanCons = mean(cmk_cons); stdCons = std(cmk_cons);
plot(linspace(24,112,100),cmk_2m','r'); hold on;
plot(linspace(24,112,100),cmk_cons','b'); 
errorbar(YPS1CMK2_mean(1,i_f:end-i_e).*60 , cmk_2mDATA, YPS1CMK2_std(5,i_f:end-i_e),'k.','markersize',20,'linewidth',2); hold on;
errorbar(YPS1CMK2_mean(1,j_f:end-j_e).*60 - 64, cmk_consDATA , YPS1CMK2_std(4,j_f:end-j_e), 'k.','markersize',20,'linewidth',2);

xlabel('time (min)'); ylabel('norm Protein');
axis tight; box off; title('pCMK2 model 2m6m vs continous 321');
