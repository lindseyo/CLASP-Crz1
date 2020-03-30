%% Plot the Np/Nc vs slope ratio plot for pYPS1 or pCMK2 - switch between folders
function [flag] = FigureSM1C()
flag = 1;

%% pYPS
% load in parameters
load(['20200212_outputs_g1vskonTFkoff_plane_pYPS1_N.mat'],'allsGen','allsGen_N')
X2 = allsGen(:,1:12);
Z2 = allsGen_N(:,1:12);
clear allsGen; clear allsGen_N;
X_all=[X2];
Z_all=[Z2];
allsGen=[X_all];
allsGen_N = [Z_all];

%  model end point
val = 60; % 5hrs  % could try 180
mod2 = {allsGen(:,1),allsGen(:,2),allsGen(:,3),allsGen(:,4),allsGen(:,5)};
mod40 = {allsGen(:,1),allsGen(:,6),allsGen(:,7),allsGen(:,8),allsGen(:,9)};
   
% define Xaxis
load('Input_idealizedLight.mat','lightInputs')
kfact = 1000;
orig_lightInputs = lightInputs;
m2_xaxis=[orig_lightInputs(1).auc,orig_lightInputs(2).auc./kfact,orig_lightInputs(3).auc./kfact,orig_lightInputs(4).auc./kfact,orig_lightInputs(5).auc./kfact];
single_xaxis=[orig_lightInputs(1).auc,orig_lightInputs(6).auc./kfact,orig_lightInputs(7).auc./kfact,orig_lightInputs(8).auc./kfact,orig_lightInputs(9).auc./kfact];

%%%%%%%%%%%%%%%%%%
% slope ratio
slopeRatio = [];
for i = 1:length(mod2{1})
    % trajectory
    x = [mod2{1}{i}(end-val),mod2{2}{i}(end-val),mod2{3}{i}(end-val),mod2{4}{i}(end-val),mod2{5}{i}(end-val)];
    y = [mod40{1}{i}(end-val),mod40{2}{i}(end-val),mod40{3}{i}(end-val),mod40{4}{i}(end-val),mod40{5}{i}(end-val)];
    % interpolate between points
    y1 = interp1(single_xaxis, y, m2_xaxis);
    % find slope ratio
    slopeRatio(i) = (x(5) - x(1))./(y1(5) - y1(1));
end

%%%%%%%%%%%%%%%%%%
%%
% Np/Nc
val = 60; % 5hrs  % could try 180
mod2_N = {allsGen_N(:,1),allsGen_N(:,2),allsGen_N(:,3),allsGen_N(:,4),allsGen_N(:,5)};
mod40_N = {allsGen_N(:,1),allsGen_N(:,6),allsGen_N(:,7),allsGen_N(:,8),allsGen_N(:,9)};

NpNc = [];
for i = 1:length(mod2_N{1})
    % trajectory
    x_N = [mod2_N{1}{i}(end-val),mod2_N{2}{i}(end-val),mod2_N{3}{i}(end-val),mod2_N{4}{i}(end-val),mod2_N{5}{i}(end-val)];
    y_N = [mod40_N{1}{i}(end-val),mod40_N{2}{i}(end-val),mod40_N{3}{i}(end-val),mod40_N{4}{i}(end-val),mod40_N{5}{i}(end-val)];
    % interpolate between points
    y1_N = interp1(single_xaxis, y_N, m2_xaxis);
    % find Np/Nc
    NpNc(i) = (x_N(5) - x_N(1))./(y1_N(5) - y1_N(1));
end

%%%%%%%%%%%%%%%%%%
%% PLOTTING YPS1 or CMK2
figure(11); 
plot(slopeRatio, NpNc,'m.'); hold on;
axis tight; box off; axis on;
title(['pYPS1 NpNc Plot: ', num2str(numel(mod2{1})),' params']); 
xlabel('slopeRatio'); ylabel('NpNc');
%xlim([1.1 1.6]); ylim([1.1 1.6]);
p = polyfit(slopeRatio, NpNc, 1);
x1 = linspace(min(slopeRatio),max(slopeRatio),30);
y1 = polyval(p,x1);
plot(x1,y1,'k','linewidth',.5);
corrC = corrcoef(slopeRatio, NpNc);
text(1.34,1.55,['y = ', num2str(round(p(1),3)), 'x - ', num2str(abs(round(p(2),3)))]);
text(1.34,1.53, ['corr coeff = ', num2str(corrC(1,2))]);
%text(1.2,1.3,['y = ', num2str(round(p(1),3)), 'x - ', num2str(abs(round(p(2),3)))]);
%text(1.2,1.27, ['corr coeff = ', num2str(corrC(1,2))]);

%% CMK
load(['20200212_outputs_g1vskonTFkoff_plane_pCMK2_N.mat'],'allsGen','allsGen_N')
X2 = allsGen(:,1:12);
Z2 = allsGen_N(:,1:12);
clear allsGen; clear allsGen_N;
X_all=[X2];
Z_all=[Z2];
allsGen=[X_all];
allsGen_N = [Z_all];

%  model end point
val = 60; % 5hrs  % could try 180
mod2 = {allsGen(:,1),allsGen(:,2),allsGen(:,3),allsGen(:,4),allsGen(:,5)};
mod40 = {allsGen(:,1),allsGen(:,6),allsGen(:,7),allsGen(:,8),allsGen(:,9)};
   
% define Xaxis
load('Input_idealizedLight.mat', 'lightInputs')
kfact = 1000;
orig_lightInputs = lightInputs;
m2_xaxis=[orig_lightInputs(1).auc,orig_lightInputs(2).auc./kfact,orig_lightInputs(3).auc./kfact,orig_lightInputs(4).auc./kfact,orig_lightInputs(5).auc./kfact];
single_xaxis=[orig_lightInputs(1).auc,orig_lightInputs(6).auc./kfact,orig_lightInputs(7).auc./kfact,orig_lightInputs(8).auc./kfact,orig_lightInputs(9).auc./kfact];

%%%%%%%%%%%%%%%%%%
% slope ratio
slopeRatio = [];
for i = 1:length(mod2{1})
    % trajectory
    x = [mod2{1}{i}(end-val),mod2{2}{i}(end-val),mod2{3}{i}(end-val),mod2{4}{i}(end-val),mod2{5}{i}(end-val)];
    y = [mod40{1}{i}(end-val),mod40{2}{i}(end-val),mod40{3}{i}(end-val),mod40{4}{i}(end-val),mod40{5}{i}(end-val)];
    % interpolate between points
    y1 = interp1(single_xaxis, y, m2_xaxis);
    % find slope ratio
    slopeRatio(i) = (x(5) - x(1))./(y1(5) - y1(1));
end

%%%%%%%%%%%%%%%%%%
%%
% Np/Nc
val = 60; % 5hrs  % could try 180
mod2_N = {allsGen_N(:,1),allsGen_N(:,2),allsGen_N(:,3),allsGen_N(:,4),allsGen_N(:,5)};
mod40_N = {allsGen_N(:,1),allsGen_N(:,6),allsGen_N(:,7),allsGen_N(:,8),allsGen_N(:,9)};

NpNc = [];
for i = 1:length(mod2_N{1})
    % trajectory
    x_N = [mod2_N{1}{i}(end-val),mod2_N{2}{i}(end-val),mod2_N{3}{i}(end-val),mod2_N{4}{i}(end-val),mod2_N{5}{i}(end-val)];
    y_N = [mod40_N{1}{i}(end-val),mod40_N{2}{i}(end-val),mod40_N{3}{i}(end-val),mod40_N{4}{i}(end-val),mod40_N{5}{i}(end-val)];
    % interpolate between points
    y1_N = interp1(single_xaxis, y_N, m2_xaxis);
    % find Np/Nc
    NpNc(i) = (x_N(5) - x_N(1))./(y1_N(5) - y1_N(1));
end

%%%%%%%%%%%%%%%%%%
%% PLOTTING CMK2
figure(12); 
plot(slopeRatio, NpNc,'b.'); hold on;
axis tight; box off; axis on;
title(['pCMK2 NpNc Plot: ', num2str(numel(mod2{1})),' params']); 
xlabel('slopeRatio'); ylabel('NpNc');
%xlim([1.1 1.6]); ylim([1.1 1.6]);
p = polyfit(slopeRatio, NpNc, 1);
x1 = linspace(min(slopeRatio),max(slopeRatio),30);
y1 = polyval(p,x1);
plot(x1,y1,'k','linewidth',.5);
corrC = corrcoef(slopeRatio, NpNc);
text(1.16,1.12,['y = ', num2str(round(p(1),3)), 'x - ', num2str(abs(round(p(2),3)))]);
text(1.16,1.11, ['corr coeff = ', num2str(corrC(1,2))]);
%text(1.2,1.3,['y = ', num2str(round(p(1),3)), 'x - ', num2str(abs(round(p(2),3)))]);
%text(1.2,1.27, ['corr coeff = ', num2str(corrC(1,2))]);