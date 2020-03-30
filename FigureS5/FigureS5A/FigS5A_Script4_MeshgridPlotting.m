function [flag] = FigS5A_Script4_MeshgridPlotting()

flag = 1; 
%%
%% load model outputs
load('20190616_outputs_rerunFIG4_pYPS1_9000.mat','allsGen','allsGen_mrna')
X3 = allsGen;
Z3 = allsGen_mrna(1:9000,13);

allsGen=[X3(1:9000,:)];
allsGen_mrna = Z3(1:9000);

%% convert model outputs to form of data
% model Output-Occupancy
mod2 = {allsGen(:,1),allsGen(:,2),allsGen(:,3),allsGen(:,4),allsGen(:,5)};
mod40 = {allsGen(:,1),allsGen(:,6),allsGen(:,7),allsGen(:,8),allsGen(:,9)};

%% load parameter inputs
load('20190616_params_kineticMod_pYPS1_g1pt06_rerunFIG4','paramset')
param10000 = paramset(:,1:9000);

%% load all experimental data
% load output-occupancy data
load('20180605_allgenes_ypsgypcmk_2mTPvsdurTrajs_adhcrz19Azl_ENDPT','GE2m_yps','GEdur_yps','sGE2m_yps','sGEdur_yps')
GE2m = GE2m_yps'; GEdur = GEdur_yps'; sGE2m = sGE2m_yps'; sGEdur = sGEdur_yps';
% define Xaxis for Output-Occupancy data
load('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\LANS OPTOGENETICS PROJECT\SYC Modeling\20180919_kinMod_paramsearch\lightInput_files_20180924\20180924_idealizedlightInputs.mat',...
    'lightInputs')
kfact = 1000;
orig_lightInputs = lightInputs;
m2_xaxis=[orig_lightInputs(1).auc,orig_lightInputs(2).auc./kfact,orig_lightInputs(3).auc./kfact,orig_lightInputs(4).auc./kfact,orig_lightInputs(5).auc./kfact];
single_xaxis=[orig_lightInputs(1).auc,orig_lightInputs(6).auc./kfact,orig_lightInputs(7).auc./kfact,orig_lightInputs(8).auc./kfact,orig_lightInputs(9).auc./kfact];
% define Xaxis for Dose response data
actiTF = [0,orig_lightInputs(11).max./kfact,orig_lightInputs(10).max./kfact,orig_lightInputs(12).max./kfact];


%% scale model data to experimental values
% best fit line through the data FOR ENDPOINT and DOSE RESPONSE data
% Output-Occupancy
val = 60;

p = polyfit(m2_xaxis,GE2m,1); % 2min
x2 = polyval(p,m2_xaxis);
figure(10); errorbar(m2_xaxis,GE2m,sGE2m,'.r','markersize',20,'linewidth',2);...
    hold on; plot(m2_xaxis,x2,'r'); xlabel('cum nuc occpuancy'); ylabel('FITC/SSC');

p = polyfit(single_xaxis,GEdur,1); % 40min
xdur = polyval(p,single_xaxis);
errorbar(single_xaxis,GEdur,sGEdur,'.b','markersize',20,'linewidth',2);...
    plot(single_xaxis,xdur,'b'); axis tight; title('best fit line to data'); box off; axis tight;

% Output-Occupancy & Dose response: scale to max and min of best fit line
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

%% indx_Rats -- sort by slope ratio
% compute slope ratio
rat_slopes = [];
for i = 1:length(mod2{1})
    x = scaled_2_endpt(i,:);
    y = scaled_C_endpt(i,:);
    rat_slopes(i) = ((x(5)-x(1))./(m2_xaxis(5)-m2_xaxis(1)))./((y(5)-y(1))./(single_xaxis(5)-single_xaxis(1)));
end
[sortedRats, indx_Rats] = sort(rat_slopes);

% figure(1); plot(sortedRats); xlabel('parameters'); ylabel('slope ratios'); title('Ratios of Slopes');
% axis tight; box off;

% get allsGen in the order of rat_slopes
% re-organize all of the parameters and trajectories by indx_Rats
allsGen_bR = allsGen(indx_Rats,:);
scaled_2_endpt_bR = scaled_2_endpt(indx_Rats,:);
scaled_C_endpt_bR = scaled_C_endpt(indx_Rats,:);

param10000_bR = param10000(:,indx_Rats);

% set up kon koff values
kon = param10000_bR(1,:); koff = param10000_bR(2,:); 

%% plotting GRID of parameters (by slope ratio) ID regions
% Figure 5B
% grid
LL = min(kon):.1:max(kon); %min(konkoff)
MM = min(koff):.1:max(koff);%max(ronroff);
[x,y] = meshgrid(LL,MM);
AA = kon;
BB = koff;
F = griddata(AA,BB,sortedRats,x,y);

% plot
figure(7); sanePColor(LL,MM,F);  hold on; box off; axis tight; colorbar;
ylim([0.1 5]); xlim([0.1 5]);
title('Slope Ratio Meshgrid'); xlabel('kon'); ylabel('koff');
