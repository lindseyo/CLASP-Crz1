function [flag] = Figure6E_left_heatmap()
flag = 1;

%% Parameter Relationships
%load set roffTHR at 0.5
%load('20181205_fastronslowkon_3state_2RANGE_roffTHR.mat','paramset');
load('20181205_slowronkon_3state_2RANGE_roffTHR.mat','paramset');

%Y3 = paramset(:,1:9500);
Y3 = paramset(:,1:20000);


param10000 = [Y3];

%% DON'T CHANGE STUFF -- this is the script to generate PARTS OF FIGURE 5
% load trajectories - set roffTHR at 0.5
load('20181205_setTHRpt5_RoffTHR_srsk2RANGE_output20000.mat','allsGen')
%load('20181205_setTHRpt5_RoffTHR_frsk2RANGE_output9500.mat','allsGen')

%X3 = allsGen(1:9500,:);
X3 = allsGen(1:20000,:);

allsGen = [X3];

% model end point
val = 60; % 5hrs  % could try 180
mod2 = {allsGen(:,1),allsGen(:,2),allsGen(:,3),allsGen(:,4),allsGen(:,5)};
mod40 = {allsGen(:,1),allsGen(:,6),allsGen(:,7),allsGen(:,8),allsGen(:,9)};
% model transferfunction
mod_transffn= {allsGen(:,1),allsGen(:,11),allsGen(:,10),allsGen(:,12)};

%%
% load data (transfer function)
load('20180921_ypsgypcmk_HillFn_4hrs','GEg_4','sGEg_4');
expG4 = GEg_4'; sexpG4 = sGEg_4';
% load data (end point)
load('20181127_gyp_pulsedNconst_data','constGE','sconstGE','pulsedGE','spulsedGE') % newer data set
GE2m = pulsedGE; GEdur = constGE; sGE2m = spulsedGE; sGEdur = sconstGE;
% define Xaxis
load('Input_idealizedLight.mat','lightInputs')
kfact = 1000;
orig_lightInputs = lightInputs;
m2_xaxis=[orig_lightInputs(1).auc,orig_lightInputs(2).auc./kfact,orig_lightInputs(3).auc./kfact,orig_lightInputs(4).auc./kfact,orig_lightInputs(5).auc./kfact];
single_xaxis=[orig_lightInputs(1).auc,orig_lightInputs(6).auc./kfact,orig_lightInputs(7).auc./kfact,orig_lightInputs(8).auc./kfact,orig_lightInputs(9).auc./kfact];
% define for transfer fn X axis
actiTF = [0,orig_lightInputs(11).max./kfact,orig_lightInputs(10).max./kfact,orig_lightInputs(12).max./kfact];

%%
% scale all trajs to experimental values
% best fit line through the data
% 2min
p = polyfit(m2_xaxis,GE2m,1);
x2 = polyval(p,m2_xaxis);
figure(10); errorbar(m2_xaxis,GE2m,sGE2m,'.r','markersize',20,'linewidth',2);...
    hold on; plot(m2_xaxis,x2,'r'); xlabel('cum nuc occpuancy'); ylabel('FITC/SSC');
% 40min
p = polyfit(single_xaxis,GEdur,1);
xdur = polyval(p,single_xaxis);
errorbar(single_xaxis,GEdur,sGEdur,'.b','markersize',20,'linewidth',2);...
    plot(single_xaxis,xdur,'b'); axis tight; title('best fit line to data');
% dose response
p = polyfit(actiTF,expG4,1);
doseresp = polyval(p,actiTF);
figure(11); errorbar(actiTF,expG4,sexpG4,'b.','markersize',20,'linewidth',2);...
    hold on; plot(actiTF,doseresp,'b'); axis tight; title('best fit line to data');...
    xlabel('cum max nuc occpuancy'); ylabel('FITC/SSC');
% test out how far it deviates from linearity (BASICALLY LINEAR)

% scale to max and min of best fit line
a = min(xdur); b = max(xdur);
e = min(doseresp); f = max(doseresp);
scaled_2_endpt = [];
scaled_C_endpt = [];
mod_trans = [];
for i = 1:length(mod2{1}) 
    i;
    % find error for endpoint
    x = [mod2{1}{i}(end-val),mod2{2}{i}(end-val),mod2{3}{i}(end-val),mod2{4}{i}(end-val),mod2{5}{i}(end-val)];
    y = [mod40{1}{i}(end-val),mod40{2}{i}(end-val),mod40{3}{i}(end-val),mod40{4}{i}(end-val),mod40{5}{i}(end-val)];
    ix = (b-a).*(x-min(y))./(max(y)-min(y))+a;
    iy = (b-a).*(y-min(y))./(max(y)-min(y))+a;
    scaled_2_endpt(i,:) = ix;
    scaled_C_endpt(i,:) = iy; 
    % find error for transf fn
    temp_mod = [mod_transffn{1}{i}(end-val),mod_transffn{2}{i}(end-val),mod_transffn{3}{i}(end-val),mod_transffn{4}{i}(end-val)];
    ix = (f-e).*(temp_mod-min(temp_mod))./(max(temp_mod)-min(temp_mod))+e;
    mod_trans(i,:) = ix;
end

%% indx_Rats -- feed original allsGen through this
% compute slope ratio
rat_slopes = [];
for i = 1:length(mod2{1})
%     x = [mod2{1}{i}(end-val),mod2{2}{i}(end-val),mod2{3}{i}(end-val),mod2{4}{i}(end-val),mod2{5}{i}(end-val)];
%     y = [mod40{1}{i}(end-val),mod40{2}{i}(end-val),mod40{3}{i}(end-val),mod40{4}{i}(end-val),mod40{5}{i}(end-val)];
    x = scaled_2_endpt(i,:);
    y = scaled_C_endpt(i,:);
    rat_slopes(i) = ((x(5)-x(1))./(m2_xaxis(5)-m2_xaxis(1)))./((y(5)-y(1))./(single_xaxis(5)-single_xaxis(1)));
end
[sortedRats, indx_Rats] = sort(rat_slopes);
% find slope ratio < 1
indx_slopeLT1 = find(sortedRats<1);
indx_slopeGT1 = find(sortedRats>1);

% figure(1); plot(sortedRats); xlabel('parameters'); ylabel('slope ratios'); title('Ratios of Slopes');
% axis tight; box off;
%% get allsGen in the order of rat_slopes
% re-organize all of the parameters and trajectories by indx_Rats
allsGen_bR = allsGen(indx_Rats,:);
scaled_2_endpt_bR = scaled_2_endpt(indx_Rats,:);
scaled_C_endpt_bR = scaled_C_endpt(indx_Rats,:);
mod_trans_bR = mod_trans(indx_Rats,:);

param10000_bR = param10000(:,indx_Rats);

% set up kon koff values
kon = param10000_bR(1,:); koff = param10000_bR(2,:); ron = param10000_bR(8,:); roff = param10000_bR(9,:);

%% THIS IS ANOTHER ERROR METRIC - FIT WITHIN ERROR BARS
upp2m = GE2m + sGE2m;
low2m = GE2m - sGE2m;
uppDur = GEdur + sGEdur;
lowDur = GEdur - sGEdur;
% compute least square error in transfer function + endpt
indxFits = []; %err_transf = [];
for i = 1:length(scaled_2_endpt_bR)
%for i = 1:length(Nmod_trans);
    % find error for endpoint
    %ix = Nscaled_2_endpt(i,:);
    ix = scaled_2_endpt_bR(i,:);
    %iy = Nscaled_C_endpt(i,:);
    iy = scaled_C_endpt_bR(i,:);
    
    uppx = ix <= upp2m;
    lowx = ix >= low2m;
    
    if sum(uppx+lowx) >= 10;
        indxFits = [indxFits,i];
    end

    % find error for transf fn
    %ix = mod_trans(i,:);
    %ix = Nmod_trans(i,:);
    
    %err_transf(i) = sum((ix - expG4).^2);

end

%%
% plot out quantitative fits of model to data
% plot the data with model [top50] (endpoint)
paraNum = 1:50;
%modTemp2 = scaled_2_endpt_bR(indx_trans(paraNum),:); mm2 = mean(modTemp2); sm2 = std(modTemp2); 
%modTempS = scaled_C_endpt_bR(indx_trans(paraNum),:); mSing = mean(modTempS); sSing = std(modTempS);
%modTemp2 = scaled_2_endpt_bR(indx_endpt(paraNum),:); mm2 = mean(modTemp2); sm2 = std(modTemp2); 
%modTempS = scaled_C_endpt_bR(indx_endpt(paraNum),:); mSing = mean(modTempS); sSing = std(modTempS);
modTemp2 = scaled_2_endpt_bR(indxFits,:); mm2 = mean(modTemp2); sm2 = std(modTemp2); 
modTempS = scaled_C_endpt_bR(indxFits,:); mSing = mean(modTempS); sSing = std(modTempS);
figure(17);
%plot(m2_xaxis,modTemp2,'r'); hold on; 
%plot(single_xaxis,modTempS,'b');
shadedErrorBar(m2_xaxis, mm2,sm2,'r'); hold on;
shadedErrorBar(single_xaxis, mSing,sSing,'b'); 
%plot(m2_xaxis, mm2,'r*','markersize',20); hold on;
%plot(single_xaxis, mSing,'b*','markersize',20); 
errorbar(m2_xaxis,GE2m,sGE2m,'.r','markersize',20,'linewidth',2); hold on;
errorbar(single_xaxis,GEdur,sGEdur,'.b','markersize',20,'linewidth',2);
axis tight; box off
xlabel('TF AUC'); ylabel('FITC/SSC');
%title(['30.5k roffTHR sequ, top',num2str(paraNum(end)),'Param'])
title(['30.5k roffTHR sequ, top',num2str(numel(indxFits)),'Param'])
%%
%%
% run parameters that fit through a further fit
% transfer function also
% compute least square error in transfer function
%nmod_trans = mod_trans_bR(indx_endpt(paraNum),:);
nmod_trans = mod_trans_bR(indxFits,:);
err_transf_1 = [];
for i = 1:size(nmod_trans,1)
    % find error for transf fn
    ix = nmod_trans(i,:);

    err_transf_1(i) = sum((ix - expG4).^2);
end
[D,indx_trans_1]=sort(err_transf_1);


indyyy = find(err_transf_1 < nanmean(err_transf_1) - .8*nanstd(err_transf_1));
meanTransffn=mean(nmod_trans(indyyy,:)); stdTransffn=std(nmod_trans(indyyy,:));
%meanTransffn=mean(nmod_trans(indx_trans_1(1:10),:)); stdTransffn=std(nmod_trans(indx_trans_1(1:10),:));
figure(3); %plot(actiTF,nmod_trans(indx_trans_1(1:5),:),'b'); hold on; 
shadedErrorBar(actiTF,meanTransffn,stdTransffn,'k'); hold on; 
errorbar(actiTF,expG4,sexpG4,'k.','markersize',20,'linewidth',2); axis tight; 
title(['Dose Resp ErrTHRpt9, top',num2str(numel(indyyy)),'Param'])
ylabel('FITC/SSC'); xlabel('max TF'); 
legend('data','model'); box off;

%% plotting scatter of parameters (by slope ratio) ID regions
% grid
% search for max of sortedRats
gindx = find(log10(kon./koff) >= -2 & log10(kon./koff) <= 0.1 & log10(ron./roff) >= -2 & log10(ron./roff) <= 0.5);
maxColor = max(sortedRats);
minColor = min(sortedRats);
colorV = linspace(minColor,maxColor,256); %256
cMap =[zeros(sum(colorV <= 1),3);parula(256-sum(colorV <= 1))];

LL = min(log10(kon./koff)):.05:max(log10(kon./koff)); %min(konkoff)
MM = min(log10(ron./roff)):.05:max(log10(ron./roff));%max(ronroff);
[x,y] = meshgrid(LL,MM);
AA = log10(kon./koff);
BB = log10(ron./roff);
F = griddata(AA,BB,sortedRats,x,y);
figure(8); sanePColor(LL,MM,F);  hold on; box off; axis tight; colorbar; colormap(cMap); %DEFINITELY CONSIDER THIS!
title('SR w Linear Region + SR < 1 + QualFits'); xlabel('log10(kon/koff)'); ylabel('log10(ron/roff)');
axis([-2 .5 -2 1]);
% axis([-1 0 -1 0.3])
plot(log10(param10000_bR(1,indxFits)./param10000_bR(2,indxFits)),...
    log10(param10000_bR(8,indxFits)./param10000_bR(9,indxFits)),'m.'); hold on;
plot(log10(param10000_bR(1,indxFits(indyyy))./param10000_bR(2,indxFits(indyyy))),...
    log10(param10000_bR(8,indxFits(indyyy))./param10000_bR(9,indxFits(indyyy))),'g.');
