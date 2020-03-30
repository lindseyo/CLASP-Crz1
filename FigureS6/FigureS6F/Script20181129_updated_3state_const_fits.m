function [flag] = Script20181129_updated_3state_const_fits()
flag = 1;
%% DON'T CHANGE STUFF -- this is the script to generate PARTS OF FIGURE 5
% load trajectories - set roffTHR at 0.5
load('20181129_3state_const_output10000.mat','allsGen')
X1=allsGen;
load('20181129_3state_const_2pt7_output10000.mat','allsGen')
X2 = allsGen;

allsGen = [X1;X2];

% model end point
val = 60; % 5hrs  % could try 180
mod2 = {allsGen(:,1),allsGen(:,2),allsGen(:,3),allsGen(:,4),allsGen(:,5)};
mod40 = {allsGen(:,1),allsGen(:,6),allsGen(:,7),allsGen(:,8),allsGen(:,9)};
% model transferfunction
mod_transffn= {allsGen(:,1),allsGen(:,11),allsGen(:,10),allsGen(:,12)};

%% Parameter Relationships
%load set roffTHR at 0.5
load('20181112_setTHRpt5_3state_roffTHR_30k.mat','paramset');
Y1 = paramset(:,1:10000);
load('20181130_params30k_genUse_roffTHR0to2pt7.mat','paramset');
Y2 = paramset(:,1:10000);
param10000 = [Y1,Y2];
%param10000 = [Y1,Y2,Y3,Y4];
%param10000 = [Y1,Y3,Y4];
%param10000= [Y1];

%%
% load data (transfer function)
load('20180921_ypsgypcmk_HillFn_4hrs','GEg_4','sGEg_4');
expG4 = GEg_4'; sexpG4 = sGEg_4';
% load data (end point)
load('20181127_gyp_pulsedNconst_data','constGE','sconstGE','pulsedGE','spulsedGE') % newer data set
%GE2m = GE2m_gyp; GEdur = GEdur_gyp; sGE2m = sGE2m_gyp; sGEdur = sGEdur_gyp;
GE2m = pulsedGE; GEdur = constGE; sGE2m = spulsedGE; sGEdur = sconstGE;
% define Xaxis
load('Input_idealizedLight.mat', 'lightInputs')
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
%xx=[.001:7];
%yy = xx.^1.1./(xx.^1.1+100^1.1); % this works, shows it's almost linear
%(only this)
%iyy = (max(doseresp)-min(doseresp))*(yy-min(yy))./(max(yy)-min(yy))+min(doseresp);
%plot(xx,iyy,'m');

% scale to max and min of best fit line
a = min(xdur); b = max(xdur);
e = min(doseresp); f = max(doseresp);
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
    % find error for transf fn
    temp_mod = [mod_transffn{1}{i}(end-val),mod_transffn{2}{i}(end-val),mod_transffn{3}{i}(end-val),mod_transffn{4}{i}(end-val)];
    ix = (f-e).*(temp_mod-min(temp_mod))./(max(temp_mod)-min(temp_mod))+e;
    mod_trans(i,:) = ix;
end

% plot out to check (just as expected, scaling this way b0 needs to be set,
% also need to exactly match traj & dose response) - could see if you
% actually get any fits (scale to max worse than scale to max and min)
modTemp2 = scaled_2_endpt(1,:);
modTempS = scaled_C_endpt(1,:); 
figure(17);
plot(m2_xaxis,modTemp2,'r'); hold on; 
plot(single_xaxis,modTempS,'b');
errorbar(m2_xaxis,GE2m,sGE2m,'.r','markersize',20,'linewidth',2); hold on;
errorbar(single_xaxis,GEdur,sGEdur,'.b','markersize',20,'linewidth',2);
axis tight; box off
xlabel('TF AUC'); ylabel('FITC/SSC');


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

figure(1); plot(sortedRats); xlabel('parameters'); ylabel('slope ratios'); title('Ratios of Slopes');
axis tight; box off;
%% get allsGen in the order of rat_slopes
% re-organize all of the parameters and trajectories by indx_Rats
allsGen_bR = allsGen(indx_Rats,:);
scaled_2_endpt_bR = scaled_2_endpt(indx_Rats,:);
scaled_C_endpt_bR = scaled_C_endpt(indx_Rats,:);
mod_trans_bR = mod_trans(indx_Rats,:);

param10000_bR = param10000(:,indx_Rats);

% set up kon koff values
kon = param10000_bR(1,:); koff = param10000_bR(2,:); ron = param10000_bR(8,:); roff = param10000_bR(9,:);

%% ERROR METRIC
% compute least square error in transfer function + endpt
err_endpt = []; %err_transf = [];
for i = 1:length(scaled_2_endpt_bR)
%for i = 1:length(Nmod_trans);
    % find error for endpoint
    %ix = Nscaled_2_endpt(i,:);
    ix = scaled_2_endpt_bR(i,:);
    %iy = Nscaled_C_endpt(i,:);
    iy = scaled_C_endpt_bR(i,:);
    endpt2(i) = sum((ix - GE2m).^2);
    endptC(i) = sum((iy - GEdur).^2);

    err_endpt(i) = endpt2(i) + endptC(i);

    % find error for transf fn
    %ix = mod_trans(i,:);
    %ix = Nmod_trans(i,:);
    
    %err_transf(i) = sum((ix - expG4).^2);

end
[C,indx_endpt]=sort(err_endpt);
%figure(1); plot(m2_xaxis,scaled_2_endpt(indx_endpt(1),:),'r'); hold on; plot(single_xaxis,scaled_C_endpt(indx_endpt(1),:)); % test out 
%figure(2); plot(C); % see error cut off % currently 30
%[D,indx_trans]=sort(err_transf);
%figure(1); plot(actiTF,mod_trans(indx_trans(80),:),'r');
%figure(4); plot(D); % currently 30

% transfer function also
% compute least square error in transfer function
err_transf = [];
for i = 1:size(mod_trans_bR,1)
    % find error for transf fn
    ix = mod_trans_bR(i,:);

    err_transf(i) = sum((ix - expG4).^2);
end
[D,indx_trans]=sort(err_transf);

%%
% plot out quantitative fits of model to data
% plot the data with model [top50] (endpoint)
paraNum = 1:100;

% metric for fits
indxxx=find(log(err_endpt) < nanmean(log(err_endpt)) - 1.3*nanstd(log(err_endpt))); % 1.3
modTemp2 = scaled_2_endpt_bR(indxxx,:); mm2 = mean(modTemp2); sm2 = std(modTemp2); 
modTempS = scaled_C_endpt_bR(indxxx,:); mSing = mean(modTempS); sSing = std(modTempS);

%modTemp2 = scaled_2_endpt_bR(indx_trans(paraNum),:); mm2 = mean(modTemp2); sm2 = std(modTemp2); 
%modTempS = scaled_C_endpt_bR(indx_trans(paraNum),:); mSing = mean(modTempS); sSing = std(modTempS);
%modTemp2 = scaled_2_endpt_bR(indx_endpt(paraNum),:); mm2 = mean(modTemp2); sm2 = std(modTemp2); 
%modTempS = scaled_C_endpt_bR(indx_endpt(paraNum),:); mSing = mean(modTempS); sSing = std(modTempS);
%modTemp2 = scaled_2_endpt_bR(indxFits,:); mm2 = mean(modTemp2); sm2 = std(modTemp2); 
%modTempS = scaled_C_endpt_bR(indxFits,:); mSing = mean(modTempS); sSing = std(modTempS);
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
title(['pGYP7 10k, top',num2str(numel(indxxx)),'Param'])
%title(['30.5k roffTHR sequ, top',num2str(numel(indxFits)),'Param'])

%% PLOT TRANSFER FUNCTION
% run parameters that fit through a further fit
% transfer function also
% compute least square error in transfer function
nmod_trans = mod_trans_bR(indx_endpt(paraNum),:);
%nmod_trans = mod_trans_bR(indxFits,:);
err_transf_1 = [];
for i = 1:size(nmod_trans,1)
    % find error for transf fn
    ix = nmod_trans(i,:);

    err_transf_1(i) = sum((ix - expG4).^2);
end
[D,indx_trans_1]=sort(err_transf_1);

indyyy = find(err_transf_1 < nanmean(err_transf_1) - 0.8*nanstd(err_transf_1));
meanTransffn=mean(nmod_trans(indyyy,:)); stdTransffn=std(nmod_trans(indyyy,:));
%meanTransffn=mean(nmod_trans(indx_trans_1(1:10),:)); stdTransffn=std(nmod_trans(indx_trans_1(1:10),:));
figure(3); %plot(actiTF,nmod_trans(indx_trans_1(1:5),:),'b'); hold on; 
shadedErrorBar(actiTF,meanTransffn,stdTransffn,'k'); hold on; 
errorbar(actiTF,expG4,sexpG4,'k.','markersize',20,'linewidth',2); axis tight; ylabel('FITC/SSC'); xlabel('max TF'); 
legend('data','model'); box off;
title(['pGYP7 10k, top',num2str(numel(indyyy)),'Param'])

%% MODEL PREDICTION - different promoters
% load data
load('20180430_sum3exp_4genes_4tfproms.mat','mGyp','sGyp','mPun','sPun','mYps','sYps','mCmk','sCmk');
gyp= [mGyp(1:2,:);mGyp(4,:)]; st_gyp=[sGyp(1:2,:);sGyp(4,:)];
%gyp = gyp(:,3:4);
%st_gyp = st_gyp(:,3:4);
gyp = gyp(:,4);
st_gyp = st_gyp(:,4);

% interpolate 
% gyp2a =interp1(m2_xaxis([1,4]),[gyp(1:2,1)],single_xaxis([1,3]));
% gyp40a = [gyp(1,1),gyp(3,1)];
gyp2t =interp1(m2_xaxis([1,4]),[gyp(1:2,1)],single_xaxis([1,3]));
gyp40t = [gyp(1,1),gyp(3,1)];
%gyp = [[gyp2a,gyp40a(2)]',[gyp2t,gyp40t(2)]'];
gyp = [gyp2t,gyp40t(2)]';
% % plot data
% figure(4);
% ydt = [];
% ctr = 1:2; hBar = bar(ctr,gyp',1); 
% for k1=1:size(gyp,1); 
%     k1
%     ctr(k1,:) = bsxfun(@plus,hBar(1).XData,[hBar(k1).XOffset]');
%     ydt(k1,:) = hBar(k1).YData; 
% end
% hold on; errorbar(ctr,ydt,st_gyp,'k.'); hold off;
% axis tight; xticklabels({'adh','tef'}); xtickangle(gca,90); ylabel('GE');
% title('gyp7 data'); ylim([0 0.5]); box off; legend('0m','2m12mX20','40mMID');

% 1=zeros; 2=adh2m20m; 3=adh2m15m; 4=adh2m12m; 5=adh2m6m; 6=adh20_full;
% 7=adh40_full; 8=adh80_full; 9=adh120_full; 10=adh240_full;
% 11=rpl240_full; 12=tef240_full; 13=beg; 14=end; 15=rpl 2m; 16=rpl single;
% 17=tef 2m; 18=tef single; 19= SS

% model end point
val = 60; % 5hrs 
tefGyp = {allsGen_bR(indx_endpt(indx_trans_1(1:10)),1),allsGen_bR(indx_endpt(indx_trans_1(1:10)),17),allsGen_bR(indx_endpt(indx_trans_1(1:10)),18)};

% interpolate
gypa = []; gypt = [];
for gg = 1:length(indx_trans_1(1:10))
    %gyp2a =interp1(m2_xaxis([1,4]),[adhGyp{1}{gg}(end-val),adhGyp{2}{gg}(end-val)],single_xaxis([1,3]));
    %gyp40a = [adhGyp{1}{gg}(end-val),adhGyp{3}{gg}(end-val)];
    gyp2t =interp1(m2_xaxis([1,4]),[tefGyp{1}{gg}(end-val),tefGyp{2}{gg}(end-val)],single_xaxis([1,3]));
    gyp40t = [tefGyp{1}{gg}(end-val),tefGyp{3}{gg}(end-val)];
    %gypa(:,gg) = [gyp2a,gyp40a(2)]';
    gypt(:,gg) = [gyp2t,gyp40t(2)]';
end

% normalize to the data
B = max(gyp(:)); A = min(gyp(:));
iy=[]; iz=[];
for i = 1:length(indx_trans_1(1:10))
    %b = gypa(:,i);
    c = gypt(:,i);
    %iy(:,i) = (B-A).*(b-min(c))./(max(c)-min(c))+A;
    iz(:,i) = (B-A).*(c-min(c))./(max(c)-min(c))+A;
end

% produce error bars
mc = mean(iz'); sc=std(iz');
%abc = [ma',mb',mc']; sabc = [sa',sb',sc'];
%abc = [mb',mc']; sabc = [sb',sc'];
abc = [mc']; sabc = [sc'];
ronTHRabc = abc; ronTHRsabc =sabc;

% combine data and model
totC = [abc,gyp];
% plotting
figure(20);
ydt=[];
ctr = 1:2; hBar = bar(ctr,totC',1); 
for k1=1:size(totC,1); 
    ctr(k1,:) = bsxfun(@plus,hBar(1).XData,[hBar(k1).XOffset]');
    ydt(k1,:) = hBar(k1).YData; 
end
hold on; errorbar(ctr,totC,[sabc,st_gyp],'k.'); hold off;
axis tight; xticklabels({'model','data'}); xtickangle(gca,90); ylabel('GE');
title('Gyp7 Diff Proms - Model top10');  ylim([0 0.5]); box off; %ylim([0 0.7]); 
axis tight; box off; legend('0m','2m12mX20','40mMID');
