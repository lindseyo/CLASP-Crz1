function [flag] = FigureS4IJ();
flag = 1;
% load gene expression
load('syc20180605_ypsgypcmk_2mTPvs40mtrajs_adhcrz19Azl.mat')

% plate format:
% B1-3 0min yps
% B4-6 0min gyp
% B7-9 0min cmk

% C1-3 2m20mX12 = 20 yps
% C4-6 2m15mX16 = 32
% C7-9 2m12mX20 = 40
% C10-12 2m6mX40 = 80
% D1-3 20m
% D4-6 40m
% D7-9 80m
% D10-12 120m
% E1-3 2m20mX12 = 20 gyp
% E4-6 2m15mX16 = 32
% E7-9 2m12mX20 = 40
% E10-12 2m6mX40 = 80
% F1-3 20m
% F4-6 40m
% F7-9 80m
% F10-12 120m
% G1-3 2m20mX12 = 20 cmk
% G4-6 2m15mX16 = 32
% G7-9 2m12mX20 = 40
% G10-12 2m6mX40 = 80
% H1-3 20m
% H4-6 40m
% H7-9 80m
% H10-12 120m

% define Xaxis
load('Input_idealizedLight.mat', 'lightInputs')
kfact = 1000;
orig_lightInputs = lightInputs;
m2_xaxis=[orig_lightInputs(1).auc,orig_lightInputs(2).auc./kfact,orig_lightInputs(3).auc./kfact,orig_lightInputs(4).auc./kfact,orig_lightInputs(5).auc./kfact];
single_xaxis=[orig_lightInputs(1).auc,orig_lightInputs(6).auc./kfact,orig_lightInputs(7).auc./kfact,orig_lightInputs(8).auc./kfact,orig_lightInputs(9).auc./kfact];

plate1_20180322 = plate1_20180605_ypsgypcmk;

%%
% get the means
allgenesMean = [];
allgenesStd = [];
for i = 1:7
    for j = 1:12
        allgenesMean(i,j) = nanmean(plate1_20180322{i,j});
        allgenesStd(i,j) = nanstd(plate1_20180322{i,j});
    end
end
% individual samples (yps)
y_0 = mean(allgenesMean(1,1:3)); sy_0 = std(allgenesMean(1,1:3));
y2m20m = mean(allgenesMean(2,1:3)); sy2m20m = std(allgenesMean(2,1:3)); 
y2m15m = mean(allgenesMean(2,4:6)); sy2m15m = std(allgenesMean(2,4:6));
y2m12m = mean(allgenesMean(2,7:9)); sy2m12m = std(allgenesMean(2,7:9));
y2m6m = mean(allgenesMean(2,10:12)); sy2m6m = std(allgenesMean(2,10:12));
y20m = mean([allgenesMean(3,1:3)]); sy20m = std([allgenesMean(3,1:3)]); 
y40m = mean(allgenesMean(3,4:6)); sy40m = std(allgenesMean(3,4:6));
y80m = mean([allgenesMean(3,7:9)]); sy80m = std([allgenesMean(3,7:9)]);
y120m = mean(allgenesMean(3,10:12)); sy120m = std(allgenesMean(3,10:12));
% (gyp)
g_0 = mean(allgenesMean(1,4:6)); sg_0 = std(allgenesMean(1,4:6));
g2m20m = mean(allgenesMean(4,1:3)); sg2m20m = std(allgenesMean(4,1:3)); 
g2m15m = mean(allgenesMean(4,4:6)); sg2m15m = std(allgenesMean(4,4:6));
g2m12m = mean(allgenesMean(4,7:9)); sg2m12m = std(allgenesMean(4,7:9));
g2m6m = mean(allgenesMean(4,10:12)); sg2m6m = std(allgenesMean(4,10:12));
g20m = mean([allgenesMean(5,1:3)]); sg20m = std([allgenesMean(5,1:3)]); 
g40m = mean(allgenesMean(5,4:6)); sg40m = std(allgenesMean(5,4:6));
g80m = mean([allgenesMean(5,7:9)]); sg80m = std([allgenesMean(5,7:9)]);
g120m = mean(allgenesMean(5,10:12)); sg120m = std(allgenesMean(5,10:12));
% (cmk)
c_0 = mean(allgenesMean(1,7:9)); sc_0 = std(allgenesMean(1,7:9));
c2m20m = mean(allgenesMean(6,1:3)); sc2m20m = std(allgenesMean(6,1:3)); 
c2m15m = mean(allgenesMean(6,4:6)); sc2m15m = std(allgenesMean(6,4:6));
c2m12m = mean(allgenesMean(6,7:9)); sc2m12m = std(allgenesMean(6,7:9));
c2m6m = mean(allgenesMean(6,10:12)); sc2m6m = std(allgenesMean(6,10:12));
c20m = mean([allgenesMean(7,1:3)]); sc20m = std([allgenesMean(7,1:3)]); 
c40m = mean(allgenesMean(7,4:6)); sc40m = std(allgenesMean(7,4:6));
c80m = mean([allgenesMean(7,7:9)]); sc80m = std([allgenesMean(7,7:9)]);
c120m = mean(allgenesMean(7,10:12)); sc120m = std(allgenesMean(7,10:12));
%means
GE2m_yps =[y_0;y2m20m;y2m15m;y2m12m;y2m6m]; sGE2m_yps =[sy_0;sy2m20m;sy2m15m;sy2m12m;sy2m6m];
GEdur_yps =[y_0;y20m;y40m;y80m;y120m]; sGEdur_yps =[sy_0;sy20m;sy40m;sy80m;sy120m]; 
GE2m_gyp =[g_0;g2m20m;g2m15m;g2m12m;g2m6m]; sGE2m_gyp =[sg_0;sg2m20m;sg2m15m;sg2m12m;sg2m6m];
GEdur_gyp =[g_0;g20m;g40m;g80m;g120m]; sGEdur_gyp =[sg_0;sg20m;sg40m;sg80m;sg120m];
GE2m_cmk =[c_0;c2m20m;c2m15m;c2m12m;c2m6m]; sGE2m_cmk =[sc_0;sc2m20m;sc2m15m;sc2m12m;sc2m6m];
GEdur_cmk =[c_0;c20m;c40m;c80m;c120m]; sGEdur_cmk =[sc_0;sc20m;sc40m;sc80m;sc120m];
%save genes
%save('20180605_allgenes_ypsgypcmk_2mTPvsdurTrajs_adhcrz19Azl.mat','GE2m_yps','sGE2m_yps','GEdur_yps','sGEdur_yps',...
%    'GE2m_gyp','sGE2m_gyp','GEdur_gyp','sGEdur_gyp','GE2m_cmk','sGE2m_cmk','GEdur_cmk','sGEdur_cmk');

% % plot all data
% figure(1); %plot(AUC2traj,GE2m_yps,'r.','linewidth',2); hold on; 
% errorbar(m2_xaxis,GE2m_yps,sGE2m_yps,'r.','markersize',20,'linewidth',2); hold on;
% %errorbar(AUCdurtraj,GEdur_yps,sGEdur_yps,'b.','markersize',20,'linewidth',2); 
% errorbar(single_xaxis,GEdur_yps,sGEdur_yps,'b.','markersize',20,'linewidth',2); 
% axis tight; box off; title('yps1'); ylabel('FITC/SSC'); xlabel('TF AUC');

% Plot pCMK2
figure(2); %plot(AUC2traj,GE2m_gyp,'r.','linewidth',2); hold on; 
errorbar(m2_xaxis,GE2m_gyp,sGE2m_gyp,'r.','markersize',20,'linewidth',2); hold on;
errorbar(single_xaxis,GEdur_gyp,sGEdur_gyp,'b.','markersize',20,'linewidth',2); 
axis tight; box off; title('gyp7'); ylabel('FITC/SSC'); xlabel('TF AUC');
legend('2minTP','durs')

% plot pGYP7
figure(3); %plot(AUC2traj,GE2m_cmk,'r.','linewidth',2); hold on; 
errorbar(m2_xaxis,GE2m_cmk,sGE2m_cmk,'r.','markersize',20,'linewidth',2); hold on;
errorbar(single_xaxis,GEdur_cmk,sGEdur_cmk,'b.','markersize',20,'linewidth',2); 
axis tight; box off; title('cmk2'); ylabel('FITC/SSC'); xlabel('TF AUC');
legend('2minTP','durs')
%%
%% plot the new reproducible set from LO
% pCMK2 from Lindsey
load('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\LANS OPTOGENETICS PROJECT\SYC Data\syc20180605_ypsgypcmk_2mTPvs40mtrajs_adhcrz19Azl\20181015_cmk2_newLO.mat');
figure(4);
errorbar(m2_xaxis,new_mCmk2_2,new_sCmk2_2,'r.','markersize',20,'linewidth',2); hold on;
errorbar(single_xaxis,new_mCmk2_const,new_sCmk2_const,'b.','markersize',20,'linewidth',2); 
axis tight; box off; title('LO cmk2'); ylabel('FITC/SSC'); xlabel('TF AUC');
legend('2minTP','durs')

%% directly compare the 2 experiments (pCMK2)
figure(5); errorbar(m2_xaxis,GE2m_cmk,sGE2m_cmk,'r.','markersize',20,'linewidth',2); hold on;
errorbar(single_xaxis,GEdur_cmk,sGEdur_cmk,'b.','markersize',20,'linewidth',2);
errorbar(m2_xaxis,new_mCmk2_2,new_sCmk2_2,'m>','markersize',10,'linewidth',2); hold on;
errorbar(single_xaxis,new_mCmk2_const,new_sCmk2_const,'g>','markersize',10,'linewidth',2); axis tight; box off; title('cmk2 comparison'); ylabel('FITC/SSC'); xlabel('TF AUC');
legend('2minTP-O','durs-O','2minTP-N','durs-N')