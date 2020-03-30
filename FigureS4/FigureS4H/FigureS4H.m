function [flag] = FigureS4H()
flag = 1;
% load gene expression
load('syc20180519_yps1_pAdh-crz19A_2minVS10VSdur.mat')

% plate format:
% all yps adh1
% A1-12 0min
% B1-3 2m20mX12 = 20
% B4-6 2m15mX16 = 32
% B7-9 2m12mX20 = 40
% B10-12 2m6mX40 = 80
% C1-3 10m120mX2 = 20
% C4-6 10m60mX4 = 40
% C7-9 10m30mX8 = 80
% C10-12 10m20mX12 = 120
% D1-3 20m
% D4-6 32m
% D7-9 48m
% D10-12 80m
% E1-3 20m
% E4-6 40m
% E7-9 80m
% E10-12 120m

% define Xaxis
load('Input_idealizedLight.mat', 'lightInputs')
kfact = 1000;
orig_lightInputs = lightInputs;
m2_xaxis=[orig_lightInputs(1).auc,orig_lightInputs(2).auc./kfact,orig_lightInputs(3).auc./kfact,orig_lightInputs(4).auc./kfact,orig_lightInputs(5).auc./kfact];
single_xaxis=[orig_lightInputs(1).auc,orig_lightInputs(6).auc./kfact,orig_lightInputs(7).auc./kfact,orig_lightInputs(8).auc./kfact,orig_lightInputs(9).auc./kfact];


plate1_20180322 = plate1_20180519;

%%
% get the means
allgenesMean = [];
allgenesStd = [];
for i = 1:5
    for j = 1:12
        allgenesMean(i,j) = nanmean(plate1_20180322{i,j});
        allgenesStd(i,j) = nanstd(plate1_20180322{i,j});
    end
end
% individual samples
yps_0 = mean(allgenesMean(1,1:12)); syps_0 = std(allgenesMean(1,1:12));
yps2m20m = mean(allgenesMean(2,1:3)); syps2m20m = std(allgenesMean(2,1:3)); 
yps2m15m = mean(allgenesMean(2,4:6)); syps2m15m = std(allgenesMean(2,4:6));
yps2m12m = mean(allgenesMean(2,7:9)); syps2m12m = std(allgenesMean(2,7:9));
yps2m6m = mean(allgenesMean(2,10:12)); syps2m6m = std(allgenesMean(2,10:12));
yps10m120m = mean(allgenesMean(3,1:3)); syps10m120m = std(allgenesMean(3,1:3)); 
yps10m60m = mean(allgenesMean(3,4:6)); syps10m60m = std(allgenesMean(3,4:6));
yps10m30m = mean(allgenesMean(3,7:9)); syps10m30m = std(allgenesMean(3,7:9));
yps10m20m = mean(allgenesMean(3,10:12)); syps10m20m = std(allgenesMean(3,10:12));
yps20m = mean([allgenesMean(4,1:3),allgenesMean(5,1:3)]); syps20m = std([allgenesMean(4,1:3),allgenesMean(5,1:3)]); 
yps32m = mean(allgenesMean(4,4:6)); syps32m = std(allgenesMean(4,4:6));
yps40m = mean(allgenesMean(5,4:6)); syps40m = std(allgenesMean(5,4:6));
yps48m = mean(allgenesMean(4,7:9)); syps48m = std(allgenesMean(4,7:9));
yps80m = mean([allgenesMean(4,10:12),allgenesMean(5,7:9)]); syps80m = std([allgenesMean(4,10:12),allgenesMean(5,7:9)]);
yps120m = mean(allgenesMean(5,10:12)); syps120m = std(allgenesMean(5,10:12));
% means
GE2m_yps =[yps_0;yps2m20m;yps2m15m;yps2m12m;yps2m6m]; sGE2m_yps =[syps_0;syps2m20m;syps2m15m;syps2m12m;syps2m6m];
GE10m_yps =[yps_0;yps10m120m;yps10m60m;yps10m30m;yps10m20m]; sGE10m_yps =[syps_0;syps10m120m;syps10m60m;syps10m30m;syps10m20m];
GEdur_yps =[yps_0;yps20m;yps32m;yps40m;yps48m;yps80m;yps120m]; sGEdur_yps =[syps_0;syps20m;syps32m;syps40m;syps48m;syps80m;syps120m];
%save genes
%save('allgenes20180519_yps1_pAdh1-crz19A_2minVS10VSdur.mat','GE2m_yps','sGE2m_yps','GE10m_yps','sGE10m_yps','GEdur_yps','sGEdur_yps');
%%
% plot all data
figure(1); 
errorbar(m2_xaxis,GE2m_yps,sGE2m_yps,'r.','markersize',20,'linewidth',2); hold on; 
errorbar(single_xaxis,GEdur_yps([1,2,4,6,7]),sGEdur_yps([1,2,4,6,7]),'b.','markersize',20,'linewidth',2); 
axis tight; box off; title('pAdh1-yps1'); ylabel('FITC/SSC'); xlabel('TF AUC');
legend('2minTP','Const')

