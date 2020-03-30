function [flag] = ReplottingTFStressPanel_GlucdropSorb_20170712()
flag = 1;
% The goal here is to plot the mean and standard error of the mean for a
% panel of TR and stresses. 

% The experiments are done by stresses across transcription factors, and
% are listed here:
%        1)20150116ZeroGlucPanel15TFs_1-8; 20150116ZeroGlucPanel15TFs_9-16
%        2)20150709_SorbitolStress_ASOE_TRpanel
%        3)20150710_PhosphateDeplet_ASOE_TRpanel

% The values we care about are the mean, the number of cells per time
% point, the standard deviation (and eventual standard error), and the
% requisite processed data that we can extract this information from is
% here:
%        1) nfMeanStd_expname.mat
%        2) expname.mat

% The strategy is to plot each TF in a single plot
%%
% 1)20150116ZeroGlucPanel15TFs_1-8; 20150116ZeroGlucPanel15TFs_9-16
cd('20150116ZeroGlucPanel15TFs_9-16_processed_data\')
load nfMeanStd_ZeroGlucPanel_9-16_experiment.mat % load for mean and std info
load ZeroGlucPanel_9-16_experiment.mat % load for single cell count info
load strainOrder.mat
exp916mean(:,:)=chanPosMean(2,:,:); % 7x90 matrix. 7 positions, 90 time points
exp916std(:,:)=chanPosStd(2,:,:);
exp916time(:,:)=chanPosTim(2,:,:); exp916time = exp916time./2;
exp916cell_info = cell_info;
expstrainOrder = strainOrder;

cell_info_cell = struct2cell(exp916cell_info); % convert structure to interface with above variables
cell_info_cell1(:,:) = cell_info_cell(:,1,:); % info: 1. struct with cell images 2. channel 3. pos 4. tim
cellNumsPerPosTim = []; %zeros(7,90);
for j = 1:numel(cell_info_cell1(1,:)) % construct a 3x90 matrix of number of cells per pos per time
    timeInfo = cell_info_cell1{4,j}/2; % divide the times by 2
    cellNumsPerPosTim(cell_info_cell1{3,j},round(timeInfo)) = numel(cell_info_cell1{1,j});
end
exp916stdErrMean = exp916std./sqrt(cellNumsPerPosTim);% define the exp916stderrmean

%for i = 1:7; figure(i+8); shadedErrorBar(exp916time(i,:), exp916mean(i,:),exp916stdErrMean(i,:)); axis([0,60, 1,3.5]); title(expstrainOrder{i+8}); end; % samples 9-16

cd ..
% 1) 20150116ZeroGlucPanel15TFs_1-8
cd('20150116ZeroGlucPanel15TFs_1-8_processed_data\')
load nfMeanStd_ZeroGlucPanel_1-8_experiment.mat % load for mean and std info
load ZeroGlucPanel_1-8_experiment.mat % load for single cell count info

exp18mean(:,:)=chanPosMean(2,:,:); % 7x90 matrix. 7 positions, 90 time points
exp18std(:,:)=chanPosStd(2,:,:);
exp18time(:,:)=chanPosTim(2,:,:);
exp18cell_info = cell_info;

clear cell_info_cell; clear cell_info_cell1;
cell_info_cell = struct2cell(exp18cell_info); % convert structure to interface with above variables
cell_info_cell1(:,:) = cell_info_cell(:,1,:); % info: 1. struct with cell images 2. channel 3. pos 4. tim
cellNumsPerPosTim = [];%zeros(8,90);
for j = 1:numel(cell_info_cell1(1,:)) % construct a 3x90 matrix of number of cells per pos per time
    timeInfo = cell_info_cell1{4,j}/2; % divide the times by 2
    cellNumsPerPosTim(cell_info_cell1{3,j},round(timeInfo)) = numel(cell_info_cell1{1,j});
end
exp18stdErrMean = exp18std./sqrt(cellNumsPerPosTim);% define the exp916stderrmean

%for i=1:8; figure(i); shadedErrorBar(exp18time(i,:), exp18mean(i,:), exp18stdErrMean(i,:)); axis([0,60, 1,3.5]); title(expstrainOrder{i}); end;

%%
cd ..
% 2)20150709_SorbitolStress_ASOE_TRpanel
cd('20150709_SorbitolStress_ASOE_TRpanel_processed_data\')
load nfMeanStd_SorbStress.mat % load for mean and std info
load SorbStress.mat % load for single cell count info
load strainOrder.mat
% check the structure of "chanPosMean" and "cell_info"
sorbmean(:,:)=chanPosMean(1,:,:); % RFP only % note that there are 30 positions for only 15 samples
sorbstd(:,:)=chanPosStd(1,:,:); 
sorbtime(:,:)=chanPosTim(1,:,:); % no need to divide time by 2 since imaged only in 1 channel
sorbcell_info = cell_info;
sorbstrainOrder = strainOrder;

clear cell_info_cell; clear cell_info_cell1;
cell_info_cell = struct2cell(sorbcell_info); % convert structure to interface with above variables
cell_info_cell1(:,:) = cell_info_cell(:,1,:); % info: 1. struct with cell images 2. channel 3. pos 4. tim
cellNumsPerPosTim = [];%zeros(8,90);
for j = 1:numel(cell_info_cell1(1,:)) % construct a 3x90 matrix of number of cells per pos per time
    %timeInfo = cell_info_cell1{4,j}/2; % no need to divide the times by 2
    timeInfo = cell_info_cell1{4,j};
    cellNumsPerPosTim(cell_info_cell1{3,j},round(timeInfo)) = numel(cell_info_cell1{1,j});
end
sorbstdErrMean = sorbstd./sqrt(cellNumsPerPosTim);% define the stderrmean

% only plot the odd ones (1,3,5,..29)
%for i=1:15; figure(i); shadedErrorBar(sorbtime(2*i-1,:), sorbmean(2*i-1,:), sorbstdErrMean(2*i-1,:)); axis([0,60, 1,3.5]); title(sorbstrainOrder{i}); end;

%%
cd ..
% 3)20150710_PhosphateDeplet_ASOE_TRpanel
cd('20150710_PhosphateDeplet_ASOE_TRpanel_processed_data')
load nfMeanStd_phosphateDeplet.mat % load for mean and std info
load phosphateDeplet.mat % load for single cell count info
load strainOrder.mat

% check the structure of "chanPosMean" and "cell_info"
phosmean(:,:)=chanPosMean(1,:,:); % RFP only % note that there are 30 positions for only 15 samples
phosstd(:,:)=chanPosStd(1,:,:); 
phostime(:,:)=chanPosTim(1,:,:); % no need to divide time by 2 since imaged only in 1 channel
phoscell_info = cell_info;
phosstrainOrder = strainOrder;

clear cell_info_cell; clear cell_info_cell1;
cell_info_cell = struct2cell(phoscell_info); % convert structure to interface with above variables
cell_info_cell1(:,:) = cell_info_cell(:,1,:); % info: 1. struct with cell images 2. channel 3. pos 4. tim
cellNumsPerPosTim = [];%zeros(8,90);
for j = 1:numel(cell_info_cell1(1,:)) % construct a 3x90 matrix of number of cells per pos per time
    %timeInfo = cell_info_cell1{4,j}/2; % no need to divide the times by 2
    timeInfo = cell_info_cell1{4,j};
    cellNumsPerPosTim(cell_info_cell1{3,j},round(timeInfo)) = numel(cell_info_cell1{1,j});
end
phosstdErrMean = phosstd./sqrt(cellNumsPerPosTim);% define the stderrmean

% only plot the odd ones (1,3,5,..29)
%for i=1:15; figure(i); shadedErrorBar(phostime(2*i-1,:), phosmean(2*i-1,:), phosstdErrMean(2*i-1,:)); axis([0,60, 1,3.5]); title(phosstrainOrder{i}); end;

%% Plot Stress Panel
% % currently it is a 3 X 5, but once I do the CaCl2, it will be a 4 X 5.
% % TR order -- Msn2, Msn4, Stb3, Dot6, Crz1, Rtg3
% % Stress order -- Gluc Drop, Osmo shock, Phosphate deplet
% figure(31); 
% subplot(3,5,1); shadedErrorBar(exp18time(1,:), exp18mean(1,:), exp18stdErrMean(1,:)); axis([0,60, 1,3]); title(expstrainOrder{1});
% subplot(3,5,2); shadedErrorBar(exp18time(2,:), exp18mean(2,:), exp18stdErrMean(2,:)); axis([0,60, 1,3]); title(expstrainOrder{2});
% subplot(3,5,3); shadedErrorBar(exp18time(3,:), exp18mean(3,:), exp18stdErrMean(3,:)); axis([0,60, 1,3]); title(expstrainOrder{3});
% subplot(3,5,4); shadedErrorBar(exp18time(5,:), exp18mean(5,:), exp18stdErrMean(5,:)); axis([0,60, 1,3]); title(expstrainOrder{5});
% subplot(3,5,5); shadedErrorBar(exp18time(7,:), exp18mean(7,:), exp18stdErrMean(7,:)); axis([0,60, 1,3]); title(expstrainOrder{7});
% subplot(3,5,6); shadedErrorBar(sorbtime(15,:), sorbmean(15,:), sorbstdErrMean(15,:)); axis([0,60, 1,3]); title(sorbstrainOrder{8});
% subplot(3,5,7); shadedErrorBar(sorbtime(29,:), sorbmean(29,:), sorbstdErrMean(29,:)); axis([0,60, 1,3]); title(sorbstrainOrder{15});
% subplot(3,5,8); shadedErrorBar(sorbtime(19,:), sorbmean(19,:), sorbstdErrMean(19,:)); axis([0,60, 1,3]); title(sorbstrainOrder{10});
% subplot(3,5,9); shadedErrorBar(sorbtime(9,:), sorbmean(9,:), sorbstdErrMean(9,:)); axis([0,60, 1,3]); title(sorbstrainOrder{5});
% subplot(3,5,10); shadedErrorBar(sorbtime(7,:), sorbmean(7,:), sorbstdErrMean(7,:)); axis([0,60, 1,3]); title(sorbstrainOrder{4});
% subplot(3,5,11); shadedErrorBar(phostime(15,:), phosmean(15,:), phosstdErrMean(15,:)); axis([0,60, 1,3]); title(phosstrainOrder{8});
% subplot(3,5,12); shadedErrorBar(phostime(29,:), phosmean(29,:), phosstdErrMean(29,:)); axis([0,60, 1,3]); title(phosstrainOrder{15});
% subplot(3,5,13); shadedErrorBar(phostime(19,:), phosmean(19,:), phosstdErrMean(19,:)); axis([0,60, 1,3]); title(phosstrainOrder{10});
% subplot(3,5,14); shadedErrorBar(phostime(9,:), phosmean(9,:), phosstdErrMean(9,:)); axis([0,60, 1,3]); title(phosstrainOrder{5});
% subplot(3,5,15); shadedErrorBar(phostime(7,:), phosmean(7,:), phosstdErrMean(7,:)); axis([0,60, 1,3]); title(phosstrainOrder{4});

% plot stress panel normalized to lowest value
minExp18mean = min(exp18mean');
minSorbmean = min(sorbmean');
minPhosmean = min(phosmean');
repminExp18mean= repmat(minExp18mean,size(exp18mean,2),1)';
repminSorbmean= repmat(minSorbmean,size(sorbmean,2),1)';
repminPhosmean= repmat(minPhosmean,size(phosmean,2),1)';

normExp18mean = exp18mean./repminExp18mean;
normSorbmean = sorbmean./repminSorbmean;
normPhosmean = phosmean./repminPhosmean;

% plot again, normalized to lowest value
figure(32);
m = 2; n = 5;
subplot(m,n,1); shadedErrorBar(exp18time(1,:), normExp18mean(1,:), exp18stdErrMean(1,:)); axis([0,60, 1,2.1]); title(expstrainOrder{1});
subplot(m,n,2); shadedErrorBar(exp18time(2,:), normExp18mean(2,:), exp18stdErrMean(2,:)); axis([0,60, 1,2.1]); title(expstrainOrder{2});
subplot(m,n,3); shadedErrorBar(exp18time(3,:), normExp18mean(3,:), exp18stdErrMean(3,:)); axis([0,60, 1,2.1]); title(expstrainOrder{3});
subplot(m,n,4); shadedErrorBar(exp18time(5,:), normExp18mean(5,:), exp18stdErrMean(5,:)); axis([0,60, 1,2.1]); title(expstrainOrder{5});
subplot(m,n,5); shadedErrorBar(exp18time(7,:), normExp18mean(7,:), exp18stdErrMean(7,:)); axis([0,60, 1,2.1]); title(expstrainOrder{7});
subplot(m,n,6); shadedErrorBar(sorbtime(15,:), normSorbmean(15,:), sorbstdErrMean(15,:)); axis([0,60, 1,2.1]); title(sorbstrainOrder{8});
subplot(m,n,7); shadedErrorBar(sorbtime(29,:), normSorbmean(29,:), sorbstdErrMean(29,:)); axis([0,60, 1,2.1]); title(sorbstrainOrder{15});
subplot(m,n,8); shadedErrorBar(sorbtime(19,:), normSorbmean(19,:), sorbstdErrMean(19,:)); axis([0,60, 1,2.1]); title(sorbstrainOrder{10});
subplot(m,n,9); shadedErrorBar(sorbtime(9,:), normSorbmean(9,:), sorbstdErrMean(9,:)); axis([0,60, 1,2.1]); title(sorbstrainOrder{5});
subplot(m,n,10); shadedErrorBar(sorbtime(7,:), normSorbmean(7,:), sorbstdErrMean(7,:)); axis([0,60, 1,2.1]); title(sorbstrainOrder{4});
 
end