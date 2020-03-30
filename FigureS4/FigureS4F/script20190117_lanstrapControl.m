function [flag] = script20190117_lanstrapControl()
flag = 1;
% load gene expression
load('syc20190117_lanstrap_trapControl.mat','plate1_20180731');
% plate format: 

plate1_20180322 = plate1_20180731;

%%
% get the means
allgenesMean1 = [];
allgenesStd1 = [];
for i = 1:3
    for j = 1:4
        allgenesMean1(i,j) = nanmean(plate1_20180322{i,j});
        allgenesStd1(i,j) = nanstd(plate1_20180322{i,j});
    end
end

%% Get the mean and std of each sample (bar plots)
one=mean(allgenesMean1(:,1)); sOne=std(allgenesMean1(:,1));
two=mean(allgenesMean1(:,2)); sTwo=std(allgenesMean1(:,2));
three=mean(allgenesMean1(:,3)); sThree=std(allgenesMean1(:,3));
four=mean(allgenesMean1(:,4)); sFour=std(allgenesMean1(:,4));
figure(1); bar([one,two,three,four]); hold on; errorbar([one,two,three,four],[sOne,sTwo,sThree,sFour],'k.');
sampleNames = {'1=pPun1','2=pADH1-Crz1','3=pADH1-Crz1*','4=pADH1-Crz1*-Lanstrap'};
xtickangle(90)
xticklabels({'1','2','3','4'}); ylabel('FITC/SSC');
text(1,1.5,sampleNames);
title('pPUN1-YFP'); axis tight; box off; 


%% 
