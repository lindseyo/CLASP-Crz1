function [flag] = script20170908_growthExp_24HRmanual512LightExp_LShark()
flag = 1;
%%
di = '20170908_24HRmanual512LightExp';
timePts = dir([di,'\*endpoint*']);

catMatrix=[];
for i=1:length(timePts);
    filename = timePts(i).name;
    [finalMatrix] = endPointOD20170907(di, filename);
    catMatrix(i,:,:) = finalMatrix;
end

avgMatrixWT = nanmean(catMatrix(:,1:3,:),2); odWT(:,:) = avgMatrixWT(:,1,:);
avgMatrixAla = nanmean(catMatrix(:,4:6,:),2); odALA(:,:) = avgMatrixAla(:,1,:);
catMatrix1=permute(catMatrix,[2,1,3]);
stdMatrixWT = nanstd(catMatrix1(1:3,:,:)); odstdWT(:,:) = stdMatrixWT(1,:,:);
stdMatrixAla = nanstd(catMatrix1(4:6,:,:)); odstdALA(:,:) = stdMatrixAla(1,:,:);

%% sample keys
% row 1:3 --> WT, row 4:6 --> Ala
% col 1:5 --> no light, col 6:10 --> light
% samples row 1:3 1. CrzWT 2. CrzWT-LANS 3. zdk-CrzWT-LANS(265) 4.
% zdk-CrzWT-LANS(264) 5. zdk-CrzWT-LANS(263)
% samples row 4:6 1. Ala 2. Ala-LANS 3. zdk-Ala-LANS(265) 4.
% zdk-Ala-LANS(263) 5. zdk-Ala(263)
% rows 7:8 are WT:Ala with 0.4M CaCl2. discarded this data.
%% PLOTTING AVERAGE ODS
timeVector = [0 1 2 3 4 5 6 7 8 9 11 13 14 15 16 22 23 25];
samplesWT = [];
samplesALA = [];
%% Comparing Crz light/nolight - not too much difference between light and no light
% compare WT light/nolight
figure(1); 
A=shadedErrorBar(timeVector, odWT(:,1), odstdWT(:,1),'r',1); hold on; B=shadedErrorBar(timeVector, odWT(:,6), odstdWT(:,6),'b',1);
legend([A.mainLine, B.mainLine], 'no light','light')
xlabel('time (hours)'); ylabel('OD600'); title('pADH1-CRZ1'); box off; axis tight;
