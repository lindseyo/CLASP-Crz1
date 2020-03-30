function [flag] = Script20181018_plotting_slopeRatios_forFIG3C()
flag = 1;
%% load TF nuclear localization
% load light inputs to calculate TF AUC
load('adhTFnuctraces_20180516.mat','adh2_p1_trace20180516','adh40_p3_trace20180516');
adh2= adh2_p1_trace20180516; adh40 = adh40_p3_trace20180516; % raw traces that need to be background subtracted or leave as is!
% define adh initials
adh2= adh2(25:37); %figure(1); plot(adh2);
adh40= adh40(1:end-8); %figure(2); plot(adh40);
%base shift if needed
baseshift = adh40(1)-adh2(1);
% zero integral
AUC_adh_0 = 0;%trapz(adh2(1).*ones(1,480))./2; % zero - done
% create trajs for these light inputs
% 2min
% single
adh2m20m = repmat([adh2(1).*ones(1,20*2-length(adh2)),adh2],[1,12]);  adh2m20m=adh2m20m-adh2m20m(1); AUC2_1 = trapz(adh2m20m)./2;
adh2m15m = repmat([adh2(1).*ones(1,15*2-length(adh2)),adh2],[1,16]); adh2m15m=adh2m15m-adh2m15m(1); AUC2_2 = trapz(adh2m15m)./2;
adh2m12m = repmat([adh2(1).*ones(1,12*2-length(adh2)),adh2],[1,20]); adh2m12m=adh2m12m-adh2m12m(1); AUC2_3 = trapz(adh2m12m)./2;
adh2m6m = repmat([adh2(1:12)],[1,40]); adh2m6m=adh2m6m-adh2m6m(1); AUC2_4 = trapz(adh2m6m)./2; 
AUC2traj = [AUC_adh_0,AUC2_1,AUC2_2,AUC2_3,AUC2_4];
% durs % note when the light shuts off it takes 5.5min to shut off, note
% when the light turns on it take 2.5 min to reach max
adh20= [adh40(1:40),adh40(end-11:end)]; 
%figure(10); plot(adh20);
adh40= adh40;
adh80= [adh40(1:end-11),adh40(40:80),adh40(40:77),adh40(end-11:end)];
adh120= [adh40(1:end-11),adh40(40:80),adh40(40:80),adh40(40:80),adh40(40:75), adh40(end-11:end)];
% ful trajs
adh20_full = [adh20(1).*ones(1,110*2),adh20,adh20(1).*ones(1,480-length(adh20)-110*2)]; adh20_full=adh20_full-adh20_full(1); AUC_1 = trapz(adh20_full)./2;
adh40_full = [adh40(1).*ones(1,100*2),adh40,adh40(1).*ones(1,480-length(adh40)-100*2)]; adh40_full=adh40_full-adh40_full(1); AUC_3 = trapz(adh40_full)./2;
adh80_full = [adh80(1).*ones(1,80*2),adh80,adh80(1).*ones(1,480-length(adh80)-80*2)]; adh80_full=adh80_full-adh80_full(1); AUC_5 = trapz(adh80_full)./2;
adh120_full = [adh120(1).*ones(1,60*2),adh120,adh120(1).*ones(1,480-length(adh120)-60*2)]; adh120_full=adh120_full-adh120_full(1); AUC_6 = trapz(adh120_full)./2;
AUCdurtraj = [AUC_adh_0,AUC_1,AUC_3,AUC_5,AUC_6];
%% load in 3 gene dataset
load('syc20180605_ypsgypcmk_2mTPvs40mtrajs_adhcrz19Azl.mat')
plate1_20180322 = plate1_20180605_ypsgypcmk;
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
y_0 = allgenesMean(1,1:3); 
y2m20m = allgenesMean(2,1:3); 
y2m15m = allgenesMean(2,4:6); 
y2m12m = allgenesMean(2,7:9); 
y2m6m = allgenesMean(2,10:12); 
y20m = [allgenesMean(3,1:3)];  
y40m = allgenesMean(3,4:6); 
y80m = [allgenesMean(3,7:9)];
y120m = allgenesMean(3,10:12); 
yps_2all = [y_0;y2m20m;y2m15m;y2m12m;y2m6m];
yps_constall = [y_0;y20m;y40m;y80m;y120m];
% (gyp)
g_0 = allgenesMean(1,4:6);
g2m20m = allgenesMean(4,1:3); 
g2m15m = allgenesMean(4,4:6); 
g2m12m = allgenesMean(4,7:9); 
g2m6m = allgenesMean(4,10:12); 
g20m = allgenesMean(5,1:3);
g40m = allgenesMean(5,4:6); 
g80m = allgenesMean(5,7:9);
g120m = allgenesMean(5,10:12);
gyp_2all = [g_0;g2m20m;g2m15m;g2m12m;g2m6m];
gyp_constall = [g_0;g20m;g40m;g80m;g120m];
% (cmk)
c_0 = allgenesMean(1,7:9); 
c2m20m = allgenesMean(6,1:3); 
c2m15m = allgenesMean(6,4:6); 
c2m12m = allgenesMean(6,7:9); 
c2m6m = allgenesMean(6,10:12); 
c20m = [allgenesMean(7,1:3)]; 
c40m = allgenesMean(7,4:6); 
c80m = [allgenesMean(7,7:9)]; 
c120m = allgenesMean(7,10:12);
cmk_2all = [c_0;c2m20m;c2m15m;c2m12m;c2m6m];
cmk_constall = [c_0;c20m;c40m;c80m;c120m];
%% load in 8 gene dataset
load('syc20180605_8targetgenes_2m20mX20vs40m_adhcrz19Azl.mat')
plate1_20180322 = plate1_20180605_8genes;
% get the means
allgenesMean = [];
allgenesStd = [];
for i = 1:6
    for j = 1:12
        allgenesMean(i,j) = nanmean(plate1_20180322{i,j});
        allgenesStd(i,j) = nanstd(plate1_20180322{i,j});
    end
end
% individual samples yps1
p_0 = allgenesMean(1,4:6); 
p2m12m = allgenesMean(2,4:6); 
p40m = allgenesMean(3,4:6); 
pun_2all = [p_0; p2m12m];
pun_constall = [p_0; p40m];

g_0 = allgenesMean(1,7:9); 
g2m12m = allgenesMean(2,7:9); 
g40m = allgenesMean(3,7:9); 
gyp_2all1 = [g_0; g2m12m];
gyp_constall1 = [g_0; g40m];

pu_0 = allgenesMean(1,10:12); 
pu2m12m = allgenesMean(2,10:12); 
pu40m = allgenesMean(3,10:12);
put_2all = [pu_0; pu2m12m];
put_constall = [pu_0; pu40m];

c_0 = allgenesMean(4,1:3);
c2m12m = allgenesMean(5,1:3);
c40m = allgenesMean(6,1:3);
cmk_2all1 = [c_0; c2m12m];
cmk_constall1 = [c_0; c40m];

e_0 = allgenesMean(4,4:6); 
e2m12m = allgenesMean(5,4:6); 
e40m = allgenesMean(6,4:6); 
ena_2all = [e_0; e2m12m];
ena_constall = [e_0; e40m];

m_0 = allgenesMean(4,7:9); 
m2m12m = allgenesMean(5,7:9);  
m40m = allgenesMean(6,7:9); 
mep_2all = [m_0; m2m12m];
mep_constall = [m_0; m40m];

pmc_0 = allgenesMean(4,10:12); 
pmc2m12m = allgenesMean(5,10:12);
pmc40m = allgenesMean(6,10:12); 
pmc_2all = [pmc_0; pmc2m12m];
pmc_constall = [pmc_0; pmc40m];


%% calculate slope ratio
% for 3 genes - via best fit line
fyps_2all = []; fyps_constall = []; yps_slope = [];
for i = 1:3
coeffs = polyfit(AUC2traj,yps_2all(:,i)',1);
fyps_2all(:,i) = polyval(coeffs, AUC2traj);
coeffs = polyfit(AUCdurtraj,yps_constall(:,i)',1);
fyps_constall(:,i) = polyval(coeffs, AUCdurtraj);
yps_slope(i) = ((fyps_2all(5,i)-fyps_2all(1,i))./(AUC2traj(5) - AUC2traj(1)))./((fyps_constall(5,i)-fyps_constall(1,i))./(AUCdurtraj(5)-AUCdurtraj(1)))
end

fgyp_2all = []; fgyp_constall = []; gyp_slope = [];
for i = 1:3
coeffs = polyfit(AUC2traj,gyp_2all(:,i)',1);
fgyp_2all(:,i) = polyval(coeffs, AUC2traj);
coeffs = polyfit(AUCdurtraj,gyp_constall(:,i)',1);
fgyp_constall(:,i) = polyval(coeffs, AUCdurtraj);
gyp_slope(i) = ((fgyp_2all(5,i)-fgyp_2all(1,i))./(AUC2traj(5) - AUC2traj(1)))./((fgyp_constall(5,i)-fgyp_constall(1,i))./(AUCdurtraj(5)-AUCdurtraj(1)))
end

fcmk_2all = []; fcmk_constall = []; cmk_slope = [];
for i = 1:3
coeffs = polyfit(AUC2traj,cmk_2all(:,i)',1);
fcmk_2all(:,i) = polyval(coeffs, AUC2traj);
coeffs = polyfit(AUCdurtraj,cmk_constall(:,i)',1);
fcmk_constall(:,i) = polyval(coeffs, AUCdurtraj);
cmk_slope(i) = ((fcmk_2all(5,i)-fcmk_2all(1,i))./(AUC2traj(5) - AUC2traj(1)))./((fcmk_constall(5,i)-fcmk_constall(1,i))./(AUCdurtraj(5)-AUCdurtraj(1)))
end

% for 8 genes
fpun_2all = []; fpun_constall = []; pun_slope = [];
for i = 1:3
coeffs = polyfit(AUC2traj([1,3]),pun_2all(:,i)',1);
fpun_2all(:,i) = polyval(coeffs, AUC2traj([1,3]));
coeffs = polyfit(AUCdurtraj([1,3]),pun_constall(:,i)',1);
fpun_constall(:,i) = polyval(coeffs, AUCdurtraj([1,3]));
pun_slope(i) = ((fpun_2all(2,i)-fpun_2all(1,i))./(AUC2traj(3) - AUC2traj(1)))./((fpun_constall(2,i)-fpun_constall(1,i))./(AUCdurtraj(3)-AUCdurtraj(1)))
end

fena_2all = []; fena_constall = []; ena_slope = [];
for i = 1:3
coeffs = polyfit(AUC2traj([1,3]),ena_2all(:,i)',1);
fena_2all(:,i) = polyval(coeffs, AUC2traj([1,3]));
coeffs = polyfit(AUCdurtraj([1,3]),ena_constall(:,i)',1);
fena_constall(:,i) = polyval(coeffs, AUCdurtraj([1,3]));
ena_slope(i) = ((fena_2all(2,i)-fena_2all(1,i))./(AUC2traj(3) - AUC2traj(1)))./((fena_constall(2,i)-fena_constall(1,i))./(AUCdurtraj(3)-AUCdurtraj(1)))
end

fput_2all = []; fput_constall = []; put_slope = [];
for i = 1:3
coeffs = polyfit(AUC2traj([1,3]),put_2all(:,i)',1);
fput_2all(:,i) = polyval(coeffs, AUC2traj([1,3]));
coeffs = polyfit(AUCdurtraj([1,3]),put_constall(:,i)',1);
fput_constall(:,i) = polyval(coeffs, AUCdurtraj([1,3]));
put_slope(i) = ((fput_2all(2,i)-fput_2all(1,i))./(AUC2traj(3) - AUC2traj(1)))./((fput_constall(2,i)-fput_constall(1,i))./(AUCdurtraj(3)-AUCdurtraj(1)))
end

fmep_2all = []; fmep_constall = []; mep_slope = [];
for i = 1:3
coeffs = polyfit(AUC2traj([1,3]),mep_2all(:,i)',1);
fmep_2all(:,i) = polyval(coeffs, AUC2traj([1,3]));
coeffs = polyfit(AUCdurtraj([1,3]),mep_constall(:,i)',1);
fmep_constall(:,i) = polyval(coeffs, AUCdurtraj([1,3]));
mep_slope(i) = ((fmep_2all(2,i)-fmep_2all(1,i))./(AUC2traj(3) - AUC2traj(1)))./((fmep_constall(2,i)-fmep_constall(1,i))./(AUCdurtraj(3)-AUCdurtraj(1)))
end

fpmc_2all = []; fpmc_constall = []; pmc_slope = [];
for i = 1:3
coeffs = polyfit(AUC2traj([1,3]),pmc_2all(:,i)',1);
fpmc_2all(:,i) = polyval(coeffs, AUC2traj([1,3]));
coeffs = polyfit(AUCdurtraj([1,3]),pmc_constall(:,i)',1);
fpmc_constall(:,i) = polyval(coeffs, AUCdurtraj([1,3]));
pmc_slope(i) = ((fpmc_2all(2,i)-fpmc_2all(1,i))./(AUC2traj(3) - AUC2traj(1)))./((fpmc_constall(2,i)-fpmc_constall(1,i))./(AUCdurtraj(3)-AUCdurtraj(1)))
end

fgyp_2all1 = []; fgyp_constall1 = []; gyp_slope1 = [];
for i = 1:3
coeffs = polyfit(AUC2traj([1,3]),gyp_2all1(:,i)',1);
fgyp_2all1(:,i) = polyval(coeffs, AUC2traj([1,3]));
coeffs = polyfit(AUCdurtraj([1,3]),gyp_constall1(:,i)',1);
fgyp_constall1(:,i) = polyval(coeffs, AUCdurtraj([1,3]));
gyp_slope1(i) = ((fgyp_2all1(2,i)-fgyp_2all1(1,i))./(AUC2traj(3) - AUC2traj(1)))./((fgyp_constall1(2,i)-fgyp_constall1(1,i))./(AUCdurtraj(3)-AUCdurtraj(1)))
end

fcmk_2all1 = []; fcmk_constall1 = []; cmk_slope1 = [];
for i = 1:3
coeffs = polyfit(AUC2traj([1,3]),cmk_2all1(:,i)',1);
fcmk_2all1(:,i) = polyval(coeffs, AUC2traj([1,3]));
coeffs = polyfit(AUCdurtraj([1,3]),cmk_constall1(:,i)',1);
fcmk_constall1(:,i) = polyval(coeffs, AUCdurtraj([1,3]));
cmk_slope1(i) = ((fcmk_2all1(2,i)-fcmk_2all1(1,i))./(AUC2traj(3) - AUC2traj(1)))./((fcmk_constall1(2,i)-fcmk_constall1(1,i))./(AUCdurtraj(3)-AUCdurtraj(1)))
end


%% plot point spread plot
percDiff = [[gyp_slope,gyp_slope1];[cmk_slope,nan(1,3)];[yps_slope,nan(1,3)];...
    [ena_slope,nan(1,3)];[put_slope(1),nan(1,1),put_slope(3),nan(1,3)];[mep_slope,nan(1,3)]]'; 
[mpercDiff,indx] = sort(nanmean(percDiff),'descend');
spercDiff = percDiff(:,indx);
%genes1 = {'pGyp7','pCmk2','pYps1','pPun1','pEna1','pPut1','pMep1'};
genes1 = {'pGyp7','pCmk2','pYps1','pEna1','pPut1','pMep1'};
sgenes1 = genes1(indx);
% plot swarm + box plot (crz1 only)
figure(9); boxplot(spercDiff); hold on; 
xticks([1:6]); xticklabels(sgenes1); set(gca,'xticklabelrotation',90);
ylabel('Slope ratio'); title('Sensitivity to frequency'); box off
%plotSpread(spercDiff,'distributionColors','g'); 
%set(findall(gca,'color','g'),'markerSize',10); axis tight; box off;
end