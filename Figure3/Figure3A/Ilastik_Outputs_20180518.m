% Script - 
% what it does: 1) loads in mean intensity value data 2) calculates mean
% and stderr 3) subtracts background 4) saves - sample order, means,
% stderrs, distributions 5) plots outputs

%function [adh_40_mVals,adh_40_sVals,adh_2_mVals,adh_2_sVals] = Ilastik_Outputs_20180518()
function [flag] = Ilastik_Outputs_20180518()
%%%%%%%%%%%%%%%%%%%%%%%%%
% 20180516
% adh1 (5/16) p1
basdir1 ='';
folder1 = 'adh1_1_csv\';
listOfFiles = [1:2:199]; listOfFiles = num2str(listOfFiles); listOfFiles = strsplit(listOfFiles,' ');
colnum = 7;
wellnum = '1';
% basdir2 = '\\elsamad.ucsf.edu\Data\Instrumentation\microscope\SYC\20180517_comp2mVS4m_SYC526_yps1strain\';
% folder2 = 'adh1\';
basdir2 = '';
folder2 = '';

[adh_2_15_mIntensity,adh_2_15_mVals,adh_2_15_sVals,adh_2_15_mbVals,adh_2_15_numVals] = findMeanStderrs(basdir1,folder1,basdir2,folder2,listOfFiles,colnum,wellnum);
% adh1 (5/16) p2
basdir1 ='';
folder1 = 'adh1_1_csv\';
listOfFiles = [1:2:199]; listOfFiles = num2str(listOfFiles); listOfFiles = strsplit(listOfFiles,' ');
colnum = 7;
wellnum = '2';
basdir2 = '';
folder2 = '';

[adh_2_15_p2_mIntensity,adh_2_15_p2_mVals,adh_2_15_p2_sVals,adh_2_15_p2_mbVals,adh_2_15_p2_numVals] = findMeanStderrs(basdir1,folder1,basdir2,folder2,listOfFiles,colnum,wellnum);
% adh1 (5/16) p3 (40min)
basdir1 ='';
folder1 = 'adh1_1_csv\';
listOfFiles = [1:2:199]; listOfFiles = num2str(listOfFiles); listOfFiles = strsplit(listOfFiles,' ');
colnum = 7;
wellnum = '3';
basdir2 = '';
folder2 = '';

[adh_40_15_mIntensity,adh_40_15_mVals,adh_40_15_sVals,adh_40_15_mbVals,adh_40_15_numVals] = findMeanStderrs(basdir1,folder1,basdir2,folder2,listOfFiles,colnum,wellnum);
% adh1 (5/16) p3 (40min)
basdir1 ='C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\LANS OPTOGENETICS PROJECT\SYC Data\20180518_3uscopy_datasets_robust nuclearQuant_Adh1_rpltef\20180516_irfpMasks\';
folder1 = 'adh1_1_csv\';
listOfFiles = [1:2:199]; listOfFiles = num2str(listOfFiles); listOfFiles = strsplit(listOfFiles,' ');
colnum = 7;
wellnum = '4';
basdir2 = '';
folder2 = '';

[adh_40_15_p2_mIntensity,adh_40_15_p2_mVals,adh_40_15_p2_sVals,adh_40_15_p2_mbVals,adh_40_15_p2_numVals] = findMeanStderrs(basdir1,folder1,basdir2,folder2,listOfFiles,colnum,wellnum);
% plot!
% % 2&10min for each promoter together
% adh1 5/16 (2 and 40min)
% figure(11); 
% time =[0.5:0.5:numel([1:2:199])./2]; shadedErrorBar(time,adh_2_15_mbVals,adh_2_15_sVals,{'color','r'},1); hold on;
% title('adh1-20180516 p1'); xlabel('time(min)'); ylabel('Fluor(bgsub)'); axis tight;
% time =[0.5:0.5:numel([1:2:199])./2]; shadedErrorBar(time,adh_2_15_p2_mbVals,adh_2_15_p2_sVals,{'color','r'},1); hold on;
% title('adh1-20180516 p2'); xlabel('time(min)'); ylabel('Fluor(bgsub)'); axis tight;
% time =[0.5:0.5:numel([1:2:199])./2]; shadedErrorBar(time,adh_40_15_mbVals,adh_40_15_sVals,{'color','b'},1); hold on;
% title('adh1-20180516 p3'); xlabel('time(min)'); ylabel('Fluor(bgsub)'); axis tight; 
% time =[0.5:0.5:numel([1:2:199])./2]; shadedErrorBar(time,adh_40_15_p2_mbVals,adh_40_15_p2_sVals,{'color','b'},1); hold on;
% title('adh1-20180516 p4'); xlabel('time(min)'); ylabel('Fluor(bgsub)'); axis tight;

% correct for the slant downwards
first=adh_40_15_mbVals(1);
last=min(adh_40_15_mbVals);
line1 = first.*ones(1,100);
coeffLine2 = polyfit([0, 100], [first, last], 1); a = coeffLine2(1); b = coeffLine2(2);
line2 = a.*[1:100] + b; %figure(14); plot(line1,'r'); hold on; plot(line2,'b');
diffs = line1 - line2; %figure(15); plot(diffs);
%figure(15); plot(adh_40_15_mbVals); hold on; plot(adh_40_15_mbVals+diffs,'r');
adh_40_15_mbVals = adh_40_15_mbVals+diffs;
first=adh_2_15_mbVals(1);
last=min(adh_2_15_mbVals);
line1= first.*ones(1,100);
coeffLine2 = polyfit([0, 100], [first, last], 1); a = coeffLine2(1); b = coeffLine2(2);
line2 = a.*[1:100] + b; %figure(16); plot(line1,'r'); hold on; plot(line2,'b');
diffs = line1 - line2; %figure(17); plot(diffs);
%figure(16); plot(adh_2_15_mbVals); hold on; plot(adh_2_15_mbVals+diffs,'r');
adh_2_15_mbVals = adh_2_15_mbVals+diffs;

% save the 20150516 data
adh2_p1_trace20180516 = adh_2_15_mbVals;
adh40_p3_trace20180516 = adh_40_15_mbVals;
adh2S_p1_trace20180516 = adh_2_15_sVals;
adh40S_p3_trace20180516 = adh_40_15_sVals
adh2_p1_FC_trace20180516 = adh_2_15_mbVals./adh_2_15_mbVals(1);
adh40_p3_FC_trace20180516 = adh_40_15_mbVals./adh_40_15_mbVals(1);
%basedir = 'C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\LANS OPTOGENETICS PROJECT\SYC Data\20180518_3uscopy_datasets_robust nuclearQuant_Adh1_rpltef\';
%save(strcat(basedir,'adhTFnuctraces_20180516.mat'),'adh2_p1_trace20180516','adh40_p3_trace20180516','adh2_p1_FC_trace20180516','adh40_p3_FC_trace20180516','adh2S_p1_trace20180516','adh40S_p3_trace20180516');

%% FOR FIGURE 3 WE WILL PLOT THE NUCLEAR LOCALIZATION WITH 95% CI WITH THE 20180516 SET OF DATA
figure(21); subplot(2,1,1)
% 95% CI calc
ci2 = 1.96.*((adh_2_15_sVals)./sqrt(adh_2_15_numVals));
size(ci2)
time =[0.5:0.5:numel([1:2:199])./2]; shadedErrorBar(time,adh_2_15_mbVals./1000,ci2./1000,{'color','r'},1); hold on;
%title('adh1-20180516 p1'); ylabel('Fluor(bgsub)'); xlabel('time(min)'); axis tight; 
%time =[0.5:0.5:numel([1:2:199])./2]; shadedErrorBar(time,adh_2_15_p2_mbVals,adh_2_15_p2_sVals,{'color','r'},1); hold on;
%title('adh1-20180516 p2'); xlabel('time(min)'); ylabel('Fluor(bgsub)'); axis tight;
axis tight; xlim([0 45]); box off; 
subplot(2,1,2)
% 95% CI calc
ci40 = 1.96.*((adh_40_15_sVals)./sqrt(adh_40_15_numVals));
size(ci40)
time =[0.5:0.5:numel([1:2:199])./2]; shadedErrorBar(time, adh_40_15_mbVals./1000,ci40./1000,{'color','b'},1); hold on;
%title('adh1-20180516 p3'); ylabel('Fluor(bgsub)'); v
%time =[0.5:0.5:numel([1:2:199])./2]; shadedErrorBar(time,adh_40_15_p2_mbVals,adh_40_15_p2_sVals,{'color','b'},1); hold on;
%title('adh1-20180516 p4'); xlabel('time(min)'); ylabel('Fluor(bgsub)'); axis tight;
axis tight; xlim([0 45]); box off; 

flag = 1;
end

function [mIntensity,mVals,sVals,mbVals,numVals] = findMeanStderrs(basdir1,folder1,basdir2,folder2,listOfFiles,colnum,wellnum)
% means and stderrs
mVals = []; sVals = [];
for i = 1:numel(listOfFiles);
    i;
    filename = strcat(basdir1,folder1,'exported_data-RFP-wheel_p',wellnum,'_t',listOfFiles{i},'_table.csv');
    A=csvread(filename,1,colnum);
    mIntensity{i} = A(:,1);
    mVals(i) = mean(A(:,1));
    sVals(i) = std(A(:,1))./sqrt(numel(A(:,1)));
    numVals = numel(A(:,1));
end
% subtract background
% imfile = strcat(basdir2,folder2,'RFP-wheel_p1_t1.tiff');
% im1 = imread(imfile);
% mBG = quantile(double(im1(:)),.95);
% mbVals = []; mbVals = mVals - mBG;
mbVals = [];
for j = 1:numel(listOfFiles);
    imfile = strcat(basdir2,folder2,'RFP-wheel_p',wellnum,'_t',listOfFiles{j},'.tiff');
    mBG = mean(double(imfile(:)));
    mbVals(j) = mVals(j) - mBG;
end
end