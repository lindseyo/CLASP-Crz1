% GRID kd vs slowness
function [flag] = Figure4B()
flag = 1;

% define x axis
load('Input_idealizedLight.mat','lightInputs')
kfact = 1000;
orig_lightInputs = lightInputs;
m2_xaxis=[orig_lightInputs(1).auc,orig_lightInputs(2).auc./kfact,orig_lightInputs(3).auc./kfact,orig_lightInputs(4).auc./kfact,orig_lightInputs(5).auc./kfact];
single_xaxis=[orig_lightInputs(1).auc,orig_lightInputs(6).auc./kfact,orig_lightInputs(7).auc./kfact,orig_lightInputs(8).auc./kfact,orig_lightInputs(9).auc./kfact];

%% Figure 4B
%% g1 = 0.025
load('20191001_model_valsGen_GRIDg1pt025.mat')
g1 = 0.025; 
% PLOTTING GRID -- code for plotting slope ratio for PROTEIN
ratioSlope = [];
for i = 1:20
val = 60;
JJ = i;
mod2_S = [allsGen{JJ,1}{1}.prot(end-val),allsGen{JJ,2}{1}.prot(end-val),allsGen{JJ,3}{1}.prot(end-val),allsGen{JJ,4}{1}.prot(end-val),allsGen{JJ,5}{1}.prot(end-val)];
mod40_S = [allsGen{JJ,1}{1}.prot(end-val),allsGen{JJ,6}{1}.prot(end-val),allsGen{JJ,7}{1}.prot(end-val),allsGen{JJ,8}{1}.prot(end-val),allsGen{JJ,9}{1}.prot(end-val)];
a = 0; b = 1;
mod40S = (mod40_S-min(mod40_S))./(max(mod40_S)-min(mod40_S));
mod2S = (mod2_S-min(mod40_S))./(max(mod40_S)-min(mod40_S));
ratioSlope(i) = ((mod2S(5) - mod2S(1))./(m2_xaxis(5)-m2_xaxis(1)))./((mod40S(5)-mod40S(1))./(single_xaxis(5)-single_xaxis(1)));
end
kdVSslowGRID = [ratioSlope(1:4);ratioSlope(5:8);ratioSlope(9:12);ratioSlope(13:16);ratioSlope(17:20)];

figure(913); imagesc(kdVSslowGRID'); yticklabels({'1 x','1/5 x','1/10 x','1/20 x'}); xticklabels({'0.5','2.3','10','20','46'}); title(['Slope ratio with decreasing kon, koff, & kd']);
colorbar; xlabel('kd'); caxis([1.1 2.1]);

end