function [flag] = FigureS7D()
flag = 1;
%% plot time course 
load('20180930_idealizedlightInputs_longer_preSS.mat','lightInputs')
orig_lightInputs = lightInputs;
kfact = 1000;
lightInput=[];
lightInput{1} = orig_lightInputs(1).lightInput./kfact; 
lightInput{2} = orig_lightInputs(2).lightInput./kfact;
lightInput{3} = orig_lightInputs(3).lightInput./kfact;
lightInput{4} = orig_lightInputs(4).lightInput./kfact;
lightInput{5} = orig_lightInputs(5).lightInput./kfact;
lightInput{6} = orig_lightInputs(6).lightInput./kfact;
lightInput{7} = orig_lightInputs(7).lightInput./kfact;
lightInput{8} = orig_lightInputs(8).lightInput./kfact;
lightInput{9} = orig_lightInputs(9).lightInput./kfact; 
lightInput{10} = orig_lightInputs(10).lightInput./kfact;
lightInput{11} = orig_lightInputs(11).lightInput./kfact; 
lightInput{12} = orig_lightInputs(12).lightInput./kfact; % rpl
lightInput{13} = orig_lightInputs(13).lightInput./kfact; % tef
lightInput{14} = orig_lightInputs(14).lightInput./kfact;
lightInput{15} = orig_lightInputs(15).lightInput./kfact; 
lightInput{16} = orig_lightInputs(16).lightInput./kfact; 
lightInput{17} = orig_lightInputs(17).lightInput./kfact; 
lightInput{18} = orig_lightInputs(18).lightInput./kfact; 

lightInputTimes = [];
lightInputTimes{1} = orig_lightInputs(1).times;
lightInputTimes{2} = orig_lightInputs(1).times;
lightInputTimes{3} = orig_lightInputs(1).times;
lightInputTimes{4} = orig_lightInputs(1).times;
lightInputTimes{5} = orig_lightInputs(1).times;
lightInputTimes{6} = orig_lightInputs(1).times;
lightInputTimes{7} = orig_lightInputs(1).times;
lightInputTimes{8} = orig_lightInputs(1).times;
lightInputTimes{9} = orig_lightInputs(1).times;
lightInputTimes{10} = orig_lightInputs(1).times;
lightInputTimes{11} = orig_lightInputs(1).times;
lightInputTimes{12} = orig_lightInputs(1).times;
lightInputTimes{13} = orig_lightInputs(1).times;
lightInputTimes{14} = orig_lightInputs(1).times;
lightInputTimes{15} = orig_lightInputs(1).times;
lightInputTimes{16} = orig_lightInputs(1).times;
lightInputTimes{17} = orig_lightInputs(1).times;
lightInputTimes{18} = orig_lightInputs(1).times;

%% parameters
% % region 2 (increase ron) [Region 3]
% kon = 0.0037*2; 
%     koff = 0.664; 
%     b1 = 7.93; g1 = 0.05; b2 = 0.06; g2 = 0.0083; b0 = 0.00316;...
%     ron = 0.2255*8.8692;  % 0.25; 
%     roff = 0.908; %0.25*(.908/.2255);
%     thr = 0.5;

    
% % region 3 (increase kon) [Region 4]
kon = (.0037*2)*13.514;%(.0037*2)*270.2703; % 0.0037*2;%
    koff = (.664/5);%0.664/270.2703; 
    b1 = 7.93; g1 = 0.05; b2 = 0.06; g2 = 0.0083; b0 = 0.00316;...
    ron = 0.2255;  % 0.25; 
    roff = 0.908; %0.25*(.908/.2255);
    thr = 0.5;
    
paramset = [kon;koff;b1;g1;b2;g2;b0;ron;roff;thr];
%%
% RUN through ODE
numParams = size(paramset,2); % number of params looping through
numLight = length(lightInput);
tspan = linspace(0,360,36e1);
%tspan = linspace(0,1200,120e1);

% storage
allsGen=[];
allsGenL = [];
model_valsGen=[];
for mm = 1:numParams
    % display
    display(strcat('current loop: ', num2str(mm)))

    % run ode for gene expression
    % parameters
    M = []; 
    M(1) = paramset(1,mm); % kon
    M(2) = paramset(2,mm); % koff
    M(3) = paramset(3,mm);% b1
    M(4) = paramset(4,mm);% g1
    M(5) = paramset(5,mm);% b2
    M(6) = paramset(6,mm); % g2
    M(7) = paramset(7,mm); % b0
    M(8) = paramset(8,mm);% ron
    M(9) = paramset(9,mm);% roff
    M(10) = paramset(10,mm);% thron
    
    % init Vals - xinitvim 
    SSvals = [];
    SSvals(1).p0SS = 1;
    SSvals(1).poffSS = 0;
    SSvals(1).ponSS = 0;
    SSvals(1).mrnaSS = 0;
    SSvals(1).protSS = 0;
    
    for nn = 1:numLight
    %%%%%%%%%%%%%%%%%%DUR%%%%%%%%%%%%%%%%%%%%%%%
    % light input (10min)
    A = struct;
    A.times = lightInputTimes{nn};
    A.nucloc = lightInput{nn};
    % run ODE
    alls = [];
    [alls] = run_Ode_2state_mult(M,A,tspan,SSvals);
    allsGen{mm,nn} = alls;
    % compare MODEL and EXP trace
    prot_end = [];
    N = length(alls);
    for k = 1:N
        %display(strcat('k=',num2str(k)));
        s = alls{k};
        prot_end(k) = s.all(end-3,end); % yaxis
    end
    model_valsGen{mm,nn} = prot_end;
    allsGen{mm,nn} = s.all(:,end);
    allsGenL{mm,nn} = s.all;
    end

end

model_valsGen = cell2mat(model_valsGen);
%% scale to experimental data
% load data (end point)
load('20181127_gyp_pulsedNconst_data','constGE','sconstGE','pulsedGE','spulsedGE') % newer data set
GE2m = pulsedGE; GEdur = constGE; sGE2m = spulsedGE; sGEdur = sconstGE;
% define Xaxis
load('Input_idealizedLight.mat', 'lightInputs')
kfact = 1000;
orig_lightInputs = lightInputs;
m2_xaxis=[orig_lightInputs(1).auc,orig_lightInputs(2).auc./kfact,orig_lightInputs(3).auc./kfact,orig_lightInputs(4).auc./kfact,orig_lightInputs(5).auc./kfact];
single_xaxis=[orig_lightInputs(1).auc,orig_lightInputs(6).auc./kfact,orig_lightInputs(7).auc./kfact,orig_lightInputs(8).auc./kfact,orig_lightInputs(9).auc./kfact];
% Get range
p = polyfit(single_xaxis,GEdur,1);
xdur = polyval(p,single_xaxis);
% errorbar(single_xaxis,GEdur,sGEdur,'.b','markersize',20,'linewidth',2);...
%     plot(single_xaxis,xdur,'b'); axis tight; title('best fit line to data');
%%
% for the freshly run
val = 60;
a = min(xdur); b = max(xdur);
JJ = 1;
mod2 = [allsGen{JJ,1}(end-val),allsGen{JJ,2}(end-val),allsGen{JJ,3}(end-val),allsGen{JJ,4}(end-val),allsGen{JJ,5}(end-val)];
mod40 = [allsGen{JJ,1}(end-val),allsGen{JJ,6}(end-val),allsGen{JJ,7}(end-val),allsGen{JJ,8}(end-val),allsGen{JJ,9}(end-val)];
    ix = (b-a).*(mod2-min(mod40))./(max(mod40)-min(mod40))+a; % scale to experimental data
    iy = (b-a).*(mod40-min(mod40))./(max(mod40)-min(mod40))+a;
% define Xaxis
m2_xaxis = [orig_lightInputs(1).auc,orig_lightInputs(2).auc/kfact,orig_lightInputs(3).auc/kfact,orig_lightInputs(4).auc/kfact,orig_lightInputs(5).auc/kfact];
%m2_xaxis=[orig_lightInputs(1).auc,orig_lightInputs(2).auc,orig_lightInputs(3).auc,orig_lightInputs(4).auc,orig_lightInputs(5).auc];
single_xaxis=[orig_lightInputs(1).auc,orig_lightInputs(6).auc./kfact,orig_lightInputs(7).auc./kfact,orig_lightInputs(8).auc./kfact,orig_lightInputs(9).auc./kfact];
%  plot protein
figure(1); 
plot(m2_xaxis, ix,'r.-','markersize',20); hold on; plot(single_xaxis,iy,'b.-','markersize',20);% model
axis tight; box off; ylabel('Protein'); xlabel('TF nuclear occupancy'); title(['Region 1: acutal data fit']);
% slope ratio
ratioSlope = ((mod2(5) - mod2(1))./(m2_xaxis(5)-m2_xaxis(1)))./((mod40(5)-mod40(1))./(single_xaxis(5)-single_xaxis(1)))
themParams = paramset;
konkoff = log10(paramset(1)./paramset(2));
ronroff = log10(paramset(8)./paramset(9));
text(10,.3,[num2str(themParams)])
text(200,.23,['log10konkoff = ',num2str(konkoff)])
text(200,.22,['log10ronroff = ',num2str(ronroff)])

%% try plotting the promoter p0, poff, pon
JJ=1;
A=allsGenL(JJ,:);
figure(13); %subplot(2,1,1); 

yyaxis left; plot([0:300-.5]-4,A{JJ,3}(1:end-60,3),'r-'); ylabel('pon'); xlabel('time(min)'); title('pon region1 roffTHR'); hold on; %title('intuitive Pon');
 plot([0:300-.5]+1,A{JJ,7}(1:end-60,3),'b-'); axis tight; xlim([100 150]); box off; %ylim([0 .9]);%ylim([0 .00095]);
  yyaxis right; plot(lightInputTimes{3}-5,lightInput{3},'r-'); hold on; % gotta plot that NUC LOC
 plot(lightInputTimes{7},lightInput{7},'b-'); 

 figure(14);
 yyaxis left; plot([0:300-.5]-4,A{JJ,3}(1:end-60,1),'r-'); ylabel('p0');  xlabel('time(min)'); title('p0 region1 roffTHR');hold on; %title('intuitive P0');
 plot([0:300-.5]+1,A{JJ,7}(1:end-60,1),'b-'); axis tight; xlim([100 150]); box off; %ylim([0 1]);
 yyaxis right; plot(lightInputTimes{3}-5,lightInput{3},'r-'); hold on; % gotta plot that NUC LOC
 plot(lightInputTimes{7},lightInput{7},'b-'); 

%  figure(15);
%  yyaxis left; plot([0:300-.5]-4,A{JJ,3}(1:end-60,2),'r-'); ylabel('poff'); xlabel('time(min)'); title('poff region1 roffTHR'); hold on; %title('intuitive Poff');
%  plot([0:300-.5]+1,A{JJ,7}(1:end-60,2),'b-'); axis tight; xlim([100 150]); box off; %ylim([0 1]);
%  yyaxis right; plot(lightInputTimes{3}-5,lightInput{3},'r-'); hold on; % gotta plot that NUC LOC
%  plot(lightInputTimes{7},lightInput{7},'b-'); 