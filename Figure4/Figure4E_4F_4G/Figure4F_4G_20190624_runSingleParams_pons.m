function [flag]= Figure4F_4G_20190624_runSingleParams_pons()
flag = 1;

%% Script plots out the 3-4 regions for FIGURE 5
load('Input_idealizedLight.mat', 'lightInputs')
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
lightInput{16} = orig_lightInputs(16).lightInput./kfact; % rpl
lightInput{17} = orig_lightInputs(17).lightInput./kfact; % tef
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

%% parameters for pYPS1
load('20190616_params_kineticMod_pYPS1_g1pt06_rerunFIG4','paramset')
param10000 = paramset(:,1:9000);

% get region1: kon/koff 
vals = param10000(1,:)./param10000(2,:);
[v,i]=sort(vals);

% make GRID kd vs slowness of kon & koff
paramset1 = [.1*20; .23*20 ;param10000(3:end,i(1))]; 
paramset3 = [.1*20/10; .23*20/10 ;param10000(3:end,i(1))]; 

% get all parameters
paramset = [paramset1, paramset3];

%%
% RUN through ODE
numParams = size(paramset,2); % number of params looping through
numLight = length(lightInput);
tspan = linspace(0,360,36e1);
% storage
%allsGen=[];
model_valsGen=[];
for mm = 1:numParams
    % display
    display(strcat('current loop: ', num2str(mm)))

    % run ode for gene expression
    % parameters
    M = []; 
    M(1) = paramset(1,mm); % kact
    M(2) = paramset(2,mm); % kinact
    M(3) = paramset(3,mm);% b1
    M(4) = paramset(4,mm);% g1
    M(5) = paramset(5,mm);% b2
    M(6) = paramset(6,mm); % g2
    M(7) = paramset(7,mm); % b0
    
    for nn = 1:numLight
    %%%%%%%%%%%%%%%%%%DUR%%%%%%%%%%%%%%%%%%%%%%%
    % light input (10min)
    A = struct;
    A.times = lightInputTimes{nn};
    A.nucloc = lightInput{nn};
    % run ODE
    alls = [];
    [alls] = run_Ode(M,A,tspan);
    %[alls] = run_Ode_linear(M,A,tspan);
    allsGen{mm,nn} = alls;
    % compare MODEL and EXP trace
    prot_end = [];
    N = length(allsGen{mm,nn});
    for k = 1:N
        %display(strcat('k=',num2str(k)));
        s = allsGen{mm,nn}{k};
        prot_end(k) = s.prot(end-3); % yaxis
    end
    model_valsGen{mm,nn} = prot_end;
    end

end

model_valsGen = cell2mat(model_valsGen);

%%
val = 60; % 5hrs

JJ = 1;
mod2_S = [allsGen{JJ,1}{1}.prot(end-val),allsGen{JJ,2}{1}.prot(end-val),allsGen{JJ,3}{1}.prot(end-val),allsGen{JJ,4}{1}.prot(end-val),allsGen{JJ,5}{1}.prot(end-val)];
mod40_S = [allsGen{JJ,1}{1}.prot(end-val),allsGen{JJ,6}{1}.prot(end-val),allsGen{JJ,7}{1}.prot(end-val),allsGen{JJ,8}{1}.prot(end-val),allsGen{JJ,9}{1}.prot(end-val)];

mod40S = (mod40_S-min(mod40_S))./(max(mod40_S)-min(mod40_S));
mod2S = (mod2_S-min(mod40_S))./(max(mod40_S)-min(mod40_S));

JJ = 2; 
mod2_S2 = [allsGen{JJ,1}{1}.prot(end-val),allsGen{JJ,2}{1}.prot(end-val),allsGen{JJ,3}{1}.prot(end-val),allsGen{JJ,4}{1}.prot(end-val),allsGen{JJ,5}{1}.prot(end-val)];
mod40_S2 = [allsGen{JJ,1}{1}.prot(end-val),allsGen{JJ,6}{1}.prot(end-val),allsGen{JJ,7}{1}.prot(end-val),allsGen{JJ,8}{1}.prot(end-val),allsGen{JJ,9}{1}.prot(end-val)];

mod40S2 = (mod40_S2-min(mod40_S2))./(max(mod40_S2)-min(mod40_S2));
mod2S2 = (mod2_S2-min(mod40_S2))./(max(mod40_S2)-min(mod40_S2));

% define Xaxis
m2_xaxis = [orig_lightInputs(1).auc,orig_lightInputs(2).auc/kfact,orig_lightInputs(3).auc/kfact,orig_lightInputs(4).auc/kfact,orig_lightInputs(5).auc/kfact];
single_xaxis=[orig_lightInputs(1).auc,orig_lightInputs(6).auc./kfact,orig_lightInputs(7).auc./kfact,orig_lightInputs(8).auc./kfact,orig_lightInputs(9).auc./kfact];

% %  plot protein
% figure(19);  
% % fast parameter set
% plot(m2_xaxis, mod2S,'r.-','markersize',20); hold on; 
% plot(single_xaxis,mod40S,'b.-','markersize',20);
% % slow parameter set
% plot(m2_xaxis, mod2S2,'r--','markersize',20); hold on; 
% plot(single_xaxis,mod40S2,'b--','markersize',20);
% 
% axis tight; box off; ylabel('Normalized Protein'); xlabel('TF nuclear occupancy');
% title('kd 2.3 O-O'); 
% % slope ratio (it's here!)
% ratioSlope = ((mod2S(5) - mod2S(1))./(m2_xaxis(5)-m2_xaxis(1)))./((mod40S(5)-mod40S(1))./(single_xaxis(5)-single_xaxis(1)))
% themParams = paramset;
% konkoff = paramset(2,JJ)./paramset(1,JJ);
% text(10,.75,[num2str(themParams(:,JJ))])
% text(200,.2,['koff/kon = ',num2str(konkoff)])
% text(200,.3,['Np/Nc = ', num2str(ratioSlope)])

%% Figure 4F plotting Pon and TF [2m20m and 40] 20191001
JJ = 2;
mod2 = allsGen{JJ,3}{1}.pon; % 2 or 3
mod40 = allsGen{JJ,13}{1}.pon;

figure(2);  title(['kinMod: region 2']);
yyaxis left; 
plot(allsGen{JJ,3}{1}.t,mod2,'r-'); hold on; xlim([0 10])
plot(allsGen{JJ,13}{1}.t,mod40,'b-'); ylabel('pon'); xlim([0 50]); ylim([0 1]); 
yyaxis right;  
plot(allsGen{JJ,3}{1}.nuclocT,allsGen{JJ,3}{1}.nucloc,'r-'); hold on; 
plot(allsGen{JJ,13}{1}.nuclocT,allsGen{JJ,13}{1}.nucloc,'b-');  ylabel('light');
xlabel('time(min)'); box off
legend('pulsed','const');

%% Figure 4G plotting Pon and TF [2m20m and 40] 20191001
JJ = 1;
mod2 = allsGen{JJ,3}{1}.pon; % 2 or 3
mod40 = allsGen{JJ,13}{1}.pon;

figure(1);  title(['kinMod: region 1']);
yyaxis left; 
plot(allsGen{JJ,3}{1}.t,mod2,'r-'); hold on; xlim([0 10])
plot(allsGen{JJ,13}{1}.t,mod40,'b-'); ylabel('pon'); xlim([0 50]); ylim([0 1]); 
yyaxis right;  
plot(allsGen{JJ,3}{1}.nuclocT,allsGen{JJ,3}{1}.nucloc,'r-'); hold on; 
plot(allsGen{JJ,13}{1}.nuclocT,allsGen{JJ,13}{1}.nucloc,'b-');  ylabel('light');
xlabel('time(min)'); box off
legend('pulsed','const');

