%% Script plots out the 3-4 regions for FIGURE 5
function [flag] = FigureS5E_kinMod_nofit_singleParamRUns()
flag = 1;
%%
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
% make GRID kd vs slowness of kon & koff
paramset1 = [1; 10; 0.1; 0.05; 0.1; 0.0083; 0.000001]; % 1
paramset2 = [paramset1(1,:)./2; paramset1(2:end,:)]; % 1/2
paramset3 = [paramset1(1,:)./4 ;paramset1(2:end,:)]; % 1/4
paramset4 = [paramset1(1,:)./8; paramset1(2:end,:)]; % 1/8
paramset5 = [paramset1(1,:)./16; paramset1(2:end,:)]; % 1/16
paramset6 = [paramset1(1,:)./32; paramset1(2:end,:)]; %1/32

% get all parameters
paramset = [paramset1, paramset2, paramset3, paramset4, paramset5, paramset6];

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

%% Figure S5E
val = 60; % 5hrs

JJ = 1;
mod2_S = [allsGen{JJ,1}{1}.prot(end-val),allsGen{JJ,2}{1}.prot(end-val),allsGen{JJ,3}{1}.prot(end-val),allsGen{JJ,4}{1}.prot(end-val),allsGen{JJ,5}{1}.prot(end-val)];
mod40_S = [allsGen{JJ,1}{1}.prot(end-val),allsGen{JJ,6}{1}.prot(end-val),allsGen{JJ,7}{1}.prot(end-val),allsGen{JJ,8}{1}.prot(end-val),allsGen{JJ,9}{1}.prot(end-val)];

JJ = 2; 
mod2_S2 = [allsGen{JJ,1}{1}.prot(end-val),allsGen{JJ,2}{1}.prot(end-val),allsGen{JJ,3}{1}.prot(end-val),allsGen{JJ,4}{1}.prot(end-val),allsGen{JJ,5}{1}.prot(end-val)];
mod40_S2 = [allsGen{JJ,1}{1}.prot(end-val),allsGen{JJ,6}{1}.prot(end-val),allsGen{JJ,7}{1}.prot(end-val),allsGen{JJ,8}{1}.prot(end-val),allsGen{JJ,9}{1}.prot(end-val)];

JJ = 3; 
mod2_S3 = [allsGen{JJ,1}{1}.prot(end-val),allsGen{JJ,2}{1}.prot(end-val),allsGen{JJ,3}{1}.prot(end-val),allsGen{JJ,4}{1}.prot(end-val),allsGen{JJ,5}{1}.prot(end-val)];
mod40_S3 = [allsGen{JJ,1}{1}.prot(end-val),allsGen{JJ,6}{1}.prot(end-val),allsGen{JJ,7}{1}.prot(end-val),allsGen{JJ,8}{1}.prot(end-val),allsGen{JJ,9}{1}.prot(end-val)];

JJ = 4; 
mod2_S4 = [allsGen{JJ,1}{1}.prot(end-val),allsGen{JJ,2}{1}.prot(end-val),allsGen{JJ,3}{1}.prot(end-val),allsGen{JJ,4}{1}.prot(end-val),allsGen{JJ,5}{1}.prot(end-val)];
mod40_S4 = [allsGen{JJ,1}{1}.prot(end-val),allsGen{JJ,6}{1}.prot(end-val),allsGen{JJ,7}{1}.prot(end-val),allsGen{JJ,8}{1}.prot(end-val),allsGen{JJ,9}{1}.prot(end-val)];

JJ = 5; 
mod2_S5 = [allsGen{JJ,1}{1}.prot(end-val),allsGen{JJ,2}{1}.prot(end-val),allsGen{JJ,3}{1}.prot(end-val),allsGen{JJ,4}{1}.prot(end-val),allsGen{JJ,5}{1}.prot(end-val)];
mod40_S5 = [allsGen{JJ,1}{1}.prot(end-val),allsGen{JJ,6}{1}.prot(end-val),allsGen{JJ,7}{1}.prot(end-val),allsGen{JJ,8}{1}.prot(end-val),allsGen{JJ,9}{1}.prot(end-val)];

JJ = 6; 
mod2_S6 = [allsGen{JJ,1}{1}.prot(end-val),allsGen{JJ,2}{1}.prot(end-val),allsGen{JJ,3}{1}.prot(end-val),allsGen{JJ,4}{1}.prot(end-val),allsGen{JJ,5}{1}.prot(end-val)];
mod40_S6 = [allsGen{JJ,1}{1}.prot(end-val),allsGen{JJ,6}{1}.prot(end-val),allsGen{JJ,7}{1}.prot(end-val),allsGen{JJ,8}{1}.prot(end-val),allsGen{JJ,9}{1}.prot(end-val)];

% define Xaxis
m2_xaxis = [orig_lightInputs(1).auc,orig_lightInputs(2).auc/kfact,orig_lightInputs(3).auc/kfact,orig_lightInputs(4).auc/kfact,orig_lightInputs(5).auc/kfact];
single_xaxis=[orig_lightInputs(1).auc,orig_lightInputs(6).auc./kfact,orig_lightInputs(7).auc./kfact,orig_lightInputs(8).auc./kfact,orig_lightInputs(9).auc./kfact];
%% PLOT FigureS5E
%  plot protein
figure(19);  
% 1
plot(m2_xaxis, mod2_S,'r-','markersize',20); hold on; 
plot(single_xaxis,mod40_S,'b-','markersize',20);
% 1/2
plot(m2_xaxis, mod2_S2,'r-','markersize',20); hold on; 
plot(single_xaxis,mod40_S2,'b-','markersize',20);
% 1/4
plot(m2_xaxis, mod2_S3,'r-','markersize',20); hold on; 
plot(single_xaxis,mod40_S3,'b-','markersize',20);
% 1/8
plot(m2_xaxis, mod2_S4,'r-','markersize',20); hold on; 
plot(single_xaxis,mod40_S4,'b-','markersize',20);
% 1/16
plot(m2_xaxis, mod2_S5,'r-','markersize',20); hold on; 
plot(single_xaxis,mod40_S5,'b-','markersize',20);
% 1/32
plot(m2_xaxis, mod2_S6,'r-','markersize',20); hold on; 
plot(single_xaxis,mod40_S6,'b-','markersize',20);

axis tight; box off; ylabel('Normalized Protein'); xlabel('TF nuclear occupancy');
title('no fit'); 