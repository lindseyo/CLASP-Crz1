%% Figure 4D - PLOT time courses with different KD!
function [flag] = runSingleParams20190624_pons_extended()
flag = 1; 

% define light input
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


% rest of function

kd = [2.3, 4.6*10];
TF2 = lightInput{3};
TFC = lightInput{13};

ponV2 = []; iponV2 = [];
ponVC = []; iponVC = [];
a = 0; b = 1;
for hh = 1:length(kd)
ponV2(hh,:) = TF2./(TF2 + kd(hh));
iponV2(hh,:) = (b-a).*(ponV2(hh,:)-min(ponV2(hh,:)))./(max(ponV2(hh,:)) - min(ponV2(hh,:))) + a;
ponVC(hh,:) = TFC./(TFC + kd(hh));
iponVC(hh,:) = (b-a).*(ponVC(hh,:)-min(ponVC(hh,:)))./(max(ponVC(hh,:)) - min(ponVC(hh,:))) + a;
end

% new way of plotting 20191010
% figure(1); plot(lightInputTimes{1}, ponV2(1,:),'r'); axis tight; box off; ...
%     xlabel('time (min)'); ylabel('pon_ss'); title('pulsed'); xlim([0 40]);% xlim([0 50]); %
% figure(2); plot(lightInputTimes{1}, ponVC(1,:),'b'); axis tight; box off; ...
%     xlabel('time (min)'); ylabel('pon_ss'); title('continuous'); xlim([0 9]); %xlim([0 50]); %

figure(3); plot(lightInputTimes{1}, ponV2(2,:),'r'); axis tight; box off; ...
   xlabel('time (min)'); ylabel('pon_ss'); title('pulsed'); xlim([0 40]); %xlim([0 50]); %
figure(4); plot(lightInputTimes{1}, ponVC(2,:),'b'); axis tight; box off; ...
   xlabel('time (min)'); ylabel('pon_ss'); title('continuous'); xlim([0 9]); %xlim([0 50]); %

figure(5); plot(lightInputTimes{1}, lightInput{3},'r'); axis tight; box off; ... % light input
    xlabel('time (min)'); ylabel('TF'); title('pulsed'); xlim([0 50]); %xlim([0 40]);...
figure(6); plot(lightInputTimes{1}, lightInput{13},'b'); axis tight; box off; ...
    xlabel('time (min)'); ylabel('TF'); title('continuous'); xlim([0 50]); %xlim([0 9]);
end