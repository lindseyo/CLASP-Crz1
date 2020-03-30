function [flag]= Figure4C()
flag = 1;
%%
% 20190930 CODE to generate dose response for Figures 4C-D
TF = linspace(0.01,10,100);
kd1 = 2.3;
kd2 = 46;
% pon as a fn of TF 
% pon = TF/(TF + kd)
pon1 = TF./(TF + kd1);
pon2 = TF./(TF + kd2);
figure(1); plot(TF,pon1,'b'); hold on;
plot(TF,pon2,'k'); box off; ylabel('pon'); xlabel('TF');
legend('2.3','46');
title('diff kds - pon');
% protein as a fn of TF
% protein = (beta2.*beta1)./(gamma1.*gamma2).*TF./(TF+kd)
beta2 = 0.06; beta1 = 0.1; gamma1 = 0.04; gamma2 = 0.0083;
protein1 = (beta2.*beta1)./(gamma1.*gamma2).*TF./(TF +kd1);
protein2 = (beta2.*beta1)./(gamma1.*gamma2).*TF./(TF +kd2);
% figure(2); plot(TF,protein1,'b'); hold on;
% plot(TF,protein2,'k'); box off; ylabel('protein'); xlabel('TF');
% legend('2.3','46');
% title('diff kds - prot');
% quasi-linear plotting of pon as a fn of TF
TFmax = 2.6;
fact = 0.6; %0.6;
b1 = fact./TFmax;
pon1L = b1.*TF;
% figure(3); plot(TF,pon1,'b'); hold on;
% plot(TF,pon1L,'k'); box off; ylabel('pon'); xlabel('TF');
% legend('non','lin'); xlim([0 2.6]); 
% title('linear vs non - pon');
% quasi-linear plotting of protein as a fn of TF
fact = 0.033; %0.6;
b1 = fact./TFmax;
protein1L = (beta2.*beta1)./(gamma1.*gamma2).*b1.*TF;
% figure(4); plot(TF,pon1,'b'); hold on;
% plot(TF,protein1L,'k'); box off; ylabel('protein'); xlabel('TF');
% legend('non','lin'); xlim([0 2.6]); 
% title('linear vs non - prot');

