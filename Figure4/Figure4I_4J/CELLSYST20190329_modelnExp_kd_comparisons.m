function [flag] = CELLSYST20190329_modelnExp_kd_comparisons()
flag = 1;
%%
% Fitting Experimental kds from pYPS1 and pCMK2
pYPS_expDose = [0.1887 0.5444	0.9636	1.1441];
stdpYPS_expDose = [0.0074	0.0749	0.0541	0.0579];
pCMK_expDose = [0.0747 0.3625	0.7733	1.3826];
stdpCMK_expDose = [0.0046	0.0076	0.0371	0.066];
TFvals = [0 0.97 2.64 6.03];

n_pYPS_expDose = (pYPS_expDose-min(pYPS_expDose))./(max(pYPS_expDose)-min(pYPS_expDose));
n_pCMK_expDose = (pCMK_expDose-min(pCMK_expDose))./(max(pCMK_expDose)-min(pCMK_expDose));

% % plot
% figure(1); plot(TFvals, n_pYPS_expDose,'b.-'); hold on; 
% plot(TFvals, n_pCMK_expDose,'m.-'); xlim([0 6]); ylim([0 1])

% equation form = 
%output = (C.*TF)./(TF + kd); % which is translated to: yf = @(b,x) (b(1).*x)./(x + b(2));
% pYPS1
x = TFvals;
y = n_pYPS_expDose;
yf = @(b,x) (b(1).*x)./(x + b(2));
B0 = [20, 40];
[Bm,normresm] = fminsearch(@(b) norm(y - yf(b,x)), B0);
SSE = sum((y - yf(Bm, x)).^2);

figure(1); 
%subplot(1,2,1); 
%plot(TFvals, n_pYPS_expDose,'k.','markersize',20); hold on;
errorbar(TFvals, n_pYPS_expDose, stdpYPS_expDose, 'k.', 'markersize',20, 'linewidth',2); hold on;
plot(TFvals, yf(Bm,x),'k');
axis tight; box off; title('pYPS1'); xlabel('TF amplitude'); ylabel('norm Protein');
text(2,0.3,{['exp kd = ', num2str(Bm(2))], ['SSE = ', num2str(SSE)],['model kd = ','2.3']})

% pCMK2
x = TFvals;
y = n_pCMK_expDose;
yf = @(b,x) (b(1).*x)./(x + b(2));
B0 = [20, 40];
[Bm,normresm] = fminsearch(@(b) norm(y - yf(b,x)), B0);
SSE = sum((y - yf(Bm, x)).^2);

%figure(3); 
%subplot(1,2,2); 
figure(2);
%plot(TFvals, n_pCMK_expDose,'k.','markersize',20); hold on;
errorbar(TFvals, n_pCMK_expDose, stdpCMK_expDose, 'k.', 'markersize',20, 'linewidth',2); hold on;
plot(TFvals, yf(Bm,x),'k'); 
axis tight; box off; title('pCMK2'); xlabel('TF amplitude'); ylabel('norm Protein');
text(2,0.3,{['exp kd = ', num2str(Bm(2))], ['SSE = ', num2str(SSE)], ['model kd = ','12.8']})
