function [flag] = FigureS4D(); 
flag = 1;
% plot all
% sample order: 
% 1 = 526 pPun1 0min light
% 2 = 526 pPun1 30min light
% 3 = 524 pPun1 0min light
% 4 = 524 pPun1 30min light
% 5 = 526 pPun1 5min light
% 6 = 526 pPun1 10min light
% 7 = 526 pPun1 20min light

load('SYC20171128_datafile.csv');

traj0 = SYC20171128_datafile(1:7,:);
time0 = SYC20171128_datafile(8:14,:)./360000 - 35./60;

% compare 526 and 524
traj1 = traj0(1:4,3:end); traj1 = 2.^traj1; %unlog
time1 = time0(1:4,3:end);

% normalized to 0 light
figure(1); plot(time1(1,:)',traj1(2,:)./traj1(1,:), 'm.-', 'linewidth',3, 'markersize',20); hold on;
plot(time1(3,:)',traj1(4,:)./traj1(3,:), 'b.-', 'linewidth',3, 'markersize',20)
xlabel('time (hrs)'); ylabel('FoldChange (FITC)'); title('Pun1 - WT vs 19A 30min');
legend('19A 30min','WT 30min'); axis tight; box off;



end