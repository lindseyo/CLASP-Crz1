function [flag] = Figure6H()
flag = 1;
% plotting of final data -- want 
load('20190518_resultsPlate1.mat','meanPlate1','sdPlate1','rowsPlate1','colsPlate1',...);
'gyp7_samples','gyp7_labels','gyp7_indices','yps1_samples','yps1_labels','yps1_indices');

% comparing promoter fusions pGYP7 and pYPS1
y_promFus = meanPlate1(1,[1:6]);
y_spromFus = sdPlate1(1,[1:6]);
y_xcoords = yps1_indices([3:5,8,6:7],1:2);

yps1_p1_x = [y_xcoords{1,1:2}(1):y_xcoords{1,1:2}(2)];
yps1_p1 = y_promFus{1,1}.*ones(1,length(yps1_p1_x));
yps1_p1_s = y_spromFus{1,1}.*ones(1,length(yps1_p1_x));

yps1_p2_x = [y_xcoords{2,1:2}(1):y_xcoords{2,1:2}(2)];
yps1_p2 = y_promFus{1,2}.*ones(1,length(yps1_p2_x));
yps1_p2_s = y_spromFus{1,2}.*ones(1,length(yps1_p2_x));

yps1_p3_x = [y_xcoords{3,1:2}(1):y_xcoords{3,1:2}(2)];
yps1_p3 = y_promFus{1,3}.*ones(1,length(yps1_p3_x));
yps1_p3_s = y_spromFus{1,3}.*ones(1,length(yps1_p3_x));

yps1_p7_x = [y_xcoords{4,1:2}(1):y_xcoords{4,1:2}(2)];
yps1_p7 = y_promFus{1,4}.*ones(1,length(yps1_p7_x));
yps1_p7_s = y_spromFus{1,4}.*ones(1,length(yps1_p7_x));

yps1_p4_x = [y_xcoords{5,1:2}(1):y_xcoords{5,1:2}(2)];
yps1_p4 = y_promFus{1,5}.*ones(1,length(yps1_p4_x));
yps1_p4_s = y_spromFus{1,5}.*ones(1,length(yps1_p4_x));

yps1_p5_x = [y_xcoords{6,1:2}(1):y_xcoords{6,1:2}(2)];
yps1_p5 = y_promFus{1,6}.*ones(1,length(yps1_p5_x));
yps1_p5_s = y_spromFus{1,6}.*ones(1,length(yps1_p5_x));

g_promFus = meanPlate1(2,[7:11]);
g_spromFus = sdPlate1(2,[7:11]);
g_xcoords = gyp7_indices([3:4,10,5:6],1:2);

gyp7_p2_x = [g_xcoords(1,1):g_xcoords(1,2)];
gyp7_p2 = g_promFus{1,1}.*ones(1,length(gyp7_p2_x));
gyp7_p2_s = g_spromFus{1,1}.*ones(1,length(gyp7_p2_x));

gyp7_p3_x = [g_xcoords(2,1):g_xcoords(2,2)];
gyp7_p3 = g_promFus{1,2}.*ones(1,length(gyp7_p3_x));
gyp7_p3_s = g_spromFus{1,2}.*ones(1,length(gyp7_p3_x));

gyp7_p9_x = [g_xcoords(3,1):g_xcoords(3,2)];
gyp7_p9 = g_promFus{1,3}.*ones(1,length(gyp7_p9_x));
gyp7_p9_s = g_spromFus{1,3}.*ones(1,length(gyp7_p9_x));

gyp7_p4_x = [g_xcoords(4,1):g_xcoords(4,2)];
gyp7_p4 = g_promFus{1,4}.*ones(1,length(gyp7_p4_x));
gyp7_p4_s = g_spromFus{1,4}.*ones(1,length(gyp7_p4_x));

gyp7_p5_x = [g_xcoords(5,1):g_xcoords(5,2)];
gyp7_p5 = g_promFus{1,5}.*ones(1,length(gyp7_p5_x));
gyp7_p5_s = g_spromFus{1,5}.*ones(1,length(gyp7_p5_x));

figure(1); 
shadedErrorBar(yps1_p1_x, yps1_p1, yps1_p1_s,'b',1); hold on;
shadedErrorBar(yps1_p2_x, yps1_p2, yps1_p2_s,'b',1);
shadedErrorBar(yps1_p3_x, yps1_p3, yps1_p3_s,'b',1);
shadedErrorBar(yps1_p7_x, yps1_p7, yps1_p7_s,'b',1);
shadedErrorBar(yps1_p4_x, yps1_p4, yps1_p4_s,'b',1);
shadedErrorBar(yps1_p5_x, yps1_p5, yps1_p5_s,'b',1);

shadedErrorBar(gyp7_p2_x, gyp7_p2, gyp7_p2_s,'k',1); hold on;
shadedErrorBar(gyp7_p3_x, gyp7_p3, gyp7_p3_s,'k',1);
shadedErrorBar(gyp7_p9_x, gyp7_p9, gyp7_p9_s,'k',1);
shadedErrorBar(gyp7_p4_x, gyp7_p4, gyp7_p4_s,'k',1);
shadedErrorBar(gyp7_p5_x, gyp7_p5, gyp7_p5_s,'k',1);
box off; xlabel('bp from TSS'); ylabel('norm H3 Occupancy (% IP/ % IP actin)'); %axis tight; 
title('Prom Fusion H3');
%%
% comparing endogenous pGYP7 and pYPS1
y_promFus = meanPlate1(3,[1:6]);
y_spromFus = sdPlate1(3,[1:6]);
y_xcoords = yps1_indices([3:5,8,6:7],1:2);

yps1_p1_x = [y_xcoords{1,1:2}(1):y_xcoords{1,1:2}(2)];
yps1_p1 = y_promFus{1,1}.*ones(1,length(yps1_p1_x));
yps1_p1_s = y_spromFus{1,1}.*ones(1,length(yps1_p1_x));

yps1_p2_x = [y_xcoords{2,1:2}(1):y_xcoords{2,1:2}(2)];
yps1_p2 = y_promFus{1,2}.*ones(1,length(yps1_p2_x));
yps1_p2_s = y_spromFus{1,2}.*ones(1,length(yps1_p2_x));

yps1_p3_x = [y_xcoords{3,1:2}(1):y_xcoords{3,1:2}(2)];
yps1_p3 = y_promFus{1,3}.*ones(1,length(yps1_p3_x));
yps1_p3_s = y_spromFus{1,3}.*ones(1,length(yps1_p3_x));

yps1_p7_x = [y_xcoords{4,1:2}(1):y_xcoords{4,1:2}(2)];
yps1_p7 = y_promFus{1,4}.*ones(1,length(yps1_p7_x));
yps1_p7_s = y_spromFus{1,4}.*ones(1,length(yps1_p7_x));

yps1_p4_x = [y_xcoords{5,1:2}(1):y_xcoords{5,1:2}(2)];
yps1_p4 = y_promFus{1,5}.*ones(1,length(yps1_p4_x));
yps1_p4_s = y_spromFus{1,5}.*ones(1,length(yps1_p4_x));

yps1_p5_x = [y_xcoords{6,1:2}(1):y_xcoords{6,1:2}(2)];
yps1_p5 = y_promFus{1,6}.*ones(1,length(yps1_p5_x));
yps1_p5_s = y_spromFus{1,6}.*ones(1,length(yps1_p5_x));

g_promFus = meanPlate1(3,[7:11]);
g_spromFus = sdPlate1(3,[7:11]);
g_xcoords = gyp7_indices([3:4,10,5:6],1:2);

gyp7_p2_x = [g_xcoords(1,1):g_xcoords(1,2)];
gyp7_p2 = g_promFus{1,1}.*ones(1,length(gyp7_p2_x));
gyp7_p2_s = g_spromFus{1,1}.*ones(1,length(gyp7_p2_x));

gyp7_p3_x = [g_xcoords(2,1):g_xcoords(2,2)];
gyp7_p3 = g_promFus{1,2}.*ones(1,length(gyp7_p3_x));
gyp7_p3_s = g_spromFus{1,2}.*ones(1,length(gyp7_p3_x));

gyp7_p9_x = [g_xcoords(3,1):g_xcoords(3,2)];
gyp7_p9 = g_promFus{1,3}.*ones(1,length(gyp7_p9_x));
gyp7_p9_s = g_spromFus{1,3}.*ones(1,length(gyp7_p9_x));

gyp7_p4_x = [g_xcoords(4,1):g_xcoords(4,2)];
gyp7_p4 = g_promFus{1,4}.*ones(1,length(gyp7_p4_x));
gyp7_p4_s = g_spromFus{1,4}.*ones(1,length(gyp7_p4_x));

gyp7_p5_x = [g_xcoords(5,1):g_xcoords(5,2)];
gyp7_p5 = g_promFus{1,5}.*ones(1,length(gyp7_p5_x));
gyp7_p5_s = g_spromFus{1,5}.*ones(1,length(gyp7_p5_x));

figure(2); 
shadedErrorBar(yps1_p1_x, yps1_p1, yps1_p1_s,'b',1); hold on;
shadedErrorBar(yps1_p2_x, yps1_p2, yps1_p2_s,'b',1);
shadedErrorBar(yps1_p3_x, yps1_p3, yps1_p3_s,'b',1);
shadedErrorBar(yps1_p7_x, yps1_p7, yps1_p7_s,'b',1);
shadedErrorBar(yps1_p4_x, yps1_p4, yps1_p4_s,'b',1);
shadedErrorBar(yps1_p5_x, yps1_p5, yps1_p5_s,'b',1);

shadedErrorBar(gyp7_p2_x, gyp7_p2, gyp7_p2_s,'k',1); hold on;
shadedErrorBar(gyp7_p3_x, gyp7_p3, gyp7_p3_s,'k',1);
shadedErrorBar(gyp7_p9_x, gyp7_p9, gyp7_p9_s,'k',1);
shadedErrorBar(gyp7_p4_x, gyp7_p4, gyp7_p4_s,'k',1);
shadedErrorBar(gyp7_p5_x, gyp7_p5, gyp7_p5_s,'k',1);
box off; xlabel('bp from TSS'); ylabel('norm H3 Occupancy (% IP/ % IP actin)'); %axis tight; 
title('Endogenous H3');
end