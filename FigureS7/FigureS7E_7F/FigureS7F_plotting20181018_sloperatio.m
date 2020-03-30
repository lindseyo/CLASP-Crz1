function [flag] = FigureS7F_plotting20181018_sloperatio()
flag = 1;
%%
nucOcc = [.98
.93
.42
.31

.56
.38
]; 

sensi = [0.56
1.09
1.33
1.15

1.22
0.83
];

genes = {
'Gyp7'
'Cmk2'
'Yps1'
'Put1'
'Mep1'
'Pun1'};

% plot 
figure(1); plot(nucOcc, sensi,'.','markersize',20); box off; axis tight;
R = corrcoef(nucOcc,sensi); R_corr = R(2)^2

title('H3 occ vs sensitivity to short pulses');
xlabel('H3 score');
ylabel('slope ratio');

dx=1; dy=1; 
text(nucOcc+dx,sensi+dy,genes);
text(.8,80,['R^2=',num2str(R_corr)])
