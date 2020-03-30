function [flag] = FigureS4A()
flag = 1;
load('imgAnalysis20180201_SYC511_pt2MCaCl2_3hrs_filtered.mat')
%load('imgAnalysis20180201_SYC511_pt2MCaCl2_3hrs.mat')

count = 0;
C = NaN(109,1000);
for i = 1:length(all_tracks_filt_vec)
    A = all_tracks_filt_vec{i}.pt2MCaCl2;
    for j = 1:length(A)
        B = A(j).nf.RFP;
        count = count + 1;
        C(count,1:length(B)) = B;
    end
end

% figure(1); 
% jj = 0;
% for ii = 78:89
%     jj = jj + 1;
%     subplot(3,4,jj); plot(C(ii,:)); ylim([1 7]); xlim([0 300])
% end

figure(2); plot(C(12,:));
figure(3); plot(C(23,:));
figure(4); plot(C(74,:));
%45
%54
%58
%74