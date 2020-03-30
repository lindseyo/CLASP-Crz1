function [flag] = FigureS1A()
flag = 1;
%%
load('workspaceForCyclopTR20161026.mat')

[val, ai] = sort(mWT123_nucleusTRs);

figure(1); plot(val, '.','markersize',30)
box off; axis tight;
title('TF Nuclear Localization State');
xlabel('Transcription Factors'); ylabel('Nuclear Metric')

end