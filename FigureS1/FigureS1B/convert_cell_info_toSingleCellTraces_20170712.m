function [singleCell_traces] = convert_cell_info_toSingleCellTraces_20170712(cell_info, chan)
% The goal of this script is to deconstruct the "cell_info" structure into
% a structure with single cell traces
%
% INPUT:
%     cell_info -- structure with 4 fields: image/content, channel, position, time
% OUTPUT:
%     singleCell_traces -- structure/array with single cell traces for each
%     sample

% NOT DONE CAN'T GO FORWARD BC DON'T HAVE THE RIGHT PROCESSED DATA INFO

if chan == 'YFP'
    j = 2;
end

clear cell_info_cell; clear cell_info_cell1;
cell_info_cell = struct2cell(cell_info);
cell_info_cell1(:,:) = cell_info_cell(:,j,:);
Pos = cell2mat(cell_info_cell1(3,:));
uniqPos = unique(Pos);

for i = 1:numel(uniqPos) % loop through all samples
end