function [flag] = rna_seq_data_analysis_20170927_pilotCrzALA36samples()
flag = 1;
%% Script: This script run a series of functions that analyzes and plots RNAseq data
% 
% Data Pre-processing: fastq files are run through Kallisto software and
% the input to functions is an "abundance.tsv" file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loading datasets
% load in matrix
load('countmatrix_files.mat')
% plot out the mean and std
%figure(1); scatter(meanLogCounts, stdLogCounts);
% visually use a cut off of 1
log2NormMatrixFiltered = log2NormMatrix;
indx = find(meanLogCounts > 2);
log2NormMatrixFiltered = log2NormMatrix(indx,:);
log2NormMatrix = log2NormMatrixFiltered;

%%
% make a mapping of keys to values to keep the samples straight
keys = cell(1,length(samples)); values = [];
for i = 1:length(samples)
    keys{i} = samples{i}{1};
    values(i) = i;
end
mapSamples = containers.Map(keys,values);

% get the common/scientific gene names
gene_names = cell(1,length(genes_comscinames));
for j = 1:length(genes_comscinames);
    gene_names{j} = genes_comscinames{j}{1};
end
gene_names = gene_names(indx);
% get the scientific gene names
gene_names_sci = cell(1,length(genes_scinames));
for k = 1:length(genes_scinames);
    gene_names_sci{k} = genes_scinames{k}{1};
end
gene_names_sci = gene_names_sci(indx);

%% hierarchical clustering dendrogram
% samples 1, 13, 17,20
temp=log2NormMatrix(:,[1,13,17,20]);
cgo = clustergram(temp,'Colormap',redbluecmap,'Standardize','Row');
end
