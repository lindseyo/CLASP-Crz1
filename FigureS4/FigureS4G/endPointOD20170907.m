function [finalMatrix] = endPointOD20170907(di,filename)
A=xlsread(strcat(di,'\',filename));
%%
% size(A)
B=A;
B(~any(~isnan(B),2),:)=[];
% size(B)
C=B;
C(:,~any(~isnan(C), 1))=[];
% size(C)

finalMatrix = C(5:12,:);