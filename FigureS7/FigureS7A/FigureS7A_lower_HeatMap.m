function [flag] = FigureS7A_upper_HeatMap()
flag = 1;

%% Parameter Relationships
%load set roffTHR at 0.5
%load('20181217_pt02ronpt6kon_3state_1RANGE_roffTHR','paramset');
load('20181217_2ronpt6kon_3state_1RANGE_roffTHR','paramset');
Y3 = paramset(:,1:1000);

param10000 = [Y3];

%% DON'T CHANGE STUFF -- this is the script to generate PARTS OF FIGURE 5
% load trajectories - set roffTHR at 0.5
%load('20181217_3state_ronTHR_pt02ronpt6kon_output1000.mat','allsGen');
load('20181217_3state_ronTHR_2ronpt6kon_output1000.mat','allsGen');

X3 = allsGen(1:1000,:);

allsGen = [X3];

% model end point
val = 60; % 5hrs  % could try 180
mod2 = {allsGen(:,1),allsGen(:,2),allsGen(:,3),allsGen(:,4),allsGen(:,5)};
mod40 = {allsGen(:,1),allsGen(:,6),allsGen(:,7),allsGen(:,8),allsGen(:,9)};
% model transferfunction
mod_transffn= {allsGen(:,1),allsGen(:,11),allsGen(:,10),allsGen(:,12)};

%%
load('20181127_gyp_pulsedNconst_data','constGE','sconstGE','pulsedGE','spulsedGE') % newer data set
GE2m = pulsedGE; GEdur = constGE; sGE2m = spulsedGE; sGEdur = sconstGE;
% define Xaxis
load('Input_idealizedLight.mat', 'lightInputs')
kfact = 1000;
orig_lightInputs = lightInputs;
m2_xaxis=[orig_lightInputs(1).auc,orig_lightInputs(2).auc./kfact,orig_lightInputs(3).auc./kfact,orig_lightInputs(4).auc./kfact,orig_lightInputs(5).auc./kfact];
single_xaxis=[orig_lightInputs(1).auc,orig_lightInputs(6).auc./kfact,orig_lightInputs(7).auc./kfact,orig_lightInputs(8).auc./kfact,orig_lightInputs(9).auc./kfact];
% define for transfer fn X axis
actiTF = [0,orig_lightInputs(11).max./kfact,orig_lightInputs(10).max./kfact,orig_lightInputs(12).max./kfact];

%%
% scale all trajs to experimental values
% best fit line through the data
% 2min
p = polyfit(m2_xaxis,GE2m,1);
x2 = polyval(p,m2_xaxis);
figure(10); errorbar(m2_xaxis,GE2m,sGE2m,'.r','markersize',20,'linewidth',2);...
    hold on; plot(m2_xaxis,x2,'r'); xlabel('cum nuc occpuancy'); ylabel('FITC/SSC');
% 40min
p = polyfit(single_xaxis,GEdur,1);
xdur = polyval(p,single_xaxis);
errorbar(single_xaxis,GEdur,sGEdur,'.b','markersize',20,'linewidth',2);...
    plot(single_xaxis,xdur,'b'); axis tight; title('best fit line to data');

% scale to max and min of best fit line
a = min(xdur); b = max(xdur);
scaled_2_endpt = [];
scaled_C_endpt = [];
mod_trans = [];
for i = 1:length(mod2{1}) 
    % find error for endpoint
    x = [mod2{1}{i}(end-val),mod2{2}{i}(end-val),mod2{3}{i}(end-val),mod2{4}{i}(end-val),mod2{5}{i}(end-val)];
    y = [mod40{1}{i}(end-val),mod40{2}{i}(end-val),mod40{3}{i}(end-val),mod40{4}{i}(end-val),mod40{5}{i}(end-val)];
    ix = (b-a).*(x-min(y))./(max(y)-min(y))+a;
    iy = (b-a).*(y-min(y))./(max(y)-min(y))+a;
    scaled_2_endpt(i,:) = ix;
    scaled_C_endpt(i,:) = iy; 

end


%% indx_Rats -- feed original allsGen through this
% compute slope ratio
rat_slopes = [];
for i = 1:length(mod2{1})
    x = scaled_2_endpt(i,:);
    y = scaled_C_endpt(i,:);
    rat_slopes(i) = ((x(5)-x(1))./(m2_xaxis(5)-m2_xaxis(1)))./((y(5)-y(1))./(single_xaxis(5)-single_xaxis(1)));
end
[sortedRats, indx_Rats] = sort(rat_slopes);
% find slope ratio < 1
indx_slopeLT1 = find(sortedRats<1);
indx_slopeGT1 = find(sortedRats>1);

figure(1); plot(sortedRats); xlabel('parameters'); ylabel('slope ratios'); title('Ratios of Slopes');
axis tight; box off;
%% get allsGen in the order of rat_slopes
% re-organize all of the parameters and trajectories by indx_Rats
allsGen_bR = allsGen(indx_Rats,:);
scaled_2_endpt_bR = scaled_2_endpt(indx_Rats,:);
scaled_C_endpt_bR = scaled_C_endpt(indx_Rats,:);

param10000_bR = param10000(:,indx_Rats);

% set up kon koff values
kon = param10000_bR(1,:); koff = param10000_bR(2,:); ron = param10000_bR(8,:); roff = param10000_bR(9,:);


%% plotting scatter of parameters (by slope ratio) ID regions
% grid
% search for max of sortedRats
%gindx = find(log10(kon./koff) >= -2 & log10(kon./koff) <= 0.1 & log10(ron./roff) >= -2 & log10(ron./roff) <= 0.5);
    maxColor = max(sortedRats);
    minColor = min(sortedRats);
    colorV = linspace(minColor,maxColor,256); %256
    cMap =[zeros(sum(colorV <= 1),3);parula(256-sum(colorV <= 1))];

    LL = min(log10(kon./koff)):.05:max(log10(kon./koff)); %min(konkoff)
    MM = min(log10(ron./roff)):.05:max(log10(ron./roff));%max(ronroff);
    [x,y] = meshgrid(LL,MM);
    AA = log10(kon./koff);
    BB = log10(ron./roff);
    F = griddata(AA,BB,sortedRats,x,y);
    figure(8); sanePColor(LL,MM,F);  hold on; box off; axis tight; colorbar; colormap(cMap); %DEFINITELY CONSIDER THIS!
    title('SR w Linear Region + SR < 1 + QualFits'); xlabel('log10(kon/koff)'); ylabel('log10(ron/roff)');
    axis([-1 1 -2 0])
%axis([-1 1 -4 -2]);
%plot(log10(param10000_bR(1,indxFits)./param10000_bR(2,indxFits)),...
%    log10(param10000_bR(8,indxFits)./param10000_bR(9,indxFits)),'m.');