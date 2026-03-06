%% GT_Figure3
% This code reproduces CPP results in Figure 3, S3, S4, S5 and S6 

%% Pick channel within posterior cluster with max amplitude across conditions
allBest = []; 
RIDE = []; 
%Load dummy EEG to have scalp locations
EEG = pop_loadset([exp.filepath '/csd_preproc_P1_' exp.thisExpCode '_P2.set']);

%% Set plotting parameters
plt.opti = 0;            %Select best channel for each subject?
ds = 1;                  %Set to 1 if data are NOT downsampled, set to (200/512) if ds to 200;
plt.RIDE = 0;            %Use RIDE to remove SL potentials? (Not used in paper)
sr = EEG.srate;          %Original sampling rate
chans = exp.cppChans;    %Cluster of posterior channels to plot CPP;
plt.epoch = 'P1'; e = 1; %Data aligned to P1, epoch = 1 to find correct plotting pars. 
plt.baseline = 0;        %Baseline data? (set to 0 because they are baselined already )
plt.font = 11;
plt.ylims = [-5 25];
lineType = '-';

%% Load data with all sortings of interest 
load(['C:\Users\elisa\Desktop\GapsTask\' exp.thisExpName '\Data\plot data\' exp.thisExpCode '_P1_allP_ERP_averageData.mat'])

%% PLOT CPP 
close all 

sortBy.acc = 0;                     %Sort trials by accuracy? Set to 0 for main figures (plot all trials, pooling error-correct), set to 1 for Fig. S6
plt.acc = 1;                        %If sorting by accuracy, which trials to plot? 1 Error, 2 Correct

exp.cppChans = exp.cppChannels;     %[4,32,31,30,21,17,18,5,19,20];
par.minus_longGap = 1;              %Subtract longest gap condition from plots? Set to 1 for Fig. S4 panels B/D
plt.onePulse = par.minus_longGap;   %Old name, same function

% CPP-P1
if sortBy.acc == 1
    [data_allC_p1,dataTopo_p1] = GT_Fig3_CPP_P1_acc(exp, plt, EEG,RIDE,allBest,allP);
else
    [data_allC_p1,dataTopo_p1] = GT_Fig3_CPP_P1(exp, plt, EEG,RIDE,allBest,allP);
end

% CPP-P2
if sortBy.acc == 1
    [data_allC_p2,dataTopo_p2] = GT_Fig3_CPP_P2_acc(exp, plt, EEG,RIDE,allBest,allP);
else
    [data_allC_p2,dataTopo_p2] = GT_Fig3_CPP_P2(exp, plt, EEG,RIDE,allBest,allP);
end

chans = exp.cppChans;

%% Figure 3
f = figure; 
f.Units = 'centimeters'; f.OuterPosition = [0.5,0.5,13,21];
lineType = '-';
plt.font = 10;
par.testLatency = [0.3,0.7];

gapIdx = 2:length(exp.gapDur{exp.thisExp});
plt.box = 1;
ylims = [-2.5,20];

subplot(3,2,1)
fill([0.6*sr, 0.6*sr,0.8*sr,0.8*sr],[ylims, flip(ylims)], 'k', 'LineStyle','none','facealpha', 0.05); hold on; %Plot pulses in time
for p1 = 1:2
    thisCol = flip(plt.satPalettes{1,p1*p1});
    stdshade(squeeze(nanmean(data_allC_p1(:,p1,round(1.3*sr):round(2.3*sr),:)))', plt.alpha, thisCol, lineType,[],plt.smooth) %Plot ERP traces
    hold on;
end

set(gca, 'fontsize', plt.font)
xticklabels({'P1', '0.5', '1'})
xticks([[0.2,0.7,1.2]*512])
xline(0.2*512)
ylabel(plt.ERP.ylabel);
ylim(ylims)

subplot(3,2,3)
fill([0.6*sr, 0.6*sr,0.8*sr,0.8*sr],[ylims, flip(ylims)], 'k', 'LineStyle','none','facealpha', 0.05); hold on; 

ylims = [-2.5,10];
c = 5;
for p1 = 1:2
    for p2 = 1:2
        c = c-1;
        thisCol = hex2rgb((plt.Rcolors{c}));
        baselines =  squeeze(nanmean(nanmean(data_allC_p2(gapIdx,p1,p2,0.15*sr:0.2*sr,:),4),1));
        stdshade(squeeze(nanmean(data_allC_p2(gapIdx,p1,p2,:,:),1))'-baselines, plt.alpha, thisCol, lineType,[],plt.smooth) %Plot ERP traces
        hold on;
    end
end
ylim([ylims]); ylabel(plt.ERP.ylabel);
set(gca, 'fontsize', plt.font)
xticklabels({'P2', '0.5', '1'})
xticks([[0.2,0.7,1.2]*512])
xline(0.2*512)
f.Units = 'centimeters';
xlabel('Time (s)'); 

subplot(3,2,5)
fill([0.6*sr, 0.6*sr,0.8*sr,0.8*sr],[ylims, flip(ylims)], 'k', 'LineStyle','none','facealpha', 0.05); hold on; 

for p1 = 1:2
        thisCol = flip(plt.satPalettes{1,p1*p1});
        baselines =  squeeze(nanmean(nanmean(nanmean(data_allC_p2(gapIdx,p1,:,0.15*sr:0.2*sr,:),4),3),1));
        stdshade(squeeze(nanmean(nanmean(data_allC_p2(gapIdx,p1,:,:,:),1),3))'-baselines, plt.alpha, thisCol, lineType,[],plt.smooth) %Plot ERP traces
        hold on;
end

ylabel(plt.ERP.ylabel);
xlabel('Time (s)'); 
ylim(ylims)
set(gca, 'fontsize', plt.font)
xticklabels({'P2', '0.5', '1'})
xticks([[0.2,0.7,1.2]*512])
xline(0.2*512)

subplot(3,2,6)
fill([0.6*sr, 0.6*sr,0.8*sr,0.8*sr],[ylims, flip(ylims)], 'k', 'LineStyle','none','facealpha', 0.05); hold on; %Plot topo time
for p2 = 1:2
        thisCol = flip(plt.satPalettes{1,p2*p2});
        baselines =  squeeze(nanmean(nanmean(nanmean(data_allC_p2(gapIdx,:,p2,0.15*sr:0.2*sr,:),4),2),1));
        stdshade(squeeze(nanmean(nanmean(data_allC_p2(gapIdx,:,p2,:,:),1),2))'-baselines, plt.alpha, thisCol, ':',[],plt.smooth) %Plot ERP traces
        hold on;
end
ylim([ylims]); ylabel(plt.ERP.ylabel);

set(gca, 'fontsize', plt.font)
xticklabels({'P2', '0.5', '1'})
xticks([[0.2,0.7,1.2]*512])
xline(0.2*512)
f.Units = 'centimeters';
xlabel('Time (s)'); 


% Plot CPP-P1 & CPP-P2 regressions, with cluster stats against zero 
allModels = {'CPP_P1', 'CPP_P2'};
for m = 1:length(allModels)
par.modelName = allModels{m}; %For CPP-P2 evoked activities 
par.baseline = 1;
%Regressions are fit using the GT_fitRegressions function, which produces 
%this file. 
load([exp.plotpath, '/' exp.thisExpCode '_' par.modelName '_behRegression_b' num2str(par.baseline) '.mat' ]);

for i = 1:size(GA.p1coh.data,1)
    permData.p1{i} = squeeze(GA.p1coh.data(i,:)); 
    permData.p2{i} = squeeze(GA.p2coh.data(i,:)); 
end

cluster = GT_cluster_CPPcoh(exp, par,permData); 

subplot(3,2,2+((m-1)*2)); 
stdshade(squeeze(GA.p1coh.data), 0.05,'k', '-', [],1) % plt.satPalettes{1,3}

%Add shaded area to restrict testing time 
ylims = [-0.1,0.1];
% if m == 1; 
    hold on; 
    fill([par.testLatency(1)+0.2,par.testLatency(1)+0.2, par.testLatency(2)+0.2,par.testLatency(2)+0.2]*512,[ylims, flip(ylims)], 'k', 'LineStyle','none','facealpha', 0.05); %Plot pulses in time
% end
hold on; stdshade(squeeze(GA.p2coh.data), 0.05, [0.4 0.4 0.4], ':', [],1)
hold on
yline(0); xline(0.2*512);
xticks([0.2,0.7,1.2]*512); xticklabels({par.modelName(end-1:end), '0.5', '1'});
if m == 2; xlabel('Time (s)'); end

lt = {'-', ':'}; 
statLabs = {'stat1', 'stat2'};
colors = {'k', [0.4 0.4 0.4]};
for stats = 1:length(fields(cluster));
    thisCol = colors{stats};% 'k'; %flip(plt.satPalettes{1,3});
    lineType = lt{stats};
    thisStat = statLabs{stats};
    idxPos = []; idxNeg = [];
    if isfield(cluster.(thisStat), 'posclusters')
        for i = 1:size(cluster.(thisStat).posclusters,2)
            if i >0
                if cluster.(thisStat).posclusters(i).prob < 0.05
                    idxPos = [idxPos,find(cluster.(thisStat).posclusterslabelmat == i)];
                    
                end
            end
        end
    end
    if isfield(cluster.(thisStat), 'negclusters')
        for i = 1:size(cluster.(thisStat).negclusters,2)
            if i >0
                if cluster.(thisStat).negclusters(i).prob < 0.05
                    idxNeg = [idxNeg,find(cluster.(thisStat).negclusterslabelmat == i)];
                end
            end
        end
    end
    hold on;     
    x = unique([idxPos, idxNeg])+(par.testLatency(1)+0.2)*512; 
    ypos = -0.085-0.0075*(stats-1);

    plot(x,repmat(ypos,1,length(x)), '*', 'color', thisCol, 'markersize',1.5); hold on;
    clear x; clear idxPos; clear idxNeg; 

end
set(gca, 'fontsize', plt.font)
xticklabels({['P' num2str(m)], '0.5', '1'})
xticks([[0.2,0.7,1.2]*512])
xline(0.2*512)
ylabel('Regression coef. (a.u.)');
ylim(ylims); yticks([-0.1,0,0.1]);

end

cLim = [-3 3];
axes('position',[.625 .835 .10 .10]);
g = 1:length(exp.gapDur{exp.thisExp}); %Include all gaps for P1 topos
diffTopo_p1 = squeeze(nanmean(dataTopo_p1(g,1,:,:,:),1) - nanmean(dataTopo_p1(:,2,:,:,:),1)); 


GT_topoplot(nanmean(diffTopo_p1,2),EEG.chanlocs, ...
    'numcontour', 0, ...
    'plotrad',0.65, ...
    'headrad', 0.5,...
    'emarker2', {chans,'*','k',1,1}); % To plot location of the selected channel
axes('position',[.625 .405 .10 .10]);

g = 2:length(exp.gapDur{exp.thisExp}); %Include only non-zero gaps for P2 topos
diffTopo_p1 = squeeze(nanmean(nanmean(dataTopo_p2(g,1,:,:,:,:),1),3) - nanmean(nanmean(dataTopo_p2(g,2,:,:,:,:),1),3)); 
diffTopo_p2 = squeeze(nanmean(nanmean(dataTopo_p2(g,:,1,:,:,:),1),2) - nanmean(nanmean(dataTopo_p2(g,:,2,:,:,:),1),2)); 

GT_topoplot(nanmean(diffTopo_p1,2),EEG.chanlocs, ...
    'numcontour', 0, ...
    'plotrad',0.65, ...
    'headrad', 0.5,...
    'emarker2', {chans,'*','k',1,1}); % To plot location of the selected channel

axes('position',[0.625 .535 .10 .10]);
GT_topoplot(nanmean(diffTopo_p2,2),EEG.chanlocs, ...
    'numcontour', 0, ...
    'plotrad',0.65, ...
    'headrad', 0.5,...
    'emarker2', {chans,'*','k',1,1}); % To plot location of the selected channel

exportgraphics(f, [exp.figPath 'GT_Fig3_Exp' num2str(exp.thisExp) '.tiff'], 'resolution', 600); 
exportgraphics(f, [exp.figPath 'GT_Fig3_Exp' num2str(exp.thisExp) '.pdf'], 'resolution', 600); 

%% For Figure S5, re-run code above setting par.minus_longGap = 0 and then run this 
f = figure; 
f.Units = 'centimeters';
f.OuterPosition = [0.5,0.5,28,9];

gapIdx = 2:length(exp.gapDur{exp.thisExp});
plt.box = 1;
ylims = [-2.5,25];
lineType = '-';
subplot(1,4,1)
fill([0.6*sr, 0.6*sr,0.8*sr,0.8*sr],[ylims, flip(ylims)], 'k', 'LineStyle','none','facealpha', 0.05); hold on; %Plot pulses in time
for p1 = 1:2
    thisCol = flip(plt.satPalettes{1,p1*p1});
    stdshade(squeeze(nanmean(data_allC_p1(:,p1,round(1.3*sr):round(2.3*sr),:)))', plt.alpha, thisCol, lineType,[],plt.smooth) %Plot ERP traces
    hold on;
end

set(gca, 'fontsize', plt.font)
xticklabels({'P1', '0.5', '1'})
xticks([[0.2,0.7,1.2]*512])
xline(0.2*512)
ylabel(plt.ERP.ylabel);
xlabel('Time (s)'); 
ylim(ylims)

subplot(1,4,2)
fill([0.6*sr, 0.6*sr,0.8*sr,0.8*sr],[ylims, flip(ylims)], 'k', 'LineStyle','none','facealpha', 0.05); hold on; 

c = 5;
for p1 = 1:2
    for p2 = 1:2
        c = c-1;
        thisCol = hex2rgb((plt.Rcolors{c}));
        baselines =  squeeze(nanmean(nanmean(data_allC_p2(gapIdx,p1,p2,0.15*sr:0.2*sr,:),4),1));
        stdshade(squeeze(nanmean(data_allC_p2(gapIdx,p1,p2,:,:),1))'-baselines, plt.alpha, thisCol, lineType,[],plt.smooth) %Plot ERP traces
        hold on;
    end
end
ylim([ylims]); 
set(gca, 'fontsize', plt.font)
xticklabels({'P2', '0.5', '1'})
xticks([[0.2,0.7,1.2]*512])
xline(0.2*512)
f.Units = 'centimeters';
xlabel('Time (s)'); 

subplot(1,4,3)
fill([0.6*sr, 0.6*sr,0.8*sr,0.8*sr],[ylims, flip(ylims)], 'k', 'LineStyle','none','facealpha', 0.05); hold on; 

for p1 = 1:2
        thisCol = flip(plt.satPalettes{1,p1*p1});
        baselines =  squeeze(nanmean(nanmean(nanmean(data_allC_p2(gapIdx,p1,:,0.15*sr:0.2*sr,:),4),3),1));
        stdshade(squeeze(nanmean(nanmean(data_allC_p2(gapIdx,p1,:,:,:),1),3))'-baselines, plt.alpha, thisCol, lineType,[],plt.smooth) %Plot ERP traces
        hold on;
end

ylim([ylims]); yticklabels('');
set(gca, 'fontsize', plt.font)
xticklabels({'P2', '0.5', '1'})
xticks([[0.2,0.7,1.2]*512])
xline(0.2*512)
f.Units = 'centimeters';
xlabel('Time (s)'); 

subplot(1,4,4)
fill([0.6*sr, 0.6*sr,0.8*sr,0.8*sr],[ylims, flip(ylims)], 'k', 'LineStyle','none','facealpha', 0.05); hold on; %Plot topo time
for p2 = 1:2
        thisCol = flip(plt.satPalettes{1,p2*p2});
        baselines =  squeeze(nanmean(nanmean(nanmean(data_allC_p2(gapIdx,:,p2,0.15*sr:0.2*sr,:),4),2),1));
        stdshade(squeeze(nanmean(nanmean(data_allC_p2(gapIdx,:,p2,:,:),1),2))'-baselines, plt.alpha, thisCol, ':',[],plt.smooth) %Plot ERP traces
        hold on;
end
ylim([ylims]); yticklabels('');
set(gca, 'fontsize', plt.font)
xticklabels({'P2', '0.5', '1'})
xticks([[0.2,0.7,1.2]*512])
xline(0.2*512)
f.Units = 'centimeters';
xlabel('Time (s)'); 

cLim = [-3 3];
axes('position',[.054 .65 .23 .23]);
g = 1:length(exp.gapDur{exp.thisExp});
% diffTopo_p1 = squeeze(nanmean(dataTopo_p1(g,1,:,:),1)) - nanmean(dataTopo_p1(g,2,:,:),1); 
diffTopo_p1 = squeeze(nanmean(dataTopo_p1(:,1,:,:,:),1) - nanmean(dataTopo_p1(:,2,:,:,:),1)); 

GT_topoplot(nanmean(diffTopo_p1,2),EEG.chanlocs, ...
    'numcontour', 0, ...
    'plotrad',0.65, ...
    'headrad', 0.5,...
    'emarker2', {chans,'*','k',1,1}); % To plot location of the selected channel

axes('position',[.464 .65 .23 .23]);
g = 1:length(exp.gapDur{exp.thisExp});
diffTopo_p1 = squeeze(nanmean(nanmean(dataTopo_p2(g,1,:,:,:,:),1),3) - nanmean(nanmean(dataTopo_p2(g,2,:,:,:,:),1),3)); 
diffTopo_p2 = squeeze(nanmean(nanmean(dataTopo_p2(g,:,1,:,:,:),1),2) - nanmean(nanmean(dataTopo_p2(g,:,2,:,:,:),1),2)); 

GT_topoplot(nanmean(diffTopo_p1,2),EEG.chanlocs, ...
    'numcontour', 0, ...
    'plotrad',0.65, ...
    'headrad', 0.5,...
    'emarker2', {chans,'*','k',1,1}); % To plot location of the selected channel
% caxis(cLim)


axes('position',[.676 .65 .23 .23]);
GT_topoplot(nanmean(diffTopo_p2,2),EEG.chanlocs, ...
    'numcontour', 0, ...
    'plotrad',0.65, ...
    'headrad', 0.5,...
    'emarker2', {chans,'*','k',1,1}); % To plot location of the selected channel

exportgraphics(f, [exp.figPath 'GT_FiS5_Exp' num2str(exp.thisExp) '.tiff'], 'resolution', 600); 
exportgraphics(f, [exp.figPath 'GT_FigS5_Exp' num2str(exp.thisExp) '.pdf'], 'resolution', 600); 