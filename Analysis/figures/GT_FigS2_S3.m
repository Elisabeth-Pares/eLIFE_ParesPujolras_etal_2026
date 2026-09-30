%% Figure S2 & S3 - plot MBL results sorted by condition
GT_FigS2_S3

%% Set plotting parameters
ds = 1              %Set to 1 if data not downsampled, set to (200/512) if ds to 200;
%Load dummy EEG to have scalp locations
EEG = pop_loadset([exp.filepath '/csd_preproc_P1_' exp.thisExpCode '_P2.set']);
% load(['C:\Users\elisa\Desktop\GapsTask\' exp.thisExpName '\Data\plot data\' exp.thisExpCode '_P1_allP_TFA_average_beta.mat'])
load(['C:\Users\elisa\Desktop\GapsTask\' exp.thisExpName '\Data\plot data\' exp.thisExpCode '_P1_allP_TFA_averageDatabeta_short.mat'])

sr = EEG.srate; %Original sampling rate
plt.epoch = 'P1'; e = 1; %Data aligned to P1, epoch = 1 to find correct plotting pars.
plt.ylims = [-5 25];
lineType = '-';
plt.chansel = [3];  %As in main manuscript; plot only C3/C4

%% MBL-P1, with boxplots
[MBL_allC_p1,dataTopo_p1] = GT_Fig3_MBL_P1(exp, plt, EEG,allP);

%% MBL-P2, with boxplots
[MBL_allC_p2,dataTopo_p2]= GT_Fig3_MBL_P2(exp, plt, EEG,allP)

%% Plot averaged across gaps
allCols = [flip(plt.satPalettes{1,1*1}); flip(plt.satPalettes{1,2*2})];

f = figure;
f.Units = 'centimeters';
f.OuterPosition = [0.5,0.5,28,9];
sr = 200;
clear x; clear y; clear excursion
gapIdx = 2:length(exp.gapDur{exp.thisExp});
plt.box = 1;
ylims = [-1,1];
xlims = [1,1.2*200];

subplot(1,4,1)
fill([0.8*sr, 0.8*sr,1*sr,1*sr],[ylims, flip(ylims)], 'k', 'LineStyle','none','facealpha', 0.05); hold on; %Plot pulses in time
for p1 = 1:2
    thisCol = flip(plt.satPalettes{1,p1*p1});
    %Compute excursion
    preP2 = squeeze(nanmean(nanmean(MBL_allC_p1(gapIdx,1,p1,1:0.2*sr,:),4),1));
    endOfTrial = squeeze(nanmean(nanmean(MBL_allC_p1(gapIdx,1,p1,0.8*sr:1*sr,:),4),1));
    excursion(p1,:) = (endOfTrial-preP2);
    
    baselines = squeeze(nanmean(nanmean(MBL_allC_p1(gapIdx,1,p1,1:0.2*sr,:),4),1));
    stdshade(squeeze(nanmean(MBL_allC_p1(:,1,p1,0.2*sr:1.4*sr,:)))'-baselines, plt.alpha, thisCol, lineType,[],plt.smooth) %Plot ERP traces
    hold on;
end

set(gca, 'fontsize', plt.font)
xticklabels({'P1', '0.5', '1'})
xticks([[0.2,0.7,1.2]*200])
xline(0.2*200)
ylabel('MBL (dB)');
xlabel('Time (s)'); xlim(xlims);
ylim(ylims)

%Add bar graph with excursion
axes('position',[.215 .65 .05 .2]);
p1cols = [flip(plt.satPalettes{1,1});flip(plt.satPalettes{1,4})];
x = 1:2; y = (nanmean(excursion')); error = nanstd(excursion')/sqrt(size(excursion,2));
b = bar([1:2], y)
b.FaceColor = 'flat';
b.CData(:,:) = p1cols ;(allCols);
xticks([1:2]); xticklabels('');
ylabel(['\Delta MBL']); box off;
set(gca, 'fontsize', 8); hold on;
errorbar(x, y, error, '.k'); box off;


clear excursion

subplot(1,4,2)
fill([0.8*sr, 0.8*sr,1*sr,1*sr],[ylims, flip(ylims)], 'k', 'LineStyle','none','facealpha', 0.05); hold on;
allCols = [];
c = 5;
for p1 = 1:2
    for p2 = 1:2
        c = c-1;
        thisCol = hex2rgb((plt.Rcolors{c}));
        allCols = [allCols; thisCol];
        
        baselines = 0;% squeeze(nanmean(nanmean(MBL_allC_p2(gapIdx,p1,p2,1:0.2*sr,:),4),1));
        %
        %Compute excursion
        preP2 = squeeze(nanmean(nanmean(MBL_allC_p2(gapIdx,p1,p2,1:0.2*sr,:),4),1));
        endOfTrial = squeeze(nanmean(nanmean(MBL_allC_p2(gapIdx,p1,p2,0.8*sr:1*sr,:),4),1));
        excursion(c,:) = (endOfTrial-preP2);
        
        stdshade(squeeze(nanmean(MBL_allC_p2(gapIdx,p1,p2,:,:),1))'-baselines, plt.alpha, thisCol, '-',[],plt.smooth) %Plot ERP traces
        
    end
end
ylim([ylims]);
set(gca, 'fontsize', plt.font)
xticklabels({'P2', '0.5', '1'})
xticks([[0.2,0.7,1.2]*200])
xline(0.2*200)
f.Units = 'centimeters';
xlabel('Time (s)'); xlim(xlims);

%Add bar graph with excursion
axes('position',[.43 .65 .05 .2]);

if c == 1
    x = 1:4; y = (nanmean(excursion')); error = nanstd(excursion')/sqrt(size(excursion,2));
    b = bar(x, y)
    b.FaceColor = 'flat';
    b.CData(:,:) = flip(allCols);
    xticks(x); xticklabels('');
    ylabel(['\Delta MBL']); box off;
    set(gca, 'fontsize', 8); hold on;
    errorbar(x, y, error, '.k'); box off;
    
end
clear excursion


subplot(1,4,3)
fill([0.8*sr, 0.8*sr,1*sr,1*sr],[ylims, flip(ylims)], 'k', 'LineStyle','none','facealpha', 0.05); hold on;

for p1 = 1:2
    thisCol = flip(plt.satPalettes{1,p1*p1});
    baselines =  0%squeeze(nanmean(nanmean(nanmean(MBL_allC_p2(gapIdx,p1,:,1:0.2*sr,:),4),3),1));
    
    %Compute excursion
    preP2 = squeeze(nanmean(nanmean(nanmean(MBL_allC_p2(gapIdx,p1,:,1:0.2*sr,:),4),3),1));
    endOfTrial = squeeze(nanmean(nanmean(nanmean(MBL_allC_p2(gapIdx,p1,:,0.8*sr:1*sr,:),4),3),1));
    excursion(p1,:) = (endOfTrial-preP2);
    
    stdshade(squeeze(nanmean(nanmean(MBL_allC_p2(gapIdx,p1,:,:,:),1),3))'-baselines, plt.alpha, thisCol, lineType,[],plt.smooth) %Plot ERP traces
    hold on;
    
end

ylim([ylims]); yticklabels('');
set(gca, 'fontsize', plt.font)
xticklabels({'P2', '0.5', '1'})
xticks([[0.2,0.7,1.2]*200])
xline(0.2*200)
f.Units = 'centimeters';
xlabel('Time (s)'); xlim(xlims);

%Add bar graph with excursion
axes('position',[.6325 .65 .05 .2]);

p1cols = [flip(plt.satPalettes{1,1});flip(plt.satPalettes{1,4})];
if c == 1
    x = 1:2; y = (nanmean(excursion')); error = nanstd(excursion')/sqrt(size(excursion,2));
    b = bar([1:2], y)
    b.FaceColor = 'flat';
    b.CData(:,:) = p1cols ;(allCols);
    xticks([1:2]); xticklabels('');
    ylabel(['\Delta MBL']); box off;
    set(gca, 'fontsize', 8); hold on;
    errorbar(x, y, error, '.k'); box off;
    
end
clear excursion



subplot(1,4,4)
fill([0.8*sr, 0.8*sr,1*sr,1*sr],[ylims, flip(ylims)], 'k', 'LineStyle','none','facealpha', 0.05); hold on; %Plot topo time
for p2 = 1:2
    thisCol = flip(plt.satPalettes{1,p2*p2});
    
    %Compute excursion
    preP2 = squeeze(nanmean(nanmean(nanmean(MBL_allC_p2(gapIdx,:,p2,1:0.2*sr,:),4),2),1));
    endOfTrial = squeeze(nanmean(nanmean(nanmean(MBL_allC_p2(gapIdx,:,p2,0.8*sr:1*sr,:),4),2),1));
    excursion(p2,:) = (endOfTrial-preP2);
    
    baselines = 0% squeeze(nanmean(nanmean(nanmean(MBL_allC_p2(gapIdx,:,p2,1:0.2*sr,:),4),2),1));
    stdshade(squeeze(nanmean(nanmean(MBL_allC_p2(gapIdx,:,p2,:,:),1),2))'-baselines, plt.alpha, thisCol, ':',[],plt.smooth) %Plot ERP traces
    
    hold on;
end

ylim([ylims]); yticklabels('');
set(gca, 'fontsize', plt.font)
xticklabels({'P2', '0.5', '1'})
xticks([[0.2,0.7,1.2]*200])
xline(0.2*200)
f.Units = 'centimeters';
xlabel('Time (s)'); xlim(xlims);

%Add bar graph with excursion
axes('position',[.84 .65 .05 .2]);
p1cols = [flip(plt.satPalettes{1,1});flip(plt.satPalettes{1,4})];
if c == 1
    x = 1:2; y = (nanmean(excursion')); error = nanstd(excursion')/sqrt(size(excursion,2));
    b = bar([1:2], y)
    b.FaceColor = 'flat';
    b.CData(:,:) = p1cols ;(allCols);
    xticks([1:2]); xticklabels('');
    ylabel(['\Delta MBL']); box off;
    set(gca, 'fontsize', 8); hold on;
    errorbar(x, y, error, '.k'); box off;
end

exportgraphics(f, [exp.figPath 'GT_FigS2_Exp' num2str(exp.thisExp) '.tiff'], 'resolution', 600);
exportgraphics(f, [exp.figPath 'GT_FigS2_Exp' num2str(exp.thisExp) '.pdf'], 'resolution', 600);