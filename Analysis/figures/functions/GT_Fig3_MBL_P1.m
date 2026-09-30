function [alldiffData, dataTopo, baselinesSL] = GT_Fig3_MBL_P1(exp, plt, EEG,allP)

%% Figure 1: plot beta lateralisation stimulus locked
factor = 2; e= 1;

ds = (200/512);
clear diffData;
xlims = plt.xlims(e,1)*200/512:plt.xlims(e,2)*ds;
plt.epoch = 'P1'; e = 1;
plt.baseline = 1;
par.db = 1;
ylims = [-1.5,1.5]

f = figure

subplot(5,4,[9:11,13:15,17:19]);

xline(200, 'k');
xline([round(0.8*200),round(0.9*200)], '--k');

allData = [];
allcontraData = []; allipsiData = [];
allContra = []; allIpsi = [];
allDiff = []; dataTopo = [];

plt.rl = 0;
if plt.rl == 1; accs = [0,1];
else, accs = 0; end

for accs = accs
    plt.acc = accs;
    f = figure
    
    for g = 1:length(exp.gapDur{exp.thisExp})
        for p1 = 1:2
            for d = 1:length(exp.dirLabs)
                
                thisLinetype = '-';
                
                thisGap = exp.gapLabs{exp.thisExp}{g};
                thisP1 = exp.cohLabs{p1};
                thisDir = exp.dirLabs{d};
                
                if plt.acc == 1
                    thisGrouping= 'np_gap_p1coh_dir_acc_';
                    dataTag = ['twoPulse' thisGap thisP1 thisDir 'Correct'];
                else
                    thisGrouping= 'np_gap_p1coh_dir_';
                    dataTag = ['twoPulse' thisGap thisP1 thisDir];
                end
                
                thisCol = (plt.satPalettes{1,p1*p1});
                
                %Compute lateralisation based on cluster of left/right hemisphere electrodes
                if strcmp(thisDir, 'dotLeft'); ipsiChan = 1; contraChan =2; otherDir= 'dotRight';
                else, ipsiChan = 2; contraChan = 1; otherDir = 'dotLeft'; end
                
                if par.db == 1; unitTag = '(dB)'; if plt.rl == 1, ylims = [-2 1]; else, ylims = [-2 0.5]; end %[-1.5 1];
                    baselines = nanmean(allP.(plt.epoch).(thisGrouping).(dataTag).sl(:,[0.2*200:0.4*200],:).^factor,2);
                    thisData = 10*log10(allP.(plt.epoch).(thisGrouping).(dataTag).sl.^factor./baselines);
                    ipsiData = nanmean(thisData(ipsiChan,xlims,:),1)%.^factor; %IPSI channels
                    contraData = nanmean(thisData(contraChan,xlims,:),1)%.^factor; %IPSI channels
                end
                if plt.baseline == 1
                    baseTag = 'base';
                    allipsiData(d,:,:) = squeeze(ipsiData - nanmean(nanmean(thisData(ipsiChan,260:300,:),1),2));
                    allcontraData(d,:,:) = squeeze(contraData - nanmean(nanmean(thisData(contraChan,260:300,:),1),2));
                end
                
            end
            
            plotcontraData = squeeze(nanmean(allcontraData));
            plotipsiData= squeeze(nanmean(allipsiData));
           
            lat = [0.4 + exp.gapDur{exp.thisExp}(g) + 0.75]; %dur of two pulses, plus dur of gap, plus minimum lead out duration
            idx = find(EEG.times(round(xlims/ds)) > (lat*1000));
            endPoint = idx(1);
            
            allContra(g,accs+1,p1,:,:) = plotcontraData;
            allIpsi(g,accs+1,p1,:,:) = plotipsiData;
            subplot(1,length(exp.gapDur{exp.thisExp}),g)
            
            fill(plt.gapShades{exp.thisExp}(g+1,:)*ds,[ylims, flip(ylims)], thisCol, 'LineStyle','none','facealpha', 0.2/4); %Plot pulses in time
            hold on
            
            diffData = squeeze(nanmean(allContra(g,accs+1,p1,:,:),2) - nanmean(allIpsi(g,accs+1,p1,:,:),2));
            alldiffData(g,accs+1,p1,:,:) = diffData;
            
            stdshade(diffData(1:endPoint,:)', plt.alpha, thisCol, thisLinetype,[],plt.smooth) %Plot ERP traces

            hold on;
            title(['Gap: ' num2str(exp.gapDur{exp.thisExp}(g)) 's'])
            
            thisCol = plt.colors(1);
            if p1 == 1
                fill(plt.gapShades{exp.thisExp}(1,:)*ds,[ylims, flip(ylims)], thisCol{1}, 'LineStyle','none','facealpha', 0.2); %Plot pulses in time
            end
            ylim(ylims);
            hold on;
            xline(plt.zero(e)*exp.srate*ds);
            yline(0, 'LineStyle', '--');
            
            xticks(round(plt.xticks(e,:)*ds));
            xticklabels(plt.xlabs(e,:)); set(gca, 'fontsize', plt.font)
            xlabel(plt.xlabel); if g == 1; ylabel(['MBL ' unitTag]); end
            clear thisData;
            xlim([1,200*2.2])
        end
    end
    f.Units = 'centimeters';
    f.OuterPosition = [0.25 0.35 25  9];
end

if plt.rl == 0
    exportgraphics(f, [exp.figPath 'GT_FigS3_MBL_Exp' num2str(exp.thisExp) '.tiff'], 'resolution', 600);
    exportgraphics(f, [exp.figPath 'GT_FigS3_MBL_Exp' num2str(exp.thisExp) '.pdf'], 'resolution', 600);
end

%% Plot average of all gaps
% DB transform data if necessary
if plt.rl == 1
    f = figure;
    excursion = []; preP1 = []; postP1 = []; allCols = [];
    gIdx = 1:length(exp.gapDur{exp.thisExp});
    for p1 = 1:2
        thisCol = (plt.satPalettes{1,p1*p1});
        allCols = [allCols; thisCol];
        subplot(1,2,1); %Beta lateralisation after P2, sorted by gap and condition
        
        if plt.ci == 0
            plotData = squeeze(nanmean(alldiffData(:,1,p1,:,:),1));
            stdshade(plotData', plt.alpha, thisCol, thisLinetype,[],plt.smooth) %Plot ERP traces
            ylabel(['MBL ' unitTag]);
        else
            stdshade(squeeze(nanmean(allIpsi(:,1,p1,:,:),1))', plt.alpha, thisCol, ':',[],plt.smooth) %Plot ERP traces
            hold on;
            stdshade(squeeze(nanmean(allContra(:,1,p1,:,:),1))', plt.alpha, thisCol, '-',[],plt.smooth) %Plot ERP traces
            hold on;
            ylabel(['Motor beta ' unitTag]);
        end
        hold on;
        
        xline(plt.zero(1)*exp.srate*ds);
        yline(0, 'LineStyle', '--');
        
        xticks(round(plt.xticks(1,1:2)*ds));
        xticklabels(plt.xlabs(1,1:2)); set(gca, 'fontsize', plt.font)
        xlabel(plt.xlabel);
        xlim([1,1.2*200]); ylim([-2.5 1])
        set(gca, 'fontsize', 11)
        box off
        % Load and plot response-aligned
        [allIpsi_RL, allContra_RL, allMBL_RL] = GT_R1_MBL_RL(exp, plt, EEG,allP, baselinesSL);
        
        subplot(1,2,2);
        hold on;
        if plt.ci == 0
            plotData = squeeze(nanmean(allMBL_RL(:,p1,:,:),1));
            stdshade(plotData', plt.alpha, thisCol, thisLinetype,[],plt.smooth) %Plot ERP traces
            ylabel('');
        else
            stdshade(squeeze(nanmean(allIpsi_RL(:,p1,:,:),1))', plt.alpha, thisCol, ':',[],plt.smooth) %Plot ERP traces
            hold on;
            stdshade(squeeze(nanmean(allContra_RL(:,p1,:,:),1))', plt.alpha, thisCol, '-',[],plt.smooth) %Plot ERP traces
            hold on;
            ylabel('');
        end
        
        
        hold on;
        xline(plt.zero(2)*exp.srate*ds);
        yline(0, 'LineStyle', '--');
        
        %     xlim(plt.xlims(2,:)*ds)
        xticks(round(plt.xticks(2,:)*ds));
        xticklabels(plt.xlabs(2,:));
        yticklabels([]);
        set(gca, 'fontsize', plt.font)
        xlabel(plt.xlabel); %ylabel([{'MBL ' unitTag}]);
        %     xlim([1,1.2*200]);
        ylim([-2.5 1])
        set(gca, 'fontsize', 11)
        box off;
        ax1 = gca;
        yruler = ax1.YRuler;
        yruler.Axle.Visible = 'off';
        
    end
    f.Units = 'centimeters';
    f.OuterPosition = [5 5 14 8];
end
