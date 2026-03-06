function [alldiffData_P2,allsumData_P2, dataTopo] = GT_Fig4_MBL_P2(exp, plt, EEG,RIDE,allBest,allP)

%% Figure 1: plot beta lateralisation stimulus locked
factor = 2; e= 1;


ds = (200/512);
clear diffData;
xlims = plt.xlims(e,1)*200/512:plt.xlims(e,2)*ds;
plt.epoch = 'P1'; e = 1;
plt.baseline = 0;
par.db = 1;

f = figure

subplot(5,4,[9:11,13:15,17:19]);

plt.font = 11
xline(200, 'k');
xline([round(0.8*200),round(0.9*200)], '--k');

chansel = plt.chansel;

allData = [];
allcontraData = []; allipsiData = [];
allContra = []; allIpsi = [];
allDiff = []; dataTopo = [];
alldiffData = [];
for g = 1:length(exp.gapDur{exp.thisExp})
    c = 5;
    for p1 = 1:2
        for p2 = 1:2
            c = c-1;
            for d = 1:length(exp.dirLabs)
                ylims = [-15 5];%
                
                thisLinetype = '-';
                thisGrouping= 'np_gap_p1p2coh_dir_';
                thisGap = exp.gapLabs{exp.thisExp}{g};
                thisP1 = exp.cohLabs{p1};
                thisP2 = exp.cohLabs{p2};
                thisDir = exp.dirLabs{d};
                
                thisCol = hex2rgb(plt.Rcolors{c});
                
                
                %Compute lateralisation based on cluster of left/right hemisphere electrodes
%                 if strcmp(thisDir, 'dotLeft'); ipsiChan = exp.movCluster(1,chansel); contraChan = exp.movCluster(2,chansel); otherDir= 'dotRight';
%                 else, ipsiChan = exp.movCluster(2,chansel); contraChan = exp.movCluster(1,chansel); otherDir = 'dotLeft'; end
%                 
            if strcmp(thisDir, 'dotLeft'); ipsiChan = 1; contraChan =2; otherDir= 'dotRight';
            else, ipsiChan = 2; contraChan = 1; otherDir = 'dotLeft'; end

                dataTag = ['twoPulse' thisGap thisP1 '_' thisP2 thisDir]; %IPSILATERAL
                if par.db == 1; unitTag = '(dB)'; ylims = [-1.5 1];
                    baselines = nanmean(allP.(plt.epoch).(thisGrouping).(dataTag).sl(:,[1.3*200:1.5*200],:).^factor,2);
                    thisData = 10*log10(allP.(plt.epoch).(thisGrouping).(dataTag).sl.^factor./baselines);
                    ipsiData = nanmean(thisData(ipsiChan,xlims,:),1)%.^factor; %IPSI channels
                    contraData = nanmean(thisData(contraChan,xlims,:),1)%.^factor; %IPSI channels
                else, unitTag = '(muV/m^{2})';ylims = [-15 5];
                    thisData = allP.(plt.epoch).(thisGrouping).(dataTag).sl;
                    ipsiData = nanmean(thisData(ipsiChan,xlims,:),1).^factor; %IPSI channels
                    contraData = nanmean(thisData(contraChan,xlims,:),1).^factor; %CONTRA channels
                end
                if plt.baseline == 1
                    baseTag = 'base';
                    allipsiData(d,:,:) = squeeze(ipsiData - nanmean(nanmean(thisData(ipsiChan,260:300,:),1),2));
                    allcontraData(d,:,:) = squeeze(contraData - nanmean(nanmean(thisData(contraChan,260:300,:),1),2));
                elseif plt.baseline == 2 %as a function of preResponse amplitude
                    %                         thisRespThresh = nanmean(allP.(plt.epoch).(thisGrouping).(dataTag).rl(:,0.7*200:0.8*200,:),2);
                    %                         thisRespThresh= permute(thisRespThresh,[2,1,3]);
                    
                    respThresh = [(nanmean(allP.P1.allCond_rl(:,0.7*200:0.8*200,:),2))];
                    respThresh = permute(respThresh,[2,1,3]);
                elseif plt.baseline == 0
                    allipsiData(d,:,:) = squeeze(ipsiData);
                    allcontraData(d,:,:) = squeeze(contraData);
                end
                
                
                if plt.onePulse ==1
                    indices = round(1.3*200:2.9*200); %2.9*512
                    subtTag = ['twoPulse' exp.gapLabs{exp.thisExp}{end} thisP1 thisDir]
                    if par.db == 1; unitTag = '(dB)'; ylims = [-1.5 1];
                        baselines = nanmean(allP.(plt.epoch).np_gap_p1coh_dir_.(subtTag).sl(:,[1.3*200:1.5*200],:).^factor,2);
                        toSubtract = 10*log10(allP.(plt.epoch).np_gap_p1coh_dir_.(subtTag).sl.^factor./baselines);
                    else
                        toSubtract = allP.(plt.epoch).np_gap_p1coh_dir_.(subtTag).sl;
                    end
                    for id = 1:size(allP.P1.np_gap_p1coh_dir_.(subtTag).sl,3)
                        thisStim_ipsi(id,:) = nanmean(toSubtract(ipsiChan,indices,id),1);
                        thisStim_contra(id,:) = nanmean(toSubtract(contraChan,indices,id),1);
                    end
                    idx =[1:321];
                    allipsiData(d,idx,:) = squeeze(allipsiData(d,idx,:)) - thisStim_ipsi';
                    allcontraData(d,idx,:) = squeeze(allcontraData(d,idx,:)) - thisStim_contra';
                end
            
            end
%             respData(d,:,:,:) = squeeze(allP.(plt.epoch).(thisGrouping).(dataTag).rl(:,:,:) - nanmean(allP.(plt.epoch).(thisGrouping).(dataTag).sl(:,1.3*200:1.5*200,:),2));
            
            plotcontraData = squeeze(nanmean(allcontraData));
            plotipsiData= squeeze(nanmean(allipsiData));
%             plotrespData = squeeze(nanmean(respData));
            
%             dataTopo(g,p1,p2,:,:) = squeeze(nanmean(respData(2,:,0.1*200:0.2*200,:) - respData(1,:,0.1*200:0.2*200,:),3));
            
            lat = [0.4 + exp.gapDur{exp.thisExp}(g) + 0.75]; %dur of two pulses, plus dur of gap, plus minimum lead out duration
            idx = find(EEG.times(round(xlims/ds)) > (lat*1000));
            endPoint = idx(1)
            
            allContra(g,p1,p2,:,:) = plotcontraData;
            allIpsi(g,p1,p2,:,:) = plotipsiData;
            
            alldiffData(g,p1,p2,:,:) = plotcontraData - plotipsiData;
            
            %Align to P2
            p2Idx = [0.2+exp.gapDur{exp.thisExp}(g)+0.2]*200;
            alldiffData_P2(g,p1,p2,:,:) = plotcontraData(p2Idx-0.2*200:p2Idx+200,:) - plotipsiData(p2Idx-0.2*200:p2Idx+200,:);
            allsumData_P2(g,p1,p2,:,:) = plotcontraData(p2Idx-0.2*200:p2Idx+200,:) + plotipsiData(p2Idx-0.2*200:p2Idx+200,:); %Sum hemispheres

            diffData = squeeze(nanmean(allContra(g,p1,p2,:,:),2) - nanmean(allIpsi(g,p1,p2,:,:),2));
            
            %         %Low pass filter at 6Hz
            %         [bfilt,afilt] = butter(3, 6*2/512, 'low');   % lo-pass
            %         for s = 1:size(diffData,2) %Subject
            %             diffData(:,s) = filtfilt(bfilt,afilt, diffData(:,s));
            %         end
            
            subplot(1,4,g)
            fill(plt.gapShades{exp.thisExp}(g+1,:)*ds,[ylims, flip(ylims)], thisCol, 'LineStyle','none','facealpha', 0.2/4); %Plot pulses in time
            hold on
            stdshade(diffData(1:endPoint,:)', plt.alpha, thisCol, thisLinetype,[],plt.smooth) %Plot ERP traces
            hold on;
            
            thisCol = plt.colors(1)% plt.colors(7,2);
            if p1 == 1 & p2 == 1; fill(plt.gapShades{exp.thisExp}(1,:)*ds,[ylims, flip(ylims)], thisCol{1}, 'LineStyle','none','facealpha', 0.2);end  %Plot pulses in time
            ylim(ylims);
            hold on;
            xline(plt.zero(e)*exp.srate*ds);
            yline(0, 'LineStyle', '--');
            title(['Gap: ' num2str(exp.gapDur{exp.thisExp}(g)) 's'])
            
            %                 xline((plt.topoTimes+250)*exp.srate/1000)
            xticks(round(plt.xticks(e,:)*ds));
            xticklabels(plt.xlabs(e,:)); set(gca, 'fontsize', plt.font)
            xlabel(plt.xlabel); ylabel([{'MBL (\muV/m^{2})'}]);
        end
    end
end


f.Units = 'centimeters';
f.OuterPosition = [0.25 0.35 25  9];
% exportgraphics(f, [exp.figPath 'GT_FigS3_MBL_Exp' num2str(exp.thisExp) '.tiff'], 'resolution', 600); 
% exportgraphics(f, [exp.figPath 'GT_FigS3_MBL_Exp' num2str(exp.thisExp) '.pdf'], 'resolution', 600); 
%% Plot average of all gaps
% DB transform data if necessary
f = figure
c = 5;
allCols = [];
gIdx = 1:length(exp.gapDur{exp.thisExp});
for p1 = 1:2
    for p2 = 1:2
        c = c-1;
        thisCol = hex2rgb((plt.Rcolors{c}));
        allCols = [allCols; thisCol];
        plotData = squeeze(nanmean(alldiffData_P2(gIdx,p1,p2,:,:),1));
        plotData  = squeeze(plotData)% - nanmean(plotData(0.1*200:0.2*200,:),1));
        subplot(2,1,1)
        stdshade(plotData', plt.alpha, thisCol, '-',[],plt.smooth) %Plot ERP traces
        hold on;
        hold on;
        %        if p1 == 1;  fill(plt.gapShades{exp.thisExp}(1,:)*ds,[ylims, flip(ylims)], thisCol{1}, 'LineStyle','none','facealpha', 0.2); end %Plot pulses in time
        %         ylim(ylims);
        hold on;
        xline(plt.zero(e)*exp.srate*ds);
        yline(0, 'LineStyle', '--');
        
        %                 xline((plt.topoTimes+250)*exp.srate/1000)
        xticks(round(plt.xticks(e,1:2)*ds));
        xticklabels(plt.xlabs(e,1:2)); set(gca, 'fontsize', plt.font)
        xlabel(plt.xlabel); ylabel([{'MBL (\muV/m^{2})'}]);
        xlim([1,1.2*200]);
        set(gca, 'fontsize', 11)
        
        preP2(gIdx,:) =  squeeze(nanmean(alldiffData_P2(gIdx,p1,p2,0.1*200:0.2*200,:),4));
        for g = 1:length(exp.gapDur{exp.thisExp})
            postP3(g,:) = min(squeeze(alldiffData_P2(g,p1,p2,0.2*200:1*200,:)));
        end
        excursion(c,:) = nanmean(preP2(gIdx,:) - postP3(gIdx,:));
        
        if c == 1
            subplot(2,1,2); %Beta lateralisation after P2, sorted by gap and condition
            b = bar([1:4], flip(nanmean(excursion')))
            b.FaceColor = 'flat';
            b.CData(:,:) = (allCols);
            xticks([1:4]); xticklabels({'HH', 'HL', 'LH', 'LL'});
            ylabel(['\Delta MBL'])
            set(gca, 'fontsize', 11)
        end
    end
end
f.Units = 'centimeters';
f.OuterPosition = [5 5 8 14];
