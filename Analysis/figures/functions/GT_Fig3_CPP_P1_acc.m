function [data_allC, dataTopo] = GT_Fig4_CPP_P1_acc(exp, plt, EEG,RIDE,allBest,allP)
sr = exp.srate;
ds = 1 %Set to 1 if data not downsampled, set to (200/512) if ds to 200;
ylims = plt.ylims; 
lineType = '-';
e = 1;
chans = exp.cppChans;

f = figure
data_allC = []; data = []; dataTopo = [];
for g = 1:length(exp.gapDur{exp.thisExp}) %Longest gap only!
    subplot(1,length(exp.gapDur{exp.thisExp}), g)
    for p1 = 1:2 %Longest gap only!
            thisGrouping= 'np_gap_p1coh_acc_';
%             thisGrouping= 'np_gap_p1coh_';

            thisGap = exp.gapLabs{exp.thisExp}{g};
            thisP1 = exp.cohLabs{p1};
            thisAcc = exp.accLabs{plt.acc}; 
            thisCol = flip(plt.satPalettes{1,p1*p1});
            
            dataTag = ['twoPulse' thisGap thisP1 thisAcc];
            incTag = ['twoPulse' thisGap thisP1 thisAcc]
%             include = allP.(plt.epoch).(thisGrouping).(incTag).nTrials > 5;
            include = logical(repmat(1,1,size(allP.(plt.epoch).(thisGrouping).(incTag).sl,3))) %allP.(plt.epoch).(thisGrouping).(incTag).nTrials > 5;

                        
            rideTag = [thisGap, '_p1' thisP1];
            baselines = nanmean(allP.(plt.epoch).(thisGrouping).(dataTag).sl(:,round(1.3*sr):round(1.5*sr),:),2);
            data =((allP.(plt.epoch).(thisGrouping).(dataTag).sl))-baselines;
            
            %From 0.4 to 0.6s post P1
            dataTopo(g,p1,:,:) = squeeze(nanmean(data(:,round(1.9*sr):round(2.1*sr),:),2)); 
            
            for s=  1:size(data,3)
                
                    if plt.opti == 1; bestCh = allBest(s,:);
                    else bestCh = exp.cppChans; end
                newData(:,s) = squeeze(nanmean(data([bestCh],:,s),1));
            end
            data = []; data = newData; 
            if plt.RIDE ==1
                indices = round(1.3*512:4.499*512);
                
                for id = 1:size(RIDE.(rideTag).s1,1)
                     if plt.opti == 1; bestCh = allBest(id,:);
                    else bestCh = exp.cppChans; end
                    thisStim(id,:) = nanmean(RIDE.(rideTag).s1{id}(:,bestCh),2);
                end
                data(indices,:) = data(indices,:) - thisStim';
            end
                        
           %Low pass filter at 6Hz
            [bfilt,afilt] = butter(3, 6*2/512, 'low');   % lo-pass
            for s = 1:size(data,2) %Subject
                nonan = find(isnan(data(:,s)), 1, 'first'); 
                if isempty(nonan);
                    data(:,s) = filtfilt(bfilt,afilt, data(:,s));
                elseif nonan ~= 1
                    data(1:nonan-1,s) = filtfilt(bfilt,afilt, data(1:nonan-1,s));
                    data(nonan:end,s) = NaN;
                elseif nonan == 1
                    data(nonan:end,s) = NaN;
                end
            end
            
            data_allC(g,p1,:,:) = data; 
             %Remove conditions with <5 trials
            data_allC(g,p1,:,include == 0) = NaN; 
            
            %Plot only up to the earliest possible response time (lead-out dot offset)
            endPoint = round([0.4 + exp.gapDur{exp.thisExp}(g) + 0.75+1.5]*sr);
            
            fill(plt.gapShades{exp.thisExp}(g+1,:)*ds,[ylims, flip(ylims)], thisCol, 'LineStyle','none','facealpha', 0.2/4); %Plot pulses in time
            hold on
            stdshade(data(1.3*sr:endPoint,:)', plt.alpha, thisCol, lineType,[],plt.smooth) %Plot ERP traces
            hold on;
            
            thisCol = plt.colors(1);
            fill(plt.gapShades{exp.thisExp}(1,:)*ds,[ylims, flip(ylims)], thisCol{1}, 'LineStyle','none','facealpha', 0.2/4);  %Plot pulses in time
            ylim(ylims);hold on;xline(plt.zero(e)*exp.srate*ds);yline(0, 'LineStyle', '--');
            
            xticks(round(plt.xticks(e,:)*ds));
            xticklabels(plt.xlabs(e,:)); set(gca, 'fontsize', plt.font)
            clear data;
            xlim([1, endPoint-1.3*sr])
            if g == 1; ylabel('CPP (\muV/m^{2})'); end
    end
    title(['Gap: ' num2str(exp.gapDur{exp.thisExp}(g)) 's'])
    xlabel('Time (s)')
    xlim([1,512*2.2])

end

f.Units = 'centimeters';
f.OuterPosition = [0.25 0.35 25 9];

exportgraphics(f,[exp.figPath, 'GT_FigS3_acc_Exp' num2str(exp.thisExp) '.tiff'], 'Resolution', 600) %Save CPP panels of Fig S3
exportgraphics(f,[exp.figPath, 'GT_FigS3_acc_Exp' num2str(exp.thisExp) '.pdf'], 'Resolution', 600) %Save CPP panels of Fig S3

%% Plot averaged across gaps 
f = figure; 
f.Units = 'centimeters';
f.OuterPosition = [0.5,0.5,8,14];

c = 5;
slope = []; peak = []; 

meanAmp = []; slope = []; 
gapIdx = 1:length(exp.gapDur{exp.thisExp});
plt.box = 1;
ylims = [-5,15];
if plt.acc == 2; lineType = '-'; else, lineType = ':'; end

for p1 = 1:2
    c = c-1;
    thisCol = flip(plt.satPalettes{1,p1*p1});
    
    baselines =  squeeze(nanmean(nanmean(data_allC(gapIdx,p1,1.45*sr:1.5*sr,:),3),1));
    hold on; 
    stdshade(squeeze(nanmean(data_allC(gapIdx,p1,1.3*sr:2.3*sr,:),1))'-baselines, plt.alpha, thisCol, lineType,[],plt.smooth) %Plot ERP traces

    
end
ylim([ylims])
set(gca, 'fontsize', plt.font)
xticklabels({'P2', '0.5', '1'})
xticks([[0.2,0.7,1.2]*512])
xline(0.2*512)
ylabel(plt.ERP.ylabel);
f.Units = 'centimeters';
xlabel('Time (s)'); 
box on
f.Units = 'centimeters';
f.OuterPosition = [0.5,0.5,6.5,9];
exportgraphics(f,[exp.figPath, 'GT_FigS6_P1_Exp' num2str(exp.thisExp) '_Acc' num2str(plt.acc) '.tiff'], 'Resolution', 600) %Save CPP panels of Fig S3
exportgraphics(f,[exp.figPath, 'GT_FigS6_P1_Exp' num2str(exp.thisExp) '_Acc' num2str(plt.acc) '.pdf'], 'Resolution', 600) %Save CPP panels of Fig S3

