function [data_allC, dataTopo] = GT_Fig4_CPP_P1(exp, plt, EEG,RIDE,allBest,allP)
sr = exp.srate;
ds = 1 %Set to 1 if data not downsampled, set to (200/512) if ds to 200;
ylims = plt.ylims; 
lineType = '-';
e = 1;
chans = exp.cppChans;

f = figure
data_allC = []; data = []; dataTopo = [];
for g = 1:length(exp.gapDur{exp.thisExp}) 
    subplot(1,length(exp.gapDur{exp.thisExp}), g)
    for p1 = 1:2 
            thisGrouping= 'np_gap_p1coh_';

            thisGap = exp.gapLabs{exp.thisExp}{g};
            thisP1 = exp.cohLabs{p1};
            
            thisCol = flip(plt.satPalettes{1,p1*p1});
            
            dataTag = ['twoPulse' thisGap thisP1]
            incTag = ['twoPulse' thisGap thisP1]
            include = logical(repmat(1,1,size(allP.(plt.epoch).(thisGrouping).(incTag).sl ,3))) 

                        
            rideTag = [thisGap, '_p1' thisP1];
            baselines = nanmean(allP.(plt.epoch).(thisGrouping).(dataTag).sl (:,round(1.3*sr):round(1.5*sr),:),2);
            data =((allP.(plt.epoch).(thisGrouping).(dataTag).sl ))-baselines;
            
            %From 0.4 to 0.6s post P1
            dataTopo(g,p1,:,:) = squeeze(nanmean(data(:,round(1.9*sr):round(2.1*sr),:),2)); 
            
            for s=  1:size(data,3)
                
                    if plt.opti == 1; bestCh = allBest(s,:);
                    else bestCh = exp.cppChans; end
                newData(:,s) = squeeze(nanmean(data([bestCh],:,s),1));
            end
            data = []; data = newData; 
            ufdata = data; 
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
            
              ignoredP(g) = NaN; 
        if any(isnan(data(1,:))), warning(['Oops! ' sum(isnan(data(1,:))) 'P excluded in G' num2str(g)]); 
        ignoredP(g) =  sum(isnan(data(1,:)));  end
    
            data_allC(g,p1,:,:) = data; 
             %Remove conditions with <5 trials
            data_allC(g,p1,:,include == 0) = NaN; 
            
            %Plot only up to the earliest possible response time (lead-out dot offset)
            endPoint = round([0.4 + exp.gapDur{exp.thisExp}(g) + 0.75+1.5]*sr);
            
            fill(plt.gapShades{exp.thisExp}(g+1,:)*ds,[ylims, flip(ylims)], thisCol, 'LineStyle','none','facealpha', 0.2/4); %Plot pulses in time
            hold on
            stdshade(ufdata(1.3*sr:endPoint,:)', plt.alpha, thisCol, lineType,[],plt.smooth) %Plot ERP traces
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

exportgraphics(f,[exp.figPath, 'GT_FigS3_Exp' num2str(exp.thisExp) '.tiff'], 'Resolution', 600) %Save CPP panels of Fig S3
exportgraphics(f,[exp.figPath, 'GT_FigS3_Exp' num2str(exp.thisExp) '.pdf'], 'Resolution', 600) %Save CPP panels of Fig S3

%% Additional plots (R1)
% %% CPP-P1
f = figure; 
f.Units = 'centimeters';
f.OuterPosition = [0.5,0.5,8,14];
tiledlayout(2,2); 
nexttile([1,2]);
peak = []; slope = []; 
clear idx; 
plt.box = 1;
lineType = '-';
for p1 = 1:2
    thisCol = flip(plt.satPalettes{1,p1*p1});
    stdshade(squeeze(nanmean(data_allC(:,p1,round(1.3*sr):round(2.3*sr),:)))', plt.alpha, thisCol, lineType,[],plt.smooth) %Plot ERP traces
    hold on;
    %Compute slopes & max amp. by subject
    for i = 1:size(data_allC,4)
        onset = round(0.195*sr); %Start of fit
        p1peak_idx= round(0.4*sr):round(0.6*sr); %GA peak amplitude to search for individual subject peaks
        offset = p1peak_idx(1);
        y = squeeze(nanmean(data_allC(:,p1,round(1.5*sr):round(2.3*sr), i)));
        x = 1:length(y);
        [peak(p1,i), idx(p1,i)] = max(y(p1peak_idx));
        B = glmfit(x(onset:idx(p1,i)+offset),y(onset:idx(p1,i)+offset));
        slope(p1,i) = B(2);
    end
end
set(gca, 'fontsize', plt.font)
xticklabels({'P1', '0.5', '1'})
xticks([[0.2,0.7,1.2]*512])
xline(0.2*512)
ylabel(plt.ERP.ylabel);
xlabel('Time (s)'); 

allCols = [flip(plt.satPalettes{1,1*1}); flip(plt.satPalettes{1,2*2})];

%Plot SLOPES 
nexttile(3);
x = [1,2]; y = [nanmean(slope(1,:)), nanmean(slope(2,:))]; 
error = [nanstd(slope(1,:))/sqrt(length(slope(1,:))), nanstd(slope(2,:))/sqrt(length(slope(2,:)))]; 
if plt.box == 1
    b = boxplot(slope(:,:)', 'BoxStyle','filled', 'Colors', allCols);
    hold on; scaling = 1;
    plot([1,2],(slope(:,:)'), 'Color', [0.7 0.7 0.7]);
else
    b = bar(x,y);
    b.FaceColor = 'flat';
    b.CData(:,:) = allCols;
    hold on; scaling = 1.75;
    errorbar(x, y, error, '.k');
end

xticks(x); xticklabels({'H', 'L'});
set(gca, 'FontSize', plt.font)
ylabel('Slope')
[h,p,CI,slopeStats] = ttest(slope(1,:), slope(2,:)); 

yl = get(gca(),'Ylim');
yl(2) = yl(2) +10*[yl(2)-yl(1)]/100;
ylim([yl])
xl = get(gca(),'Xlim')+scaling;
yl = get(gca(),'Ylim');
yl = yl -5*[yl(2)-yl(1)]/100;

textp = num2str(p);
if p < 0.001; thisTitle = '***'; fw = 'bold';
elseif p<0.01; thisTitle = '**'; fw = 'bold';
elseif p < 0.05; thisTitle = '*'; fw = 'bold';
else; thisTitle = textp(1:4); fw = 'normal'; end
t = text(xl(1),yl(2),thisTitle, ...
    'VerticalAlignment','top', ...
    'HorizontalAlignment', 'center',...
    'FontSize',12);
% 
%Plot AMP
nexttile(4);
x = [1,2]; y = [nanmean(peak(1,:)), nanmean(peak(2,:))]; 
error = [std(peak(1,:))/sqrt(length(peak(1,:))), std(peak(2,:))/sqrt(length(peak(1,:)))]; 
if plt.box == 1
    b = boxplot(peak(:,:)', 'BoxStyle','filled', 'Colors', allCols);
    hold on;
        scaling = 1;

    plot([1:2],(peak(:,:)'), 'Color', [0.7 0.7 0.7]);
else
    b = bar(x,y);
    b.FaceColor = 'flat';
    b.CData(:,:) = allCols;
    hold on;
    scaling = 1.75;
    errorbar(x, y, error, '.k');
end

xticks(x); xticklabels({'H', 'L'});
set(gca, 'FontSize', plt.font)
[h,p,CI,ampStats] = ttest(peak(1,:), peak(2,:)); 
ylabel('Max. amp (\muV/m^{2})')
% yticks('');
yl = get(gca(),'Ylim');
yl(2) = yl(2) +10*[yl(2)-yl(1)]/100;
ylim([yl])
xl = get(gca(),'Xlim')+scaling;
yl = get(gca(),'Ylim');
yl = yl -5*[yl(2)-yl(1)]/100;
  textp = num2str(p);
    if p < 0.001; thisTitle = '***'; fw = 'bold'; 
    elseif p<0.01; thisTitle = '**'; fw = 'bold'; 
    elseif p < 0.05; thisTitle = '*'; fw = 'bold'; 
    else; thisTitle = textp(1:4); fw = 'normal'; end
    
    t = text(xl(1),yl(2),thisTitle, ...
        'VerticalAlignment','top', ...
        'HorizontalAlignment', 'center',...
        'FontSize',12);
set(gca, 'FontSize', plt.font)

axes('position',[.25 .75 .2 .2]);
diffTopo = squeeze(nanmean(dataTopo(:,1,:,:,:),1) - nanmean(dataTopo(:,2,:,:,:),1)); 

GT_topoplot(nanmean(diffTopo,2),EEG.chanlocs, ...
    'numcontour', 0, ...
    'plotrad',0.65, ...
    'headrad', 0.5,...
    'emarker2', {chans,'*','k',1,1}); % To plot location of the selected channel

exportgraphics(f,[exp.figPath, 'GT_FigR1_Exp' num2str(exp.thisExp) '.tiff'], 'Resolution', 600) %Save CPP panels of Fig S3
exportgraphics(f,[exp.figPath, 'GT_FigR1_Exp' num2str(exp.thisExp) '.pdf'], 'Resolution', 600) %Save CPP panels of Fig S3

% %% Plot averaged across gaps 
% f = figure; 
% f.Units = 'centimeters';
% f.OuterPosition = [0.5,0.5,6.5,9];
% peak = []; slope = []; 
% clear idx; 
% plt.box = 1;
% lineType = '-';
% 
% for p1 = 1:2
%     thisCol = flip(plt.satPalettes{1,p1*p1});
%     stdshade(squeeze(nanmean(data_allC(:,p1,round(1.3*sr):round(2.3*sr),:)))', plt.alpha, thisCol, lineType,[],plt.smooth) %Plot ERP traces
%     hold on;
%     %Compute slopes & max amp. by subject
%     for i = 1:size(data_allC,4)
%         y = squeeze(nanmean(data_allC(:,p1,round(1.5*sr):round(2.3*sr), i)));
%         x = 1:length(y);
%         [peak(p1,i), idx(p1,i)] = max(y);
%         B = glmfit(x(1:idx(p1,i)),y(1:idx(p1,i)));
%         slope(p1,i) = B(2);
%     end
% end
% set(gca, 'fontsize', plt.font)
% xticklabels({'P1', '0.5', '1'})
% xticks([[0.2,0.7,1.2]*512])
% xline(0.2*512)
% ylabel(plt.ERP.ylabel);
% xlabel('Time (s)'); 
