function [data_allC_p2, dataTopo,include] = GT_Fig4_CPP_P2(exp, plt, EEG,RIDE,allBest,allP)
sr = exp.srate;
ds = 1 %Set to 1 if data not downsampled, set to (200/512) if ds to 200;
ylims = plt.ylims; 
lineType = '-';
e = 1;
chans = exp.cppChans;

f = figure
f.Units = 'centimeters';
f.OuterPosition = [0.25 0.35 25 9];
data_allC = []; dataTopo = [];


for g = 1:length(exp.gapDur{exp.thisExp}) %All gaps
    subplot(1,length(exp.gapDur{exp.thisExp}), g)
    c = 5;
    for p1 = 1:2 
        for p2 = 1:2
            c = c-1;
            thisGrouping= 'np_gap_p1p2coh_';

            thisGap = exp.gapLabs{exp.thisExp}{g};
            thisP1 = exp.cohLabs{p1};
            thisP2 = exp.cohLabs{p2};

            lineType = '-';
            thisCol = hex2rgb((plt.Rcolors{c}));
            dataTag = ['twoPulse' thisGap thisP1 '_' thisP2]
            incTag = ['twoPulse' thisGap thisP1 '_' thisP2]
            
            include = logical(repmat(1,1,size(allP.(plt.epoch).(thisGrouping).(incTag).sl ,3))) %allP.(plt.epoch).(thisGrouping).(incTag).nTrials > 5;

            disp(['Included participants: ' num2str(sum(include))]);
            if sum(include == 0) > 0
                warning('oops!')
            end
            rideTag = [thisGap, '_p1' thisP1 '_p2' thisP2];
            baselines = 0%nanmean(allP.(plt.epoch).(thisGrouping).(dataTag).sl(:,round(1.3*sr):round(1.5*sr),:),2);
            data =((allP.(plt.epoch).(thisGrouping).(dataTag).sl ))-baselines;
            
            p2Idx = (1.5 + exp.gapDur{exp.thisExp}(g))*sr;
            %From 0.4 to 0.6s post P2
            dataTopo(g,p1,p2,:,:,:) = squeeze(nanmean(data(:,p2Idx+0.7*sr:p2Idx+0.9*sr,:),2));% -  squeeze(nanmean(data(:,p2Idx:p2Idx+0.2*sr,:),2));
            for t = 1:5
                sPoint = p2Idx-0.2*(t)*sr; 
                ePoint = sPoint + 0.1*sr;
                topoPreP2(g,p1,p2,:,:,:,t) =squeeze(nanmean(data(:,sPoint:ePoint,:),2));% -  squeeze(nanmean(data(:,p2Idx:p2Idx+0.2*sr,:),2));
            end
    
            for s = 1:size(data,3)
                if plt.opti == 1; bestCh = allBest(s,:); 
                else bestCh = exp.cppChans; end
                newData(:,s) = squeeze(nanmean(data([bestCh],:,s),1));
            end
            data = []; data = newData;
            
                                  clear thisStim; 

            if plt.RIDE ==1
                indices = round(1.3*512:4.499*512);
                
                for id = 1:size(RIDE.(rideTag).s1,1)
                    if plt.opti == 1; bestCh = allBest(id,:);
                    else bestCh = exp.cppChans; end
                    thisStim(id,:) = nanmean(RIDE.(rideTag).s1{id}(:,[bestCh]),2);
                end
                data(indices,:) = data(indices,:) - thisStim';
            end
            clear thisStim; 
            
             if plt.onePulse ==1
                indices = round(1.3*512:2.9*512); %2.9*512
                subtTag = ['twoPulse' exp.gapLabs{exp.thisExp}{end} thisP1]
                toSubtract = allP.(plt.epoch).np_gap_p1coh_.(subtTag).sl ; 
               for id = 1:size(toSubtract,3)
                    if plt.opti == 1; bestCh = allBest(id,:);
                    else bestCh = exp.cppChans; end
                    thisStim(id,:) = nanmean(toSubtract(bestCh,indices,id),1);
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

            data_allC(g,p1,p2,:,:) = data; 
            data_allC_p2(g,p1,p2,:,:) = data(p2Idx:p2Idx+1.2*sr,:);

            %Remove conditions with <5 trials
            data_allC(g,p1,p2,:,include == 0) = NaN; 
            data_allC_p2(g,p1,p2,:,include == 0) = NaN;

            endPoint = round([0.4 + exp.gapDur{exp.thisExp}(g) + 0.75+1.5]*sr);
            
            fill(plt.gapShades{exp.thisExp}(g+1,:)*ds,[ylims, flip(ylims)], thisCol, 'LineStyle','none','facealpha', 0.2/4); %Plot pulses in time
            hold on
            stdshade(data(1.3*sr:endPoint,include)', plt.alpha, thisCol, lineType,[],plt.smooth) %Plot ERP traces
            hold on;
            
            hold on
            thisCol = plt.colors(1);
            fill(plt.gapShades{exp.thisExp}(1,:)*ds,[ylims, flip(ylims)], thisCol{1}, 'LineStyle','none','facealpha', 0.2/4);  %Plot pulses in time
            ylim(ylims);
            hold on;
            xline(plt.zero(e)*exp.srate*ds);
            yline(0, 'LineStyle', '--');
            
            xticks(round(plt.xticks(e,:)*ds));
            xticklabels(plt.xlabs(e,:)); set(gca, 'fontsize', plt.font)
            clear data;
            xlim([1, endPoint-1.3*sr])
            if g == 1; ylabel('CPP (\muV/m^{2})'); end
        end
    end
    xlabel('Time (s)')
    xlim([1,512*2.2])
    ylim([-5,25])
    title(['Gap: ' num2str(exp.gapDur{exp.thisExp}(g)) 's'])
end
exportgraphics(f,[exp.figPath, 'GT_FigS4_Exp' num2str(exp.thisExp) '_Clean' num2str(plt.onePulse) '.tiff'], 'Resolution', 600) %Save CPP panels of Fig S3
exportgraphics(f,[exp.figPath, 'GT_FigS4_Exp' num2str(exp.thisExp) '_Clean' num2str(plt.onePulse) '.pdf'], 'Resolution', 600) %Save CPP panels of Fig S3



%% Plot averaged across gaps 
f = figure; 
f.Units = 'centimeters';
f.OuterPosition = [0.5,0.5,8,14];

c = 5;
slope = []; peak = []; 
tiledlayout(2,2); 
nexttile([1,2]);
meanAmp = []; slope = []; 
gapIdx = 2:length(exp.gapDur{exp.thisExp});
plt.box = 1;
ylims = [-5,35];
lineType = '-';

for p1 = 1:2
    for p2 = 1:2
        c = c-1;
        thisCol = hex2rgb((plt.Rcolors{c}));
        baselines =  squeeze(nanmean(nanmean(data_allC_p2(gapIdx,p1,p2,0.15*sr:0.2*sr,:),4),1));

        stdshade(squeeze(nanmean(data_allC_p2(gapIdx,p1,p2,:,:),1))'-baselines, plt.alpha, thisCol, lineType,[],plt.smooth) %Plot ERP traces
        hold on;
        
          %Compute slopes & max amp. by subject 
        for i = 1:size(data_allC,5)
            y = squeeze(nanmean(data_allC_p2(gapIdx,p1,p2,:, i),1)-baselines(i));
            x = 1:length(y); 
            [peak(p1,p2,i), idx(p1,p2,i)] = max(y(0.6*sr:0.8*sr));
            B = glmfit(x,y); 
            meanAmp(p1,p2,i) = peak(p1,p2,i);%nanmean(y(0.6*sr:0.8*sr)); %

            slope(p1,p2,i) = B(2);
        end
        
    end
end
ylim([ylims])
set(gca, 'fontsize', plt.font)
xticklabels({'P2', '0.5', '1'})
xticks([[0.2,0.7,1.2]*512])
xline(0.2*512)
ylabel(plt.ERP.ylabel);
f.Units = 'centimeters';
xlabel('Time (s)'); 


%Plot SLOPES 
nexttile(3)
x = [1:4];
y = []; error = []; 
c =0;
for p1 = 1:2
    for p2 = 1:2
        c = c+1; 
        allCols(c,1:3) = hex2rgb((plt.Rcolors{c})); 
        y = [y; nanmean(slope(p1, p2,:))];
        allSlopes(c,:) = slope(p1,p2,:);

        error = [error; std(slope(p1, p2,:))/sqrt(length(slope(p1, p2,:)))];
    end
end

if plt.box == 1
    b = boxplot(allSlopes(:,:)', 'BoxStyle','filled', 'Colors', flip(allCols));
      hold on;
else
    b = bar(x,y);
    hold on;
    errorbar(x, y, error, '.k');
    xticks(x); xticklabels({}); %{'HH', 'HL', 'LH', 'LL'});
    b.FaceColor = 'flat';
    b.CData(:,:) = flip(allCols);
end
xticks(x); xticklabels({}); %{'HH', 'HL', 'LH', 'LL'});
set(gca, 'FontSize', plt.font)
ylabel('Slope')
% yticks([0,0.01,0.02]);
    yl = get(gca(),'Ylim');
    yl(2) = yl(2) +20*[yl(2)-yl(1)]/100;
    ylim([yl])
    xl = get(gca(),'Xlim')+1;
    yl = get(gca(),'Ylim');
    yl = yl -5*[yl(2)-yl(1)]/100;
    

% Test for main effects
tSlopes = array2table(allSlopes');
tSlopes.Properties.VariableNames = {'HH','HL','LH','LL'};
WithinStructure = table([1 1 2 2]',[1 2 1 2]','VariableNames',{'P1','P2'});
WithinStructure.P1 = categorical(WithinStructure.P1);
WithinStructure.P2 = categorical(WithinStructure.P2);

rm = fitrm(tSlopes, 'HH,HL,LH,LL~1','WithinDesign',WithinStructure);
ranovatable = ranova(rm,'WithinModel','P1*P2');
ranovatable

effects = {' p1:', ' p2:'}
for pulse = 1:2
    p = ranovatable(3+2*(pulse-1),5).pValue;
    textp = num2str(p);
    if p < 0.001; thisTitle = '***'; fw = 'bold'; 
    elseif p<0.01; thisTitle = '**'; fw = 'bold'; 
    elseif p < 0.05; thisTitle = '*'; fw = 'bold'; 
    else; thisTitle = textp(1:4); fw = 'normal'; end
  
    yl = get(gca(),'Ylim');
    xl = get(gca(),'Xlim');

    t = text(xl(1),yl(2)-10*(pulse-1)*yl(2)/100,[effects{pulse}, thisTitle], ...
        'VerticalAlignment','top', ...
        'HorizontalAlignment', 'left',...
        'FontSize',10, 'FontWeight', fw);

    hold on
end


%Plot AMP
nexttile(4)
x = [1:4];
y = []; error = []; 
c = 0; 
allCols = [];
for p1 = 1:2
    for p2 = 1:2
        c = c+1; 
        allCols(c,1:3) = hex2rgb((plt.Rcolors{c})); 
        y = [y; nanmean(meanAmp(p1, p2,:),3)];
        allAmp(c,:) = meanAmp(p1, p2,:);
        error = [error; std(meanAmp(p1, p2,:))/sqrt(length(meanAmp(p1, p2,:)))];
    end
end

if plt.box == 1
    b = boxplot(allAmp', 'BoxStyle','filled', 'Colors', flip(allCols));
    hold on;
else
    b = bar(x,y);
    hold on;
    errorbar(x, y, error, '.k');
    xticks(x); xticklabels({}); %{'HH', 'HL', 'LH', 'LL'});
    b.FaceColor = 'flat';
    b.CData(:,:) = flip(allCols);

end

% Test for main effects
tAmp = array2table(allAmp');
tAmp.Properties.VariableNames = {'HH','HL','LH','LL'};
WithinStructure = table([1 1 2 2]',[1 2 1 2]','VariableNames',{'P1','P2'});
WithinStructure.P1 = categorical(WithinStructure.P1);
WithinStructure.P2 = categorical(WithinStructure.P2);

rm = fitrm(tAmp, 'HH,HL,LH,LL~1','WithinDesign',WithinStructure);
tAmp_ranovatable = ranova(rm,'WithinModel','P1*P2');


xticks(x); xticklabels({}); %{'HH', 'HL', 'LH', 'LL'});
yl = get(gca(),'Ylim');
yl(2) = yl(2) +20*[yl(2)-yl(1)]/100;
ylim([yl])
xl = get(gca(),'Xlim')+1;
yl = get(gca(),'Ylim');
yl = yl -5*[yl(2)-yl(1)]/100;
ylabel('Max. amp (\muV/m^{2})')


effects = {' p1:', ' p2:'}
for pulse = 1:2
    p = tAmp_ranovatable(3+2*(pulse-1),5).pValue;
      textp = num2str(p);

    if p < 0.001; thisTitle = '***'; fw = 'bold'; 
    elseif p<0.01; thisTitle = '**'; fw = 'bold'; 
    elseif p < 0.05; thisTitle = '*'; fw = 'bold'; 
    else; thisTitle = textp(1:4); fw = 'normal'; end
  
    yl = get(gca(),'Ylim');
    xl = get(gca(),'Xlim');

    t = text(xl(1),yl(2)-10*(pulse-1)*yl(2)/100,[effects{pulse}, thisTitle], ...
        'VerticalAlignment','top', ...
        'HorizontalAlignment', 'left',...
        'FontSize',10, 'FontWeight', fw);

    hold on
end

set(gca, 'FontSize', plt.font)

axes('position',[.25 .75 .2 .2]);
g = 1:length(exp.gapDur{exp.thisExp});
diffTopo_p1 = squeeze(nanmean(nanmean(dataTopo(g,1,:,:,:,:),1),3) - nanmean(nanmean(dataTopo(g,2,:,:,:,:),1),3)); 
diffTopo_p2 = squeeze(nanmean(nanmean(dataTopo(g,:,1,:,:,:),1),2) - nanmean(nanmean(dataTopo(g,:,2,:,:,:),1),2)); 

GT_topoplot(nanmean(diffTopo_p1,2),EEG.chanlocs, ...
    'numcontour', 0, ...
    'plotrad',0.65, ...
    'headrad', 0.5,...
    'emarker2', {chans,'*','k',1,1}); % To plot location of the selected channel


axes('position',[.45 .75 .2 .2]);
GT_topoplot(nanmean(diffTopo_p2,2),EEG.chanlocs, ...
    'numcontour', 0, ...
    'plotrad',0.65, ...
    'headrad', 0.5,...
    'emarker2', {chans,'*','k',1,1}); % To plot location of the selected channel

set(gca, 'FontSize', plt.font)

%% Plot high-low topo for long gap, just before P2
figure; 
g = length(exp.gapDur{exp.thisExp});
timeLabs = {'[-1,-0.9]','[-0.8,-0.7]','[-0.6,-0.5]','[-0.4,-0.3]', '[-0.2,-0.1]'};
ii = 0;
for idx = [5:-1:1]; 
   diffTopo_p1 = squeeze(nanmean(nanmean(topoPreP2(g,1,:,:,:,:,idx),1),3) - nanmean(nanmean(topoPreP2(g,2,:,:,:,:,idx),1),3)); 
ii = ii+1;
subplot(1,5,ii)
title(timeLabs{ii})
GT_topoplot(nanmean(diffTopo_p1,2),EEG.chanlocs, ...
    'numcontour', 0, ...
    'plotrad',0.65, ...
    'headrad', 0.5,...
    'emarker2', {chans,'*','k',1,1}); % To plot location of the selected channel
caxis([-5,5])
end

%% Plot without topos, to match sims width 
%% Plot averaged across gaps 
f = figure; 
f.Units = 'centimeters';
f.OuterPosition = [0.5,0.5,8,14];

c = 5;
slope = []; peak = []; 

meanAmp = []; slope = []; 
gapIdx = 2:length(exp.gapDur{exp.thisExp});
plt.box = 1;
ylims = [-5,15];
lineType = '-';

for p1 = 1:2
    for p2 = 1:2
        c = c-1;
        thisCol = hex2rgb((plt.Rcolors{c}));
        baselines =  squeeze(nanmean(nanmean(data_allC_p2(gapIdx,p1,p2,0.15*sr:0.2*sr,:),4),1));
        if p1 ~= p2
        stdshade(squeeze(nanmean(data_allC_p2(gapIdx,p1,p2,:,:),1))'-baselines, plt.alpha, thisCol, lineType,[],plt.smooth) %Plot ERP traces
        hold on;
        end
             
    end
end
ylim([ylims])
set(gca, 'fontsize', plt.font)
xticklabels({'P2', '0.5', '1'})
xticks([[0.2,0.7,1.2]*512])
xline(0.2*512)
ylabel(plt.ERP.ylabel);
f.Units = 'centimeters';
xlabel('Time (s)'); 
f.Units = 'centimeters';
f.OuterPosition = [0.5,0.5,6.5,9];

