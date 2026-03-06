%% Plot simulated accuracies
function [] = GT_S_plotAccuracies(allSims,allPars,simPars,exp,plt)

%% New version, for paper
%Load real data
f = figure;
for e = 1:2 %For both experiments
    
    exp.thisExp = e; %Experiment to analyse
    [exp,  pars]  = GT_S_setup(exp); %Get correct paths/values

    load(['C:\Users\elisa\Desktop\GapsTask\Modelling & simulations\GT_Sonic\GT_gaDat_' exp.thisExpCode '.mat']);
    sims = allSims{e};
    pars = allPars{e};
    
    pars.accType = 'bin';
    
    acc = []; p1acc = []; accLab = {}; noBacc = [];
    for c =length(pars.cLab):-1:3
        thisCond = sims.simdat(:,1) == c;
        if strcmp(pars.accType, 'bin'), thisAcc = 2; elseif strcmp(pars.accType, 'erf'), thisAcc = 3; end
        sims.simdat(sims.simdat(:,2) == -1,2) = 0;
        acc(c-2) = mean(sims.simdat(thisCond,thisAcc));
        accLab{c-2} = pars.cLab{c};
    end
    
    acc = reshape(acc, [4,length(exp.gapDur{e})]);
    accLab = reshape(accLab, [4, length(exp.gapDur{e})]);
    
    
    for c = 1:2%length(pars.cLab)
        thisCond = sims.simdat(:,1) == c;
        if strcmp(pars.accType, 'bin'), thisAcc = 2; elseif strcmp(pars.accType, 'erf'), thisAcc = 3; end
        sims.simdat(sims.simdat(:,2) == -1,2) = 0;
        p1acc(c) = mean(sims.simdat(thisCond,thisAcc));
        p1accLab{c} = pars.cLab{c};
    end
    
    acc = [repmat(p1acc(2),1,length(exp.gapDur{e}));repmat(p1acc(1),1,length(exp.gapDur{e})); flip(acc)];
    
    subplot(2,3,[1,2]+3*(e-1))
    r1 = plot([-0.3:0.01:-0.2],repmat(acc(1,1)',1,11), 'LineStyle', ':','LineWidth', 2); hold on;
    rr1 = plot([-0.3:0.01:-0.2],repmat(squeeze(ga.p(1,1,1)),1,11),  'LineStyle','-','LineWidth', 2); hold on; %', 'none', 'Marker', '*', 'MarkerSize', 5); hold on;
    
    r2 = plot([-0.2:0.01:-0.1],repmat(acc(2,1)',1,11), 'LineStyle', ':','LineWidth', 2); hold on;
    rr2 = plot([-0.2:0.01:-0.1],repmat(squeeze(ga.p(1,1,2)),1,11),  'LineStyle', '-','LineWidth', 2); hold on; %', 'none', 'Marker', '*', 'MarkerSize', 5); hold on;
    
    hold on;
    
    l = plot(exp.gapDur{e},acc(3:end,:)', 'LineStyle', ':','LineWidth', 2); hold on;
    ll = plot(exp.gapDur{e},squeeze(ga.p(1,:,3:end)), 'LineWidth', 2); hold on; %', 'none', 'Marker', '*', 'MarkerSize', 5); hold on;
    
    flippedCols = (plt.Rcolors);
    
    for k = 1:4
        l(k).Color = hex2rgb(flippedCols{k});
        ll(k).Color = hex2rgb(flippedCols{k});
        hold on;
    end
    r1.Color = hex2rgb(flippedCols{1});
    rr1.Color= hex2rgb(flippedCols{1});
    
    r2.Color = hex2rgb(flippedCols{4});
    rr2.Color= hex2rgb(flippedCols{4});
    
    xticks([-0.25,-0.15,exp.gapDur{e}]);
    gapLabels = []; 
    for g = 1:length(exp.gapDur{e}); 
        gapLabels = [gapLabels, {num2str(exp.gapDur{e}(g))}];
    end
    xticklabels({'L-', 'H-', gapLabels{:}})
    set(gca, 'fontsize', 11);
    ylim([0.6,0.9]);
    xlim([-0.4,1.08]);

    if e == 2; xlabel('Gap duration (s)', 'fontsize', 11); else xlabel(''); end
    %     ylabel('P(correct)', 'fontsize', 11);
    
    xline(-0.05, '-k')
    box off;
        if e == 2; ylabel('                                               P(correct)'); end
    
    % Early terminations
    subplot(2,3,[3]+3*(e-1))
    % Panel C, %early terminations
   
    flippedCols = flip(plt.Rcolors);
    for conds = 3:6%length(pars.cLab)/4
        thisLab = pars.cLab{conds}(6:end);
        indices{conds} = contains(pars.cLab, thisLab);
        condList{conds-2} = thisLab;
        disp(['Ncells:' num2str(sum(indices{conds}))]);
    end
    
    for c = 1:length(pars.cLab)
        thisP1prop_z1(c) = nanmean(sims.bndTP1(c,:));
    end
    
    for conds = 3:6
        thisP1prop_z1xgap(conds) = nanmean(thisP1prop_z1(indices{conds}));
    end
   
    
    flippedCols = flip(plt.Rcolors);
    i = 0;
    for c = [6:-1:3]
        i = i+1;
        l = bar(i,mean(thisP1prop_z1xgap(c))); hold on;
        l.FaceColor = hex2rgb(flippedCols{c-2});
        yline(0.5, '--k')
%         ylabel({exp.name{e}});
        if e == 2; ylabel('                                               P(Early termination)'); end
        if e == 2;
            xticks([1,2,3,4]); xticklabels({'LL','LH','HL', 'HH'});
        
        end
        ylim([0.2,0.9])
        set(gca, 'fontsize', 11);

    box off; 

    end
    
end

f.Units = 'centimeters';
f.OuterPosition = [2 2 18 14]














% % {[plt.modelTitle]},...
% % sgtitle([{['EXP 1: b:' num2str(simPars.b(1)) ', d:' num2str(simPars.d(1))]}]),...', sz:' num2str(simPars.sz(1))]},...
% %     {['EXP 2: b:' num2str(simPars.b(2)) ', d:' num2str(simPars.d(2)) ', sz:' num2str(simPars.sz(2))]}]);


% %% New version, for paper
% close all
% f = figure;
%
% for e = 1:2 %For both experiments
%
%     sims = allSims{e};
%     pars = allPars{e};
%       exp.thisExp = e; %Experiment to analyse
%             [exp,  pars]  = GT_S_setup(exp); %Get correct paths/values
%     %Load real data
%     load(['C:\Users\elisa\Desktop\GapsTask\Modelling & simulations\GT_Sonic\GT_gaDat_' exp.thisExpCode '.mat']);
%     pars.accType = 'bin';
%
%     acc = []; p1acc = []; accLab = {}; noBacc = [];
%     for c =length(pars.cLab):-1:3
%         thisCond = sims.simdat(:,1) == c;
%         if strcmp(pars.accType, 'bin'), thisAcc = 2; elseif strcmp(pars.accType, 'erf'), thisAcc = 3; end
%         sims.simdat(sims.simdat(:,2) == -1,2) = 0;
%         acc(c-2) = mean(sims.simdat(thisCond,thisAcc));
%         accLab{c-2} = pars.cLab{c};
%     end
%
%     acc = reshape(acc, [4,length(exp.gapDur{e})]);
%     accLab = reshape(accLab, [4, length(exp.gapDur{e})]);
%
%
%     for c = 1:2%length(pars.cLab)
%         thisCond = sims.simdat(:,1) == c;
%         if strcmp(pars.accType, 'bin'), thisAcc = 2; elseif strcmp(pars.accType, 'erf'), thisAcc = 3; end
%         sims.simdat(sims.simdat(:,2) == -1,2) = 0;
%         p1acc(c) = mean(sims.simdat(thisCond,thisAcc));
%         p1accLab{c} = pars.cLab{c};
%     end
%
%     acc = [repmat(p1acc(2),1,length(exp.gapDur{e}));repmat(p1acc(1),1,length(exp.gapDur{e})); flip(acc)];
%
%     %% Panel B, 2-Pulse pooled accuracies
%     subplot(2,2,1+(e-1)*2)
% % figure
%     flippedCols = (plt.Rcolors);
%     flippedCols = {flippedCols{1}, flippedCols{1:2}, flippedCols{4}, flippedCols{3:4}};
%
%     cIdx = [1,4,2,3,5,6];
%     i = 0;
%     for c = 1:6
%         i = i+1;
%         realAcc(cIdx(i)) = nanmean(ga.p(1,:,c),2);
%     end
%
%     i = 0;
%     for c = [1,4]%1:6
%         i = i+1;
%         l = bar(i,mean(realAcc(c))); hold on;
%         if c == 1 || c == 4
%             hatchfill2(l(1),'single','HatchAngle',45,'hatchcolor','w', 'HatchLineWidth', 2);
%         end
%
%         l.FaceColor = [1 1 1];
%         l.EdgeColor = hex2rgb(flippedCols{c}); l.LineWidth = 2;
%         if e == 2; xticks([1,2,3,4,5,6]); xticklabels({'L-','H-'});
%         else xticklabels (''); end
%
%     end
%
%     i = 0;
%     for c = [1,4]
%         i = i+1;
%         plot(i,mean(acc(c,:)), '*', 'color','k', 'markersize', 5); hold on;  %Simulated data
%     end
%
%
%     ylim([0.6,0.9]);
%     set(gca, 'fontsize', 11);
%     if e == 1, title('2 pulse'); end
%     box off;
%         ylabel('P(correct)');
%
%
%     %% Panel D, Accuracies x gap line graph
%     subplot(2,2,[2+(e-1)*2, 2+(e-1)*2])
%     l = plot(exp.gapDur{e},acc(3:end,:)', 'LineStyle', '--'); hold on;
%     ll = plot(exp.gapDur{e},squeeze(ga.p(1,:,3:end)), 'LineWidth', 2); hold on; %', 'none', 'Marker', '*', 'MarkerSize', 5); hold on;
%
%     flippedCols = (plt.Rcolors);
%
%     for k = 1:4
%         l(k).Color = hex2rgb(flippedCols{k});
%         ll(k).Color = hex2rgb(flippedCols{k});
%     end
%     xticks(exp.gapDur{e});
%     set(gca, 'fontsize', 11);
%     ylim([0.6,0.9]);
%     xlim([0,1.08]);
%     %     yticks([]);
%     if e == 2; xlabel('Gap duration (s)', 'fontsize', 11); else xlabel(''); end
%     %     ylabel('P(correct)', 'fontsize', 11);
%     if e == 1, title('2 pulse x gap'); end
%
% %     for c = 1:2
% %         thisCol = hex2rgb(flippedCols{c*c});
% %         yline(p1acc(c), 'linestyle', ':', 'color', thisCol, 'linewidth', 2); %Simulated
% %         yline(squeeze(ga.p(1,:,c)), 'linestyle', '--', 'color', plt.Rcolors{c*c}, 'linewidth', 1); %Real
% %
% %         hold on;
% %     end
%
%     box off;
% end
%
%
%
% f.Units = 'centimeters';
% f.OuterPosition = [2 2 15 14]
% % {[plt.modelTitle]},...
% % sgtitle([{['EXP 1: b:' num2str(simPars.b(1)) ', d:' num2str(simPars.d(1))]}]),...', sz:' num2str(simPars.sz(1))]},...
% %     {['EXP 2: b:' num2str(simPars.b(2)) ', d:' num2str(simPars.d(2)) ', sz:' num2str(simPars.sz(2))]}]);
%
% sgtitle([{['EXP 1: b:' num2str(simPars.b(1)) ', d:' num2str(simPars.d(1)) ', z:' num2str(simPars.z(1))]},...
%     {['EXP 2: b:' num2str(simPars.b(2)) ', d:' num2str(simPars.d(2)) ', z:' num2str(simPars.z(2))]}]);
%
%
%% Old versions
% %%
% f = figure
%
%
% for e = 1:2 %For both experiments
%
%     sims = allSims{e};
%     pars = allPars{e};
%       exp.thisExp = e; %Experiment to analyse
%             [exp,  pars]  = GT_S_setup(exp); %Get correct paths/values
%     %Load real data
%     load(['C:\Users\elisa\Desktop\GapsTask\Modelling & simulations\GT_Sonic\GT_gaDat_' exp.thisExpCode '.mat']);
%     pars.accType = 'bin';
%
%     acc = []; p1acc = []; accLab = {}; noBacc = [];
%     for c = 3:length(pars.cLab)
%         thisCond = sims.simdat(:,1) == c;
%         if strcmp(pars.accType, 'bin'), thisAcc = 2; elseif strcmp(pars.accType, 'erf'), thisAcc = 3; end
%         sims.simdat(sims.simdat(:,2) == -1,2) = 0;
%         acc(c-2) = mean(sims.simdat(thisCond,thisAcc));
%         accLab{c-2} = pars.cLab{c};
%     end
%     for c = 1:2
%         p1acc(c) = mean(sims.simdat(sims.simdat(:,1) == c,2)); %1p coherence
%     end
%
%     acc = reshape(acc, [4,length(exp.gapDur{e})]);
%     accLab = reshape(accLab, [4, length(exp.gapDur{e})]);
%
%     %% Panel B, 2-Pulse pooled accuracies
%     subplot(2,2,1+(e-1)*2)
%     flippedCols = flip(plt.Rcolors);
%
%     i = 0;
%     for c = [6:-1:3]
%         i = i+1;
%         realAcc(i) = nanmean(ga.p(1,:,c),2);
%     end
%
%     i = 0;
%     for c = [4:-1:1]
%         i = i+1;
%         l = bar(i,mean(realAcc(c))); hold on;
%         l.FaceColor = hex2rgb(flippedCols{c});
%         if e == 2; xticks([1,2,3,4]); xticklabels({'LL','LH', 'HL', 'HH'});
%         else xticklabels (''); end
%         plot(i,mean(acc(c,:)), '*', 'color','k'); hold on;  %Simulated data
%     end
%     ylim([0.5,0.9]);
%     %     ylabel('P(correct)');
%     set(gca, 'fontsize', 11);
%     if e == 1, title('2 pulse'); end
%     box off;
%
%
%     %% Panel D, Accuracies x gap line graph
%     subplot(2,2,[2+(e-1)*2, 2+(e-1)*2])
%     l = plot(exp.gapDur{e},acc', 'LineStyle', 'none', 'Marker', '*', 'MarkerSize', 5); hold on;
%     ll = plot(exp.gapDur{e},squeeze(ga.p(1,:,c+2:end)), 'LineWidth', 2)%', 'none', 'Marker', '*', 'MarkerSize', 5); hold on;
%
%     %     la = plot(exp.gapDur{e},noBacc', 'LineWidth', 2, 'LineStyle', ':', 'Marker', 'diamond', 'MarkerSize', 5);
%
%     for k = 1:size(acc,1)
%         l(k).Color = hex2rgb(flippedCols{k});
%         ll(k).Color = hex2rgb(plt.Rcolors{k});
%     end
%     xticks(exp.gapDur{e});
%     set(gca, 'fontsize', 11);
%     ylim([0.6,0.9]);
%     xlim([0,1.08]);
%     %     yticks([]);
%     if e == 2; xlabel('Gap duration (s)', 'fontsize', 11); else xlabel(''); end
%     %     ylabel('P(correct)', 'fontsize', 11);
%     if e == 1, title('2 pulse x gap'); end
%
%     for c = 1:2
%         thisCol = hex2rgb(flippedCols{c*c});
%         yline(p1acc(c), 'linestyle', ':', 'color', thisCol, 'linewidth', 2); %Simulated
%         yline(squeeze(ga.p(1,:,c)), 'linestyle', '--', 'color', plt.Rcolors{c*c}, 'linewidth', 1); %Real
%
%         hold on;
%     end
%
%     box off;
%
%
% end
%
% f.Units = 'centimeters';
% f.OuterPosition = [2 2 11 12]
% % {[plt.modelTitle]},...
% sgtitle([{['EXP 1: b:' num2str(simPars.b(1)) ', d:' num2str(simPars.d(1)) ', sz:' num2str(simPars.sz(1))]},...
%     {['EXP 2: b:' num2str(simPars.b(2)) ', d:' num2str(simPars.d(2)) ', sz:' num2str(simPars.sz(2))]}]);
%
%
% %% As main figure
% close all
% f = figure;
%
% for e = 1:2 %For both experiments
%
%     sims = allSims{e};
%     pars = allPars{e};
%       exp.thisExp = e; %Experiment to analyse
%             [exp,  pars]  = GT_S_setup(exp); %Get correct paths/values
%     %Load real data
%     load(['C:\Users\elisa\Desktop\GapsTask\Modelling & simulations\GT_Sonic\GT_gaDat_' exp.thisExpCode '.mat']);
%     pars.accType = 'bin';
%
%     acc = []; p1acc = []; accLab = {}; noBacc = [];
%     for c =length(pars.cLab):-1:3
%         thisCond = sims.simdat(:,1) == c;
%         if strcmp(pars.accType, 'bin'), thisAcc = 2; elseif strcmp(pars.accType, 'erf'), thisAcc = 3; end
%         sims.simdat(sims.simdat(:,2) == -1,2) = 0;
%         acc(c-2) = mean(sims.simdat(thisCond,thisAcc));
%         accLab{c-2} = pars.cLab{c};
%     end
%
%     acc = reshape(acc, [4,length(exp.gapDur{e})]);
%     accLab = reshape(accLab, [4, length(exp.gapDur{e})]);
%
%
%     for c = 1:2%length(pars.cLab)
%         thisCond = sims.simdat(:,1) == c;
%         if strcmp(pars.accType, 'bin'), thisAcc = 2; elseif strcmp(pars.accType, 'erf'), thisAcc = 3; end
%         sims.simdat(sims.simdat(:,2) == -1,2) = 0;
%         p1acc(c) = mean(sims.simdat(thisCond,thisAcc));
%         p1accLab{c} = pars.cLab{c};
%     end
%
%     acc = [repmat(p1acc(2),1,length(exp.gapDur{e}));repmat(p1acc(1),1,length(exp.gapDur{e})); flip(acc)];
%
%     %% Panel B, 2-Pulse pooled accuracies
%     subplot(2,2,1+(e-1)*2)
% % figure
%     flippedCols = (plt.Rcolors);
%     flippedCols = {flippedCols{1}, flippedCols{1:2}, flippedCols{4}, flippedCols{3:4}};
%
%     cIdx = [1,4,2,3,5,6];
%     i = 0;
%     for c = 1:6
%         i = i+1;
%         realAcc(cIdx(i)) = nanmean(ga.p(1,:,c),2);
%     end
%
%     i = 0;
%     for c = 1:6
%         i = i+1;
%         l = bar(i,mean(realAcc(c))); hold on;
%         if c == 1 || c == 4
%             hatchfill2(l(1),'single','HatchAngle',45,'hatchcolor','w', 'HatchLineWidth', 2);
%         end
%
%         l.FaceColor = [1 1 1];
%         l.EdgeColor = hex2rgb(flippedCols{c}); l.LineWidth = 2;
%         if e == 2; xticks([1,2,3,4,5,6]); xticklabels({'L-','LL','LH','H-', 'HL', 'HH'});
%         else xticklabels (''); end
%
%     end
%
%     i = 0;
%     for c = 1:6
%         i = i+1;
%         plot(cIdx(i),mean(acc(c,:)), '*', 'color','k', 'markersize', 5); hold on;  %Simulated data
%     end
%
%
%     ylim([0.6,0.9]);
%     set(gca, 'fontsize', 11);
%     if e == 1, title('2 pulse'); end
%     box off;
%         ylabel('P(correct)');
%
%
%     %% Panel D, Accuracies x gap line graph
%     subplot(2,2,[2+(e-1)*2, 2+(e-1)*2])
%     l = plot(exp.gapDur{e},acc(3:end,:)', 'LineStyle', 'none', 'Marker', '*', 'MarkerSize', 5); hold on;
%     ll = plot(exp.gapDur{e},squeeze(ga.p(1,:,3:end)), 'LineWidth', 2); hold on; %', 'none', 'Marker', '*', 'MarkerSize', 5); hold on;
%
%     flippedCols = (plt.Rcolors);
%
%     for k = 1:4
%         l(k).Color = hex2rgb(flippedCols{k});
%         ll(k).Color = hex2rgb(flippedCols{k});
%     end
%     xticks(exp.gapDur{e});
%     set(gca, 'fontsize', 11);
%     ylim([0.6,0.9]);
%     xlim([0,1.08]);
%     %     yticks([]);
%     if e == 2; xlabel('Gap duration (s)', 'fontsize', 11); else xlabel(''); end
%     %     ylabel('P(correct)', 'fontsize', 11);
%     if e == 1, title('2 pulse x gap'); end
%
% %     for c = 1:2
% %         thisCol = hex2rgb(flippedCols{c*c});
% %         yline(p1acc(c), 'linestyle', ':', 'color', thisCol, 'linewidth', 2); %Simulated
% %         yline(squeeze(ga.p(1,:,c)), 'linestyle', '--', 'color', plt.Rcolors{c*c}, 'linewidth', 1); %Real
% %
% %         hold on;
% %     end
%
%     box off;
%
%
% end
%
% f.Units = 'centimeters';
% f.OuterPosition = [2 2 15 14]
% % {[plt.modelTitle]},...
% sgtitle([{['EXP 1: b:' num2str(simPars.b(1)) ', d:' num2str(simPars.d(1)) ', sz:' num2str(simPars.sz(1))]},...
%     {['EXP 2: b:' num2str(simPars.b(2)) ', d:' num2str(simPars.d(2)) ', sz:' num2str(simPars.sz(2))]}]);
%
% sgtitle([{['EXP 1: b:' num2str(simPars.b(1)) ', d:' num2str(simPars.d(1)) ', z:' num2str(simPars.z(1))]},...
%     {['EXP 2: b:' num2str(simPars.b(2)) ', d:' num2str(simPars.d(2)) ', z:' num2str(simPars.z(2))]}]);
%
%
