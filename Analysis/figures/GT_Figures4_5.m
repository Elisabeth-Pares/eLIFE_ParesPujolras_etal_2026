function [] = GT_Figures4_5(allSims, allPars,  simPars,exp, plt)
% This code simulates a diffusion process using the best-fitting parameters
% for each experiment. The default parameters can be used to reproduce:
% Fig. 4 - behavioural simulations
% Fig. 5 - CPP & MBL simulations  in the main manuscript. 

% Changing specific parameters of the code as indicated below can be used
% to reproduce Figures S6, S7 and S8 in the supplementary material. 

% Fig. S6 - simulate CPP/MBL traces in correct/error trials separately 

% Fig. S7 - explores parameters capture assumptions about the time for the signal (CPP) to
% decay given no bound crossing (pars.stepAfterOffset), and the duration of post-decisional 
% evidence accumulation (pars.delayDecay)

% Fig. S8 - illustrates what the CPP would look like under the assumption
% of no signal decay during the inter-pulse gap.(pars.decayBlue)


%% INITIALISE EXP PARAMETERS
clc
clear

exp.thisExp = 1; %Experiment to analyse

[exp, plt, pars]  = GT_setup(exp);

exp.sub_id{1}            = [1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]; % P8 excluded because only 4 blocks were recorded.
exp.sub_id{2}            = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,23]; % P1 & P21 excluded because they could not perform task well

%% Plot SONIC fitted models
close all

modelNames = {'FixedBound_1DV_2DR', 'FixedBound_1DV_2DR_DRvar'}; %Model names 
Selvecs = {...
    [1 0 1 1 0 0 0 0 0 0 0 0],...;%FixedBound_1DV_2DR - Model in the main manuscript, with Bound, DR1 and DR2
    [1 0 1 1 0 0 0 0 0 1 0 0],...;   %FixedBound_2DR_DRvar - Supplementary model with additional trial-by-trial variability in DR2 
    };

N = 1000; %Number of simulations
i =0;
for m = 1
            i = i+1;
            for e = 1:2 %For the two experiments
                modelNames = {'FixedBound_1DV_2DR', 'FixedBound_1DV_2DR_DRvar'};
                modelName = modelNames{m};

                exp.thisExp = e; %Experiment to analyse
                [exp,  pars]  = GT_S_setup(exp); %Get correct paths/values
                pars.nSim = N;
                
                pars.stepAfterbndT = 20; % Time for signal to go back to zero after hitting bound: 200ms (Fig. 5)
                pars.stepAfterOffset = 20; %Time for signal to go back to zero after dot offset if no bound is hit: 200 ms (Fig. 5). For Fig. S7, set to [10,50,100]
                pars.decayDelay = 0; %Time until signal starts falling off - post-decisional accumulation after bound is hit: assume no post-bound acc. For Fig. S8, set to [5,10]
                
                %Load real data
                load(['C:\Users\elisa\Desktop\GapsTask\Modelling & simulations\GT_Sonic\GT_gaDat_' exp.thisExpCode '.mat']);
                
                %Load model results
                load(['C:\Users\elisa\Desktop\GapsTask\Modelling & simulations\GT_Sonic\GT_S_ModelFits\FixedBound_1DV_2DR_fitType0_n10000.mat']);

                pm = squeeze(bestpm(e,3,:));
                Sel = pm ~= 0; pvar = pm(pm ~= 0); pfix = pm(pm == 0);
                
                pars = GT_S_getSimPars(exp, pars, modelName);
                pars.simpleFit = 0; %Always simulate all conditions
                [allSims{e}, allPars{e}] = GT_S_simul_DDM_clean(pars, pvar, pfix, Sel);
                
            end
            plt.modelTitle = ['Model ' num2str(i), 'fitType = 0'];
            pars.accType = 'bin';
end

% save('C:\Users\elisa\Desktop\GapsTask\Modelling & simulations\GT_Sonic\GT_allSims_n50000.mat', 'allSims', 'allPars', '-v7.3')

% ===========================================================================
%% PLOT SIMULATIONS, with real data overlaid
% ===========================================================================
%% Figure 4
%Load real data
f = figure;
for e = 1:2 %For both experiments
    
    exp.thisExp = e; %Experiment to analyse
    [exp,  pars]  = GT_S_setup(exp); %Get correct paths/values
    
    load(['C:\Users\elisa\Desktop\GapsTask\Modelling & simulations\GT_Sonic\GT_gaDat_' exp.thisExpCode '.mat']);
    sims = allSims{e};
    sims.simdat(sims.simdat(:,2) == -1,2) = 0;
    
    pars = allPars{e};
    pars.accType = 'bin';
    
    acc = []; p1acc = []; accLab = {}; noBacc = [];
    
    for c =length(pars.cLab):-1:3
        thisCond = sims.simdat(:,1) == c;
        
        if strcmp(pars.accType, 'bin'), thisAcc = 2; elseif strcmp(pars.accType, 'erf'), thisAcc = 3; end
        acc(c-2) = mean(sims.simdat(thisCond,thisAcc));
        accLab{c-2} = pars.cLab{c};
    end
    
    acc = reshape(acc, [4,length(exp.gapDur{e})]);
    accLab = reshape(accLab, [4, length(exp.gapDur{e})]);
    box on
    
    for c = 1:2%length(pars.cLab)
        thisCond = sims.simdat(:,1) == c;
        if strcmp(pars.accType, 'bin'), thisAcc = 2; elseif strcmp(pars.accType, 'erf'), thisAcc = 3; end
        sims.simdat(sims.simdat(:,2) == -1,2) = 0;
        p1acc(c) = mean(sims.simdat(thisCond,thisAcc));
        p1accLab{c} = pars.cLab{c};
    end
    
    acc = [repmat(p1acc(2),1,length(exp.gapDur{e}));repmat(p1acc(1),1,length(exp.gapDur{e})); flip(acc)];
    
    subplot(2,4,[1,2]+4*(e-1))
    r1 = plot([-0.3:0.01:-0.2],repmat(acc(1,1)',1,11), 'LineStyle', '-','LineWidth', 2); hold on;
    rr1 = plot([-0.3:0.01:-0.2],repmat(squeeze(ga.p(1,1,1)),1,11),  'LineStyle',':','LineWidth', 2); hold on; %', 'none', 'Marker', '*', 'MarkerSize', 5); hold on;
    
    r2 = plot([-0.2:0.01:-0.1],repmat(acc(2,1)',1,11), 'LineStyle', '-','LineWidth', 2); hold on;
    rr2 = plot([-0.2:0.01:-0.1],repmat(squeeze(ga.p(1,1,2)),1,11),  'LineStyle', ':','LineWidth', 2); hold on; %', 'none', 'Marker', '*', 'MarkerSize', 5); hold on;
    
    hold on;
    
    l = plot(exp.gapDur{e},acc(3:end,:)', 'LineStyle', '-','LineWidth', 2); hold on;
    ll = plot(exp.gapDur{e},squeeze(ga.p(1,:,3:end)), 'LineWidth', 2, 'LineStyle', ':'); hold on; %', 'none', 'Marker', '*', 'MarkerSize', 5); hold on;
    
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
    
    xline(-0.05, '-k')
    box on;
    if e == 2; ylabel('                                               P(correct)'); end
    
    
    
    %% Early terminations
    subplot(2,4,[3]+4*(e-1))
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
        if e == 2; ylabel('                                               P(Early termination)'); end
        xticks([1,2,3,4]); xticklabels({'LL','LH','HL', 'HH'});
        
        ylim([0.2,0.9])
        set(gca, 'fontsize', 11);
        
        box on;
        
    end
    
    %% Plot DV at P1
    subplot(2,4,[4]+4*(e-1))
    
    for c = 1:length(pars.cLab)
        idx = sims.bndTP1(c,:) == 0;
        thisP1_DV(c) = nanmean(abs(sims.p1DV_end(c,idx)),2);
    end
    
    for conds = 3:6
        thisP1_DV_xgap(conds) = nanmean(thisP1_DV(indices{conds}));
    end
    
    
    flippedCols = flip(plt.Rcolors);
    i = 0;
    for c = [6:-1:3]
        i = i+1;
        l = bar(i,mean(thisP1_DV_xgap(c))); hold on;
        l.FaceColor = hex2rgb(flippedCols{c-2});
        if e == 2; ylabel('                                         |DV| at P1 if B not hit'); end
        xticks([1,2,3,4]); xticklabels({'LL','LH','HL', 'HH'});
        
                if e == 2; ylim([0.013,0.016]); elseif e == 1; ylim([0.015,0.02]); end
        set(gca, 'fontsize', 11);
        
        box on;
        
    end
    
end

f.Units = 'centimeters';
f.OuterPosition = [2 2 21 14]
exportgraphics(f, [exp.figPath, 'GT_Fig4_noLine.tiff'], 'Resolution', 600)
exportgraphics(f, [exp.figPath, 'GT_Fig4_noLine.pdf'], 'Resolution', 600)

% ====================================================================
%% PLOT EEG SIMULATIONS
% ====================================================================
% Plot MBL & CPP for each experiment separately 
ylims = ([0 0.035]);
sortBy.acc = 1; 
decayBlue = 1; %Does the CPP go down when dots turn blue? 1 - yes, 0 - no

for e = 1:2 %For both experiments
    f = figure
    sims = allSims{e};
    pars = allPars{e};
    
    sims.wfDV(:,pars.nSim+1:end,:) = [];
    sims.decayDV(:,pars.nSim+1:end,:) = [];

    pars.sr = 100; %Simulations sampling rate
    ds = pars.sr/512; %To scale plotting limits
       
    acc_byCond=[];
    acc_byCond = reshape(sims.simdat(:,2), [size(sims.decayDV,2), size(sims.decayDV,1)])';
    
    for gap = 1:length(exp.gapLabs{e});
        cols = (plt.satPalettes);
        thisCol = cols(2,gap);

        cidx = contains(pars.cLab, [exp.gapLabs{e}{gap}]);
        if e == 1 & gap == 1,         cidx = contains(pars.cLab, [exp.gapLabs{e}{gap}]) & ~contains(pars.cLab, [exp.gapLabs{e}{gap+1}]) ; end
        lt = '-';
        
        if decayBlue == 1; thisCPP = squeeze(nanmean(abs(sims.decayDV(cidx,:,:)),2));
        else, thisCPP = squeeze(nanmean(abs(sims.wfDV(cidx,:,:)),2)); end

        subplot(2,1,2) %PLOT CPP 
        ylims = [0 0.035]
        baseline = 0;
        plot(squeeze(nanmean(thisCPP(:,55:end))-baseline),'LineWidth', 2, 'Color',thisCol{:}, 'LineStyle', lt);
        xlim([1,2.2*pars.sr])
        hold on;
        fill(plt.gapShades{e}(gap+1,:)*ds,[ylims, flip(ylims)], thisCol{:}, 'LineStyle','none','facealpha', 0.2); %Plot pulses in time
        hold on;yline(0,'--k'); hold on;
        if gap == 1;  fill(plt.gapShades{e}(1,:)*ds,[ylims, flip(ylims)], 'k', 'LineStyle','none','facealpha', 0.2);end %Plot pulses in time
        ylabel('Simulated CPP');
        xticks([0.2*pars.sr, pars.sr+0.2*pars.sr, pars.sr*2.2]); xticklabels({'0', '1',  '2'});xline(0.2*pars.sr); set(gca, 'fontsize', 12);
        
        %Compute %CPP-P1/CPP-P2 across gaps 
        p1Peak = [0.75*100,0.95*100]; p2Peak =[round((0.75+0.2+exp.gapDur{e}(gap))*100): round((0.75+exp.gapDur{e}(gap)+0.2+0.2)*100)];
        cpp1 = nanmean(nanmean(thisCPP(:,p1Peak),2)); cpp2 = nanmean(nanmean(thisCPP(:,p2Peak),2))-nanmean(nanmean(thisCPP(:,p2Peak(1)),2));
        percCPP(e,gap,:) = cpp2/cpp1*100; clear cpp1; clear cpp2; 
        
        subplot(2,1,1) %PLOT MBL (signed DV) -- consider plotting correct only, so bound can be illustrated.
        baseline = 0;        ylims = [0 0.05]
        thisMBL = squeeze(nanmean((sims.MBL(cidx,:,:)),2));
        plot(squeeze(nanmean(thisMBL(:,55:end))-baseline),'LineWidth', 2, 'Color',thisCol{:}, 'LineStyle', lt);
        xlim([1,2.2*pars.sr])
        hold on;
        fill(plt.gapShades{e}(gap+1,:)*ds,[ylims, flip(ylims)], thisCol{:}, 'LineStyle','none','facealpha', 0.2); %Plot pulses in time
        hold on;yline(0,'--k'); hold on;
        if gap == 1;  fill(plt.gapShades{e}(1,:)*ds,[ylims, flip(ylims)], 'k', 'LineStyle','none','facealpha', 0.2);end %Plot pulses in time
        ylabel('Simulated MBL');    xticks([0.2*pars.sr, pars.sr+0.2*pars.sr, pars.sr*2.2]); xticklabels({'0', '1',  '2'});xline(0.2*pars.sr); set(gca, 'fontsize', 12);
    end

    f.Units = 'centimeters';
    f.OuterPosition = [0.25 0.35 11 15]; 
    
end
disp(['Exp. 1 ratio:' num2str(nanmean(percCPP(1,2:3)))]) %Average over non-zero gaps
disp(['Exp. 2 ratio:' num2str(nanmean(percCPP(2,2:4)))]) %Average over non-zero gaps


%% Sort & plot DV for each pulse separately
f = figure
sortBy.acc = 0;
decayBlue = 1;

if sortBy.acc == 1, ceIdxx = [1:2]; else ceIdxx = 1; end

if decayBlue == 1; varp1 = 'decayP1'; varp2 = 'decayP2';
else varp1 = 'simp1', varp2 = 'simp2'; end
clear all_p1mbl; clear all_p1sim; clear all_p2mbl; clear all_p2sim; 
lineTypes = {'-', ':'};

for e = 1:2 %For both experiments
    
    sims = allSims{e};
    pars = allPars{e};
    
    t0 = 0.75*pars.sr;
    
    acc_byCond=[];
    acc_byCond = reshape(sims.simdat(:,2), [size(sims.p1DV,2), size(sims.p1DV,1)])';
    
    for ce = ceIdxx(1:end)
        sortBy.acc = 1;
        accVals = [1,-1];
        thisAcc = accVals(ce);
        c = 2;
        
        for g = 1:length(exp.gapDur{e})
            subplot(2,length(exp.gapDur{e}),g+(length(exp.gapDur{e})*(e-1)));
            thisGap = exp.gapLabs{e}{g};
                    cc = 5;

            for p1 = 1:2
                thisP1 = ['P1' exp.cohLabs{p1}];
                for p2 = 1:2
                    c = c+1;
                    cc = cc-1;
                    thisP2 = ['P2', exp.cohLabs{p2}];
                    
                    if sortBy.acc == 1, ceIdx = acc_byCond(c,:) == thisAcc; else ceIdx = logical(ones(1, size(sims.p1DV,2))); end
                    
                    cTag = [thisGap '_' thisP1 '_' thisP2];
                    allLabs{cc} = cTag;
                    
                    p1sim = squeeze(nanmean(abs(sims.(varp1)(c,ceIdx,:)-sims.(varp1)(c,ceIdx,1)),2));
                    p2sim = squeeze(nanmean(abs(sims.(varp2)(c,ceIdx,:))-(sims.(varp2)(c,ceIdx,1)),2));
                        
                    if decayBlue == 1
                        cppsim = squeeze(nanmean(abs(sims.decayDV(c,ceIdx,:)),2));        
                    else
                        cppsim = squeeze(nanmean(abs(sims.wfDV(c,ceIdx,:)),2));
                    end
                    
                    p1mbl = squeeze(nanmean(sims.MBL(c,ceIdx,:),2));%-sims.MBL(c,ceIdx,1)),2));
                    p2mbl= squeeze(nanmean(sims.MBL(c,ceIdx,:),2));%-(sims.MBL(c,ceIdx,1)),2));
                    hold on;
                    plot(1:length(p2sim(1:2.75*pars.sr)),p2sim(1:2.75*pars.sr), 'LineWidth', 2, 'Color', plt.Rcolors{cc}, 'LineStyle', ':');hold on;
                    
                    all_p1sim(e,ce,cc,g,:) = (p1sim(t0-0.2*pars.sr:t0-0.2*pars.sr+2*pars.sr));
                    all_p2sim(e,ce,cc,g,:) = (p2sim(t0+(exp.gapDur{e}(g)*pars.sr):t0+(exp.gapDur{e}(g)*pars.sr)+1.2*pars.sr));
                    
                    all_cppsim(e,ce,cc,g,:) = (cppsim(t0-0.2*pars.sr:t0-0.2*pars.sr+2*pars.sr));
                                        
                    all_p1mbl(e,ce,cc,g,:) = (p1mbl(t0-0.2*pars.sr:t0-0.2*pars.sr+2*pars.sr));
                    all_p2mbl(e,ce,cc,g,:) = (p1mbl(t0+(exp.gapDur{e}(g)*pars.sr):t0+(exp.gapDur{e}(g)*pars.sr)+1.2*pars.sr));
                    
                    
                end
                title(thisGap);
                xticks([t0,t0+pars.sr, t0+pars.sr*2]); xticklabels({'0', '1', '2'});
                xlabel('Time (s)');
                if g == 1; ylabel([{exp.thisExpName}, {'Simulated DV'}]), else ylabel('Simulated DV');end
                xline([t0,t0+0.2*pars.sr+(exp.gapDur{e}(g)*pars.sr)], '-k')
                
            end
        end
    end
end

%% CPP-P1
sortBy.acc = 0; 
if sortBy.acc == 1, ceIdxx = [1:2]; else ceIdxx = 1; end

ylims = ([0 0.035]); xlims = [1,1.2*100];
for ce = ceIdxx(1:end)
    ls = lineTypes{ce};
    for e = 1:2
        f = figure;
        gIdx = 1:length(exp.gapDur{e}); %Include all gaps in plot
        for p1 = 1:2
            if p1 == 1; cc = find(contains(allLabs, 'P1high'));
            else cc = find(contains(allLabs, 'P1low')); end
            thisCol = flip(plt.satPalettes{1,p1*p1});
            
            subplot(2,2,3) %Plot CPP-P1
            ylims = [0,0.035+pars.decayDelay/1000]; ylim(ylims);
            if sortBy.acc == 1
                plot(squeeze(nanmean(nanmean(all_p1sim(e,ce,cc,gIdx,1:end),4),3)),'LineWidth', 2, 'Color', thisCol, 'LineStyle', ls);
            else
                plot(squeeze(nanmean(nanmean(nanmean(all_p1sim(e,:,cc,gIdx,1:end),4),3),2)),'LineWidth', 2, 'Color', thisCol, 'LineStyle', ls);
            end
            hold on;
            %         ylabel('a.u.');
            title('CPP-P1')
            set(gca, 'fontsize', 12);
            hold on;xlabel('Time (s)');
            xticks([0.2*pars.sr, pars.sr/2+0.2*pars.sr, pars.sr+0.2*pars.sr]); xticklabels({'P1', '0.5', '1'});
            xline(0.2*pars.sr)
            ylim(ylims); xlim(xlims);
            set(gca, 'fontsize', 12)
            
            subplot(2,2,1) %Plot MBL-P1
            ylims = [0,0.05]; ylim(ylims);
            if sortBy.acc == 1
                plot(squeeze(nanmean(nanmean(all_p1mbl(e,ce,cc,gIdx,1:end),4),3)),'LineWidth', 2, 'Color', thisCol, 'LineStyle', ls);
            else
                plot(squeeze(nanmean(nanmean(nanmean(all_p1mbl(e,:,cc,gIdx,1:end),4),3),2)),'LineWidth', 2, 'Color', thisCol, 'LineStyle', ls);
            end
            hold on;
            title('MBL-P1'); xlim(xlims);
            
            set(gca, 'fontsize', 12, 'ydir', 'reverse')
            xticks([0.2*pars.sr, pars.sr/2+0.2*pars.sr, pars.sr+0.2*pars.sr]); xticklabels({'P1', '0.5', '1'});
            xline(0.2*pars.sr)
            ylim(ylims);
            
        end
        
        gIdx = 1:length(exp.gapDur{e});
        cc = 5;
        c = 0;
        for p1 = 1:2
            for p2 = 1:2
                c = c+1;
                cc = cc-1;
                subplot(2,2,4) %CPP-P2
                ylims = [0,0.015]; ylim(ylims);
                if sortBy.acc == 1
                    plot(squeeze(nanmean(nanmean(all_p2sim(e,ce,cc,gIdx,1:end),4),3)),'LineWidth', 2, 'Color', plt.Rcolors{cc}, 'LineStyle', ls);
                else
                    plot(squeeze(nanmean(nanmean(nanmean(all_p2sim(e,:,cc,gIdx,1:end),4),3),2)),'LineWidth', 2, 'Color', plt.Rcolors{cc}, 'LineStyle', ls);
                end
                %             plot(squeeze(nanmean(all_p2sim(e,ce,cc,gIdx,1:end),4)),'LineWidth', 2, 'Color', plt.Rcolors{cc}, 'LineStyle', ls);
                xlim([1,1.2*pars.sr])
                hold on;xlabel('Time (s)');
                xticks([0.2*pars.sr, pars.sr/2+0.2*pars.sr, pars.sr+0.2*pars.sr]); xticklabels({'P2', '0.5', '1'});
                xline(0.2*pars.sr);     set(gca, 'fontsize', 12)
                title('CPP-P2')
                
                
                subplot(2,2,2) %MBL-P2
                ylims = [0,0.05]; ylim(ylims);
                if sortBy.acc == 1
                    plot(squeeze(nanmean(nanmean(all_p2mbl(e,ce,cc,gIdx,1:end),4),3)),'LineWidth', 2, 'Color', plt.Rcolors{cc}, 'LineStyle', ls);
                else
                    plot(squeeze(nanmean(nanmean(nanmean(all_p2mbl(e,:,cc,gIdx,1:end),4),3),2)),'LineWidth', 2, 'Color', plt.Rcolors{cc}, 'LineStyle', ls);
                end
                xlim([1,1.2*pars.sr])
                hold on;
                xticks([0.2*pars.sr, pars.sr/2+0.2*pars.sr, pars.sr+0.2*pars.sr]); xticklabels({'P2', '0.5', '1'});
                xline(0.2*pars.sr);      set(gca, 'fontsize', 12, 'ydir', 'reverse')
                title('MBL-P2')
                
                
            end
        end
        xticks([0.2*pars.sr, pars.sr/2+0.2*pars.sr, pars.sr+0.2*pars.sr]); xticklabels({'P2', '0.5', '1'});
        xline(0.2*pars.sr)
        %     ylabel('a.u.')
        
        set(gca, 'fontsize', 12)
        f.Units = 'centimeters';
        f.OuterPosition = [0.5,0.5,14,15];
        exportgraphics(f, [exp.figPath, 'GT_Fig5B_Exp' num2str(e) '.tiff'], 'Resolution', 600)
        exportgraphics(f, [exp.figPath, 'GT_Fig5B_Exp' num2str(e) '.pdf'], 'Resolution', 600)
    end
end

%% Plot averages over gaps 
clear thisCPP;  
for e = 1:2
     f=figure
     lineTypes = {'-', '--'}
     for ce = ceIdxx(1:end) %Loop through accuracies
         lt = lineTypes{ce};
         for gap = 1:length(exp.gapDur{e})
             cols = (plt.satPalettes);
             thisCol = cols(2,gap);
             
             thisCPP(ce,gap,:,:) = squeeze(nanmean(all_cppsim(e,ce,:,gap,:),3));
             
             subplot(2,1,2) %PLOT CPP (|deltaDV|)
             ylims = [0 0.035]
             baseline = 0;
             if ce == ceIdxx(end)
                 plot(squeeze(nanmean(thisCPP(:,gap,:,:),1)-baseline),'LineWidth', 2, 'Color',thisCol{:}, 'LineStyle', '-');
             end
             
             xlim([1,2.2*pars.sr])
             hold on;
             fill(plt.gapShades{e}(gap+1,:)*ds,[ylims, flip(ylims)], thisCol{:}, 'LineStyle','none','facealpha', 0.1); %Plot pulses in time
             hold on;yline(0,'--k'); hold on;
             if gap == 1;  fill(plt.gapShades{e}(1,:)*ds,[ylims, flip(ylims)], 'k', 'LineStyle','none','facealpha', 0.1);end %Plot pulses in time
             ylabel('Simulated CPP');
             xticks([0.2*pars.sr, pars.sr+0.2*pars.sr, pars.sr*2.2]); xticklabels({'0', '1',  '2'});xline(0.2*pars.sr); set(gca, 'fontsize', 12);
                     ylim(ylims);    xlabel('Time (s)');
                     box on
             
             thisMBL = squeeze(nanmean(all_p1mbl(e,ce,:,gap,:),3));
             
             subplot(2,1,1) 
             baseline = 0;        ylims = [-0.05 0.05]
             plot(squeeze((thisMBL(:,1:end))-baseline),'LineWidth', 2, 'Color',thisCol{:}, 'LineStyle', lt);
             xlim([1,2.2*pars.sr])
             hold on;
             fill(plt.gapShades{e}(gap+1,:)*ds,[ylims, flip(ylims)], thisCol{:}, 'LineStyle','none','facealpha', 0.1); %Plot pulses in time
             hold on;yline(0,'--k'); hold on;
             if gap == 1;  fill(plt.gapShades{e}(1,:)*ds,[ylims, flip(ylims)], 'k', 'LineStyle','none','facealpha', 0.1);end %Plot pulses in time
             ylabel('Simulated MBL');    xticks([0.2*pars.sr, pars.sr+0.2*pars.sr, pars.sr*2.2]); xticklabels({'0', '1',  '2'});xline(0.2*pars.sr); set(gca, 'fontsize', 12);
             ylim(ylims);set(gca, 'ydir', 'reverse'); %So that errors go up, as in fig. 2
             yline([-bestpm(e,3,1),bestpm(e,3,1)], ':r', 'linewidth', 1) %Plot bound 
             box on
         end
     end
     f.Units = 'centimeters';
     f.OuterPosition = [0.25 0.35 11 15];
     exportgraphics(f, [exp.figPath, 'GT_Fig5A_Exp' num2str(e) '.tiff'], 'Resolution', 600)
     exportgraphics(f, [exp.figPath, 'GT_Fig5A_Exp' num2str(e) '.pdf'], 'Resolution', 600)    
end

%% Plot Figure S8 - simulate CPP with no decay during gap 
%Fot these plots, set decayBlue to = 0, sortBy.acc = 0 in simulations
%above

clear thisCPP;      
f=figure

for e = 1:2
     lineTypes = {'-', '--'}
         for gap = 1:length(exp.gapDur{e})
             cols = (plt.satPalettes);
             thisCol = cols(2,gap);
             
             thisCPP(gap,:,:) = squeeze(nanmean(nanmean(all_cppsim(e,:,:,gap,:),3),2));
             
             subplot(2,1,e) %PLOT CPP (|deltaDV|)
             ylims = [0 0.055]
             baseline = 0;
            
             plot(squeeze(thisCPP(gap,:,:)-baseline),'LineWidth', 2, 'Color',thisCol{:}, 'LineStyle', '-');
             
             xlim([1,2.2*pars.sr])
             hold on;
             fill(plt.gapShades{e}(gap+1,:)*ds,[ylims, flip(ylims)], thisCol{:}, 'LineStyle','none','facealpha', 0.1); %Plot pulses in time
             hold on;yline(0,'--k'); hold on;
             if gap == 1;  fill(plt.gapShades{e}(1,:)*ds,[ylims, flip(ylims)], 'k', 'LineStyle','none','facealpha', 0.1);end %Plot pulses in time
             ylabel('Simulated CPP');
             xticks([0.2*pars.sr, pars.sr+0.2*pars.sr, pars.sr*2.2]); xticklabels({'0', '1',  '2'});xline(0.2*pars.sr); set(gca, 'fontsize', 12);
                     ylim(ylims);    xlabel('Time (s)');
                     box on
     end
     f.Units = 'centimeters';
     f.OuterPosition = [0.25 0.35 11 15];

end

%% Plot Figure S6  - simulated CPP sorted by error/correct
% Plot 
ylims = ([0 0.035]); xlims = [1,1.2*100];
for ce = ceIdxx(1:end)
    f = figure;
    ls = lineTypes{ce};
    for e = 1:2
        gIdx = 1:length(exp.gapDur{e}); %Include all gaps in plot
        for p1 = 1:2
            if p1 == 1; cc = find(contains(allLabs, 'P1high'));
            else cc = find(contains(allLabs, 'P1low')); end
            thisCol = flip(plt.satPalettes{1,p1*p1});
            
            subplot(2,2,1+(2*(e-1))) %Plot CPP-P1
            ylims = [0,0.06]; ylim(ylims);
            if sortBy.acc == 1
                plot(squeeze(nanmean(nanmean(all_p1sim(e,ce,cc,gIdx,1:end),4),3)),'LineWidth', 2, 'Color', thisCol, 'LineStyle', ls);
            else
                plot(squeeze(nanmean(nanmean(nanmean(all_p1sim(e,:,cc,gIdx,1:end),4),3),2)),'LineWidth', 2, 'Color', thisCol, 'LineStyle', ls);
            end
            hold on;
            %         ylabel('a.u.');
                if e == 1; title('CPP-P1'); end
            set(gca, 'fontsize', 12);
            hold on;xlabel('Time (s)');
            xticks([0.2*pars.sr, pars.sr/2+0.2*pars.sr, pars.sr+0.2*pars.sr]); xticklabels({'P1', '0.5', '1'});
            xline(0.2*pars.sr)
            ylim(ylims); xlim(xlims);
            %     ylabel('a.u.');
            set(gca, 'fontsize', 12)
            
        end
        
        gIdx = 1:length(exp.gapDur{e});
        cc = 5;
        c = 0;
        for p1 = 1:2
            for p2 = 1:2
                c = c+1;
                cc = cc-1; hold on;
            subplot(2,2,2+(2*(e-1))) %Plot CPP-P2
                ylims = [0,0.03]; ylim(ylims);
                if sortBy.acc == 1
                    plot(squeeze(nanmean(nanmean(all_p2sim(e,ce,cc,gIdx,1:end),4),3)),'LineWidth', 2, 'Color', plt.Rcolors{cc}, 'LineStyle', ls);
                else
                    plot(squeeze(nanmean(nanmean(nanmean(all_p2sim(e,:,cc,gIdx,1:end),4),3),2)),'LineWidth', 2, 'Color', plt.Rcolors{cc}, 'LineStyle', ls);
                end
                xlim([1,1.2*pars.sr])
                hold on;xlabel('Time (s)');
                xticks([0.2*pars.sr, pars.sr/2+0.2*pars.sr, pars.sr+0.2*pars.sr]); xticklabels({'P2', '0.5', '1'});
                xline(0.2*pars.sr);     set(gca, 'fontsize', 12)
                if e == 1; title('CPP-P2'); end
                
               
                              
            end
        end
        xticks([0.2*pars.sr, pars.sr/2+0.2*pars.sr, pars.sr+0.2*pars.sr]); xticklabels({'P2', '0.5', '1'});
        xline(0.2*pars.sr)
        %     ylabel('a.u.')
        
        set(gca, 'fontsize', 12)
        f.Units = 'centimeters';
        f.OuterPosition = [0.5,0.5,13,15];
        
    end
end