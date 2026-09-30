%% GT MODELLING FOR SONIC

% Specify which model to run
N = 10000; 

%parList = {'b', 'bc', 'dlow', 'dhigh', 'z', 'accT', 'accTvar', 'eta', 'sz', 'sd', 'N', 'seed'};
parList = {'b', 'bc', 'dlow', 'dhigh', 'dp2', 'accT', 'accTvar', 'eta', 'sz', 'sd', 'N', 'seed'};

modelNames = {'FixedBound_1DV_2DR',...
 'FixedBound_1DV_2DR_lowDRP2', 'FixedBound_1DV_2DR_invDRP2',...
 'NoBound_1DV_2DR_lowDRP2', 'NoBound_1DV_2DR_invDRP2'};

Selvecs = {...
    [1 0 1 1 0 0 0 0 0 0 0 0],...;%FixedBound_1DV_2DR (3) - FINAL 
    [1 0 1 1 1 0 0 0 0 0 0 0],...%FixedBound_1DV_2DR_lowDRP2(9) %R2: scale P2 drift rate 
    [1 0 1 1 1 0 0 0 0 0 0 0],...%FixedBound_1DV_2DR_invDRP2(10) %R2: scale P2 drift rate inversely to P1
    [0 0 1 1 1 0 0 0 0 0 0 0],...%NoBound_1DV_2DR_lowDRP2(11) %R2: scale P2 drift rate 
    [0 0 1 1 1 0 0 0 0 0 0 0]}...%NoBound_1DV_2DR_invDRP2(12) %R2: scale P2 drift rate inversely to P1 

% Fit for both experiments
for e = 1:2 
    % New seed!
    rng('shuffle'); seed = sum(clock); rng(seed);
    
    %% Load specific experiment settings 
    exp.thisExp = e; %Experiment to analyse
    [exp, pars]  = GT_S_setup(exp); %Get correct experiment parameters
    
    pars.nSim = N;
    pars.simpleFit = 0; % 0: all, 1: 1p and 0gap, 2: 1p, 0 and long gap; 3 = 1p and long only; 

    %Load grand-averaged data 
    thisFile = fullfile(".", ['GT_gaDat_' exp.thisExpCode '.mat']); 
    load(thisFile)
  
    %% Parameter specifications
    par_range = [0.01, 0.1;...%0.1;...%0.01 0.1;...%0 1;...    % Bound - before 1
        0.01, 0.1;...       % Rate of bound collapse/increase [0.000125 to 0.002] - allow it to be 0 - would indicate fixed bound
        0.001 0.25;...         % Drift rate [Restricted from .001 -1 to .001 0.25 to help fitting]
        0.001 0.25;...         % Drift rate bias - option to increase DR if bound not hit after DV1
        0 2;...            % dp2: allow drift rate at P2 to be systematically lower/higher than in P1 
        0 1;...             % AccT
        0 1;...             % AccTvar
        0 1;                % eta
        -0.4 0.4;...        % sz
        -0.4 0.4;...        % sd
        N N;...
        seed seed];

%     par_range = [par_range; ...
%         N N;...
%         seed seed];
    
    pm_if_fixed = [Inf 0 0.15 0 1 0 0 0 0 0 N seed];
    
    %% Fit model
    p = [];
    for m=1
        Sel = Selvecs{m};
        assert(length(Sel) == size(par_range,1))
        numfreepm=sum(Sel);
        disp(['Fitting ' modelNames{m} ' model, for Experiment ' num2str(e), ', fitType' num2str(pars.simpleFit)])
        
        options = optimoptions('particleswarm','UseParallel',true, 'SwarmSize', 10*sum(Sel),'MaxStallIterations', 40); %Edited on 7/10/24 - tried to extend max stall iterations to try avoiding local minima
        
        [p, g2(e,m)] = particleswarm(@(p) GT_S_GSquare(p,exp,pars,pm_if_fixed(Sel == 0),Sel,ga, modelNames{m}), sum(Sel), par_range((Sel==1),1)',par_range((Sel==1),2)', options);
        bestpm(e,m,Sel==1) = [p]; 
        bestpm(e,m,[Sel==0]) = pm_if_fixed(Sel==0);
        
        disp(['best parameters: ' num2str(bestpm(e,m,:))])
        AIC(e,m) = g2(e,m) + 2*numfreepm;
        BIC(e,m) = g2(e,m) + numfreepm*log(sum(ga.n(:)));
    end
    
end
% Save model fit
savingFolder = fullfile(".","GT_S_ModelFits");
if ~isfolder(savingFolder)
    mkdir(savingFolder)
end
fitType = ['_fitType' num2str(pars.simpleFit)];
save(fullfile(savingFolder, [modelNames{m}, fitType]), 'bestpm', 'AIC', 'BIC', 'modelNames', 'Selvecs', 'par_range')

