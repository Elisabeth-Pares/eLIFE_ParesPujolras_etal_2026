%% GT MODELLING FOR SONIC
%Fit simple bounded drift diffusion model

N = 10000; % Number of simulations to fit model 

%Parameter list [not all possible free parameters were used for the paper]
parList = {'b', 'bc', 'd', 'db', 'z', 'accT', 'accTvar', 'eta', 'sz', 'sd', 'N', 'seed'};

%Names of the models to fit. 
%Note: Names are important, because they are used in
%the GT_getSimPars function later on to determine some aspects of the
%accumulation process. For example, 1DV indicates that decisions are based
%on one single DV that extends all the way from P1 to P2, rather than the 
%assuming that choices are based on the sum of two independent decision
%processes where P1 + P2 are fully integrated. 

modelNames = {'FixedBound_1DV_2DR', 'FixedBound_1DV_2DR_DRvar'};
% Vectors
Selvecs = {...
    [1 0 1 1 0 0 0 0 0 0 0 0],...;%FixedBound_1DV_2DR 
    [1 0 1 1 0 0 0 0 0 1 0 0]}    %FixedBound_2DR_DRvar   

% Fit for both experiments
for e = 1:2 
    % New seed!
    rng('shuffle'); seed = sum(clock); rng(seed);
    
    %% Load specific experiment settings 
    exp.thisExp = e; %Experiment to analyse
    [exp, pars]  = GT_S_setup(exp); %Get correct experiment parameters
    
    pars.nSim = N;
    pars.simpleFit = 0; %  Option to fit model to subsets of trials. 0: all, 1: 1p and 0gap, 2: 1p, 0 and long gap; 3 = 1p and long only; 

    %Load grand-averaged data 
    thisFile = fullfile(".", ['GT_gaDat_' exp.thisExpCode '.mat']); 
    load(thisFile)
  
    %% Parameter specifications
    par_range = [0.01, 0.1;...%0.1;...%0.01 0.1;...%0 1;...    % Bound - before 1
        0.01, 0.1;...       % Rate of bound collapse/increase [0.000125 to 0.002] - allow it to be 0 - would indicate fixed bound
        0.001 1;...         % Drift rate 1 
        0.001 1;...         % Drift rate 2
        -2 2;...            % Bias: starting point shift
        0 1;...             % AccT: time of accumulation onset
        0 1;...             % AccTvar; variability in accT
        0 1;                % eta
        -0.4 0.4;...        % sz: variabiliyt in starting point shift
        -0.4 0.4;...        % sd
        N N;...
        seed seed];
    
    pm_if_fixed = [1 0 0 0 0 0 0 0 0 0 N seed]; %Which default values to use if parameters are not free
    
    %% Fit model
    p = [];
    for m=1
        Sel = Selvecs{m};
        assert(length(Sel) == size(par_range,1))
        numfreepm=sum(Sel);
        disp(['Fitting ' modelNames{m} ' model, for Experiment ' num2str(e), ', fitType' num2str(pars.simpleFit)])
        
        options = optimoptions('particleswarm','UseParallel',false, 'SwarmSize', 10*sum(Sel),'MaxStallIterations', 40); %Edited on 7/10/24 - tried to extend max stall iterations to try avoiding local minima
        pars.acc = 'erf';
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

