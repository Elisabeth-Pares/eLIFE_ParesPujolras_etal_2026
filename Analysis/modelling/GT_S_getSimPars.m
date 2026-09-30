function pars = GT_S_getSimPars(exp, pars, modelName)

%NOTE: not all options are implemented in the paper or the subsequent
%simulation code 
%% General
if contains(modelName, '1DV')
    pars.biased = 1; %Assume one single decision, where accumulation at P2 starts off from the final state of the DV at the end of P1. 
else
    pars.biased = 0; %Assume two independent decisions;
end

% Do NOT accumulate noise over gap --> would predict decreasing accuracy over time 
% If set to = 1, the model would continue accumulating during the
% incoherent motion gap. 
pars.accGap = 0; 

%Sampling rate for the simulations
pars.sr = 100; 

%% Bound  - To allow bound to step down after P1, or linearly decay over time [not implemented]
if contains(modelName, 'StepBound'); pars.boundDecay = 'step'; 
elseif contains(modelName, 'LinBound'); pars.boundDecay = 'linear';
else pars.boundDecay = 'NA';
end

%% Bias -  To introduce a bias to the accumulation at P2 [not implemented]
if contains(modelName, 'StepBias'); pars.bias = 'zBias'; 
elseif contains(modelName, 'LinBias'); pars.bias = 'zBias_linear'; 
elseif contains(modelName, 'RevLBias');  pars.bias = 'zBias_linearReverse'; 
else pars.bias = 'NA'; 
end

%% DR options - To estimate 2 drift rates (high/low coherence separately). If set to NA, force linear scaling.
if contains(modelName, '2DR')
    pars.drModulation = '2DR';

end

%% DR 
% if contains(modelName, 'StepDR'); pars.drModulation = 'step'; 
% elseif contains(modelName, 'LinDR'); pars.drModulation = 'linear'; 
% elseif contains(modelName, 'BiasDR');  pars.drModulation = 'bias'; 
% elseif contains(modelName, 'invprob_global'); pars.drModulation = 'invprob'; 
%     pars.globalgain = 1; 
%     elseif contains(modelName, 'invprob_DR'); pars.drModulation = 'invprob'; 
%     pars.globalgain = 0; 
% elseif contains(modelName, 'alphaDR')
%     pars.drModulation = 'alpha';
% elseif contains(modelName, '2DR')
%     pars.drModulation = '2DR';
% elseif contains(modelName, 'invDRP2')
%     pars.drModulation = 'invDRP2';
% else pars.drModulation = 'NA'; 
%     
% end
if contains(modelName, 'invDRP2')
    pars.drModulation = 'invDRP2';
else pars.drModulation = 'NA'; 
    
end



