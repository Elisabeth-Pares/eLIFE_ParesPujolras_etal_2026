function [exp, pars] = GT_S_setup(exp)

%% EXPERIMENT & PARTICIPANT INFO 
exp.code            = {'EP1', 'EP2'};
exp.name            = {'Experiment 1', 'Experiment 2'};
exp.numbers         = [1,2];

exp.sub_id{1}            = [1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]; % P8 excluded because only 4 blocks were recorded.
exp.sub_id{2}            = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,23]; % P1 & P21 excluded because they could not perform task well


exp.thisExpName= exp.name{exp.thisExp};
exp.thisExpCode= exp.code{exp.thisExp};

%% Experimental parameters 
exp.pulseDur = 0.2; % in s

%Gap duration is the only parameter varying in the two tasks
exp.gapDur{1} = [0,0.5,1]; %in seconds
exp.gapDur{2} = [0,0.120, 0.360, 1.080]; %in seconds

exp.gapLabs{1} = {'gap0', 'gap05', 'gap1'};
exp.gapLabs{2} = {'gap0', 'gap120', 'gap360', 'gap1080'};

exp.respKeys = [67,77]; %Left (c), right (m)
exp.accuracies = [0,1]; %Wrong, Right
exp.dotDirs = [0,-180]; %0: - Right, -180: left 

exp.cohLabs = {'high', 'low'};
exp.p1Labs = {'p1coh_high', 'p1coh_low'};
exp.p2Labs = {'p2coh_high', 'p2coh_low'};
exp.p1cohLabs = {'p1high', 'p1low'};
exp.p2cohLabs = {'p2high', 'p2low'};
exp.dirLabs = {'dotRight', 'dotLeft'};
exp.movLabs = {'Left', 'Right'};
exp.accLabs = {'Incorrect', 'Correct'};
exp.npulseLabs = {'onePulse', 'twoPulse'};
exp.cppLabs = {'Early', 'Late'};
exp.alphaLabs = {'Low', 'High'}; %PreStimulus power 
exp.repLabs = {'Alt', 'Rep'}; %Alternation/Repetition response

exp.blockLength = 70; 

%% SIMULATION PARAMETERS 
pars.sr = 100;      %Simulation sampling rate
pars.nSim= 1000;    %Number of simulations
pars.maxIter =100;  %Max iterations for fitting 
pars.nRepMat=100;   
pars.stepAfterbndT = 0.1*pars.sr; %%When is activity reset back to 0 after hitting bound? 0.1s*sr. Originally: 50 (for 500 sr)
pars.accType = 'erf'; %How to compute accuracy 

%Get condition labels for this experiment 
c = 0;
% Do for 1pulse
for p1 = 1:2
    c = c+1;
    
    tag = ['onePulse_P1' exp.cohLabs{p1}];
    %Update data structure
    pars.cLab{c} = tag;
    pars.gDur(c) = -0.2;%1000;%exp.gapDur{exp.thisExp}(g);
end

%Sort by conditions (p1coh x p2coh x gap)
for g = 1:length(exp.gapDur{exp.thisExp})
    for p1 = 1:2
        for p2 = 1:2
            c = c+1;
            
            tag = [exp.gapLabs{exp.thisExp}{g}, '_P1' exp.cohLabs{p1}, '_P2' exp.cohLabs{p2}];
            pars.cLab{c} = tag;
            pars.gDur(c) = exp.gapDur{exp.thisExp}(g);
            clear t; clear idx; clear tag;
        end
    end
end

