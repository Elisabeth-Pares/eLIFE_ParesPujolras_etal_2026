function [g2] = GT_S_GSquare(pvar, exp,pars,pfix, Sel,dat, modelName)
%% Observed data:
p_obs_all = dat.p;
p_obs = [];
% Put in format compatible with simulations
for c = 2:-1:1 %OnePulse
     if c <3
        p_obs = [p_obs; p_obs_all(:,1,c)];
     end
end
for g = 1:length(exp.gapDur{exp.thisExp})
    for c = 6:-1:3
        p_obs = [p_obs; p_obs_all(:,g,c)];
    end
end

p_obs = p_obs'; 
p_obs(2,:) = 1-p_obs(1,:);

%% Get analysis parameters based on modelName
% Specifies how to run the accumulation process 
pars = GT_S_getSimPars(exp, pars, modelName); 

%% Simulate data 
[simdata, pars] = GT_S_simul_DDM_clean(pars, pvar, pfix, Sel);

cond = simdata.tt(:,1); 
outcome = simdata.RT(:,outcomeIdx);

% Now we're going to make a matrix called 'pred' that computes the proportion of simulated trials of a certain condition that fall
% in each of outcome categories of the real data
% pred will be outcome x condition
for c=1:size(p_obs,2)   %  number of conditions

    numtr_thiscond = length(find(cond==c));
    % order of outcomes in each column of dat and pred: corr, err, miss
    if strcmp(pars.acc, 'bin'); %If fitting based on binary accuracy (-1 error, 1 correct)
        p_pred(2,c) = length(find(cond==c & outcome==-1))./numtr_thiscond; % errors
        p_pred(1,c) = length(find(cond==c & outcome==1))./numtr_thiscond; % corrects
    elseif strcmp(pars.acc, 'erf') %If fitting based on probabilistic accuracy 
        p_pred(2,c) = length(find(cond==c & outcome<=0.5))./numtr_thiscond; % errors
        p_pred(1,c) = length(find(cond==c & outcome>0.5))./numtr_thiscond; % corrects
    end
    % correction for any p_pred=0 (to avoid div by zero)
    p_pred(find(p_pred(:,c)<0.0001 & p_obs(:,c)>0),c) = 0.0001; % only apply correction where needed
    p_pred(:,c) = p_pred(:,c)./sum(p_pred(:,c)); % re-normalise so all predicted proportions sum to 1 again
end

% compute G-square
g2vals=[];
for c=1:size(p_pred,2) % condition
    for o=1:size(p_pred,1)   % error/correct
        if p_obs(o,c)>0
            g2vals = [g2vals sum(dat.n(:,c))*p_obs(o,c)*log(p_obs(o,c)/p_pred(o,c))]; 
        end
    end
end

g2=nansum(g2vals)*2; 
