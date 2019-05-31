function Out = RLdecayTmaze1(num_trial,RLtype,RLparas,Rews,decay_rate,DAdep_paras)

%-----
% This file is associated with the following article, which has been provisionally accepted for publication in PLOS Computational Biology
% (initially submitted on May 11, 2016, and provisionally accepted on Sep 14, 2016):
% Authors: Ayaka Kato (1) & Kenji Morita (2)
% Affiliations:
%  (1) Department of Biological Sciences, Graduate School of Science, The University of Tokyo, Tokyo, Japan
%  (2) Physical and Health Education, Graduate School of Education, The University of Tokyo, Tokyo Japan
% Title: Forgetting in Reinforcement Learning Links Sustained Dopamine Signals to Motivation
% Short title: Dynamic Equilibrium in Reinforcement Learning
% Correspondence: Kenji Morita (morita@p.u-tokyo.ac.jp)
%-----
% Out = RLdecayTmaze1(num_trial,RLtype,RLparas,Rews,decay_rate,DAdep_paras)
%
%   DAdep_paras: [DAdep_factor (mulitiplicative), DAdep_start_trial], e.g., [0.25, 501]
%
% <states>	1       2       3       4           5/6     7/8     9
% <actions> 1S2G	4S5G	7S8G	10S11/12	13S14G	19S20G      Arm1
%                                               16S17G	22S23G      Arm2
%

% parameters
p_alpha = RLparas(1);
p_beta = RLparas(2);
p_gamma = RLparas(3);
num_tstep = num_trial * 100;

% initialization of the variables
ArmChoices = zeros(num_trial,1);
States = zeros(num_tstep,1);
Choices = zeros(num_tstep,1);
Qs = zeros(num_trial,23);
TDs = zeros(num_tstep,1);
endsteps = zeros(num_trial,1); % time steps for reaching the trial end

% main loop
Qnow = zeros(1,23);
CurrS = 1;
prevA = [];
k_trial = 1;
for k_tstep = 1:num_tstep
    
    % save the state
    States(k_tstep) = CurrS;

    % set DA depletion
    if k_trial >= DAdep_paras(2)
        DAfactor = DAdep_paras(1);
    else
        DAfactor = 1;
    end
    
    % main
    if isempty(prevA)
        Choices(k_tstep) = 0 + actselect(Qnow([1,2]),p_beta);
        if RLtype == 'Q'
            TDs(k_tstep) = 0 + p_gamma*max(Qnow([1,2])) - 0;
        elseif RLtype == 'S'
            TDs(k_tstep) = 0 + p_gamma*Qnow(Choices(k_tstep)) - 0;
        end
        Qnow = Qnow * (1 - decay_rate);
        prevA = Choices(k_tstep);
        if Choices(k_tstep) == 2
            CurrS = 2;
        end
    elseif CurrS == 4
        Choices(k_tstep) = 9 + actselect(Qnow([10,11,12]),p_beta);
        if RLtype == 'Q'
            TDs(k_tstep) = Rews(prevA) + p_gamma*max(Qnow([10,11,12])) - Qnow(prevA);
        elseif RLtype == 'S'
            TDs(k_tstep) = Rews(prevA) + p_gamma*Qnow(Choices(k_tstep)) - Qnow(prevA);
        end
        Qnow(prevA) = Qnow(prevA) + DAfactor*p_alpha*TDs(k_tstep);
        Qnow = Qnow * (1 - decay_rate);
        prevA = Choices(k_tstep);
        if Choices(k_tstep) == 11
            ArmChoices(k_trial) = 1;
            CurrS = 5;
        elseif Choices(k_tstep) == 12
            ArmChoices(k_trial) = 2;
            CurrS = 6;
        end
    elseif CurrS == 9 % trial end
        TDs(k_tstep) = Rews(prevA) + p_gamma*0 - Qnow(prevA);
        Qnow(prevA) = Qnow(prevA) + DAfactor*p_alpha*TDs(k_tstep);
        Qnow = Qnow * (1 - decay_rate);
        Qs(k_trial,:) = Qnow;
        endsteps(k_trial) = k_tstep;
        if k_trial == num_trial
            States = States(1:k_tstep);
            Choices = Choices(1:k_tstep);
            TDs = TDs(1:k_tstep);
            break;
        else
            k_trial = k_trial + 1;
            CurrS = 1;
            prevA = [];
        end
    else
        Choices(k_tstep) = 3*(CurrS-1) + actselect(Qnow(3*(CurrS-1)+[1,2]),p_beta);
        if RLtype == 'Q'
            TDs(k_tstep) = Rews(prevA) + p_gamma*max(Qnow(3*(CurrS-1)+[1,2])) - Qnow(prevA);
        elseif RLtype == 'S'
            TDs(k_tstep) = Rews(prevA) + p_gamma*Qnow(Choices(k_tstep)) - Qnow(prevA);
        end
        Qnow(prevA) = Qnow(prevA) + DAfactor*p_alpha*TDs(k_tstep);
        Qnow = Qnow * (1 - decay_rate);
        prevA = Choices(k_tstep);
        if Choices(k_tstep) == 3*(CurrS-1) + 2 % Go
            if (CurrS == 1) || (CurrS == 2) || (CurrS == 3) || (CurrS == 8)
                CurrS = CurrS + 1;
            elseif (CurrS == 5) || (CurrS == 6) || (CurrS == 7)
                CurrS = CurrS + 2;
            end
        end
    end
end
if k_tstep == num_tstep
    error('number of time steps is not enough');
end

% output variables
Out.ArmChoices = ArmChoices;
Out.States = States;
Out.Choices = Choices;
Out.Qs = Qs;
Out.TDs = TDs;
Out.endsteps = endsteps;


% innter-file function

function chosen_option_index = actselect(optionQs,p_beta)

Pchoose = zeros(1,length(optionQs));
for k = 1:length(optionQs)
    Pchoose(k) = exp(p_beta*optionQs(k))/sum(exp(p_beta*optionQs));
end

chosen_option_index = [];
tmp_rand = rand;
tmp = 0;
for k = 1:length(optionQs)
    tmp = tmp + Pchoose(k);
    if tmp_rand <= tmp
        chosen_option_index = k;
        break;
    end
end
if isempty(chosen_option_index)
    error('choice is not made');
end
