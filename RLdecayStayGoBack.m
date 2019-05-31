function Out = RLdecayStayGoBack(num_trial,RLtype,RLparas,decay_rate,DAdep_paras)

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
% Out = RLdecayStayGoBack(num_trial,RLtype,RLparas,decay_rate,DAdep_paras);
%
% <input variables>
%   num_trial: number of trials
%   RLtype: 'Q':Q-learning, 'S':SARSA
%	RLparas: [p_alpha, p_beta, p_gamma]
%   decay_rate: rate of the decay of the learned values per time-step, e.g., 0.01
%   DAdep_paras: [DAdep_factor (mulitiplicative), DAdep_start_trial], e.g., [0.25, 251]
% <output variable>
%   Out.TDs: TD error at each time step
%   Out.States: tansitions of states
%   Out.goalsteps: index of time step when reaching the goal
%   Out.Vs_whole: learned values of all the actions evaluated at the end of each trial

% number of states
num_state = 7;

% maximum number of time steps to simulate
maxstep = num_trial * num_state * 20;

% set random numbers
rands_for_action = rand(1,maxstep);

% RL parameters
p_alpha = RLparas(1);
p_beta = RLparas(2);
p_gamma = RLparas(3);

% DA depletion parameters
DAdep_factor = DAdep_paras(1); % the size of value update becomes DAdep_factor * p_alpha * TDerror
DAdep_start_trial = DAdep_paras(2); % index of trial from which DA depletion starts

% reward
Rs = [zeros(1,num_state-1), 1];  % reward at each state (time step) for all the trials

% variables to save
TDs = zeros(1,maxstep); % TD error at each time step for all the trials, initialization
States = zeros(1,maxstep); % tansitions of states, initialization
goalsteps = []; % index of time step when reaching the goal, initialization
Vs_whole = []; % learned values of all the actions evaluated at the end of each trial, initialization

% run simulation
current_trial = 1; % initialization
current_state = 1; % initialization
previous_action = []; % initialization
Vs_latest = zeros(3,num_state); % latest learned values of all the actions (row: [1-Stay;2-Go;3-Back], column: at each state), initialization
for k_tstep = 1:maxstep
    States(k_tstep) = current_state;
    if current_state ~= num_state
        tmp_chooseStay = exp(p_beta * Vs_latest(1,current_state));
        tmp_chooseGo = exp(p_beta * Vs_latest(2,current_state));
        if current_state == 1
            tmp_chooseBack = 0;
        else
            tmp_chooseBack = exp(p_beta * Vs_latest(3,current_state));
        end
        if rands_for_action(k_tstep) <= tmp_chooseStay / (tmp_chooseStay + tmp_chooseGo + tmp_chooseBack)
            current_action = 1;
        elseif rands_for_action(k_tstep) <= (tmp_chooseStay + tmp_chooseGo) / (tmp_chooseStay + tmp_chooseGo + tmp_chooseBack)
            current_action = 2;
        else
            current_action = 3;
        end
        if RLtype == 'Q'
            upcoming_value = max(Vs_latest(:,current_state));
        elseif RLtype == 'S'
            upcoming_value = Vs_latest(current_action,current_state);
        end
        if isempty(previous_action) % It is assumed that there is no "previous action" at the beginning of a trial
            TDs(k_tstep) = Rs(current_state) + p_gamma * upcoming_value - 0;
        else
            TDs(k_tstep) = Rs(current_state) + p_gamma * upcoming_value - Vs_latest(previous_action,States(k_tstep-1));
            if current_trial >= DAdep_start_trial
            	tmpTDeffective = TDs(k_tstep) * DAdep_factor;
            else
                tmpTDeffective = TDs(k_tstep);
            end
            Vs_latest(previous_action,States(k_tstep-1)) = Vs_latest(previous_action,States(k_tstep-1)) + p_alpha * tmpTDeffective;
        end
        Vs_latest = Vs_latest * (1 - decay_rate); % decay of learned values
        if current_action == 2 % Go
            current_state = current_state + 1;
        elseif current_action == 3 % Back
            current_state = current_state - 1;
            if current_state == 0
                error('state is 0');
            end
        end
        previous_action = current_action;
    else % when reaching the goal
        goalsteps = [goalsteps; k_tstep];
        TDs(k_tstep) = Rs(current_state) + p_gamma * 0 - Vs_latest(previous_action,States(k_tstep-1)); % It is assumed that there is no "upcoming value" at the goal
        if current_trial >= DAdep_start_trial
            tmpTDeffective = TDs(k_tstep) * DAdep_factor;
        else
            tmpTDeffective = TDs(k_tstep);
        end
        Vs_latest(previous_action,States(k_tstep-1)) = Vs_latest(previous_action,States(k_tstep-1)) + p_alpha * tmpTDeffective;
        Vs_latest = Vs_latest * (1 - decay_rate); % decay of learned values
        Vs_whole = [Vs_whole; Vs_latest(1,:) Vs_latest(2,:) Vs_latest(3,:)];
        if current_trial == num_trial
            TDs = TDs(1:k_tstep);
            States = States(1:k_tstep);
            break;
        end
        current_trial = current_trial + 1;
        current_state = 1;
        previous_action = [];
    end
end
if length(goalsteps) < num_trial
    error('planned number of trials have not been completed');
else
    Out.TDs = TDs;
    Out.States = States;
    Out.goalsteps = goalsteps;
    Out.Vs_whole = Vs_whole;
end
