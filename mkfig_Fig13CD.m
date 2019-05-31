function Out = mkfig_Fig13CD

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

% to use the same random numbers as used in the simulations presented in the figures in the paper
load used_rand_twister_for_Fig13CD
rand('twister',used_rand_twister_for_Fig13CD);

% parameters
maxstep = 1000;
num_trial = 500;
RLtype = 'Q';
RLparas = [0.5,5,1];
decay_rate = 0.01;
DAdep_paras = [1 1001];
middlerew = 0.1;

% number of states
num_state = 7;

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
Rs = [zeros(1,num_state-1), 1]; % reward at each state (time step) for all the trials
Rs(ceil(num_state/2)) = middlerew; % small reward at the middle

% variables to save
TDs = zeros(1,maxstep); % TD error at each time step for all the trials, initialization
States = zeros(1,maxstep); % tansitions of states, initialization
goalsteps = []; % index of time step when reaching the goal, initialization
Vs_whole = []; % learned values of all the actions evaluated at the end of each trial, initialization
Vs_everystep = []; % learned values at every time step

% run simulation
current_state = 1; % initialization
previous_action = []; % initialization
Vs_latest = zeros(1,2*num_state); % latest (i.e., updated at each time step) learned values of all the actions, initialization
for k_tstep = 1:maxstep
    current_trial = size(Vs_whole,1)+1;
    States(k_tstep) = current_state;
    if current_state ~= num_state
        prob_chooseNoGo = 1 / (1 + exp(p_beta * (Vs_latest(2*current_state) - Vs_latest(2*current_state-1))));
        if isnan(prob_chooseNoGo)
            break;
        end
        current_action = 2*current_state - (rands_for_action(k_tstep) <= prob_chooseNoGo);
        if RLtype == 'Q'
            upcoming_value = max(Vs_latest([2*current_state-1,2*current_state]));
        elseif RLtype == 'S'
            upcoming_value = Vs_latest(current_action);
        end
        if isempty(previous_action) % It is assumed that there is no "previous action" at the beginning of a trial
            TDs(k_tstep) = Rs(current_state) + p_gamma * upcoming_value - 0;
        else
            TDs(k_tstep) = Rs(current_state) + p_gamma * upcoming_value - Vs_latest(previous_action);
            if current_trial >= DAdep_start_trial
            	tmpTDeffective = TDs(k_tstep) * DAdep_factor;
            else
                tmpTDeffective = TDs(k_tstep);
            end
            Vs_latest(previous_action) = Vs_latest(previous_action) + p_alpha * tmpTDeffective;
        end
        Vs_latest = Vs_latest * (1 - decay_rate); % decay of learned values
        Vs_everystep = [Vs_everystep; Vs_latest];
        if current_action == 2*current_state % Go
            current_state = current_state + 1;
        end
        previous_action = current_action;
    else % when reaching the goal
        goalsteps = [goalsteps; k_tstep];
        TDs(k_tstep) = Rs(current_state) + p_gamma * 0 - Vs_latest(previous_action); % It is assumed that there is no "upcoming value" at the goal
        if current_trial >= DAdep_start_trial
            tmpTDeffective = TDs(k_tstep) * DAdep_factor;
        else
            tmpTDeffective = TDs(k_tstep);
        end
        Vs_latest(previous_action) = Vs_latest(previous_action) + p_alpha * tmpTDeffective;
        Vs_latest = Vs_latest * (1 - decay_rate); % decay of learned values
        Vs_everystep = [Vs_everystep; Vs_latest];
        Vs_whole = [Vs_whole; Vs_latest];
        if current_trial == num_trial
            TDs = TDs(1:k_tstep);
            States = States(1:k_tstep);
            break;
        end
        current_state = 1;
        previous_action = [];
    end
end
if length(goalsteps) < num_trial
    if isnan(prob_chooseNoGo)
        error('probability of choice became NaN');
    end
end

% output
Out.TDs = TDs;
Out.States = States;
Out.goalsteps = goalsteps;
Out.Vs_whole = Vs_whole;
Out.Vs_everystep = Vs_everystep;
save OutFig13CD Out

% plot
save_fig = 1;
critical_step = find(Vs_everystep(:,7)>1,1,'first');
if ~isempty(critical_step)
    if (critical_step > maxstep/2) && (critical_step < maxstep*(3/4))
        % figure for action values
        F = figure;
        A = axes;
        hold on;
        P = image(1+63*flipud(Vs_everystep(:,1:12)));
        C = colorbar; set(C,'YTick',1+63*[0:0.1:1],'YTickLabel',[0:0.1:1]);
        axis([0.5 12.5 0.5 maxstep+0.5]);
        set(A,'Box','off');
        %set(A,'PlotBoxAspectRatio',[1 1 1]);
        set(A,'FontName','Ariel','FontSize',20);
        set(A,'XTick',[1:12],'XTickLabel',[1:12]);
        set(A,'YTick',0.5+[0:100:maxstep-100],'YTickLabel',maxstep-[0:100:maxstep-100]);
        if save_fig
            print(F,'-depsc','Fig13C');
        end
        % figure for action values, without colorbar
        F = figure;
        A = axes;
        hold on;
        P = image(1+63*flipud(Vs_everystep(:,1:12)));
        %C = colorbar; set(C,'YTick',1+63*[0:0.1:1],'YTickLabel',[0:0.1:1]);
        axis([0.5 12.5 0.5 maxstep+0.5]);
        set(A,'Box','off');
        %set(A,'PlotBoxAspectRatio',[1 1 1]);
        set(A,'FontName','Ariel','FontSize',20);
        set(A,'XTick',[1:12],'XTickLabel',[1:12]);
        set(A,'YTick',0.5+[0:100:maxstep-100],'YTickLabel',maxstep-[0:100:maxstep-100]);
        if save_fig
            print(F,'-depsc','Fig13C_wo_colorbar');
        end
        % figure for state transitions
        F = figure;
        A = axes;
        hold on;
        P = plot(States(critical_step-80:critical_step+40),[120:-1:0],'k.-');
        axis([0.5 7.5 0 120]);
        set(A,'Box','off');
        %set(A,'PlotBoxAspectRatio',[1 1 1]);
        set(A,'FontName','Ariel','FontSize',20);
        set(A,'XTick',[1:7],'XTickLabel',[1:7]);
        tmpY = (critical_step+40) - floor((critical_step+40)/10)*10;
        set(A,'YTick',[tmpY:10:120],'YTickLabel',floor((critical_step+40)/10)*10:-10:ceil((critical_step-80)/10));
        if save_fig
            print(F,'-depsc','Fig13D');
        end
    end
end
