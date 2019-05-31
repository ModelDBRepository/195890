% mkfig_Fig11EF

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

% average ntspt (number of time steps per trial) and TD error (DA) in the case of [alpha0,beta0,gamma0] and decay_rate=0.01
num_sim = 20;
num_trial = 500;
rewsizevar.ntspt = zeros(5,num_sim);
rewsizevar.TD_per_tstep = zeros(5,num_sim);
rewsizevar.TD_per_trial = zeros(5,num_sim);
%
k_reward_size = 1;
load D_rewardsize0p50_alpha
for k3 = 1:num_sim
    rewsizevar.ntspt(k_reward_size,k3) = length(Dalpha.Out{6}{6}{k3}.TDs)/num_trial;
    rewsizevar.TD_per_tstep(k_reward_size,k3) = sum(Dalpha.Out{6}{6}{k3}.TDs)/length(Dalpha.Out{6}{6}{k3}.TDs);
    rewsizevar.TD_per_trial(k_reward_size,k3) = sum(Dalpha.Out{6}{6}{k3}.TDs)/num_trial;
end
clear Dalpha
%
k_reward_size = 2;
load D_rewardsize0p75_alpha
for k3 = 1:num_sim
    rewsizevar.ntspt(k_reward_size,k3) = length(Dalpha.Out{6}{6}{k3}.TDs)/num_trial;
    rewsizevar.TD_per_tstep(k_reward_size,k3) = sum(Dalpha.Out{6}{6}{k3}.TDs)/length(Dalpha.Out{6}{6}{k3}.TDs);
    rewsizevar.TD_per_trial(k_reward_size,k3) = sum(Dalpha.Out{6}{6}{k3}.TDs)/num_trial;
end
clear Dalpha
%
k_reward_size = 3;
load Dalpha
for k3 = 1:num_sim
    rewsizevar.ntspt(k_reward_size,k3) = length(Dalpha.Out{6}{6}{k3}.TDs)/num_trial;
    rewsizevar.TD_per_tstep(k_reward_size,k3) = sum(Dalpha.Out{6}{6}{k3}.TDs)/length(Dalpha.Out{6}{6}{k3}.TDs);
    rewsizevar.TD_per_trial(k_reward_size,k3) = sum(Dalpha.Out{6}{6}{k3}.TDs)/num_trial;
end
clear Dalpha
%
k_reward_size = 4;
load D_rewardsize1p25_alpha
for k3 = 1:num_sim
    rewsizevar.ntspt(k_reward_size,k3) = length(Dalpha.Out{6}{6}{k3}.TDs)/num_trial;
    rewsizevar.TD_per_tstep(k_reward_size,k3) = sum(Dalpha.Out{6}{6}{k3}.TDs)/length(Dalpha.Out{6}{6}{k3}.TDs);
    rewsizevar.TD_per_trial(k_reward_size,k3) = sum(Dalpha.Out{6}{6}{k3}.TDs)/num_trial;
end
clear Dalpha
%
k_reward_size = 5;
load D_rewardsize1p50_alpha
for k3 = 1:num_sim
    rewsizevar.ntspt(k_reward_size,k3) = length(Dalpha.Out{6}{6}{k3}.TDs)/num_trial;
    rewsizevar.TD_per_tstep(k_reward_size,k3) = sum(Dalpha.Out{6}{6}{k3}.TDs)/length(Dalpha.Out{6}{6}{k3}.TDs);
    rewsizevar.TD_per_trial(k_reward_size,k3) = sum(Dalpha.Out{6}{6}{k3}.TDs)/num_trial;
end
clear Dalpha

% average ntspt (number of time steps per trial) and TD error (DA) in the case of [alpha0,beta0,gamma0] and decay_rate=0
num_sim = 20;
num_trial = 500;
rewsizevar_d0.ntspt = zeros(5,num_sim);
rewsizevar_d0.TD_per_tstep = zeros(5,num_sim);
rewsizevar_d0.TD_per_trial = zeros(5,num_sim);
%
k_reward_size = 1;
load D_rewardsize0p50_alpha
for k3 = 1:num_sim
    rewsizevar_d0.ntspt(k_reward_size,k3) = length(Dalpha.Out{1}{6}{k3}.TDs)/num_trial;
    rewsizevar_d0.TD_per_tstep(k_reward_size,k3) = sum(Dalpha.Out{1}{6}{k3}.TDs)/length(Dalpha.Out{1}{6}{k3}.TDs);
    rewsizevar_d0.TD_per_trial(k_reward_size,k3) = sum(Dalpha.Out{1}{6}{k3}.TDs)/num_trial;
end
clear Dalpha
%
k_reward_size = 2;
load D_rewardsize0p75_alpha
for k3 = 1:num_sim
    rewsizevar_d0.ntspt(k_reward_size,k3) = length(Dalpha.Out{1}{6}{k3}.TDs)/num_trial;
    rewsizevar_d0.TD_per_tstep(k_reward_size,k3) = sum(Dalpha.Out{1}{6}{k3}.TDs)/length(Dalpha.Out{1}{6}{k3}.TDs);
    rewsizevar_d0.TD_per_trial(k_reward_size,k3) = sum(Dalpha.Out{1}{6}{k3}.TDs)/num_trial;
end
clear Dalpha
%
k_reward_size = 3;
load Dalpha
for k3 = 1:num_sim
    rewsizevar_d0.ntspt(k_reward_size,k3) = length(Dalpha.Out{1}{6}{k3}.TDs)/num_trial;
    rewsizevar_d0.TD_per_tstep(k_reward_size,k3) = sum(Dalpha.Out{1}{6}{k3}.TDs)/length(Dalpha.Out{1}{6}{k3}.TDs);
    rewsizevar_d0.TD_per_trial(k_reward_size,k3) = sum(Dalpha.Out{1}{6}{k3}.TDs)/num_trial;
end
clear Dalpha
%
k_reward_size = 4;
load D_rewardsize1p25_alpha
for k3 = 1:num_sim
    rewsizevar_d0.ntspt(k_reward_size,k3) = length(Dalpha.Out{1}{6}{k3}.TDs)/num_trial;
    rewsizevar_d0.TD_per_tstep(k_reward_size,k3) = sum(Dalpha.Out{1}{6}{k3}.TDs)/length(Dalpha.Out{1}{6}{k3}.TDs);
    rewsizevar_d0.TD_per_trial(k_reward_size,k3) = sum(Dalpha.Out{1}{6}{k3}.TDs)/num_trial;
end
clear Dalpha
%
k_reward_size = 5;
load D_rewardsize1p50_alpha
for k3 = 1:num_sim
    rewsizevar_d0.ntspt(k_reward_size,k3) = length(Dalpha.Out{1}{6}{k3}.TDs)/num_trial;
    rewsizevar_d0.TD_per_tstep(k_reward_size,k3) = sum(Dalpha.Out{1}{6}{k3}.TDs)/length(Dalpha.Out{1}{6}{k3}.TDs);
    rewsizevar_d0.TD_per_trial(k_reward_size,k3) = sum(Dalpha.Out{1}{6}{k3}.TDs)/num_trial;
end
clear Dalpha

% average ntspt (number of time steps per trial) and TD error (DA) in the case of [alpha0,beta0,gamma0] and decay_rate=0 or 0.01 for 1-100 and 401-500 trials
num_sim = 20;
num_trial = 500;
for j1 = [1 6] % 1:decay_rate=0, 6:decay_rate=0.01
    for j2 = 1:5 % 1:1-100 trials, 2:101-200 trials, ..., 5:401-500 trials
        rewsizevar_more.ntspt{j1}{j2} = zeros(5,num_sim);
        rewsizevar_more.TD_per_tstep{j1}{j2} = zeros(5,num_sim);
    end
end
%
k_reward_size = 1;
load D_rewardsize0p50_alpha
for j1 = [1 6] % 1:decay_rate=0, 6:decay_rate=0.01
    for j2 = 1:5 % 1:1-100 trials, 2:101-200 trials, ..., 5:401-500 trials
        for k3 = 1:num_sim
            tmp_goalsteps_100s = [0;Dalpha.Out{j1}{6}{k3}.goalsteps([100:100:500])];
            rewsizevar_more.ntspt{j1}{j2}(k_reward_size,k3) = (tmp_goalsteps_100s(j2+1) - tmp_goalsteps_100s(j2))/100;
            rewsizevar_more.TD_per_tstep{j1}{j2}(k_reward_size,k3) = ...
                sum(Dalpha.Out{j1}{6}{k3}.TDs([tmp_goalsteps_100s(j2)+1:tmp_goalsteps_100s(j2+1)])) / (tmp_goalsteps_100s(j2+1) - tmp_goalsteps_100s(j2));
        end
    end
end
clear Dalpha
%
k_reward_size = 2;
load D_rewardsize0p75_alpha
for j1 = [1 6] % 1:decay_rate=0, 6:decay_rate=0.01
    for j2 = 1:5 % 1:1-100 trials, 2:101-200 trials, ..., 5:401-500 trials
        for k3 = 1:num_sim
            tmp_goalsteps_100s = [0;Dalpha.Out{j1}{6}{k3}.goalsteps([100:100:500])];
            rewsizevar_more.ntspt{j1}{j2}(k_reward_size,k3) = (tmp_goalsteps_100s(j2+1) - tmp_goalsteps_100s(j2))/100;
            rewsizevar_more.TD_per_tstep{j1}{j2}(k_reward_size,k3) = ...
                sum(Dalpha.Out{j1}{6}{k3}.TDs([tmp_goalsteps_100s(j2)+1:tmp_goalsteps_100s(j2+1)])) / (tmp_goalsteps_100s(j2+1) - tmp_goalsteps_100s(j2));
        end
    end
end
clear Dalpha
%
k_reward_size = 3;
load Dalpha
for j1 = [1 6] % 1:decay_rate=0, 6:decay_rate=0.01
    for j2 = 1:5 % 1:1-100 trials, 2:101-200 trials, ..., 5:401-500 trials
        for k3 = 1:num_sim
            tmp_goalsteps_100s = [0;Dalpha.Out{j1}{6}{k3}.goalsteps([100:100:500])];
            rewsizevar_more.ntspt{j1}{j2}(k_reward_size,k3) = (tmp_goalsteps_100s(j2+1) - tmp_goalsteps_100s(j2))/100;
            rewsizevar_more.TD_per_tstep{j1}{j2}(k_reward_size,k3) = ...
                sum(Dalpha.Out{j1}{6}{k3}.TDs([tmp_goalsteps_100s(j2)+1:tmp_goalsteps_100s(j2+1)])) / (tmp_goalsteps_100s(j2+1) - tmp_goalsteps_100s(j2));
        end
    end
end
clear Dalpha
%
k_reward_size = 4;
load D_rewardsize1p25_alpha
for j1 = [1 6] % 1:decay_rate=0, 6:decay_rate=0.01
    for j2 = 1:5 % 1:1-100 trials, 2:101-200 trials, ..., 5:401-500 trials
        for k3 = 1:num_sim
            tmp_goalsteps_100s = [0;Dalpha.Out{j1}{6}{k3}.goalsteps([100:100:500])];
            rewsizevar_more.ntspt{j1}{j2}(k_reward_size,k3) = (tmp_goalsteps_100s(j2+1) - tmp_goalsteps_100s(j2))/100;
            rewsizevar_more.TD_per_tstep{j1}{j2}(k_reward_size,k3) = ...
                sum(Dalpha.Out{j1}{6}{k3}.TDs([tmp_goalsteps_100s(j2)+1:tmp_goalsteps_100s(j2+1)])) / (tmp_goalsteps_100s(j2+1) - tmp_goalsteps_100s(j2));
        end
    end
end
clear Dalpha
%
k_reward_size = 5;
load D_rewardsize1p50_alpha
for j1 = [1 6] % 1:decay_rate=0, 6:decay_rate=0.01
    for j2 = 1:5 % 1:1-100 trials, 2:101-200 trials, ..., 5:401-500 trials
        for k3 = 1:num_sim
            tmp_goalsteps_100s = [0;Dalpha.Out{j1}{6}{k3}.goalsteps([100:100:500])];
            rewsizevar_more.ntspt{j1}{j2}(k_reward_size,k3) = (tmp_goalsteps_100s(j2+1) - tmp_goalsteps_100s(j2))/100;
            rewsizevar_more.TD_per_tstep{j1}{j2}(k_reward_size,k3) = ...
                sum(Dalpha.Out{j1}{6}{k3}.TDs([tmp_goalsteps_100s(j2)+1:tmp_goalsteps_100s(j2+1)])) / (tmp_goalsteps_100s(j2+1) - tmp_goalsteps_100s(j2));
        end
    end
end
clear Dalpha


% plot
save_fig = 1;

% chance level
chance_step7 = 0; % initialization
for k = 0:Inf
    tmp = k * (nchoosek(6-1+k,k) / (2^k));
    if 1/(2^k) == 0
        break;
    end
    chance_step7 = chance_step7 + tmp;
end
chance_step7 = 7 + (chance_step7/(2^6));

% Fig. 11E
reward_size_set = [0.5:0.25:1.5];
F = figure;
A = axes;
hold on;
P = plot([0 2],[7 7],'k');
P = plot([0 2],chance_step7*[1 1],'k:');
P = plot(reward_size_set,mean(rewsizevar_more.ntspt{1}{1},2),'bo'); set(P,'MarkerSize',20);
P = plot(reward_size_set,mean(rewsizevar_more.ntspt{1}{5},2),'bx'); set(P,'MarkerSize',20);
P = plot(reward_size_set,mean(rewsizevar_more.ntspt{6}{1},2),'ko'); set(P,'MarkerSize',20);
P = plot(reward_size_set,mean(rewsizevar_more.ntspt{6}{5},2),'kx'); set(P,'MarkerSize',20);
P = errorbar(reward_size_set,mean(rewsizevar_d0.ntspt,2),std(rewsizevar_d0.ntspt,1,2)/sqrt(num_sim),'b-'); set(P,'MarkerSize',20);
P = plot(reward_size_set,mean(rewsizevar_d0.ntspt,2),'b--');
P = errorbar(reward_size_set,mean(rewsizevar.ntspt,2),std(rewsizevar.ntspt,1,2)/sqrt(num_sim),'k-'); set(P,'MarkerSize',20);
P = plot(reward_size_set,mean(rewsizevar.ntspt,2),'k--');
axis([0 2 6 14]);
set(A,'Box','off');
%set(A,'PlotBoxAspectRatio',[1 1 1]);
set(A,'FontName','Ariel','FontSize',20);
set(A,'XTick',[0:0.5:2],'XTickLabel',[0:0.5:2]);
set(A,'YTick',[6:1:14],'YTickLabel',[6:1:14]);
if save_fig
    print(F,'-depsc','Fig11E');
end

% Fig. 11F
reward_size_set = [0.5:0.25:1.5];
F = figure;
A = axes;
hold on;
P = plot(reward_size_set,mean(rewsizevar_more.TD_per_tstep{1}{1},2),'bo'); set(P,'MarkerSize',20);
P = plot(reward_size_set,mean(rewsizevar_more.TD_per_tstep{1}{5},2),'bx'); set(P,'MarkerSize',20);
P = plot(reward_size_set,mean(rewsizevar_more.TD_per_tstep{6}{1},2),'ko'); set(P,'MarkerSize',20);
P = plot(reward_size_set,mean(rewsizevar_more.TD_per_tstep{6}{5},2),'kx'); set(P,'MarkerSize',20);
P = errorbar(reward_size_set,mean(rewsizevar_d0.TD_per_tstep,2),std(rewsizevar_d0.TD_per_tstep,1,2),'b-'); set(P,'MarkerSize',20);
P = plot(reward_size_set,mean(rewsizevar_d0.TD_per_tstep,2),'b--');
P = errorbar(reward_size_set,mean(rewsizevar.TD_per_tstep,2),std(rewsizevar.TD_per_tstep,1,2),'k-'); set(P,'MarkerSize',20);
P = plot(reward_size_set,mean(rewsizevar.TD_per_tstep,2),'k--');
axis([0 2 0 0.3]);
set(A,'Box','off');
%set(A,'PlotBoxAspectRatio',[1 1 1]);
set(A,'FontName','Ariel','FontSize',20);
set(A,'XTick',[0:0.5:2],'XTickLabel',[0:0.5:2]);
set(A,'YTick',[0:0.05:0.3],'YTickLabel',[0:0.05:0.3]);
if save_fig
    print(F,'-depsc','Fig11F');
end
