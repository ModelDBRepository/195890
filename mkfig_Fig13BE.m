% mkfig_Fig13BE

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
load used_rand_twister_for_Fig13BE

% number of simulations and trials
num_sim = 20;
num_trial = 500;

% type of reinforcement learning algorithm
RLtype = 'Q';

% default parameter values
alpha0 = 0.5;
beta0 = 5;
gamma0 = 1;
decay_rate = 0.01;
DAdep_paras = [1,1001];
middlerew_set = [0:0.01:0.1]; % small reward at the middle

% sim varying middlerew, with decay
%rand('twister',sum(100*clock));
%Dsim.rand_twister = rand('twister');
Dsim.rand_twister = used_rand_twister_for_Fig13BE;
rand('twister',Dsim.rand_twister);
Dsim.ntspt = NaN(length(middlerew_set),num_sim); % number of time steps per trial
for k1 = 1:length(middlerew_set)
    for k2 = 1:num_sim
        fprintf('gamma %d-%d\n',k1,k2);
        Dsim.Out{k1}{k2} = RLdecayStayGo5(num_trial, RLtype, [alpha0,beta0,gamma0], decay_rate, DAdep_paras, middlerew_set(k1));
        if length(Dsim.Out{k1}{k2}.goalsteps) == num_trial
            Dsim.ntspt(k1,k2) = length(Dsim.Out{k1}{k2}.States)/num_trial;
        end
    end
end
Dsim.ntspt_mean = mean(Dsim.ntspt,3);
Dsim.ntspt_std = std(Dsim.ntspt,1,3);
D_decay = Dsim;
save D13BE_decay D_decay
clear Dsim


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

% Fig. 13B
F = figure;
A = axes;
hold on;
P = plot(middlerew_set,sum(~isnan(D_decay.ntspt),2)/num_sim,'k.-'); set(P,'MarkerSize',20);
axis([0 0.1 0 1]);
set(A,'Box','off');
%set(A,'PlotBoxAspectRatio',[1 1 1]);
set(A,'FontName','Ariel','FontSize',20);
set(A,'XTick',[0:0.01:0.1],'XTickLabel',[0:0.01:0.1]);
set(A,'YTick',[0:0.1:1],'YTickLabel',[0:10:100]);
if save_fig
    print(F,'-depsc','Fig13B');
end

% Fig. 13E
tmp_SE = NaN(length(middlerew_set),1);
for k = 1:length(middlerew_set)
    tmp_N = sum(~isnan(D_decay.ntspt(k,:)));
    if tmp_N > 0
        tmp_SE(k) = std2(D_decay.ntspt(k,:),1,2)/sqrt(tmp_N);
    end
end
F = figure;
A = axes;
hold on;
P = plot([0 0.1],7*[1 1],'k--');
P = plot([0 0.1],chance_step7*[1 1],'k:');
P = errorbar(middlerew_set,mean2(D_decay.ntspt,2),tmp_SE,'k-');
P = plot(middlerew_set,mean2(D_decay.ntspt,2),'k-');
axis([0 0.1 6 15]);
set(A,'Box','off');
%set(A,'PlotBoxAspectRatio',[1 1 1]);
set(A,'FontName','Ariel','FontSize',20);
set(A,'XTick',[0:0.01:0.1],'XTickLabel',[0:0.01:0.1]);
set(A,'YTick',[6:1:15],'YTickLabel',[6:1:15]);
if save_fig
    print(F,'-depsc','Fig13E');
end
