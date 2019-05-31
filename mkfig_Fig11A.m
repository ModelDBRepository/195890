% mkfig_Fig11A

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
load used_rand_twister_for_Fig11A

% reward size
reward_size = 0.5;

% number of simulations and trials
num_sim = 20;
num_trial = 500;

% type of reinforcement learning algorithm
RLtype = 'Q';

% default parameter values
num_state = 7;
alpha0 = 0.5;
beta0 = 5;
gamma0 = 1;
DAdep_paras = [1,1001];

% varying parameter values
decay_rate_set = [0:0.002:0.02];
alpha_set = [0:0.1:1];
beta_set = [0:1:10];
gamma_set = [0:0.1:1];

% sim varying alpha
%rand('twister',sum(100*clock));
%Dsim.rand_twister = rand('twister');
Dsim.rand_twister = used_rand_twister_for_Fig11A.alpha;
rand('twister',Dsim.rand_twister);
Dsim.ntspt = zeros(length(decay_rate_set),length(alpha_set),num_sim); % number of time steps per trial
for k1 = 1:length(decay_rate_set)
    for k2 = 1:length(alpha_set)
        Dsim.Vend{k1}{k2} = zeros(num_sim,num_state*2);
        Dsim.Vs_whole_ave{k1}{k2} = zeros(num_trial,num_state*2);
        Dsim.ntsptAllbin5{k1}{k2} = zeros(num_sim,100);
        for k3 = 1:num_sim
            fprintf('alpha %d-%d-%d\n',k1,k2,k3);
            Dsim.Out{k1}{k2}{k3} = RLdecayStayGo4(num_trial, RLtype, [alpha_set(k2),beta0,gamma0], reward_size, decay_rate_set(k1), DAdep_paras);
            Dsim.ntspt(k1,k2,k3) = length(Dsim.Out{k1}{k2}{k3}.States)/num_trial;
            Dsim.Vend{k1}{k2}(k3,:) = Dsim.Out{k1}{k2}{k3}.Vs_whole(end,:);
            Dsim.Vs_whole_ave{k1}{k2} = Dsim.Vs_whole_ave{k1}{k2} + Dsim.Out{k1}{k2}{k3}.Vs_whole/num_sim;
            Dsim.ntsptAllbin5{k1}{k2}(k3,:) = mean(reshape(diff([0;Dsim.Out{k1}{k2}{k3}.goalsteps]),5,100),1);
        end
        Dsim.ntsptAllbin5_mean{k1}{k2} = mean(Dsim.ntsptAllbin5{k1}{k2},1);
        Dsim.ntsptAllbin5_std{k1}{k2} = std(Dsim.ntsptAllbin5{k1}{k2},1,1);
    end
end
Dsim.ntspt_mean = mean(Dsim.ntspt,3);
Dsim.ntspt_std = std(Dsim.ntspt,1,3);
Dalpha = Dsim;
save D_rewardsize0p50_alpha Dalpha
clear Dsim

% sim varying beta
%rand('twister',sum(100*clock));
%Dsim.rand_twister = rand('twister');
Dsim.rand_twister = used_rand_twister_for_Fig11A.beta;
rand('twister',Dsim.rand_twister);
Dsim.ntspt = zeros(length(decay_rate_set),length(beta_set),num_sim); % number of time steps per trial
for k1 = 1:length(decay_rate_set)
    for k2 = 1:length(beta_set)
        Dsim.Vend{k1}{k2} = zeros(num_sim,num_state*2);
        Dsim.Vs_whole_ave{k1}{k2} = zeros(num_trial,num_state*2);
        Dsim.ntsptAllbin5{k1}{k2} = zeros(num_sim,100);
        for k3 = 1:num_sim
            fprintf('beta %d-%d-%d\n',k1,k2,k3);
            Dsim.Out{k1}{k2}{k3} = RLdecayStayGo4(num_trial, RLtype, [alpha0,beta_set(k2),gamma0], reward_size, decay_rate_set(k1), DAdep_paras);
            Dsim.ntspt(k1,k2,k3) = length(Dsim.Out{k1}{k2}{k3}.States)/num_trial;
            Dsim.Vend{k1}{k2}(k3,:) = Dsim.Out{k1}{k2}{k3}.Vs_whole(end,:);
            Dsim.Vs_whole_ave{k1}{k2} = Dsim.Vs_whole_ave{k1}{k2} + Dsim.Out{k1}{k2}{k3}.Vs_whole/num_sim;
            Dsim.ntsptAllbin5{k1}{k2}(k3,:) = mean(reshape(diff([0;Dsim.Out{k1}{k2}{k3}.goalsteps]),5,100),1);
        end
        Dsim.ntsptAllbin5_mean{k1}{k2} = mean(Dsim.ntsptAllbin5{k1}{k2},1);
        Dsim.ntsptAllbin5_std{k1}{k2} = std(Dsim.ntsptAllbin5{k1}{k2},1,1);
    end
end
Dsim.ntspt_mean = mean(Dsim.ntspt,3);
Dsim.ntspt_std = std(Dsim.ntspt,1,3);
Dbeta = Dsim;
save D_rewardsize0p50_beta Dbeta
clear Dsim

% sim varying gamma
%rand('twister',sum(100*clock));
%Dsim.rand_twister = rand('twister');
Dsim.rand_twister = used_rand_twister_for_Fig11A.gamma;
rand('twister',Dsim.rand_twister);
Dsim.ntspt = zeros(length(decay_rate_set),length(gamma_set),num_sim); % number of time steps per trial
for k1 = 1:length(decay_rate_set)
    for k2 = 1:length(gamma_set)
        Dsim.Vend{k1}{k2} = zeros(num_sim,num_state*2);
        Dsim.Vs_whole_ave{k1}{k2} = zeros(num_trial,num_state*2);
        Dsim.ntsptAllbin5{k1}{k2} = zeros(num_sim,100);
        for k3 = 1:num_sim
            fprintf('gamma %d-%d-%d\n',k1,k2,k3);
            Dsim.Out{k1}{k2}{k3} = RLdecayStayGo4(num_trial, RLtype, [alpha0,beta0,gamma_set(k2)], reward_size, decay_rate_set(k1), DAdep_paras);
            Dsim.ntspt(k1,k2,k3) = length(Dsim.Out{k1}{k2}{k3}.States)/num_trial;
            Dsim.Vend{k1}{k2}(k3,:) = Dsim.Out{k1}{k2}{k3}.Vs_whole(end,:);
            Dsim.Vs_whole_ave{k1}{k2} = Dsim.Vs_whole_ave{k1}{k2} + Dsim.Out{k1}{k2}{k3}.Vs_whole/num_sim;
            Dsim.ntsptAllbin5{k1}{k2}(k3,:) = mean(reshape(diff([0;Dsim.Out{k1}{k2}{k3}.goalsteps]),5,100),1);
        end
        Dsim.ntsptAllbin5_mean{k1}{k2} = mean(Dsim.ntsptAllbin5{k1}{k2},1);
        Dsim.ntsptAllbin5_std{k1}{k2} = std(Dsim.ntsptAllbin5{k1}{k2},1,1);
    end
end
Dsim.ntspt_mean = mean(Dsim.ntspt,3);
Dsim.ntspt_std = std(Dsim.ntspt,1,3);
Dgamma = Dsim;
save D_rewardsize0p50_gamma Dgamma
clear Dsim


% plot
save_fig = 1;
savename = 'Fig11A';

% left
fig11max = ceil(max([max(max(Dalpha.ntspt_mean)), max(max(Dbeta.ntspt_mean)), max(max(Dgamma.ntspt_mean))]));
fig11min = num_state;
F = figure;
A = axes;
hold on;
P = image(1+63*(Dalpha.ntspt_mean' - fig11min)/(fig11max - fig11min));
C = colorbar; set(C,'YTick',1+63*[0:1:fig11max-fig11min]/(fig11max-fig11min),'YTickLabel',[fig11min:1:fig11max]);
axis([0.5 length(decay_rate_set)+0.5 0.5 length(alpha_set)+0.5]);
set(A,'Box','off');
%set(A,'PlotBoxAspectRatio',[1 1 1]);
set(A,'FontName','Ariel','FontSize',20);
set(A,'XTick',[1:length(decay_rate_set)],'XTickLabel',decay_rate_set);
set(A,'YTick',[1:length(alpha_set)],'YTickLabel',alpha_set);
if save_fig
    print(F,'-depsc',[savename '-left']);
end

% middle
F = figure;
A = axes;
hold on;
P = image(1+63*(Dbeta.ntspt_mean' - fig11min)/(fig11max - fig11min));
C = colorbar; set(C,'YTick',1+63*[0:1:fig11max-fig11min]/(fig11max-fig11min),'YTickLabel',[fig11min:1:fig11max]);
axis([0.5 length(decay_rate_set)+0.5 0.5 length(beta_set)+0.5]);
set(A,'Box','off');
%set(A,'PlotBoxAspectRatio',[1 1 1]);
set(A,'FontName','Ariel','FontSize',20);
set(A,'XTick',[1:length(decay_rate_set)],'XTickLabel',decay_rate_set);
set(A,'YTick',[1:length(beta_set)],'YTickLabel',beta_set);
if save_fig
    print(F,'-depsc',[savename '-middle']);
end

% right
F = figure;
A = axes;
hold on;
P = image(1+63*(Dgamma.ntspt_mean' - fig11min)/(fig11max - fig11min));
C = colorbar; set(C,'YTick',1+63*[0:1:fig11max-fig11min]/(fig11max-fig11min),'YTickLabel',[fig11min:1:fig11max]);
axis([0.5 length(decay_rate_set)+0.5 0.5 length(gamma_set)+0.5]);
set(A,'Box','off');
%set(A,'PlotBoxAspectRatio',[1 1 1]);
set(A,'FontName','Ariel','FontSize',20);
set(A,'XTick',[1:length(decay_rate_set)],'XTickLabel',decay_rate_set);
set(A,'YTick',[1:length(gamma_set)],'YTickLabel',gamma_set);
if save_fig
    print(F,'-depsc',[savename '-right']);
end
