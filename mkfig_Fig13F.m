% mkfig_Fig13F

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
load used_rand_twister_for_Fig13F

% number of simulations and trials
num_sim = 20;
num_trial = 500;

% type of reinforcement learning algorithm
RLtype = 'Q';

% default parameter values
alpha0 = 0.5;
beta0 = 5;
gamma0 = 1;
DAdep_paras = [1,1001];
middlerew = 0.1; % small reward at the middle

% varying parameter values
decay_rate_set = [0:0.002:0.02];
alpha_set = [0:0.1:1];
beta_set = [0:1:10];
gamma_set = [0:0.1:1];

% sim varying gamma
%rand('twister',sum(100*clock));
%Dsim.rand_twister = rand('twister');
Dsim.rand_twister = used_rand_twister_for_Fig13F;
rand('twister',Dsim.rand_twister);
Dsim.ntspt = NaN(length(decay_rate_set),length(gamma_set),num_sim); % number of time steps per trial
for k1 = 1:length(decay_rate_set)
    for k2 = 1:length(gamma_set)
        Dsim.ntsptAllbin5{k1}{k2} = NaN(num_sim,100);
        for k3 = 1:num_sim
            fprintf('gamma %d-%d-%d\n',k1,k2,k3);
            Dsim.Out{k1}{k2}{k3} = RLdecayStayGo5(num_trial, RLtype, [alpha0,beta0,gamma_set(k2)], decay_rate_set(k1), DAdep_paras, middlerew);
            if length(Dsim.Out{k1}{k2}{k3}.goalsteps) == num_trial
                Dsim.ntspt(k1,k2,k3) = length(Dsim.Out{k1}{k2}{k3}.States)/num_trial;
                Dsim.ntsptAllbin5{k1}{k2}(k3,:) = mean(reshape(diff([0;Dsim.Out{k1}{k2}{k3}.goalsteps]),5,100),1);
            end
        end
        Dsim.ntsptAllbin5_mean{k1}{k2} = mean(Dsim.ntsptAllbin5{k1}{k2},1);
        Dsim.ntsptAllbin5_std{k1}{k2} = std(Dsim.ntsptAllbin5{k1}{k2},1,1);
    end
end
Dsim.ntspt_mean = mean(Dsim.ntspt,3);
Dsim.ntspt_std = std(Dsim.ntspt,1,3);
Dgamma = Dsim;
save D13Fgamma Dgamma
clear Dsim


% check if 500 trials were completed in all the simulations
for k1 = 1:length(decay_rate_set)
    for k2 = 1:length(gamma_set)-1
        for k3 = 1:num_sim
            if length(Dgamma.Out{k1}{k2}{k3}.goalsteps) ~= num_trial
                error('500 trials were not completed');
            end
        end
    end
end


% plot
save_fig = 1;

% number of simulation runs in which subject settled at S4 before completing 500 trials
num_sim_not_completed_500 = zeros(length(decay_rate_set),length(gamma_set));
for k1 = 1:length(decay_rate_set)
    for k2 = 1:length(gamma_set)
        for k3 = 1:num_sim
            if length(Dgamma.Out{k1}{k2}{k3}.goalsteps) < num_trial
                num_sim_not_completed_500(k1,k2) = num_sim_not_completed_500(k1,k2) + 1;
            end
        end
    end
end

% Fig. 13F
fig13Fmax = ceil(max(max(Dgamma.ntspt_mean)));
fig13Fmin = 7;
F = figure;
A = axes;
hold on;
P = image(1+63*(Dgamma.ntspt_mean' - fig13Fmin)/(fig13Fmax - fig13Fmin));
C = colorbar; set(C,'YTick',1+63*[0:1:fig13Fmax-fig13Fmin]/(fig13Fmax-fig13Fmin),'YTickLabel',[fig13Fmin:1:fig13Fmax]);
axis([0.5 length(decay_rate_set)+0.5 0.5 length(gamma_set)+0.5]);
set(A,'Box','off');
%set(A,'PlotBoxAspectRatio',[1 1 1]);
set(A,'FontName','Ariel','FontSize',20);
set(A,'XTick',[1:length(decay_rate_set)],'XTickLabel',decay_rate_set);
set(A,'YTick',[1:length(gamma_set)],'YTickLabel',gamma_set);
if save_fig
    print(F,'-depsc','Fig13F');
end
