% mkfig_Fig9

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
load used_rand_twister_for_Fig9

% number of simulations and trials
num_sim = 20;
num_trial = 500;

% default parameter values
alpha0 = 0.5;
beta0 = 5;
gamma0 = 1;

% varying parameter values
dfactor_set = [1:-0.002:0.98]; % 1 - decay_rate
alpha_set = [0:0.1:1];
beta_set = [0:1:10];
gamma_set = [0:0.1:1];

% sim varying alpha
%rand('twister',sum(100*clock));
%Dsim.rand_twister = rand('twister');
Dsim.rand_twister = used_rand_twister_for_Fig9;
rand('twister',Dsim.rand_twister);
Dsim.ntspt = zeros(length(dfactor_set),length(alpha_set),num_sim); % number of time steps per trial
for k1 = 1:length(dfactor_set)
    for k2 = 1:length(alpha_set)
        Dsim.Vend{k1}{k2} = zeros(num_sim,7);
        Dsim.Vs_whole_ave{k1}{k2} = zeros(num_trial,7);
        Dsim.ntsptAllbin5{k1}{k2} = zeros(num_sim,100);
        for k3 = 1:num_sim
            fprintf('%d-%d-%d\n',k1,k2,k3);
            Dsim.Out{k1}{k2}{k3} = RLdecayStayGoStateValues(num_trial, [alpha_set(k2),beta0,gamma0], 1-dfactor_set(k1), [1, 1001]);
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
save Dalpha_StateValues Dalpha
clear Dsim
clear Dalpha


% plot
save_fig = 1;
load Dalpha_StateValues

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

% Fig. 9A
k_alpha = 6;
F = figure;
A = axes;
hold on;
P = plot([0 1-dfactor_set(end)],[7,7],'k:');
P = plot([0 1-dfactor_set(end)],chance_step7*[1,1],'k--');
P = errorbar(1-dfactor_set,Dalpha.ntspt_mean(:,k_alpha),Dalpha.ntspt_std(:,k_alpha)/sqrt(num_sim),'kx-');
set(P,'MarkerSize',20);
P = plot(1-dfactor_set,Dalpha.ntspt_mean(:,k_alpha),'k:');
YMax = chance_step7 + 1;
axis([0 1-dfactor_set(end) 6 YMax]);
set(A,'Box','off');
%set(A,'PlotBoxAspectRatio',[1 1 1]);
set(A,'FontName','Ariel','FontSize',20);
set(A,'XTick',[0:0.002:1-dfactor_set(end)],'XTickLabel',[0:0.002:1-dfactor_set(end)]);
set(A,'YTick',[6:1:YMax],'YTickLabel',[6:1:YMax]);
if save_fig
    print(F,'-depsc','Fig9A');
end

% Fig. 9B bottom panels
k_decay_set = [6 11];
k_alpha = 6;
tmp = 0;
for k_decay = k_decay_set
    tmp = max(tmp, max(mean(Dalpha.Vend{k_decay}{k_alpha},1) + std(Dalpha.Vend{k_decay}{k_alpha},1,1)/sqrt(num_sim)));
end
Fig9BbYmax = max(1, ceil(tmp*10)/10);
for k_decay = k_decay_set
    F = figure;
    A = axes;
    hold on;
    P = plot([0.5 7.5],[1 1],'k:');
    P = errorbar([1:7],mean(Dalpha.Vend{k_decay}{k_alpha},1),std(Dalpha.Vend{k_decay}{k_alpha},1,1)/sqrt(num_sim),'kx-');
    set(P,'MarkerSize',20);
    P = plot([1:7],mean(Dalpha.Vend{k_decay}{k_alpha},1),'k--');
    axis([0.5 7.5 0 Fig9BbYmax]);
    set(A,'Box','off');
    %set(A,'PlotBoxAspectRatio',[1 1 1]);
    set(A,'FontName','Ariel','FontSize',20);
    set(A,'XTick',[1:1:7],'XTickLabel',[1:1:7]);
    set(A,'YTick',[0:0.1:Fig9BbYmax],'YTickLabel',[0:0.1:Fig9BbYmax]);
    if save_fig
        print(F,'-depsc',['Fig9B-bottom_' num2str(k_decay)]);
    end
end

% Fig. 9B top panels
k_decay_set = [6 11];
k_alpha = 6;
for k_decay = k_decay_set
    F = figure;
    A = axes;
    hold on;
    P = image(1+63*flipud(Dalpha.Vs_whole_ave{k_decay}{k_alpha}(:,1:7)));
    C = colorbar; set(C,'YTick',1+63*[0:0.1:1],'YTickLabel',[0:0.1:1]);
    axis([0.5 7.5 0.5 500.5]);
    set(A,'Box','off');
    %set(A,'PlotBoxAspectRatio',[1 1 1]);
    set(A,'FontName','Ariel','FontSize',20);
    set(A,'XTick',[1:7],'XTickLabel',[1:7]);
    set(A,'YTick',0.5+[0:100:400],'YTickLabel',500-[0:100:400]);
    if save_fig
        print(F,'-depsc',['Fig9B-top_' num2str(k_decay)]);
    end
end
% without color bar
k_decay_set = [6 11];
k_alpha = 6;
for k_decay = k_decay_set
    F = figure;
    A = axes;
    hold on;
    P = image(1+63*flipud(Dalpha.Vs_whole_ave{k_decay}{k_alpha}(:,1:7)));
    %C = colorbar; set(C,'YTick',1+63*[0:0.1:1],'YTickLabel',[0:0.1:1]);
    axis([0.5 7.5 0.5 500.5]);
    set(A,'Box','off');
    %set(A,'PlotBoxAspectRatio',[1 1 1]);
    set(A,'FontName','Ariel','FontSize',20);
    set(A,'XTick',[1:7],'XTickLabel',[1:7]);
    set(A,'YTick',0.5+[0:100:400],'YTickLabel',500-[0:100:400]);
    if save_fig
        print(F,'-depsc',['Fig9B-top_wocolorbar_' num2str(k_decay)]);
    end
end
