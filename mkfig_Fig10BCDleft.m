% mkfig_Fig10BCDleft

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
load used_rand_twister_for_Fig10BCDleft

% number of simulations and trials
num_sim = 20;
num_trial = 500;

% type of reinforcement learning algorithm
RLtype = 'S';

% default parameter values
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
Dsim.rand_twister = used_rand_twister_for_Fig10BCDleft.alpha;
rand('twister',Dsim.rand_twister);
Dsim.ntspt = zeros(length(decay_rate_set),length(alpha_set),num_sim); % number of time steps per trial
for k1 = 1:length(decay_rate_set)
    for k2 = 1:length(alpha_set)
        Dsim.Vend{k1}{k2} = zeros(num_sim,7*2);
        Dsim.Vs_whole_ave{k1}{k2} = zeros(num_trial,7*2);
        Dsim.ntsptAllbin5{k1}{k2} = zeros(num_sim,100);
        for k3 = 1:num_sim
            fprintf('alpha %d-%d-%d\n',k1,k2,k3);
            Dsim.Out{k1}{k2}{k3} = RLdecayStayGo3(num_trial, RLtype, [alpha_set(k2),beta0,gamma0], decay_rate_set(k1), DAdep_paras);
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
save DSARSAalpha Dalpha
clear Dsim
clear Dalpha

% sim varying beta
%rand('twister',sum(100*clock));
%Dsim.rand_twister = rand('twister');
Dsim.rand_twister = used_rand_twister_for_Fig10BCDleft.beta;
rand('twister',Dsim.rand_twister);
Dsim.ntspt = zeros(length(decay_rate_set),length(beta_set),num_sim); % number of time steps per trial
for k1 = 1:length(decay_rate_set)
    for k2 = 1:length(beta_set)
        Dsim.Vend{k1}{k2} = zeros(num_sim,7*2);
        Dsim.Vs_whole_ave{k1}{k2} = zeros(num_trial,7*2);
        Dsim.ntsptAllbin5{k1}{k2} = zeros(num_sim,100);
        for k3 = 1:num_sim
            fprintf('beta %d-%d-%d\n',k1,k2,k3);
            Dsim.Out{k1}{k2}{k3} = RLdecayStayGo3(num_trial, RLtype, [alpha0,beta_set(k2),gamma0], decay_rate_set(k1), DAdep_paras);
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
save DSARSAbeta Dbeta
clear Dsim
clear Dbeta

% sim varying gamma
%rand('twister',sum(100*clock));
%Dsim.rand_twister = rand('twister');
Dsim.rand_twister = used_rand_twister_for_Fig10BCDleft.gamma;
rand('twister',Dsim.rand_twister);
Dsim.ntspt = zeros(length(decay_rate_set),length(gamma_set),num_sim); % number of time steps per trial
for k1 = 1:length(decay_rate_set)
    for k2 = 1:length(gamma_set)
        Dsim.Vend{k1}{k2} = zeros(num_sim,7*2);
        Dsim.Vs_whole_ave{k1}{k2} = zeros(num_trial,7*2);
        Dsim.ntsptAllbin5{k1}{k2} = zeros(num_sim,100);
        for k3 = 1:num_sim
            fprintf('gamma %d-%d-%d\n',k1,k2,k3);
            Dsim.Out{k1}{k2}{k3} = RLdecayStayGo3(num_trial, RLtype, [alpha0,beta0,gamma_set(k2)], decay_rate_set(k1), DAdep_paras);
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
save DSARSAgamma Dgamma
clear Dsim
clear Dgamma


% plot
save_fig = 1;
load DSARSAalpha
load DSARSAbeta
load DSARSAgamma

% Fig. 10B-left
fig6max = ceil(max([max(max(Dalpha.ntspt_mean)), max(max(Dbeta.ntspt_mean)), max(max(Dgamma.ntspt_mean))]));
fig6min = 7;
F = figure;
A = axes;
hold on;
P = image(1+63*(Dalpha.ntspt_mean' - fig6min)/(fig6max - fig6min));
C = colorbar; set(C,'YTick',1+63*[0:1:fig6max-fig6min]/(fig6max-fig6min),'YTickLabel',[fig6min:1:fig6max]);
axis([0.5 length(decay_rate_set)+0.5 0.5 length(alpha_set)+0.5]);
set(A,'Box','off');
%set(A,'PlotBoxAspectRatio',[1 1 1]);
set(A,'FontName','Ariel','FontSize',20);
set(A,'XTick',[1:length(decay_rate_set)],'XTickLabel',decay_rate_set);
set(A,'YTick',[1:length(alpha_set)],'YTickLabel',alpha_set);
if save_fig
    print(F,'-depsc','Fig10B-left');
end

% Fig. 10B-middle
F = figure;
A = axes;
hold on;
P = image(1+63*(Dbeta.ntspt_mean' - fig6min)/(fig6max - fig6min));
C = colorbar; set(C,'YTick',1+63*[0:1:fig6max-fig6min]/(fig6max-fig6min),'YTickLabel',[fig6min:1:fig6max]);
axis([0.5 length(decay_rate_set)+0.5 0.5 length(beta_set)+0.5]);
set(A,'Box','off');
%set(A,'PlotBoxAspectRatio',[1 1 1]);
set(A,'FontName','Ariel','FontSize',20);
set(A,'XTick',[1:length(decay_rate_set)],'XTickLabel',decay_rate_set);
set(A,'YTick',[1:length(beta_set)],'YTickLabel',beta_set);
if save_fig
    print(F,'-depsc','Fig10B-middle');
end

% Fig. 10B-right
F = figure;
A = axes;
hold on;
P = image(1+63*(Dgamma.ntspt_mean' - fig6min)/(fig6max - fig6min));
C = colorbar; set(C,'YTick',1+63*[0:1:fig6max-fig6min]/(fig6max-fig6min),'YTickLabel',[fig6min:1:fig6max]);
axis([0.5 length(decay_rate_set)+0.5 0.5 length(gamma_set)+0.5]);
set(A,'Box','off');
%set(A,'PlotBoxAspectRatio',[1 1 1]);
set(A,'FontName','Ariel','FontSize',20);
set(A,'XTick',[1:length(decay_rate_set)],'XTickLabel',decay_rate_set);
set(A,'YTick',[1:length(gamma_set)],'YTickLabel',gamma_set);
if save_fig
    print(F,'-depsc','Fig10B-right');
end

% Fig. 10C
k_decay_set = [6];
k_alpha = 6;
tmp = 0;
for k_decay = k_decay_set
    tmp = max(tmp, max(mean(Dalpha.Vend{k_decay}{k_alpha}(:,2:2:12),1) + std(Dalpha.Vend{k_decay}{k_alpha}(:,2:2:12),1,1)/sqrt(num_sim)));
end
tmp = max(tmp, max(mean(Dgamma.Vend{1}{8}(:,2:2:12),1) + std(Dgamma.Vend{1}{8}(:,2:2:12),1,1)/sqrt(num_sim)));
Fig10CYmax = max(1, ceil(tmp*10)/10);
for k_decay = k_decay_set
    F = figure;
    A = axes;
    hold on;
    P = plot([0.5 12.5],[1 1],'k:');
    P = errorbar([2:2:12],mean(Dalpha.Vend{k_decay}{k_alpha}(:,2:2:12),1),std(Dalpha.Vend{k_decay}{k_alpha}(:,2:2:12),1,1)/sqrt(num_sim),'kx-');
    set(P,'MarkerSize',20);
    P = plot([2:2:12],mean(Dalpha.Vend{k_decay}{k_alpha}(:,2:2:12),1),'k--');
    P = errorbar([1:2:11],mean(Dalpha.Vend{k_decay}{k_alpha}(:,1:2:11),1),std(Dalpha.Vend{k_decay}{k_alpha}(:,1:2:11),1,1)/sqrt(num_sim),'rx-');
    set(P,'MarkerSize',20,'Color',0.5*[1 1 1]);
    P = plot([1:2:11],mean(Dalpha.Vend{k_decay}{k_alpha}(:,1:2:11),1),'r--');
    set(P,'Color',0.5*[1 1 1]);
    axis([0.5 12.5 0 Fig10CYmax]);
    set(A,'Box','off');
    %set(A,'PlotBoxAspectRatio',[1 1 1]);
    set(A,'FontName','Ariel','FontSize',20);
    set(A,'XTick',[1:1:12],'XTickLabel',[1:1:12]);
    set(A,'YTick',[0:0.1:Fig10CYmax],'YTickLabel',[0:0.1:Fig10CYmax]);
    if save_fig
        print(F,'-depsc',['Fig10C_' num2str(k_decay)]);
    end
end

% Fig. 10D-left
k_decay = 6;
k_alpha = 6;
TDs_Stay_Go = zeros(2,num_sim); % 1st row: Stay, 2nd row: Go
for k3 = 1:num_sim
    tmpTDs_Stay = [];
    tmpTDs_Go = [];
    for k_tstep = 1:length(Dalpha.Out{k_decay}{k_alpha}{k3}.TDs)-1
        if Dalpha.Out{k_decay}{k_alpha}{k3}.States(k_tstep) == Dalpha.Out{k_decay}{k_alpha}{k3}.States(k_tstep+1) % Stay
            tmpTDs_Stay = [tmpTDs_Stay, Dalpha.Out{k_decay}{k_alpha}{k3}.TDs(k_tstep)];
        elseif Dalpha.Out{k_decay}{k_alpha}{k3}.States(k_tstep) + 1 == Dalpha.Out{k_decay}{k_alpha}{k3}.States(k_tstep+1) % Go
            tmpTDs_Go = [tmpTDs_Go, Dalpha.Out{k_decay}{k_alpha}{k3}.TDs(k_tstep)];
        end
    end
    TDs_Stay_Go(1,k3) = mean(tmpTDs_Stay);
    TDs_Stay_Go(2,k3) = mean(tmpTDs_Go);
end
F = figure;
A = axes;
hold on;
P = plot([0.5 2.5],[0 0],'k--');
P = bar([1 2],mean(TDs_Stay_Go,2));
shading flat; set(P,'BarWidth',0.75,'FaceColor',0.5*[1 1 1]);
P = errorbar([1 2],mean(TDs_Stay_Go,2),std(TDs_Stay_Go,1,2)/sqrt(num_sim),'kx');
axis([0.5 2.5 -0.05 0.2]);
set(A,'Box','off');
set(A,'PlotBoxAspectRatio',[1 2 1]);
set(A,'FontName','Ariel','FontSize',20);
set(A,'XTick',[1 2],'XTickLabel',[]);
set(A,'YTick',[-0.05:0.05:0.2],'YTickLabel',[-0.05:0.05:0.2]);
if save_fig
    print(F,'-depsc','Fig10D-left');
end
