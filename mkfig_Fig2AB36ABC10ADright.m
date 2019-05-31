% mkfig_Fig2AB36ABC10ADright

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
load used_rand_twister_for_Fig2AB36ABC10ADr

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

% varying parameter values
decay_rate_set = [0:0.002:0.02];
alpha_set = [0:0.1:1];
beta_set = [0:1:10];
gamma_set = [0:0.1:1];

% sim varying alpha
%rand('twister',sum(100*clock));
%Dsim.rand_twister = rand('twister');
Dsim.rand_twister = used_rand_twister_for_Fig2AB36ABC10ADr.alpha;
rand('twister',Dsim.rand_twister);
Dsim.ntspt = zeros(length(decay_rate_set),length(alpha_set),num_sim); % number of time steps per trial
for k1 = 1:length(decay_rate_set)
    for k2 = 1:length(alpha_set)
        Dsim.Vend{k1}{k2} = zeros(num_sim,7*2);
        Dsim.Vs_whole_ave{k1}{k2} = zeros(num_trial,7*2);
        Dsim.ntsptAllbin5{k1}{k2} = zeros(num_sim,100);
        for k3 = 1:num_sim
            fprintf('alpha %d-%d-%d\n',k1,k2,k3);
            Dsim.Out{k1}{k2}{k3} = RLdecayStayGo1(num_trial, RLtype, [alpha_set(k2),beta0,gamma0], decay_rate_set(k1), DAdep_paras);
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
save Dalpha Dalpha
clear Dsim
clear Dalpha

% sim varying beta
%rand('twister',sum(100*clock));
%Dsim.rand_twister = rand('twister');
Dsim.rand_twister = used_rand_twister_for_Fig2AB36ABC10ADr.beta;
rand('twister',Dsim.rand_twister);
Dsim.ntspt = zeros(length(decay_rate_set),length(beta_set),num_sim); % number of time steps per trial
for k1 = 1:length(decay_rate_set)
    for k2 = 1:length(beta_set)
        Dsim.Vend{k1}{k2} = zeros(num_sim,7*2);
        Dsim.Vs_whole_ave{k1}{k2} = zeros(num_trial,7*2);
        Dsim.ntsptAllbin5{k1}{k2} = zeros(num_sim,100);
        for k3 = 1:num_sim
            fprintf('beta %d-%d-%d\n',k1,k2,k3);
            Dsim.Out{k1}{k2}{k3} = RLdecayStayGo1(num_trial, RLtype, [alpha0,beta_set(k2),gamma0], decay_rate_set(k1), DAdep_paras);
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
save Dbeta Dbeta
clear Dsim
clear Dbeta

% sim varying gamma
%rand('twister',sum(100*clock));
%Dsim.rand_twister = rand('twister');
Dsim.rand_twister = used_rand_twister_for_Fig2AB36ABC10ADr.gamma;
rand('twister',Dsim.rand_twister);
Dsim.ntspt = zeros(length(decay_rate_set),length(gamma_set),num_sim); % number of time steps per trial
for k1 = 1:length(decay_rate_set)
    for k2 = 1:length(gamma_set)
        Dsim.Vend{k1}{k2} = zeros(num_sim,7*2);
        Dsim.Vs_whole_ave{k1}{k2} = zeros(num_trial,7*2);
        Dsim.ntsptAllbin5{k1}{k2} = zeros(num_sim,100);
        for k3 = 1:num_sim
            fprintf('gamma %d-%d-%d\n',k1,k2,k3);
            Dsim.Out{k1}{k2}{k3} = RLdecayStayGo1(num_trial, RLtype, [alpha0,beta0,gamma_set(k2)], decay_rate_set(k1), DAdep_paras);
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
save Dgamma Dgamma
clear Dsim
clear Dgamma


% plot
save_fig = 1;
load Dalpha
load Dbeta
load Dgamma

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

% Fig. 2A
k_alpha = 6;
F = figure;
A = axes;
hold on;
P = plot(decay_rate_set([1,end]),[7,7],'k:');
P = plot(decay_rate_set([1,end]),chance_step7*[1,1],'k--');
P = errorbar(decay_rate_set,Dalpha.ntspt_mean(:,k_alpha),Dalpha.ntspt_std(:,k_alpha)/sqrt(num_sim),'kx-');
set(P,'MarkerSize',20);
P = plot(decay_rate_set,Dalpha.ntspt_mean(:,k_alpha),'k:');
YMax = chance_step7 + 1;
axis([0 decay_rate_set(end) 6 YMax]);
set(A,'Box','off');
%set(A,'PlotBoxAspectRatio',[1 1 1]);
set(A,'FontName','Ariel','FontSize',20);
set(A,'XTick',decay_rate_set,'XTickLabel',decay_rate_set);
set(A,'YTick',[6:1:YMax],'YTickLabel',[6:1:YMax]);
if save_fig
    print(F,'-depsc','Fig2A');
end

% Fig. 2B
k_decay_set = [1 6 11];
k_alpha = 6;
tmp = 0;
for k_decay = k_decay_set
    tmp = max(tmp, max(Dalpha.ntsptAllbin5_mean{k_decay}{k_alpha}+Dalpha.ntsptAllbin5_std{k_decay}{k_alpha}/sqrt(num_sim)));
end
YMax = ceil(tmp);
for k_decay = k_decay_set
    F = figure;
    A = axes;
    hold on;
    P = plot([0 100],[7,7],'k:');
    P = plot([0 100],chance_step7*[1,1],'k:');
    P = plot([1:100]-0.5,Dalpha.ntsptAllbin5_mean{k_decay}{k_alpha},'k');
    P = plot([1:100]-0.5,Dalpha.ntsptAllbin5_mean{k_decay}{k_alpha}+Dalpha.ntsptAllbin5_std{k_decay}{k_alpha}/sqrt(num_sim),'k--');
    P = plot([1:100]-0.5,Dalpha.ntsptAllbin5_mean{k_decay}{k_alpha}-Dalpha.ntsptAllbin5_std{k_decay}{k_alpha}/sqrt(num_sim),'k--');
    axis([0 100 6 YMax]);
    set(A,'Box','off');
    %set(A,'PlotBoxAspectRatio',[1 1 1]);
    set(A,'FontName','Ariel','FontSize',20);
    set(A,'XTick',[0:10:100],'XTickLabel',[0:50:500]);
    set(A,'YTick',[6:1:YMax],'YTickLabel',[6:1:YMax]);
    if save_fig
        print(F,'-depsc',['Fig2B_' num2str(k_decay)]);
    end
end

% Fig. 3A
k_decay_set = [1 6 11];
k_alpha = 6;
tmp = 0;
for k_decay = k_decay_set
    tmp = max(tmp, max(mean(Dalpha.Vend{k_decay}{k_alpha}(:,2:2:12),1) + std(Dalpha.Vend{k_decay}{k_alpha}(:,2:2:12),1,1)/sqrt(num_sim)));
end
tmp = max(tmp, max(mean(Dgamma.Vend{1}{8}(:,2:2:12),1) + std(Dgamma.Vend{1}{8}(:,2:2:12),1,1)/sqrt(num_sim)));
Fig3AYmax = max(1, ceil(tmp*10)/10);
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
    axis([0.5 12.5 0 Fig3AYmax]);
    set(A,'Box','off');
    %set(A,'PlotBoxAspectRatio',[1 1 1]);
    set(A,'FontName','Ariel','FontSize',20);
    set(A,'XTick',[1:1:12],'XTickLabel',[1:1:12]);
    set(A,'YTick',[0:0.1:Fig3AYmax],'YTickLabel',[0:0.1:Fig3AYmax]);
    if save_fig
        print(F,'-depsc',['Fig3A_' num2str(k_decay)]);
    end
end

% Fig. 3B
k_decay_set = [1 6 11];
k_alpha = 6;
for k_decay = k_decay_set
    F = figure;
    A = axes;
    hold on;
    P = image(1+63*flipud(Dalpha.Vs_whole_ave{k_decay}{k_alpha}(:,1:12)));
    C = colorbar; set(C,'YTick',1+63*[0:0.1:1],'YTickLabel',[0:0.1:1]);
    axis([0.5 12.5 0.5 500.5]);
    set(A,'Box','off');
    %set(A,'PlotBoxAspectRatio',[1 1 1]);
    set(A,'FontName','Ariel','FontSize',20);
    set(A,'XTick',[1:12],'XTickLabel',[1:12]);
    set(A,'YTick',0.5+[0:100:400],'YTickLabel',500-[0:100:400]);
    if save_fig
        print(F,'-depsc',['Fig3B_' num2str(k_decay)]);
    end
end
% without colorbar
k_decay_set = [1 6 11];
k_alpha = 6;
for k_decay = k_decay_set
    F = figure;
    A = axes;
    hold on;
    P = image(1+63*flipud(Dalpha.Vs_whole_ave{k_decay}{k_alpha}(:,1:12)));
    %C = colorbar; set(C,'YTick',1+63*[0:0.1:1],'YTickLabel',[0:0.1:1]);
    axis([0.5 12.5 0.5 500.5]);
    set(A,'Box','off');
    %set(A,'PlotBoxAspectRatio',[1 1 1]);
    set(A,'FontName','Ariel','FontSize',20);
    set(A,'XTick',[1:12],'XTickLabel',[1:12]);
    set(A,'YTick',0.5+[0:100:400],'YTickLabel',500-[0:100:400]);
    if save_fig
        print(F,'-depsc',['Fig3B_wocolorbar_' num2str(k_decay)]);
    end
end

% Fig. 3C
k_decay_set = [1 6 11];
k_alpha = 6;
k_sim = 11;
k_trial_start = 401;
for k_decay = k_decay_set
    F = figure;
    A = axes;
    hold on;
    for k_trial = k_trial_start:k_trial_start+9
        P = plot((Dalpha.Out{k_decay}{k_alpha}{k_sim}.goalsteps(k_trial-1)+0.5)*[1 1],[0 1],'k:');
        tmp_tsteps = [Dalpha.Out{k_decay}{k_alpha}{k_sim}.goalsteps(k_trial-1)+1:Dalpha.Out{k_decay}{k_alpha}{k_sim}.goalsteps(k_trial)];
        P = plot(tmp_tsteps,Dalpha.Out{k_decay}{k_alpha}{k_sim}.TDs(tmp_tsteps),'k');
    end
    axis([Dalpha.Out{k_decay}{k_alpha}{k_sim}.goalsteps(k_trial_start-1)+[0 50] 0 1]);
    set(A,'Box','off');
    set(A,'PlotBoxAspectRatio',[2 1 1]);
    set(A,'FontName','Ariel','FontSize',20);
    tmp = ceil(Dalpha.Out{k_decay}{k_alpha}{k_sim}.goalsteps(k_trial_start-1)/10)*10;
    set(A,'XTick',tmp+[0:10:50],'XTickLabel',tmp+[0:10:50]);
    set(A,'YTick',[0:0.1:1],'YTickLabel',[0:0.1:1]);
    if save_fig
        print(F,'-depsc',['Fig3C_' num2str(k_decay)]);
    end
end

% Fig. 6A
k_decay = 1;
k_gamma = 9;
F = figure;
A = axes;
hold on;
P = plot([0.5 12.5],[1 1],'k:');
P = errorbar([2:2:12],mean(Dgamma.Vend{k_decay}{k_gamma}(:,2:2:12),1),std(Dgamma.Vend{k_decay}{k_gamma}(:,2:2:12),1,1)/sqrt(num_sim),'kx-');
set(P,'MarkerSize',20);
P = plot([2:2:12],mean(Dgamma.Vend{k_decay}{k_gamma}(:,2:2:12),1),'k--');
P = errorbar([1:2:11],mean(Dgamma.Vend{k_decay}{k_gamma}(:,1:2:11),1),std(Dgamma.Vend{k_decay}{k_gamma}(:,1:2:11),1,1)/sqrt(num_sim),'rx-');
set(P,'MarkerSize',20,'Color',0.5*[1 1 1]);
P = plot([1:2:11],mean(Dgamma.Vend{k_decay}{k_gamma}(:,1:2:11),1),'r--');
set(P,'Color',0.5*[1 1 1]);
axis([0.5 12.5 0 Fig3AYmax]);
set(A,'Box','off');
%set(A,'PlotBoxAspectRatio',[1 1 1]);
set(A,'FontName','Ariel','FontSize',20);
set(A,'XTick',[1:1:12],'XTickLabel',[1:1:12]);
set(A,'YTick',[0:0.1:Fig3AYmax],'YTickLabel',[0:0.1:Fig3AYmax]);
if save_fig
    print(F,'-depsc','Fig6A');
end

% Fig. 6B
k_decay = 1;
k_gamma = 9;
F = figure;
A = axes;
hold on;
P = image(1+63*flipud(Dgamma.Vs_whole_ave{k_decay}{k_gamma}(:,1:12)));
C = colorbar; set(C,'YTick',1+63*[0:0.1:1],'YTickLabel',[0:0.1:1]);
axis([0.5 12.5 0.5 500.5]);
set(A,'Box','off');
%set(A,'PlotBoxAspectRatio',[1 1 1]);
set(A,'FontName','Ariel','FontSize',20);
set(A,'XTick',[1:12],'XTickLabel',[1:12]);
set(A,'YTick',0.5+[0:100:400],'YTickLabel',500-[0:100:400]);
if save_fig
    print(F,'-depsc','Fig6B');
end
% without colorbar
k_decay = 1;
k_gamma = 9;
F = figure;
A = axes;
hold on;
P = image(1+63*flipud(Dgamma.Vs_whole_ave{k_decay}{k_gamma}(:,1:12)));
%C = colorbar; set(C,'YTick',1+63*[0:0.1:1],'YTickLabel',[0:0.1:1]);
axis([0.5 12.5 0.5 500.5]);
set(A,'Box','off');
%set(A,'PlotBoxAspectRatio',[1 1 1]);
set(A,'FontName','Ariel','FontSize',20);
set(A,'XTick',[1:12],'XTickLabel',[1:12]);
set(A,'YTick',0.5+[0:100:400],'YTickLabel',500-[0:100:400]);
if save_fig
    print(F,'-depsc','Fig6B_wocolorbar');
end

% Fig. 6C
k_decay = 1;
k_gamma = 9;
k_sim = 1;
k_trial_start = 401;
F = figure;
A = axes;
hold on;
for k_trial = k_trial_start:k_trial_start+9
    P = plot((Dgamma.Out{k_decay}{k_gamma}{k_sim}.goalsteps(k_trial-1)+0.5)*[1 1],[0 1],'k:');
    tmp_tsteps = [Dgamma.Out{k_decay}{k_gamma}{k_sim}.goalsteps(k_trial-1)+1:Dgamma.Out{k_decay}{k_gamma}{k_sim}.goalsteps(k_trial)];
    P = plot(tmp_tsteps,Dgamma.Out{k_decay}{k_gamma}{k_sim}.TDs(tmp_tsteps),'k');
end
axis([Dgamma.Out{k_decay}{k_gamma}{k_sim}.goalsteps(k_trial_start-1)+[0 50] 0 1]);
set(A,'Box','off');
set(A,'PlotBoxAspectRatio',[2 1 1]);
set(A,'FontName','Ariel','FontSize',20);
tmp = ceil(Dgamma.Out{k_decay}{k_gamma}{k_sim}.goalsteps(k_trial_start-1)/10)*10;
set(A,'XTick',tmp+[0:10:50],'XTickLabel',tmp+[0:10:50]);
set(A,'YTick',[0:0.1:1],'YTickLabel',[0:0.1:1]);
if save_fig
    print(F,'-depsc','Fig6C');
end

% Fig. 10A-left
fig10Amax = ceil(max([max(max(Dalpha.ntspt_mean)), max(max(Dbeta.ntspt_mean)), max(max(Dgamma.ntspt_mean))]));
fig10Amin = 7;
F = figure;
A = axes;
hold on;
P = image(1+63*(Dalpha.ntspt_mean' - fig10Amin)/(fig10Amax - fig10Amin));
C = colorbar; set(C,'YTick',1+63*[0:1:fig10Amax-fig10Amin]/(fig10Amax-fig10Amin),'YTickLabel',[fig10Amin:1:fig10Amax]);
axis([0.5 length(decay_rate_set)+0.5 0.5 length(alpha_set)+0.5]);
set(A,'Box','off');
%set(A,'PlotBoxAspectRatio',[1 1 1]);
set(A,'FontName','Ariel','FontSize',20);
set(A,'XTick',[1:length(decay_rate_set)],'XTickLabel',decay_rate_set);
set(A,'YTick',[1:length(alpha_set)],'YTickLabel',alpha_set);
if save_fig
    print(F,'-depsc','Fig10A-left');
end

% Fig. 10A-middle
F = figure;
A = axes;
hold on;
P = image(1+63*(Dbeta.ntspt_mean' - fig10Amin)/(fig10Amax - fig10Amin));
C = colorbar; set(C,'YTick',1+63*[0:1:fig10Amax-fig10Amin]/(fig10Amax-fig10Amin),'YTickLabel',[fig10Amin:1:fig10Amax]);
axis([0.5 length(decay_rate_set)+0.5 0.5 length(beta_set)+0.5]);
set(A,'Box','off');
%set(A,'PlotBoxAspectRatio',[1 1 1]);
set(A,'FontName','Ariel','FontSize',20);
set(A,'XTick',[1:length(decay_rate_set)],'XTickLabel',decay_rate_set);
set(A,'YTick',[1:length(beta_set)],'YTickLabel',beta_set);
if save_fig
    print(F,'-depsc','Fig10A-middle');
end

% Fig. 10A-right
F = figure;
A = axes;
hold on;
P = image(1+63*(Dgamma.ntspt_mean' - fig10Amin)/(fig10Amax - fig10Amin));
C = colorbar; set(C,'YTick',1+63*[0:1:fig10Amax-fig10Amin]/(fig10Amax-fig10Amin),'YTickLabel',[fig10Amin:1:fig10Amax]);
axis([0.5 length(decay_rate_set)+0.5 0.5 length(gamma_set)+0.5]);
set(A,'Box','off');
%set(A,'PlotBoxAspectRatio',[1 1 1]);
set(A,'FontName','Ariel','FontSize',20);
set(A,'XTick',[1:length(decay_rate_set)],'XTickLabel',decay_rate_set);
set(A,'YTick',[1:length(gamma_set)],'YTickLabel',gamma_set);
if save_fig
    print(F,'-depsc','Fig10A-right');
end

% Fig. 10D-right
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
    print(F,'-depsc','Fig10D-right');
end
