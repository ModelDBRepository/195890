% mkfig_Fig5

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
load used_rand_twister_for_Fig5

% parameters
num_sim = 20;
num_trial = 1000;
RLtype = 'Q';
p_alpha = 0.5;
p_beta = 5;
p_gamma = 1;
decay_rate = 0.01;
DAdep_paras = [0.25,501]; % depletion to 1/4 after 500 trials

% simulations
%rand('twister',sum(100*clock));
%Dsim.rand_twister = rand('twister');
Dsim.rand_twister = used_rand_twister_for_Fig5;
rand('twister',Dsim.rand_twister);

% reward
Rews{1} = zeros(1,23); Rews{1}([11,17]) = [0.5 1];
Rews{2} = zeros(1,23); Rews{2}([11,12]) = [0.5 1];

% main
for k1 = 1:2
    for k2 = 1:num_sim
        fprintf('%d-%d\n',k1,k2);
        Dsim.Out{k1}{k2} = RLdecayTmaze1(num_trial,RLtype,[p_alpha,p_beta,p_gamma],Rews{k1},decay_rate,DAdep_paras);
    end
end

% choose 2 ratio
bin = 10; % num_trial/bin should be an integer
for k1 = 1:2
    Dsim.choose2ratio{k1} = NaN(num_sim,num_trial/bin);
    for k2 = 1:num_sim
        Dsim.choose2ratio{k1}(k2,:) = sum(reshape(Dsim.Out{k1}{k2}.ArmChoices-1,bin,num_trial/bin),1)/bin;
    end
end

% time to reach state 4
bin = 10; % num_trial/bin should be an integer
for k1 = 1:2
    Dsim.avetime{k1} = NaN(num_sim,num_trial/bin); % all, to reach state 4
    for k2 = 1:num_sim
        tmp_times = diff([0;Dsim.Out{k1}{k2}.endsteps]);
        tmp_times3 = NaN(num_trial,1);
        for k_trial = 1:num_trial
            tmp = [0;Dsim.Out{k1}{k2}.endsteps];
            tmp_tsteps = [tmp(k_trial)+1:tmp(k_trial+1)]; % time steps (from start to end) for k_trial
            tmp_times3(k_trial) = find(Dsim.Out{k1}{k2}.States(tmp_tsteps)==4,1,'first');
        end
        Dsim.avetime{k1}(k2,:) = mean2(reshape(tmp_times3,bin,num_trial/bin),1);
    end
end

% save
save('DTmaze','Dsim');

% plot
save_fig = 1;

% Fig. 5B,F
tmpletters = 'BF';
for k1 = 1:2
    F = figure;
    A = axes;
    hold on;
    P = plot([0 num_trial/bin],[0 0],'k:');
    P = plot([0 num_trial/bin],[0.5 0.5],'k:');
    P = plot([0 num_trial/bin],[1 1],'k:');
    P = plot((1/2)*(num_trial/bin)*[1 1],[-0.1 1.1],'k.-');
    P = plot([0.5:1:num_trial/bin-0.5],mean2(Dsim.choose2ratio{k1},1)+sem2(Dsim.choose2ratio{k1},1,1,0),'k--');
    P = plot([0.5:1:num_trial/bin-0.5],mean2(Dsim.choose2ratio{k1},1)-sem2(Dsim.choose2ratio{k1},1,1,0),'k--');
    P = plot([0.5:1:num_trial/bin-0.5],mean2(Dsim.choose2ratio{k1},1),'k');
    axis([0 num_trial/bin -0.1 1.1]);
    set(A,'Box','off');
    %set(A,'PlotBoxAspectRatio',[1 1 1]);
    set(A,'FontName','Ariel','FontSize',20);
    set(A,'XTick',[0:100/bin:num_trial/bin],'XTickLabel',[0:100:num_trial]);
    set(A,'YTick',[0:0.1:1],'YTickLabel',[0:0.1:1]);
    if save_fig
        print(F,'-depsc',['Fig5' tmpletters(k1)]);
    end
end

% Fig. 5C,G
tmpletters = 'CG';
for k1 = 1:2
    F = figure;
    A = axes;
    hold on;
    P = plot([0 num_trial/bin],[4 4],'k:');
    P = plot((1/2)*(num_trial/bin)*[1 1],[3 8],'k.-');
    P = plot([0.5:1:num_trial/bin-0.5],mean2(Dsim.avetime{k1},1)+sem2(Dsim.avetime{k1},1,1,0),'k--');
    P = plot([0.5:1:num_trial/bin-0.5],mean2(Dsim.avetime{k1},1)-sem2(Dsim.avetime{k1},1,1,0),'k--');
    P = plot([0.5:1:num_trial/bin-0.5],mean2(Dsim.avetime{k1},1),'k');
    axis([0 num_trial/bin 3 8]);
    set(A,'Box','off');
    %set(A,'PlotBoxAspectRatio',[1 1 1]);
    set(A,'FontName','Ariel','FontSize',20);
    set(A,'XTick',[0:100/bin:num_trial/bin],'XTickLabel',[0:100:num_trial]);
    set(A,'YTick',[3:1:8],'YTickLabel',[3:1:8]);
    if save_fig
        print(F,'-depsc',['Fig5' tmpletters(k1)]);
    end
end

% Fig. 5D,H
for k1 = 1:2
    aveQs{k1} = zeros(num_trial,23);
    for k2 = 1:num_sim
        aveQs{k1} = aveQs{k1} + Dsim.Out{k1}{k2}.Qs/num_sim;
    end
end
Tmaze_value_plot(mean(aveQs{1}(251:500,:),1),jet,[6,0.8,60],'Fig5D-left');
Tmaze_value_plot(mean(aveQs{1}(751:1000,:),1),jet,[6,0.8,60],'Fig5D-right');
Tmaze_value_plot(mean(aveQs{2}(251:500,:),1),jet,[6,0.8,60],'Fig5H-left');
Tmaze_value_plot(mean(aveQs{2}(751:1000,:),1),jet,[6,0.8,60],'Fig5H-right');
