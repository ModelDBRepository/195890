% mkfig_Fig14

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
load used_rand_twister_for_Fig14

% parameters
num_sim = 20;
num_trial = 1000;
RLtype = 'Q';
p_alpha = 0.5;
p_beta = 5;
p_gamma = 1;
decay_rate = 0.01;
DAdep_paras = [0.25,501]; % depletion to 1/4 after 500 trials
velo_Stay_factor = 0.5;
rewarded_state = [5,8; 5,6];

% simulations
%rand('twister',sum(100*clock));
%Dsim.rand_twister = rand('twister');
Dsim.rand_twister = used_rand_twister_for_Fig14;
rand('twister',Dsim.rand_twister);

% reward
Rews{1} = zeros(1,20); Rews{1}([11,17]) = [0.5 1];
Rews{2} = zeros(1,20); Rews{2}([11,12]) = [0.5 1];

% main
for k1 = 1:2
    for k2 = 1:num_sim
        fprintf('%d-%d\n',k1,k2);
        Dsim.Out{k1}{k2} = RLdecayTmaze2(num_trial,RLtype,[p_alpha,p_beta,p_gamma],Rews{k1},decay_rate,DAdep_paras,velo_Stay_factor);
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
save('Dvelo','Dsim');

% plot
save_fig = 1;

% Fig. 14C,G
tmpletters = 'CG';
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
        print(F,'-depsc',['Fig14' tmpletters(k1)]);
    end
end

% Fig. 14D,H
tmpletters = 'DH';
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
        print(F,'-depsc',['Fig14' tmpletters(k1)]);
    end
end

% Fig. 14EI
tmpletters = 'EI';
for k1 = 1:2
    avelos{1}{k1} = NaN(num_sim,10); % all, before DA depletion
    avelos{2}{k1} = NaN(num_sim,10); % all, after DA depletion
    avelos{3}{k1} = NaN(num_sim,10); % after stay, before DA depletion
    avelos{4}{k1} = NaN(num_sim,10); % after stay, after DA depletion
    for k2 = 1:num_sim
        for k_state = 1:10
            tmp1 = zeros(length(Dsim.Out{k1}{k2}.States),1);
            tmp2 = zeros(length(Dsim.Out{k1}{k2}.States),1);
            tmp1(Dsim.Out{k1}{k2}.endsteps(250)+1:Dsim.Out{k1}{k2}.endsteps(500)) = 1;
            tmp2(Dsim.Out{k1}{k2}.endsteps(750)+1:Dsim.Out{k1}{k2}.endsteps(1000)) = 1;
            tmp_velo1 = Dsim.Out{k1}{k2}.Velocities((Dsim.Out{k1}{k2}.States==k_state)&tmp1);
            tmp_velo2 = Dsim.Out{k1}{k2}.Velocities((Dsim.Out{k1}{k2}.States==k_state)&tmp2);
            tmp_afterstay = (Dsim.Out{k1}{k2}.States == [0;Dsim.Out{k1}{k2}.States(1:end-1)]);
            tmp_velo3 = Dsim.Out{k1}{k2}.Velocities((Dsim.Out{k1}{k2}.States==k_state)&tmp1&tmp_afterstay);
            tmp_velo4 = Dsim.Out{k1}{k2}.Velocities((Dsim.Out{k1}{k2}.States==k_state)&tmp2&tmp_afterstay);
            avelos{1}{k1}(k2,k_state) = mean2(tmp_velo1,1);
            avelos{2}{k1}(k2,k_state) = mean2(tmp_velo2,1);
            avelos{3}{k1}(k2,k_state) = mean2(tmp_velo3,1);
            avelos{4}{k1}(k2,k_state) = mean2(tmp_velo4,1);
        end
    end
    F = figure;
    A = axes;
    hold on;
    P = plot([0.5 7.5],[0 0],'b');
    P = plot([0.5 7.5],[1 1],'b');
    if k1 == 1
        % all
        P = errorbar([1:4 4.85 5.85 6.85],mean2(avelos{1}{k1}(:,[1:4 6 8 10]),1),sem2(avelos{1}{k1}(:,[1:4 6 8 10]),1,1,0),'k');
        P = plot([1:4 4.85 5.85 6.85],mean2(avelos{1}{k1}(:,[1:4 6 8 10]),1),'k--');
        P = errorbar([4 5.15 6.15],mean2(avelos{1}{k1}(:,[4 5 9]),1),sem2(avelos{1}{k1}(:,[4 5 9]),1,1,0),'k');
        P = plot([4 5.15 6.15],mean2(avelos{1}{k1}(:,[4 5 9]),1),'k:');
        P = errorbar([1:4 4.85 5.85 6.85],mean2(avelos{2}{k1}(:,[1:4 6 8 10]),1),sem2(avelos{2}{k1}(:,[1:4 6 8 10]),1,1,0),'r');
        P = plot([1:4 4.85 5.85 6.85],mean2(avelos{2}{k1}(:,[1:4 6 8 10]),1),'r--');
        P = errorbar([4 5.15 6.15],mean2(avelos{2}{k1}(:,[4 5 9]),1),sem2(avelos{2}{k1}(:,[4 5 9]),1,1,0),'r');
        P = plot([4 5.15 6.15],mean2(avelos{2}{k1}(:,[4 5 9]),1),'r:');
        % after stay
        P = errorbar([1:4 4.85 5.85 6.85],mean2(avelos{3}{k1}(:,[1:4 6 8 10]),1),sem2(avelos{3}{k1}(:,[1:4 6 8 10]),1,1,0),'c');
        P = plot([1:4 4.85 5.85 6.85],mean2(avelos{3}{k1}(:,[1:4 6 8 10]),1),'c--');
        P = errorbar([4 5.15 6.15],mean2(avelos{3}{k1}(:,[4 5 9]),1),sem2(avelos{3}{k1}(:,[4 5 9]),1,1,0),'c');
        P = plot([4 5.15 6.15],mean2(avelos{3}{k1}(:,[4 5 9]),1),'c:');
        P = errorbar([1:4 4.85 5.85 6.85],mean2(avelos{4}{k1}(:,[1:4 6 8 10]),1),sem2(avelos{4}{k1}(:,[1:4 6 8 10]),1,1,0),'m');
        P = plot([1:4 4.85 5.85 6.85],mean2(avelos{4}{k1}(:,[1:4 6 8 10]),1),'m--');
        P = errorbar([4 5.15 6.15],mean2(avelos{4}{k1}(:,[4 5 9]),1),sem2(avelos{4}{k1}(:,[4 5 9]),1,1,0),'m');
        P = plot([4 5.15 6.15],mean2(avelos{4}{k1}(:,[4 5 9]),1),'m:');
    elseif k1 == 2
        % all
        P = errorbar([1:4 4.85 5.85],mean2(avelos{1}{k1}(:,[1:4 6 10]),1),sem2(avelos{1}{k1}(:,[1:4 6 10]),1,1,0),'k');
        P = plot([1:4 4.85 5.85],mean2(avelos{1}{k1}(:,[1:4 6 10]),1),'k--');
        P = errorbar([4 5.15 6.15],mean2(avelos{1}{k1}(:,[4 5 9]),1),sem2(avelos{1}{k1}(:,[4 5 9]),1,1,0),'k');
        P = plot([4 5.15 6.15],mean2(avelos{1}{k1}(:,[4 5 9]),1),'k:');
        P = errorbar([1:4 4.85 5.85],mean2(avelos{2}{k1}(:,[1:4 6 10]),1),sem2(avelos{2}{k1}(:,[1:4 6 10]),1,1,0),'r');
        P = plot([1:4 4.85 5.85],mean2(avelos{2}{k1}(:,[1:4 6 10]),1),'r--');
        P = errorbar([4 5.15 6.15],mean2(avelos{2}{k1}(:,[4 5 9]),1),sem2(avelos{2}{k1}(:,[4 5 9]),1,1,0),'r');
        P = plot([4 5.15 6.15],mean2(avelos{2}{k1}(:,[4 5 9]),1),'r:');
        % after stay
        P = errorbar([1:4 4.85 5.85],mean2(avelos{3}{k1}(:,[1:4 6 10]),1),sem2(avelos{3}{k1}(:,[1:4 6 10]),1,1,0),'c');
        P = plot([1:4 4.85 5.85],mean2(avelos{3}{k1}(:,[1:4 6 10]),1),'c--');
        P = errorbar([4 5.15 6.15],mean2(avelos{3}{k1}(:,[4 5 9]),1),sem2(avelos{3}{k1}(:,[4 5 9]),1,1,0),'c');
        P = plot([4 5.15 6.15],mean2(avelos{3}{k1}(:,[4 5 9]),1),'c:');
        P = errorbar([1:4 4.85 5.85],mean2(avelos{4}{k1}(:,[1:4 6 10]),1),sem2(avelos{4}{k1}(:,[1:4 6 10]),1,1,0),'m');
        P = plot([1:4 4.85 5.85],mean2(avelos{4}{k1}(:,[1:4 6 10]),1),'m--');
        P = errorbar([4 5.15 6.15],mean2(avelos{4}{k1}(:,[4 5 9]),1),sem2(avelos{4}{k1}(:,[4 5 9]),1,1,0),'m');
        P = plot([4 5.15 6.15],mean2(avelos{4}{k1}(:,[4 5 9]),1),'m:');
    end
    axis([0.5 7.5 -0.1 1.1]);
    set(A,'Box','off');
    %set(A,'PlotBoxAspectRatio',[1 1 1]);
    set(A,'FontName','Ariel','FontSize',20);
    set(A,'XTick',[1:4 4.85 5.15 5.85 6.15 6.85 7.15],'XTickLabel',[]);
    set(A,'YTick',[0:0.1:1],'YTickLabel',[0:0.1:1]);
    if save_fig
        print(F,'-depsc',['Fig14' tmpletters(k1)]);
    end
end
