% mkfig_Fig2C4

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
load used_rand_twister_for_Fig2C4

% parameters
num_sim = 20;
num_trial = 500;
RLtype = 'Q';
p_alpha = 0.5;
p_beta = 5;
p_gamma = 1;
decay_rate_set = [0:0.002:0.02];
DAdep_factor_set = [0 0.1 0.2 0.25];
DAdep_start_trial = 251;

% simulations
%rand('twister',sum(100*clock));
%Dsim.rand_twister = rand('twister');
Dsim.rand_twister = used_rand_twister_for_Fig2C4;
rand('twister',Dsim.rand_twister);
for k1 = 1:length(decay_rate_set)
    for k2 = 1:length(DAdep_factor_set)
        Dsim.ntsptAllbin5{k1}{k2} = zeros(num_sim,100);
        for k3 = 1:num_sim
            fprintf('%d-%d-%d\n',k1,k2,k3);
            Dsim.Out{k1}{k2}{k3} = RLdecayStayGo2(num_trial,RLtype,[p_alpha,p_beta,p_gamma],decay_rate_set(k1),[DAdep_factor_set(k2),DAdep_start_trial]);
            Dsim.ntsptAllbin5{k1}{k2}(k3,:) = mean(reshape(diff([0;Dsim.Out{k1}{k2}{k3}.goalsteps]),5,100),1);
        end
        Dsim.ntsptAllbin5_mean{k1}{k2} = mean(Dsim.ntsptAllbin5{k1}{k2},1);
        Dsim.ntsptAllbin5_std{k1}{k2} = std(Dsim.ntsptAllbin5{k1}{k2},1,1);
    end
end
Ddep_decay = Dsim;
save Ddep_decay Ddep_decay

% plot
save_fig = 1;
load Ddep_decay
Dsim = Ddep_decay;

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

% Fig. 2C
tmp = 0;
for k1 = [1 6 11] %1:length(decay_rate_set)
    for k2 = [1 4] %1:length(DAdep_factor_set)
        tmp = max(tmp, max(Dsim.ntsptAllbin5_mean{k1}{k2}+Dsim.ntsptAllbin5_std{k1}{k2}/sqrt(num_sim)));
    end
end
YMax = ceil(tmp);
for k1 =  [1 6 11] %1:length(decay_rate_set)
    for k2 = [1 4] %1:length(DAdep_factor_set)
        F = figure;
        A = axes;
        hold on;
        P = plot([0 100],[7,7],'k:');
        P = plot([0 100],chance_step7*[1,1],'k:');
        P = plot(50*[1 1],[6 YMax],'k-.');
        P = plot([1:100]-0.5,Dsim.ntsptAllbin5_mean{k1}{k2},'k');
        P = plot([1:100]-0.5,Dsim.ntsptAllbin5_mean{k1}{k2}+Dsim.ntsptAllbin5_std{k1}{k2}/sqrt(num_sim),'k--');
        P = plot([1:100]-0.5,Dsim.ntsptAllbin5_mean{k1}{k2}-Dsim.ntsptAllbin5_std{k1}{k2}/sqrt(num_sim),'k--');
        axis([0 100 6 YMax]);
        set(A,'Box','off');
        %set(A,'PlotBoxAspectRatio',[1 1 1]);
        set(A,'FontName','Ariel','FontSize',20);
        set(A,'XTick',[0:10:100],'XTickLabel',[0:50:500]);
        set(A,'YTick',[6:1:YMax],'YTickLabel',[6:1:YMax]);
        if save_fig
            print(F,'-depsc',['Fig2C_' num2str(k1) '_' num2str(k2)]);
        end
    end
end

% Fig. 4
for k1 = 1:length(decay_rate_set)
    for k2 = 1:length(DAdep_factor_set)
        Dsim.Vend{k1}{k2} = zeros(num_sim,7*2);
        Dsim.Vs_whole_ave{k1}{k2} = zeros(num_trial,7*2);
        for k3 = 1:num_sim
            Dsim.Vend{k1}{k2}(k3,:) = Dsim.Out{k1}{k2}{k3}.Vs_whole(end,:);
            Dsim.Vs_whole_ave{k1}{k2} = Dsim.Vs_whole_ave{k1}{k2} + Dsim.Out{k1}{k2}{k3}.Vs_whole/num_sim;
        end
    end
end
% bottom panels
tmp = 0;
for k1 =  [1 6 11] %1:length(decay_rate_set)
    for k2 = [1 4] %1:length(DAdep_factor_set)
        tmp = max(tmp, max(mean(Dsim.Vend{k1}{k2}(:,2:2:12),1) + std(Dsim.Vend{k1}{k2}(:,2:2:12),1,1)/sqrt(num_sim)));
    end
end
Fig4bottomYmax = max(1, ceil(tmp*10)/10);
for k1 =  [1 6 11] %1:length(decay_rate_set)
    for k2 = [1 4] %1:length(DAdep_factor_set)
        F = figure;
        A = axes;
        hold on;
        P = plot([0.5 12.5],[1 1],'k:');
        P = errorbar([2:2:12],mean(Dsim.Vend{k1}{k2}(:,2:2:12),1),std(Dsim.Vend{k1}{k2}(:,2:2:12),1,1)/sqrt(num_sim),'kx-');
        set(P,'MarkerSize',20);
        P = plot([2:2:12],mean(Dsim.Vend{k1}{k2}(:,2:2:12),1),'k--');
        P = errorbar([1:2:11],mean(Dsim.Vend{k1}{k2}(:,1:2:11),1),std(Dsim.Vend{k1}{k2}(:,1:2:11),1,1)/sqrt(num_sim),'rx-');
        set(P,'MarkerSize',20,'Color',0.5*[1 1 1]);
        P = plot([1:2:11],mean(Dsim.Vend{k1}{k2}(:,1:2:11),1),'r--');
        set(P,'Color',0.5*[1 1 1]);
        axis([0.5 12.5 0 Fig4bottomYmax]);
        set(A,'Box','off');
        %set(A,'PlotBoxAspectRatio',[1 1 1]);
        set(A,'FontName','Ariel','FontSize',20);
        set(A,'XTick',[1:1:12],'XTickLabel',[1:1:12]);
        set(A,'YTick',[0:0.1:Fig4bottomYmax],'YTickLabel',[0:0.1:Fig4bottomYmax]);
        if save_fig
            print(F,'-depsc',['Fig4_bottom_panels_' num2str(k1) '_' num2str(k2)]);
        end
    end
end
% top panels
for k1 =  [1 6 11] %1:length(decay_rate_set)
    for k2 = [1 4] %1:length(DAdep_factor_set)
        F = figure;
        A = axes;
        hold on;
        P = image(1+63*flipud(Dsim.Vs_whole_ave{k1}{k2}(:,1:12)));
        C = colorbar; set(C,'YTick',1+63*[0:0.1:1],'YTickLabel',[0:0.1:1]);
        P = plot([0.5 12.5],250.5*[1 1],'k:');
        axis([0.5 12.5 0.5 500.5]);
        set(A,'Box','off');
        %set(A,'PlotBoxAspectRatio',[1 1 1]);
        set(A,'FontName','Ariel','FontSize',20);
        set(A,'XTick',[1:12],'XTickLabel',[1:12]);
        set(A,'YTick',0.5+[0:100:400],'YTickLabel',500-[0:100:400]);
        if save_fig
            print(F,'-depsc',['Fig4_top_pannels_' num2str(k1) '_' num2str(k2)]);
        end
    end
end
% without colorbar
for k1 =  [1 6 11] %1:length(decay_rate_set)
    for k2 = [1 4] %1:length(DAdep_factor_set)
        F = figure;
        A = axes;
        hold on;
        P = image(1+63*flipud(Dsim.Vs_whole_ave{k1}{k2}(:,1:12)));
        %C = colorbar; set(C,'YTick',1+63*[0:0.1:1],'YTickLabel',[0:0.1:1]);
        P = plot([0.5 12.5],250.5*[1 1],'k:');
        axis([0.5 12.5 0.5 500.5]);
        set(A,'Box','off');
        %set(A,'PlotBoxAspectRatio',[1 1 1]);
        set(A,'FontName','Ariel','FontSize',20);
        set(A,'XTick',[1:12],'XTickLabel',[1:12]);
        set(A,'YTick',0.5+[0:100:400],'YTickLabel',500-[0:100:400]);
        if save_fig
            print(F,'-depsc',['Fig4_top_pannels_wocolorbar_' num2str(k1) '_' num2str(k2)]);
        end
    end
end
