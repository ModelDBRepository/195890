function plot_nullclines(RLparas,d,savename)

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
% plot_nullclines(RLparas,d,savename);
%	RLparas: [alpha (learning rate), beta (inverse temperature), gamma (time discount factor)]
%	d: decay-degree
%	savename: name of the figure that is created by this function; if this is [], figure is not saved

% RL parameters
a = RLparas(1); % alpha (learning rate)
b = RLparas(2); % beta (inverse temperature
g = RLparas(3); % gamma (time discount factor)

% nullclines (calculate Q2 given Q1)
Q1s = [0:0.001:1];
Q2s_Q1nullcline = NaN(1,length(Q1s));
Q2s_Q2nullcline = NaN(1,length(Q1s));
for k = 1:length(Q1s)
    Q2s_Q1nullcline(k) = fzero(@(Q2) Q1nullcline_Q2top(Q2, a, b, g, d, Q1s(k)),1);
    Q2s_Q2nullcline(k) = fzero(@(Q2) Q2nullcline_Q2top(Q2, a, b, g, d, Q1s(k)),0.5);
end

% nullclines (calculate Q1 given Q2)
Q2s = [0:0.001:1];
Q1s_Q1nullcline = NaN(1,length(Q2s));
for k = 1:length(Q2s)
    Q1s_Q1nullcline(k) = fzero(@(Q1) Q1nullcline_Q1top(Q1, a, b, g, d, Q2s(k)),1);
end

% vector flow
flow_magnification = 0.5;
flowpointsX = [0:0.1:0.8 0.85 0.9 0.95 1];
flowpointsY = [0:0.2:0.6 0.75 0.85 0.9 0.95 1];
flowX = zeros(length(flowpointsX),length(flowpointsY));
flowY = zeros(length(flowpointsX),length(flowpointsY));
for k1 = 1:length(flowpointsX)
    for k2 = find(flowpointsX(k1) < flowpointsY, 1, 'first'):length(flowpointsY)
        Q1 = flowpointsX(k1);
        Q2 = flowpointsY(k2);
        P1 = exp(b*Q1) / (exp(b*Q1) + exp(b*Q2));
        flowX(k1,k2) = (P1/(1-P1)) * a * (g*Q2 - Q1) - d*Q1;
        flowY(k1,k2) = a * (1 - Q2) - d*Q2;
    end
end

% plot
F = figure;
A = axes;
hold on;
axis([0 1 0 1]);
P = plot([0 1],[0 1],'k:');
P = plot(Q1s,Q2s_Q1nullcline,'r.'); set(P,'MarkerSize',8);
P = plot(Q1s_Q1nullcline,Q2s,'r.'); set(P,'MarkerSize',8);
P = plot(Q1s,Q2s_Q2nullcline,'b.'); set(P,'MarkerSize',8);
for k1 = 1:length(flowpointsX)
    for k2 = find(flowpointsX(k1) < flowpointsY, 1, 'first'):length(flowpointsY)
        P = plot([flowpointsX(k1) flowpointsX(k1)+flow_magnification*flowX(k1,k2)],...
            [flowpointsY(k2) flowpointsY(k2)+flow_magnification*flowY(k1,k2)],'k');
    end
end
set(A,'PlotBoxAspectRatio',[1 1 1]);
set(A,'Box','off');
set(A,'FontSize',28);
set(A,'XTick',[0:0.1:1],'XTickLabel',[0:0.1:1]);
set(A,'YTick',[0:0.1:1],'YTickLabel',[0:0.1:1]);
if ~isempty(savename)
    print(F,'-depsc',savename);
end
