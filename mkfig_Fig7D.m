% mkfig_Fig7D

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
load used_rand_twister_for_Fig7D

% simulation
%rand('twister',sum(100*clock));
%Dsim.rand_twister = rand('twister');
Dsim.rand_twister = used_rand_twister_for_Fig7D;
rand('twister',Dsim.rand_twister);
Dsim.Out = RLdecayStayGo1(2000,'Q',[0.5,5,1],0.0045,[1 2001]);

% save
Dbi = Dsim;
save Dbi Dbi

% Fig. 7D
save_fig = 1;
F = figure;
A = axes;
hold on;
P = image(1+63*flipud(Dsim.Out.Vs_whole(:,1:12)));
C = colorbar; set(C,'YTick',1+63*[0:0.1:1],'YTickLabel',[0:0.1:1]);
axis([0.5 12.5 0.5 2000.5]);
set(A,'Box','off');
%set(A,'PlotBoxAspectRatio',[1 1 1]);
set(A,'FontName','Ariel','FontSize',20);
set(A,'XTick',[1:12],'XTickLabel',[1:12]);
set(A,'YTick',0.5+[0:500:1500],'YTickLabel',2000-[0:500:1500]);
if save_fig
    print(F,'-depsc','Fig7D');
end
% without colorbar
F = figure;
A = axes;
hold on;
P = image(1+63*flipud(Dsim.Out.Vs_whole(:,1:12)));
%C = colorbar; set(C,'YTick',1+63*[0:0.1:1],'YTickLabel',[0:0.1:1]);
axis([0.5 12.5 0.5 2000.5]);
set(A,'Box','off');
%set(A,'PlotBoxAspectRatio',[1 1 1]);
set(A,'FontName','Ariel','FontSize',20);
set(A,'XTick',[1:12],'XTickLabel',[1:12]);
set(A,'YTick',0.5+[0:500:1500],'YTickLabel',2000-[0:500:1500]);
if save_fig
    print(F,'-depsc','Fig7D_wocolorbar');
end
