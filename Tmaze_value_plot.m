function Tmaze_value_plot(Vs,color_matrix,plotparas,savename)

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
% Tmaze_value_plot(Vs,color_matrix,plotparas,savename);
%color_matrix = jet; % get the colormap matrix "JET"

% set the colormap
colormap(color_matrix);

% colorbar
F = figure;
A = axes;
hold on;
C = colorbar; set(C,'YTick',1+63*[0:0.1:1],'YTickLabel',[0:0.1:1]);
set(A,'FontName','Ariel','FontSize',20);
if ~isempty(savename)
    print(gcf,'-depsc',[savename '_colorbar']);
end

% main figure
F = figure;
A = axes;
hold on;

% Go
arrow_width = plotparas(1); % 6;
arrow_length = plotparas(2); % 0.75;
P = plot([0 arrow_length],[0 0],'Line','-','LineWidth',arrow_width,'Color',color_matrix(round(1 + 63*Vs(2)),:));
P = plot(1+[0 arrow_length],[0 0],'Line','-','LineWidth',arrow_width,'Color',color_matrix(round(1 + 63*Vs(5)),:));
P = plot(2+[0 arrow_length],[0 0],'Line','-','LineWidth',arrow_width,'Color',color_matrix(round(1 + 63*Vs(8)),:));
P = plot([3 3],[0 arrow_length],'Line','-','LineWidth',arrow_width,'Color',color_matrix(round(1 + 63*Vs(12)),:));
P = plot([3 3],1+[0 arrow_length],'Line','-','LineWidth',arrow_width,'Color',color_matrix(round(1 + 63*Vs(17)),:));
P = plot([3 3],2+[0 arrow_length],'Line','-','LineWidth',arrow_width,'Color',color_matrix(round(1 + 63*Vs(23)),:));
P = plot([3 3],[0 -arrow_length],'Line','-','LineWidth',arrow_width,'Color',color_matrix(round(1 + 63*Vs(11)),:));
P = plot([3 3],-1+[0 -arrow_length],'Line','-','LineWidth',arrow_width,'Color',color_matrix(round(1 + 63*Vs(14)),:));
P = plot([3 3],-2+[0 -arrow_length],'Line','-','LineWidth',arrow_width,'Color',color_matrix(round(1 + 63*Vs(20)),:));

% Stay
circle_size = plotparas(3); % 50;
P = plot(0,0,'Marker','.','MarkerSize',circle_size,'Color',color_matrix(round(1 + 63*Vs(1)),:));
P = plot(1,0,'Marker','.','MarkerSize',circle_size,'Color',color_matrix(round(1 + 63*Vs(4)),:));
P = plot(2,0,'Marker','.','MarkerSize',circle_size,'Color',color_matrix(round(1 + 63*Vs(7)),:));
P = plot(3,0,'Marker','.','MarkerSize',circle_size,'Color',color_matrix(round(1 + 63*Vs(10)),:));
P = plot(3,1,'Marker','.','MarkerSize',circle_size,'Color',color_matrix(round(1 + 63*Vs(16)),:));
P = plot(3,2,'Marker','.','MarkerSize',circle_size,'Color',color_matrix(round(1 + 63*Vs(22)),:));
P = plot(3,-1,'Marker','.','MarkerSize',circle_size,'Color',color_matrix(round(1 + 63*Vs(13)),:));
P = plot(3,-2,'Marker','.','MarkerSize',circle_size,'Color',color_matrix(round(1 + 63*Vs(19)),:));

axis([-1 5 -3 3]);
set(A,'Box','off');
set(A,'PlotBoxAspectRatio',[1 1 1]);
set(A,'FontName','Ariel','FontSize',20);
set(A,'XTick',[],'XTickLabel',[]);
set(A,'YTick',[],'YTickLabel',[]);
if ~isempty(savename)
    print(F,'-depsc',[savename '_main']);
end
