% mkfig_Fig8B-right

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

save_fig = 1;

% parameters to vary
RLparas = [0.5, 6, 1]; % [alpha (learning rate), beta (inverse temperature), gamma (time discount factor)]
ds = [0:0.0001:0.0195 0.01951:0.00001:0.0205 0.0206:0.0001:0.0410 0.04101:0.00001:0.0420 0.0421:0.0001:0.1]; % set of decay-degree
fzero_inits = [0.01:0.01:0.99]; % set of initial values used for "fzero" to solve the equation for Q1 ("formula_Q1 = 0")
difference_threshold = 0.0001; % threshold for distinguishing different equibrium points

% RL parameters
a = RLparas(1); % alpha (learning rate)
b = RLparas(2); % beta (inverse temperature
g = RLparas(3); % gamma (time discount factor)

% determine the values of Q1 at equibrium points
Q1s = zeros(length(ds),length(fzero_inits)); % value of "Stay" (Q(A_11) in Fig.1)
for k1 = 1:length(ds)
    fprintf('%d\n',k1);
    for k2 = 1:length(fzero_inits)
        Q1s(k1,k2) = fzero(@(Q1) formula_Q1(Q1, a, b, g, ds(k1)), fzero_inits(k2));
    end
end
dsM = ds' * ones(1,length(fzero_inits)); % transform "ds" to the matrix form (the same form as Q1s)
Q2s = a ./ (a + dsM); % value of "Go" (Q(A_12) in Fig.1)
P1s = 1 ./ (1 + exp(b*(Q2s - Q1s)));

% Jacobian matrix and eigen values
D11s = a * exp(b*(Q1s-Q2s)) .* (b*g*Q2s - b*Q1s - 1) - dsM;
D12s = a * exp(b*(Q1s-Q2s)) .* (b*Q1s - b*g*Q2s + g);
D21s = zeros(length(ds),length(fzero_inits));
D22s = -a - dsM;
lamda1 = (1/2)*(D11s + D22s + sqrt((D11s+D22s).^2 - 4*(D11s.*D22s - D12s.*D21s)));
lamda2 = (1/2)*(D11s + D22s - sqrt((D11s+D22s).^2 - 4*(D11s.*D22s - D12s.*D21s)));

% Types of equibrium points
Types_of_equibriums = zeros(length(ds),length(fzero_inits));
for k1 = 1:length(ds)
    for k2 = 1:length(fzero_inits)
        if isreal(lamda1(k1,k2)) % if eigen values are real numbers
            if max(lamda1(k1,k2),lamda2(k1,k2)) < 0 % asymptotically stable
                Types_of_equibriums(k1,k2) = 1;
            elseif max(lamda1(k1,k2),lamda2(k1,k2)) > 0 % unstable (for at least one direction)
                Types_of_equibriums(k1,k2) = 2;
            else
                Types_of_equibriums(k1,k2) = 3;
            end
        else
            if (D11s(k1,k2)+D22s(k1,k2))/2 < 0 % asymptotically stable
                Types_of_equibriums(k1,k2) = 4;
            elseif (D11s(k1,k2)+D22s(k1,k2))/2 > 0 % unstable
                Types_of_equibriums(k1,k2) = 5;
            else
                Types_of_equibriums(k1,k2) = 6;
            end
        end
    end
end
Symbols_for_plot{1} = '.';
Symbols_for_plot{2} = 'o';
Symbols_for_plot{3} = 'x';
Symbols_for_plot{4} = '.';
Symbols_for_plot{5} = 'o';
Symbols_for_plot{6} = 'x';

% extract solutions that are different with each other
num_equibriums = zeros(1,length(ds));
for k1 = 1:length(ds)
    fzero_init_indice{k1} = [];
    for k2 = 1:length(fzero_inits)
        if isempty(fzero_init_indice{k1})
            fzero_init_indice{k1} = k2;
        elseif min(abs(Q1s(k1,k2) - Q1s(k1,fzero_init_indice{k1}))) > difference_threshold
            fzero_init_indice{k1} = [fzero_init_indice{k1}, k2];
        end
    end
    num_equibriums(k1) = length(fzero_init_indice{k1});
end

% Fig. 8B-right
F = figure;
A = axes;
hold on;
for k1 = 1:length(ds)
    for kk = 1:length(fzero_init_indice{k1})
        k2 = fzero_init_indice{k1}(kk);
        %P = plot(ds(k1),Q1s(k1,k2),['b' Symbols_for_plot{Types_of_equibriums(k1,k2)}]);
        %P = plot(ds(k1),Q2s(k1,k2),['r' Symbols_for_plot{Types_of_equibriums(k1,k2)}]);
        if (Types_of_equibriums(k1,k2) == 1) || (Types_of_equibriums(k1,k2) == 4)
            tmpcolor1 = 'r';
            tmpcolor2 = 'b';
            tmpsize = 15;
        else
            tmpcolor1 = [1 0.7 0.7];
            tmpcolor2 = [0.7 0.7 1];
            tmpsize = 7.5;
        end
        P = plot(ds(k1),Q1s(k1,k2),'.','color',tmpcolor1,'MarkerSize',tmpsize);
        P = plot(ds(k1),Q2s(k1,k2),'.','color',tmpcolor2,'MarkerSize',tmpsize);
    end
end
axis([0 0.1 0 1]);
%set(A,'PlotBoxAspectRatio',[1 1 1]);
set(A,'Box','off');
set(A,'FontSize',16);
set(A,'XTick',[0:0.01:0.1],'XTickLabel',[0:0.01:0.1]);
set(A,'YTick',[0:0.1:1],'YTickLabel',[0:0.1:1]);
if save_fig
    print(F,'-depsc','Fig8B-right');
end
