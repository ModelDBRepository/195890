function out = Q1nullcline_Q1top(Q1, a, b, g, d, Q2)

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
% out = Q1nullcline_Q1top(Q1, a, b, g, d, Q2);
% "Q1nullcline_Q1top = 0" is the equation for the Q1 nullcline
%
%   Q1: value of "Stay"
%   Q2: value of "Go"
%   a: alpha (learning rate)
%   b: beta (inverse temperature)
%   g: gamma (time discount factor)
%   d: decay-degree (per trial on average, d=0 means no decay)

P1 = 1 / (1 + exp(b*(Q2 - Q1)));
out = (P1/(1-P1)) * a*(g*Q2 - Q1) - d*Q1;