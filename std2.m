function out = std2(inp,flag,dim)

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
% taking std for numbers omitting NaN

m = size(inp,1);
n = size(inp,2);

if dim == 1
    out = NaN(1,n);
    for k1 = 1:n
        tmp = [];
        for k2 = 1:m
            if ~isnan(inp(k2,k1))
                tmp = [tmp, inp(k2,k1)];
            end
        end
        if ~isempty(tmp)
            out(k1) = std(tmp,flag);
        end
    end
elseif dim == 2
    out = NaN(m,1);
    for k1 = 1:m
        tmp = [];
        for k2 = 1:n
            if ~isnan(inp(k1,k2))
                tmp = [tmp, inp(k1,k2)];
            end
        end
        if ~isempty(tmp)
            out(k1) = std(tmp,flag);
        end
    end
end
