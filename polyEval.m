function At = polyEval(A,alpha)
% function out = polyEval(A,alpha)
%
% Evaluation the polytopic system for a given alpha vector
%
% input:  A={A1,A2,...,An} -> vertices of polytope
%         alpha=[a1 a2 ... an], simplex vector
%        
% output: At             -> evaluation system in alpha
%
% E.g.
% A={randn(2),randn(2)};
% At=polyEval(A,[0.2 0.8])
%
%
% Date: 23/09/2017
% Author: Marcos Rogério Fernandes
% Email: eng.marofe@hotmail.com
%% Eval
if sum(alpha)==1
[n m]=size(A{1});
At=zeros(n,m);
for i=1:length(A)
    At = At+alpha(i)*A{i};
end
else
    error('Alpha does not belong to the unit simplex!')
end
end