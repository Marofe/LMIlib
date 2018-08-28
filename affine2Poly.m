function [A V]=affine2Poly(Am,t)
% function [A V] = affine2Poly(A,t)
%
% Transform Affine Uncertain model into Polytopic Uncertain Model 
%
% input:  A={A0,A1,A2,...,Am} ->  matrices of the system
%         t={[mint1 maxt1],[mint2 maxt2],...,[mintm maxtm]}->  parameters interval          
%        
% output: out.A -> 2^m vertices of the Polytope model
%         out.V -> Cartesian Product of t1,t2,...,tm
%
% E.g.
% t={[-1 1],[-1 1]};
% At={randn(3),randn(3),randn(3)}
% out=affine2Poly(A,t)
%
% Author: Marcos Rogério Fernandes
% E-mail: eng.marofe@hotmail.com
% Date: 23/09/2017

V=cartesianProduct(t{:});
n=size(Am{1},1);
m=size(Am{1},2);
out.dim=n;
out.dim2=m;
out.N=size(V,1);
for i=1:size(V,1)
    Ac=zeros(n,m);
    for j=1:size(V,2)
        Ac=Ac+V(i,j)*Am{j+1};
    end
    A{i}=Am{1}+Ac;
end
end