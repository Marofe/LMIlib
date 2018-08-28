function out=makePoly(dim,N,d,param)
% function out = makePoly(dim,N,d,param)
%
% Make a polytope with random stable vertices A1,A2,...,An and the greatest
% eigvalue at d
%
% input:  A={A1,A2,...AN} -> state-space vertices matrices
%         param.dt -> (optional) used for ajust the grid (default dt=0.01)
%         param.n  -> (optional) used for ajust the number of
%         random points in the polytope
%        
% output: out.N             -> number of vertices
%         out.dim           -> dimension of the system
%         out.clock         -> times took for create the polytope
%         out.eig           -> eigvalues of the convex combination of the vertices
%         out.eigV          -> eigvalues of the vertices
%	  out.alpha         -> values of simplex 
%	  out.maxEig        -> gretest eigvalue in the polytope
%         out.alphaMaxEig   -> simplex value of the greatest eigvalue
%
% E.g.
% poly=makePoly(3,2,-1)
%
%
% Date: 23/09/2017
% Author: Marcos Rogério Fernandes 
% Email: eng.marofe@hotmail.com

if nargin==3
    param=0;
end
%% generate vertices stable and randomness
for i=1:N
    A{i}=randn(dim,dim);
    k=max(real(eig(A{i})));
    if k>0
        A{i}=A{i}-k*eye(dim);
    end
end
poly=checkPoly(A);
k=-d+poly.maxEig;
A=gsubtract(A,k*eye(dim));
out=checkPoly(A,param);
end