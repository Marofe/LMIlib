function out=makePolyD(dim,N,d,param)
% function out = makePolyD(dim,N,d,param)
%
% Make a polytope with random stable vertices A1,A2,...,An for discrete
% system and the greatest eigvalue radius at 1-d
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
%         out.r           -> radius of the convex combination of the vertices
%         out.rV          -> radius of the vertices
%	  out.alpha         -> values of simplex 
%	  out.maxEig        -> greatest eigvalue in the polytope
%     out.maxR          -> greatest radius in the polytope
%         out.alphaMaxEig   -> simplex value of the greatest eigvalue
%
% E.g.
% poly=makePolyD(3,2,-1)
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
    k=max(abs(eig(A{i})));
    if k>=1
        A{i}=A{i}/k;
    end
end
poly=checkPolyD(A);
k=poly.maxR/d;
A=gmultiply(A,1/k);
out=checkPolyD(A,param);
end