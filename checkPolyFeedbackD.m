function out = checkPolyFeedbackD(sys,K,param)
% function out = checkPolyFeedback(sys,param)
%
% Check the robust state feedback for a polytopic system (A+BK).
%
% input:  poly=makePolyABCD() -> state-space polytopic system
%         K      -> Robust gain 
%         param.dt -> (optional) used for ajust the grid (default dt=0.01)
%         param.n  -> (optional) used for ajust the number of
%         random points in the polytope
%        
% output: out.N             -> number of vertices
%         out.dim           -> dimension of the system
%         out.clock         -> times took for create the polytope
%         out.eig           -> eigvalues of the convex combination of the vertices
%         out.eigV          -> eigvalues of the vertices
%	      out.alpha         -> values of simplex 
%	      out.maxEig        -> gretest eigvalue in the polytope
%         out.alphaMaxEig   -> simplex point of the greatest eigvalue
%
% E.g.
% poly=makePolyABCD(3,3,2,2,3,0.1);
% K=[1 1 1] 
% poly_cl=checkPolyFeedBack(poly,K)
%
%
% Date: 23/09/2017
% Author: Marcos Rogério Fernandes
% Email: eng.marofe@hotmail.com
%% takes the parameters
out.N=length(sys.A);
if isfield(sys,'B')
    sys.Bu=sys.B;
end
%% take closed loop vertices
for i=1:out.N
Acl{i}=sys.A{i}+sys.Bu{i}*K;
end
if nargin==2
    param=0;
end
%% check polytope
out=checkPolyD(Acl,param);
end