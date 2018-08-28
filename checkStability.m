function out = checkStability(poly,type)
% function out = checkH2(sys,type,param)
%
% Check robust Stability for polytopic system.
%
% input:  sys=makePoly() -> state-space polytopic system 
%         type  -> (optional) used for ajust the LMI test
%        
% output: out.stable        -> flag 0: unfeasible; 1:stable
%         out.r             -> residual primal
%         out.time          -> yalmit+solver time
%
% E.g.
% sys=makePoly(3,2,3,-1);
% s=checkStability(sys)
%
%
% Date: 23/09/2017
% Author: Marcos Rogério Fernandes
% Email: eng.marofe@hotmail.com
%% takes the parameters
out.N=poly.N; %number of vertices
out.dim=poly.dim; %number of states
if nargin==1
 type='quadratic';
end
%% Test Quadratic Stability P(alpha)=P
if type=='quadratic'
    P=sdpvar(out.dim,out.dim,'symmetric');
    lmis=[P>=0];
    for i=1:out.N
        lmis=[lmis (poly.V{i}'*P+P*poly.V{i})<=0];
    end
    setup=sdpsettings('verbose',0,'solver','sedumi');
    sol=solvesdp(lmis,[],setup);
    out.time=sol.yalmiptime+sol.solvertime;
    out.r=min(checkset(lmis));
    if (out.r>0)&&(sol.problem==0)
        out.stable=1;
    else
        out.stable=0;
    end
end