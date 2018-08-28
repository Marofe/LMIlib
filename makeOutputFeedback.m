function out = makeOutputFeedback(sys,param)
% function out = makeOutputFeedback(sys,param)
%
% Synthesis Output feedback Control for Continuous (A,B_u,B_w,C,D_u) Polytopic System.
%
% input:  sys= -> state-space (A,B_u,B_w,C,D_u) polytopic system
%         param.teorem -> teorem for derivate the LMIs
%                         -> 1 : Quadratic Stability
%                         -> 2 : Lyapunov Affine
%
% output: out.N             -> number of vertices
%         out.dim           -> dimension of the system
%         out.m_w            -> dimension of disturbance vector
%         out.m_u           -> dimension of control input vector
%         out.p        -> dimension of controlled output vector
%         out.clock         -> time to mount and solve the problem
%         out.L             -> output feedback gain
%         out.feas          -> feasibility of LMI set
%         out.r             -> minimum primal residual
%         out.sol           -> solvesdp object result
%
% E.g.
% sys=makePolyABuBwCDu(3,2,2,3,1.1);
% ctrl=makeOutputFeedback(sys)
%
%
% Date: 5/11/2017
% Author: Marcos Rogério Fernandes
% Email: eng.marofe@hotmail.com
%% setup
out.N=length(sys.A);
out.dim=size(sys.A{1},1);
takeToc=1;
if isfield(sys,'B')
    sys.Bu=sys.B;
end
out.m_u=size(sys.Bu{1},2); %control input vector size
out.p=size(sys.C{1},1);
out.teorem=1; %default teorem for derivate LMIs conditions
if nargin == 2
    if isfield(param,'toc')
        takeToc=param.toc;
    end
    if isfield(param,'tol')
        out.tol=param.tol;
    end
    if isfield(param,'teorem')
        out.teorem=param.teorem;
    end
end
if takeToc==1
    tic
end
%% Robust  State Feedback Control Synthesis
%mount LMI set for each vertice
switch out.teorem
    case 1
        %% Teorem 1 C=[I | 0]
        W11=sdpvar(out.p,out.p,'symmetric');
        W22=sdpvar(out.dim-out.p,out.dim-out.p,'symmetric');
        Z11=sdpvar(out.m_u,out.p,'full');
        W=[W11 zeros(out.p,out.dim-out.p);zeros(out.dim-out.p,out.p) W22];
        Z=[Z11 zeros(out.m_u,out.dim-out.p)];
        lmis=[W>=0];
        for i=1:out.N
            lmis=[lmis (sys.A{i}*W+W*sys.A{i}'+sys.Bu{i}*Z+Z'*sys.Bu{i}'<=0)];
        end
otherwise
    error('Teorem not found!')
end
%solve LMI
setup=sdpsettings('verbose',0,'solver','sedumi');
out.sol=solvesdp(lmis,[],setup);

%capture the result
out.r=min(checkset(lmis));
if (out.r>0)&&(out.sol.problem==0)
    switch out.teorem
        case 1
            out.L=double(Z11)*inv(double(W11));
            out.W=double(W); 
    end
    out.feas=1;
else
    warning('LMI infeasible')
    out.feas=0;
    out.info=out.sol.info;
end
if takeToc==1
    out.clock=toc;
else
    out.clock=-1;
end
end