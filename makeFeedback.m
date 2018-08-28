function out = makeFeedback(sys,param)
% function out = makeFeedback(sys,param)
%
% Synthesis State-feedback Control for Continuous (A,B_u,B_w,C,D_u) Polytopic System.
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
%         out.output        -> dimension of controlled output vector
%         out.clock         -> time to mount and solve the problem
%         out.K             -> robust gain that minimize H2 norm
%         out.feas          -> feasibility of LMI set
%         out.r             -> minimum primal residual
%         out.sol           -> solvesdp object result
%
% E.g.
% sys=makePolyABuBwCDu(3,2,2,3,1.1);
% ctrl=makeFeedback(sys)
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
    if isfield(param,'sigma')&&out.teorem~=3
        out.teorem=3;
    end
end
if takeToc==1
    tic
end
%% Robust  State Feedback Control Synthesis
%mount LMI set for each vertice
switch out.teorem
    case 1
        %% Teorem 1 (Quadratic Stability State Feedback )
        W=sdpvar(out.dim,out.dim,'symmetric');
        Z=sdpvar(out.m_u,out.dim,'full');
        lmis=[W>=0];
        for i=1:out.N
            lmis=[lmis (sys.A{i}*W+W*sys.A{i}'+sys.Bu{i}*Z+Z'*sys.Bu{i}'<=0)];
        end
    case 2
        %% Teorem 2 (Lyapunov Affine)
        G=sdpvar(out.dim,out.dim,'symmetric');
        Z=sdpvar(out.m2,out.dim,'full');
        lmis=[];
        for i=1:out.N
            W{i}=sdpvar(out.dim,out.dim,'symmetric');
            L11=W{i};
            L12=G'*sys.A{i}'+Z'*sys.Bu{i}';
            L21=L12';
            L22=G+G'-W{i};
            L=[L11 L12;L21 L22];
            lmis=[lmis W{i}>=0 L>=0];
        end
    case 3
        %% Greatest autovalue
        W=sdpvar(out.dim,out.dim,'symmetric');
        Z=sdpvar(out.m_u,out.dim,'full');
        lmis=[W>=0];
        for i=1:out.N
            lmis=[lmis (sys.A{i}*W+W*sys.A{i}'+sys.Bu{i}*Z+Z'*sys.Bu{i}'+2*param.sigma*W<=0)];
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
            out.K=double(Z)*inv(double(W));
            out.W=double(W);
        case 2
            out.K=double(Z)*inv(double(G));
            out.W=W;
        case 3
            out.K=double(Z)*inv(double(W));
            out.W=W;   
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