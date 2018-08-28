function out = makeH2ctrl(sys,param)
% function out = makeH2ctrl(sys,param)
%
% Synthesis H2 control for (A,B_u,B_w,C,D_u,D_w) polytopic system.
%
% input:  sys= -> state-space (A,B_u,B_w,C,D_u,D_w) polytopic system 
%        
% output: out.N             -> number of vertices
%         out.dim           -> dimension of the system
%         out.clock         -> times took for create the polytope
%         out.K             -> robust gain with minimize H2 norm
%         out.feas          -> feasibility of LMIs
%         out.r             -> minimum primal residual
%         out.sol           -> solvesdp result
%
% E.g.
% sys=makePolyABCD(3,2,3,-1);
% K=makeH2ctrl(sys)
%
%
% Date: 5/11/2017
% Author: Marcos Rogério Fernandes
% Email: eng.marofe@hotmail.com
%% setup
out.N=length(sys.A);
out.dim=size(sys.A{1},1);
takeToc=1;
out.m1=size(sys.Bw{1},2); %disturbe vector size
m_w=out.m1;
out.m2=size(sys.Bu{1},2); %control input vector size
m_u=out.m2;
p=size(sys.C{1},1); %output vector size
if nargin == 2
    if isfield(param,'toc')
        takeToc=param.toc;
    end
end
if takeToc==1
    tic
end

%% Robust  H2 Control Synthesis
%solver variables
X=sdpvar(out.m2,out.m2,'symmetric');
Z=sdpvar(out.m2,out.dim,'full');
W=sdpvar(out.dim,out.dim,'symmetric');

%mount LMI set for each vertice
lmis=[W>=0];
for i=1:out.N
%     M=[X sys.Bw{i}';sys.Bw{i} W];
%     L11=sys.A{i}*W+W*sys.A{i}'+Z'*sys.Bu{i}'+sys.Bu{i}*Z;
%     L12=W*sys.C{i}'+Z'*sys.D{i}';
%     L22=-eye(p);
%     L=[L11 L12;L12' L22];
%     lmis=[lmis M>=0 -L>=0];
M=[X sys.C{i}*W+sys.Du{i}*Z;W*sys.C{i}'+Z'*sys.Du{i} W];
L11=W*sys.A{i}'+sys.A{i}*W+Z'*sys.Bu{i}'+sys.Bu{i}*Z;
L12=sys.Bw{i};
L21=L12';
L22=-eye(m_w);
L=[L11 L12;L21 L22];
lmis=[lmis M>=0 L<=0];
end

%solve LMI
out.sol=solvesdp(lmis,trace(X));

%acess
out.r=min(checkset(lmis));
if (out.r>0||abs(out.r)<1e-7)&&(out.sol.problem==0)
    out.K=double(Z)*inv(double(W));
    out.H2=sqrt(trace(double(X)));
    out.W=double(W);
    out.feas=1;
else
    warning('LMI infactivel')
    out.feas=0;
end

if takeToc==1
    out.clock=toc;
else
    out.clock=-1;
end
end