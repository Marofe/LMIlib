function out = makeH2Filter(sys,param)
% function out = makeH2Filter(sys,param)
%
% Synthesis H2 filter for (A,B_w,C1,C2,D21) polytopic system.
%
% input:  sys= -> state-space (A,B_w,C1,C2,D21) polytopic system 
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
% h2filter=makeH2Filter(sys)
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
    sys.Bw=sys.B;
end
out.m_w=size(sys.Bw{1},2); %disturbe vector size
out.p=size(sys.C1{1},1); %controlled output vector size
out.q=size(sys.C2{1},1); %measurement output vector size
out.U=eye(out.dim); %default state space representation of filter
if nargin == 2
    if isfield(param,'toc')
        takeToc=param.toc;
    end
    if isfield(param,'U')
        out.U=param.U;
    end
end
if takeToc==1
    tic
end

%% Robust  H2 Filter Synthesis
%solver variables
X=sdpvar(out.dim,out.dim,'symmetric');
W=sdpvar(out.dim,out.dim,'symmetric');
F=sdpvar(out.p,out.dim,'full');
L=sdpvar(out.dim,out.q,'full');
G=sdpvar(out.dim,out.dim,'full');
M=sdpvar(out.m_w,out.m_w,'full')

%mount LMI set for each vertice
lmis=[X-W>=0];
for i=1:out.N
    T11=X;
    T12=W;
    T13=X*sys.Bw{i}+L*sys.D21{i};
    T22=W;
    T23=W*sys.Bw{i}+L*sys.D21{i};
    T33=M;
    lmis=[lmis [T11 T12 T13;T12' T22 T23;T13' T23' T33]>=0];
    T11=X*sys.A{i}+sys.A{i}'*X+sys.C2{i}'*L'+L*sys.C2{i};
    T12=sys.A{i}'*W+G+sys.C2{i}'*L';
    T13=sys.C1{i}';
    T22=G+G';
    T23=-F';
    T33=-eye(out.p);
    lmis=[lmis [T11 T12 T13;T12' T22 T23;T13' T23' T33]<=0];
end

%solve LMI
setup=sdpsettings('verbose',0,'solver','sedumi');
out.sol=solvesdp(lmis,trace(M),setup);

%acess solution
out.r=min(checkset(lmis));
if (out.r>0||abs(out.r)<1e-7)&&(out.sol.problem==0)
    out.feas=1;
    out.Af=inv(out.U')*double(G)*inv(double(W))*out.U';
    out.Bf=inv(out.U')*double(L);
    out.Cf=double(F)*inv(double(W))*out.U';
    out.Df=zeros(out.p,out.q);
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