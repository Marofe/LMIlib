function out = makeHinfFilter(sys,param)
% function out = makeH2Filter(sys,param)
%
% Synthesis Hinf filter for (A,B_w,C1,C2,D21) polytopic system.
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
Z=sdpvar(out.dim,out.dim,'symmetric');
F=sdpvar(out.p,out.dim,'full');
L=sdpvar(out.dim,out.q,'full');
G=sdpvar(out.dim,out.dim,'full');
Df=sdpvar(out.p,out.q,'full')
mu=sdpvar;
%mount LMI set for each vertice
lmis=[X>=0 Z>=0 [Z Z;Z X]>=0];
for i=1:out.N
    T11=Z*sys.A{i}+sys.A{i}'*Z;
    T12=Z*sys.A{i}+sys.A{i}'*X+sys.C2{i}'*L'+G';
    T13=sys.C1{i}'-sys.C2{i}'*Df'-F';
    T14=Z*sys.Bw{i};
    T22=sys.A{i}'*X+X*sys.A{i}+sys.C2{i}'*L'+L*sys.C2{i};
    T23=sys.C1{i}'-sys.C2{i}'*Df';
    T24=X*sys.Bw{i}+L*sys.D21{i};
    T33=-eye(out.p);
    T34=sys.D11{i}-Df*sys.D21{i};
    T44=-mu*eye(out.m_w);
    lmis=[lmis [T11 T12 T13 T14;T12' T22 T23 T24;T13' T23' T33 T34;T14' T24' T34' T44]<=0];
end

%solve LMI
setup=sdpsettings('verbose',0,'solver','sedumi');
out.sol=solvesdp(lmis,mu,setup);

%acess solution
out.r=min(checkset(lmis));
if (out.r>0||abs(out.r)<1e-7)&&(out.sol.problem==0)
    out.feas=1;
    out.V=inv(out.U')*(eye(out.dim)-double(X)*inv(double(Z)));
    out.Af=inv(out.U')*double(G)*inv(out.V*double(Z));
    out.Bf=inv(out.U')*double(L);
    out.Cf=double(F)*inv(out.V*double(Z));
    out.Df=double(Df);
    out.Hinf=sqrt(double(mu));
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