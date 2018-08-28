function out = makeHinfCtrlD(sys,param)
% function out = makeHinfCtrlD(sys,param)
%
% Synthesis Hinf Control for Discrete (A,B_u,B_w,C,D_u) Polytopic System.
%
% input:  sys= -> state-space (A,B_u,B_w,C,D_u) polytopic system 
%         param.tol -> tolerance for primal residual (default:1e-7)
%         param.teorem -> teorem for derivate the LMIs
%                         -> 0 : Riccati approach
%                         -> 1 : Bounded Real Lemma
%         param.lyapOrder -> order for Lyapunov matrix
%                            -> 0 : Quadratic Stability (default)
%                            -> 1 : Affine Lyapunov (W=a1W1+a2W2)
%                            -> n>1 : Polynomial Lyapunov order n
%        
% output: out.N             -> number of vertices
%         out.dim           -> dimension of the system
%         out.m1            -> dimension of disturbance vector
%         out.m2            -> dimension of control input vector
%         out.output        -> dimension of controlled output vector
%         out.clock         -> time to mount and solve the problem
%         out.K             -> robust gain that minimize H2 norm
%         out.Hinf          -> Hinf norm value
%         out.feas          -> feasibility of LMI set
%         out.r             -> minimum primal residual
%         out.sol           -> solvesdp object result
%
% E.g.
% sys=makePolyABuBwCDu(3,2,2,3,1.1);
% hInfctrl=makeHinfCtrlD(sys)
%
%
% Date: 5/11/2017
% Author: Marcos Rogério Fernandes
% Email: eng.marofe@hotmail.com
%% setup
out.N=length(sys.A);
out.dim=size(sys.A{1},1);
takeToc=1;
out.m_w=size(sys.Bw{1},2); %disturbe vector size
if isfield(sys,'B')
    sys.Bu=sys.B;
end
out.m_u=size(sys.Bu{1},2); %control input vector size
out.tol = 1e-7; %default tolerance for primal residual
out.teorem=1; %default teorem for derivate LMIs conditions
out.p=size(sys.C{1},1); %output vector size
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
%% Robust  Hinf Control Synthesis
%solver variables
mu=sdpvar;
Z=sdpvar(out.m_u,out.dim,'full');
W=sdpvar(out.dim,out.dim,'symmetric');

%mount LMI set for each vertice
switch out.teorem
    case 0
        %% Riccati approach
        if out.N~=1
            error('Riccati approach is only for unique vertice polytope!')
        end
    case 1
        %% Teorem 1 (Bounded Real Lemma)
        lmis=[W>=0];
        for i=1:out.N
            L11=W;
            L12=sys.A{i}*W+sys.Bu{i}*Z;
            L13=zeros(out.dim,out.p);
            L14=sys.Bw{i};
            L22=W;
            L23=W*sys.C{i}'+Z'*sys.Du{i}';
            L24=zeros(out.dim,out.m_w);
            L33=eye(out.p);
            L34=zeros(out.p,out.m_w);
            L44=mu*eye(out.m_w);
            lmis=[lmis [L11 L12 L13 L14;L12' L22 L23 L24;L13' L23' L33 L34;L14' L24' L34' L44]>=0];
        end  
    case 2
        %% Teorem 2 (Extend version using Finsler)
        G=sdpvar(out.dim,out.dim,'full');
        lmis=[W>=0];
        for i=1:out.N
            L11=W;
            L12=sys.A{i}*G+sys.Bu{i}*Z;
            L13=zeros(out.dim,out.p);
            L14=sys.Bw{i};
            L22=G+G'-W;
            L23=G'*sys.C{i}'+Z'*sys.Du{i}';
            L24=zeros(out.dim,out.m_w);
            L33=eye(out.p);
            L34=zeros(out.p,out.m_w);
            L44=mu*eye(out.m_w);
            lmis=[lmis [L11 L12 L13 L14;L12' L22 L23 L24;L13' L23' L33 L34;L14' L24' L34' L44]>=0];
        end 
    otherwise
        error('Teorem not found!')
end
if out.teorem > 0
%solve LMI
setup=sdpsettings('verbose',0,'solver','sedumi');
obj_min=mu; %default objetive function for minimization 
if nargin==2
    if isfield(param,'obj_min')
        switch param.obj_min
            case 0
                obj_min=[];
                lmis=[lmis mu==param.Hinf^2];
            case 1
                epsilon=1e-5;
                if isfield(param,'epsilon')
                    epsilon=param.epsilon;
                end
                obj_min=mu+epsilon*trace(W);
        end
    end
end
out.sol=solvesdp(lmis,obj_min,setup);

%capture the result
out.W=double(W);
out.r=min(checkset(lmis));
if ((abs(out.r)<out.tol)||(out.r>0))&&(out.sol.problem==0)
    switch out.teorem
        case 1
            out.K=double(Z)*inv(double(W));  
        case 2
            out.K=double(Z)*inv(double(G));
    end
    out.Hinf=sqrt(double(mu));
    out.feas=1;
else
    warning('LMI infeasible')
    out.feas=0;
    out.info=[num2str(out.sol.problem) ' - ' out.sol.info];
end
else if out.N==1
    %% Riccati solution
    % F'XF - X -F'XB(B'XB + R)  B'XF + Q - SR  S'  with  F:=A-BR  S'.
    %A'PA-P+C'C+(A'PB2+C'D)inv(D'D+B2'PB2)(B2'PA+D'C)+C'C=0
    %S= 
    %K=-inv(D'D+B2'PB2)(B2'PA+D'C)
    [out.W out.eig out.K out.sol]=dare(sys.A{1},sys.Bu{1},sys.C{1}'*sys.C{1},sys.Du{1}'*sys.Du{1},sys.C{1}'*sys.Du{1});
    out.r=min(eig(out.W));
    if out.r>0
    out.feas=1;
    out.K=-inv(sys.Du{1}'*sys.Du{1}+sys.Bu{1}'*out.W*sys.Bu{1})*(sys.Bu{1}'*out.W*sys.A{1}+sys.Du{1}'*sys.C{1});
    out.H2=sqrt(trace(sys.Bw{1}'*out.W*sys.Bw{1}))+trace(sys.Dw{1}'*sys.Dw{1});
    else
        warning('Infeasible!')
        out.feas=0;
        out.info='No positive solution for Riccati Equation';
    end
end

if takeToc==1
    out.clock=toc;
else
    out.clock=-1;
end
end