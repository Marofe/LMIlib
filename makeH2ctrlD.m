function out = makeH2ctrlD(sys,param)
% function out = makeH2ctrlD(sys,param)
%
% Synthesis H2 Control for Discrete (A,B_u,B_w,C,D_u,D_w) Polytopic System.
%
% input:  sys= -> state-space (A,B_u,B_w,C,D_u,D_w) polytopic system 
%         param.tol -> tolerance for primal residual (default:1e-7)
%         param.teorem -> teorem for derivate the LMIs
%                         -> 1 : Primal (A'PA-P+C'C<0)
%                         -> 2 : Dual (AWA'-W+BB'<0)
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
%         out.feas          -> feasibility of LMI set
%         out.r             -> minimum primal residual
%         out.sol           -> solvesdp object result
%
% E.g.
% sys=makePolyDABuBwCDuDw(3,2,2,3,1.1);
% h2ctrl=makeH2ctrlD(sys)
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
out.m_w=size(sys.Bw{1},2); %disturbe vector size
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
%% Robust  H2 Control Synthesis
%solver variables
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
        %% Teorem 1 (A'PA-P+C'C<0)
        X=sdpvar(out.m_w,out.m_w,'symmetric');
        lmis=[W>=0];
        for i=1:out.N
            L11=W;
            L12=W*sys.A{i}'+Z'*sys.Bu{i}';
            L13=W*sys.C{i}'+Z'*sys.Du{i}';
            L22=W;
            L23=zeros(out.dim,out.p);
            L33=eye(out.p);
            lmis=[lmis [X sys.Bw{i}';sys.Bw{i} W]>=0 [L11 L12 L13;L12' L22 L23;L13' L23' L33]>=0];
        end
    case 2
        %% Teorem 2 - Dual (AWA'-W+BB'<0)
        X=sdpvar(out.p,out.p,'symmetric');
        lmis=[W>=0];
        for i=1:out.N
            L11=W;
            L12=sys.A{i}*W+sys.Bu{i}*Z;
            L13=sys.Bw{i};
            L22=W;
            L23=zeros(out.dim,out.m_w);
            L33=eye(out.m_w);
            lmis=[lmis [X sys.C{i}*W+sys.Du{i}*Z;W*sys.C{i}'+Z'*sys.Du{i}' W]>=0 [L11 L12 L13;L12' L22 L23;L13' L23' L33]>=0];
        end
    case 3
        %% Teorem 3 - Extension Finsler Primal
        G=sdpvar(out.dim,out.dim,'full');
        X=sdpvar(out.m_w,out.m_w,'symmetric');
        lmis=[W>=0];
        for i=1:out.N
            L11=G+G'-W;
            L12=G'*sys.A{i}'+Z'*sys.Bu{i}';
            L13=G'*sys.C{i}'+Z'*sys.Du{i}';
            L22=W;
            L23=zeros(out.dim,out.p);
            L33=eye(out.p);
            lmis=[lmis [X sys.Bw{i}';sys.Bw{i} W]>=0 [L11 L12 L13;L12' L22 L23;L13' L23' L33]>=0];
        end
   case 4
        %% Teorem 4 - Extension Finsler Dual
        G=sdpvar(out.dim,out.dim,'full');
        X=sdpvar(out.p,out.p,'symmetric');
        lmis=[W>=0];
        for i=1:out.N
            L11=W;
            L12=sys.A{i}*G+sys.Bu{i}*Z;
            L13=sys.Bw{i};
            L21=L12';
            L22=G+G'-W;
            L23=zeros(out.dim,out.m_w);
            L31=L13';
            L32=L23';
            L33=eye(out.m_w);
            lmis=[lmis [X sys.C{i}*G+sys.Du{i}*Z;G'*sys.C{i}'+Z'*sys.Du{i}' G+G'-W]>=0 [L11 L12 L13;L21 L22 L23;L31 L32 L33]>=0];
        end      
    otherwise
        error('Teorem not found!')
end
if out.teorem > 0
%solve LMI
setup=sdpsettings('verbose',0,'solver','sedumi');
obj_min=trace(X);%default objetive function for minimization
if nargin==2
     if isfield(param,'obj_min')
        switch param.obj_min            
            case 0 
                %% without minimization
                obj_min=[];
                if isfield(param,'H2')
                    lmis=[lmis trace(X)==param.H2^2];
                end
            case 1
                %% with regularization of W
                obj_min=trace(X)+param.epsilon*trace(W);
        end
     end
end
out.sol=solvesdp(lmis,obj_min,setup);

%capture the result
out.W=double(W);
out.r=min(checkset(lmis));
if (abs(out.r)<out.tol||out.r>0)&&out.sol.problem==0
    switch out.teorem
        case 1
            out.K=double(Z)*inv(double(W));
        case 2
            out.K=double(Z)*inv(double(W));
        case 3
            out.K=double(Z)*inv(double(G));
        case 4
            out.K=double(Z)*inv(double(G));        
    end
    out.H2=sqrt(trace(double(X)))+trace(sys.Dw{1}'*sys.Dw{1});
    out.feas=1;
else
    warning('LMI infeasible')
    out.feas=0;
    out.info=[num2str(out.sol.problem) ' - ' out.sol.info];
end
else if out.N==1
    %% Riccati solution
    %A'PA-P-(A'PB2+C'D)inv(D'D+B2'PB2)(B2'PA+D'C)+C'C=0
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