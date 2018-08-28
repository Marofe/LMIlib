function out = H2LMID(A,B,C,D,param)
% function out = H2LMID(A,B,C,D,param)
%
% Compute the H2 Norm for Discrete System by LMI
%
% input:  A,B,C,D ->  matrices of the system (without uncertainty)
%         param.type -> type of LMI: 1-Controlability Gramian
%                                    2-Controlability Extended (Finsler)
% 
%        
% output: out.feas -> feasibility of the LMI
%         out.H2 -> norm H2
%  
%
%
%
% Author: Marcos Rogério Fernandes
% E-mail: eng.marofe@hotmail.com
% Date: 23/09/2017

%% parameters

p=1e-7;
n=size(A,1); %number of states
nm=size(B,2); %number of inputs
np=size(C,1); %number of outputs
P=sdpvar(n);
type=1; %default type is Bounded Real Lemma
if nargin == 5
    if isfield(param,'type')
        type=param.type;
    end
end
tic
%% Controlability Gramian
if type==1
lmis=[(P>=0) (A*P*A'-P+B*B'<=0)];
end
%% Extended (Finsler)
if type==2
    F=sdpvar(n,n,'full');
    G=sdpvar(n,n,'full');
    
    L11=P-A*F'-F*A';
    L12=F-A*G;
    L13=B;
    L21=L12';
    L22=-P+G+G';
    L23=zeros(n,nm);
    L31=L13';
    L32=L23';
    L33=eye(nm);
    L=[L11 L12 L13;L21 L22 L23;L31 L32 L33];
    lmis=[(P>=0) (L>=0)];
end

setup=sdpsettings('verbose',0,'solver','sedumi');

%% solve LMI
sol=solvesdp(lmis,trace(C*P*C'),setup);
out.time=toc;
%% Test feasibility
out.l=min(checkset(lmis));
if out.l > -p
    out.feas=1;
    out.H2=sqrt(trace(C*double(P)*C'));
else
    out.feas=0;
    out.H2=-1;
end
out.type=type;
end