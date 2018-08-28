function out = HinfLMID(A,B,C,D,param)
% function out = HinfLMID(A,B,C,D,param)
%
% Compute the Hinf Norm for Discrete System by LMI
%
% input:  A,B,C,D ->  matrices of the system (without uncertainty)
%         param.type -> type of LMI: 1-Bouded Real Lema
%                                    2-Bouded Real Lema Extended (Finsler)
% 
%        
% output: out.feas -> feasibility of the LMI
%         out.Hinf -> norm Hinf
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
m=sdpvar(1);
type=1; %default type is Bounded Real Lemma
if nargin == 5
    if isfield(param,'type')
        type=param.type;
    end
end
tic
%% Bounded real lemma
if type==1
L11=A'*P*A-P+C'*C;
L12=A'*P*B+C'*D;
L22=B'*P*B+D'*D-eye(size(B,2))*m;
L=[L11 L12;L12' L22];
lmis=[(P>=0) (L<=0)];
end
%% Bouded real lemma extended (Finsler)
if type==2
    F1=sdpvar(n,n,'full');
    F2=sdpvar(n,n,'full');
    G1=sdpvar(n,np,'full');
    G2=sdpvar(n,np,'full');
    F3=sdpvar(np,n,'full');
    F4=sdpvar(nm,n,'full');
    G3=sdpvar(np,np,'full');
    G4=sdpvar(nm,np,'full');
    
    L11=P+F1+F1';
    L12=-F1*A-G1*C+F2';
    L13=G1+F3';
    L14=-F1*B-G1*D+F4';
    L21=L12';
    L22=-P-F2*A-A'*F2'-G2*C-C'*G2';
    L23=G2-A'*F3'-C'*G3';
    L24=-F2*B-G2*D-A'*F4'-C'*G4';
    L31=L13';
    L32=L23';
    L33=G3+G3'+eye(np);
    L34=-F3*B-G3*D+G4';
    L41=L14';
    L42=L24';
    L43=L34';
    L44=-F4*B-B'*F4'-G4*D-D'*G4'-m*eye(nm);
    L=[L11 L12 L13 L14;L21 L22 L23 L24;L31 L32 L33 L34;L41 L42 L43 L44];
    lmis=[(P>=0) (L<=0)];
end

setup=sdpsettings('verbose',0,'solver','sedumi');

%% solve LMI
sol=solvesdp(lmis,m,setup);
out.time=toc;
%% Test feasibility
out.l=min(checkset(lmis));
if out.l > -p
    out.feas=1;
    out.Hinf=sqrt(double(m));
else
    out.feas=0;
    out.Hinf=-1;
end
out.type=type;
end