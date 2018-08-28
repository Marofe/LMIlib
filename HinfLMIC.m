function out = HinfLMIC(A,B,C,D,param)
% function out = HinfLMIC(y,param)
%
% Compute the Hinf Norm for Continuous System by LMI
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
L11=A'*P+P*A+C'*C;
L12=P*B+C'*D;
L22=D'*D-eye(size(B,2))*m;
L=[L11 L12;L12' L22];
lmis=[(P>=0) (L<=0)];
end
%% Bouded real lemma extended (Finsler)
if type==2
    H=sdpvar(2*n+np+nm,np+n,'full');
    DD = [zeros(n) P zeros(n,np) zeros(n,nm);
        P zeros(n) zeros(n,np) zeros(n,nm);...
        zeros(np,n) zeros(np,n) eye(np,np) zeros(np,nm);...
        zeros(nm,n) zeros(nm,n) zeros(nm,np) -m*eye(nm)];
    BB = [eye(n) -A zeros(np) -B;...
        zeros(np,n) -C eye(np) -D];
    PP = [P zeros(n) zeros(n,np) zeros(n,nm);
        zeros(n) eye(n) zeros(n,np) zeros(n,nm);...
        zeros(np,n) zeros(np,n) eye(np) zeros(np,nm);...
        zeros(nm,n) zeros(nm,n) zeros(nm,np) eye(nm)];
    lmis=[(PP>=0) (DD+H*BB+BB'*H'<=0)];
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