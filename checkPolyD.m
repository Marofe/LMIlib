function out = checkPolyD(A,param)
% function out = checkPolyD(A,param)
%
% Check the polytope with gives vertices A1,A2,...,An for Discrete System.
%
% input:  A={A1,A2,...AN} -> state-space vertices matrices
%         param.dt -> (optional) used for ajust the grid (default dt=0.01)
%         param.n  -> (optional) used for ajust the number of
%         random points in the polytope
%        
% output: out.N             -> number of vertices
%         out.dim           -> dimension of the system
%         out.clock         -> times took for create the polytope
%         out.r           -> radius of the convex combination of the vertices
%         out.eig         -> eigvalues 
%         out.rV          -> radius of the vertices
%         out.eigV        -> eigvalues of the vertices
%	      out.alpha         -> values of simplex 
%	      out.maxR        -> gretest radiu in the polytope
%         out.alphaMaxR   -> simplex point of the greatest radiu
%
% E.g.
% A={randn(3),randn(3),randn(3)}
% poly=checkPolyD(A)
%
%
% Date: 23/09/2017
% Author: Marcos Rogério Fernandes
% Email: eng.marofe@hotmail.com
%% takes the parameters
out.N=length(A);
% take vertices eig values
for i=1:out.N
    out.eigV(i,:)=eig(A{i});
    out.rV(i,:)=abs(out.eigV(i,:));
end
out.V=A;
out.dim=size(A{1},1);
out.dt=0.01; %default grid
out.n=1000; %default number of randomness point in polytope
out.subpoly=1; %default take subpolytopes for N>3
out.links=1; %default take points on the links between vertices
takeToc=1;
if nargin == 2
    if isfield(param,'dt')
        out.dt=param.dt;
    end
    if isfield(param,'n')
        out.n=param.n;
    end
    if isfield(param,'subpoly')
        out.subpoly=param.subpoly;
    end
    if isfield(param,'links')
        out.links=param.links;
    end
    if isfield(param,'toc')
        takeToc=param.toc;
    end
end
if takeToc==1
    tic
end
a=0:out.dt:1;
n=length(a);
%% create polytope with N vertices
if out.N==1
    out.alpha=1;
    out.eig=eig(A{1});
    out.r=abs(out.eig);
elseif out.N==2
    for i=1:n
        out.eig(i,:)=eig(a(i)*A{1}+(1-a(i))*A{2});
        out.r(i,:)=abs(out.eig(i,:));
        out.alpha(i,1)=a(i);
        out.alpha(i,2)=1-a(i);
    end
elseif out.N>2
    k=1;
    %% polytope greater than 2
    if out.links==1
    L=combnk(1:out.N,2);
    for j=1:size(L,1)
        for i=1:n
            out.eig(k,:)=eig(a(i)*A{L(j,1)}+(1-a(i))*A{L(j,2)});
            out.r(k,:)=abs(out.eig(k,:));
            out.alpha(k,L(j,1))=a(i);
            out.alpha(k,L(j,2))=1-a(i);
            k=k+1;
        end
    end
    end
    %% one thousand of point equaly likely inside the polytope
    for n=1:out.n
    %generate alpha vector: alpha=(alpha1,alpha2,...alphaN)
    out.alpha(k,1)=1-rand^(1/(out.N-1));
    for j=2:out.N-1
        out.alpha(k,j)=(1-sum(out.alpha(k,1:j-1)))*(1-rand^(1/(out.N-j)));
    end
    out.alpha(k,out.N)=1-sum(out.alpha(k,1:out.N-1));
    %make convex combination using alpha
    Ac=zeros(out.dim,out.dim);
    for i=1:out.N
        Ac=Ac+out.alpha(k,i)*A{i};
    end
    out.eig(k,:)=eig(Ac);
    out.r(k,:)=abs(out.eig(k,:));
    k=k+1;
    end
    %% take subpolytopes
    if (out.N > 3) && (out.subpoly==1)
        L=combnk(1:out.N,3);
        param2.links=0;
        param2.n=150;
        for i=1:size(L,1)
            poly = checkPolyD({A{L(i,1)},A{L(i,2)},A{L(i,3)}},param2);
            out.r = [out.r; poly.r];
            out.eig = [out.eig; poly.eig];
            alpha = zeros(size(poly.alpha,1),out.N);
            alpha(:,L(i,1))=poly.alpha(:,1);
            alpha(:,L(i,2))=poly.alpha(:,2);
            alpha(:,L(i,3))=poly.alpha(:,3);
            out.alpha = [out.alpha; alpha];
        end
    end
end
v=max(out.r,[],2);
[out.maxR c]=max(v);
if out.N>1
out.alphaMaxR=out.alpha(c,:);   
end
out.maxEig=out.eig(c,:);
if takeToc==1
    out.clock=toc;
else
    out.clock=-1;
end
end