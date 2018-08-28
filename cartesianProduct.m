function out=cartesianProduct(varargin)
% function out = cartesianProduct(v1,v2,...,vm)
%
% Compute the Cartesian Product of vectors v1,v2,...,vm 
%
% input:  v1,v2,...,vm ->  vectors 
%        
% output: out -> cartesian product result
%
% E.g.
% x1=[-1 1];
% x2=[-1 1];
% x1x2=cartesianProduct(x1,x2)
%
%
[X{1:nargin}] = ndgrid(varargin{:});
    for i=1:nargin
        T(:,i) = X{i}(:);
    end
    %out = T;
    out = unique(T , 'rows');
end