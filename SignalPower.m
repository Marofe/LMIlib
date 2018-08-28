function out = SignalPower(y,dt)
% function out = SignalPower(y,param)
%
% Compute the power of signal y
%
% input:  y ->  signal
%         dt -> (optional) used for ajust the sample time (default dt=0.01)
%        
% output: out -> power of signal y
%  
%
% E.g.
% y=randn(1,1000);
% yPower=SignalPower(y)
%
%
% Date: 23/09/2017
% Author: Marcos Rogério Fernandes
% E-mail: eng.marofe@hotmail.com

if nargin == 1
    dt = 0.01;
end
t = 0:dt:(length(y)-1)*dt;
out=sqrt(trapz(y.^2)*dt);
end