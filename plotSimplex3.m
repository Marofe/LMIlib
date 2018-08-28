function out = plotSimplex3(poly,t)
% function out = plotSimplex3(poly)
%
% Plot the simplex visited values for polytopes with three vertices (only).
%
% input:  poly -> polytope generate by makePoly or createPoly function
%         t    -> if specified will be used for pause time, thus generating
%         an animation
%        
% output: none
%
% E.g.
% A={randn(3),randn(3),randn(3)};
% poly=createPoly(A);
% plotSimplex3(poly)
%
%
% Date: 23/09/2017
% Author: eng.marofe@hotmail.com
dt=0;
if nargin==2
    dt=t;
end
if poly.N == 3
figure
plot3(poly.alpha(1,1),poly.alpha(1,2),poly.alpha(1,3),'*b')
hold on
grid on
title('Simplex')
xlabel('\alpha_!')
ylabel('\alpha_2')
zlabel('\alpha_3')
axis([0 1.2 0 1.2 0 1.2])
view(60,30)
if dt>0
    pause(dt)
end
for i=2:length(poly.alpha)
    plot3(poly.alpha(i,1),poly.alpha(i,2),poly.alpha(i,3),'*b');
    if dt>0
    pause(dt)
    end
end
else
    disp('Number of vertices is not trhee!')
end
end