function out = polyRootLocus(poly,animation,param)
% function out = polyRootLocus(poly,animation)
%
% Plot eigvalues cloud of the polytope
%
% input:  poly=makePoly()/checkPoly() -> polytope object
%         animation -> (optional) if specified, plot the eigvalues with
%         pause time, generating an animation. The value of animation will
%         be used for define the time. ex: animation=0.01, then
%         0.01s between each frame will be used.
%        
%        
% output: none

%
% E.g.
% A={randn(3),randn(3),randn(3)}
% poly=createPoly(A)
% polyRootLocus(poly)
%
%
% Date: 23/09/2017
% Author: eng.marofe@hotmail.com

if nargin >=2
        figure
        %plot(real(poly.eig),imag(poly.eig),'k.')
        hold on
        grid on
        title('Robust Root Locus')
        xlabel('Real axis (\sigma)')
        ylabel('Imaginary axis (j\omega)')
        if isfield(param,'sigma')
            line(-[param.sigma param.sigma],[-5 5],'LineStyle','--','Color','red');
        end
        %% plot eigvalues of vertices
        Inst=max(real(poly.eigV),[],2)>=0;
        plot(real(poly.eigV(~Inst,:)),imag(poly.eigV(~Inst,:)),'g*')
        if sum(Inst)>0
            plot(real(poly.eigV(Inst,:)),imag(poly.eigV(Inst,:)),'r*')
        end
        %% plot eigvalues of convex combinations
        for k = 1:size(poly.eig,1)
            if max(real(poly.eig(k,:)))<0
                plot(real(poly.eig(k,:)),imag(poly.eig(k,:)),'*g')
            else
                plot(real(poly.eig(k,:)),imag(poly.eig(k,:)),'*r')
            end
            pause(animation)
        end
else
    figure
    plot(real(poly.eig),imag(poly.eig),'k.')
    hold on
    Inst=max(real(poly.eig),[],2)>=0;
    plot(real(poly.eig(~Inst,:)),imag(poly.eig(~Inst,:)),'g*')
    if sum(Inst)>0
        plot(real(poly.eig(Inst,:)),imag(poly.eig(Inst,:)),'r*')
    end
    grid on
    title('Robust Root Locus')
    xlabel('Real axis (\sigma)')
    ylabel('Imaginary axis (j\omega)')
end
end