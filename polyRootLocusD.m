function out = polyRootLocusD(poly,animation,param)
% function out = polyRootLocusD(poly,animation)
%
% Plot the polytope root locus for discrete system with vertices A1,A2,...,An.
%
% input:  poly=makePolyD()/checkPolyD() -> polytope object
%         animation -> (optional) if specified, plot the eigvalues with
%         pause time, generating an animation. The value of animation will
%         be used for define the time. ex: animation=0.01, then
%         0.01s between each frame will be used.
%        
% output: none
%
% E.g.
% A={randn(3),randn(3),randn(3)}
% poly=createPoly(A)
% polyRootLocus(poly,0.01)
%
% Date: 23/09/2017
% Author: eng.marofe@hotmail.com
if nargin >=2
        figure
        plot(real(poly.eig),imag(poly.eig),'k.')
        hold on
        grid on
        title('Robust Root Locus for Discrete Systems')
        xlabel('Real axis (\sigma)')
        ylabel('Imaginary axis (j\omega)')
        viscircles([0 0],1,'Color','k','LineStyle','--')
        if isfield(param,'pho')
            viscircles([0 0],param.pho,'Color','red','LineStyle','--')    
        end
        axis('equal')
        %% plot eigvalues of vertices
        plot(real(poly.eigV),imag(poly.eigV),'g*')
        hold on
        Inst=max(abs(poly.eigV),[],2)>=1;
        if sum(Inst)>0
            plot(real(poly.eigV(Inst)),imag(poly.eigV(Inst)),'r*')
        end
        %% plot eigvalues of convex combinations
        for k = 1:size(poly.eig,1)
            if max(abs(poly.eig(k,:)))<1
              plot(real(poly.eig(k,:)),imag(poly.eig(k,:)),'*g')
            else
              plot(real(poly.eig(k,:)),imag(poly.eig(k,:)),'*r')
            end
            pause(animation)
        end
else
    figure
    Inst=max(abs(poly.eig),[],2)>=1;
    plot(real(poly.eig(~Inst,:)),imag(poly.eig(~Inst,:)),'g*')
    hold on
    if sum(Inst)>0
        plot(real(poly.eig(Inst,:)),imag(poly.eig(Inst,:)),'r*')
    end   
    grid on
    title('Robust Root Locus for Discrete Systems')
    xlabel('Real axis (\sigma)')
    ylabel('Imaginary axis (j\omega)')
    viscircles([0 0],1,'Color','k','LineStyle','--')
    axis('equal')
end
end