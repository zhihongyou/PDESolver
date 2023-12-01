
dire0='../data/';

t0=0;
dt=1;
ts=1000;
M=128+6;
N=M;
dn=2;
h=1;
plotType=2;
tPaus=0.01;

% figure('Renderer', 'painters', 'Position', [500 300 1600 800]);
% figure('Renderer', 'painters', 'Position', [500 300 1000 500]);
figure('Renderer', 'painters', 'Position', [500 300 900 400]);
% figure('Renderer', 'painters', 'Position', [500 300 400 200]);
rhot=[];

for t=t0:dt:ts
    tStr=num2str(t);
    if mod(t,1)==0
        tStr=[tStr '.0'];
    end
    
    
    % figure;
    if isfile([dire0 'phivx_' num2str(t) '.dat'])
        clf;
        vx=reshape(importdata([dire0 'phivx_' num2str(t) '.dat']),[N,M])';
        vy=reshape(importdata([dire0 'phivy_' num2str(t) '.dat']),[N,M])';
        Qxx=reshape(importdata([dire0 'sigxx_' num2str(t) '.dat']),[N,M])';
        Qxy=reshape(importdata([dire0 'sigxy_' num2str(t) '.dat']),[N,M])';
        phi=reshape(importdata([dire0 'phi_' num2str(t) '.dat']),[N,M])';
        S2=reshape(importdata([dire0 'phi_' num2str(t) '.dat']),[N,M])';
        
        vx=vx(4:N-3,4:N-3);
        vy=vy(4:N-3,4:N-3);
        Qxx=Qxx(4:N-3,4:N-3);
        Qxy=Qxy(4:N-3,4:N-3);
        phi=phi(4:N-3,4:N-3);
        S2=S2(4:N-3,4:N-3);    
        S=2*sqrt(Qxx.^2+Qxy.^2);                       

        % clf;
        if plotType==0        
            subplot(2,3,1)        
            imagesc(Qxx);
            title(['t=' num2str(t)]);
            subplot(2,3,2)
            imagesc(Qxy);
            subplot(2,3,3)
            imagesc(S2);
            subplot(2,3,4)
            imagesc(S21);
            subplot(2,3,5)
            surf(S2-S21);
            % surf(dtQxy-dtQxy1);
            % subplot(2,3,6)
            % surf(dtomega-dtomega1);        
        elseif plotType==1        
            subplot(2,3,1)
            surf(dtQxx);
            title(['t=' num2str(t)]);
            subplot(2,3,2)
            surf(dtQxy);
            subplot(2,3,3)
            surf(dtomega);
            subplot(2,3,4)
            surf(dtQxx-dtQxx1);
            subplot(2,3,5)
            surf(dtQxy-dtQxy1);
            subplot(2,3,6)
            surf(dtomega-dtomega1);
        elseif plotType==2            
            theta=0.5*atan2(Qxy,Qxx);
            v=sqrt(vx.^2+vy.^2);
            vortex=(d1xO2(vy,h)-d1yO2(vx,h));
            vx1=vx(dn:dn:end,dn:dn:end);
            vy1=vy(dn:dn:end,dn:dn:end);
            Qxx1=Qxx(dn:dn:end,dn:dn:end);
            Qxy1=Qxy(dn:dn:end,dn:dn:end);
            v1=v(dn:dn:end,dn:dn:end);
            S1=S(dn:dn:end,dn:dn:end);    
            xc=[1,(N-6)/dn];
            yc=[1,(M-6)/dn];
            vortex1=vortex(dn:dn:end,dn:dn:end);
            theta1=theta(dn:dn:end,dn:dn:end);  
            
            subplot(1,3,1)
            imagesc(xc,yc,phi,[0 1]);
            colormap('Gray');
            colorbar;
            hold on;
            S1=S1./(S1+1e-10);
            quiver(S1.*cos(theta1),S1.*sin(theta1),0.4,'Color','k', ...
                   'ShowArrowHead','off','AutoScale','off');
            quiver(-S1.*cos(theta1),-S1.*sin(theta1),0.4,'Color','k', ...
                   'ShowArrowHead','off','AutoScale','off');
            set(gca, 'YDir','normal')    
            title(['\fontsize{32} Q, t=',num2str(t)]);
            axis off;
            axis square;
            subplot(1,3,2)
            imagesc(xc,yc,S,[min(min(S)),max(max(S))+1e-12]);
            colormap(gca,'summer');
            colorbar;
            hold on;
            S1=S1./(S1+1e-10);
            quiver(S1.*cos(theta1),S1.*sin(theta1),0.4,'Color','k', ...
                   'ShowArrowHead','off','AutoScale','off');
            quiver(-S1.*cos(theta1),-S1.*sin(theta1),0.4,'Color','k', ...
                   'ShowArrowHead','off','AutoScale','off');
            set(gca, 'YDir','normal')    
            title(['\fontsize{32} Q, t=',num2str(t)]);
            axis off;
            axis square;
            subplot(1,3,3)
            maxw=max(max(abs(vortex)));
            imagesc(xc,yc,vortex,[-maxw,maxw]);
            colormap(gca,'jet');
            colorbar;
            hold on;
            quiver(vx1,vy1,0.8,'Color','k');
            set(gca, 'YDir','normal')    
            title(['\fontsize{32} v, t=',num2str(t)]);
            axis off;
            axis square;
        end    

        pause(tPaus);
    end
end