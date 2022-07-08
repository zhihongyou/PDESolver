
dire0='/home/you/Dropbox/Apps/github/PDESolver/data/';

t0=0;
dt=1;
ts=200;
M=64+6;
N=M;
dn=1;
h=1;
plotType=2;
tPaus=0.01;

% figure('Renderer', 'painters', 'Position', [500 300 1600 800]);
% figure('Renderer', 'painters', 'Position', [500 300 1000 500]);
figure('Renderer', 'painters', 'Position', [500 300 800 900]);
% figure('Renderer', 'painters', 'Position', [500 300 400 200]);
rhot=[];

for t=t0:dt:ts
    tStr=num2str(t);
    if mod(t,1)==0
        tStr=[tStr '.0'];
    end
    
    clf;
    % figure;
    vx=reshape(importdata([dire0 'incomFlow.vx_' num2str(t) '.dat']),[N,M])';
    vy=reshape(importdata([dire0 'incomFlow.vy_' num2str(t) '.dat']),[N,M])';
    omega=reshape(importdata([dire0 'incomFlow.omega_' num2str(t) '.dat']),[N,M])';
    Qxx=reshape(importdata([dire0 'Qxx_' num2str(t) '.dat']),[N,M])';
    Qxy=reshape(importdata([dire0 'Qxy_' num2str(t) '.dat']),[N,M])';
    S2=reshape(importdata([dire0 'S2_' num2str(t) '.dat']),[N,M])';
    dtomega=reshape(importdata([dire0 'incomFlow.omega_' num2str(t) '.dat']),[N,M])';
    dtQxx=reshape(importdata([dire0 'Qxx_' num2str(t) '.dat']),[N,M])';
    dtQxy=reshape(importdata([dire0 'Qxy_' num2str(t) '.dat']),[N,M])';
    
    vx=vx(4:67,4:67);
    vy=vy(4:67,4:67);
    Qxx=Qxx(4:67,4:67);
    Qxy=Qxy(4:67,4:67);
    S2=S2(4:67,4:67);    
    S=2*sqrt(Qxx.^2+Qxy.^2);
    omega=omega(4:67,4:67);    
    dtQxx=dtQxx(4:67,4:67);
    dtQxy=dtQxy(4:67,4:67);
    dtomega=dtomega(4:67,4:67);
    
    gamma=1;
    K=1;
    A2=-0.25;
    A4=1;
    lambda=0.7;
    eta=10;
    Gamma=0.0;
    alpha=-4;
    
    isFE=1;
    isCov=1;
    isFL=1;
    iswQ=1;
        
    S21=4*(Qxx.^2+Qxy.^2);
    dxvx=d1xO4(vx,h);
    dxvy=d1xO4(vy,h);
    dyvx=d1yO4(vx,h);
    dtQxx1=isFE/gamma*(K*(d2xO4(Qxx,1)+d2yO4(Qxx,1))-A2*Qxx-A4/2*S21.*Qxx);
    dtQxx1=dtQxx1-isCov*(vx.*d1xO4(Qxx,h)+vy.*d1yO4(Qxx,h))+isFL*lambda*dxvx-iswQ*(dxvy-dyvx).*Qxy;
    dtQxy1=isFE/gamma*(K*(d2xO4(Qxy,1)+d2yO4(Qxy,1))-A2*Qxy-A4/2*S21.*Qxy);
    dtQxy1=dtQxy1-isCov*(vx.*d1xO4(Qxy,h)+vy.*d1yO4(Qxy,h))+isFL*lambda/2*(dxvy+dyvx)+iswQ*(dxvy-dyvx).*Qxx;
    dtomega1=eta*(d2xO4(omega,1)+d2yO4(omega,1))-Gamma*omega;
    dtomega1=dtomega1+alpha*(d2xO4(Qxy,h)-d2yO4(Qxy,h)-2*d1x1yO4(Qxx,h));

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
        
        imagesc(xc,yc,S,[min(min(S)),max(max(S))+1e-12]);
        % imagesc(xc,yc,v,[min(min(v)),max(max(v))+1e-12]);
        colormap(gca,'summer');
        hold on;
        % kk=find(S1<1e-2);
        % S1(kk)=0;
        S1=S1./(S1+1e-10);
        quiver(S1.*cos(theta1),S1.*sin(theta1),0.4,'Color','k', ...
               'ShowArrowHead','off','AutoScale','off');
        quiver(-S1.*cos(theta1),-S1.*sin(theta1),0.4,'Color','k', ...
               'ShowArrowHead','off','AutoScale','off');
        set(gca, 'YDir','normal')    
        title(['\fontsize{32} \phi, t=',num2str(t)]);
        axis off;
        axis square;

    end    

    pause(tPaus);
end