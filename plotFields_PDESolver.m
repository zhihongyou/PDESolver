
dire0='/home/you/Dropbox/Apps/github/PDESolver/data/';

t0=0;
dt=2;
ts=200;
M=64+6;
N=M;
dn=1;
plotType=10;
tPaus=0.01;

% figure('Renderer', 'painters', 'Position', [500 300 1600 800]);
% figure('Renderer', 'painters', 'Position', [500 300 1000 500]);
figure('Renderer', 'painters', 'Position', [500 300 800 400]);
% figure('Renderer', 'painters', 'Position', [500 300 400 200]);
rhot=[];

for t=t0:dt:ts
    tStr=num2str(t);
    if mod(t,1)==0
        tStr=[tStr '.0'];
    end
    
    clf;
    % figure;
    phia=reshape(importdata([dire0 'phia_' num2str(t) '.dat']),[N,M])';
    % mua=reshape(importdata([dire0 'mua_' num2str(t) '.dat']),[N,M])';
    laplacea=reshape(importdata([dire0 'phia_' num2str(t) '.dat']),[N,M])';
    mua=reshape(importdata([dire0 'phia_' num2str(t) '.dat']),[N,M])';
    phib=reshape(importdata([dire0 'phib_' num2str(t) '.dat']),[N,M])';
    mub=reshape(importdata([dire0 'phib_' num2str(t) '.dat']),[N,M])';
    laplaceb=reshape(importdata([dire0 'phib_' num2str(t) '.dat']),[N,M])';
    phia=phia(4:67,4:67);
    mua=mua(4:67,4:67);
    phib=phib(4:67,4:67);
    mub=mub(4:67,4:67);
    laplacea=laplacea(4:67,4:67);
    laplaceb=laplaceb(4:67,4:67);
    laplacea1=d2xO4(phia,1)+d2yO4(phia,1);
    laplaceb1=d2xO4(phib,1)+d2yO4(phib,1);
    mua1=laplacea1-0.5*laplaceb1;
    mub1=laplaceb1+0.5*laplacea1;
    % mua1=sin(phia);
    
    subplot(2,3,1)
    imagesc(phia);
    title(num2str(t));
    subplot(2,3,4)
    imagesc(phib);
    % subplot(2,3,2)
    % surf(mua1);
    % subplot(2,3,5)
    % surf(laplaceb-laplaceb1);
    % subplot(2,3,3)
    % surf(mua-mua1);
    % subplot(2,3,6)
    % surf(mub-mub1);
    

    % clf;
    if plotType==0
        subplot(1,2,1)
        imagesc(phi);
        % colormap(gca,'hsv');
        set(gca, 'YDir','normal')
        title(['\fontsize{32} \phi, t=',num2str(t)]);
        axis off;
        axis equal;
        
        subplot(1,2,2)
        imagesc(rho);
        % colormap(gca,'hsv');
        set(gca, 'YDir','normal')
        title(['\fontsize{32} \rho, t=',num2str(t)]);
        axis off;
        axis equal;
        
    elseif plotType==1
        subplot(1,2,1)
        surf(phi);
        title(['\fontsize{32} \phi, t=',num2str(t)]);
        
        subplot(1,2,2)
        surf(rho);

    end

    pause(tPaus);
end
% figure;
% % lapl=0*phi(4:67,4:67)+2.5*(d2xO4(phi(4:67,4:67),1)+d2yO4(phi(4:67,4:67),1));
% rhs1=0*phi+2.5*(d2xO4(phi,1)+d2yO4(phi,1));
% lapl1=0*phi+1*(d2xO4(phi,1)+d2yO4(phi,1));
% subplot(2,3,1)
% surf(lapl);
% subplot(2,3,2)
% surf(lapl1);
% subplot(2,3,3)
% surf(lapl-lapl1);
% subplot(2,3,4)
% surf(rhs);
% subplot(2,3,5)
% surf(rhs1);
% subplot(2,3,6)
% surf(rhs-rhs1)

% subplot(1,3,1)
% surf(phi);
% subplot(1,3,2)
% surf(rho);
% subplot(1,3,3)
% surf(phi-rho);