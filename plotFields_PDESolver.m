
dire0='/home/you/Dropbox/Apps/github/PDESolver/data/';

t0=1;
dt=1;
ts=100;
M=64+6;
N=M;
dn=1;
plotType=10;
tPaus=0.01;

% figure('Renderer', 'painters', 'Position', [500 300 1600 800]);
figure('Renderer', 'painters', 'Position', [500 300 1000 500]);
% figure('Renderer', 'painters', 'Position', [500 300 800 800]);
% figure('Renderer', 'painters', 'Position', [500 300 400 200]);
rhot=[];

for t=t0:dt:ts
    tStr=num2str(t);
    if mod(t,1)==0
        tStr=[tStr '.0'];
    end
    
    clf;
    % figure;
    phi=reshape(importdata([dire0 'phi_' num2str(t) '.dat']),[N,M])';
    rhs=reshape(importdata([dire0 'rhs_' num2str(t) '.dat']),[N,M])';
    lapl=reshape(importdata([dire0 'laplace_' num2str(t) '.dat']),[N,M])';
    subplot(1,3,1)
    imagesc(phi);
    title(num2str(t));
    subplot(1,3,2)
    imagesc(rhs);
    subplot(1,3,3)
    imagesc(lapl);
    rhs1=0*phi+2.5*(d2xO4(phi,1)+d2yO4(phi,1));
    lapl1=0*phi+1*(d2xO4(phi,1)+d2yO4(phi,1));
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