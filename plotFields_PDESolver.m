
dire0='/home/you/Dropbox/Apps/github/PDESolver/data/';

t0=0;
dt=5;
ts=100;
M=64;
N=M;
dn=1;
plotType=0;
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
    
    phi0=importdata([dire0 'phi_' num2str(t) '.dat']);
    rho0=importdata([dire0 'phib_' num2str(t) '.dat']);    
    phi=reshape(phi0,[N M])';
    rho=reshape(rho0,[N M])';    
    rhot=[rhot;t sum(sum(rho))];

    clf;
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
