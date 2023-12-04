
dire0='../data/';

t0=0;
dt=2;
ts=1000;
M=64+6;
N=M;
dn=1;
h=1;
plotType=0;
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
    
    % figure;
    if isfile([dire0 'phi_' num2str(t) '.dat'])
        clf;
        mu=reshape(importdata([dire0 'mu_' num2str(t) '.dat']),[N,M])';
        phi=reshape(importdata([dire0 'phi_' num2str(t) '.dat']),[N,M])';
        mu=mu(4:67,4:67);
        phi=phi(4:67,4:67);

        % clf;
        if plotType==0    
            subplot(1,2,1)
            imagesc(mu);
            % colro
            set(gca, 'YDir','normal')    
            title(['\fontsize{32} \mu, t=',num2str(t)]);
            axis off;
            axis square;
            subplot(1,2,2)            
            imagesc(phi);
            % imagesc(phi,[-0.01,0.01]);
            % colro
            set(gca, 'YDir','normal')    
            title(['\fontsize{32} \phi, t=',num2str(t)]);
            axis off;
            axis square;
        end    

        pause(tPaus);
    end
end