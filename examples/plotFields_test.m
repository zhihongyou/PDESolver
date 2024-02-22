
dire0='data/';

Nx=128;
Ny=128;
h=1;
figure;

for t=1:10
    % figure;
    f=reshape(importdata([dire0 'f_' num2str(t) '.dat']),[Ny,Nx]);
    phi=reshape(importdata([dire0 'phi_' num2str(t) '.dat']),[Ny,Nx]);
    % phirhs0=reshape(importdata([dire0 'phi_rhs_3.dat']),[Ny,Nx]);
    % f=f0(4:end-3,4:end-3);
    % phi=phi0(4:end-3,4:end-3);
    % phirhs=phirhs0(4:end-3,4:end-3);
    f1=2*phi+10*LaplO4(phi,h);
    clf;

    subplot(2,2,1)
    surf(phi);
    title(['t=' num2str(t)]);

    subplot(2,2,2)
    surf(f);

    subplot(2,2,3)
    surf(f1);

    subplot(2,2,4)
    surf(abs(f1-f));
    pause(0.1);
end