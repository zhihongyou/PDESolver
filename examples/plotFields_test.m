
dire0='data/';

Nx=128+6;
Ny=128+6;
h=1;
figure;

for t=1:10
    f0=reshape(importdata([dire0 'f_' num2str(t) '.dat']),[Ny,Nx]);
    phi0=reshape(importdata([dire0 'phi_' num2str(t) '.dat']),[Ny,Nx]);
    phirhs0=reshape(importdata([dire0 'phi_rhs_3.dat']),[Ny,Nx]);
    f=f0(4:end-3,4:end-3);
    phi=phi0(4:end-3,4:end-3);
    phirhs=phirhs0(4:end-3,4:end-3);
    f1=LaplO4(phi,h);
    clf;

    subplot(2,2,1)
    surf(f);
    title(['t=' num2str(t)]);

    subplot(2,2,2)
    surf(phi);

    subplot(2,2,3)
    surf(phirhs);

    subplot(2,2,4)
    surf(abs(f1-f));
    pause(0.1);
end