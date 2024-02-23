
dire0='data/';

Nx=128;
Ny=128;
h=1;
figure;

for t=1:10
    % figure;
    f=reshape(importdata([dire0 'f_' num2str(t) '.dat']),[Ny,Nx])';
    omega=reshape(importdata([dire0 'incompFlow.omega_' num2str(t) '.dat']),[Ny,Nx])';
    vx=reshape(importdata([dire0 'incompFlow.vx_' num2str(t) '.dat']),[Ny,Nx])';
    vy=reshape(importdata([dire0 'incompFlow.vy_' num2str(t) '.dat']),[Ny,Nx])';
    phi=reshape(importdata([dire0 'incompFlow.phi_' num2str(t) '.dat']),[Ny,Nx])';
    % phirhs0=reshape(importdata([dire0 'phi_rhs_3.dat']),[Ny,Nx]);
    % f=f0(4:end-3,4:end-3);
    % phi=phi0(4:end-3,4:end-3);
    % phirhs=phirhs0(4:end-3,4:end-3);
    f1=2*omega-10*LaplO4(omega,h);
    % f1=omega;
    omega1=-LaplO4(phi,h);
    omega2=d1xO4(vy,h)-d1yO4(vx,h);
    clf;

    subplot(2,3,1)
    surf(f);
    title(['t=' num2str(t)]);

    subplot(2,3,2)
    surf(f1);

    subplot(2,3,3)
    surf(f-f1);

    subplot(2,3,4)
    surf(omega);
    subplot(2,3,5)
    surf(omega2);
    subplot(2,3,6)
    surf(omega-omega2);
    pause(0.1);
end