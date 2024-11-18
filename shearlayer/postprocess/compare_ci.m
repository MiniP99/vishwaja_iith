clear,clf
inputdir='/home/vsj/Codes/plot/';
outputdir='/home/vsj/Codes/plot/myplots/';
outputdir1='/home/vsj/Codes/plot/myplots/compare_ci/';
mkdir(outputdir1);
Lx = 80; Ly = 80; Lz = 40;
nx = 180; ny = 240; nz = 84;
dx = Lx/nx; x = linspace(0,Lx-dx,nx); 
dy = Ly/ny; y = linspace(-Ly/2,Ly/2,ny); 
dz = Lz/nz; z = linspace(0,Lz-dz,nz); 
c=sqrt(1.4);Mc=2.0;du=Mc*2*c;rho_0=1;Re_theta=1000;
%delta_0=Re_theta/(rho_0*Re*du);
delta_0=1;
figure(1)
fname1='/home/vsj/Codes/runs_amd/amd_data/vel_tke_flux_tavg_nor_2.0_amd_ci2.dat';
fname2='/home/vsj/Codes/runs_amd/amd_data/vel_tke_flux_tavg_nor_2.0_amd_ci3.dat';
f1 = readmatrix(fname1,'ConsecutiveDelimitersRule','join');  
f2 = readmatrix(fname2,'ConsecutiveDelimitersRule','join');  
del_ci2=3.42586246;del_ci3=3.40021219 ;
figure(1),clf
plot(y/del_ci2,f1(:,6),'-b','LineWidth',2,'DisplayName','c_{I}=0.03');hold on
plot(y/del_ci3,f2(:,6),'--r','LineWidth',2,'DisplayName','c_{I}=0.003');hold on
%plot(y/del_ci2,f1(:,6),'-b','LineWidth',2,'HandleVisibility','off');hold on
%plot(y/del_ci3,f2(:,6),'--r','LineWidth',2,'HandleVisibility','off');hold on
xlabel('y/{\delta_\theta(\tau)}');
xlim([-10 10])
h=ylabel('$\sqrt{|R_{12}|}/{\Delta{u}}$');set(h,'Interpreter','latex')
legend show
ax=gca;  ax.YRuler.Exponent = 0;
grid on;hAx=gca;set(hAx,'xminorgrid','on','yminorgrid','on')
screen2jpeg(strcat(outputdir1,'R12.png'));
