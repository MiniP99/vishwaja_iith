clear
inputdir='/home/vsj/Codes/shearlayer_les/check/';
dir=strcat('/home/vsj/Codes/shearlayer_les/check/postprocess/');
dir1=strcat('/home/vsj/Codes/shearlayer_les/runs180/mgm_dyn_data/');
outputdir='/home/vsj/Codes/shearlayer_les/check/contours/';
mkdir(outputdir)
caseid="mgm_dyn";          %% change mc, domain

%%------------- Read Coordinates x, y, z-----------------------------
fname_coords=strcat(inputdir,'shearlayer1mat_coords.h5');
X=h5read(fname_coords,'//X');Y=h5read(fname_coords,'//Y');Z=h5read(fname_coords,'//Z');
x=X(:,1,1)';y=Y(1,:,1);z=squeeze(Z(1,1,:))';
Lx=x(end);     Ly=2*y(end);  Lz=z(end);
nx=size(x,2);  ny=size(y,2); nz=size(z,2);
dx = Lx/nx;    dy= Ly/ny;    dz= Lz/nz; 
kx = fftshift((-nx/2:nx/2-1)*(2*pi/Lx));
ky = fftshift((-ny/2:ny/2-1)*(2*pi/Ly));
kz = fftshift((-nz/2:nz/2-1)*(2*pi/Lz));
[kx3D,ky3D,kz3D] = ndgrid(kx,ky,kz);

% inputs
c=sqrt(1.4);Mc=2.0;du=Mc*2*c;rho_0=1;Re_theta=1000;
%delta_0=Re_theta/(rho_0*Re*du);
delta_0=1;



clf
time=0:250;

file_video='/home/vsj/Codes/shearlayer_les/check/contours/uvel.mp4';
obj=VideoWriter(file_video);
obj.Quality=100;
obj.FrameRate=6;
open(obj);

ll=0;

for i=1:size(time,2)
counter=time(i)/1;
tau = time(i)*du/delta_0;

file=strcat(inputdir,'shearlayer1mat_',sprintf('%04d',counter),'.h5')
u=h5read(file,'//u');v=h5read(file,'//v');w=h5read(file,'//w');
rho=h5read(file,'//rho');p=h5read(file,'//p');
mu=h5read(file,'//mu');
bulk=h5read(file,'//bulk');
T=h5read(file,'//T');kap=h5read(file,'//kap');
t11sgs=h5read(file,'//t11');t12sgs=h5read(file,'//t12');t13sgs=h5read(file,'//t13');
t22sgs=h5read(file,'//t22');t23sgs=h5read(file,'//t23');t33sgs=h5read(file,'//t33');

%calculating derivatives-----------------------------
dudx=ddx_compact(u,dx,nx,ny,nz);dudy=ddy_compact(u,dy,nx,ny,nz);dudz=ddz_compact(u,dz,nx,ny,nz);
dvdx=ddx_compact(v,dx,nx,ny,nz);dvdy=ddy_compact(v,dy,nx,ny,nz);dvdz=ddz_compact(v,dz,nx,ny,nz);
dwdx=ddx_compact(w,dx,nx,ny,nz);dwdy=ddy_compact(w,dy,nx,ny,nz);dwdz=ddz_compact(w,dz,nx,ny,nz);
% calculating averages,fluctuations------------------
rho_bar=avgxz(rho,nx,nz);
u_bar=avgxz(u,nx,nz);[u_tilde,u_pp]=favre_avg_fluct(rho,u,nx,ny,nz);
v_bar=avgxz(v,nx,nz);[v_tilde,v_pp]=favre_avg_fluct(rho,v,nx,ny,nz);
w_bar=avgxz(w,nx,nz);[w_tilde,w_pp]=favre_avg_fluct(rho,w,nx,ny,nz);
% tausgs average
t11sgs_bar=avgxz(t11sgs,nx,nz);t12sgs_bar=avgxz(t12sgs,nx,nz);t13sgs_bar=avgxz(t13sgs,nx,nz);
t22sgs_bar=avgxz(t22sgs,nx,nz);t23sgs_bar=avgxz(t23sgs,nx,nz);t33sgs_bar=avgxz(t33sgs,nx,nz);
% calculating derivatives of favre_avg
du_tildx=zeros(1,ny);du_tildy=ddy1D_compact(u_tilde,dy,nx,ny,nz);du_tildz=zeros(1,ny);
dv_tildx=zeros(1,ny);dv_tildy=ddy1D_compact(v_tilde,dy,nx,ny,nz);dv_tildz=zeros(1,ny);
dw_tildx=zeros(1,ny);dw_tildy=ddy1D_compact(w_tilde,dy,nx,ny,nz);dw_tildz=zeros(1,ny);
% calculating  derivatives of favre_fluct--------------
du_ppdx=ddx_compact(u_pp,dx,nx,ny,nz);du_ppdy=ddy_compact(u_pp,dy,nx,ny,nz);du_ppdz=ddz_compact(u_pp,dz,nx,ny,nz);
dv_ppdx=ddx_compact(v_pp,dx,nx,ny,nz);dv_ppdy=ddy_compact(v_pp,dy,nx,ny,nz);dv_ppdz=ddz_compact(v_pp,dz,nx,ny,nz);
dw_ppdx=ddx_compact(w_pp,dx,nx,ny,nz);dw_ppdy=ddy_compact(w_pp,dy,nx,ny,nz);dw_ppdz=ddz_compact(w_pp,dz,nx,ny,nz);
%--------------------------------- calculating momentum thickness------------------------------------
[R11,R12,R13,R22,R23,R33]=reynold_stress(rho,u_pp,v_pp,w_pp,nx,nz);
R11 = rho_bar.*R11;R12 = rho_bar.*R12;R13 = rho_bar.*R13;R22 = rho_bar.*R22;R23 = rho_bar.*R23;R33 = rho_bar.*R33;
R11=R11-t11sgs_bar;R12=R12-t12sgs_bar;R13=R13-t13sgs_bar;R22=R22-t22sgs_bar;R23=R23-t23sgs_bar;R33=R33-t33sgs_bar;

func1=rho_bar.*(0.5*du-u_tilde).*(0.5*du+u_tilde);
delta_theta=(1/(rho_0*(du^2)))*sum(func1)*dy;

%Dissipation------------------------------
[tau11,tau12,tau13,tau22,tau23,tau33]=shearstress(mu,bulk,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,t11sgs,t12sgs,t13sgs,t22sgs,t23sgs,t33sgs);               %%tausgs included
[tau11_bar,tau11_rey_fluct,tau12_bar,tau12_rey_fluct,tau13_bar,tau13_rey_fluct,tau22_bar,tau22_rey_fluct,tau23_bar,tau23_rey_fluct,tau33_bar,tau33_rey_fluct]=tau_avg_fluct(tau11,tau12,tau13,tau22,tau23,tau33,nx,ny,nz);

%func_tau=tau11.*du_ppdx+tau12.*du_ppdy+tau13.*du_ppdz+tau12.*dv_ppdx+tau22.*dv_ppdy+tau23.*dv_ppdz+tau13.*dw_ppdx+tau23.*dw_ppdy+tau33.*dw_ppdz;
func_tau=tau11_rey_fluct.*du_ppdx+tau12_rey_fluct.*du_ppdy+tau13_rey_fluct.*du_ppdz+tau12_rey_fluct.*dv_ppdx+tau22_rey_fluct.*dv_ppdy+tau23_rey_fluct.*dv_ppdz+tau13_rey_fluct.*dw_ppdx+tau23_rey_fluct.*dw_ppdy+tau33_rey_fluct.*dw_ppdz;

diss=avgxz(func_tau,nx,nz);

%Integrated Reynold stresses
R11_int=(sum(R11)*dy)/(rho_0*(du^2)*delta_theta);R12_int=(sum(R12)*dy)/(rho_0*(du^2)*delta_theta);
R13_int=(sum(R13)*dy)/(rho_0*(du^2)*delta_theta);R22_int=(sum(R22)*dy)/(rho_0*(du^2)*delta_theta);
diss_int = (sum(diss)*dy)/(du^3);

u_draw=squeeze(u(:,:,nz/2))';
n=256;
t1=linspace(min(min(u_draw)),max(max(u_draw)),n+1);
t1=t1(1:n);


figure(1)
subplot(3,3,[1,2,4,5])
b=contourf(x, y, u_draw,'LineColor','None','LevelList',t1);
zmin=min(u_draw,[],'all')
zmax=max(u_draw,[],'all')
pbaspect([2 2 1])
ax=gca;ax.YAxis.FontSize = 6;ax.XAxis.FontSize = 6;
%tit=title(strcat('\boldmath$\tau=t\Delta{u}/\delta_{\theta,0}=',sprintf('%4.0f',tau),'$'));set(tit,'Fontsize',8,'interpreter','latex');
tit=title(strcat('\boldmath$\tau=',sprintf('%4.0f',tau),'$'));set(tit,'Fontsize',8,'interpreter','latex');
xlab=xlabel('\boldmath$x/\delta_{\theta,0}$');set(xlab,'Fontsize',8,'interpreter','latex');
ylab=ylabel('\boldmath$y/\delta_{\theta,0}$');set(ylab,'Fontsize',8,'interpreter','latex');
c=colorbar;c.FontSize=6;colormap(jet(n));caxis([-2.3664 2.3664]);c.Ticks=([-3,-2,-1,0,1,2,3]);
c.Position(1) = c.Position(1)+0.05;c.Position(3) = c.Position(3)-0.02;

subplot(3,3,[9])
plot(tau,delta_theta,'-or','HandleVisibility','off');hold on
if(i==size(time,2))
fname = strcat(dir,'delta_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
f=readmatrix(fname);
idx1=find(f(:,1)==110);idx2=find(f(:,1)==250);idx3=find(f(:,1)==250);
tau1=(f(idx1:idx2,1)*du)/delta_0;
del1=f(idx1:idx2,2)/delta_0;
eq=polyfit(tau1,del1,1);
delta1 = polyval(eq,tau1);
plot(tau1,delta1,':k','LineWidth',1,'Displayname',strcat('\boldmath$\dot{\delta}_{\theta}=',sprintf('%1.3f',eq(1,1)),'$'));
end
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
leg=legend('Location','northwest','box','off','Fontsize',6,'interpreter','latex');legend show
leg.ItemTokenSize = [10 10];
ax=gca;ax.YAxis.FontSize = 6;ax.XAxis.FontSize = 6;
xlab=xlabel('\boldmath$\tau$');set(xlab,'Fontsize',8,'interpreter','latex');
ylab=ylabel('\boldmath$\delta_{\theta}/ \delta_{\theta,0}$');set(ylab,'Fontsize',8,'interpreter','latex');
xticks([0 400 800 1200])
xlim([0 1200]);ylim([1 6]);


subplot(3,3,[7,8])
yyaxis left
if(i==1)
plot(tau,R11_int,'ob','DisplayName','\boldmath$R_{11,int}$');hold on
plot(tau,-R12_int,'<m','DisplayName','\boldmath$R_{12,int}$');hold on
plot(tau,R22_int,'^g','DisplayName','\boldmath$R_{22,int}$');hold on
else
plot(tau,R11_int,'ob','HandleVisibility','off');hold on
plot(tau,-R12_int,'<m','HandleVisibility','off');hold on
plot(tau,R22_int,'^g','HandleVisibility','off');hold on
end
ax.YColor='k';
ax=gca;ax.YAxis(1).Exponent = -2;
ax=gca;ax.YAxis(1).FontSize = 6;ax.XAxis.FontSize = 6;ax.YAxis(1).Color = 'k';
xlim([0 1200]);ylim([0 0.09]);
ylab=ylabel('\boldmath$$R_{ij,int}$');set(ylab,'Fontsize',8,'interpreter','latex');

yyaxis right
if(i==1)
plot(tau,diss_int,'>r','DisplayName','\boldmath$\epsilon_{int}$');hold on;ax=gca;ax.YAxis(2).Exponent = 0;
else
plot(tau,diss_int,'>r','HandleVisibility','off');hold on;ax=gca;ax.YAxis(2).Exponent = 0;
end
if(i==size(time,2))
p=patch([tau1(1) tau1(1) tau1(end) tau1(end)],[min(ylim) max(ylim) max(ylim) min(ylim)],[0.9 0.9 0.9],'EdgeColor','None','HandleVisibility','off');alpha(p,0.5)
end
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
ax.YColor='r';
ax=gca;ax.YAxis(2).FontSize = 6;ax.XAxis.FontSize = 6;
ax=gca;ax.YAxis(2).Exponent = -3;
%ylab=ylabel('\boldmath$\epsilon_{int}$');set(ylab,'Fontsize',10,'interpreter','latex');
xlim([0 1200]);ylim([0 0.003]);
xlab=xlabel('\boldmath$\tau$');set(xlab,'Fontsize',8,'interpreter','latex');
leg=legend('Location','northwest','box','off','Fontsize',7.5,'interpreter','latex');legend show
leg.ItemTokenSize = [5 5];


subplot(3,3,[3,6])
all_marks = {'o','+','*','x','d','^','v','>','<','p','h'};

if(rem(time(i),25)==0)
ll=ll+1;
p=plot(y/delta_theta,sqrt(abs(R11))/du,'Linewidth',1,'DisplayName',strcat('\boldmath$',sprintf('%4.0f',tau),'$'));hold on
txt1 = '\boldmath$\tau$';
t = text(-11,0.135,txt1,'Fontsize',7.5,'FontWeight','bold','Interpreter','Latex');

if(time(i)<109)
p.Color(4)=0.5
end

if(i==size(time,2))
%frnam = strcat(dir1,'vel_tke_flux_tavg_nor_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
%fr = readmatrix(frnam,'ConsecutiveDelimitersRule','join'); 
%plot(y/delta_theta,fr(:,5),'k','LineStyle','-','LineWidth',1,'HandleVisibility','off');hold on
end

end

grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
ax=gca;ax.YAxis.FontSize = 6;ax.XAxis.FontSize = 6;
ax=gca;ax.YAxis.Exponent = -2;
xlab=xlabel('\boldmath$y/\delta_{\theta}$');set(xlab,'Fontsize',8,'interpreter','latex');
ylab=ylabel('\boldmath$\sqrt{R_{11}}/\Delta{u}$');set(ylab,'Fontsize',8,'interpreter','latex');
leg=legend('Location','northwest','box','off','Fontsize',6,'interpreter','latex');legend show
leg.ItemTokenSize = [5 5];
xlim([-15 15]);ylim([0 0.14]);

%%-----final-------------
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
set(findall(gcf,'-property','Markersize'),'Markersize',1)
%set(findall(gcf,'-property','Linestyle'),'Linestyle','-')
file_png= strcat(outputdir,'u_contxy_',sprintf('%04d',time(i)),'_',sprintf('%1.1f',Mc),'_',caseid,'.png');
screen2jpeg(file_png)

%%-----create video-----
f= imread(file_png);
writeVideo(obj,f);
end
obj.close();
