clear all, clf,clc

inputdir = "/home/vsj/Codes/vortexring/";
outputdir = "/home/vsj/Codes/vortexring/radial_plot/";
mkdir(outputdir)

caseid = ["pr10_285mm_256","pr10_285mm_512"];
caseidname = ["256x256x256","512x512x512"];
%caseid = ["ctest","cvltest"];
%caseidname = ["cgrid.F90","cvlgrid.F90"];

myg = [0.4660 0.6740 0.1880]; myb = [0.3010 0.7450 0.9330]; myv = [0.4940 0.1840 0.5560]; myy = [0.9290 0.6940 0.1250];
col={'b','r',myg,myy,'c'};
all_marks = {'o','+','v','s','p'};

counter = 520;                
for m=1:size(counter,2)
figure(1),clf
for p=1:size(caseid,2)-1
fname_coords=strcat(inputdir,caseid(p),'/vortexring_coords.h5');
x=h5read(fname_coords,'//X');y=h5read(fname_coords,'//Y');z=h5read(fname_coords,'//Z');
Lx=x(end,1,1);  Ly=2*y(1,end,1);    Lz=2*z(1,1,end);
[nx,ny,nz] = size(x);
dxi = Lx/(nx-1);  deta = Ly/(ny-1);  dzeta= Lz/(nz-1); 
x1 = 0;     xn = Lx;
y1 = -Ly/2; yn = Ly/2;
z1 = -Lz/2; zn = Lz/2;
for i=1:nx
    for j=1:ny
        for k=1:nz
            xi(i,j,k)    = x1 + ((i-1)*dxi); 
            eta(i,j,k)   = y1 + ((j-1)*deta); 
            zeta(i,j,k)  = z1 + ((k-1)*dzeta); 
        end
    end
end

%Find metric multipliers to find derivatives in physical space
dxidxij = get_metrics(x,y,z,dxi,deta,dzeta);

file=strcat(inputdir,caseid(p),'/vortexring_',sprintf('%04d',counter(m)),'.h5');
time = round(h5readatt(file,'/','Time')*10^6) 
u=h5read(file,strcat('//u'));  v=h5read(file,strcat('//v'));   w=h5read(file,strcat('//w'));

dudxi   = ddxi_compact10_gen(u,dxi);
dudeta  = ddeta_compact10_sym(u,deta);
dudzeta = ddzeta_compact10_sym(u,dzeta);
dudx = dxidxij(:,:,:,1).*dudxi + dxidxij(:,:,:,4).*dudeta + dxidxij(:,:,:,7).*dudzeta;
dudy = dxidxij(:,:,:,2).*dudxi + dxidxij(:,:,:,5).*dudeta + dxidxij(:,:,:,8).*dudzeta;
dudz = dxidxij(:,:,:,3).*dudxi + dxidxij(:,:,:,6).*dudeta + dxidxij(:,:,:,9).*dudzeta;

dvdxi   = ddxi_compact10_gen(v,dxi);
dvdeta  = ddeta_compact10_sym(v,deta);
dvdzeta = ddzeta_compact10_sym(v,dzeta);
dvdx = dxidxij(:,:,:,1).*dvdxi + dxidxij(:,:,:,4).*dvdeta + dxidxij(:,:,:,7).*dvdzeta;
dvdy = dxidxij(:,:,:,2).*dvdxi + dxidxij(:,:,:,5).*dvdeta + dxidxij(:,:,:,8).*dvdzeta;
dvdz = dxidxij(:,:,:,3).*dvdxi + dxidxij(:,:,:,6).*dvdeta + dxidxij(:,:,:,9).*dvdzeta;

dwdxi   = ddxi_compact10_gen(w,dxi);
dwdeta  = ddeta_compact10_sym(w,deta);
dwdzeta = ddzeta_compact10_sym(w,dzeta);
dwdx = dxidxij(:,:,:,1).*dwdxi + dxidxij(:,:,:,4).*dwdeta + dxidxij(:,:,:,7).*dwdzeta;
dwdy = dxidxij(:,:,:,2).*dwdxi + dxidxij(:,:,:,5).*dwdeta + dxidxij(:,:,:,8).*dwdzeta;
dwdz = dxidxij(:,:,:,3).*dwdxi + dxidxij(:,:,:,6).*dwdeta + dxidxij(:,:,:,9).*dwdzeta;

ome_x=  0.5*(dwdy-dvdz);
ome_y=  0.5*(dudz-dwdx);
ome_z=  0.5*(dvdx-dudy);

mag_omega = sqrt(ome_x.^2 + ome_y.^2 + ome_z.^2);
u_draw_xy =  squeeze(mag_omega(:,:,round(nz/2)));
xmax1 = max(u_draw_xy(5:end,150:end),[],'all')
[idx_ome,idy_ome]=find(abs(u_draw_xy-xmax1)<10^-5);

idx1 = 37;idx2=156 ; %pr10_165mm_256_t=520
idx3 = 38;idx4=99 ;  %pr10_165mm_256_t=520
u_y1 = v(:,idx2,round(nz/2));
u_y2 = -v(:,idx4,round(nz/2));
u_y0 = (u_y1+u_y2)/2;

idx_x1 = 29; idx_x2 = 46; idx_x3 = 48; idx_x4 = 51;


subplot(2,7,[1,2,3,8,9,10])
n=256;
t1=linspace(min(min(u_draw_xy)),max(max(u_draw_xy)),n+1);
t1=t1(1:n);
xstr=x(:,1,1)';ystr=y(1,:,1);zstr=squeeze(z(1,1,:))';
contourf(xstr, ystr, u_draw_xy','LineColor','None','LevelList',t1);hold on
xlab=xlabel('\boldmath$x(m)$');set(xlab,'Fontsize',10,'interpreter','latex');
ylab=ylabel('\boldmath$y(m)$');set(ylab,'Fontsize',10,'interpreter','latex');
colormap(jet(n));
c=colorbar;c.FontSize=8;c.Location='northoutside';
xlim([0 0.2]);ylim([-0.1 0.1])
xline(xstr(idx_x1),'--k'); xline(xstr(idx_x2),'--k'); xline(xstr(idx_x3),'--k'); xline(xstr(idx_x4),'--k');
yline(ystr(idx2(1)),'--k'); yline(ystr(idx4(1)),'--k');

subplot(2,7,[5,6,7,12,13,14])
frad = strcat(inputdir,'radial_pr10_285_520.csv');
radial_vel = readmatrix(frad,'ConsecutiveDelimiters','join');
txt=strcat('PIV (t=520\mu{s})');
plot(radial_vel(:,1)-0.0072,radial_vel(:,2),'ok','LineWidth',1,'DisplayName',txt,'Markerfacecolor','k','Markersize',5);hold on

txt=strcat(caseidname(p));
plot(xstr,u_y0,'Marker',all_marks{p},'color',col{p},'LineWidth',1, 'Linestyle','None','Displayname',txt,'Markersize',5);hold on
xlim([0 0.15]);
xline(xstr(idx_x1),'--k'); xline(xstr(idx_x2),'--k'); xline(xstr(idx_x3),'--k'); xline(xstr(idx_x4),'--k');
legend show;leg=legend('Location','northwest','box','off','Fontsize',8);
leg.ItemTokenSize = [3 3];
xlab=xlabel('\boldmath$x(m)$');set(xlab,'Fontsize',10,'interpreter','latex');
ylab=ylabel('\boldmath$v(m/s)$');set(ylab,'Fontsize',10,'interpreter','latex');
title(strcat('\boldmath$v(y=',sprintf('%1.3f',ystr(idx2(1))),'m):t=',sprintf('%04d',round((time))),'\mu{s}$'),'Fontsize',10,'interpreter','latex');
end

file_png= strcat(outputdir,'radial_',sprintf('%04d',round((time))),'.png');
screen2jpeg(file_png)
end
