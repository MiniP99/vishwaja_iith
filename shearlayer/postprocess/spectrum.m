clear,clf
mainid = "runs180";
runid  = "run20";
caseid = ["base","lad","amd_dyn","sigma_dyn","mgm_dyn"];
caseid_name = ["Base","LAD","AMD-Dyn","Sigma-Dyn","MGM-Dyn"];
%t_i=[120,160,150,150,120];t_f=[300,340,330,330,300] %Mc=1.2
t_i=[0,110,100,130,110];t_f=[0,250,240,280,250] %Mc=2.0


outputdir1 =strcat('/home/vsj/Codes/plot/myplots/',mainid,'/');
c=sqrt(1.4);Mc=2.0;du=Mc*2*c;rho_0=1;Re_theta=1000;delta_0=1;
delta_0=1;      %delta_0=Re_theta/(rho_0*Re*du);
myg = [0.4660 0.6740 0.1880]; myb = [0.3010 0.7450 0.9330]; myv = [0.4940 0.1840 0.5560]; myy = [0.9290 0.6940 0.1250];
col={myv,myy,myg,'r','b'};
col1={'k','k',myg,'r','b'};
all_marks = {'p','s','d','^','v'};
figure(1),clf
for k=2:size(caseid,2)
inputdir =strcat('/home/vsj/Codes/',mainid,'/',runid,'/',caseid(k),'/');
outputdir=strcat('/home/vsj/Codes/',mainid,'/',runid,'/',caseid(k),'/postprocess_new/');

%%------------- Read Coordinates x, y, z-----------------------------
fname_coords=strcat(inputdir,'shearlayer1mat_coords.h5');
X=h5read(fname_coords,'//X');Y=h5read(fname_coords,'//Y');Z=h5read(fname_coords,'//Z');
x=X(:,1,1)';y=Y(1,:,1);z=squeeze(Z(1,1,:))';
Lx=x(end);Ly=2*y(end);Lz=z(end);
nx=size(x,2);ny=size(y,2);nz=size(z,2);
dx = Lx/nx; dy = Ly/ny; dz = Lz/nz; 
kx = fftshift((-nx/2:nx/2-1)*(2*pi/Lx));
ky = fftshift((-ny/2:ny/2-1)*(2*pi/Ly));
kz = fftshift((-nz/2:nz/2-1)*(2*pi/Lz));

fname = strcat(outputdir,'delta_',sprintf('%1.1f',Mc),'_',caseid(k),'.dat');
delta1= readmatrix(fname);t1=delta1(:,1);d1=delta1(:,2);d2=delta1(:,4);y1=delta1(:,9);
fname1 = strcat(outputdir,'Rey_kol_',sprintf('%1.1f',Mc),'_',caseid(k),'.dat');
rey_kol= readmatrix(fname);kol=rey_kol(:,5);


for i=40
counter=i;
time=counter*5
file=strcat(inputdir,'shearlayer1mat_',sprintf('%04d',counter),'.h5')
u=h5read(file,'//u');

idx=find(t1==time);    
del=d1(idx);y_c=y1(idx);
eeta=kol(idx);
[val_yc,idx_yc]=min(abs(y_c-y));

Ek_all = zeros(nx,nz);

for j=1:nz
u_yc= u(:,idx_yc,j);
uh = fft(u_yc);
Ek = uh.*conj(uh);%Ek = (0.5/(nx*ny*nz))*Ek;
Ek_all(:,j) = Ek;
end
Ek_zavg = mean(Ek_all');
Ek_zavg = Ek_zavg./((du^2).*del);

Ek2_zavg = zeros(1,nx/2+1);
Ek2_zavg(1) = Ek_zavg(1);
for ii=2:nx/2+1
  Ek2_zavg(ii) = Ek_zavg(ii) + Ek_zavg(nx-ii+2);
end

%loglog(kx.*del,Ek_zavg,'o-', kx(1:nx/2+1)*del, Ek2_zavg, '-x')
loglog(kx(1:nx/2+1).*del,Ek2_zavg,'Color',col{k},'Marker',all_marks{k},'Markerindices',1:5:length(kx(1:nx/2+1)),'LineWidth',1.5,'DisplayName',caseid_name(k));
%loglog(kx(1:nx/2+1).*del,Ek2_zavg./kx(1:nx/2+1).^(-5/3),'Color',col{k},'Marker',all_marks{k},'Markerindices',1:5:length(x),'LineWidth',1,'DisplayName',caseid_name(k));
%loglog(kx.*del,Ek_zavg,'Color',col{k},'LineWidth',1,'DisplayName',caseid_name(k));
hold on
end
end
x_axis = kx(1:nx/2+1).*del;
y_axis = x_axis.^(-5/3);
xlim([min(abs(kx(1:nx/2+1)*del)) max(abs(kx(1:nx/2+1)*del))]);
loglog(x_axis,y_axis*50,'--k','LineWidth',2,'DisplayName','Slope=-5/3');hold on
ylim([min(Ek2_zavg) max(y_axis(2:nx/2+1)*50)]);
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
legend show;legend('Location','southwest','box','off','Fontname','Timesnewroman','Fontsize',13);
ax=gca;ax.YAxis.FontSize = 13;ax.XAxis.FontSize = 13;
xlab=xlabel('\boldmath$k_{x}\delta_{\theta}$','Interpreter','latex');set(xlab,'Fontsize',17);
ylab=ylabel('\boldmath$E_{k}/\Delta{u^2}\delta_\theta$', 'Interpreter','latex');set(ylab,'Fontsize',17);
set(findall(gcf,'-property','Markersize'),'Markersize',6)
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
screen2jpeg(strcat(outputdir1,'spectrum_',sprintf('%1.1f',Mc),'_.jpeg'));

