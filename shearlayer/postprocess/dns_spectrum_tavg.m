clear,clf
mainid = "dns_data_csl";
runid  = "mc2p0";
caseid = ["dns"];
caseid_name = ["DNS"];

time_step =10;
t_i=[50];t_f=[110] 

outputdir1 =strcat('/home/vsj/Codes/plot/myplots/',mainid,'/spectrum_fil/');
mkdir(outputdir1)
c=sqrt(1.4);Mc=2.0;du=Mc*2*c;rho_0=1;Re_theta=1000;delta_0=1;
delta_0=1;      %delta_0=Re_theta/(rho_0*Re*du);
myg = [0.4660 0.6740 0.1880]; myb = [0.3010 0.7450 0.9330]; myv = [0.4940 0.1840 0.5560]; myy = [0.9290 0.6940 0.1250];

col={myv,myy,myg,'r','b'};
col1={'k','k',myg,'r','b'};
all_marks = {'p','s','d','^','v'};


figure(1),clf
for k=1:size(caseid,2)
inputdir =strcat('/home/vsj/',mainid,'/',runid,'/','my_fil/');
outputdir=strcat('/home/vsj/',mainid,'/',runid,'/postprocess_fil/');
outputdirfil=strcat('/home/vsj/',mainid,'/',runid,'/postprocess_fil/');

%%------------- Read Coordinates x, y, z-----------------------------
%fname_coords=strcat(inputdir,'shearlayer1mat_coords.h5');
Lx = 80; Ly = 80; Lz = 40;
%nx = 1024; ny = 1448; nz = 512;
nx = 171 ; ny = 242; nz = 86;
dx = Lx/nx; x = linspace(0,Lx-dx,nx); 
dy = Ly/ny; y = linspace(-Ly/2,Ly/2,ny); 
dz = Lz/nz; z = linspace(0,Lz-dz,nz); 
kx = fftshift((-nx/2:nx/2-1)*(2*pi/Lx));
ky = fftshift((-ny/2:ny/2-1)*(2*pi/Ly));
kz = fftshift((-nz/2:nz/2-1)*(2*pi/Lz));

fname = strcat(outputdir,'delta_',sprintf('%1.1f',Mc),'_',caseid(k),'.dat');
delta1= readmatrix(fname);t1=delta1(:,1);d1=delta1(:,2);%d1 is delta_theta

fname = strcat(outputdirfil,'delta_',sprintf('%1.1f',Mc),'_',caseid(k),'.dat');
delta2= readmatrix(fname);y1=delta2(:,8); t1_new=delta2(:,1); %y1=yc
%time,delta_theta,delta_dot,delta_omega,D_w,delta_99,delta_y,y_c,U_del/du);

time=t_i(k):time_step:t_f(k);

Ek2_zavg = zeros(1,round(nx/2+1));
Ek_time = zeros(size(time,2),size(Ek2_zavg,2));
del_time=zeros(size(time,2),1);

size(Ek_time)
for i=1:size(time,2)
t=time(i);
counter=t/time_step;                  %% Change if necessary--------
file=strcat(inputdir,'sl_fil',sprintf('%04d',counter),'.h5');
u=h5read(file,'//u');

idx=find(t1==t)
idx_new = find(t1_new==t)
del=d1(idx);  y_c=y1(idx_new);
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

Ek2_zavg = zeros(1,round(nx/2+1));
Ek2_zavg(1) = Ek_zavg(1);
for ii=2:nx/2+1
  Ek2_zavg(ii) = Ek_zavg(ii) + Ek_zavg(nx-ii+2);
end

Ek_time(i,:) = Ek2_zavg;
del_time(i) = del;

end
Ek_tavg = mean(Ek_time);
del_tavg = mean(del_time);

fname = strcat(outputdir1,'spectrum_tavg_',sprintf('%1.1f',Mc),'_',caseid(k),'.dat');
fileID = fopen(fname,'w');   
for p=1:size(Ek_time,2)
fprintf(fileID,'%3.5e\r\n',Ek_tavg(p));
end
fclose(fileID);


loglog(kx(1:round(nx/2+1)).*del_tavg,Ek_tavg,'Color',col{k},'Marker',all_marks{k},'Markerindices',1:5:length(kx(1:nx/2+1)),'LineWidth',1.5,'DisplayName',caseid_name(k));
hold on


if(k==size(caseid,2))
x_axis = kx(1:nx/2+1).*del_tavg;
y_axis = x_axis.^(-5/3);
xlim([min(abs(kx(1:nx/2+1)*del_tavg)) max(abs(kx(1:nx/2+1)*del_tavg))]);
loglog(x_axis,y_axis*50,'--k','LineWidth',2,'DisplayName','Slope=-5/3');hold on
ylim([min(Ek_tavg) max(y_axis(2:nx/2+1)*50)]);
end


end

grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
legend show;legend('Location','southwest','box','off','Fontname','Times','Fontsize',17);
ax=gca;ax.YAxis.FontSize = 20;ax.XAxis.FontSize = 20;
xlab=xlabel('\boldmath$k_{x}\delta_{\theta}$','Interpreter','latex');set(xlab,'Fontsize',22);
ylab=ylabel('\boldmath$E_{k}/\Delta{u^2}\delta_\theta$', 'Interpreter','latex');set(ylab,'Fontsize',22);
set(findall(gcf,'-property','Markersize'),'Markersize',5)
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
set(findall(gcf,'-property','FontName'),'FontName','Timesnewroman');
screen2jpeg(strcat(outputdir1,'spectrum_',sprintf('%1.1f',Mc),'_.jpeg'));

