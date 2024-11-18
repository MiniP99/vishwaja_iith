clear,clf
inputdir='/home/vsj/Codes/shearlayer_les/runs180/';
inputdir1='/home/vsj/Codes/shearlayer_les/runs180/smag_dyn_data/';
outputdir='/home/vsj/Codes/shearlayer_les/plot/myplots/runs180/predil_plots/smag_dyn_180/';
mkdir(outputdir);
caseid='smag_dyn';
figure(1),clf
myg = [0.4660 0.6740 0.1880]; myb = [0.3010 0.7450 0.9330]; myv = [0.4940 0.1840 0.5560]; myy = [0.9290 0.6940 0.1250];
col={'r',myg,'b',myv,myy,'k'};
all_marks = {'o','v','d','^','s','>','x','v','>','<','p','h'};
%mc=["02","04","08","12","16","20"];
%Mc=[0.2,0.4,0.8,1.2,1.6,2.0];
%Lx_mat=[150,150,100,100,80,80];Ly_mat=[200,200,100,100,80,80];Lz_mat=[75,75,50,50,40,40];

mc=["02","08","12","16","20"];
Mc=[0.2,0.8,1.2,1.6,2.0];
Lx_mat=[150,100,100,80,80];Ly_mat=[200,100,100,80,80];Lz_mat=[75,50,50,40,40];

du=2*Mc*sqrt(1.4);
nx=180;ny=240;nz=84;
%nx=128;ny=180;nz=64;
%nx=256;ny=360;nz=128;

t_i=[560,300,210,140,130];t_f=[1000,500,380,300,280]; %Smag_dyn
%t_i=[660,350,300,0,450];t_f=[1000,550,480,0,600];  %Smag_sta
time_step = [20,5,5,5,5];


%t_i = [400,280,150,120,0,0];t_f=[900,600,350,300,0,0];    %base
%t_i = [400,280,150,160,140,110];t_f=[900,600,350,340,280,250];    %lad
%t_i = [560,380,200,150,110,100];t_f=[1000,800,400,330,270,240];   %AMD_dyn
%t_i = [560,400,200,150,110,100];t_f=[900,800,400,330,250,240];   %AMD_sta
%t_i = [560,480,200,150,110,130];t_f=[1000,800,400,330,270,280];   %sigma
%t_i = [460,380,150,120,110,110];t_f=[900,700,350,300,250,250];    %mgm
%t_i = [460,280,150,120,110,110];t_f=[900,600,350,300,250,250];    %mgm_sta
%time_step = [20,20,5,5,5,5];
%%% ---------------- Plot of averaged profiles along y at different Mc --------------------- %%
for i=2%1:2:size(mc,2)
step = t_i(i):time_step(i):t_f(i);
Lx=Lx_mat(i);dx = Lx/nx; x = linspace(-Lx/2,Lx/2,nx)'; 
Ly=Ly_mat(i);dy = Ly/ny; y = linspace(-Ly/2,Ly/2,ny)'; 
Lz=Lz_mat(i);dz = Lz/nz; z = linspace(-Lz/2,Lz/2,nz)'; 

fname = strcat(inputdir,'run',mc(i),'/',caseid,'/postprocess','/delta_',sprintf('%1.1f',Mc(i)),'_',caseid,'.dat');
%fname = strcat(inputdir,'/postprocess','/delta_',sprintf('%1.1f',Mc(i)),'_',caseid,'.dat');
delta1= readmatrix(fname);t1=delta1(:,1);d1=delta1(:,2);d2=delta1(:,4);


for j= 1: size(step,2)

counter=step(j)/time_step(i);
file=strcat(inputdir,'run',mc(i),'/',caseid,'/shearlayer1mat_',sprintf('%04d',counter),'.h5');
%file=strcat(inputdir,'shearlayer1mat_',sprintf('%04d',counter),'.h5');

idx=find(t1==step(j));    
del=d1(idx);

%----------------Read fields data-------------------------------------
u=h5read(file,'//u');v=h5read(file,'//v');w=h5read(file,'//w');
rho=h5read(file,'//rho');p=h5read(file,'//p');
mu=h5read(file,'//mu');T=h5read(file,'//T');kap=h5read(file,'//kap');bulk=h5read(file,'//bulk');
%t11sgs=h5read(file,'//t11');t12sgs=h5read(file,'//t12');t13sgs=h5read(file,'//t13');
%t22sgs=h5read(file,'//t22');t23sgs=h5read(file,'//t23');t33sgs=h5read(file,'//t33');
%q1sgs=h5read(file,'//q1');q2sgs=h5read(file,'//q2');q3sgs=h5read(file,'//q3');

% calculating averages,fluctuations------------------
rho_bar=avgxz(rho,nx,nz);
u_bar=avgxz(u,nx,nz);[u_tilde,u_pp]=favre_avg_fluct(rho,u,nx,ny,nz);
v_bar=avgxz(v,nx,nz);[v_tilde,v_pp]=favre_avg_fluct(rho,v,nx,ny,nz);
% calculating  derivatives of favre_fluct--------------
du_ppdx=ddx_compact(u_pp,dx,nx,ny,nz);du_ppdy=ddy_compact(u_pp,dy,nx,ny,nz);du_ppdz=ddz_compact(u_pp,dz,nx,ny,nz);
dv_ppdx=ddx_compact(v_pp,dx,nx,ny,nz);dv_ppdy=ddy_compact(v_pp,dy,nx,ny,nz);dv_ppdz=ddz_compact(v_pp,dz,nx,ny,nz);

[p_bar,p_reynold_fluct]    = reynolds_avg_fluct(p,nx,ny,nz);
[rho_bar,rho_reynold_fluct]= reynolds_avg_fluct(rho,nx,ny,nz);
[u_bar,u_reynold_fluct]    = reynolds_avg_fluct(u,nx,ny,nz);
[v_bar,v_reynold_fluct]    = reynolds_avg_fluct(v,nx,ny,nz);
[w_bar,w_reynold_fluct]    = reynolds_avg_fluct(w,nx,ny,nz);
%Pressure-strain terms----------------------
pi11 = avgxz(p_reynold_fluct.*2.*du_ppdx,nx,ny);         
pi12 = avgxz(p_reynold_fluct.*(du_ppdy+dv_ppdx),nx,ny);   
pi22 = avgxz(p_reynold_fluct.*2.*dv_ppdy,nx,ny);         

pi11_bar(j,:) = pi11.*(del/(du(i)^3));
pi12_bar(j,:) = pi12.*(del/(du(i)^3));
pi22_bar(j,:) = pi22.*(del/(du(i)^3));

end
%%% Time average
pi11_tavg= mean(pi11_bar);
pi12_tavg= mean(pi12_bar);
pi22_tavg= mean(pi22_bar);
size(pi11_bar)

fname = strcat(outputdir,'precorr_tavg_',sprintf('%1.1f',Mc(i)),'_',caseid,'.dat');
fileID = fopen(fname,'w');   
for p=1:size(pi11_tavg,2)
fprintf(fileID,'%3.5e %3.5e %3.5e\r\n',pi11_tavg(p),pi12_tavg(p),pi22_tavg(p));
end
fclose(fileID);

fname2 = strcat(inputdir1,'Rij_KE_P12_dis_pi_int_tavg_',caseid,'.dat');
f2 = readmatrix(fname2,'ConsecutiveDelimitersRule','join'); 
%%f2 = importdata(fname2); 
idx=find(abs(f2(:,1)-Mc(i)) < 0.001);del_tavg=f2(idx,2);
%%1-Mc,2-del_avg,3-R11,4-R12,5-R13,6-R22,7-R23,8-R33,9-ke,10-P12,11-diss,12-pi11,13-pi12,14-pi22
plot(y/del_tavg,pi11_tavg,'r','DisplayName','\pi 11','Linewidth',1.5);hold on
plot(y/del_tavg,pi12_tavg,'b','DisplayName','\pi 12','Linewidth',1.5);hold on
plot(y/del_tavg,pi22_tavg,'g','DisplayName','\pi 22','Linewidth',1.5);hold on
end

ax=gca;ax.YAxis.FontSize = 15;ax.XAxis.FontSize = 15;
xlab=xlabel('\boldmath$y/{\delta_\theta(\tau)}$','Interpreter','latex');set(xlab,'Fontsize',20);
ylab=ylabel('\boldmath$\pi \delta_\theta/\Delta{u}^3$','Interpreter','latex');set(ylab,'Fontsize',20);
xlim([-10 10])
legend show;legend('Location','Northeast','box','off','Fontsize',14,'Fontname','Timesnewroman');
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
set(findall(gcf,'-property','Markersize'),'Markersize',6)
%ax=gca;  ax.YRuler.Exponent = 0;
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
screen2jpeg(strcat(outputdir,'pi_',caseid,'.jpeg'));
