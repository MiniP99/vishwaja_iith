clear,clf
inputdir='/home/vsj/Codes/runs180/';
inputdir1='/home/vsj/Codes/runs180/mgm_dyn_data/';
outputdir='/home/vsj/Codes/plot/myplots/runs180/sgs_prod_plots/';
mkdir(outputdir);
caseid='mgm_dyn';
figure(1),clf
myg = [0.4660 0.6740 0.1880]; myb = [0.3010 0.7450 0.9330]; myv = [0.4940 0.1840 0.5560]; myy = [0.9290 0.6940 0.1250];
col={'r',myg,'b',myv,myy,'k'};
all_marks = {'o','v','d','^','s','>','x','v','>','<','p','h'};
mc=["02","04","08","12","16","20"];
Mc=[0.2,0.4,0.8,1.2,1.6,2.0];
Lx_mat=[150,150,100,100,80,80];Ly_mat=[200,200,100,100,80,80];Lz_mat=[75,75,50,50,40,40];
du=2*Mc*sqrt(1.4);
nx=180;ny=240;nz=84;


%t_i = [560,380,200,150,110,100];t_f=[1000,800,400,330,270,240];   %AMD
%t_i = [560,480,200,150,110,130];t_f=[1000,800,400,330,270,280];   %sigma
t_i = [460,380,150,120,110,110];t_f=[900,700,350,300,250,250];    %mgm
time_step = [20,20,5,5,5,5];
%%% ---------------- Plot of averaged profiles along y at different Mc --------------------- %%
for i=1:size(mc,2)
step = t_i(i):time_step(i):t_f(i);
Lx=Lx_mat(i);dx = Lx/nx; x = linspace(-Lx/2,Lx/2,nx)'; 
Ly=Ly_mat(i);dy = Ly/ny; y = linspace(-Ly/2,Ly/2,ny)'; 
Lz=Lz_mat(i);dz = Lz/nz; z = linspace(-Lz/2,Lz/2,nz)'; 

fname = strcat(inputdir,'run',mc(i),'/',caseid,'/postprocess','/delta_',sprintf('%1.1f',Mc(i)),'_',caseid,'.dat');
delta1= readmatrix(fname);t1=delta1(:,1);d1=delta1(:,2);d2=delta1(:,4);

diss_sgs_bar = zeros(size(step,2),ny);
diss_sgs_tavg = zeros(1,ny);

for j= 1: size(step,2)

counter=step(j)/time_step(i);
file=strcat(inputdir,'run',mc(i),'/',caseid,'/shearlayer1mat_',sprintf('%04d',counter),'.h5');

idx=find(t1==step(j));    
del=d1(idx);

%----------------Read fields data-------------------------------------
u=h5read(file,'//u');v=h5read(file,'//v');w=h5read(file,'//w');
rho=h5read(file,'//rho');p=h5read(file,'//p');
mu=h5read(file,'//mu');T=h5read(file,'//T');kap=h5read(file,'//kap');bulk=h5read(file,'//bulk');
t11sgs=h5read(file,'//t11');t12sgs=h5read(file,'//t12');t13sgs=h5read(file,'//t13');
t22sgs=h5read(file,'//t22');t23sgs=h5read(file,'//t23');t33sgs=h5read(file,'//t33');
q1sgs=h5read(file,'//q1');q2sgs=h5read(file,'//q2');q3sgs=h5read(file,'//q3');

% calculating averages,fluctuations------------------
%rho_bar=avgxz(rho,nx,nz);
%u_bar=avgxz(u,nx,nz);[u_tilde,u_pp]=favre_avg_fluct(rho,u,nx,ny,nz);
%v_bar=avgxz(v,nx,nz);[v_tilde,v_pp]=favre_avg_fluct(rho,v,nx,ny,nz);
%w_bar=avgxz(w,nx,nz);[w_tilde,w_pp]=favre_avg_fluct(rho,w,nx,ny,nz);
%T_bar=avgxz(T,nx,nz);[T_tilde,T_pp]=favre_avg_fluct(rho,T,nx,ny,nz);
%calculating derivatives-----------------------------
dudx=ddx_compact(u,dx,nx,ny,nz);dudy=ddy_compact(u,dy,nx,ny,nz);dudz=ddz_compact(u,dz,nx,ny,nz);
dvdx=ddx_compact(v,dx,nx,ny,nz);dvdy=ddy_compact(v,dy,nx,ny,nz);dvdz=ddz_compact(v,dz,nx,ny,nz);
dwdx=ddx_compact(w,dx,nx,ny,nz);dwdy=ddy_compact(w,dy,nx,ny,nz);dwdz=ddz_compact(w,dz,nx,ny,nz);
dTdx=ddx_compact(T,dx,nx,ny,nz);dTdy=ddy_compact(T,dy,nx,ny,nz);dTdz=ddz_compact(T,dz,nx,ny,nz);
% calculating  derivatives of favre_fluct--------------
%du_ppdx=ddx_compact(u_pp,dx,nx,ny,nz);du_ppdy=ddy_compact(u_pp,dy,nx,ny,nz);du_ppdz=ddz_compact(u_pp,dz,nx,ny,nz);
%dv_ppdx=ddx_compact(v_pp,dx,nx,ny,nz);dv_ppdy=ddy_compact(v_pp,dy,nx,ny,nz);dv_ppdz=ddz_compact(v_pp,dz,nx,ny,nz);
%dw_ppdx=ddx_compact(w_pp,dx,nx,ny,nz);dw_ppdy=ddy_compact(w_pp,dy,nx,ny,nz);dw_ppdz=ddz_compact(w_pp,dz,nx,ny,nz);
%dT_ppdx=ddx_compact(T_pp,dx,nx,ny,nz);dT_ppdy=ddy_compact(T_pp,dy,nx,ny,nz);dT_ppdz=ddz_compact(T_pp,dz,nx,ny,nz);



%Subgrid Production------------------------------
dx1 = (dx^2)/12;dy1 = (dy^2)/12;dz1 = (dz^2)/12;
G1T = dx1.*dudx.*dTdx +  dy1.*dudy.*dTdy + dz1.*dudz.*dTdz ;
G2T = dx1.*dvdx.*dTdx +  dy1.*dvdy.*dTdy + dz1.*dvdz.*dTdz ;
G3T = dx1.*dwdx.*dTdx +  dy1.*dwdy.*dTdy + dz1.*dwdz.*dTdz ;

prod_sgs_T= avgxz(-(G1T.*dTdx + G2T.*dTdy +G3T.*dTdz),nx,nz);
prod_sgs_bar_T(j,:) = prod_sgs_T.*(del/(du(i)^3));



%Subgrid Dissipation------------------------------

%func_tau=t11sgs.*du_ppdx+t12sgs.*du_ppdy+t13sgs.*du_ppdz+t12sgs.*dv_ppdx+t22sgs.*dv_ppdy+t23sgs.*dv_ppdz+t13sgs.*dw_ppdx+t23sgs.*dw_ppdy+t33sgs.*dw_ppdz;
%func_tau2=q1sgs.*dT_ppdx+q2sgs.*dT_ppdy+q3sgs.*dT_ppdz;
%diss_sgs=avgxz(func_tau,nx,nz);
%diss_sgs_bar(j,:) = diss_sgs.*(del/(du(i)^3));
%
%diss_sgs_T=avgxz(func_tau2,nx,nz);
%diss_sgs_bar_T(j,:) = diss_sgs_T.*(del/(du(i)^3));
end
%%% Time average
%diss_sgs_tavg= mean(diss_sgs_bar);
%size(diss_sgs_bar)

prod_sgs_tavg_T= mean(prod_sgs_bar_T);
size(prod_sgs_bar_T)
%diss_sgs_tavg_T= mean(diss_sgs_bar_T);
%size(diss_sgs_bar_T)

%fname = strcat(outputdir,'sgs_diss_tavg_',sprintf('%1.1f',Mc(i)),'_',caseid,'.dat');
%fileID = fopen(fname,'w');   
%for p=1:size(diss_sgs_tavg,2)
%fprintf(fileID,'%3.5e %3.5e\r\n',diss_sgs_tavg(p),diss_sgs_tavg_T(p));
%end
%fclose(fileID);


fname2 = strcat(inputdir1,'Rij_KE_P12_dis_pi_int_tavg_',caseid,'.dat');
f2 = readmatrix(fname2,'ConsecutiveDelimitersRule','join'); 
idx=find(abs(f2(:,1)-Mc(i)) < 0.001);del_tavg=f2(idx,2);
%1-Mc,2-del_avg,3-R11,4-R12,5-R13,6-R22,7-R23,8-R33,9-ke,10-P12,11-diss,12-pi11,13-pi12,14-pi22
%plot(y(8:end-8)/del_tavg,-diss_sgs_tavg(8:end-8),'LineWidth',1.5,'Marker',all_marks{i},'Markerindices',8:15:length(y)-16,'Color',col{i},'DisplayName',sprintf('%1.1f',Mc(i)));hold on
plot(y/del_tavg,prod_sgs_tavg_T,'LineWidth',1.5,'Marker',all_marks{i},'Markerindices',1:15:length(y),'Color',col{i},'DisplayName',sprintf('%1.1f',Mc(i)));hold on

if(i==6)
ylim([min(abs(prod_sgs_tavg_T)) max(abs(prod_sgs_tavg_T))])
end

end

ax=gca;ax.YAxis.FontSize = 15;ax.XAxis.FontSize = 15;
xlab=xlabel('\boldmath$y/{\delta_\theta(\tau)}$','Interpreter','latex');set(xlab,'Fontsize',20);
ylab=ylabel('\boldmath$P_{sgs,T}$','Interpreter','latex');set(ylab,'Fontsize',20);
xlim([-21 21])
legend show;legend('Location','Northeast','box','off','Fontsize',14,'Fontname','Timesnewroman');
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
set(findall(gcf,'-property','Markersize'),'Markersize',6)
%ax=gca;  ax.YRuler.Exponent = 0;
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
screen2jpeg(strcat(outputdir,'sgs_prod_',caseid,'.jpeg'));
