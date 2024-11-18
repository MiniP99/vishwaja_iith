clear,clf    %%% mainid,runid,caseid,mc,SGS,time_step,time
mainid = "runs180";
runid  = "run02";
caseid = "mgm_dyn";
caseid_name = "mgm-dyn";

c=sqrt(1.4);Mc=0.2;du=Mc*2*c;rho_0=1;Re_theta=1000;delta_0=1;
delta_0=1;      %delta_0=Re_theta/(rho_0*Re*du);

inputdir =strcat('/home/vsj/Codes/',mainid,'/',runid,'/',caseid,'/');
inputdir1=strcat('/home/vsj/Codes/',mainid,'/');
outputdir=strcat('/home/vsj/Codes/',mainid,'/',runid,'/',caseid,'/postprocess/');
mkdir(outputdir);
outputdir1=strcat('/home/vsj/Codes/plot/myplots/',mainid,'/sgs_sij/');
mkdir(outputdir1);

myg = [0.4660 0.6740 0.1880]; myb = [0.3010 0.7450 0.9330]; myv = [0.4940 0.1840 0.5560]; myy = [0.9290 0.6940 0.1250];
col={'r','b',myg,'m','y'};
all_marks = {'o','s','d','^','v'};
%%------------- Read Coordinates x, y, z-----------------------------
fname_coords=strcat(inputdir,'shearlayer1mat_coords.h5');
X=h5read(fname_coords,'//X');Y=h5read(fname_coords,'//Y');Z=h5read(fname_coords,'//Z');
x=X(:,1,1)';y=Y(1,:,1);z=squeeze(Z(1,1,:))';
Lx=x(end);Ly=2*y(end);Lz=z(end);
nx=size(x,2);ny=size(y,2);nz=size(z,2);
%dx = Lx/nx; dy = Ly/ny; dz = Lz/nz; 
dx = X(2,1,1)-X(1,1,1); dy = Y(1,2,1)-Y(1,1,1);dz = Z(1,1,2)-Z(1,1,1);

fname = strcat(outputdir,'delta_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
delta1= readmatrix(fname);t1=delta1(:,1);d1=delta1(:,2);d2=delta1(:,4);

time_step=20;
time = [460:time_step:900];
%%-------------Start loop-------------------------------------------
for i=1:size(time,2)
counter=time(i)/time_step;
file=strcat(inputdir,'shearlayer1mat_',sprintf('%04d',counter),'.h5');

idx=find(t1==time(i));    
del=d1(idx);

%----------------Read fields data-------------------------------------
u=h5read(file,'//u');v=h5read(file,'//v');w=h5read(file,'//w');
rho=h5read(file,'//rho');p=h5read(file,'//p');
mu=h5read(file,'//mu');T=h5read(file,'//T');kap=h5read(file,'//kap');bulk=h5read(file,'//bulk');
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
% calculating  derivatives of favre_fluct--------------
du_ppdx=ddx_compact(u_pp,dx,nx,ny,nz);du_ppdy=ddy_compact(u_pp,dy,nx,ny,nz);du_ppdz=ddz_compact(u_pp,dz,nx,ny,nz);
dv_ppdx=ddx_compact(v_pp,dx,nx,ny,nz);dv_ppdy=ddy_compact(v_pp,dy,nx,ny,nz);dv_ppdz=ddz_compact(v_pp,dz,nx,ny,nz);
dw_ppdx=ddx_compact(w_pp,dx,nx,ny,nz);dw_ppdy=ddy_compact(w_pp,dy,nx,ny,nz);dw_ppdz=ddz_compact(w_pp,dz,nx,ny,nz);
% tausgs average
t11sgs_bar=avgxz(t11sgs,nx,nz);t12sgs_bar=avgxz(t12sgs,nx,nz);t13sgs_bar=avgxz(t13sgs,nx,nz);
t22sgs_bar=avgxz(t22sgs,nx,nz);t23sgs_bar=avgxz(t23sgs,nx,nz);t33sgs_bar=avgxz(t33sgs,nx,nz);

%Strain rate
S11 = dudx; S12 = 0.5*(dudy+dvdx); S13 = 0.5*(dudz+dwdx);
S22 = dvdy; S23 = 0.5*(dvdz+dwdy); S33 = dwdz;

%Subgris Dissipation------------------------------

func_tau=t11sgs.*du_ppdx+t12sgs.*du_ppdy+t13sgs.*du_ppdz+t12sgs.*dv_ppdx+t22sgs.*dv_ppdy+t23sgs.*dv_ppdz+t13sgs.*dw_ppdx+t23sgs.*dw_ppdy+t33sgs.*dw_ppdz;

diss_sgs=avgxz(func_tau,nx,nz);
diss_sgs_bar(i,:) = diss_sgs.*(del/(du^3));

%t11s11 = t11sgs.*S11;t12s12 = t12sgs.*S12;t22s22 = t22sgs.*S22;  %% Non-dimensionalise
%t11s11_bar(i,:) = avgxz(t11s11,nx,nz).*(del/(du^3));
%t12s12_bar(i,:) = avgxz(t12s12,nx,nz).*(del/(du^3));
%t22s22_bar(i,:) = avgxz(t22s22,nx,nz).*(del/(du^3));
end
%% Time average

diss_sgs_tavg= mean(diss_sgs_bar);


fname2 = strcat(inputdir1,caseid,'_data/','Rij_KE_P12_dis_pi_int_tavg_',caseid,'.dat');
f2 = readmatrix(fname2,'ConsecutiveDelimitersRule','join'); 
idx=find(abs(f2(:,1)-Mc) < 0.001);del_tavg=f2(idx,2);
%dis_text = ["\boldmath$-\tau_{11,sgs}S_{11}\delta_{\theta}/\Delta{u}^3$","\boldmath$-\tau_{12,sgs}S_{12}\delta_{\theta}/\Delta{u}^3$","\boldmath$-\tau_{22,sgs}S_{22}\delta_{\theta}/\Delta{u}^3$"];
%1-Mc,2-del_avg,3-R11,4-R12,5-R13,6-R22,7-R23,8-R33,9-ke,10-P12,11-diss,12-pi11,13-pi12,14-pi22
plot(y/del_tavg,-diss_sgs_tavg,'LineWidth',1.5,'Marker',all_marks{1},'Markerindices',1:15:length(y),'Color',col{1},'DisplayName',sprintf('1.1f',Mc));hold on
%ax=gca;  ax.YRuler.Exponent = 0; 
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
%ylim([min(abs(-t22s22_tavg)) max(abs(-t22s22_tavg))])
xlab=xlabel('\boldmath$y/{\delta_\theta(\tau)}$','Interpreter','latex');set(xlab,'Fontsize',15);
legend show;legend('Location','northeast','box','off','Fontsize',10,'interpreter','latex');
set(findall(gcf,'-property','Markersize'),'Markersize',6)
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
screen2jpeg(strcat(outputdir1,'sgs_sij_',sprintf('%1.1f',Mc),'_',caseid,'.jpeg'));
