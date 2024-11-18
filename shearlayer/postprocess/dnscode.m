%%============ Code for finding momentum thickness, del_99, TKE Budgets ====================%%

clear,clf    %%% mainid,runid,caseid,mc,SGS,time_step,time
caseid = "dns";

inputdir =strcat('/home/vsj/dns_data_csl/mc0p8/');
outputdir=strcat('/home/vsj/dns_data_csl/mc0p8/postprocess/');
mkdir(outputdir);

%%------------- Read Coordinates x, y, z-----------------------------
Lx = 100; Ly = 100; Lz = 50;
nx = 1024; ny = 1448; nz = 512;
dx = Lx/nx; x = linspace(0,Lx-dx,nx); 
dy = Ly/ny; y = linspace(-Ly/2,Ly/2,ny); 
dz = Lz/nz; z = linspace(0,Lz-dz,nz); 
c=sqrt(1.4);Mc=0.8;du=Mc*2*c;rho_0=1;Re_theta=1000;
delta_0=1;

%%-------------Start loop-------------------------------------------
for i=14:16
counter=i
time=counter*10;
file=strcat(inputdir,'diss/sl_shear',sprintf('%04d',counter),'.h5');
file1=strcat(inputdir,'sl',sprintf('%04d',counter),'.h5');

%----------------Read fields data-------------------------------------

u=h5read(file1,'//u');%v=h5read(file1,'//v');w=h5read(file1,'//w');
rho=h5read(file1,'//rho');
%p=h5read(file1,'//p');
%mu=h5read(file1,'//mu');bulk=h5read(file1,'//bulk');

tau11=h5read(file,'//tau11');tau12=h5read(file,'//tau12');tau13=h5read(file,'//tau13');
tau22=h5read(file,'//tau22');tau23=h5read(file,'//tau23');tau33=h5read(file,'//tau33');
du_ppdx=h5read(file,'//duppdx');du_ppdy=h5read(file,'//duppdy');du_ppdz=h5read(file,'//duppdz');
dv_ppdx=h5read(file,'//dvppdx');dv_ppdy=h5read(file,'//dvppdy');dv_ppdz=h5read(file,'//dvppdz');
dw_ppdx=h5read(file,'//dwppdx');dw_ppdy=h5read(file,'//dwppdy');dw_ppdz=h5read(file,'//dwppdz');

%[p_bar,p_reynold_fluct]    = reynolds_avg_fluct(p,nx,ny,nz);
[u_tilde,u_pp]=favre_avg_fluct(rho,u,nx,ny,nz);
%[v_tilde,v_pp]=favre_avg_fluct(rho,v,nx,ny,nz);
%[w_tilde,w_pp]=favre_avg_fluct(rho,w,nx,ny,nz);

% calculating  derivatives of favre_fluct--------------
%du_ppdx=ddx_compact(u_pp,dx,nx,ny,nz);du_ppdy=ddy_compact(u_pp,dy,nx,ny,nz);%du_ppdz=ddz_compact(u_pp,dz,nx,ny,nz);
%dv_ppdx=ddx_compact(v_pp,dx,nx,ny,nz);dv_ppdy=ddy_compact(v_pp,dy,nx,ny,nz);%dv_ppdz=ddz_compact(v_pp,dz,nx,ny,nz);
%R22=favre(rho,v_pp.*v_pp,nx,nz);
%R12=favre(rho,u_pp.*v_pp,nx,nz);R23=favre(rho,v_pp.*w_pp,nx,nz);

% calculating derivatives of favre_avg
%du_tildy=ddy1D_compact(u_tilde,dy,nx,ny,nz);dv_tildy=ddy1D_compact(v_tilde,dy,nx,ny,nz);dw_tildy=ddy1D_compact(w_tilde,dy,nx,ny,nz);

%u_tilde=avgxz(rho.*u,nx,nz)./avgxz(rho,nx,nz);
func1=avgxz(rho,nx,nz).*(0.5*du-u_tilde).*(0.5*du+u_tilde);
delta_theta=(1/(rho_0*(du^2)))*sum(func1)*dy;

%Pij = ((-(R12.*du_tildy+R22.*dv_tildy+R23.*dw_tildy).*avgxz(rho,nx,nz))*delta_theta)/(du^3);
%P12 = ((-(R12.*du_tildy).*avgxz(rho,nx,nz))*delta_theta)/(du^3);


%%calculating derivatives-----------------------------
%dudx=ddx_compact(u,dx,nx,ny,nz);dudy=ddy_compact(u,dy,nx,ny,nz);dudz=ddz_compact(u,dz,nx,ny,nz);
%dvdx=ddx_compact(v,dx,nx,ny,nz);dvdy=ddy_compact(v,dy,nx,ny,nz);dvdz=ddz_compact(v,dz,nx,ny,nz);
%dwdx=ddx_compact(w,dx,nx,ny,nz);dwdy=ddy_compact(w,dy,nx,ny,nz);dwdz=ddz_compact(w,dz,nx,ny,nz);
%% calculating  derivatives of favre_fluct--------------
%du_ppdx=ddx_compact(u_pp,dx,nx,ny,nz);du_ppdy=ddy_compact(u_pp,dy,nx,ny,nz);%du_ppdz=ddz_compact(u_pp,dz,nx,ny,nz);
%dv_ppdx=ddx_compact(v_pp,dx,nx,ny,nz);dv_ppdy=ddy_compact(v_pp,dy,nx,ny,nz);%dv_ppdz=ddz_compact(v_pp,dz,nx,ny,nz);
%dw_ppdx=ddx_compact(w_pp,dx,nx,ny,nz);dw_ppdy=ddy_compact(w_pp,dy,nx,ny,nz);dw_ppdz=ddz_compact(w_pp,dz,nx,ny,nz);
%
%tau11 =  ((4/3)*mu + bulk).*dudx + (bulk-(2/3)*mu).*(dvdy+dwdz) ;
%tau12 = mu.*(dudy+dvdx);
%tau13 = mu.*(dudz+dwdx);
%tau22 = ((4/3)*mu + bulk).*dvdy + (bulk-(2/3)*mu).*(dudx+dwdz) ;
%tau23 = mu.*(dvdz+dwdy);
%tau33 = ((4/3)*mu + bulk).*dwdz + (bulk-(2/3)*mu).*(dvdy+dudx) ;
%
diss=avgxz(tau11.*du_ppdx+tau12.*du_ppdy+tau13.*du_ppdz+tau12.*dv_ppdx+tau22.*dv_ppdy+tau23.*dv_ppdz+tau13.*dw_ppdx+tau23.*dw_ppdy+tau33.*dw_ppdz,nx,nz).*(delta_theta/(du^3));

%Pressure-strain terms----------------------
%pi_11 = avgxz(p_reynold_fluct.*2.*du_ppdx,nx,ny).*(delta_theta/(du^3));
%pi_12 = avgxz(p_reynold_fluct.*(du_ppdy+dv_ppdx),nx,ny).*(delta_theta/(du^3));
%pi_22 = avgxz(p_reynold_fluct.*2.*dv_ppdy,nx,ny).*(delta_theta/(du^3));

%%To find centreline
%S = ddy1D_compact(u_tilde,dy,nx,ny,nz);
%[value_yc,idx_yc] = max(abs(S)); 
%y_c  = y(idx_yc);    %[val_yc,idx_yc]=min(abs(y_c(i)-y));


%%%---------------------------Non-dimensional terms-------------------------------
%fname=strcat(outputdir,'precorr_nor_',sprintf('%1.1f',Mc),'_',caseid,'_',sprintf('%04d',time),'.dat');
%fileID=fopen(fname,'w');
%for j=1:ny
%fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f\n',y(j),pi_11(1,j),pi_12(1,j),pi_22(1,j));
%end
%fclose(fileID);

fname=strcat(outputdir,'diss_nor_',sprintf('%1.1f',Mc),'_',caseid,'_',sprintf('%04d',time),'.dat');
fileID=fopen(fname,'w');
for j=1:ny
fprintf(fileID,'%12.8f %12.8f \n',y(j),diss(1,j));
end
fclose(fileID);

%fname=strcat(outputdir,'y_c_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
%fileID=fopen(fname,'a');
%fprintf(fileID,'%12.8f %12.8f\n',idx_yc,y_c);
%fclose(fileID);

%fname = strcat(outputdir,'delta_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
%fileID = fopen(fname,'a');   %%% U_del/du
%fprintf(fileID,'%04d %2.8f %2.8f %2.8f %2.8f %2.8f %2.8f %2.8f %2.8f\r\n',time,delta_theta,delta_dot,delta_omega,D_w,delta_99,delta_y,y_c,U_del/du);
%fclose(fileID);
end
