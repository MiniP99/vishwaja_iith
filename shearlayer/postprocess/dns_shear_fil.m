clear,clf    %%% mainid,runid,caseid,mc,SGS,time_step,time
caseid = "dns";

inputdir =strcat('/home/vsj/dns_data_csl/mc0p2/my_fil/');
inputdir1 =strcat('/home/vsj/dns_data_csl/mc0p2/diss_fil/');
outputdir=strcat('/home/vsj/dns_data_csl/mc0p2/postprocess/');
mkdir(outputdir);mkdir(inputdir1);

%%------------- Read Coordinates x, y, z-----------------------------
Lx = 150; Ly = 200; Lz = 75;
nx = 1024; ny = 1448; nz = 512;
dx = Lx/nx; x = linspace(0,Lx-dx,nx); 
dy = Ly/ny; y = linspace(-Ly/2,Ly/2,ny); 
dz = Lz/nz; z = linspace(0,Lz-dz,nz); 
c=sqrt(1.4);Mc=0.2;du=Mc*2*c;rho_0=1;Re_theta=1000;
delta_0=1;

%%-------------Start loop-------------------------------------------
for i=101:113
counter=i;
time=counter*10;
%file=strcat(inputdir,'shearlayer1mat_',sprintf('%04d',counter),'.h5');
file=strcat(inputdir,'sl_fil',sprintf('%04d',counter),'.h5');

%----------------Read fields data-------------------------------------

u=h5read(file,'//u');v=h5read(file,'//v');w=h5read(file,'//w');
rho=h5read(file,'//rho');
%p=h5read(file,'//p');
%mu=h5read(file,'//mu');bulk=h5read(file,'//bulk');
%calculating derivatives-----------------------------
%dudx=ddx_compact(u,Lx,nx,ny,nz);dudy=ddy_compact(u,Ly,nx,ny,nz);dudz=ddz_compact(u,Lz,nx,ny,nz);
%dvdx=ddx_compact(v,Lx,nx,ny,nz);dvdy=ddy_compact(v,Ly,nx,ny,nz);dvdz=ddz_compact(v,Lz,nx,ny,nz);
%dwdx=ddx_compact(w,Lx,nx,ny,nz);dwdy=ddy_compact(w,Ly,nx,ny,nz);dwdz=ddz_compact(w,Lz,nx,ny,nz);
%
%tau11 =  ((4/3)*mu + bulk).*dudx + (bulk-(2/3)*mu).*(dvdy+dwdz) ;
%tau12 = mu.*(dudy+dvdx);
%tau13 = mu.*(dudz+dwdx);
%tau22 = ((4/3)*mu + bulk).*dvdy + (bulk-(2/3)*mu).*(dudx+dwdz) ;
%tau23 = mu.*(dvdz+dwdy);
%tau33 = ((4/3)*mu + bulk).*dwdz + (bulk-(2/3)*mu).*(dvdy+dudx) ;
%
%file2=strcat(inputdir1,'sl_shear_fil',sprintf('%04d',counter),'.h5');
%h5create(file2,'//tau11',[nx ny nz]);h5create(file2,'//tau12',[nx ny nz]);h5create(file2,'//tau13',[nx ny nz]);
%h5create(file2,'//tau22',[nx ny nz]);h5create(file2,'//tau23',[nx ny nz]);h5create(file2,'//tau33',[nx ny nz]);
%
%h5write(file2,'//tau11',tau11);h5write(file2,'//tau12',tau12);h5write(file2,'//tau13',tau13);
%h5write(file2,'//tau22',tau22);h5write(file2,'//tau23',tau23);h5write(file2,'//tau33',tau33);
%
% calculating averages,fluctuations------------------
[u_tilde,u_pp]=favre_avg_fluct(rho,u,nx,ny,nz);
[v_tilde,v_pp]=favre_avg_fluct(rho,v,nx,ny,nz);
[w_tilde,w_pp]=favre_avg_fluct(rho,w,nx,ny,nz);
% calculating  derivatives of favre_fluct--------------
du_ppdx=ddx_compact(u_pp,Lx,nx,ny,nz);du_ppdy=ddy_compact(u_pp,Ly,nx,ny,nz);du_ppdz=ddz_compact(u_pp,Lz,nx,ny,nz);
dv_ppdx=ddx_compact(v_pp,Lx,nx,ny,nz);dv_ppdy=ddy_compact(v_pp,Ly,nx,ny,nz);dv_ppdz=ddz_compact(v_pp,Lz,nx,ny,nz);
dw_ppdx=ddx_compact(w_pp,Lx,nx,ny,nz);dw_ppdy=ddy_compact(w_pp,Ly,nx,ny,nz);dw_ppdz=ddz_compact(w_pp,Lz,nx,ny,nz);

file2=strcat(inputdir1,'sl_shear_fil',sprintf('%04d',counter),'.h5');
h5create(file2,'//duppdx',[nx ny nz]);h5create(file2,'//duppdy',[nx ny nz]);h5create(file2,'//duppdz',[nx ny nz]);
h5create(file2,'//dvppdx',[nx ny nz]);h5create(file2,'//dvppdy',[nx ny nz]);h5create(file2,'//dvppdz',[nx ny nz]);
h5create(file2,'//dwppdx',[nx ny nz]);h5create(file2,'//dwppdy',[nx ny nz]);h5create(file2,'//dwppdz',[nx ny nz]);

h5write(file2,'//duppdx',du_ppdx);h5write(file2,'//duppdy',du_ppdy);h5write(file2,'//duppdz',du_ppdz);
h5write(file2,'//dvppdx',dv_ppdx);h5write(file2,'//dvppdy',dv_ppdy);h5write(file2,'//dvppdz',dv_ppdz);
h5write(file2,'//dwppdx',dw_ppdx);h5write(file2,'//dwppdy',dw_ppdy);h5write(file2,'//dwppdz',dw_ppdz);

end