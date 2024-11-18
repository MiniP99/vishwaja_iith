clear,clf    %%% mainid,runid,caseid,mc,SGS,time_step,time
caseid = "dns";

inputdir =strcat('/home/vsj/dns_data_csl/mc2p0/');
inputdir1 =strcat('/home/vsj/dns_data_csl/mc2p0/my_fil/');
outputdir=strcat('/home/vsj/dns_data_csl/mc2p0/postprocess/');
mkdir(outputdir);
mkdir(inputdir1);

%%------------- Read Coordinates x, y, z-----------------------------
Lx = 80; Ly = 80; Lz = 40;
nx = 1024; ny = 1448; nz = 512;
dx = Lx/nx; x = linspace(0,Lx-dx,nx); 
dy = Ly/ny; y = linspace(-Ly/2,Ly/2,ny); 
dz = Lz/nz; z = linspace(0,Lz-dz,nz); 
c=sqrt(1.4);Mc=2.0;du=Mc*2*c;rho_0=1;Re_theta=1000;

%%------------- Read Coordinates x, y, z for LES-----------------------------
nx_les = 180; ny_les = 240; nz_les = 84;
dx_les = Lx/nx_les; x_les = linspace(0,Lx-dx_les,nx_les); 
dy_les = Ly/ny_les; y_les = linspace(-Ly/2,Ly/2,ny_les); 
dz_les = Lz/nz_les; z_les = linspace(0,Lz-dz_les,nz_les); 
%%-------------Start loop-------------------------------------------
for i=5:11
counter=i;
time=counter*10;
%file=strcat(inputdir,'shearlayer1mat_',sprintf('%04d',counter),'.h5');
file=strcat(inputdir,'sl',sprintf('%04d',counter),'.h5');

%----------------Read fields data-------------------------------------
u=h5read(file,'//u');v=h5read(file,'//v');w=h5read(file,'//w');
rho=h5read(file,'//rho');p=h5read(file,'//p');
mu=h5read(file,'//mu');bulk=h5read(file,'//bulk');
T=h5read(file,'//T');

Nx = round(((Lx/nx_les)/(Lx/nx))/2);
Ny = round(((Ly/ny_les)/(Ly/ny))/2);
Nz = round(((Lz/nz_les)/(Lz/nz))/2);


u_fil = box_filter(u,Nx,Ny,Nz);
v_fil = box_filter(v,Nx,Ny,Nz);
w_fil = box_filter(w,Nx,Ny,Nz);
T_fil = box_filter(T,Nx,Ny,Nz);
rho_fil = box_filter(rho,Nx,Ny,Nz);
p_fil = box_filter(p,Nx,Ny,Nz);
bulk_fil = box_filter(bulk,Nx,Ny,Nz);
mu_fil = box_filter(mu,Nx,Ny,Nz);

[nx_new ny_new nz_new] = size(u_fil);

file2=strcat(inputdir1,'sl_fil',sprintf('%04d',counter),'.h5');
h5create(file2,'//u',[nx_new ny_new nz_new]);h5create(file2,'//v',[nx_new ny_new nz_new]);
h5create(file2,'//w',[nx_new ny_new nz_new]);h5create(file2,'//T',[nx_new ny_new nz_new]);
h5create(file2,'//rho',[nx_new ny_new nz_new]);h5create(file2,'//p',[nx_new ny_new nz_new]);
h5create(file2,'//bulk',[nx_new ny_new nz_new]);h5create(file2,'//mu',[nx_new ny_new nz_new]);

h5write(file2,'//u',u_fil);h5write(file2,'//v',v_fil);h5write(file2,'//w',w_fil);
h5write(file2,'//rho',rho_fil);h5write(file2,'//p',p_fil);h5write(file2,'//T',T_fil);
h5write(file2,'//bulk',bulk_fil);h5write(file2,'//mu',mu_fil);

%dx_new = Lx/nx_new; dy_new = Ly/ny_new; dz_new = Lz/nz_new;
%
%x_new = linspace(0,Lx-dx_new,nx_new); 
%y_new = linspace(-Ly/2,Ly/2,ny_new); 
%z_new = linspace(0,Lz-dz_new,nz_new);
% 
%dil_draw = squeeze(u(:,:,nz/2))';
%dil_draw_fil = squeeze(u_fil(:,:,nz_new/2))';
%n=256;
%lev=linspace(min(min(dil_draw)),max(max(dil_draw)),n+1);
%lev=lev(1:n);
%contourf(x, y, dil_draw,'LineColor','None','LevelList',lev)
%screen2jpeg(strcat(outputdir,'uplot_',sprintf('%1.1f',Mc),'_',caseid,'.png'));
%
%lev=linspace(min(min(dil_draw_fil)),max(max(dil_draw_fil)),n+1);
%lev=lev(1:n);
%contourf(x_new, y_new, dil_draw_fil,'LineColor','None','LevelList',lev)
%screen2jpeg(strcat(outputdir,'uplot_fil_',sprintf('%1.1f',Mc),'_',caseid,'.png'));
end
