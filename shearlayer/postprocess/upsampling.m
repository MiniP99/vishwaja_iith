clear
inputdir='/home/vsj/Codes/runs256/run02/mgm/';
outputdir='/home/vsj/Codes/runs256/run02/mgm/postprocess/';
mkdir(outputdir)
file=strcat(inputdir,'restart_0001','.h5');
rho=h5read(file,'//rho');rhou=h5read(file,'//rhou');rhov=h5read(file,'//rhov');
rhow=h5read(file,'//rhow');TE=h5read(file,'//TE');

%%------------- Read Coordinates x, y, z-----------------------------
fname_coords=strcat(inputdir,'shearlayer1mat_coords.h5');
X=h5read(file,'//X');Y=h5read(file,'//Y');Z=h5read(file,'//Z');
x=X(:,1,1)';y=Y(1,:,1);z=squeeze(Z(1,1,:))';
Lx=x(end);Ly=2*y(end);Lz=z(end);
nx=size(x,2);ny=size(y,2);nz=size(z,2);
dx = Lx/nx; dy = Ly/ny; dz = Lz/nz; 
[Y,X,Z]=meshgrid(y,x,z);


nx_up=512;ny_up=720;nz_up=256;
Lx=80;Ly=80;Lz=40;
dx_up = Lx/nx_up; x_up = linspace(0,Lx-dx_up,nx_up); 
dy_up = Ly/ny_up; y_up = linspace(-Ly/2,Ly/2,ny_up); 
dz_up = Lz/nz_up; z_up = linspace(0,Lz-dz_up,nz_up);
[X_up,Y_up,Z_up]=meshgrid(y_up,x_up,z_up);

rho_up=interp3(X,Y,Z,rho,X_up,Y_up,Z_up);rhou_up=interp3(X,Y,Z,rhou,X_up,Y_up,Z_up);
rhov_up=interp3(X,Y,Z,rhov,X_up,Y_up,Z_up);rhow_up=interp3(X,Y,Z,rhow,X_up,Y_up,Z_up);
TE_up=interp3(X,Y,Z,TE,X_up,Y_up,Z_up);

file2='restart_0001.h5';
h5create(file2,'//rho',[nx_up ny_up nz_up]);
h5create(file2,'//rhou',[nx_up ny_up nz_up]);
h5create(file2,'//rhov',[nx_up ny_up nz_up]);
h5create(file2,'//rhow',[nx_up ny_up nz_up]);
h5create(file2,'//TE',[nx_up ny_up nz_up]);

h5write(file2,'//rho',rho_up);h5write(file2,'//rhou',rhou_up);
h5write(file2,'//rhov',rhov_up);h5write(file2,'//rhow',rhow_up);
h5write(file2,'//TE',TE_up);
h5writeatt(file2,'/','Time',568.7874);
h5writeatt(file2,'/','step',5000);



