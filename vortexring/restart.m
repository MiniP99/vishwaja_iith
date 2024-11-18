clear
%inputdir='/home/vsj/Codes/vortexring/pr03_165mm_512/';
%outputdir='/home/vsj/Codes/vortexring/pr03_165mm_512/';
inputdir='/home/vsj/Codes/vortexring/pr10_285mm_512/';
outputdir='/home/vsj/Codes/vortexring/pr10_285mm_512/';
mkdir(outputdir)
file=strcat(inputdir,'vortexring_0005','.h5');
k=h5info(file);
rho=h5read(file,'//rho');u=h5read(file,'//u');v=h5read(file,'//v');
w=h5read(file,'//w');p=h5read(file,'//p');
time = k.Attributes.Value;
nx = size(u,1); ny=size(u,2); nz=size(u,3);

rhou = rho.*u;
rhov = rho.*v;
rhow = rho.*w;

gam=1.4;
e = p ./ ((gam - 1)*rho);
TE = rho.*( e + 0.5*(u.*u + v.*v + w.*w));
%%------------- Read Coordinates x, y, z-----------------------------

file2=strcat(outputdir,'restart_0000.h5');
h5create(file2,'//rho',[nx ny nz]);
h5create(file2,'//rhou',[nx ny nz]);
h5create(file2,'//rhov',[nx ny nz]);
h5create(file2,'//rhow',[nx ny nz]);
h5create(file2,'//TE',[nx ny nz]);

h5write(file2,'//rho',rho);h5write(file2,'//rhou',rhou);
h5write(file2,'//rhov',rhov);h5write(file2,'//rhow',rhow);
h5write(file2,'//TE',TE);
h5writeatt(file2,'/','Time',time);
h5writeatt(file2,'/','step',1070);



