clear all, clf
inputdir = "/home/vsj/Codes/jet/";
outputdir = "/home/vsj/Codes/jet/postprocess/";
mkdir(postprocess)

fname = strcat(inputdir,'jet_0000.h5');
h5_info = h5info(fname)
%To find the name of the datasets{h5_info.Datasets(number)} Eg: h5_info.Datasets(1) = T
%h5_info.Attributes(1) = Time



%To read the datasets  '//T'= To read dataset (T)
u = h5read(fname,'//u'); v = h5read(fname,'//v');  w = h5read(fname,'//w');
p = h5read(fname,'//p'); rho = h5read(fname,'//rho');  T = h5read(fname,'//T');
mu= h5read(fname,'//mu'); bulk = h5read(fname,'//bulk');kap = h5read(fname,'//kap');

time =  h5readatt(fname,'/','Time');

fname_coords=strcat(inputdir,'jet_coords.h5');
X=h5read(fname_coords,'//X');Y=h5read(fname_coords,'//Y');Z=h5read(fname_coords,'//Z');
x=X(:,1,1)';y=Y(1,:,1);z=squeeze(Z(1,1,:))';
%Lx=x(end);     Ly=2*y(end);  Lz=2*z(end);
nx=size(x,2);  ny=size(y,2); nz=size(z,2);
%dx = Lx/nx;    dy= Ly/ny;    dz= Lz/nz; 

var=["u","v","w","p","rho","T","mu","kap","bulk"];

for j=1%:size(var,2)
u=h5read(file,strcat('//',var(1,j)));
u_draw=squeeze(u(1,:,:));

n=256;
t1=linspace(min(min(u_draw)),max(max(u_draw)),n+1);
t1=t1(1:n);

figure(1),clf
contourf(y, z, u_draw,'LineColor','None','LevelList',t1);
ax=gca;ax.YAxis.FontSize = 13;ax.XAxis.FontSize = 13;
xlab=xlabel('\boldmath$y$');set(xlab,'Fontsize',17,'interpreter','latex');
ylab=ylabel('\boldmath$z$');set(ylab,'Fontsize',17,'interpreter','latex');
c=colorbar;c.FontSize=13;colormap(jet(n));%caxis([-3 3])
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
file_png= strcat(outputdir,var(1,j),'_contyz_initial.png');
screen2jpeg(file_png)
end
