% Read the data from the .dat file
clear all,clc
inputdir = "/home/vsj/Codes/vortexring/noise/";
outputdir = "/home/vsj/Codes/vortexring/noise/";
mkdir(outputdir)

%Read coordinates from *_coords.h5 file
fname_coords=strcat(inputdir,'vortexring_coords.h5');
x=h5read(fname_coords,'//X');y=h5read(fname_coords,'//Y');z=h5read(fname_coords,'//Z');
Lx=x(end,1,1);  Ly=2*y(1,end,1);    Lz=2*z(1,1,end);
[nx,ny,nz] = size(x);
dxi = Lx/(nx-1);  deta = Ly/(ny-1);  dzeta= Lz/(nz-1); 

data = readmatrix(strcat(inputdir,'u_noise_00002.dat'));

%file=strcat(inputdir,'vortexring_',sprintf('%04d',1),'.h5');
%u=h5read(file,strcat('//u'));

% Extract the i, j, and A values
i = data(:, 1);
j = data(:, 2);
A_values = data(:, 3);

% Create an empty matrix A with appropriate dimensions
max_i = max(i);
max_j = max(j);
A = zeros(max_i, max_j);

% Fill the A matrix with the corresponding values
for k = 1:length(A_values)
    A(i(k), j(k)) = A_values(k);
end
%A = squeeze(u(1,:,:));
 
figure(1),clf
n=256;
t1=linspace(min(min(A)),max(max(A)),n+1);
t1=t1(1:n);
contourf(y(1,:,1), squeeze(z(1,1,:)), A,'LineColor','None','LevelList',t1);hold on
ax=gca;ax.YAxis.FontSize = 12;ax.XAxis.FontSize = 12;
xlab=xlabel('\boldmath$y(m)$');set(xlab,'Fontsize',15,'interpreter','latex');
ylab=ylabel('\boldmath$z(m)$');set(ylab,'Fontsize',15,'interpreter','latex');
c=colorbar;c.FontSize=10;%c.Ticks=[0 10000 20000];
colormap(jet(n));
xlim([-0.05 0.05]);ylim([-0.05 0.05])
file_png= strcat(outputdir,'noise_.png');
screen2jpeg(file_png)
