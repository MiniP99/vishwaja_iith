clear,clf
inputdir='/home/vsj/Codes/shearlayer_les/plot/';
outputdir='/home/vsj/Codes/shearlayer_les/plot/';
outputdir1='/home/vsj/Codes/shearlayer_les/plot/myplots/comp_smag/';
mkdir(outputdir1);

c=sqrt(1.4);Mc=2.0;du=Mc*2*c;rho_0=1;Re_theta=1000;
delta_0=1;
caseid = ["base","lad","amd_dyn","sigma_dyn","mgm_dyn","smag_dyn"];
caseid_name = ["Base","LAD","AMD","Sigma","MGM","Smag"];
myg = [0.4660 0.6740 0.1880]; myb = [0.3010 0.7450 0.9330]; myv = [0.4940 0.1840 0.5560]; myy = [0.9290 0.6940 0.1250];
col={myv,myy,myg,'r','b','g'};
all_marks = {'o','v','d','^','s','>','x','v','>','<','p','h'};
%%------------- Read Coordinates x, y, z-----------------------------
fname_coords=strcat('/home/vsj/Codes/shearlayer_les/runs180/run20/smag_dyn/shearlayer1mat_coords.h5');
X=h5read(fname_coords,'//X'); Y=h5read(fname_coords,'//Y'); Z=h5read(fname_coords,'//Z');
x=X(:,1,1)'; y=Y(1,:,1); z=squeeze(Z(1,1,:))';
Lx=x(end);Ly=2*y(end);Lz=z(end);
nx=size(x,2);ny=size(y,2);nz=size(z,2);
dx = Lx/nx; dy = Ly/ny; dz = Lz/nz; 
%----------------------------------------------------------------------
figure(1),clf

for i=2:size(caseid,2)
fname1=strcat('/home/vsj/Codes/shearlayer_les/runs180/',caseid(i),'_data/','vel_tke_flux_tavg_nor_',sprintf('%1.1f',Mc),'_',caseid(i),'.dat');
fname_theta=strcat('/home/vsj/Codes/shearlayer_les/runs180/',caseid(i),'_data/','delta_tavg_',caseid(i),'.dat');

f1 = readmatrix(fname1,'ConsecutiveDelimitersRule','join');
f2 = readmatrix(fname_theta,'ConsecutiveDelimitersRule','join');
idx=find(abs(f2(:,1)-Mc) < 0.001);del_tavg=f2(idx,2);

sub_txt = ["\sqrt{|R_{11}|}/{\Delta{u}}","\sqrt{|R_{12}|}/{\Delta{u}}","\sqrt{|R_{22}|}/{\Delta{u}}","\epsilon"];
sub_data= [f1(:,5),f1(:,6),f1(:,8),f1(:,12)];

for k=1:3
subplot(2,2,k)
plot(y/del_tavg,sub_data(:,k),'Color',col{i},'DisplayName',caseid_name(i));hold on

if i==size(caseid,2)
ny_dns=1448;
dy_dns = Ly/ny_dns; y_dns = linspace(-Ly/2,Ly/2,ny_dns)'; 
fname1_dns = strcat('/home/vsj/dns_data_csl/dns_pp_data/','vel_rey_nor_',sprintf('%1.1f',Mc),'_','dns.dat');
f1_dns = readmatrix(fname1_dns,'ConsecutiveDelimitersRule','join');  
fname2_dns = strcat('/home/vsj/dns_data_csl/dns_pp_data/','delta_tavg_','dns.dat');
f2_dns = readmatrix(fname2_dns,'ConsecutiveDelimitersRule','join'); 
idx=find(abs(f2_dns(:,1)-Mc) < 0.001); %Mc,delta_theta
del_tavg_dns=f2_dns(idx,2);
plot(y_dns/del_tavg_dns,f1_dns(:,k+2),'marker','o','color','k','LineWidth',1,'Linestyle','none','Displayname','DNS','markerindices',1:40:length(y_dns));hold on
end

xlab=xlabel('\boldmath$y/{\delta_\theta(\tau)}$','Interpreter','latex');set(xlab,'Fontsize',15);
xlim([-10 10]);
h=ylabel(strcat('\boldmath$', sub_txt(k),'$'));set(h,'Interpreter','latex');
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
legend show;leg=legend('Location','Northeast','Fontsize',6,'box','off');leg.ItemTokenSize=[5 5]
end


subplot(2,2,4)
%%------time-averaged TKE Budget plot-----------------------
f1(:,12) = -f1(:,12);
plot(y/del_tavg,f1(:,11),'Color',col{i},'LineStyle','-','HandleVisibility','off');hold on
plot(y/del_tavg,f1(:,12),'Color',col{i},'LineStyle','-.','HandleVisibility','off');hold on
if i==size(caseid,2)
ny_dns=1448;
dy_dns = Ly/ny_dns; y_dns = linspace(-Ly/2,Ly/2,ny_dns)';
diss_dns = strcat('/home/vsj/dns_data_csl/dns_pp_data/','diss_nor_',sprintf('%1.1f',Mc),'_','dns.dat');
prod_dns = strcat('/home/vsj/dns_data_csl/dns_pp_data/','prod_nor_',sprintf('%1.1f',Mc),'_','dns.dat');
diss_dns = readmatrix(diss_dns,'ConsecutiveDelimitersRule','join');  
prod_dns = readmatrix(prod_dns,'ConsecutiveDelimitersRule','join');  
fname2_dns = strcat('/home/vsj/dns_data_csl/dns_pp_data/','delta_tavg_','dns.dat');
f2_dns = readmatrix(fname2_dns,'ConsecutiveDelimitersRule','join'); 
idx=find(abs(f2_dns(:,1)-Mc) < 0.001); %Mc,delta_theta
del_tavg_dns=f2_dns(idx,2);
plot(y_dns/del_tavg_dns,-diss_dns(:,2),'marker','o','color','k','LineWidth',1,'Linestyle','none','Displayname','DNS-Diss','markerindices',1:40:length(y_dns));hold on
plot(y_dns/del_tavg_dns,prod_dns(:,2),'marker','+','color','k','LineWidth',1,'Linestyle','none','Displayname','DNS-Prod','markerindices',1:40:length(y_dns));hold on
end 
txt1 = ' Prod(P)\rightarrow';
t = text(-5,0.0004,txt1);
txt2 = 'Diss(\epsilon)\rightarrow';
t = text(-5,-0.0002,txt2);
%plot(y/del_tavg,f1(:,13),'Color',col{i},'LineStyle','-');hold on
legend show;leg=legend('Location','Northeast','Fontsize',6,'box','off');leg.ItemTokenSize=[10 10];
xlabel('y/{\delta_\theta(\tau)}','FontWeight','bold');
xlim([-5 5]);
ax=gca;  ax.YRuler.Exponent = 0;
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off');
end

ax=gca;  ax.YRuler.Exponent = 0;
set(findall(gcf,'-property','FontSize'),'FontSize',7);
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold');
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
set(findall(gcf,'-property','Markersize'),'Markersize',3);
tit=sgtitle(strcat('\boldmath$ M_{c}=',sprintf('%1.1f',Mc),'$'));set(tit,'Interpreter','latex');

screen2jpeg(strcat(outputdir1,'comp_smg_',sprintf('%1.1f',Mc),'.png'));

