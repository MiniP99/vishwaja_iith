%%======================== Time average ============================%%
clear,clf    %%% mainid,runid,caseid,mc,time
mainid = "runs180";
runid  = "run08";
caseid = "smag_dyn";

c=sqrt(1.4);Mc=0.8;du=Mc*2*c;rho_0=1;Re_theta=1000;delta_0=1;

inputdir =strcat('/home/vsj/Codes/shearlayer_les/',mainid,'/',runid,'/',caseid,'/');
outputdir=strcat('/home/vsj/Codes/shearlayer_les/',mainid,'/',runid,'/',caseid,'/postprocess/');
outputdir1=strcat('/home/vsj/Codes/shearlayer_les/',mainid,'/',caseid,'_data/');
%inputdir =strcat('/home/vsj/Codes/check_mgm_sta_20/',mainid,'/');
%outputdir=strcat('/home/vsj/Codes/check_mgm_sta_20/',mainid,'/postprocess/');
%outputdir1=strcat('/home/vsj/Codes/check_mgm_sta_20/',mainid,'/',caseid,'_data/');
mkdir(outputdir1);

%%------------- Read Coordinates x, y, z-----------------------------
fname_coords=strcat(inputdir,'shearlayer1mat_coords.h5');
X=h5read(fname_coords,'//X');Y=h5read(fname_coords,'//Y');Z=h5read(fname_coords,'//Z');
x=X(:,1,1)';y=Y(1,:,1);z=squeeze(Z(1,1,:))';
Lx=x(end);Ly=2*y(end);Lz=z(end);
nx=size(x,2);ny=size(y,2);nz=size(z,2);
%dx = Lx/nx; dy = Ly/ny; dz = Lz/nz;
dx = X(2,1,1)-X(1,1,1); dy = Y(1,2,1)-Y(1,1,1);dz = Z(1,1,2)-Z(1,1,1);

%%----------------- Read data and find time-avg------------------------ 
fname = strcat(outputdir,'delta_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
delta1= readmatrix(fname);
t1=delta1(:,1);  %---time
d1=delta1(:,2);  %---momentum thickness
d2=delta1(:,4);  %---vorticity thickness
d99=delta1(:,6);

time=300:5:500; time_step=(time(end)-time(1))/(size(time,2)-1);

%%----Finding Momemtum thickness rate from slope of momentum thickness vs time in self-similar state----
fname = strcat(outputdir,'delta_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
f=readmatrix(fname);
idx1=find(f(:,1)==time(1));idx2=find(f(:,1)==time(end));
tau=(f(idx1:idx2,1)*du)/delta_0;
del=f(idx1:idx2,2)/delta_0;
eq=polyfit(tau,del,1);  %%eq(1,1) gives momentum thickness rate--


%%----Finding averaged intergrated Rij,TKE Budget--------------
fname = strcat(outputdir,'Rij_P12_dis_pi_nor_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
f=readmatrix(fname,'ConsecutiveDelimitersRule','join');%%% time,R11,R12,R13,R22,R23,R33,KE,P12,dissi,Pi_11,pi_12,p_22
idx1=find(f(:,1)==time(1));idx2=find(f(:,1)==time(end));
int_avg=mean(f(idx1:idx2,:));                                          

%%-----Finding averaged profiles of Reynold stress,TKE Budget-----
for i=1:size(time,2)
t=time(i);
tau=(t*du)/delta_0;
f3=zeros(ny,17);
fname1=strcat(outputdir,'vel_tke_flux_nor_',sprintf('%1.1f',Mc),'_',caseid,'_',sprintf('%04d',t),'.dat');
f3=dlmread(fname1);                   %1-ubar,2-u_til,3-vbar,4-vtil,5-R11,6-R12,7-R13,8-R22,9-R23,10-R33,11-prod,12-eps,13-trans
idx=find(t1==t);   del=d1(idx);  del_omega=d2(idx);  del_99=d99(idx);  %14-baro,15-pre_dil,16,17-mass_flux(u,v)
del_t(1,i)=del; del_t_omega(1,i)=del_omega; del_t_99(1,i)=del_99;
post_data(:,:,i)=f3;
end
post_data_tavg = sum(post_data,3)/size(time,2);
del_tavg=sum(del_t)/size(time,2); del_tavg_omega=sum(del_t_omega)/size(time,2);del_tavg_99=sum(del_t_99)/size(time,2);


%%-----Finding shearlayer centre(S is max), decorrelation length scale(delta_y),Mach(turb,grad),U_del-----

[R12_nor,M_turb_tavg,M_turbv_tavg,M_tau_tavg,M_grad_tavg,delta_y_tavg,U_del_tavg,y_c_tavg]=get_correlation_matsuno(inputdir,du,time(1),time(end),time_step,nx,ny,nz,Ly,dy,y);

%%-----write R12_nor(interms of Rey fluct) to a file------------
fname         = strcat(outputdir1,'R12_tavg_nor_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
writematrix(R12_nor',fname,'Delimiter','tab');

%%---write mach nos to a file----------
fname = strcat(outputdir1,'mach_tavg_',caseid,'.dat');
fileID = fopen(fname,'a');
fprintf(fileID,'%1.1f %2.8f %2.8f %2.8f %2.8f \r\n',Mc,M_turb_tavg,M_turbv_tavg,M_tau_tavg,M_grad_tavg);
fclose(fileID);

%%----write del_theta,del_omega,del_99,del_y,del_dot,y_c,U_del to a file----------------
fname = strcat(outputdir1,'delta_tavg_',caseid,'.dat');   
fileID = fopen(fname,'a');
fprintf(fileID,'%1.1f %2.8f %2.8f %2.8f %2.8f %2.8f %5.8f %2.8f \r\n',Mc,del_tavg,del_tavg_omega,del_99,delta_y_tavg,eq(1,1),y_c_tavg,U_del_tavg);
fclose(fileID);

%%%---write integrated value to file----------------------
fname = strcat(outputdir1,'Rij_KE_P12_dis_pi_int_tavg_',caseid,'.dat');   
fileID = fopen(fname,'a');
fprintf(fileID,'%1.1f %2.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f \r\n',Mc,del_tavg,int_avg(1,2),int_avg(1,3),int_avg(1,4),int_avg(1,5),int_avg(1,6),int_avg(1,7),int_avg(1,8),int_avg(1,9),int_avg(1,10),int_avg(1,11),int_avg(1,12),int_avg(1,13));
fclose(fileID);

%%----Write averaged profiles to a file------------
fname = strcat(outputdir1,'vel_tke_flux_tavg_nor_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
writematrix(post_data_tavg,fname,'Delimiter','tab');

%%%%%%-----------------------------------------------------------------------------------------------
%%------time-averaged TKE Budget plot-----------------------
post_data_tavg(:,12)=(-1)*post_data_tavg(:,12);
clf
figure(3)      %%%Normalised budget
plot(y(10:end-10)/del_tavg,post_data_tavg(10:end-10,11:14),'LineWidth',2);
hold on
plot(y(10:end-10)/del_tavg,post_data_tavg(10:end-10,15),'--','LineWidth',2);
legend({'P','\epsilon','T','B','\pi'},'Location','Northeast','box','off');
xlabel('y/{\delta_\theta(\tau)}','FontWeight','bold');
%h=ylabel('$\sqrt{|R_{22}|}/{\Delta{u}}$');set(h,'Interpreter','latex')
%xlim([-2.5 2.5])
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
ax=gca;  ax.YRuler.Exponent = 0;
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
screen2jpeg(strcat(outputdir,'tke_tavg_',sprintf('%1.1f',Mc),'_',caseid,'.png'));
%%%------------------------------------------------------------------------------------------------------
%%------time-averaged Reynold stress plot-----------------------
clf
figure(7),clf
R_data=[post_data_tavg(:,5),post_data_tavg(:,6),post_data_tavg(:,7),post_data_tavg(:,8),post_data_tavg(:,9),post_data_tavg(:,10)];
Rij=["R_{11}","R_{12}","R_{13}","R_{22}","R_{23}","R_{33}"];
for i=1:size(Rij,2)
txt=strcat('\surd{',Rij(1,i),'}','/ \Delta','u');
plot(y(10:end-10)/del_tavg,R_data(10:end-10,i),'DisplayName',txt,'LineWidth',1);
hold on 
end
legend show;legend('Location','Northeast','box','off');
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
xlabel('y/{\delta_\theta(\tau)}','Fontweight','bold');
%xlim([-2.5 2.5])
ax=gca; ax.YRuler.Exponent = 0;
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
screen2jpeg(strcat(outputdir,'Rij_tavg_nor_',sprintf('%1.1f',Mc),'_',caseid,'.png'));
%%%--------------------------------------------------------------------------------------------------------------
%%------time-averaged Reynold stress,dissipation subplot and also instantaneous plot-----------------------
clf
tau=(time*du)/delta_0;
col={'m','c','b','r','g','y'};

sub_data_tavg = [post_data_tavg(10:end-10,5),post_data_tavg(10:end-10,6),post_data_tavg(10:end-10,8),-post_data_tavg(10:end-10,12)];
sub_txt = ["\sqrt{|R_{11}|}/{\Delta{u}}","\sqrt{|R_{12}|}/{\Delta{u}}","\sqrt{|R_{22}|}/{\Delta{u}}","\epsilon"];
figure(22),clf
for i=1:size(time,2)
fname=strcat(outputdir,'vel_tke_flux_nor_',sprintf('%1.1f',Mc),'_',caseid,'_',sprintf('%04d',time(i)),'.dat');
f3=dlmread(fname);  %%  1-ub,    2-util,  3-vbar,  4-vtil,  5-R11,  6-R12, 7-R13, 8-R22, 9-R23, 10-R33
idx=find(t1==time(i));    %%  11-prod,12-eps,   13-trans, 14-baro, 15-pre_dil,  16-mass_ru,  17-mass_rv
del=d1(idx);
xaxis=y/del;
txt1 =strcat('\tau=',sprintf('%4.0f',tau(i)));
sub_data= [f3(10:end-10,5),f3(10:end-10,6),f3(10:end-10,8),f3(10:end-10,12)];

for k=1:4
subplot(2,2,k)
p=plot(xaxis(10:end-10),sub_data(:,k),col{k},'LineWidth',1);p.Color(4)=0.2;hold on
plot(y(10:end-10)/del_tavg,sub_data_tavg(:,k),'-k','LineWidth',1);hold on
xlabel('y/{\delta_\theta(\tau)}');xlim([-10 10]);
h=ylabel(strcat('\boldmath$', sub_txt(k),'$'));set(h,'Interpreter','latex');
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
end

end
ax=gca; ax.YRuler.Exponent = 0;
set(findall(gcf,'-property','FontSize'),'FontSize',6);set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
%legend show;leg=legend('Location','Northeast','Fontsize',6,'box','off');
%leg.ItemTokenSize = [6 6];
screen2jpeg(strcat(outputdir,'subplot_',sprintf('%1.1f',Mc),'_',caseid,'.png'));


%for i=1:size(Rij,2)
%txt=strcat('\surd{',Rij(1,i),'}','/ \Delta','u');
%plot(y/del_tavg,R_data(1:ny,i),'DisplayName',txt,'LineWidth',1);
%hold on 
%end
%legend show
%xlabel('y/{\delta_\theta(\tau)}');
%%xlim([-2.5 2.5])
%ax=gca; ax.YRuler.Exponent = 0;
%grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')me1_dns = strcat(dir1,'dns_pp_data/','prod_nor_',sprintf('%1.1f',Mc(i)),'_','dns.dat');
