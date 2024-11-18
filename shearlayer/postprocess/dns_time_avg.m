%%======================== Time average ============================%%
clear,clf    %%% mainid,runid,caseid,mc,SGS,time_step,time
caseid = "dns";
inputdir =strcat('/home/vsj/dns_data_csl/mc0p8/postprocess/');
outputdir=strcat('/home/vsj/dns_data_csl/mc0p8/postprocess/');
outputdir1=strcat('/home/vsj/dns_data_csl/dns_pp_data/');
mkdir(outputdir1);

%%------------- Read Coordinates x, y, z-----------------------------
Lx = 100; Ly = 100; Lz = 50;
nx = 1024; ny = 1448; nz = 512;
dx = Lx/nx; x = linspace(0,Lx-dx,nx); 
dy = Ly/ny; y = linspace(-Ly/2,Ly/2,ny); 
dz = Lz/nz; z = linspace(0,Lz-dz,nz); 
c=sqrt(1.4);Mc=0.8;du=Mc*2*c;rho_0=1;Re_theta=1000;
delta_0=1;
%%----------------- Read data and find time-avg------------------------ 
fname = strcat(outputdir,'delta_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
delta1= readmatrix(fname);
t1=delta1(:,1);  %---time
d1=delta1(:,2);  %---momentum thickness


fname_pre = strcat(outputdir,'precorr_int_nor_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
pre_int=readmatrix(fname_pre);
pi_t1=pre_int(:,1);  %---time
pi11_int=pre_int(:,2);
pi12_int=pre_int(:,3);
pi22_int=pre_int(:,4);

time=140:10:160; time_step=(time(end)-time(1))/(size(time,2)-1);

%%----Finding Momemtum thickness rate from slope of momentum thickness vs time in self-similar state----
fname = strcat(outputdir,'delta_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
f=readmatrix(fname);
idx1=find(f(:,1)==time(1));idx2=find(f(:,1)==time(end));
tau=(f(idx1:idx2,1)*du)/delta_0;
del=f(idx1:idx2,2)/delta_0;
eq=polyfit(tau,del,1);  %%eq(1,1) gives momentum thickness rate--

%%-----Finding averaged profiles of Reynold stress,TKE Budget-----
for i=1:size(time,2)
t=time(i);
tau=(t*du)/delta_0;
fname1=strcat(outputdir,'vel_rey_nor_',sprintf('%1.1f',Mc),'_',caseid,'_',sprintf('%04d',t),'.dat');
fname2=strcat(outputdir,'diss_nor_',sprintf('%1.1f',Mc),'_',caseid,'_',sprintf('%04d',t),'.dat');
f4=dlmread(fname2);   %%y,prod
f3=dlmread(fname1);                   %1-ubar,2-u_til,3-vbar,4-vtil,5-R11,6-R12,7-R13,8-R22,9-R23,10-R33,11-prod,12-eps,13-trans
idx=find(t1==t);   del=d1(idx);        %14-baro,15-pre_dil,16,17-mass_flux(u,v)
del_t(1,i)=del; 
post_data(:,:,i)=f3;post_prod(:,:,i)=f4;

idxp=find(pi_t1==t);          %14-baro,15-pre_dil,16,17-mass_flux(u,v)
%pi11=pi11_int(idxp);pi12=pi12_int(idxp);pi22=pi22_int(idxp);
%pi11_t(1,i)=pi11;  pi12_t(1,i)=pi12;  pi22_t(1,i)=pi22;

end
post_data_tavg = sum(post_data,3)/size(time,2);
post_prod_tavg= sum(post_prod,3)/size(time,2);
del_tavg=sum(del_t)/size(time,2); 
%pi11_tavg=sum(pi11_t)/size(time,2); pi12_tavg=sum(pi12_t)/size(time,2); pi22_tavg=sum(pi22_t)/size(time,2); 


%%----write del_theta,del_omega,del_99,del_y,del_dot,y_c,U_del to a file----------------
%fname = strcat(outputdir1,'delta_tavg_',caseid,'.dat');   
%fileID = fopen(fname,'a');
%fprintf(fileID,'%1.1f %2.8f %2.8f \r\n',Mc,del_tavg,eq(1,1));
%fclose(fileID);

%fname = strcat(outputdir1,'precorr_nor_tavg_',caseid,'.dat');   
%fileID = fopen(fname,'a');
%fprintf(fileID,'%1.1f %12.8f %12.8f %12.8f \r\n',Mc,pi11_tavg,pi12_tavg,pi22_tavg);
%fclose(fileID);



%%----Write averaged profiles to a file------------
fname = strcat(outputdir1,'diss_nor_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
writematrix(post_prod_tavg,fname,'Delimiter','tab');
%fname = strcat(outputdir1,'vel_rey_nor_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
%writematrix(post_data_tavg,fname,'Delimiter','tab');
%fname = strcat(outputdir1,'precorr_nor_',sprintf('%1.1f',Mc),'_',caseid,'.dat');
%writematrix(post_prod_tavg,fname,'Delimiter','tab');

%%%%%%-----------------------------------------------------------------------------------------------
%%------time-averaged TKE Budget plot-----------------------
clf
figure(3)      %%%Normalised budget
plot(y/del_tavg,post_data_tavg(:,3),'LineWidth',2);hold on
plot(y/del_tavg,post_data_tavg(:,4),'LineWidth',2);hold on
plot(y/del_tavg,post_data_tavg(:,5),'LineWidth',2);hold on
%plot(y/del_tavg,post_prod_tavg(:,2),'LineWidth',2);hold on
%legend({'$\sqrt{|R_{11}|}/{\Delta{u}}$','$\sqrt{|R_{12}|}/{\Delta{u}}$','$\sqrt{|R_{22}|}/{\Delta{u}}$'},'Location','Northeast','box','off','Interpreter','latex');
xlabel('y/{\delta_\theta(\tau)}','FontWeight','bold');
%h=ylabel('$\sqrt{|R_{22}|}/{\Delta{u}}$');set(h,'Interpreter','latex')
%xlim([-2.5 2.5])
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
ax=gca;  ax.YRuler.Exponent = 0;
grid on;hAx=gca;set(hAx,'xminorgrid','off','yminorgrid','off')
screen2jpeg(strcat(outputdir,'p11_',sprintf('%1.1f',Mc),'_',caseid,'.png'));
