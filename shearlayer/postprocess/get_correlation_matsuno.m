function [R12_nor,M_turb_tavg,M_turbv_tavg,M_tau_tavg,M_grad_tavg,delta_y_tavg,U_del_tavg,y_c_tavg]=get_correlation_matsuno(inputdir,du,t1,t2,time_step,nx,ny,nz,Ly,dy,y)
time=t1:time_step:t2;
for i=1:size(time,2)
t=time(i);
counter=t/time_step;                  %% Change if necessary--------
file=strcat(inputdir,'shearlayer1mat_',sprintf('%04d',counter),'.h5');
%file=strcat(inputdir,'sl_fil_',sprintf('%04d',counter),'.h5');
%read data
u=h5read(file,'//u');v=h5read(file,'//v');w=h5read(file,'//w');
rho=h5read(file,'//rho');p=h5read(file,'//p');
mu=h5read(file,'//mu');
%T=h5read(file,'//T');kap=h5read(file,'//kap');
%bulk=h5read(file,'//bulk');
t11sgs=h5read(file,'//t11');t12sgs=h5read(file,'//t12');t13sgs=h5read(file,'//t13');
t22sgs=h5read(file,'//t22');t23sgs=h5read(file,'//t23');t33sgs=h5read(file,'//t33');
%q1sgs=h5read(file,'//q1');q2sgs=h5read(file,'//q2');q3sgs=h5read(file,'//q3');
%t11sgs=zeros(nx,ny,nz);t12sgs=zeros(nx,ny,nz);t13sgs=zeros(nx,ny,nz);
%t22sgs=zeros(nx,ny,nz);t23sgs=zeros(nx,ny,nz);t33sgs=zeros(nx,ny,nz);

rho_bar=avgxz(rho,nx,nz);
t11sgs_bar=avgxz(t11sgs,nx,nz);t12sgs_bar=avgxz(t12sgs,nx,nz);t13sgs_bar=avgxz(t13sgs,nx,nz);
t22sgs_bar=avgxz(t22sgs,nx,nz);t23sgs_bar=avgxz(t23sgs,nx,nz);t33sgs_bar=avgxz(t33sgs,nx,nz);
[u_bar,u_rey_fluct]    = reynolds_avg_fluct(u,nx,ny,nz);
[v_bar,v_rey_fluct]    = reynolds_avg_fluct(v,nx,ny,nz);
[u_tilde,u_pp]=favre_avg_fluct(rho,u,nx,ny,nz);
[v_tilde,v_pp]=favre_avg_fluct(rho,v,nx,ny,nz);
[w_tilde,w_pp]=favre_avg_fluct(rho,w,nx,ny,nz);

%% R12 in terms of reynolds fluctuations-------------
[R12_bar,R12_rey_fluct]= reynolds_avg_fluct(u_rey_fluct.*v_rey_fluct,nx,ny,nz);
R12_mat(i,:)=R12_bar;
%%---------------------------------------------------

S = ddy1D_compact(u_tilde,dy,nx,ny,nz);
C = avgxz(sqrt(1.4*(p./rho)),nx,nz);

[R11,R12,R13,R22,R23,R33]=reynold_stress(rho,u_pp,v_pp,w_pp,nx,nz);
R11 = rho_bar.*R11;R12 = rho_bar.*R12;R13 = rho_bar.*R13;R22 = rho_bar.*R22;R23 = rho_bar.*R23;R33 = rho_bar.*R33;
R11=R11-t11sgs_bar;R12=R12-t12sgs_bar;R13=R13-t13sgs_bar;R22=R22-t22sgs_bar;R23=R23-t23sgs_bar;R33=R33-t33sgs_bar;
[val_R11,idx_R11]=max(abs(R11));[val_R22,idx_R22]=max(abs(R22));[val_R12,idx_R12]=max(abs(R12));

%y_c(i)=(y(idx_R11)+y(idx_R12)+y(idx_R22))./3;
%[value_yc,idx_yc] = min(abs(u_tilde)); 
[value_yc,idx_yc] = max(abs(S)); 
y_c(i) = y(idx_yc);    %[val_yc,idx_yc]=min(abs(y_c(i)-y));
S_yc = S(idx_yc);
C_yc = C(idx_yc);

%%Turbulent Mach number---
num1 = sqrt(abs((R11+R22+R33)./rho_bar));num1_yc=num1(idx_yc);
num2 = sqrt(abs(R22./rho_bar));num2_yc=num2(idx_yc);
num3 = sqrt(abs(R12./rho_bar));num3_yc=num3(idx_yc);
M_t(i) = num1_yc/C_yc;M_tv(i) = num2_yc/C_yc;M_tau(i) = num3_yc/C_yc;

%% Gradient Mach Number---
RHS(i) = 0.1* avgxz(v_rey_fluct(:,idx_yc,:).*v_rey_fluct(:,idx_yc,:),nx,nz);
LHS=zeros(1,ny/2);
for k=1:idx_yc-1%:size(delta_y,2)
index1=idx_yc-k;index2=idx_yc+k;
if(index2<=ny)
LHS(k)=avgxz(v_rey_fluct(:,index1,:).*v_rey_fluct(:,index2,:),nx,nz);
end
end
[variable,index]=min(abs(RHS(i)-LHS));
delta_y(i) = index*2*dy;
M_g(i) = S_yc*delta_y(i)/C_yc;

%u_delta-------
y1=y_c(i)+delta_y(i)/2;y2=y_c(i)-delta_y(i)/2; [val_y1,idx_y1]=min(abs(y1-y)); [val_y2,idx_y2]=min(abs(y2-y));
U_del(i)=u_tilde(idx_y1)-u_tilde(idx_y2);
end

R12_data_tavg = sum(R12_mat,1)/size(time,2);
R12_nor       = abs(R12_data_tavg)/(du*du);

M_turb_tavg = sum(M_t)/size(time,2);M_turbv_tavg = sum(M_tv)/size(time,2);M_tau_tavg = sum(M_tau)/size(time,2);
M_grad_tavg = sum(M_g)/size(time,2);
delta_y_tavg =sum(delta_y)/size(time,2);
U_del_tavg = (sum(U_del)/size(time,2))/du;       %% normalised
y_c_tavg   = sum(y_c)/size(time,2);

