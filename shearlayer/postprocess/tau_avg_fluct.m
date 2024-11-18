function [tau11_bar,tau11_rey_fluct,tau12_bar,tau12_rey_fluct,tau13_bar,tau13_rey_fluct,tau22_bar,tau22_rey_fluct,tau23_bar,tau23_rey_fluct,tau33_bar,tau33_rey_fluct]=tau_avg_fluct(tau11,tau12,tau13,tau22,tau23,tau33,nx,ny,nz)

[tau11_bar,tau11_rey_fluct]=reynolds_avg_fluct(tau11,nx,ny,nz);
[tau12_bar,tau12_rey_fluct]=reynolds_avg_fluct(tau12,nx,ny,nz);
[tau13_bar,tau13_rey_fluct]=reynolds_avg_fluct(tau13,nx,ny,nz);
[tau22_bar,tau22_rey_fluct]=reynolds_avg_fluct(tau22,nx,ny,nz);
[tau23_bar,tau23_rey_fluct]=reynolds_avg_fluct(tau23,nx,ny,nz);
[tau33_bar,tau33_rey_fluct]=reynolds_avg_fluct(tau33,nx,ny,nz);

end
