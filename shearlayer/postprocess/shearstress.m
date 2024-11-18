function [tau11,tau12,tau13,tau22,tau23,tau33]=shearstress(mu,bulk,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,t11sgs,t12sgs,t13sgs,t22sgs,t23sgs,t33sgs)
lambda=bulk-(2/3)*mu;
bambda=(4/3)*mu + bulk;
tau11=bambda.*dudx+lambda.*(dvdy+dwdz);
tau12=mu.*(dudy+dvdx);
tau13=mu.*(dudz+dwdx);
tau22=bambda.*dvdy+lambda.*(dudx+dwdz);
tau23=mu.*(dvdz+dwdy);
tau33=bambda.*dwdz+lambda.*(dvdy+dudx);
tau11=tau11-t11sgs;
tau12=tau12-t12sgs;
tau13=tau13-t13sgs;
tau22=tau22-t22sgs;
tau23=tau23-t23sgs;
tau33=tau33-t33sgs;
end

