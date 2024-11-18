function [u2, v2, w2] = scale_uvw_TI(u, v, w, TI_reqd, U_scale)

  TI_current = sqrt(mean(u.^2 + v.^2 + w.^2,'all')/3)/U_scale;
  TI_scale = TI_reqd/TI_current;
  u2 = TI_scale*u;   v2 = TI_scale*v;    w2 = TI_scale*w;
  %sqrt(mean(u2.^2 + v2.^2 + w2.^2,'all')/3)/U_scale

end
