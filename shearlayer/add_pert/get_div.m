function maxdiv = get_div(u, v, w, kx3D, ky3D, kz3D)

  dudx = ddx(u, kx3D);   dvdy = ddy(v, ky3D);   dwdz = ddz(w, kz3D);
  dilat = dudx + dvdy + dwdz;
  maxdiv = max(abs(dilat),[],'all');

end
