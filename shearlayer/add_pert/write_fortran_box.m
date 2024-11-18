function done = write_fortran_box(fname, writearr, format)

  done = 0;
  fid = fopen(fname, 'w');
  fwrite(fid, writearr, format);
  fclose(fid);
  done = 1;

end
