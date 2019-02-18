pro write_throughput_csv

  fname = '../etc/ci_throughput-airmass_1.00.fits'

  outname = repstr(fname, '.fits', '.csv')

  if file_test(outname) then stop
  str = mrdfits(fname, 1)

  print, outname

  write_csv, outname, str, header=tag_names(str)

end
