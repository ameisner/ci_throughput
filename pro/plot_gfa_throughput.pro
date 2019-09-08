pro plot_gfa_throughput

  assemble_gfa_transmission, outstr=outstr, airmass=1.0

  xtitle = 'wavelength (nm)'
  ytitle = 'throughput'
  times = textoidl('\times')
  title = 'atmosphere ' + times + ' primary reflectance ' + times + $
          ' DESI corrector ' + times + ' vignetting ' + times + $
          ' r filter ' + times + ' e2v QE'
  plot, outstr.lambda_nm, outstr.transmission, charsize=2, xtitle=xtitle, $
      ytitle=ytitle, title='airmass = 1'

  xyouts, 510, 0.472 + 0.07, title, charsize=1.7

  bitmap = tvrd(true=1)

  write_png, 'gfa_throughput.png', bitmap

end
