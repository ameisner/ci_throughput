pro plot_ci_transmission

  assemble_ci_transmission, outstr=outstr, airmass=1.0

  xtitle = 'wavelength (nm)'
  ytitle = 'throughput'
  times = textoidl('\times')
  title = 'atmosphere ' + times + ' primary reflectance ' + times + $
      ' DESI corrector ' + times + ' astrodon filter ' + times + ' CI QE'
  plot, outstr.lambda_nm, outstr.transmission, charsize=2, xtitle=xtitle, $
      ytitle=ytitle, title='airmass = 1'

  xyouts, 513, 0.472, title, charsize=1.8

  bitmap = tvrd(true=1)

  write_png, 'ci_throughput.png', bitmap

end
