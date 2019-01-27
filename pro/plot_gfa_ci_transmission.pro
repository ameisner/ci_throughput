pro plot_gfa_ci_transmission

; I don't actuall know what the "No.1", "No.2", "No.3" GFA
; transmission measurements are -- whether they're three separate 
; measurements of the same filter of measurements of three different
; physical realizations of the same GFA filter design

  readcol, '../etc/astrodon_photometrics_gen2_sloan.csv', lambda_filter_nm_ci, $
      filter_transmission_percent_ci, F='F, F'

  readcol, '../etc/gfa_filter_transmission_DESI-1297.dat', $
      lambda_filter_nm_gfa, filter_transmission_percent_gfa_1, $
                            filter_transmission_percent_gfa_2, $
                            filter_transmission_percent_gfa_3, F='F, F, F, F'


; the three GFA filter transmission curves are virtually indistinguishable

  filter_transmission_percent_gfa = (filter_transmission_percent_gfa_1 + $
                                    filter_transmission_percent_gfa_2 + $
                                    filter_transmission_percent_gfa_3)/3.0

  plot, lambda_filter_nm_gfa, filter_transmission_percent_gfa, $
      xrange=[400, 900], /xst, charsize=2, $
      xtitle='wavelength (nm)', ytitle='% transmission', yrange=[0, 100], /yst
  oplot, lambda_filter_nm_ci, filter_transmission_percent_ci, $
      color=djs_icolor('red')

  xyouts, 440, 60, 'GFA', charsize=3
  xyouts, 440, 50, 'CI', charsize=3, color=djs_icolor('red')

;  oplot, lambda_filter_nm_gfa, filter_transmission_percent_gfa_3, $
;      color=djs_icolor('green')

  bitmap = tvrd(true=1)

  write_png, 'filter_transmission-gfa_ci.png', bitmap

end
