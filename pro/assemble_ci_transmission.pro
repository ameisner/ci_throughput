pro assemble_ci_transmission, airmass=airmass, outstr=outstr, write=write

; common grid onto which to interpolate various inputs
; only need to consider 550 nm to 700 nm based on astrodon r' filter curve
  lambda_common_nm = 500.0 + findgen(251.0)

;  atmosphere*aluminum*corrector*filter*QE

  if ~keyword_set(airmass) then airmass = 1.0

  atm_str = mrdfits('../etc/ZenithExtinction-KPNO.fits', 1)

; will be creating transmission curve on 1 nm spaced grid
; so want to smooth atmospheric extinction to be critically sampled on
; that grid

; ZenithExtinction-KPNO.fits is uniformly gridded in wavelength, with
; spacing of 0.1 angstrom

; sigma = 10 angstroms
  atm_extinction_smooth = gauss_smooth(atm_str.extinction, 100.0)

  lambda_atm_nm = atm_str.wavelength/10.0
  atm_transmission_frac = 10^(-1.0*airmass*atm_extinction_smooth/2.5)

  readcol, '../etc/primary_reflectance.dat', lambda_mirror_nm, $
      primary_reflectance_frac, F='F, F'

  readcol, '../etc/desi_corrector_throughput.dat', lambda_corrector_nm, $
      corrector_transmission_frac, F='F, F'

  readcol, '../etc/astrodon_photometrics_gen2_sloan.csv', lambda_filter_nm, $
      filter_transmission_percent, F='F, F'

  filter_transmission_frac = filter_transmission_percent/100.0

  readcol, '../etc/ci_ccd_qe_fig5.csv', lambda_qe_nm, qe_frac, F='F, F'

  window, xsize=950, ysize=650
  xtitle = 'wavelength (nm)'
  ytitle = 'throughput'
  plot, lambda_mirror_nm, primary_reflectance_frac, xrange=[350, 1055], /xst, $
      charsize=1.4, xtitle=xtitle, ytitle=ytitle
  oplot, lambda_atm_nm, atm_transmission_frac, color=djs_icolor('blue')
  oplot, lambda_corrector_nm, corrector_transmission_frac, $
      color=djs_icolor('green')
  oplot, lambda_filter_nm, filter_transmission_frac, color=djs_icolor('red')
  oplot, lambda_qe_nm, qe_frac, color=djs_icolor('yellow')

  xyouts, 775, 0.6, 'atmosphere', charsize=2, color=djs_icolor('blue')
  xyouts, 775, 0.55, 'airmass = ' + string(airmass, format='(F4.2)'), $
      charsize=2, color=djs_icolor('blue')

  xyouts, 875, 0.3, 'KAF-6303 QE', charsize=2, color=djs_icolor('yellow')

  xyouts, 412.5, 0.15, "astrodon r' filter", charsize=2, color=djs_icolor('red')

  xyouts, 370, 0.91, charsize=1.5, 'primary mirror (reflectance)'

  xyouts, 364, 0.8575, 'DESI corrector', color=djs_icolor('green'), $
      charsize=1.4

  bitmap = tvrd(true=1)
  write_png, 'ci_throughput_factors.png', bitmap

; all of these multiplicative throughput factors on a common
; wavelength grid are fraction 0->1

  atm_transmission_common = interpol(atm_transmission_frac, $
      lambda_atm_nm, lambda_common_nm)
  primary_reflectance_common = interpol(primary_reflectance_frac, $
      lambda_mirror_nm, lambda_common_nm)
  corrector_transmission_common = interpol(corrector_transmission_frac, $
      lambda_corrector_nm, lambda_common_nm)
  filter_transmission_common = interpol(filter_transmission_frac, $
      lambda_filter_nm, lambda_common_nm)
  qe_common = interpol(qe_frac, lambda_qe_nm, lambda_common_nm)

  oplot, lambda_common_nm, atm_transmission_common, color=djs_icolor('red')
  oplot, lambda_common_nm, primary_reflectance_common, color=djs_icolor('red')
  oplot, lambda_common_nm, corrector_transmission_common, $
      color=djs_icolor('red')
  oplot, lambda_common_nm, filter_transmission_common, color=djs_icolor('red')
  oplot, lambda_common_nm, qe_common, color=djs_icolor('red')

  total_transmission = atm_transmission_common*primary_reflectance_common*corrector_transmission_common*filter_transmission_common*qe_common

  plot, lambda_common_nm, total_transmission, psym=3
; create output structure

  outstr = replicate({lambda_nm: 0.0, nu_hz: 0.0, transmission: 0.0}, $
      n_elements(lambda_common_nm))

  outstr.lambda_nm = lambda_common_nm
  outstr.transmission = total_transmission

  c_air = 299792458.0/1.0003 ; m/s

  nu_hz = c_air/(lambda_common_nm*(1e-9))

  outstr.nu_hz = nu_hz

  if keyword_set(write) then begin
      outname = '../etc/ci_throughput-airmass_' + $
          string(airmass, format='(F4.2)') + '.fits'
      print, outname
      if file_test(outname) then stop
      _outstr = struct_trimtags(outstr, select=['LAMBDA_NM'])
      addstr = replicate({throughput: 0.0d}, n_elements(_outstr))
      addstr.throughput = outstr.transmission
      _outstr = struct_addtags(_outstr, addstr)
      mwrfits, _outstr, outname
  endif

end
