; don't forget vignetting !!!!!!!!!!!!!!!

pro assemble_gfa_transmission, airmass=airmass, outstr=outstr, write=write

; common grid onto which to interpolate various inputs
; GFA filter transmission remains nonzero until very close to 750 nm
; so it might be good to extend this a little bit redder, but don't worry
; about it for now ...
  lambda_common_nm = 500.0 + findgen(251.0)

;  atmosphere*aluminum*corrector*filter*QE*VIGNETTING

  vignetting_fac = 0.95910780 ; updated 2019 December 16, get_vignetting_fac.pro
  
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

    readcol, '../etc/gfa_filter_transmission_DESI-1297.dat', $
      lambda_filter_nm, filter_transmission_percent_gfa_1, $
                        filter_transmission_percent_gfa_2, $
                        filter_transmission_percent_gfa_3, F='F, F, F, F'


; the three GFA filter transmission curves are virtually indistinguishable

  filter_transmission_percent = (filter_transmission_percent_gfa_1 + $
                                filter_transmission_percent_gfa_2 + $
                                filter_transmission_percent_gfa_3)/3.0

  filter_transmission_frac = filter_transmission_percent/100.0

  readcol, '../etc/gfa_qe_midband.csv', lambda_qe_nm, qe_percent, F='F, F'

  qe_frac = qe_percent/100.0
  
  window, xsize=1200, ysize=650

    dfpsplot, 'gfa_throughput_factors.eps', /encap, $
              /color, ysize=8, xsize=15
    
  xtitle = 'wavelength (nm)'
  ytitle = 'throughput'
  plot, lambda_mirror_nm, primary_reflectance_frac, xrange=[350, 1055], /xst, $
        charsize=1.4, xtitle=xtitle, ytitle=ytitle, $
        title='DESI GFA throughput multiplicative factors', thick=2, charthick=2
  
  oplot, [-1000, 2000], vignetting_fac + [0, 0], color=djs_icolor('magenta'), $
           linestyle=1, thick=2
  oplot, lambda_atm_nm, atm_transmission_frac, color=djs_icolor('blue'), $
         thick=2
  oplot, lambda_corrector_nm, corrector_transmission_frac, $
      color=djs_icolor('green'), thick=2
  oplot, lambda_filter_nm, filter_transmission_frac, color=djs_icolor('red'), $
         thick=2
  oplot, lambda_qe_nm, qe_frac, color=djs_icolor('orange'), thick=2

  charthick = 3
  xyouts, 775 + 27.5, 0.65 + 0.12, 'KPNO', charsize=2, $
          color=djs_icolor('blue'), $
          charthick=charthick
  xyouts, 775 + 27.5, 0.6 + 0.12, 'atmosphere', charsize=2, $
          color=djs_icolor('blue'), $
          charthick=charthick
  xyouts, 775 + 27.5, 0.55 + 0.12, 'airmass = ' + $
          string(airmass, format='(F4.2)'), $
          charsize=2, color=djs_icolor('blue'), charthick=charthick

  xyouts, 880, 0.425, 'e2v CCD230-42 QE', charsize=2, $
      color=djs_icolor('orange'), charthick=charthick

  xyouts, 898, 0.375, '(midband coating)', charsize=2, $
          color=djs_icolor('orange'), charthick=charthick
  
  xyouts, 467.5+25, 0.15, 'r filter', charsize=2, $
      color=djs_icolor('red'), charthick=charthick

  xyouts, 440.0+25, 0.10, '(DESI-1297)', charsize=2, $
      color=djs_icolor('red'), charthick=charthick
  
  xyouts, 370, 0.90, charsize=1.5, 'primary mirror (reflectance)', $
          charthick=charthick

  xyouts, 390, 0.962, charsize=1.7, 'vignetting', color=djs_icolor('magenta'), $
          charthick=charthick
  
  xyouts, 364, 0.8575 - 0.0065, 'DESI corrector', color=djs_icolor('green'), $
      charsize=1.4, charthick=charthick

  dfpsclose
;  bitmap = tvrd(true=1)
;  write_png, 'gfa_throughput_factors.png', bitmap

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

  total_transmission = atm_transmission_common*primary_reflectance_common*corrector_transmission_common*filter_transmission_common*qe_common*vignetting_fac

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
      outname = '../etc/gfa_throughput-airmass_' + $
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
