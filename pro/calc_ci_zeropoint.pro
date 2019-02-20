pro calc_ci_zeropoint, airmass=airmass

  assemble_ci_transmission, airmass=airmass, outstr=outstr

  ergs_per_second_per_hertz_per_cm2 = 3.631e-20 ; mag = 0 AB

  detected_ergs_per_second_per_hertz_per_cm2 = $
      outstr.transmission*ergs_per_second_per_hertz_per_cm2

!P.multi= [0, 2, 2]
  plot, outstr.lambda_nm, detected_ergs_per_second_per_hertz_per_cm2, charsize=2

  h_ergs_by_seconds = 6.62607015e-27 ; Planck constant in erg * s

; ergs per photon
  photon_energy_ergs = outstr.nu_hz*h_ergs_by_seconds

  detected_photons_per_second_per_hertz_per_cm2 = $
      detected_ergs_per_second_per_hertz_per_cm2/photon_energy_ergs

  plot, outstr.lambda_nm, detected_photons_per_second_per_hertz_per_cm2, $
      charsize=2

  area_sq_meters = 8.658739421 ; desi.yaml
  area_sq_cm = area_sq_meters*100.0*100.0

  detected_photons_per_second_per_hertz = $
      area_sq_cm*detected_photons_per_second_per_hertz_per_cm2 

  plot, outstr.lambda_nm, detected_photons_per_second_per_hertz

  c_air = 299792458.0/1.0003 ; m/s
  d_nu_hz = (1e-9)*c_air/((outstr.lambda_nm*(1e-9))^2)

  detected_photons_per_second = $
      total(d_nu_hz*detected_photons_per_second_per_hertz)

  print, 'detected photons per second from flat-spectrum (in f_nu) 0 mag AB', $
      detected_photons_per_second

  print, 'implied zeropoint (AB mag of 1 electron/second source)', $
      2.5*alog10(detected_photons_per_second)

end
