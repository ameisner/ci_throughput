pro frac_area_aper_mask

;  http://www-kpno.kpno.noao.edu/kpno-misc/mayall_params.html

; Central obscuration		diam 1800 mm (DESI cage)

; Used clear aperture		diam 3797 mm (Measured by Phil Massey)

; Bare mirror clear aperture diam 4002.6 mm

  a_no_mask = !pi*(4.0026^2 - 1.8^2)/4
  a_with_mask = !pi*(3.797^2 - 1.8^2)/4

  print, a_no_mask/a_with_mask

end
