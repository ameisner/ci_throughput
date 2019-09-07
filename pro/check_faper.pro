pro check_faper, psf

  ; Paul's CI ETC sheet says that, for an aperture size 
  ; with diameter equal to the FWHM, the fraction of flux counted
  ; is 0.5

  ; just curious whether i can reproduce that value approximately; certainly
  ; the exact answer will depend on the shape of the PSF

  ; i will try it with a gaussian psf

  ; use a very oversampled PSF to avoid getting burned by fractional
  ; pixel issues

  psf = psf_gaussian(npixel=301, fwhm=10.0)

  psf = psf/total(psf)

end
