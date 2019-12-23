pro ingest_fiducial_spectrum, outstr

  fname = '../etc/fiducial-spectrum-25401.txt'

  readcol, fname, lambda_aa, flux, F='F, F'

  outstr = replicate({lambda_aa: 0.0, flux: 0.0}, n_elements(flux))
  outstr.lambda_aa = lambda_aa
  outstr.flux = flux

end
