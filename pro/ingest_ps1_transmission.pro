pro ingest_ps1_transmission, outstr, write=write

  skipline = 1

; does the PS1 transmission curve that i downloaded have some assumption
; about the airmass baked in?
  
  fname = '../etc/PAN-STARRS_PS1.r.dat'
  readcol, fname, lambda_aa, t, F='F, F', skipline=skipline

  outstr = replicate({lambda_nm: 0.0, throughput: 0.0}, n_elements(t))

  outstr.lambda_nm = lambda_aa/10.0
  outstr.throughput = t

  if keyword_set(write) then begin
      outname = repstr(fname, '.dat', '.fits')
      if file_test(outname) then stop
      mwrfits, outstr, outname
  endif
  
end
