function get_raw_extnum, extname, fname, image=image

; this assumes file is fpacked in the same way as the CI raw data,
; specifically the beginning / end indices of the for loop below

  if n_elements(fname) NE 1 then stop

  if ~file_test(fname) then stop

  if n_elements(extname) NE 1 then stop
  if (extname NE 'CIE') AND (extname NE 'CIN') AND $
     (extname NE 'CIC') AND (extname NE 'CIS') AND $
     (extname NE 'CIW') then stop

  fits_info, fname, n_ext=n_ext, /silent

  for i=2L, n_ext do begin
      h = headfits(fname, ex=i)
      if strtrim(sxpar(h, 'EXTNAME'), 2) EQ strtrim(extname, 2) then begin
          result = keyword_set(image) ? (i-1) : i
          return, result
      endif
  endfor

  return, -1

end

pro _test

  fname = '/project/projectdirs/desi/spectro/data/20190406/00004486/ci-00004486.fits.fz'

  fname = '/project/projectdirs/desi/spectro/data/20190417/00007577/ci-00007577.fits.fz'

  fname = '/project/projectdirs/desi/spectro/data/20190417/00007578/ci-00007578.fits.fz'

  fname = '/project/projectdirs/desi/spectro/data/20190417/00007579/ci-00007579.fits.fz'

  ex = get_raw_extnum('CIC', fname, /image)
  
  print, ex
  _ = readfits(fname, ex=ex, h)

  print, sxpar(h, 'EXTNAME')

end
