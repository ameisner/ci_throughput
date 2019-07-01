pro compare_4486_7577_images

  ; just want to see how close to perfectly aligned the original 4486
  ; CIC image and the newer 7577 image from arjun later in april 
  ; actually are

  old = readfits('/project/projectdirs/desi/users/ameisner/CI/reduced/v0001/20190406/ci-00004486/ci-00004486_reduced.fits', ex=2, h_old)

  if strtrim(sxpar(h_old, 'EXTNAME'), 2) NE 'CIC' then stop

  new = readfits('/project/projectdirs/desi/users/ameisner/CI/reduced/v0001/20190417/ci-00007577/ci-00007577_reduced.fits', ex=2, h_new)

  if strtrim(sxpar(h_new, 'EXTNAME'), 2) NE 'CIC' then stop

  print, 'got here'

  stop

end

; bright star at (x, y) = (2381, 544) in expid = 7577 is at
; (x,y) = (2209, 436), so they're shifted by approx 203.1 pixels
; which is approx 0.13*203.1 asec = 26.4 asec = 0.44 arcmin

; 7577-7581 inclusive are all very nearly perfectly lined up with one
; another (probably down to the level of any tracking imperfections
; since i'm guessing that the telescope was not intentionally moved)

pro compare_4486_7577_astrom

; does 4486 have a good astrometric solution ? i'd hope/assue so !
;  --> yes CIC of 4486 has contrast = 6.8118916

   astr_old = mrdfits('/project/projectdirs/desi/users/ameisner/CI/reduced/v0001/20190406/ci-00004486/ci-00004486_astr-a.fits', 1)

;;;  astr_old = mrdfits('/project/projectdirs/desi/users/ameisner/CI/reduced/v0001/20190417/ci-00007581/ci-00007581_astr-a.fits', 1)

  astr_old = astr_old[where(strtrim(astr_old.extname, 2) EQ 'CIC')]

  astr_new = mrdfits('/project/projectdirs/desi/users/ameisner/CI/reduced/v0001/20190417/ci-00007577/ci-00007577_astr-a.fits', 1)

  astr_new = astr_new[where(strtrim(astr_new.extname, 2) EQ 'CIC')]

  xcen = 1535.5000d
  ycen = 1023.5000d

  xy2ad, xcen, ycen, astr_old, ra_old, dec_old
  xy2ad, xcen, ycen, astr_new, ra_new, dec_new

  dangle = djs_diff_angle(ra_old, dec_old, ra_new, dec_new)
 

  print, dangle*3600.0d

end
