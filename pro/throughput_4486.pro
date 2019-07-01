; need to get best possible astrometric solution

function get_as_built_astrom, skyra, skydec, extname=extname

  if ~keyword_set(extname) then extname = 'CIC'

  astrom = mrdfits('/project/projectdirs/desi/users/ameisner/CI/ci_reduce_etc/viewer_astrom_index-as_built.bigtan.fits', 1)
  astr = astrom[where(astrom.extname EQ 'CIC')]

  if n_elements(skyra) NE 1 then stop
  if n_elements(skydec) NE 1 then stop

  if size(skyra, /type) NE 5 then stop
  if size(skydec, /type) NE 5 then stop


  astr.crval = [skyra, skydec]

  return, astr
end

function shifted_astrom, skyra, skydec, extname=extname

  astr_raw = get_as_built_astrom(skyra, skydec, extname=extname)

  astr_shifted = astr_raw

  astr_shifted.crpix[0] += 475.0d
  astr_shifted.crpix[1] -= 57.5d

  return, astr_shifted
end

function get_shifted_astrom, h

  return, shifted_astrom(sxpar(h, 'SKYRA'), sxpar(h, 'SKYDEC'))
end

pro get_cic_image_and_header, im, h, raw=raw, astr=astr

  if ~keyword_set(raw) then $
  fname = $
      '/global/cscratch1/sd/ameisner/real_data3/ci-00004486_reduced.fits.gz' $
  else $
  fname = '/project/projectdirs/desi/spectro/data/20190406/00004486/ci-00004486.fits.fz'
 
  fname_astr = '/project/projectdirs/desi/users/ameisner/CI/reduced/v0001/20190406/ci-00004486/ci-00004486_astr-a.fits'

  astr = mrdfits(fname_astr, 1)
  astr = astr[where(strtrim(astr.extname, 2) EQ 'CIC')]

  print, 'READING : ' + fname + ' @@@@@@@@@@@@@@@@'
  im = readfits(fname, h, ex=2)

  mask = readfits('/project/projectdirs/desi/users/ameisner/CI/ci_reduce_etc/CI_static_badpixels.fits', ex=1, hmask)

  if strtrim(sxpar(hmask, 'EXTNAME'), 2) NE 'CIC' then stop
  if strtrim(sxpar(h, 'EXTNAME'), 2) NE 'CIC' then stop

  intx = djs_maskinterp(im, mask NE 0, iaxis=0)
  inty = djs_maskinterp(im, mask NE 0, iaxis=1)

  im = (intx + inty)/2

end

pro get_ps1_sample, cat, raw=raw

  get_cic_image_and_header, im, h, raw=raw, astr=astr

  cat = read_ps1cat(astr.crval[0], astr.crval[1])

  rmag = reform(cat.median[1,*], n_elements(cat))

  cat = cat[where(rmag LT 17)]

  ad2xy, cat.ra, cat.dec, astr, x, y

  good = (x GT 20) AND (y GT 20) AND (x LT 3051) AND (y LT 2027)

  cat = cat[where(good)]
  x = x[where(good)]
  y = y[where(good)]

  addstr = replicate({x: 0.0d, y: 0.0d}, n_elements(cat))
  addstr.x = x
  addstr.y = y
  cat = struct_addtags(cat, addstr)

  rmag = reform(cat.median[1,*], n_elements(cat))

  addstr = replicate({rmag: 0.0}, n_elements(cat))

  addstr.rmag = rmag

  cat = struct_addtags(cat, addstr)

end

;pro optimize_astrom, astr_guess

;  get_ps1_sample, cat

  

;end

pro try_recentroid, x, y, im

  for i=0L, n_elements(x)-1 do begin
         xcen = x[i]
         ycen = y[i]
         djs_photcen, xcen, ycen, im, $
             cbox=9, cmaxshift=5.0
;    [ calg=, cbox=, cmaxiter=, cmaxshift=, fwhm=, /fixfw, ceps=, qmaxshift= ]
         x[i] = xcen
         y[i] = ycen
  endfor

end

pro try_aper_phot, cat, im=im, raw=raw

; cat arg is meant to be OUTPUT
; im also meant as optional output

; assumes x, y input are already refined with try_recentroid !!

  get_cic_image_and_header, im, h, raw=raw

  get_ps1_sample, cat, raw=raw

  x = cat.x
  y = cat.y 

  try_recentroid, x, y, im

  fluxes = []
  skyvals =[]
  for i=0L, n_elements(x)-1 do begin
      xcen = x[i]
      ycen = y[i]
      flux = djs_phot(xcen, ycen, 5.0, [50.0, 55.0], im, $
        calg='none', cbox=9, cmaxshift=5.0, $
        skyval=skyval, peakval=peakval)
      fluxes = [fluxes, flux]
      skyvals = [skyvals, skyval]
  endfor

  print, minmax(fluxes)

  plot, cat.rmag, -2.5*alog10(fluxes), psym=1

  oplot, [0,100], [0,100]-27.867

;  plothist, skyvals

  cat.x = x ; put in the recentered centroid coords
  cat.y = y

  addstr = replicate({skyval: 0.0, flux: 0.0}, n_elements(cat))

  addstr.skyval = skyvals
  addstr.flux = fluxes

  cat = struct_addtags(cat, addstr)

end

function min_edge_dist_ci, x_pix, y_pix

  min_edge_dist = 10000L

  min_edge_dist <= (x_pix + 0.5)
  min_edge_dist <= (y_pix + 0.5)
  min_edge_dist <= (3071.5 - x_pix)
  min_edge_dist <= (2047.5 - y_pix)

  return, min_edge_dist

end

function get_cutouts, im, cat, sidelen=sidelen

; think i need to go bigger eventually
  if ~keyword_set(sidelen) then sidelen = 101

; require odd
  if (sidelen MOD 2) EQ 0 then stop

  _edge_dist = min_edge_dist_ci(cat.x, cat.y)

  half = sidelen/2

  _cat = cat[where(_edge_dist GT half)]

  cube = fltarr(sidelen, sidelen, n_elements(_cat))
  for i=0L, n_elements(_cat)-1 do begin
      xmax = long(round(_cat[i].x)) + half
      xmin = long(round(_cat[i].x)) - half
      ymax = long(round(_cat[i].y)) + half
      ymin = long(round(_cat[i].y)) - half

      cutout = im[xmin:xmax, ymin:ymax]
; STILL NEED TO SSHIFT INTO ALIGNMENT !!!
      cutout -= _cat[i].skyval
      cutout = sshift2d(cutout, [long(round(_cat[i].x)) - _cat[i].x, $
                                 long(round(_cat[i].y)) - _cat[i].y])
      ; normalize using PS1 r mag, standardize to 15th mag arbitrarily
      cutout *= 10^((_cat[i].rmag-15.0)/2.5)
      cube[*, *, i] = cutout
  endfor

  return, cube
end

function get_psf, cube

  return, median(cube, dim=3)
end

function get_aperture_corr, psf

  sz = (size(psf, /dim))[0]
  if (sz[0] MOD 2) EQ 0 then stop

  xcen = double(sz[0]/2)
  ycen = xcen

  flux_aper = djs_phot(xcen, ycen, 5.0, [50.0, 55.0], psf, $
        calg='none', $
        skyval=skyval)

  flux_tot = djs_phot(xcen, ycen, 50.0, [50.0, 55.0], psf, $
        calg='none', $
        skyval=skyval)

  print, flux_aper
  print, flux_tot

  rat = flux_tot/flux_aper

  return, rat[0]
end

pro get_zp_e_per_s, raw=raw

  get_cic_image_and_header, _, h, raw=raw

  exptime = sxpar(h, 'EXPTIME')
  gain = 1.64

 ; get cat and im
  try_aper_phot, cat, im=im, raw=raw

  cube =  get_cutouts(im, cat, sidelen=101)

  psf = get_psf(cube)
  aper_corr = get_aperture_corr(psf)
  flux_tot_adu = cat.flux*aper_corr

help, aper_corr
  flux_tot_e = flux_tot_adu*gain

  rate_tot_e_per_s = flux_tot_e/exptime

  inst_mag = -2.5*alog10(rate_tot_e_per_s)

;  plot, cat.rmag, inst_mag, psym=1


  zp = median(cat.rmag - inst_mag)
  print, 'zero point is', zp

  plot, cat.rmag, inst_mag, psym=1, charsize=2.5, xtitle='r' + $
      textoidl('!Dps1!N'), $
      ytitle='-2.5'+ textoidl('\times') + 'log' + textoidl('!D10!N') + $
       '(e-/s)', title='EXPID = 4486, EXTNAME = CIC'

  oplot, [0, 100], [0, 100]-zp, color=djs_icolor('red')

  xyouts, 14.3, -10.0, 'y = x - ' + strtrim(string(zp, format='(F10.3)') ,2 ), $
      color=djs_icolor('red'), charsize=3

  xyouts, 15.0, -12.0, 'assumed gain = ' + string(gain, format='(F5.3)') + $
      ' e-/ADU', charsize=2, color=djs_icolor('green')

  bitmap = tvrd(true=1)
  write_png, 'zpt_summary_4486_CIC.png', bitmap
stop
; remember to compare to value that incorporates correct airmass !!
; also there's the aperture mask issue, which should cause real
; measurement to be deeper than forecasted by a little bit

end

; https://arxiv.org/abs/1303.3634 says saturation in PS1 r happens at
; r ~ 13.5
; seems like a good range of mags
; ww = where((gaia.phot_g_mean_mag GT 14.5) AND (gaia.phot_g_mean_mag LT 17.0))

; airmass = 1.5908490000000000

pro calc_pred_zp

  h = headfits('/project/projectdirs/desi/spectro/data/20190406/00004486/ci-00004486.fits.fz',ex=1)

  calc_ci_zeropoint, airmass=sxpar(h, 'AIRMASS')
end
