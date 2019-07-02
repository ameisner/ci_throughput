function night_from_expid, expid

  if n_elements(expid) NE 1 then stop

  if (expid NE 4486) AND (expid NE 7577) AND (expid NE 7578) AND $
     (expid NE 7579) AND (expid NE 7580) AND (expid NE 7581) then stop

  if expid EQ 4486 then return, '20190406'

  return, '20190417'

end

pro get_cic_image_and_header, im, h, raw=raw, astr=astr, h_raw=h_raw, $
                              expid=expid, bitmask=bitmask

  if ~keyword_set(expid) then expid = 4486

  if n_elements(expid) NE 1 then stop

  expid_string = string(expid, format='(I08)')

  night = night_from_expid(expid)

  fname_raw = '/project/projectdirs/desi/spectro/data/' + night + $
      '/' + expid_string + '/ci-' + expid_string + '.fits.fz'

  if ~file_test(fname_raw) then stop

; need this to get the airmass
  h_raw = headfits(fname_raw, ex=1, /silent)

  if ~keyword_set(raw) then $
  fname = '/project/projectdirs/desi/users/ameisner/CI/reduced/v0001/' + $
      night + '/ci-' + expid_string + '/ci-' + expid_string + '_reduced.fits' $
  else $
      fname = fname_raw
 
  fname_astr = $
  '/project/projectdirs/desi/users/ameisner/CI/reduced/v0001/' + night + $
  '/ci-' + expid_string + '/ci-' + expid_string + '_astr-a.fits'

  astr = mrdfits(fname_astr, 1)
  astr = astr[where(strtrim(astr.extname, 2) EQ 'CIC')]

  print, 'READING : ' + fname + ' @@@@@@@@@@@@@@@@'

  if ~keyword_set(raw) then ex = 2 else ex = get_raw_extnum('CIC', fname, $
      /image)
  im = readfits(fname, h, ex=ex)
  if keyword_set(raw) then im = float(im)

  mask = readfits('/project/projectdirs/desi/users/ameisner/CI/ci_reduce_etc/CI_static_badpixels.fits', ex=1, hmask)

  if strtrim(sxpar(hmask, 'EXTNAME'), 2) NE 'CIC' then stop
  if strtrim(sxpar(h, 'EXTNAME'), 2) NE 'CIC' then stop

  intx = djs_maskinterp(im, mask NE 0, iaxis=0)
  inty = djs_maskinterp(im, mask NE 0, iaxis=1)

  im = (intx + inty)/2

  bitmask = $
    readfits('/project/projectdirs/desi/users/ameisner/CI/reduced/v0001/' + $
    night + '/ci-' + expid_string + '/ci-' + expid_string + '_bitmask.fits.gz')

end

pro get_ps1_sample, cat, raw=raw, expid=expid

  get_cic_image_and_header, im, h, raw=raw, astr=astr, expid=expid

; use corner coordinates + center coordinates to get PS1 catalog

;           LL    UL       UR       LR       center
  x_vals = [0.0d, 0.0d   , 3071.0d, 3071.0d, 1535.5d]
  y_vals = [0.0d, 2047.0d, 2047.0d, 0.0d   , 1023.5d]

  xy2ad, x_vals, y_vals, astr, ra_vals, dec_vals

  cat = read_ps1cat(ra_vals, dec_vals)

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

pro try_aper_phot, cat, im=im, raw=raw, expid=expid

; cat arg is meant to be OUTPUT
; im also meant as optional output

; assumes x, y input are already refined with try_recentroid !!

  get_cic_image_and_header, im, h, raw=raw, expid=expid

  get_ps1_sample, cat, raw=raw, expid=expid

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

pro append_bitmask_info, cat, bitmask

; cat gets modified via the addition of another bitmask column that's
; just nearest neighbor interpolated off of the bitmask image

  ix = ((long(round(cat.x)) > 0) < 3172)
  iy = ((long(round(cat.y)) > 0) < 2047)

  bitmask_vals = bitmask[ix, iy]

  addstr = replicate({ci_bitmask: 0L}, n_elements(cat))
  addstr.ci_bitmask = bitmask_vals

  cat = struct_addtags(cat, addstr)

end

pro get_zp_e_per_s, raw=raw, expid=expid

  if ~keyword_set(expid) then expid = 4486

  if n_elements(expid) NE 1 then stop

  get_cic_image_and_header, _, h, raw=raw, astr=astr, h_raw=h_raw, $
      expid=expid, bitmask=bitmask

  exptime = sxpar(h, 'EXPTIME')
  gain = 1.64

 ; get cat and im
  try_aper_phot, cat, im=im, raw=raw, expid=expid

  append_bitmask_info, cat, bitmask

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

  calc_ci_zeropoint, airmass=sxpar(h_raw, 'AIRMASS'), zp_pred=zp_pred, /silent
  print, 'PREDICTED zero point is', zp_pred

  !p.multi = [0, 1, 1]
  plot, cat.rmag, inst_mag, psym=1, charsize=2.5, xtitle='r' + $
      textoidl('!Dps1!N'), $
      ytitle='-2.5'+ textoidl('\times') + 'log' + textoidl('!D10!N') + $
       '(e-/s)', title='EXPID = ' + strtrim(string(expid), 2) + $
       ', EXTNAME = CIC'

  oplot, [0, 100], [0, 100]-zp, color=djs_icolor('red')

  xyouts, 14.3, -10.0, 'y = x - ' + strtrim(string(zp, format='(F10.3)') ,2 ), $
      color=djs_icolor('red'), charsize=3

  xyouts, 15.0, -12.0, 'assumed gain = ' + string(gain, format='(F5.3)') + $
      ' e-/ADU', charsize=2, color=djs_icolor('green')

  bitmap = tvrd(true=1)
  write_png, 'zpt_summary_' + strtrim(string(expid), 2) + '_CIC.png', bitmap
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
