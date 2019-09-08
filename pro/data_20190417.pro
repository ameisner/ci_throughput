pro data_20190417

  expid = 7577 + lindgen(5)

  flist = '/project/projectdirs/desi/spectro/data/20190417/' + $
      string(expid, format='(I08)') + '/ci-' + $
      string(expid, format='(I08)') + '.fits.fz'

  print, n_elements(unique(flist))

  print, total(file_test(flist))

  for i=0L, n_elements(flist)-1 do begin
      h = headfits(flist[i], ex=1, /silent)
      print, sxpar(h, 'EXPID'), '  ~~  ', sxpar(h, 'EXPTIME'), ' ~~ ', $
          sxpar(h, 'AIRMASS')
  endfor

end
