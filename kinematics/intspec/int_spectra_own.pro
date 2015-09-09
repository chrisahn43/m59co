pro int_spectra_own,minradius,maxradius,xcen,ycen

  IF (NOT KEYWORD_SET(minradius)) then minradius=0.
  IF (NOT KEYWORD_SET(maxradius)) then maxradius=1.

  infile='../vucd3_combine_best8.fits'
  outfile='vucd3_combine_best8_int_fc'+STRTRIM(FIX(minradius),2)+'-'+STRTRIM(FIX(maxradius),2)+'.fits'

  cube=DOUBLE(READFITS(infile,head,ext=1,/SILENT))
  varcube=DOUBLE(READFITS(infile,head,ext=2,/SILENT))
  cubesize=SIZE(cube,/dim)
  lambda0=SXPAR(head,'CRVAL3')
  dlambda=SXPAR(head,'CD3_3')
  lambda=FINDGEN(cubesize[2])*dlambda+lambda0
  outspec=DBLARR(cubesize[2])
  outrawspec=DBLARR(cubesize[2])
  outvar=DBLARR(cubesize[2])
  outsky=DBLARR(cubesize[2])

END
