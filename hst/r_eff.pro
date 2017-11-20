PRO R_EFF
  scale=0.05
  infile='m59co_475_skysubtract.fits'
  fits_read,infile,img,head
  find_galaxy,img,m,e,a,xc,yc
  xc=316 & yc=299
  rad=findgen(100)+1
  flux=fltarr(n_elements(rad))
  for i=0,n_elements(rad)-1 do begin
     aper,img,xc,yc,maxflux,fluxerr,0.,skyerr,1,rad[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     flux[i]=maxflux
  endfor
  djs_plot,rad*scale,flux,psym=4
  stop
END
