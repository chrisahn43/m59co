PRO TOTALCOLOR

  scale=0.05
  radius=50.
  zeropoint_g=26.05923
  zeropoint_z=24.84245
  gsun=5.11
  zsun=4.52
  fits_read,'/Users/chrisahn/research/code/gemini15/m59co/data/m59_g_lucy.fits',oldimg,head
  hrotate,oldimg,head,img,newhead,1
  find_galaxy,img,m,e,a,xcg,ycg
  readcol,'m59co_mge_outputsersic.dat',gintensity,gsigmaarc,gq,gpa,format='F,F,F,F'
  mge2image,img,xcg,ycg,gintensity,gsigmaarc,gq,gpa,gimg,zeropoint=zeropoint_g,scale=scale,msun=gsun
  
  fits_read,'/Users/chrisahn/research/code/gemini15/m59co/data/m59_z_lucy.fits',oldimg,head
  hrotate,oldimg,head,img,newhead,1
  find_galaxy,img,m,e,a,xcz,ycz
  readcol,'m59co_mge_outputsersic850_free.dat',zintensity,zsigmaarc,zq,zpa,format='F,F,F,F'

  mge2image,img,xcz,ycz,zintensity,zsigmaarc,zq,zpa,zimg,zeropoint=zeropoint_z,scale=scale,msun=zsun

  aper,gimg,xcg,ycg,gflux,fluxerr,0.,skyerr,1,radius,-1,[1,1],/silent,setskyval=0.,/flux,/exact
  aper,zimg,xcz,ycz,zflux,fluxerr,0.,skyerr,1,radius,-1,[1,1],/silent,setskyval=0.,/flux,/exact

  area=!PI*(radius)^2

  gmag=zeropoint_g+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(gflux)
  zmag=zeropoint_z+5*alog10(scale)+2.5*alog10(area)-2.5*alog10(zflux)

  color=gmag-zmag
  print,color
  stop
END
