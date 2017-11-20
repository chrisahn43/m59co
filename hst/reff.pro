PRO REFF
  scale=0.05
  fits_read,'/Users/chrisahn/research/code/gemini15/m59co/hst/make_decon/serfreemodelg.fits',gimg
  radius=findgen(100)+1
  find_galaxy,gimg,m,e,a,xc,yc
  gflux=fltarr(n_elements(radius))
  for i=0,n_elements(radius)-1 do begin
     aper,gimg,xc,yc,flux,fluxerr,0.,skyerr,1,radius[i],-1,[1,1],/silent,setskyval=0.,/flux,/exact
     gflux[i]=flux
  endfor
  radius=radius*scale
  djs_plot,radius,gflux,psym=4,ytitle='counts',xtitle='Radius ["]',xran=[0,2.5],/xstyle
  stop
END
