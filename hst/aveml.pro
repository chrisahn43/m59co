PRO AVEML
  scale=0.05
  zeropoint=26.05923
  msun=5.10
  mlgin=2.84950d;2.26893d;
  mlgout=5.47835d;4.94259d;
  rad=9.497
  readcol,'m59co_mge_outputsersic.dat',intensity,sigmaarc,q,pa,format='F,F,F,F'
  fits_read,'/Users/chrisahn/research/code/gemini15/m59co/data/m59_g_lucy.fits',oldimg,head
  stop
  hrotate,oldimg,head,img,newhead,1
  find_galaxy,img,m,e,a,xc,yc
  inind=where(pa gt 20.)
  inintensity=intensity[inind]
  insigma=sigmaarc[inind]
  inq=q[inind]
  inpa=pa[inind]
  mge2image,img,xc,yc,inintensity,insigma,inq,inpa,inmodel,zeropoint=zeropoint,scale=scale,msun=msun

  outind=where(pa lt 20.)
  outintensity=intensity[outind]
  outsigma=sigmaarc[outind]
  outq=q[outind]
  outpa=pa[outind]
  mge2image,img,xc,yc,outintensity,outsigma,outq,outpa,outmodel,zeropoint=zeropoint,scale=scale,msun=msun

  aper,inmodel,xc,yc,influx,influxerr,0.,skyerr,1,50.,-1,[1,1],/silent,setskyval=0.,/flux,/exact
  aper,outmodel,xc,yc,outflux,outfluxerr,0.,skyerr,1,50.,-1,[1,1],/silent,setskyval=0.,/flux,/exact
  inmag=zeropoint-2.5*alog10(influx)
  inabsmag=inmag-5*(alog10(16.5e6)-1)
  inlum=10^((inabsmag-msun)/(-2.5))
  outmag=zeropoint-2.5*alog10(outflux)
  outabsmag=outmag-5*(alog10(16.5e6)-1)
  outlum=10^((outabsmag-msun)/(-2.5))
  aveml=((influx*mlgin)+(outflux*mlgout))/(influx+outflux)
  print,aveml
  gvcolorin=0.4216
  gvcolorout=0.4814
  avecolor=((influx*gvcolorin)+(outflux*gvcolorout))/(influx+outflux)
  stop
  
END
