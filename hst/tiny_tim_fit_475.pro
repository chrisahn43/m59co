PRO TINY_TIM_FIT_475

  fits_read,'./psf/temp_psf.fits',img,h
  scale=0.05
  ngauss=10
  minlevel=1.e-40
  
  find_galaxy,img,majoraxis,eps,ang,xc,yc,fraction=0.3,/plot

  sectors_photometry,img,eps,ang,xc,yc,radius,angle,counts,minlevel=minlevel

  mge_fit_sectors,radius,angle,counts,eps,sol=sol,ngauss=ngauss,qbound=[0.9999999,1.],scale=scale

  modelnorm=1./total(sol[0,*])
  modelpeak=modelnorm*(sol[0,*])
  modelsig=sol[1,*]
  forprint,modelpeak,modelsig,sol[2,*],format='F,F,F',textout='tinytim_fit_475.dat'
  stop
END
