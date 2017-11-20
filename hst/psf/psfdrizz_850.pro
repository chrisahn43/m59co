PRO psfdrizz_850
  fits_read,'f850lp_flt.fits',zero,header,exten_no=0
  fits_read,'f850lp_flt.fits',img,header_sci,exten_no=1
  fits_read,'f850lp_flt.fits',err,header_err,exten_no=2
  fits_read,'f850lp_flt.fits',dq,header_dq,exten_no=3
  outhead=header
  imgsize=size(img,/dim)
  fits_read,'f850lp_psf.fits',psf
  temp=img[2494:3494,0:500]
  find_galaxy,temp,majorAxis,eps,ang,xc,yc
  find_galaxy,psf,majoraxis1,eps1,ang1,xcpsf,ycpsf
  img[*,*]=0.0
  stop
  img[(xc+2494)-(xcpsf-1):(xc+2494)+xcpsf,0:yc+(ycpsf-1)]=psf[*,45:242]
  outim='psfonimg_850lp.fits'
  fits_write,outim,zero,outhead
  fits_open,outim,fcbout,/update
  fits_write,fcbout,img,header_sci,extname='SCI'
  fits_write,fcbout,err,header_err,extname='ERR'
  fits_write,fcbout,dq,header_dq,extname='DQ'
  fits_close,fcbout
    
  stop
END
