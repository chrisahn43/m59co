@sersic2image
@mge2image
PRO MAKE_DECON_IMAGE
  zeropoint=26.05923
  zeropointz=24.84245
  scale=0.05d
  Msun=5.11d
  Msunz=4.52
  fits_read,'/Users/chrisahn/research/code/gemini15/m59co/hst/galfit_output/m59co_475_galfit_model.fits',img,exten_no=1
;  fits_read,'/Users/chrisahn/research/code/gemini15/m59co/data/m59_g_lucy.fits',oldimg,head
;  hrotate,oldimg,head,img,newhead,1
  fits_read,'/Users/chrisahn/research/code/gemini15/m59co/hst/galfit_output/m59co_850_galfit_model.fits',imgz,exten_no=1
;  fits_read,'/Users/chrisahn/research/code/gemini15/m59co/data/m59_z_lucy.fits',oldimgz,headz
;  hrotate,oldimgz,headz,imgz,newheadz,1
  imsizeg=size(img,/dim)
  imsizez=size(imgz,/dim)

  find_galaxy,img,m,e,a,xcg,ycg
  find_galaxy,imgz,m,e,a,xcz,ycz

                                ;CREATE G BAND FIXED MGE MODEL IMAGE
  readcol,'/Users/chrisahn/research/code/gemini15/m59co/hst/m59co_mge_outputsersic.dat',fixlum,fixsig,fixq,fixpa,format='D,D,D,D'
  mge2image,img,xcg,ycg,fixlum,fixsig,fixq,fixpa,model,zeropoint=zeropoint,scale=scale,msun=msun,fits='mgefixmodelg.fits'

                                ;CREATE G BAND FIXED SERSIC MODEL IMAGE
  mtot=[19.56,18.32]
  re=[2.96,12.25]
  n=[1.02,1.21]
  q=[0.99,0.98]
  pa=[34.06,17.69]
  sersic2image,img,xcg,ycg,mtot,re,n,q,pa,model,zeropoint=zeropoint,fits='serfixmodelg.fits'
                                ;CREATE G BAND FREE MGE MODEL IMAGE
  
  readcol,'/Users/chrisahn/research/code/gemini15/m59co/hst/m59co_mge_outputsersic_free.dat',freelum,freesig,freeq,freepa,format='D,D,D,D'
  mge2image,img,xcg,ycg,freelum,freesig,freeq,freepa,model,zeropoint=zeropoint,scale=scale,msun=msun,fits='mgefreemodelg.fits'

                                ;CREATE G BAND FREE SERSIC MODEL IMAGE
  mtot=[19.4,18.38]
  re=[3.13,12.85]
  n=[1.06,1.09]
  q=[0.97,0.98]
  pa=[-65.24,88.42]
  sersic2image,img,xcg,ycg,mtot,re,n,q,pa,model,zeropoint=zeropoint,fits='serfreemodelg.fits'

                                ;CREATE Z BAND FREE MGE MODEL
  readcol,'/Users/chrisahn/research/code/gemini15/m59co/hst/m59co_mge_outputsersic850_free.dat',freelum,freesig,freeq,freepa,format='D,D,D,D'
  mge2image,imgz,xcz,ycz,freelum,freesig,freeq,freepa,model,zeropoint=zeropointz,scale=scale,msun=msunz,fits='mgefreemodelz.fits'
  
                                ;CREATE Z BAND FREE SERSIC MODEL
  mtot=[18.17,16.72]
  re=[2.96,12.25]
  n=[1.02,1.21]
  q=[0.99,0.98]
  pa=[34.06,17.69]
  sersic2image,imgz,xcz,ycz,mtot,re,n,q,pa,model,zeropoint=zeropointz,fits='serfreemodelz.fits'

                                ;CREATE Z BAND FIX MGE MODEL
  readcol,'/Users/chrisahn/research/code/gemini15/m59co/hst/m59co_mge_outputsersic850.dat',fixlum,fixsig,fixq,fixpa,format='D,D,D,D'
  mge2image,imgz,xcz,ycz,fixlum,fixsig,fixq,fixpa,model,zeropoint=zeropointz,scale=scale,msun=msunz,fits='mgefixmodelz.fits'

                                ;CREATE Z BAND FIX SERSIC MODEL
  mtot=[17.98,16.78]
  re=[3.13,12.85]
  n=[1.06,1.09]
  q=[0.97,0.98]
  pa=[-65.24,88.42]
  sersic2image,imgz,xcz,ycz,mtot,re,n,q,pa,model,zeropoint=zeropointz,fits='serfixmodelz.fits'
  
  
  
END
