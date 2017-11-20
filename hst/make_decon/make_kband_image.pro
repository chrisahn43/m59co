PRO MAKE_KBAND_IMAGE

  zeropoint_g_ab=26.05923
  zeropoint_z_ab=24.84245
  zeropoint_g_ve=26.16252
  zeropoint_z_ve=24.32305
  extinct=0.27
  readcol,'./johnson_ssp.dat',johnz,johnage,johnmbol,u,b,v,r,i,j,h,k,format='F,F,F,F,F,F,F,F,F,F,F'
  readcol,'./hstwfc_ssp.dat',hstz,hstage,hstmbol,f435w,f475w,f502n,f550m,f555w,f606w,f625w,f658n,f660n,f775w,f814w,f850lp,f892n,format='F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'
  kext=extinct*0.11471
  k=k+kext
  f475wext=extinct*1.19119
  f475w=f475w+f475wext
  f850ext=extinct*0.4844
  f850lp=f850lp+f850ext
  ind=where(hstage gt 1.e9)
  f475w=f475w[ind]
  f850lp=f850lp[ind]
  k=k[ind]
  f475wab=(f475w+zeropoint_g_ab-zeropoint_g_ve)
  f850lpab=(f850lp+zeropoint_z_ab-zeropoint_z_ve)
  color=f475wab-f850lpab
  kcolor=f850lpab-k
  djs_plot,color,kcolor,psym=4
  stop

  fits_read,'serfixmodelg.fits',gimg
  fits_read,'serfreemodelz.fits',zimg
  stop
  zimg=shift(zimg,0,1)
  gimg=gimg[1:*,1:*]
  zimg=zimg[1:*,1:*]
  imgsize=size(gimg,/dim)
  gmag=fltarr(imgsize[0],imgsize[1])
  zmag=fltarr(imgsize[0],imgsize[1])

  for i=0,imgsize[0]-1 do begin
     for j=0,imgsize[1]-1 do begin
        gmag[i,j]=zeropoint_g_ab-2.5*alog10(gimg[i,j])
        zmag[i,j]=zeropoint_z_ab-2.5*alog10(zimg[i,j])
     endfor
  endfor
  diff=gmag-zmag
  kimg=fltarr(imgsize[0],imgsize[1])

  for i=0,imgsize[0]-1 do begin
     for j=0,imgsize[1]-1 do begin
        cind=where(diff[i,j] gt color-0.017 and diff[i,j] lt color+0.017,count)
        if (count eq 1) then begin
           z_k=kcolor[cind]
           kimg[i,j]=zimg[i,j]*(10^(-0.4*z_k))
        endif
        if (count gt 1) then begin
           z_k=mean(kcolor[cind])
           kimg[i,j]=zimg[i,j]*(10^(-0.4*z_k))
        endif
     endfor
  endfor
  writefits,'kbandimage.fits',kimg
  stop
END
