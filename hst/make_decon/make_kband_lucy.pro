PRO MAKE_KBAND_LUCY

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
;  stop

  fits_read,'../../data/m59_g_lucy.fits',gimg
  gimg=gimg/375.
  fits_read,'../../data/m59_z_lucy.fits',zimg
  zimg=zimg/560.
  imgsize=size(gimg,/dim)
  xarr=(findgen(imgsize[0])) # (fltarr(imgsize[1])+1)
  yarr=(fltarr(imgsize[0])+1) # (findgen(imgsize[1]))
  noise=fltarr(imgsize[0],imgsize[1])
  noise[*,*]=sqrt(0.00206449)
  gmag=fltarr(imgsize[0],imgsize[1])
  zmag=fltarr(imgsize[0],imgsize[1])
  voronoi_2d_binning,xarr,yarr,gimg,noise,25,binnumber,/plot
  newgimg=fltarr(imgsize[0],imgsize[1])
  for i=0,n_elements(binnumber)-1 do begin
     ind=where(binnumber eq i,nind)
     xind=ind MOD imgsize[0]
     yind=ind / imgsize[0]
     for j=0,nind-1 do begin
        newgimg[xind[j],yind[j]]=gimg[xind[j],yind[j]]
     endfor
  endfor
  imgsize=size(zimg,/dim)
  xarr=(findgen(imgsize[0])) # (fltarr(imgsize[1])+1)
  yarr=(fltarr(imgsize[0])+1) # (findgen(imgsize[1]))
  noise=fltarr(imgsize[0],imgsize[1])
  noise[*,*]=sqrt(0.00200689)
  voronoi_2d_binning,xarr,yarr,zimg,noise,25,binnumber,/plot
  newzimg=fltarr(imgsize[0],imgsize[1])
  for i=0,n_elements(binnumber)-1 do begin
     ind=where(binnumber eq i,nind)
     xind=ind MOD imgsize[0]
     yind=ind / imgsize[0]
     for j=0,nind-1 do begin
        newzimg[xind[j],yind[j]]=gimg[xind[j],yind[j]]
     endfor
  endfor
  gimg=newgimg
  zimg=newzimg
  stop
  
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
        cind=where(diff[i,j] gt color-0.57 and diff[i,j] lt color+0.57,count)
        if (count eq 1) then begin
           z_k=kcolor[cind]
           kimg[i,j]=zimg[i,j]*(10^(-0.4*z_k))
        endif
        if (count gt 1) then begin
           z_k=mean(kcolor[cind])
           kimg[i,j]=zimg[i,j]*(10^(-0.4*z_k))
        endif
;        if (count eq 0) then stop
     endfor
  endfor
  writefits,'kbandimage_lucy.fits',kimg
  stop


END
