@myfunct.pro
pro mask_skysubtract
  infile='../data/HST_9401_09_ACS_WFC_F475W_drz.fits'
  fits_read,infile,img,h
  imgsize=size(img,/dim)
  fits_read,infile,mask,exten_no=3
  mdrizz=[39.2494828065,39.0788202159]
  expt=375.
  mdrizzcount=mdrizz/expt
  avgdrizz=mean(mdrizzcount)
  img=img+avgdrizz
  temp=img[4000:4630,3700:4296]
;  center=gauss2dfit(temp,A)
;  xcen=A[4]+4270.
;  ycen=A[5]+3870.
  xcen=4315 & ycen=3998
  r=radial_dist(imgsize[0],imgsize[1],xcen,ycen)
  for i=0,imgsize[0]-1 do begin
     for j=0,imgsize[1]-1 do begin
        if r[i,j] le 100. then mask[i,j]=0
     endfor
  endfor
  xcenstar1=4053 & ycenstar1=3833
  rstar1=radial_dist(imgsize[0],imgsize[1],xcenstar1,ycenstar1)
  for i=0,imgsize[0]-1 do begin
     for j=0,imgsize[1]-1 do begin
        if rstar1[i,j] le 10. then mask[i,j]=0
     endfor
  endfor
  xcenstar2=4338 & ycenstar2=4127
  rstar2=radial_dist(imgsize[0],imgsize[1],xcenstar2,ycenstar2)
  for i=0,imgsize[0]-1 do begin
     for j=0,imgsize[1]-1 do begin
        if rstar2[i,j] le 10. then mask[i,j]=0
     endfor
  endfor
  xcenstar3=4156 & ycenstar3=3704
  rstar3=radial_dist(imgsize[0],imgsize[1],xcenstar3,ycenstar3)
  for i=0,imgsize[0]-1 do begin
     for j=0,imgsize[1]-1 do begin
        if rstar3[i,j] le 10. then mask[i,j]=0
     endfor
  endfor
  newmask=fltarr(imgsize[0],imgsize[1])
  temperr=img*750.
  for i=0,imgsize[0]-1 do begin
     for j=0,imgsize[1]-1 do begin
        if (mask[i,j] eq 5 and temperr[i,j] gt 0.) then newmask[i,j]=1./sqrt(temperr[i,j]) else newmask[i,j]=0.
     endfor
  endfor
  newimg=img[4000:4630,3700:4296]
  newmask=newmask[4000:4630,3700:4296]

  x = (findgen(631)) # (fltarr(597)+1)
  y = (fltarr(631)+1) # (findgen(597))
  z=newimg
  p0=[0.028,1.e-6,1.e-6]
  aa=mpfit2dfun('myfunct',x,y,z,err,p0,weights=newmask,covar=covar,perror=perror)
  plane=aa[0]+aa[1]*x+aa[2]*y
  
  stop
  writefits,'sky_mask.fits',plane
  stop
END
