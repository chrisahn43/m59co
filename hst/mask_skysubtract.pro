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
  temp=img[4270:4360,3870:4110]
  center=gauss2dfit(temp,A)
  xcen=A[4]+4270.
  ycen=A[5]+3870.
  r=radial_dist(imgsize[0],imgsize[1],xcen,ycen)
  dist_ellipse,ell,[imgsize[0],imgsize[1]],xcen,ycen,2.3,25
  stop
  for i=0,imgsize[0]-1 do begin
     for j=0,imgsize[1]-1 do begin
        if r[i,j] le 100. then mask[i,j]=0
        if (ell[i,j] le 100.) and (r[i,j] le 100.) then begin
           print,i,j
           stop
        endif
        
     endfor
  endfor
  
  stop
END
