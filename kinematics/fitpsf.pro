PRO FITPSF
                                ;THIS CODE IS USED TO FIND THE
                                ;KINEMATIC PSF FOR A DOUBLE GAUSSIAN
                                ;AND A GAUSSIAN+ MOFFAT FUNCTION
                                ;PROVIDED BY ANIL. SECOND TIME USING
                                ;KBAND IMAGE

;  Gauss+gauss
;Iter     11   CHI-SQUARE =       1382.7916          DOF = 955
;    P(0) =              4.32970 inner gaussian width=0.216"
;    P(1) =              42.8391
;    P(2) =             0.398045
;    P(3) =             0.275526
;    P(4) =              18.5852 outer gaussian width=0.929"
;    P(5) =             0.904470  (f_outer=0.90/(1+0.90))
;Gauss+Moffat
;Iter      8   CHI-SQUARE =       1544.9711          DOF = 956
;    P(0) =              2.74309 inner gaussian width=0.137"
;    P(1) =              42.4860
;    P(2) =             0.393562
;    P(3) =             0.280785
;    P(4) =              21.6700 outer moffat width=1.05"
;    P(5) =              3.16352 f_moffat=3.16/(1+3.16)
;    P(6) =              4.76500

;GAUSS + MOFFAT PSF
  scaling=0.05
  p=[2.74309,42.4860,0.393562,0.280785,21.6700,3.16352,4.76500]
  r=range(0.01,2.5,100,/log)
  psfwidth1=0.216/2.35
  psfwidth2=p[4]/2.35*scaling
  psf1=exp(-(r/psfwidth1)^2/2.)/(2.*!PI*psfwidth1^2)
  psf1=psf1/TOTAL(psf1)
  psf2=1./(1+r/psfwidth2)^p[6]
  psf2=p[5]*psf2/total(psf2)
  psf=psf1+psf2
  psf=psf/total(psf)
  psf=p[1]*psf

  mge_fit_1d,r,psf2,ngauss=10,sol=sol
  forprint,sol[0,*],sol[1,*],format='F,F',textout='temp.dat'
  readcol,'temp.dat',modelweight,modelsig,format='F,F'
  spawn, 'rm -rf temp.dat'
  

  modelflux=fltarr(n_elements(r))
  modelpeak=((total(psf2)/(1+total(psf2))))*(modelweight)/(total(modelweight))
  modelpeak=[1-total(modelpeak),modelpeak]
  modelsig=[psfwidth1,modelsig]
  forprint,modelpeak,modelsig,format='F,F',textout='kinematic_psf_moffat.dat'
  modelsig=sol[1,*]
  for i=0,n_elements(r)-1 do begin
     temp=0.
     for j=0,n_elements(modelweight)-1 do begin
        temp+= (modelweight[j]*exp(-(1./(2*modelsig[j]^2))*(r[i])^2))
        modelflux[i]=temp
     endfor
  endfor
;  stop
;  djs_plot,r,psf2,psym=2,/xlog
;  djs_oplot,r,modelflux,color='blue'
;  print,mean(((psf2-modelflux)/psf2)*100)
  stop
;DOUBLE GAUSSIAN
  p=[4.32970,42.8391,0.398045,0.275526,18.5852,0.904470]
  psfwidth1=0.216/2.35
  psfwidth2=0.929/2.35
  psf1=exp(-(r/psfwidth1)^2/2.)/(2.*!PI*psfwidth1^2)
  psf1=psf1/TOTAL(psf1)
  psf2=exp(-(r/psfwidth2)^2/2.)/(2.*!PI*psfwidth2^2)
  psf2=p[5]*psf2/total(psf2)
  psf=psf1+psf2
  psf=psf/total(psf)
  lightfrac=p[5]/(1+p[5])
  modelpeak=[(1-lightfrac),lightfrac]
  modelsig=[psfwidth1,psfwidth2]
  forprint,modelpeak,modelsig,format='F,F',textout='kinematic_psf_gauss.dat'

  stop

  p=[4.27511,0.0122303,0.597095,0.2125112,19.1576,0.850073]
  psfwidth1=(p[0]*0.05)/2.35
  psfwidth2=(p[4]*0.05)/2.35
  lightfrac=p[5]/(1+p[5])
  modelpeak=[(1-lightfrac),lightfrac]
  modelsig=[psfwidth1,psfwidth2]
  forprint,modelpeak,modelsig,format='F,F',textout='kinematic_oldpsf_gauss.dat'
END
