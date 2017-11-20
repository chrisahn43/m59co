@sersic2mge.pro
pro m59co_sersic2mge
                                ;Start with F475W best-fits output from Galfit
                                ;Free parameters fit
  scale=0.05
  Msun=5.11d
  A_B=0.107d
  zeropt=26.05923
  re1=3.13d*scale
  re2=12.85d*scale
  sersic1={mag:19.4d,re:re1,n:1.06d,pa:-65.24d,q:0.97d}
  sersic2={mag:18.38d,re:re2,n:1.09d,pa:88.42d,q:0.98d}
  sersics=[sersic1,sersic2]
  mge=sersics2mge(sersics,Msun,A_B=A_B)
  forprint,mge[*,0],mge[*,1],mge[*,2],mge[*,3],textout='m59co_mge_outputsersic_free.dat',format='F,F,F,F'
                                ;Fixed paramters fit
  re1=2.96*scale
  re2=12.25*scale
  sersic1={mag:19.56d,re:re1,n:1.02d,pa:34.06d,q:0.99d}
  sersic2={mag:18.32d,re:re2,n:1.21d,pa:17.69d,q:0.98d}
  sersics=[sersic1,sersic2]
  mge=sersics2mge(sersics,Msun,A_B=A_B)
  forprint,mge[*,0],mge[*,1],mge[*,2],mge[*,3],textout='m59co_mge_outputsersic.dat',format='F,F,F,F'
                                ;Now do F850LP best-fit output from Galfit
                                ;Free parameters fit
  Msun=4.52d
  A_B=0.041d
  zeropt=24.84245
  re1=2.96d*scale
  re2=12.25d*scale
  sersic1={mag:18.17d,re:re1,n:1.02d,pa:34.06d,q:0.99d}
  sersic2={mag:16.72d,re:re2,n:1.21d,pa:17.69d,q:0.98d}
  sersics=[sersic1,sersic2]
  mge=sersics2mge(sersics,Msun,A_B=A_B)
  forprint,mge[*,0],mge[*,1],mge[*,2],mge[*,3],textout='m59co_mge_outputsersic850_free.dat',format='F,F,F,F'
                                ;Fixed parameters fit
  re1=3.13d*scale
  re2=12.85d*scale
  sersic1={mag:17.98d,re:re1,n:1.06d,pa:-65.24d,q:0.97d}
  sersic2={mag:16.78d,re:re2,n:1.09d,pa:88.42d,q:0.98d}
  sersics=[sersic1,sersic2]
  mge=sersics2mge(sersics,Msun,A_B=A_B)
  forprint,mge[*,0],mge[*,1],mge[*,2],mge[*,3],textout='m59co_mge_outputsersic850.dat',format='F,F,F,F'
  
  stop
END
