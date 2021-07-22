 
c***********************************************************************
 
	subroutine phder(fl,tl,phi,hst,nosd,notd)
 
c	routine pour EFF, Auteur J. Christensen Dalsgaard, adaptation:
c	A. Baglin, M. Auvergne, P. Morel
 
	implicit real*8 (a-h,o-z)
 
      logical nosd,notd
      dimension phi(30)
      dimension sig(10),tmn(10),ff(4),gg(4),c(48),hst(10)
      common/phdsms/ s0,sf,sg,sff,sfg,sgg,sfff,sffg,sfgg,sggg,
     .   cfg,tf,tg,tff,tfg,tgg,tfff,tffg,tfgg,tggg
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .   ct,crho,cpe,che,ca03,caa,ckh,car
      common/ln10/ amm
	integer ic,imax,i,icase,ic2,l,mio,kk,l0,k,in,im,ig,ik
 
c  logical equivalence
      equivalence(s0,sig(1)),(cfg,tmn(1))
c     ***********************************
c  coefficients
      data ic,c/4,2.315472,7.128660,7.504998,2.665350,7.837752,23.507934
     .,23.311317,7.987465,9.215560,26.834068,25.082745,8.020509,3.693280
     .,10.333176,9.168960,2.668248,2.315472,6.748104,6.564912,2.132280,
     .7.837752,21.439740,19.080088,5.478100,9.215560,23.551504,19.015888
     .,4.679944,3.693280,8.859868,6.500712,1.334124,1.157736,3.770676,
     .4.015224,1.402284,8.283420,26.184486,28.211372,10.310306,14.755480
     .,45.031658,46.909420,16.633242,7.386560,22.159680,22.438048,
     .7.664928/
c
c  number of sums
      imax=10
      if(notd) imax=6
      if(nosd) imax=3
c  f,g
      f=1.e1**fl
      t=1.e1**tl
      ts=ct*t
      wf= sqrt(1+f)
      g=ts*wf
c
c  1/(1+f), 1/(1+g) etc.
c
      vf=1/(1+f)
      vg=1+g
      fdf=vg*g
      fdf=vf*fdf* sqrt(fdf)
      vg=1/vg
      vfg=vf*vg
      ug=g*vg
      uf=f*vf
c
      ug2=ug*vg
      ug3=vg*ug2
c  powers of f and g
      ff(1)=f
      gg(1)=1
      do 10 i=2,ic
      ff(i)=f*ff(i-1)
      gg(i)=g*gg(i-1)
   10 fdf=fdf*vfg
c  test on size of f and g
      icase=1
      if(.not.notd) go to 12
      if(f.lt.1.e-4) icase=icase+1
      if(g.lt.1.e-4) icase=icase+2
c
   12 ic2=ic*ic
c
c  calculate phi* and derivatives
c
      l=1
      anu=1.5e0
      mio=ic
      an32=2.5e0-ic
      kk=-10
c
      l0=1
      do 50 k=1,3
      kk=kk+10
      if(k-2) 18,16,17
c  reset fdf for k=2 and 3
   16 anu=2.5e0
      fdf=g*fdf
      go to 18
   17 mio=mio+1
      fdf=fdf*vf
   18 do 19 i=1,imax
   19 sig(i)=0.e0
      annu=anu-1
c
c  the summation
c
      l=l0
      go to (20,25,30,35), icase
c  the general case
   20 do 23 in=1,ic
      annu=annu+1
      do 23 im=1,ic
c  phimn*(f**(m+1))*(g**n)
      cfg=c(l)*ff(im)*gg(in)
c
      tg=annu*cfg
      tf=im*cfg
      if(nosd) go to 21
c  second derivatives
      tgg=annu*tg
      tfg=im*tg
      tff=im*tf
      if(notd) go to 21
c  third derivatives
      tggg=annu*tgg
      tfgg=im*tgg
      tffg=im*tfg
      tfff=im*tff
c  summing
   21 do 22 i=1,imax
   22 sig(i)=sig(i)+tmn(i)
   23 l=l+1
c  the summation is finished
      if(nosd) go to 45
c
c  the sigma tilde (cf (22.2)) are stored in the corresponding sigma.
c  this is o.k. provided that we go backwards.
c
      s02=s0*s0
      sg2=sg*sg
      sf2=sf*sf
      if(notd) go to 24
c  third derivatives
      s03=s02*s0
      sfff=(sfff*s02-sf*(3*sff*s0-2*sf2))/s03
      sffg=(sffg*s02-sff*sg*s0-2*sf*(sfg*s0-sg*sf))/s03
      sfgg=(sfgg*s02-sgg*sf*s0-2*sg*(sfg*s0-sg*sf))/s03
      sggg=(sggg*s02-sg*(3*sgg*s0-2*sg2))/s03
c  second derivatives
   24 sff=(sff*s0-sf2)/s02
      sfg=(sfg*s0-sf*sg)/s02
      sgg=(sgg*s0-sg2)/s02
      go to 45
c  f is small
   25 do 28 in=1,ic
      annu=annu+1
      do 27 im=1,2
      cfg=c(l)*ff(im)*gg(in)
      sig(im)=sig(im)+cfg
      tg=annu*cfg
      sg=sg+tg
      if(nosd) go to 27
      sgg=sgg+annu*tg
      sfg=sfg+im*tg
   27 l=l+1
   28 l=l+2
c  the summation is finished. set precursors for sigma tilde
      sff=s0*sf
      s0=s0+sf
      sf=s0+sf
      if(nosd) go to 45
      sgg=sgg*s0-sg*sg
      go to 40
c  g is small
   30 ig=1
      do 33 in=1,2
      annu=annu+1
      do 32 im=1,4
      cfg=c(l)*ff(im)*gg(in)
      sig(ig)=sig(ig)+cfg
      tf=im*cfg
      sf=sf+tf
      if(nosd) go to 32
      sff=sff+im*tf
      sfg=sfg+annu*tf
   32 l=l+1
   33 ig=3
c  the summation is finished. set precursors for sigma tilde.
      sgg=s0*sg
      s0=s0+sg
      sg=anu*s0+sg
      if(nosd) go to 45
      sff=sff*s0-sf*sf
      go to 40
c  both f ang g are small
   35 ig=3
c  in this case we must also zero sfg
      sfg=0.e0
      do 38 in=1,2
      annu=annu+1
      do 37 im=1,2
      cfg=c(l)*ff(im)*gg(in)
      sig(im)=sig(im)+cfg
      sig(ig)=sig(ig)+cfg
   37 l=l+1
      ig=5
   38 l=l+2
c  the summation is finished. set precursors for sigma tilde.
      sff=s0*sf
      s0=s0+sf
      sf=s0+sf
      sgg=sg*sfg
      sg=anu*s0+sfg
      if(nosd) go to 45
c  set final values of the sigma tilde.
   40 s02=s0*s0
      sff=sff/s02
      sgg=sgg/s02
      if(f*g.lt.1.00001e-8) go to 42
      sfg=(sfg*s0-sf*sg)/s02
      go to 45
c  approximate expression for sfg (may need fixing up, if f = o(1)
c  or g = o(1))
   42 sfg=f*g*(c(l0+5)-c(l0+1)*c(l0+4)/c(l0))/c(l0)
c
c  phi* and first derivatives
c
   45 phi(kk+1)=fdf*s0
      pht=an32*ug+sg/s0
      phi(kk+3)=pht
      phi(kk+2)=(pht/2-mio)*uf+sf/s0
      if(nosd) go to 50
c
c  second derivatives of phi*.
c
      phtt=an32*ug2+sgg
      phi(kk+6)=phtt
      phi(kk+5)=phtt*uf/2+sfg
      phi(kk+4)=sff+uf*(sfg+vf*(pht/2-mio+f*phtt/4))
c
      if(notd) go to 50
c  third derivatives
      phttt=an32*ug3*(1-g)+sggg
      phi(kk+10)=phttt
      phi(kk+9)=sfgg+uf*phttt/2
      phfft=sffg+uf*(sfgg+vf*(phtt+f*phttt/2)/2)
      phi(kk+8)=phfft
      phi(kk+7)=sfff+uf*(sffg+phfft/2+vf*(1.5*sfg+f*sfgg/4
     .   +vf*((1-f)*(pht/2-mio)+f*phtt/2)))
   50 l0=l0+ic2
c
c  h* and its derivatives (pp 23-25)
c
      do 55 i=2,imax
      ik=20+i
   55 phi(ik)=phi(ik)-phi(i)
c
      hs=phi(21)/phi(1)
      wft1=2*g
      hst(1)=hs+wft1
c
      hf=phi(22)
      ht=phi(23)
      wft2=ts*f/wf
      hst(2)=hs*hf+wft2
      hst(3)=hs*ht+wft1
      if(nosd) go to 58
c  second derivatives
      hff=phi(24)
      hft=phi(25)
      htt=phi(26)
      wft3=uf*(1+f/2)*ts/wf
      hst(4)=hs*(hf*hf+hff)+wft3
      hst(5)=hs*(hf*ht+hft)+wft2
      hst(6)=hs*(ht*ht+htt)+wft1
      if(notd) go to 58
c  third derivatives
      hst(7)=hs*(hf*(hf*hf+3*hff)+phi(27))+uf*vf*(1+f*(2+f)/4)*ts/wf
      hst(8)=hs*(hf*(hf*ht+2*hft)+ht*hff+phi(28))+wft3
      hst(9)=hs*(ht*(ht*hf+2*hft)+hf*htt+phi(29))+wft2
      hst(10)=hs*(ht*(ht*ht+3*htt)+phi(30))+wft1
c  change to derivatives wrt log10 f and log10 t
   58 fct=amm
      do 60 i=2,imax
      if(i.eq.4.or.i.eq.7) fct=fct*amm
   60 hst(i)=hst(i)*fct
      return
 
	end
