c*******************************************************************
 
	subroutine hviona(fl,tl,x,y,z,nosd,notd,anu,ue,anur,uer)
 
c	routine pour EFF, Auteur J. Christensen Dalsgaard, adaptation:
c	A. Baglin, M. Auvergne, P. Morel
 
	implicit real*8 (a-h,o-z)
 
c      implicit real*16 (a-h,o-z)
      logical nosd,notd
      dimension anu(1),ue(1),anur(1),uer(1)
      dimension eir(10),dr(10),hr(10),gr(10),dmup(10),chir(27),ir(6),
     .   izr(6),abr(6),amr(6),iz(23),name(23),chi(374),xi(29),phi(30),
     .   hst(30),am(23)
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .   ct,crho,cpe,che,ca03,caa,ckh,car,ergev
      common/potetc/ chi,am,iz,name
      common/hvabnd/ ab(10),iab
      common/eqscnt/ anz0,ihvz
      common/dmuder/ dmup,idmu
      common/ln10/ amm,amm2

c D. Cordier, 5 octobre 2011 : on supprime "name" de la ligne suivante, il
c y a en effet un incohérence, c'est un caractère dans le common "potetc"
c	integer*4 iz,name,iab,ihvz,idmu,icount,iwrite,lj,k,izj,i,im,irab,
	integer*4 iz,iab,ihvz,idmu,icount,iwrite,lj,k,izj,i,im,irab,	
     1	ir,jfirst,j,ljr,izi,izrj,izr,idmx,izoj,lj0,imax,
     2	ii,l
c
      data icount/0/,iwrite/3/
c
      if(icount.gt.0) go to 10
      icount=1
      if(iwrite.ne.2) go to 61
      write(6,200)
      lj=0
      do 3 k=1,10
      izj=iz(k)
      do 2 i=1,izj
      lj=lj+1
    2 xi(i)=chi(lj)
      im=min(izj,13)
    3 write(6,210) name(k),izj,(xi(i),i=1,im)
  200 format(///' ionization potentials:'/)
  210 format(1x,a4,i4,13f9.2)
c  reset to restricted set of variables
   61 irab=3
      ir(1)=1
      ir(2)=3
      ir(3)=10
c  element with lowest ionization potential (here fe)
      jfirst=3
c
      j=1
      lj=0
      ljr=0
      do 5 i=1,10
      izi=iz(i)
      if(i.ne.ir(j)) go to 5
      izrj=izi
      if(j.gt.2) izrj=1
      izr(j)=izrj
      amr(j)=am(i)
      do 4 k=1,izrj
      lj=lj+1
      ljr=ljr+1
    4 chir(ljr)=chi(lj)
      izi=izi-izrj
      j=j+1
    5 lj=lj+izi
      abr(1)=1.02*ab(1)
      abr(2)=1.35*ab(3)
      abr(3)=3.71*ab(10)
c  reset anz0
      anz0=0
      do 63 j=1,irab
   63 anz0=anz0+abr(j)*izr(j)/amr(j)
c
      if(iwrite.ne.1) go to 8
c
      write(6,220)
      lj=0
      do 7 k=1,irab
      izrj=izr(k)
      do 6 i=1,izrj
      lj=lj+1
    6 xi(i)=chir(lj)
      im=min(izrj,13)
    7 write(6,210) name(ir(k)),izrj,(xi(i),i=1,im)
  220 format(///' ionization potentials after resetting:'//)
      write(6,230) anz0
  230 format(//' now anz0 =',f10.5//)
c  change to summed ionization potentials
    8 lj=0
      do 9 k=1,irab
      sum=0
      izj=izr(k)
      do 9 i=1,izj
      lj=lj+1
      sum=sum+chir(lj)
    9 chir(lj)=sum
c
   10 f=1.e1**fl
      t=1.e1**tl
c  k*t, in ev
      tk=ck1*t
      if(idmu.eq.1) go to 15
c
      call phder(fl,tl,phi,hst,nosd,notd)
c
c
c  psi and derivatives
      wf= sqrt(1+f)
      psi=2*wf+log(f/((1+wf)*(1+wf)))
      psif=amm*wf
      psff=amm2*f/2/wf
c
      zf=x+2*y+6*z
      zf2=zf*zf
      zf3=zf*zf2
c
      ak0=ckh*zf2
c  a0**3*ne/(k*t)
      aa=caa*phi(1)/(tk*zf3)
c
c  delta mu and derivatives
c
      bmu=tk+20*ak0
      dmu=aa*bmu
      dmps=dmu-psi
      dmup(1)=dmps
      ref=phi(2)
      ret=phi(3)
      dmuf=dmu*ref*amm
      dmut=aa*(bmu*ret-20*ak0)*amm
      dmux=aa*(3*tk+20*ak0)/zf
      dmup(2)=dmuf-psif
      dmup(3)=dmut
      dmup(4)=dmux
c
      if(nosd) go to 18
      dmup(5)=(dmuf*ref+dmu*phi(4)*amm)*amm-psff
      dmup(6)=(dmut*ref+dmu*phi(5)*amm)*amm
      dmup(7)=dmux*ref*amm
      dmup(8)=aa*(20*ak0*(1-2*ret)+bmu*(phi(6)+ret*ret))*amm2
      dmup(9)=aa*(ret*(3*tk+20*ak0)-20*ak0)*amm/zf
      dmup(10)=aa*(12*tk+40*ak0)/zf2
c
      go to 18
   15 dmps=dmup(1)
c
   18 idmx=10
      if(nosd) idmx=4
      do 19 i=1,10
      anu(i)=0.e0
      ue(i)=0.e0
      anur(i)=0.e0
   19 uer(i)=0.e0
c
c
      lj=0
      tki=1.e0/tk
      do 50 j=1,irab
      izoj=iz(ir(j))
      izj=izr(j)
c
      lj0=lj
      do 20 i=1,izj
      lj=lj+1
      dxp=i*dmps-chir(lj)*tki
   20 xi(i)=dxp
c  complete or no ionization?
c
      if(izj-1) 21,21,25
c  only one level
   21 phm=xi(1)
      if(phm) 22,22,23
   22 imax=0
      phmh=phm
      phm=0
      go to 30
c
   23 imax=1
      phmh=0
      go to 30
c  more than one level
   25 ii=izj+1
      xi(ii)=0
      imax=1
      phm=xi(1)
c  largest xi
      do 26 i=2,ii
      if(xi(i).le.phm) go to 26
      imax=i
      phm=xi(i)
   26 continue
c  largest but one xi
      phmh=-100
      do 27 i=1,ii
      if(i.eq.imax) go to 27
      phmh=max(phmh,xi(i))
   27 continue
c
      if(imax.eq.ii) imax=0
c  test for complete ionization
   30 dphm=phm-phmh
      dptst=19
      if(imax.lt.izj.or.dphm.le.dptst) go to 32
c  complete ionization
      fct1=abr(j)/amr(j)
      anur(j)=fct1*izj
      uer(j)=fct1*chir(lj)
      anu(1)=anu(1)+anur(j)
      ue(1)=ue(1)+uer(j)
      go to 50
c  test for no ionization
   32 if(imax.gt.0) go to 34
      if(j.eq.jfirst) dptst=160
      if(dphm.gt.dptst) go to 50
c  general case
   34 do 35 i=1,idmx
      dr(i)=0.e0
      hr(i)=0.e0
   35 gr(i)=0.e0
      if(phm.le.dptst) dr(1)=omega(0,izoj)* exp(-phm)
      do 40 i=1,izj
      cchi=chir(lj0+i)
      dxp=xi(i)-phm
      if(dxp.lt.-dptst) go to 40
      dxp=omega(i,izoj)* exp(dxp)
      eir(1)=1
      do 36 k=2,4
   36 eir(k)=i*dmup(k)
      eir(3)=eir(3)+cchi*amm*tki
      ii=4
      if(nosd) go to 38
      do 37 k=2,4
      do 37 l=k,4
      ii=ii+1
   37 eir(ii)=eir(k)*eir(l)+i*dmup(ii)
      eir(8)=eir(8)-amm2*cchi*tki
c
   38 do 39 k=1,idmx
      eeir=dxp*eir(k)
      dr(k)=dr(k)+eeir
      hr(k)=hr(k)+i*eeir
   39 gr(k)=gr(k)+cchi*eeir
   40 continue
c
      dr1i=1.e0/dr(1)
c
      fct1=abr(j)/(amr(j)*dr(1))
      hrdr=hr(1)*dr1i
      grdr=gr(1)*dr1i
      anur(j)=fct1*hr(1)
      uer(j)=fct1*gr(1)
      anu(1)=anu(1)+anur(j)
      ue(1)=ue(1)+uer(j)
c  derivatives
      ii=4
      do 48 k=2,4
      anu(k)=anu(k)+fct1*(hr(k)-dr(k)*hrdr)
      ue(k)=ue(k)+fct1*(gr(k)-dr(k)*grdr)
      if(nosd) go to 48
      do 45 l=k,4
      ii=ii+1
      anu(ii)=anu(ii)+fct1*(hr(ii)-(dr(k)*hr(l)+dr(l)*hr(k)+(dr(ii)
     .   -2*dr(k)*dr(l)*dr1i)*hr(1))*dr1i)
   45 ue(ii)=ue(ii)+fct1*(gr(ii)-(dr(k)*gr(l)+dr(l)*gr(k)+(dr(ii)
     .   -2*dr(k)*dr(l)*dr1i)*gr(1))*dr1i)
   48 continue
   50 continue
c
      return
 
	end
 
c******************************************************************
 
	blockdata bleqst
 
c	routine pour EFF, Auteur J. Christensen Dalsgaard, adaptation:
c	A. Baglin, M. Auvergne, P. Morel
 
	implicit real*8 (a-h,o-z)
 
c      implicit real*16 (a-h,o-z)

c D. Cordier, 5 octobre 2011, "name" doit être une chaîne de caractères
c	integer*4 name,iz
	integer*4 iz
	character*4 name
      common/potetc/ chi(374),am(23),iz(23),name(23)
      data iz/6,7,8,10,11,12,13,14,18,26,13*9999/
      data name/'   c','   n','   o','  ne','  na',
     .   '  mg','  al','  si','   a','  fe',13*'bido'/
      data am/12.00,14.01,16.00,20.17,22.99,24.31,26.98,28.08,
     .   39.94,55.84,13*99999./
      data chi/11.26,24.38,47.86,64.48,391.99,489.84,
     .   14.54,29.60,47.43,77.45,97.86,551.92,666.83,
     .   13.61,35.15,54.93,77.39,113.87,138.08,739.11,871.12,
     .   21.56,41.07,63.5,97.16,126.4,157.91,207.3,239.,1196.,1360.,
     .   5.14,47.29,71.65,98.88,138.60,172.36,208.44,264.15,299.78,
     .   1465.,1646.,
     .   7.64,15.03,80.12,109.29,141.23,186.86,225.31,265.96,327.90,
     .   367.36,1761.2,2085.,
     .   5.98,18.82,28.44,119.96,153.77,190.42,241.93,285.13,330.1,
     .   398.5,441.9,2085.5,2299.,
     .   8.15,16.34,33.46,45.13,166.73,205.11,246.41,303.87,
     .   351.83,401.3,476.0,523.2,2436.,2666.,
     .   15.75,27.62,40.90,59.79,75.0,91.3,124.,143.46,421.,480.,
     .   539.5,621.1,688.5,755.5,854.4,3*1.d20,
     .   7.90,16.18,30.64,56.,79.,105.,133.,151.,235.,262.,290.,321.,
     .   355.,390.,457.,11*1.d20,249*88888./
 
	end
