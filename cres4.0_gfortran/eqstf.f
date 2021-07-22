 
c**********************************************************************
 
	subroutine eqstf(fl,tl,x,y,z,nosd,notd)
 
c	routine pour EFF, Auteur J. Christensen Dalsgaard, adaptation:
c	A. Baglin, M. Auvergne, P. Morel
 
	implicit real*8 (a-h,o-z)
 
c      implicit real*16 (a-h,o-z)
 
	logical tstl,nosd,notd,cmplio
	integer ihvz,idmu,l,n,ider,iders,k,ia,i,ii,l1,m,m1,jj,j,k1,k2,
     1	kd,imx,i10,lm,nu,nb,kl,nuu,ln,lt,lf,mn
 
      dimension phi(30),hst(10),ex(3),ea(30),xi(30),dxt(4),dxt1(4),
     .  dxt2(4),ane(10),rho(20),dne(10),dph(20),ddmu(11),he(20),pe(10),
     .   hi(20),pi(10),hh(20),ph(10),hr(20),pr(10),ht(20),pt(20),cp(4),
     .   dad(4),dlt(4),xii(30),anuh(10),ueh(10),anuhr(23),uehr(23)
      common/eqscnt/ anh0,ihvz
      common/ln10/ amm,amm2,amm3
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .   ct,crho,cpe,che,ca03,caa,ckh,car,ergev
      common/dmuder/ dmu,dmuf,dmut,dmux,dmuff,dmuft,dmufx,dmutt,dmutx,
     .   dmuxx,idmu
      common/eqsout/ ea,xii,dne,dph,he,pe,hi,pi,hh,ph,hr,pr
      common/eqstd/ xii1(4),ane,rho,ht,pt,cp,dad,dlt,gm1,tprh,trhp,
     .   rhxp
      common/eqsaux/psi
c
      equivalence(ex(1),exh),(ddmu(1),dmu)
c
      fxp(a)= exp(min(8.d1,max(a,-1.d2)))
      nuu(l,n)=((l-1)*(8-l))/2+n-l+5
c
c  number of derivatives
      ider=20
      if(notd) ider=10
      if(nosd) ider=4
      iders=min(ider,10)
c
      f=1.e1**fl
      t=1.e1**tl
c
      call phder(fl,tl,phi,hst,nosd,notd)
c
c
c  psi and derivatives
      wf= sqrt(1+f)
      psi=2*wf+log(f/((1+wf)*(1+wf)))
      psif=amm*wf
      psff=amm2*f/(2*wf)
c
      zf=x+2*y+6*z
      zf2=zf*zf
      zf3=zf*zf2
c
      ak0=ckh*zf2
c  k*t, in ev
      tk=ck1*t
c  a0**3*ne/(k*t)
      aa=caa*phi(1)/(tk*zf3)
c
c  delta mu and derivatives
c
      bmu=tk+20*ak0
      dmu=aa*bmu
c
      ref=phi(2)
      ret=phi(3)
      dmuf=dmu*ref*amm
      dmut=aa*(bmu*ret-20*ak0)*amm
      dmux=aa*(3*tk+20*ak0)/zf
c
      if(nosd) go to 10
      dmuff=(dmuf*ref+dmu*phi(4)*amm)*amm-psff
      dmuft=(dmut*ref+dmu*phi(5)*amm)*amm
      dmufx=dmux*ref*amm
      dmutt=aa*(20*ak0*(1-2*ret)+bmu*(phi(6)+ret*ret))*amm2
      dmutx=aa*(ret*(3*tk+20*ak0)-20*ak0)*amm/zf
      dmuxx=aa*(12*tk+40*ak0)/zf2
   10 dmu=dmu-psi
      dmuf=dmuf-psif
      idmu=1
c  test for complete ionization of h and he
      cmplio=(dmu-54.4/tk).gt.19
      if(cmplio) go to 31
c
c  e h, e he, e he+ and derivatives
c
      k=-10
      do 25 ia=1,3
      k=k+10
      ext=ex(ia)/tk
      eea=dmu-ext
      if(eea+30) 15,15,21
c  no ionization
   15 do 16 i=1,10
   16 ea(k+i)=0.e0
      go to 25
c
   21 eea=fxp(eea)
      if(ia-2)  22,23,22
   22 eea=eea/2
      go to 24
   23 eea=2*eea
   24 ea(k+1)=eea
c  first derivatives
      ea(k+2)=dmuf*eea
      ea(k+3)=(amm*ext+dmut)*eea
      ea(k+4)=dmux*eea
      if(nosd) go to 25
c  second derivatives
      ea(k+5)=dmuff*eea+dmuf*ea(k+2)
      ea(k+6)=dmuft*eea+dmuf*ea(k+3)
      ea(k+7)=dmufx*eea+dmux*ea(k+2)
      ea(k+8)=(dmutt-amm2*ext)*eea+(amm*ext+dmut)*ea(k+3)
      ea(k+9)=dmutx*eea+dmux*ea(k+3)
      ea(k+10)=dmuxx*eea+dmux*ea(k+4)
   25 continue
c  reset e he+
      eea=ea(21)
      if(eea.eq.0) go to 35
      eea1=ea(11)
      ea(21)=eea*eea1
      ii=24
      if(nosd) go to 27
      do 26 l=1,3
      l1=l+21
      do 26 m=l,3
      ii=ii+1
      m1=m+21
   26 ea(ii)=eea*ea(ii-10)+ea(l1)*ea(m1-10)+ea(l1-10)*ea(m1)+
     .   eea1*ea(ii)
   27 do 28 i=22,24
   28 ea(i)=ea(i)*eea1+eea*ea(i-10)
   30 continue
      go to 35
c
c  x h+, x he+, x he++ and derivatives
c
c  complete ionization
   31 xi(1) =1
      xi(11)=0
      xi(21)=1
      do 33 i=2,22,10
      jj=i+8
      do 33 j=i,jj
   33 xi(j)=0
      go to 50
c  partial ionization
   35 dnm=1+ea(1)
      xi(1)=ea(1)/dnm
      ii=4
      dnm2=dnm*dnm
      dnm3=dnm*dnm2
      do 40 l=1,3
      l1=l+1
      eal=ea(l1)
      xi(l1)=eal/dnm2
      if(nosd) go to 40
      do 38 m=l,3
      m1=m+1
      ii=ii+1
   38 xi(ii)=ea(ii)/dnm2-2*eal*ea(m1)/dnm3
   40 continue
c
      dnm=1+ea(11)+ea(21)
      dnm2=dnm*dnm
      dnm3=dnm*dnm2
      k1=11
      k2=21
   45 kd=k1-k2
      ii=k1+3
      ea1=ea(k1)
      ea2=ea(k2)
      xi(k1)=ea1/dnm
      do 48 l=1,3
      l1=k1+l
      eal1=ea(l1)
      eal2=ea(k2+l)
      anm=(1+ea2)*eal1-ea1*eal2
      xi(l1)=anm/dnm2
      if(nosd) go to 48
      do 46 m=l,3
      ii=ii+1
      eam1=ea(k1+m)
      eam2=ea(k2+m)
   46 xi(ii)=(eal1*eam2+(1+ea2)*ea(ii)-eal2*eam1-ea1*ea(ii-kd))/dnm2
     .   -2*anm*(eam1+eam2)/dnm3
   48 continue
      if(k1.eq.21) go to 50
      k1=21
      k2=11
      go to 45
c  ionization of heavy elements
   50 if(ihvz) 51,51,52
   51 anuh(1)=anh0
      call zero(anuh(2),9)
      call zero(ueh,10)
      go to 53
   52 call hviona(fl,tl,x,y,z,nosd,notd,anuh,ueh,anuhr,uehr)
c
c  ne and derivatives
c
c  combine he fractions for an in xi(20+.), and for ionization
c  energy in xi(10+.)
   53 exhc=exhe+exhep
c
      imx=10
      if(nosd) imx=4
      do 55 i=1,21,10
   55 call store(xi(i),xii(i),imx)
      xia=xi(1)
      imx=imx+20
c
      do 56 i=21,imx
      i10=i-10
      xio=exhe*xi(i10)+exhc*xi(i)
      xi(i)=2*xi(i)+xi(i10)
   56 xi(i10)=xio
c  terms needed for x-derivatives
      do 57 l=1,4
   57 dxt(l)=av*(xi(l)/ah-xi(20+l)/ahe)
c
      ane(1)=av*(xi(1)*x/ah+xi(21)*y/ahe+z*anuh(1))
c
      ii=4
      do 60 l=1,3
      l1=l+1
      anel=(xi(l1)*x/ah+xi(l+21)*y/ahe+z*anuh(l1))*av
      tstl=l.eq.3
      if(tstl) anel=anel+dxt(1)
      ane(l1)=anel
      if(nosd) go to 60
      do 58 m=l,3
      m1=m+1
      ii=ii+1
      anelm=(xi(ii)*x/ah+xi(20+ii)*y/ahe+z*anuh(ii))*av
      if(tstl) anelm=anelm+dxt(m1)
      if(m.eq.3) anelm=anelm+dxt(l1)
   58 ane(ii)=anelm
   60 continue
c
c  the density and derivatives (note that rho(2) = dlog rho/dlog f,
c  and so on).
c
      anee=ane(1)
      rho(1)=crho*phi(1)/anee
      ii=4
      jj=10
      do 63 l=1,3
      l1=l+1
      rhol=-ane(l1)/amm/anee
      tstl=l.le.2
      if(tstl) rhol=rhol+phi(l1)
      rho(l1)=rhol
      if(nosd) go to 63
      do 62 m=l,3
      ii=ii+1
      m1=m+1
      lm=l+m
      rholm=(ane(l1)*ane(m1)/anee-ane(ii))/anee/amm
      if(tstl.and.m.le.2) rholm=rholm+amm*phi(2+lm)
      if(notd) go to 62
      do 61 n=m,3
      jj=jj+1
      rhd=-2*ane(l1)*ane(m1)*ane(n+1)/amm/anee**3
      if(l.lt.3.and.m.lt.3.and.n.lt.3) rhd=rhd+amm2*phi(4+lm+n)
   61 rho(jj)=rhd
   62 rho(ii)=rholm
   63 continue
      call zero(pt,20)
      call zero(ht,20)
c
c
c  delta p, delta h and derivatives.
c
c  test for complete ionization
      if(cmplio) go to 80
c
c  delta ne**2
c
c  reset xi(1) and xi(21) (note that xi(2-10) and xi(22-30) are still
c  consistent).
      xi(1)=-1/(ea(1)+1)
      xi(21)=-(2+ea(11))/(1+ea(11)+ea(21))
c
      do 65 l=1,4
      dxt1(l)=xi(l+20)/ahe-xi(l)/ah
      dxtl=-xi(l+20)*y/ahe-xi(l)*x/ah-anuh(l)*z
      if(l.eq.1) dxtl=dxtl+anh0*z
      dxtl2=ane(l)
      if(l.le.3) go to 64
      dxtl=dxtl+dxt1(1)
      dxtl2=dxtl2+avda
   64 dxt(l)=dxtl*av
   65 dxt2(l)=dxtl2
c
      ann=ane(1)+(x/ah+2*y/ahe+z*anh0)*av
      xtt=dxt(1)
      dne(1)=ann*xtt
c
      ii=4
      do 70 l=1,3
      l1=l+1
      dne(l1)=dxt2(l1)*xtt+ann*dxt(l1)
      if(nosd) go to 70
      tstl=l.eq.3
      do 68 m=l,3
      m1=m+1
      ii=ii+1
      dnlm=-xi(20+ii)*y/ahe-xi(ii)*x/ah
      if(tstl) dnlm=dnlm+dxt1(m1)
      if(m.eq.3) dnlm=dnlm+dxt1(l1)
   68 dne(ii)=ane(ii)*xtt+dxt2(l1)*dxt(m1)+dxt2(m1)*dxt(l1)
     .   +ann*dnlm*av
   70 continue
c  quantities common to delta p and delta h (dph(15-20) is used as
c  intermediate storage).
      a03=ca03/zf3/2
      call zero(dph(15),6)
      dxt1(1)=0
      c1=amm*tk
      dxt1(2)=c1
      dph(18)=amm2*tk
      dph(19)=3*c1/zf
      dnee=dne(1)
      a03=a03*rho(1)*rho(1)
      nu=2
      nb=20
      k1=1
   75 do 77 l=2,4
   77 dxt(l-1)=dne(l)+nu*amm*dnee*rho(l)
c
      dxt1(3)=(3*tk+nb*ak0)/zf
      bmu=tk+nb*ak0
      dph(k1)=a03*bmu*dnee
      ii=k1+3
      jj=4
      do 79 l=1,3
      l1=l+1
      dph(k1+l)=a03*(bmu*dxt(l)+dnee*dxt1(l))
      if(nosd) go to 79
      do 78 m=l,3
      ii=ii+1
      jj=jj+1
      m1=m+1
      dphlm=bmu*(dne(jj)+nu*amm*(dne(l1)*rho(m1)+dne(m1)*rho(l1)
     .   +dnee*(rho(jj)+nu*amm*rho(m1)*rho(l1))))+dxt(l)*dxt1(m)
     .   +dxt(m)*dxt1(l)+dnee*dph(10+jj)
      if(m.eq.3.and.l.eq.3) dphlm=dphlm+dnee*(12*tk+2*nb*ak0)/zf2
   78 dph(ii)=a03*dphlm
   79 continue
      if(nu.eq.1) go to 90
      nu=1
      nb=40
      a03=a03/rho(1)
      k1=11
      go to 75
c  complete ionization
   80 call zero(dne,10)
      call zero(dph,20)
      go to 100
   90 continue
      do 95 i=1,iders
      pt(i)=pt(i)+dph(i)
   95 ht(i)=ht(i)+dph(i+10)
c  electron pressure and enthalpy
  100 ii=4
      jj=10
      pee=cpe*phi(11)
      pt(1)=pt(1)+pee
      hsst=hst(1)
      hee=che*anee*hsst
      ht(1)=ht(1)+hee
      pe(1)=pee
      he(1)=hee
      do 110 l=1,3
      l1=l+1
      hll=0
      hel=hsst*ane(l1)
      pel=0
      tstl=l.eq.3
      if(tstl) go to 102
      pel=amm*pee*phi(11+l)
      hll=hst(l1)
      hel=hel+anee*hll
      pt(l1)=pt(l1)+pel
  102 pe(l1)=pel
      hel=che*hel
      ht(l1)=ht(l1)+hel
      he(l1)=hel
      if(nosd) go to 110
      do 108 m=l,3
      ii=ii+1
      m1=m+1
      pelm=0
      helm=ane(ii)*hsst+hll*ane(m1)
      tstl=tstl.or.m.eq.3
      if(tstl) go to 104
      lm=l+m
      pelm=amm2*pee*(phi(11+m)*phi(11+l)+phi(12+lm))
      pt(ii)=pt(ii)+pelm
      helm=helm+ane(l1)*hst(m1)+anee*hst(2+lm)
  104 pe(ii)=pelm
      helm=che*helm
      ht(ii)=ht(ii)+helm
      he(ii)=helm
      if(notd) go to 108
      do 106 n=m,3
      jj=jj+1
      helm=0
      if(tstl) go to 106
      helm=ane(n+1)*hst(2+lm)
      if(n.eq.3) go to 105
      helm=helm+ane(l1)*hst(2+m+n)+ane(m1)*hst(2+l+n)
     .   +anee*hst(4+lm+n)
      pt(jj)=pt(jj)+amm3*pee*(phi(11+l)*phi(11+m)*phi(11+n)
     .   +phi(12+lm)*phi(11+n)+phi(12+l+n)*phi(11+m)
     .   +phi(12+m+n)*phi(11+l)+phi(14+l+m+n))
  105 helm=che*helm
      ht(jj)=ht(jj)+helm
  106 he(jj)=helm
  108 continue
  110 continue
c  ionization enthalpy
c  (pi is introduced for consistency)
      call zero(pi,10)
      call zero(hi(11),10)
c
      xi(1)=xia
      averg=av*ergev
      do 112 l=1,4
  112 dxt(l)=exh*xi(l)/ah-xi(10+l)/ahe
c
      hi1=averg*(exh*xia*x/ah+xi(11)*y/ahe+z*ueh(1))
      hi(1)=hi1
      ht(1)=ht(1)+hi1
      ii=4
      do 115 l=2,4
      dhi=exh*xi(l)*x/ah+xi(l+10)*y/ahe+z*ueh(l)
      tstl=l.eq.4
      if(tstl) dhi=dhi+dxt(1)
      dhi=averg*dhi
      ht(l)=ht(l)+dhi
      hi(l)=dhi
      if(nosd) go to 115
      do 114 m=l,4
      ii=ii+1
      dhi=exh*xi(ii)*x/ah+xi(10+ii)*y/ahe+z*ueh(ii)
      if(tstl) dhi=dhi+dxt(m)
      if(m.eq.4) dhi=dhi+dxt(l)
      dhi=averg*dhi
      ht(ii)=ht(ii)+dhi
  114 hi(ii)=dhi
  115 continue
c  pressure and enthalpy of heavy particles
  120 anh=(x/ah+y/ahe+z/az)*av
c  k*t, in ergs
      tk=ck2*t
      rhh=rho(1)
      phh=anh*rhh*tk
      ph(1)=phh
      ph(2)=amm*phh*rho(2)
      drht=1+rho(3)
      ph(3)=amm*phh*drht
      ph(4)=amm*phh*rho(4)+rhh*tk*avd1
      do 121 i=1,4
  121 pt(i)=pt(i)+ph(i)
      if(nosd) go to 125
      ph(10)=amm*(phh*rho(10)+rho(4)*(ph(4)+rhh*tk*avd1))
      do 122 k=1,3
      k1=k+1
      ph(k+4)=amm*(ph(k1)*rho(2)+phh*rho(k+4))
      if(k.gt.1) ph(k+6)=amm*(ph(k1)*drht+phh*rho(k+6))
  122 continue
      do 123 i=5,10
  123 pt(i)=pt(i)+ph(i)
c
      if(notd) go to 125
      do 124 k=1,3
      k1=k+1
      do 124 l=k,3
      kl=nuu(k,l)
      pt(6+kl)=pt(6+kl)+amm*(ph(kl)*rho(2)+ph(k1)*rho(l+4)+
     .   ph(l+1)*rho(k+4)+phh*rho(kl+6))
      if(k.gt.1) pt(9+kl)=pt(9+kl)+amm*(ph(kl)*drht+ph(k1)*
     .   rho(6+l)+ph(l+1)*rho(6+k)+phh*rho(9+kl))
  124 continue
      pt(20)=pt(20)+amm*(ph(10)*rho(4)+2*ph(4)*rho(10)
     .   +phh*rho(20)+rhh*tk*avd1*(amm*rho(4)*rho(4)+rho(10)))
c
  125 hhh=2.5*anh*tk
      hh(1)=hhh
      hh(2)=0
      hh(3)=amm*hhh
      hh(4)=2.5*tk*avd1
      hh(5)=0
      hh(6)=0
      hh(7)=0
      hh(8)=amm2*hhh
      hh(9)=amm*hh(4)
      hh(10)=0
      call zero(hh(11),10)
      hh(17)=amm*hh(8)
      hh(18)=amm*hh(9)
      ht(1)=ht(1)+hhh
      ht(3)=ht(3)+amm*hhh
      hh4=2.5*tk*avd1
      ht(4)=ht(4)+hh4
      if(nosd) go to 130
      ht(8)=ht(8)+amm2*hhh
      ht(9)=ht(9)+amm*hh4
      if(notd) go to 130
      ht(17)=ht(17)+amm3*hhh
      ht(18)=ht(18)+amm2*hh4
c  pressure and enthalpy of radiation
  130 t4=t*t*t*t
      prr=car*t4/3
      pr(1)=prr
      call zero(pr(2),9)
      pr(3)=4*amm*pr(1)
      pr(8)=4*amm*pr(3)
      pt(1)=pt(1)+prr
      pt(3)=pt(3)+4*amm*prr
      if(nosd) go to 135
      pt(8)=pt(8)+16*amm2*prr
      if(notd) go to 135
      pt(17)=pt(17)+64*amm3*prr
c
  135 hrr=4*pr(1)/rhh
      hr(1)=hrr
      do 136 i=2,4
  136 dxt(i)=-rho(i)
      dxt(3)=4+dxt(3)
      ii=4
      jj=10
      do 140 l=1,3
      l1=l+1
      hr(l1)=amm*hrr*dxt(l1)
      if(nosd) go to 140
      do 138 m=l,3
      ii=ii+1
      m1=m+1
      hr(ii)=amm*(hr(l1)*dxt(m1)-hrr*rho(ii))
      if(notd) go to 138
      do 137 n=m,3
      jj=jj+1
      ln=nuu(l,n)
      mn=nuu(m,n)
  137 hr(jj)=amm*(hr(ii)*dxt(n+1)-amm*hrr*dxt(m1)*rho(ln)-hr(l1)*rho(mn)
     .   -hrr*rho(jj))
  138 continue
  140 continue
      do 145 i=1,ider
  145 ht(i)=ht(i)+hr(i)
c  change to derivatives of log p
      ptt=pt(1)
      if(notd) go to 155
      jj=10
      do 152 l=1,3
      do 152 m=l,3
      do 152 n=m,3
      lm=nuu(l,m)
      ln=nuu(l,n)
      mn=nuu(m,n)
      jj=jj+1
  152 pt(jj)=(pt(jj)+(-(pt(lm)*pt(n+1)+pt(ln)*pt(m+1)+pt(mn)*pt(l+1))
     .   +2*pt(m+1)*pt(l+1)*pt(n+1)/ptt)/ptt)/ptt/amm
  155 ii=4
      do 158 l=2,4
      ptl=pt(l)/ptt/amm
      if(nosd) go to 158
      do 156 m=l,4
      ii=ii+1
  156 pt(ii)=(pt(ii)-pt(l)*pt(m)/ptt)/ptt/amm
  158 pt(l)=ptl
c
c  cp and dad
c
  160 pf=pt(2)
      hf=ht(2)
      dxtt=ht(3)*pf-hf*pt(3)
      lt=6
      do 165 l=1,3
      lf=l+4
      if(l.gt.1) lt=6+l
  165 dxt(l+1)=ht(lt)*pf+ht(3)*pt(lf)-ht(lf)*pt(3)-hf*pt(lt)
c
      cpp=dxtt/amm/t/pf
      cp(1)=cpp
      if(nosd) go to 173
      do 170 l=2,4
      dcp=dxt(l)/amm/t/pf-cpp*pt(l+3)/pf
      if(l.eq.3) dcp=dcp-cpp*amm
  170 cp(l)=dcp
c
  173 prh=amm*ptt/rhh
      anum=pf*prh-hf
      dad(1)=anum/dxtt
      if(nosd) go to 177
      do 175 l=2,4
      lf=l+3
  175 dad(l)=(prh*(amm*(pt(l)-rho(l))*pf+pt(lf))-ht(lf)-anum*dxt(l)
     .   /dxtt)/dxtt
c  further derivatives
  177 rhf=rho(2)
      dxtt=pt(3)*rhf-pf*rho(3)
      gm1=pf/(rhf-dad(1)*dxtt)
      tprh=rhf/dxtt
      trhp=-pf/dxtt
      rhxp=amm*x*(rho(4)-rhf*pt(4)/pf)
c  delta and derivatives
      delta=-1/trhp
      dlt(1)=delta
      lt=6
      if(nosd) go to 190
      do 180 l=2,4
      lf=l+3
      if(l.gt.2) lt=l+5
  180 dlt(l)=(pt(lt)*rhf+pt(3)*rho(lf)-pt(lf)*rho(3)-pf*rho(lt))/pf
     .   -delta*pt(lf)/pf
  190 xii1(1)=xii(1)
      xii1(2)=xii(11)
      xii1(3)=xii(21)
      xii1(4)=anuh(1)/anh0
      return
 
	end
