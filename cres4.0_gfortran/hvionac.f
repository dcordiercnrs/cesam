c
c*******************************************************************************
c
	subroutine hvionac(fl,tl,x,y,z,nosd,notd,anu,ue,anur,uer)
c
c  Calculate ionization of heavy elements.
c
c  Modification 17/5/90: include possibility of including full ionization
c     of heavy elements, for ihvz = 4. Note that this involves changing
c     size of arrays in commons /hvredu/ and /xhv/
c
c  Modification 4/6/90: include average degree of ionization for
c     each element in array xhvmn(10) at the beginning of
c     common/hviond/
c
      implicit double precision (a-h,o-z)
	implicit integer(i-n)
      logical nosd,notd
      character name*5
      dimension anu(1),ue(1),anur(1),uer(1)
      dimension eir(10),dr(10),hr(10),gr(10),xi(29),phi(30),hst(30)
      common/constsc/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev,exhm
      common/potetcc/ chi(125),am(10),iz(10)
      common/hvname/ name(10)
      common/hvabndc/ ab(10),iab
      common /hvredu/ irab,jfirst,ir(10),chir(125),izr(10),abr(10),
     *  amr(10)
      common /hvcntl/ icount,iwrite,dptst0,dptst1
      common/hviond/ xhvmn(10),xhv(125)
      common/eqscntc/ anz0,anze0,ihvz,iprrad,ihmin
	common/dmudec/ dmup(10),idmu
      common/eqdpco/ frhi,bdcoh,bdcoz,idpco
      common/ln10c/ amm,amm2,amm3
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
c
      save
c
c  test for restricting set of elements
c
      if(icount.gt.0) go to 10
      icount=1
      if(iwrite.eq.1) then
        write(istdpr,200)
        write(istdpr,205)
        lj=0
        do 3 k=1,10
        izj=iz(k)
        do 2 i=1,izj
        lj=lj+1
    2   xi(i)=chi(lj)
        im=min0(izj,13)
    3   write(istdpr,210) name(k),ab(k),izj,(xi(i),i=1,im)
      end if
c
c  test for restricting set of elements or levels included.
c  for ihvz = 4 keep everything, but for consistency storec
c  full set of information in restricted arrays
c
      if(ihvz.eq.1) then
c
c  reset to restricted set of variables
c
        irab=3
        ir(1)=1
        ir(2)=3
        ir(3)=10
c
c  reset abundances, to get reasonable fit to full case
c
        abr(1)=1.02*ab(1)
        abr(2)=1.35*ab(3)
        abr(3)=3.71*ab(10)
c
c  element with lowest ionization potential (here fe)
c
        jfirst=3
c
c  number of elements treated fully
c
        jiz1=2
c
      else
c
c  use full set of elements
c
        irab=10
        do 76 i=1,iab
        abr(i)=ab(i)
   76   ir(i)=i
c
c  lowest potential is of na
c
        jfirst=5
c
c  number of elements treated fully
c
        if(ihvz.eq.2) then
          jiz1=0
        else if(ihvz.eq.3) then
          jiz1=3
        else if(ihvz.eq.4) then
          jiz1=10
        else
          write(istdpr,110) ihvz
          if(istdou.ne.istdpr) write(istdou,110) ihvz
          stop
        end if
      end if
c
c  set possibly reduced set of ionization data
c
   79 j=1
      lj0=0
      ljr=0
      do 5 i=1,10
c
      izi=iz(i)
      if(i.eq.ir(j)) then
c
c  test for inclusion of all levels
c
        if(j.le.jiz1) then
          izrj=izi
        else
          izrj=1
        end if
        izr(j)=izrj
        amr(j)=am(i)
c
c  restorec ionization potentials
c
        lj=lj0
        do 4 k=1,izrj
        lj=lj+1
        ljr=ljr+1
    4   chir(ljr)=chi(lj)
        j=j+1
      end if
c
    5 lj0=lj0+izi
c
c  reset anz0 and anze0
c
      anz0=0
      anze0=0
      do 77 j=1,irab
      anz0=anz0+abr(j)*izr(j)/amr(j)
   77 anze0=anze0+abr(j)*izr(j)
c
      if(iwrite.eq.1) then
c
        write(istdpr,220)
        write(istdpr,205)
        lj=0
        do 7 k=1,irab
        izrj=izr(k)
        do 6 i=1,izrj
        lj=lj+1
    6   xi(i)=chir(lj)
        im=min0(izrj,13)
    7   write(istdpr,210) name(ir(k)),abr(k),izrj,(xi(i),i=1,im)
        write(istdpr,230) anz0
c
      end if
c
c  change to summed ionization potentials
c
      lj=0
      do 9 k=1,irab
      sum=0
      izj=izr(k)
      do 9 i=1,izj
      lj=lj+1
      sum=sum+chir(lj)
    9 chir(lj)=sum
c
c  set total number of levels for (possibly) reduced set
c
      nlvtot=lj
c
c  ***************************************
c
c  end of initialization section
c
   10 f=1.d1**fl
      t=1.d1**tl
c  test for departure of heavy elements
      if(idpco.lt.2) bdcoz=1
c  k*t, in ev
      tk=ck1*t
c
c  test whether phderc has already been called
c
      if(idmu.eq.1) go to 15
c
      call phderc(fl,tl,phi,hst,nosd,notd)
c
c
c  psi and derivatives
      wf= sqrt(1+f)
      psi=2*wf+log(f/((1+wf)*(1+wf)))
      psif=amm*wf
      psff=amm2*f/2/wf
c
      zf=x+2*y+anze0*z
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
      tki=1.d0/tk
      bdcozl=log(bdcoz)
c
c  ***************************************
c
c  start loop over elements
c
      do 50 j=1,irab
      izoj=iz(ir(j))
      izj=izr(j)
c..      write(6,*) 'j, izj, izoj', j, izj, izoj
c
c  skip elements with zero abundance
c
      if(abr(j).le.0) then
        anur(j)=0
        uer(j)=0
        xhvmn(j)=0
        lj=lj+izj
        go to 50
      end if
c
c  set limit for no ionization
c
      dptstn=dptst1
      if(j.eq.jfirst) dptstn=dptst0
c
c  set exponents phi in saha equations into xi
c
      lj0=lj
      do 20 i=1,izj
      lj=lj+1
      xhv(lj)=0
      dxp=i*dmps-chir(lj)*tki+bdcozl
   20 xi(i)=dxp
c
c  -----------------
c
c  test for complete or no ionization
c
      if(izj-1) 21,21,25
c
c  only one level
c
   21 phm=xi(1)
      if(phm) 22,22,23
c
   22 imax=0
      phm=0
c  test for no ionization
      if(xi(1).lt.-dptstn) go to 50
c
      go to 34
c
   23 imax=1
c  test for complete ionization
      if(xi(1).gt.dptst1) go to 29
      go to 34
c
c  more than one level
c
   25 ii=izj+1
      xi(ii)=0
      imax=1
      phm=xi(1)
c
c  set phm to largest xi
c
      do 26 i=2,ii
      if(xi(i).le.phm) go to 26
      imax=i
      phm=xi(i)
   26 continue
c
      if(imax.ne.izj) go to 30
c
c  test for complete ionization
c
      izj1=izj-1
      dphm=phm-xi(1)
      if(izj1.eq.1) go to 28
      do 27 i=1,izj1
   27 dphm=min(dphm,phm-xi(i))
   28 if(dphm.le.dptst1) go to 34
c
c  complete ionization
c
   29 xhv(lj)=1
c..      write(6,*) 'complete ionization at lj =',lj
      fct1=abr(j)/amr(j)
      anur(j)=fct1*izj
      uer(j)=fct1*chir(lj)
      anu(1)=anu(1)+anur(j)
      ue(1)=ue(1)+uer(j)
      xhvmn(j)=1
      go to 50
c
   30 if(imax.ne.ii) go to 34
c
c  test for no ionization
c
      do 31 i=1,izj
      if(xi(i).gt.-dptstn) go to 34
   31 continue
c
c  no ionization. skip element completely
c
      xhvmn(j)=0
c
      go to 50
c
c  ************************************
c
c  general case
c
   34 do 35 i=1,idmx
      dr(i)=0.e0
      hr(i)=0.e0
   35 gr(i)=0.e0
      if(phm.le.dptst1) dr(1)=omegac(0,izoj)*exp(-phm)
      do 40 i=1,izj
      lji=lj0+i
      cchi=chir(lji)
      dxp=xi(i)-phm
      if(dxp.lt.-dptst1) go to 40
      dxp=omegac(i,izoj)*exp(dxp)
c..      print *,' j,i,dxp',j,i,dxp
      xhv(lji)=dxp
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
      dr1i=1.d0/dr(1)
c  scale xhv
      do 42 i=1,izj
      lji=lj0+i
      xhv(lji)=dr1i*xhv(lji)
c..      write(6,*) 'At i, j, lji =',i,j,lji,'  xhv =',xhv(lji)
   42 continue
c
      fct1=abr(j)/(amr(j)*dr(1))
      hrdr=hr(1)*dr1i
      grdr=gr(1)*dr1i
      anur(j)=fct1*hr(1)
      uer(j)=fct1*gr(1)
      anu(1)=anu(1)+anur(j)
      ue(1)=ue(1)+uer(j)
c
      xhvmn(j)=hr(1)/(izj*dr(1))
c
c..      print *,' j,anur(j),anu(1)',j,anur(j),anu(1)
c  derivatives
      ii=4
      do 48 k=2,4
      anu(k)=anu(k)+fct1*(hr(k)-dr(k)*hrdr)
      ue(k)=ue(k)+fct1*(gr(k)-dr(k)*grdr)
      if(nosd) go to 48
      do 45 l=k,4
      ii=ii+1
      anu(ii)=anu(ii)+fct1*(hr(ii)-(dr(k)*hr(l)+dr(l)*hr(k)+(dr(ii)
     .  -2*dr(k)*dr(l)*dr1i)*hr(1))*dr1i)
   45 ue(ii)=ue(ii)+fct1*(gr(ii)-(dr(k)*gr(l)+dr(l)*gr(k)+(dr(ii)
     .  -2*dr(k)*dr(l)*dr1i)*gr(1))*dr1i)
   48 continue
   50 continue
c..      if(idgeos.ge.3) write(istdpr,250) (xhv(i),i=1,nlvtot)
c
      return
  110 format(//' **** error in s/r hvionac. ihvz =',i5,
     *  ' not allowed')
  200 format(///' Original heavy element data.'/)
  205 format(' Name, abundance, Z, ionization potentials:'/)
  210 format(1x,a4,f8.5,i4,13f8.2)
  220 format(///' Heavy element data after resetting:'/)
  230 format(//' now anz0 =',f10.5//)
  250 format(/' xhv:'/(1p10e12.4))
      end
