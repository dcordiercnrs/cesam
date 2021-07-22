c
c******************************************************************************
c
	subroutine eqstfc(fl,tl,x,y,z,nosd,notd)
c
c  This version of the EFF routine includes
c  effects of Coulomb interaction in the Debye-Hueckel approximation,
c  using the F4 from the MHD package.
c
c  Modifications initiated 19/4/90.
c
c  jc-d equation of state updated to be compatible with  dog equation
c  of state routines. notice in particular that common/eqscntc/ has
c  been changed, by the inclusion of anhe0 after anh0. in addition  a
c  dummy sbroutine seteqs has been added. furthermore the dimensions
c  in common /potetc/ have been changed, to allow for only up to
c  10 heavy elements.
c
c  controls:
c  ********
c
c  in argument list:
c
c  logical nosd and notd are switches for no calculation of second and
c  third derivatives, respectively.
c
c  in common/eqscntc/:
c
c  ihvz determines treatment of heavy elements.
c  ihvz = 0: assume heavy elements to be fully ionized everywhere
c       = 1: complete treatment of ionization of c and o; first level
c            of ionization of fe. reset abundances to fit results of
c            full treatment of all 10 elements considered.
c       = 2: first level of ionization of all 10 elements.
c       = 3: complete treatment of c, n and o; first level of
c            ionization of remaining elements.
c       = 4: complete treatment of ionization of all 10 elements
c            included.
c
c  iprrad .ne. 0: include pressure and enthalpy of radiation (in
c                 diffusion approximation)
c  iprrad  =   0: do not include pressure and enthalpy of radiation.
c
c  ihmin   =   1: include h-
c         .ne. 1: do not include h-
c
c  these are initialized  in block data bleqstc to ihvz = 1,
c  iprrad = 1, ihmin = 0.
c
c  in common/eqdpco/, relating to fudge departure coefficients to
c  mimick nlte effects:
c
c  idpco = 0: do not include departure coefficients
c        = 1: determine departure coefficient for h (returned in
c             bdcoh) such that n(h+)/n(h neutr.) = frhi.
c        = 2: treat h as for idpco = 1. take departure coefficient
c             for he and heavy elements to be bdcoz.
c        = 3: take departure coefficient for h to be bdcoh, and
c             departure coefficient for he and heavy elements to be
c             bdcoz.
c
c  frhi, bdcoh and bdcoz are initialized to 1., and idpco to 0.
c
c  addition 5/1/84: flag iomfll in common/hvomcl/.
c    if iomfll = 0 statistical weight omegac for fully ionized state
c    of heavy elements is set to 15.
c    this corresponds to the situation before 5/1/84. the
c    origin of this error is unclear.
c    if iomfll .ne. 0 omegac for fully ionized state is set to 1,
c    as it should be.
c
c  Controls in common /ccoulm/
c    epsdmu: convergence criterion for iteration for Coulomb
c    effect in Saha equation. (Default: 1.e-12)
c    icoulm: determines how Coulomb terms are included:
c      icoulm = -1: do not include Coulomb effects (should correspond
c         to basic eqstf).
c      icoulm = 0: include Coulomb effects fully.
c      icoulm = 1: include Coulomb effects in Saha equations, but not
c         explicit contributions to pressure and enthalpy.
c      icoulm = 2: include explicit contributions to pressure and
c         enthalpy from Coulomb effects, but not effect
c         in Saha equations.
c      (Default: 0)
c    iclmit: determines type of Coulomb iteration
c      iclmit = 0: backsubstitution
c      iclmit = 1: Newton iteration
c      (Default: 1)
c    iclmsd; if iclmsd = 1, set include second derivatives of
c      Coulomb terms
c      (Default: 0)
c
c
c  modification 15/9/86: save statement added in all sbroutines.
c     (needed in f77, apparently)
c
c  modification 3/1/88: dimension of dxt1 increased to 20.
c
c  Modification 6/8/88: second derivative of He ionization fractions
c     corrected. Also in commom/eqsout/ is now storecd ea before rescaling,
c     for test of derivatives.
c
c  Modification 17/5/90: include possibility of including detailed
c     ionization of heavy elements, for ihvz = 4. Note that this
c     involves changing size of arrays in commons /hvredu/ and /xhv/
c
c  Modified 22/5/90 to include all levels of ionization of Ar and Fe.
c     Previously only the first 15 levels were included, the remainder
c     being forced to be unionized
c
c  Modification 4/6/90: include average degree of ionization for
c     each element in array xhvmn(10) at the beginning of
c     common/hviond/
c
c  Modified 4/6/90, adding array gmm1(4) in common/ eqstd/ containing
c     Gamma1 and derivatives.
c
c  Modified 4/6/90, to include numerical second derivatives of
c     Coulomb terms.
c
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 9/3/90
c
c            *********************************
c
c  results returned in common/eqstdc/. in addition value of fl
c  is set into common/eqstfl/.
c
      implicit double precision (a-h,o-z)
	implicit integer(i-n)
      parameter(nspe = 6, npar = 3, npar2=npar+2)
      logical tstl,nosd,notd,cmplio, secder, thrder
      dimension phi(30),hst(10),ex(3),ea(30),xi(30),dxt(4),dxt1(20),
     .  dxt2(4),anuh(10),ueh(10),anuhr(23),uehr(23),
     *  ehm(10),aneder(10), sn(nspe), dmucpr(npar),
     *  eahat(30), dmucc(npar,4), dmucp(npar,4,6),
     *  pcoulc(4), pcoulp(4,6), hcoulc(4), hcoulp(4,6)
	common/dmudec/ dmu,dmuf,dmut,dmux,dmuff,dmuft,dmufx,dmutt,dmutx,
     1	dmuxx,idmu
      common/eqscntc/ anh0,anhe0,ihvz,iprrad,ihmin
      common/ln10c/ amm,amm2,amm3
      common/constsc/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev,exhm
c
c  note: second dimension of dmuc must be at least as great as npar + 2
c
      common/dmucdr/ dmuc(npar,10)
      common/df4der/ df4(npar), d2f4i(npar,npar), d2f4f(npar,npar),
     .          d2f4x(npar),d2f4r(npar),d2f4t(npar),
     .          dp4i(npar),dh4i(npar), p4, h4, dp4dr, dp4dt, dp4dx,
     .          dh4dr, dh4dt, dh4dx
      common/eqsoutc/ east(30),xii(30),dne(10),dph(20),he(20),pe(10),
     *  hi(20),pi(10),hh(20),ph(10),hr(20),pr(10), pcoul(10), hcoul(10)
      common/eqstdc/ xii1(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     *  dlt(4),gm1,tprh,trhp,rhxp,gmm1(4)
      common/eqstfl/ flcm
      common/eqdpco/ frhi,bdcoh,bdcoz,idpco
      common/xhminc/ xhm(10)
c
c  controls for Coulomb effect.
c  epsdmu is accuracy requirement in backsubstitution
c  icoulm is used for selectively switching off parts of Coulomb effect
c
      common/ccoulm/ epsdmu, icoulm, iclmit, iclmsd, nitdmu, epssdr
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
      common/cversn/ ivrbcs,ivrrhs,ivreos,ivropc,ivreng
c
      external bleqstc, blstio
c
      equivalence(ex(1),exh)
c
      data flprev, tlprev, xprev /-1.d10, -1.d10, -1.d10/
c
      save
c
      fxp(a)=exp(min(85.d0,max(a,-85.d0)))
      nuu(l,n)=((l-1)*(8-l))/2+n-l+5
c
c  *****************************************************************
c
c  set equation of state version number, depending on icoulm
c
      if(icoulm.eq.0) then
        ivreos=1
      else if(icoulm.eq.2) then
        ivreos=2
      end if
c
c  *****************************************************************
c
c  for changing internal logics, set logical variables for setting
c  second and third derivatives
c
      secder = .not. nosd
      thrder = .not. notd
c
c  set y
      y=1-x-z
c  storec fl in common
      flcm=fl
c  number of derivatives
      if(thrder) then
        ider=20
      else if(secder) then
        ider=10
      else
        ider=4
      end if
      iders=min0(ider,10)
c
      iclsdr = 0
c
c  entry point for setting numerical second derivatives of
c  Coulomb terms
c
    5 f=1.d1**fl
      t=1.d1**tl
c
      call phderc(fl,tl,phi,hst,nosd,notd)
c
c
c  psi and derivatives
      wf= sqrt(1+f)
      psi=2*wf+log(f/((1+wf)*(1+wf)))
      psif=amm*wf
      psff=amm2*f/(2*wf)
c
      zf=x+2*y+anhe0*z
      zf2=zf*zf
      zf3=zf*zf2
c
      ak0=ckh*zf2
c
c  set k*t, in ev and ergs
c
      tkev=ck1*t
      tkergs=ck2*t
c  a0**3*ne/(k*t)
      aa=caa*phi(1)/(tkev*zf3)
c
c  delta mu and derivatives
c
      bmu=tkev+20*ak0
      dmu=aa*bmu
c
      ref=phi(2)
      ret=phi(3)
      dmuf=dmu*ref*amm
      dmut=aa*(bmu*ret-20*ak0)*amm
      dmux=aa*(3*tkev+20*ak0)/zf
c
      if(secder) then
        dmuff=(dmuf*ref+dmu*phi(4)*amm)*amm-psff
        dmuft=(dmut*ref+dmu*phi(5)*amm)*amm
        dmufx=dmux*ref*amm
        dmutt=aa*(20*ak0*(1-2*ret)+bmu*(phi(6)+ret*ret))*amm2
        dmutx=aa*(ret*(3*tkev+20*ak0)-20*ak0)*amm/zf
        dmuxx=aa*(12*tkev+40*ak0)/zf2
      end if
c
      dmu=dmu-psi
      dmuf=dmuf-psif
      idmu=1
c
c  ***************************************
c
c  ionization of heavy elements
c
      if(ihvz.le.0) then
        anuh(1)=anh0
        call zeroc(anuh(2),9)
        call zeroc(ueh,10)
c
      else
c
        call hvionac(fl,tl,x,y,z,nosd,notd,anuh,ueh,anuhr,uehr)
c
      end if
c
c  test for complete ionization of h and he
c
      cmplio=(dmu-54.4/tkev).gt.19
c  test for use of departure coefficient
      if(idpco .ge.1) cmplio=cmplio.and.frhi.gt.1.e6
      if(cmplio) then
c
c  complete ionization
c
        xi(1) =1
        xi(11)=0
        xi(21)=1
        do 10 i=2,22,10
        jj=i+8
        do 10 j=i,jj
   10   xi(j)=0
c
c  to avoid problems with logics in Coulomb term treatment,
c  initialize iclder
c
        iclder = 0
c
        go to 30
c
      end if
c
c  A note on logics of Coulomb iteration: iclder = 0 flags for
c  initial passes iterating for Coulomb correction, and iclder = 1
c  flags for final pass setting derivatives
c  Similarly, when setting second derivatives numerically, iclsdr
c  flags for modifications of variables
c
c  e h, e he, e he+ and derivatives
c
c  test for initializing Coulomb terms or using previous values
c  Note: for icoulm = 2 switch off effect of Coulomb term on ionization
c  For icoulm = -1 switch off effect of Coulomb term entirely
c  Initialization also depends on whether Newton iteration or
c  backsubstitution is used
c
      if(icoulm.eq.2.or.icoulm.eq.-1) then
        call zeroc(dmuc,npar*10)
        call zeroc(dmucpr,npar)
      else
c
c  initialize for Coulomb iteration. Switch off second derivatives
c  for iteration passes.
c
        iclder=0
        secder=.false.
        call zeroc(dmuc,npar*10)
        if(abs(fl-flprev).gt.0.1.or.abs(tl-tlprev).gt.0.1
     *    .or.abs(x-xprev).gt.0.1) then
          call zeroc(dmucpr,npar)
        end if
      end if
c
      flprev=fl
      tlprev=tl
      xprev=x
      nitdmu=0
c
c  entry point for setting Coulomb derivatives
c
   15 k=-10
c..      write(6,*) 'iclsdr, fl, tl, x', iclsdr, fl, tl, x
      do 25 ia=1,3
      k=k+10
      ext=ex(ia)/tkev
      eea=dmu-ext+dmuc(ia,1)
      if(eea+30.le.0) then
c
c  no ionization
c
        do 16 i=1,10
   16   ea(k+i)=0.e0
        go to 25
c
      end if
c
      eea=fxp(eea)
c
c  set statistical weights
c
      if(ia.ne.2) then
        eea=eea/2
      else
        eea=2*eea
      end if
c
c  test for departure coefficient for h or he
c
      if(ia.eq.1.and.idpco.gt.0) then
        if(idpco.le.2) bdcoh=frhi/eea
        eea=bdcoh*eea
      else if(ia.ge.2.and.idpco.ge.2) then
        eea=bdcoz*eea
      end if
c
      ea(k+1)=eea
c  first derivatives
      ea(k+2)=(dmuf+dmuc(ia,2))*eea
      ea(k+3)=(amm*ext+dmut+dmuc(ia,3))*eea
      ea(k+4)=(dmux+dmuc(ia,4))*eea
      if(secder) then
c  second derivatives, with Coulomb contribution
        ea(k+5)=(dmuff+dmuc(ia,5))*eea+(dmuf+dmuc(ia,2))*ea(k+2)
        ea(k+6)=(dmuft+dmuc(ia,6))*eea+(dmuf+dmuc(ia,2))*ea(k+3)
        ea(k+7)=(dmufx+dmuc(ia,7))*eea+(dmux+dmuc(ia,4))*ea(k+2)
        ea(k+8)=(dmutt+dmuc(ia,8)-amm2*ext)*eea
     *         +(amm*ext+dmut+dmuc(ia,3))*ea(k+3)
        ea(k+9)=(dmutx+dmuc(ia,9))*eea+(dmux+dmuc(ia,4))*ea(k+3)
        ea(k+10)=(dmuxx+dmuc(ia,10))*eea+(dmux+dmuc(ia,4))*ea(k+4)
      end if
c
   25 continue
c
c  test for setting up for Coulomb iteration
c
      call storec(ea,eahat,30)
c
c  in iteration pass, include dmuc from previous case in ea
c  unless ignoring Coulomb effect on ionization
c
      if(iclder.ne.1.and.icoulm.ne.-1.and.icoulm.ne.2) then
        do 26 i=1,3
        is=1+10*(i-1)
   26   ea(is)=eahat(is)*exp(dmucpr(i))
      end if
c
c  entry point for continuing Coulomb iteration by backsubstitution
c  For subsequent iterations, include Coulomb term
c
   27 if(iclder.ne.1.and.nitdmu.gt.0) then
        do 28 i=1,3
        is=1+10*(i-1)
   28   ea(is)=eahat(is)*exp(dmuc(i,1))
      end if
c
c  entry point for continuing Coulomb iteration by Newton iteration
c  for diagnostic or iteration purposes, storec current ea
c
   29 call storec(ea,east,30)
c
c  set degrees of ionization of H and He
c
      call hheion(ea(1), ea(11), ea(21), xi(1), xi(11), xi(21), secder)
c
c  inclusion of h-
c
   30 if(ihmin.eq.1) then
        call hmnion(tl, ea(1), ehm, xhm, secder)
      else
        call zeroc(xhm,10)
      end if
c
c  ne and derivatives
c
c  combine he fractions for an in xi(20+.), and for ionization
c  energy in xi(10+.)
c
      exhc=exhe+exhep
c
      if(secder) then
        imx=10
      else
        imx=4
      end if
      do 35 i=1,21,10
   35 call storec(xi(i),xii(i),imx)
      xia=xi(1)
      imx=imx+20
c
      do 36 i=21,imx
      i10=i-10
      xio=exhe*xi(i10)+exhc*xi(i)
      xi(i)=2*xi(i)+xi(i10)
   36 xi(i10)=xio
c
      ider2=min0(ider,10)
      if(ihmin.eq.1) then
c
c  in the case with H- there is a risk of negative (unphysical)
c  Ne. Write error message and reset to values without
c  combine h and h-
c
        dxt1(1)=xi(1)-xhm(1)
        ane(1)=av*(dxt1(1)*x/ah+xi(21)*y/ahe+z*anuh(1))
        if(ane(1).le.0) then
          write(istdpr,1010) xi(1), xhm(1), anuh(1)
          do 37 l=1,ider2
   37     dxt1(l)=xi(l)
          ane(1)=av*(dxt1(1)*x/ah+xi(21)*y/ahe+z*anuh(1))
        else
          do 38 l=2,ider2
   38     dxt1(l)=xi(l)-xhm(l)
        end if
      else
c
c  no H-
c
        do 39 l=1,ider2
   39   dxt1(l)=xi(l)
c
        ane(1)=av*(dxt1(1)*x/ah+xi(21)*y/ahe+z*anuh(1))
c
      end if
c
c  terms needed for x-derivatives
c
      do 40 l=1,4
   40 dxt(l)=av*(dxt1(l)/ah-xi(20+l)/ahe)
c
      ii=4
      do 44 l=1,3
      l1=l+1
      anel=(dxt1(l1)*x/ah+xi(l+21)*y/ahe+z*anuh(l1))*av
      tstl=l.eq.3
      if(tstl) anel=anel+dxt(1)
      ane(l1)=anel
      if(secder) then
        do 42 m=l,3
        m1=m+1
        ii=ii+1
        anelm=(dxt1(ii)*x/ah+xi(20+ii)*y/ahe+z*anuh(ii))*av
        if(tstl) anelm=anelm+dxt(m1)
        if(m.eq.3) anelm=anelm+dxt(l1)
   42   ane(ii)=anelm
      end if
   44 continue
c
c  as auxiliary quantities, to avoid over- and underflow, set
c  derivatives of ne divided by ne
c
      do 46 i=2,iders
   46 aneder(i)=ane(i)/ane(1)
c
c  the density and derivatives (note that rho(2) = dlog rho/dlog f,
c  and so on).
c
      anee=ane(1)
      rho(1)=phi(1)*(crho/anee)
      ii=4
      jj=10
      do 50 l=1,3
      l1=l+1
      rhol=-aneder(l1)/amm
      tstl=l.le.2
      if(tstl) rhol=rhol+phi(l1)
      rho(l1)=rhol
      if(secder) then
        do 48 m=l,3
        ii=ii+1
        m1=m+1
        lm=l+m
        rholm=(aneder(l1)*aneder(m1)-aneder(ii))/amm
        if(tstl.and.m.le.2) rholm=rholm+amm*phi(2+lm)
        if(thrder) then
          do 47 n=m,3
          jj=jj+1
          rhd=-2*aneder(l1)*aneder(m1)*aneder(n+1)/amm
          if(l.lt.3.and.m.lt.3.and.n.lt.3) rhd=rhd+amm2*phi(4+lm+n)
   47     rho(jj)=rhd
        end if
   48   rho(ii)=rholm
      end if
   50 continue
c
c  test for skipping Coulomb terms entirely
c  Also skip this bit for final pass when setting second derivatives
c
      if(icoulm.eq.-1.or.iclsdr.eq.7) go to 60
c
c  set Coulomb terms. Prepare for calling WD f4 routine
c  Note: we currently skip this in the pass for setting derivatives,
c  unless for full ionization, where Coulomb part has not been passed
c  previously
c  If not setting derivatives in f4der is implemented in future,
c  this may require further thought
c
      if(iclder.ne.1.or.cmplio.or.icoulm.eq.2) then
c
        anx=x*av/ah
        any=y*av/ahe
c
        sn(1)=(1.d0-xii(1))*anx
        sn(2)=xii(1)*anx
        sn(3)=(1.d0-xii(11)-xii(21))*any
        sn(4)=xii(11)*any
        sn(5)=xii(21)*any
c
c  to avoid rounding error problems in subtractions in sn(1) and sn(3),
c  reset them from the ea-s, or to zero for full ionization
c
        if(cmplio) then
          sn(1)=0
          sn(3)=0
        else
          sn(1)=anx/(1+east(1))
          sn(3)=any/(1+east(11)*(1+east(21)))
        end if
c
c  number of free electrons. For consistency include only those
c  coming from H and He
c
        sn(6)=xii(1)*anx+(xii(11)+2*xii(21))*any
c
        rhl=log10(rho(1))
c        if(idgeos.ge.3) write(istdpr,*) 'calling f4der at',rhl,tl,sn
c
        call f4der(rhl,tl,sn,nspe,f4,e4,p4,h4,d2f4r2,d2f4t2,
     .             d2f4rt,df4,d2f4i,d2f4f,d2f4x,d2f4r,d2f4t,
     .             dp4i,dp4dr,dp4dt,dp4dx,dh4i,dh4dr,dh4dt,dh4dx,npar)
c
        if(icoulm.ne.2) then
          dmuc(1,1)=-df4(1)/tkergs
          dmuc(2,1)=-df4(2)/tkergs
          dmuc(3,1)=-df4(3)/tkergs
c
c  test for continuing Coulomb iteration
c
          ddmuc1=dmuc(1,1)-dmucpr(1)
          ddmuc2=dmuc(2,1)-dmucpr(2)
          ddmuc3=dmuc(3,1)-dmucpr(3)
c          if(idgeos.ge.2)
c	write(6,*) ' ddmuc1-3:', ddmuc1, ddmuc2, ddmuc3
c        else
          ddmuc1=0
          ddmuc2=0
          ddmuc3=0
        end if
c
      end if
c
c..      write(6,*) 'iclder =',iclder
c
c  test for convergence
c
      if(.not.cmplio.and.icoulm.ne.2.and.iclder.eq.0.and.
     *   (abs(ddmuc1).ge.epsdmu.or.
     *    abs(ddmuc2).ge.epsdmu.or.
     *    abs(ddmuc3).ge.epsdmu)) then
c
c  test for failure to converge
c
        if(nitdmu.gt.20) then
          write(istdpr,1100) fl, tl, x, z,
     *      ddmuc1, ddmuc2, ddmuc3
          if(istdou.ne.istdpr) write(istdou,1100) fl, tl, x, z,
     *      ddmuc1, ddmuc2, ddmuc3
          stop
        end if
c
c
        nitdmu=nitdmu+1
c
c  storec previous value
c
        do 51 i=1,npar
   51   dmucpr(i)=dmuc(i,1)
c
c  test for backsubstitution or Newton iteration
c
        if(iclmit.eq.0) then
          go to 27
        else
c
c  storec derivatives in dmuc
c
c..          write(6,*) 'd2f4f(1,.):',(d2f4f(1,j),j=1,3)
c..          write(6,*) 'd2f4f(2,.):',(d2f4f(2,j),j=1,3)
c..          write(6,*) 'd2f4f(3,.):',(d2f4f(3,j),j=1,3)
          do 52 i=1,npar
          dmuc(i,2)=-d2f4r(i)/tkergs
          do 52 j=1,npar
   52     dmuc(i,j+2)=-d2f4f(i,j)/tkergs
c
          if(idgeos.ge.4) write(6,*) 'xi(1-3):',xii(1), xii(11), xii(21)
          call clmnew(east,eahat,dmuc,ane(1),x,y,ea,npar,nitdmu)
          go to 29
        end if
c
c
c  test for needing last pass of ionization section, to set
c  derivatives
c
      else if(.not.cmplio.and.icoulm.ne.2.and.iclder.eq.0) then
c
c  set derivatives of delta mu c
c
        do 54 i=1,npar
        do 53 j=2,4
        dmuc(i,j)=0
        do 53 k=1,npar
   53   dmuc(i,j)=dmuc(i,j)-d2f4i(i,k)*xii(10*(k-1)+j)
c
c  add explicit contributions from derivatives of
c  dmuc wrt log rho, log T and X at fixed xi
c  Note: last term in dmuc(i,3) corrects for division by kT.
c
        dmuc(i,2)=dmuc(i,2)-d2f4r(i)*rho(2)
        dmuc(i,3)=dmuc(i,3)-d2f4r(i)*rho(3)-d2f4t(i)+amm*df4(i)
c..        write(6,*) 'dmucx:',i,dmuc(i,4),-d2f4r(i)*rho(4),-d2f4x(i)
        dmuc(i,4)=dmuc(i,4)-d2f4r(i)*rho(4)-d2f4x(i)
c
        do 54 j=2,4
   54   dmuc(i,j)=dmuc(i,j)/tkergs
c
        iclder=1
c
c  restorec secder, unless a later pass will be made to compute
c  second derivatives
c
        secder=.not.nosd
c
        go to 15
c
      else
c
c  otherwise storec pressure and enthalpy from Coulomb effect,
c  transforming derivatives
c
        pcoul(1)=p4
        hcoul(1)=h4
c
c..        write(6,*) 'setting pcoul(1) =',pcoul(1)
c
        call zeroc(pcoul(2),9)
        call zeroc(hcoul(2),9)
c
        do 56 j=2,4
        do 56 k=1,npar
        pcoul(j)=pcoul(j)+dp4i(k)*xii(10*(k-1)+j)
   56   hcoul(j)=hcoul(j)+dh4i(k)*xii(10*(k-1)+j)
c
c  add explicit contributions from derivatives of
c  dmuc wrt log rho, log T and X at fixed xi
c
        pcoul(2)=pcoul(2)+dp4dr*rho(2)
        pcoul(3)=pcoul(3)+dp4dr*rho(3)+dp4dt
        pcoul(4)=pcoul(4)+dp4dr*rho(4)+dp4dx
        hcoul(2)=hcoul(2)+dh4dr*rho(2)
        hcoul(3)=hcoul(3)+dh4dr*rho(3)+dh4dt
        hcoul(4)=hcoul(4)+dh4dr*rho(4)+dh4dx
c..        write(6,*) 'In eqstf dh4dr, dh4dt =',dh4dr, dh4dt
      end if
c
c  test for setting second derivatives of Coulomb terms numerically
c
      if(iclmsd.eq.1.and..not.nosd) then
c
        if(iclsdr.eq.0) then
          flc=fl
          tlc=tl
          xc=x
          if(icoulm.ne.2) call storec(dmuc,dmucc,4*npar)
          call storec(pcoul,pcoulc,4)
          call storec(hcoul,hcoulc,4)
c
          fl=flc-epssdr
          iclsdr=1
          go to 5
c
        else
          if(icoulm.ne.2) call storec(dmuc,dmucp(1,1,iclsdr),4*npar)
          call storec(pcoul,pcoulp(1,iclsdr),4)
          call storec(hcoul,hcoulp(1,iclsdr),4)
c
          if(iclsdr.le.5) then
            if(iclsdr.eq.1) then
              fl=flc+epssdr
            else if(iclsdr.eq.2) then
              fl=flc
              tl=tlc-epssdr
            else if(iclsdr.eq.3) then
              tl=tlc+epssdr
            else if(iclsdr.eq.4) then
              tl=tlc
              x=xc-epssdr
              y=1-x-z
            else if(iclsdr.eq.5) then
              x=xc+epssdr
              y=1-x-z
            end if
c
            iclsdr=iclsdr+1
            go to 5
c
          else
c
            x=xc
            y=1-x-z
c
            if(icoulm.ne.2) call storec(dmucc,dmuc,4*npar)
            call storec(pcoulc,pcoul,4)
            call storec(hcoulc,hcoul,4)
c
c  set second derivatives
c
            epsdi2=1.d0/(2*epssdr)
            pcoul(5) =(pcoulp(2,2)-pcoulp(2,1))*epsdi2
            hcoul(5) =(hcoulp(2,2)-hcoulp(2,1))*epsdi2
            pcoul(6) =(pcoulp(3,2)-pcoulp(3,1))*epsdi2
            hcoul(6) =(hcoulp(3,2)-hcoulp(3,1))*epsdi2
            pcoul(7) =(pcoulp(4,2)-pcoulp(4,1))*epsdi2
            hcoul(7) =(hcoulp(4,2)-hcoulp(4,1))*epsdi2
            pcoul(8) =(pcoulp(3,4)-pcoulp(3,3))*epsdi2
            hcoul(8) =(hcoulp(3,4)-hcoulp(3,3))*epsdi2
            pcoul(9) =(pcoulp(3,6)-pcoulp(3,5))*epsdi2
            hcoul(9) =(hcoulp(3,6)-hcoulp(3,5))*epsdi2
            pcoul(10)=(pcoulp(4,6)-pcoulp(4,5))*epsdi2
            hcoul(10)=(hcoulp(4,6)-hcoulp(4,5))*epsdi2
c
            if(icoulm.ne.2) then
              do 58 i=1,npar
              dmuc(i,5) =(dmucp(i,2,2)-dmucp(i,2,1))*epsdi2
              dmuc(i,6) =(dmucp(i,3,2)-dmucp(i,3,1))*epsdi2
              dmuc(i,7) =(dmucp(i,4,2)-dmucp(i,4,1))*epsdi2
              dmuc(i,8) =(dmucp(i,3,4)-dmucp(i,3,3))*epsdi2
              dmuc(i,9) =(dmucp(i,3,6)-dmucp(i,3,5))*epsdi2
   58         dmuc(i,10)=(dmucp(i,4,6)-dmucp(i,4,5))*epsdi2
c
c  make final pass of ionization section to set second derivatives
c  correctly
c
              iclder=1
              secder=.true.
              iclsdr=7
              go to 15
            end if
c
          end if
c
        end if
c
      end if
 
c
c  ****************************************************
c
c  end skipping Coulomb terms
c
   60 continue
c
c  start setting total pressure and enthalpy
c
c  as initialization, zero arrays
c
      call zeroc(pt,20)
      call zeroc(ht,20)
c
c
c  delta p, delta h and derivatives.
c
c  test for complete ionization
      if(cmplio) go to 80
c
c  delta ne**2
c
c  note: on 3/1/84 definition of dne was changed to
c        (delta ne **2)/av**2 to avoid overflows in calculation
c        of delta p and delta h on univac (with exponents limited
c        to 38). a corresponding change was made in the definition
c        of ca03 in s/r setcnsc.
c
c
c  reset xi(1) and xi(21) (note that xi(2-10) and xi(22-30) are still
c  consistent).
      xi(1)=-1/(east(1)+1)
      xi(21)=-(2+east(11))/(1+east(11)*(1+east(21)))
c
      do 65 l=1,4
      dxt1(l)=xi(l+20)/ahe-xi(l)/ah
      dxtl=-xi(l+20)*y/ahe-xi(l)*x/ah-anuh(l)*z
      if(l.eq.1) dxtl=dxtl+anh0*z
      dxtl2=ane(l)
      if(l.le.3) go to 64
      dxtl=dxtl+dxt1(1)
      dxtl2=dxtl2+avda
   64 dxt(l)=dxtl
   65 dxt2(l)=dxtl2/av
c
c  ne0/av:
c
      ann=ane(1)/av+(x/ah+2*y/ahe+z*anh0)
      xtt=dxt(1)
      xttiav=xtt/av
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
   68 dne(ii)=ane(ii)*xttiav+dxt2(l1)*dxt(m1)+dxt2(m1)*dxt(l1)
     .  +ann*dnlm
   70 continue
c
c  quantities common to delta p and delta h (dph(15-20) is used as
c  intermediate storage).
c
c  Note: here the constant for delta p and delta H is set to C1 = ca03/2
c
      a03=ca03/zf3/2
c
      call zeroc(dph(15),6)
      dxt1(1)=0
      c1=amm*tkev
      dxt1(2)=c1
      dph(18)=amm2*tkev
      dph(19)=3*c1/zf
      dnee=dne(1)
      a03=a03*rho(1)*rho(1)
      nu=2
      nb=20
      k1=1
   75 do 77 l=2,4
   77 dxt(l-1)=dne(l)+nu*amm*dnee*rho(l)
c
      dxt1(3)=(3*tkev+nb*ak0)/zf
      bmu=tkev+nb*ak0
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
     .  +dnee*(rho(jj)+nu*amm*rho(m1)*rho(l1))))+dxt(l)*dxt1(m)
     .  +dxt(m)*dxt1(l)+dnee*dph(10+jj)
      if(m.eq.3.and.l.eq.3) dphlm=dphlm+dnee*(12*tkev+2*nb*ak0)/zf2
   78 dph(ii)=a03*dphlm
   79 continue
      if(nu.eq.1) go to 90
      nu=1
      nb=40
      a03=a03/rho(1)
      k1=11
      go to 75
c  complete ionization
   80 call zeroc(dne,10)
      call zeroc(dph,20)
      go to 100
   90 continue
      do 95 i=1,iders
      pt(i)=pt(i)+dph(i)
   95 ht(i)=ht(i)+dph(i+10)
c
  100 continue
c
c  *********************************************
c
c  electron pressure and enthalpy
c
      ii=4
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
     .  +anee*hst(4+lm+n)
      pt(jj)=pt(jj)+amm3*pee*(phi(11+l)*phi(11+m)*phi(11+n)
     .  +phi(12+lm)*phi(11+n)+phi(12+l+n)*phi(11+m)
     .  +phi(12+m+n)*phi(11+l)+phi(14+l+m+n))
  105 helm=che*helm
      ht(jj)=ht(jj)+helm
  106 he(jj)=helm
  108 continue
  110 continue
c
c  *********************************************
c
c  ionization enthalpy
c
c  (pi is introduced for consistency)
      call zeroc(pi,10)
      call zeroc(hi(11),10)
c
      xi(1)=xia
      averg=av*ergev
c  combine h and h-
c	Write (6,*)ider
      do 111 l=1,ider
  111 dxt1(l)=exh*xi(l)-exhm*xhm(l)
      do 112 l=1,4
  112 dxt(l)=dxt1(l)/ah-xi(10+l)/ahe
c
      hi1=averg*(dxt1(1)*x/ah+xi(11)*y/ahe+z*ueh(1))
      hi(1)=hi1
      ht(1)=ht(1)+hi1
      ii=4
      do 115 l=2,4
      dhi=dxt1(l)*x/ah+xi(l+10)*y/ahe+z*ueh(l)
      tstl=l.eq.4
      if(tstl) dhi=dhi+dxt(1)
      dhi=averg*dhi
      ht(l)=ht(l)+dhi
      hi(l)=dhi
      if(nosd) go to 115
      do 114 m=l,4
      ii=ii+1
      dhi=dxt1(ii)*x/ah+xi(10+ii)*y/ahe+z*ueh(ii)
      if(tstl) dhi=dhi+dxt(m)
      if(m.eq.4) dhi=dhi+dxt(l)
      dhi=averg*dhi
      ht(ii)=ht(ii)+dhi
  114 hi(ii)=dhi
  115 continue
c
c  *********************************************
c
c  pressure and enthalpy of heavy particles
c
  120 anh=(x/ah+y/ahe+z/az)*av
      rhh=rho(1)
      phh=anh*rhh*tkergs
      ph(1)=phh
      ph(2)=amm*phh*rho(2)
      drht=1+rho(3)
      ph(3)=amm*phh*drht
      ph(4)=amm*phh*rho(4)+rhh*tkergs*avd1
      do 121 i=1,4
  121 pt(i)=pt(i)+ph(i)
      if(nosd) go to 125
      ph(10)=amm*(phh*rho(10)+rho(4)*(ph(4)+rhh*tkergs*avd1))
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
     .  ph(l+1)*rho(k+4)+phh*rho(kl+6))
      if(k.gt.1) pt(9+kl)=pt(9+kl)+amm*(ph(kl)*drht+ph(k1)*
     .  rho(6+l)+ph(l+1)*rho(6+k)+phh*rho(9+kl))
  124 continue
      pt(20)=pt(20)+amm*(ph(10)*rho(4)+2*ph(4)*rho(10)
     .  +phh*rho(20)+rhh*tkergs*avd1*(amm*rho(4)*rho(4)+rho(10)))
c
  125 hhh=2.5*anh*tkergs
      hh(1)=hhh
      hh(2)=0
      hh(3)=amm*hhh
      hh(4)=2.5*tkergs*avd1
      hh(5)=0
      hh(6)=0
      hh(7)=0
      hh(8)=amm2*hhh
      hh(9)=amm*hh(4)
      hh(10)=0
      call zeroc(hh(11),10)
      hh(17)=amm*hh(8)
      hh(18)=amm*hh(9)
      ht(1)=ht(1)+hhh
      ht(3)=ht(3)+amm*hhh
      hh4=2.5*tkergs*avd1
      ht(4)=ht(4)+hh4
      if(nosd) go to 130
      ht(8)=ht(8)+amm2*hhh
      ht(9)=ht(9)+amm*hh4
      if(notd) go to 130
      ht(17)=ht(17)+amm3*hhh
      ht(18)=ht(18)+amm2*hh4
c
c  *********************************************
c
c  pressure and enthalpy of radiation (included if iprrad.ne.0)
c
  130 call zeroc(pr,10)
      call zeroc(hr,20)
      if(iprrad.eq.0) go to 145
      t4=t*t*t*t
      prr=car*t4/3
      pr(1)=prr
      call zeroc(pr(2),9)
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
     .  -hrr*rho(jj))
  138 continue
  140 continue
      do 142 i=1,ider
  142 ht(i)=ht(i)+hr(i)
c     if(notd) go to 145
c..      write(istdpr,14201) pe,ph,pr,(dph(i),i=1,10),pt
14201 format(/' pe:',1p10e12.4/' ph:',1p10e12.4/' pr:',10e12.4/
     *  ' dp:',10e12.4/' pt:',10e12.4/4x,10e12.4)
c
c  pressure and enthalpy from Coulomb effect
c
  145 if(icoulm.ne.1.and.icoulm.ne.-1) then
        do 150 i=1,ider
        pt(i)=pt(i)+pcoul(i)
  150   ht(i)=ht(i)+hcoul(i)
      end if
c
c  *********************************************
c
c  change to derivatives of log p
c
      ptt=pt(1)
c
c  first divide all derivatives by pt, to avoid over- and
c  underflow
c
      do 15010 i=2,ider
15010 pt(i)=pt(i)/ptt
c
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
     .  +2*pt(m+1)*pt(l+1)*pt(n+1)))/amm
  155 ii=4
      do 158 l=2,4
      ptl=pt(l)/amm
      if(nosd) go to 158
      do 156 m=l,4
      ii=ii+1
  156 pt(ii)=(pt(ii)-pt(l)*pt(m))/amm
  158 pt(l)=ptl
c     if(.not.notd) write(istdpr,15801) pt
15801 format(/' pt:',1p10e12.4/4x,10e12.4)
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
      fcpp=1./(amm*t*pf)
      cpp=dxtt*fcpp
      cp(1)=cpp
      if(nosd) go to 173
      do 170 l=2,4
      dcp=dxt(l)*fcpp-cpp*pt(l+3)/pf
      if(l.eq.3) dcp=dcp-cpp*amm
  170 cp(l)=dcp
c
  173 prh=amm*ptt/rhh
      anum=pf*prh-hf
      dad(1)=anum/dxtt
      if(nosd) go to 177
      do 175 l=2,4
      lf=l+3
  175 dad(l)=(prh*(amm*(pt(l)-rho(l))*pf+pt(lf))-ht(lf)-dad(1)*dxt(l)
     .  )/dxtt
c  further derivatives
  177 rhf=rho(2)
      dxtt=pt(3)*rhf-pf*rho(3)
c
c  Gamma1 and derivatives
c
      dxt(1)=rhf-dad(1)*dxtt
      gm1=pf/dxt(1)
      gmm1(1)=gm1
      if(secder) then
        dxt(2)=rho(5)-dad(2)*dxtt
     *    -dad(1)*(pt(6)*rho(2)+pt(3)*rho(5)-pt(5)*rho(3)-pt(2)*rho(6))
        dxt(3)=rho(6)-dad(3)*dxtt
     *    -dad(1)*(pt(8)*rho(2)+pt(3)*rho(6)-pt(6)*rho(3)-pt(2)*rho(8))
        dxt(4)=rho(7)-dad(4)*dxtt
     *    -dad(1)*(pt(9)*rho(2)+pt(3)*rho(7)-pt(7)*rho(3)-pt(2)*rho(9))
        gmm1(2)=(pt(5)-gm1*dxt(2))/dxt(1)
        gmm1(3)=(pt(6)-gm1*dxt(3))/dxt(1)
        gmm1(4)=(pt(7)-gm1*dxt(4))/dxt(1)
      else
        call zeroc(gmm1(2),3)
      end if
c
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
     .  -delta*pt(lf)/pf
  190 xii1(1)=xii(1)
      xii1(2)=xii(11)
      xii1(3)=xii(21)
      xii1(4)=anuh(1)/anh0
      return
 1010 format(/' **** error in eqstf. With H-, Ne .le. 0.'/
     *         '      x(H+) =',1pe11.3,'  x(H-) =',e11.3,
     *         ' x(heavies) =',e11.3/
     *         '      Ignore H-')
 1100 format(/' ***** error in s/r eqstf.',
     *  ' Iteration for Coulomb term failed to converge.'/
     *  7x,' log f, log T, X, Z =',4f10.5//
     *  ' last changes in delta mu_c(1-3):'/ 1p3e13.5)
      end
