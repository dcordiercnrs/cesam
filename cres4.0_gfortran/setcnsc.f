 
c*****************************************************************************
c
	subroutine setcnsc
c
c  sets the physical and mathematical constants used in the programme
c
c  numerical constants from CODATA report, with additional values
c  for atomic masses and ionization potentials. this is
c  the GONG set of values.
c
c  corresponds to iver = 2
c
c  modified on 8/12/87 to approach a consistent setting based on
c  minimum set of basic constants. so far keep original values,
c  but print consistent values
c
c  modified 15/3/88 to use consistent numerical constants
c  note that a version number iver has been added in common/fconst/
c  to indicate which set of values have been used.
c
c  also include constants for simple equation of state, with
c  no partial degeneracy, for use for GONG models.
c
c  This version also calls WD routine to prepare for inclusion
c  of Coulomb terms.
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 9/3/90
c
      implicit double precision(a-h,o-z)
	implicit integer(i-n)
      common/constsc/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev1,exhm
c
c  common giving fundamental constants
c
      common /fconst/ pi, ebase, amu, ame, clight, planck, boltzm,
     *  cgrav, amsun, echar, ergev, syear, iver
      common/eqphsm/ rho00, f0
      common/ln10c/ amm,amm2,amm3
      common/eqstdc/ ccc1(94)
      common/eqsoutc/ ccc2(230)
	common/dmudec/ ccc3(10),idmu
      common/cversn/ ivrbcs,ivrrhs,ivreos,ivropc,ivreng
c
      save
c
c     write(6,100)
c
c  version number for constants
c
c
c  *****************************************************************
c
c  set equation of state version number to 1 for EFF with Coulomb
c  terms (note that actual value may depend on icoulm, and hence
c  may be reset later, after call of eqstfc
c
      ivreos = 1
c
c  *****************************************************************
c
      iver = 2
c
      pi=4.d0*atan(1.d0)
      ebase=exp(1.d0)
c
c  atomic mass unit
      amu=1.6605402d-24
c  electron mass
      ame=9.1093897d-28
c  speed of light
      clight=2.99792458d10
c  planck's constant
      planck=6.6260755d-27
c  boltzman's constant, cgs
      boltzm=1.380658d-16
c  gravitational constant (chosen from planetary orbit data, assuming
c  msun = 1.989e33. consistent, but not identical, with CODATA value).
      cgrav=6.67232d-8
c  solar mass
      amsun=1.989d33
c  electron charge, in ESU
      echar=4.8032068d-10
c  number of ergs in 1 ev
      ergev=1.60217733d-12
c  number of seconds in a year
c  (should be replaced with second value)
      syear =3.155692597d7
c
c  start on derived values, largely
c
c     write(6,110)
c  avogadro's number (should be 1/amu)
      avp=6.02217d23
      av=1.0/amu
c     write(6,120) 'av',avp,av
c  atomic weights of h and he
      ah=1.007825d0
      ahe=4.002603d0
c  average atomic weight for heavy elements.
c  before 4/1/84 was given value 17.8, but was in fact
c  most often reset to the current value in calling programme.
      az=16.389d0
c  auxiliary quantities for eggleton fudge calculation
      avda=av*(1/ah-2/ahe)
      avd1=av*(1/ah-1/ahe)
c  boltzmann's constant in ev/deg
      ck1p=8.6170837d-5
      ck1=boltzm/ergev
c     write(6,120) 'ck1',ck1p,ck1
c  the same, in ergs/deg
      ck2=boltzm
c  ionization potentials
      exh=13.595d0
      exhe=24.580d0
      exhep=54.403d0
c  constants for transition to starred variables
      ctp=1.686304d-10
      crhop=1.759547d30
      cpep =1.440588d24
      chep =8.187265d-7
c
      alamc=planck/(ame*clight)
      alamc3=alamc**3
      ct=boltzm/(ame*clight**2)
      crho=8*pi/alamc3
      che=ame*clight**2
      cpe=crho*che
c     write(6,120) 'ct',ctp,ct
c     write(6,120) 'crho',crhop,crho
c     write(6,120) 'cpe',cpep,cpe
c     write(6,120) 'che',chep,che
c  constants for pressure ionization
      ca03p=2.147d-24
      caap=1.759547d30*ca03p
      ckhp=13.5d0
c
      a0=planck/(2*pi*echar)
      a0=(a0/ame)*a0
c..      write(6,*) 'a0 =',a0
c
c  for the moment, set additional fudge factor by hand at this
c  point. need to think about a more suitable way to do it later
c
      efffac=15
      ca03=efffac*a0**3
      caa=crho*ca03
      ckh=exh
c
c  change from ev to ergs and include factor av**2 to compensate for
c  redefinition of dne (see notes of 3/1/84)
c
      ca03p=av*1.602192d-12*ca03p*av
      ca03=av*ergev*ca03*av
c
c     write(6,120) 'ca03',ca03p,ca03
c     write(6,120) 'caa',caap,caa
c     write(6,120) 'ckh',ckhp,ckh
c
c  the radiation constant
c
      carp=7.5647d-15
      car=boltzm/(clight*planck)
c..      write(6,*) 'boltzm/(clight*planck)', carc
      car=8*pi**5*boltzm*car**3/15
c     write(6,120) 'car',carp,car
c  number of ergs in 1 ev
      ergev1=ergev
c  ionization potential for h-
      exhm=0.754d0
c  ln 10
      amm=log(1.d1)
      amm2=amm*amm
      amm3=amm2*amm
c
c constant for phderc, with no partial degeneracy
c
      rho00=sqrt(2*pi)*ebase**2/8
      f0=av/(rho00*crho*ct**1.5d0)
c..     write(6,*) 'f0',f0
c
c  set commons from s/r eqstfc to zero
      call zeroc(ccc1,90)
      call zeroc(ccc2,210)
      call zeroc(ccc3,10)
c
c====================== initialize mhd constants ==================
      nspe=6
      call setf4(nspe)
c=====================================================================
c
      return
  100 format(//2x,75('*')//
     *  ' Set constants from CODATA Report.'/' See Cohen & Taylor,',
     *  ' Rev. Mod. Phys., vol. 59, 1121 (1987)'/
     *  ' ++++ Double precision version'/)
  110 format(/' variable   old value    current, consistent value'/)
  120 format(a10,1p2e18.7)
      end
