 
c***********************************************************************
 
	subroutine setcns
 
c	routine pour EFF, Auteur J. Christensen Dalsgaard, adaptation:
c	A. Baglin, M. Auvergne, P. Morel
 
c	Les valeurs des constantes ont ete modifiees pour rendre
c	leur valeurs coherentes avec le ss-pg ctes du code de Pierre.
c	Michel Juillet 90.
c	les constantes modifies sont signalees par !*
c  sets the physical and mathematical constants used in the program
 
	implicit real*8 (a-h,o-z)
 
c      implicit real*16 (a-h,o-z)
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .   ct,crho,cpe,che,ca03,caa,ckh,car,ergev
      common/ln10/ amm,amm2,amm3
      common/eqstd/ ccc1(90)
      common/eqsout/ ccc2(210)
      common/dmuder/ ccc3(11)
c
c
c  avogadro's number
      av=6.0221366d23		!*
c  atomic weights of h and he
      ah=1.007825d0		!*
      ahe=4.002603d0		!*
      az=17.8d0
      avda=av*(1/ah-2/ahe)
      avd1=av*(1/ah-1/ahe)
c  boltzmann's constant in ev/deg
      ck1=8.6170837d-5
c  the same, in ergs/deg
      ck2=1.380658d-16	!*
c  ionization potentials
      exh=13.595d0
      exhe=24.580d0
      exhep=54.403d0
c  constants for transition to starred variables
      ct=1.686304d-10
      crho=1.759547d30
      cpe =1.440588d24
      che =8.187265d-7
c  constants for pressure ionization
      ca03=2.147d-24
      caa=1.759547d30*ca03
      ckh=13.5d0
c  change from ev to ergs
      ca03=1.60217733d-12*ca03	!*
c  the radiation constant
      car=7.5659122d-15		!*
c  number of ergs in 1 ev
      ergev=1.60217733d-12		!*
c  ln 10
      amm=log(1.d1)
      amm2=amm*amm
      amm3=amm2*amm
c  set commons from s/r eqstf to zero
      call zero(ccc1,90)
      call zero(ccc2,210)
      call zero(ccc3,10)
      return
 
	end
