c
c***************************************************************************
c
	subroutine setf4(nspe)
c
c  Note: commons free, map and spes have been renamed freecl,
c  mapcl and spescl to avoid conflict when calling routine
c  from MHD table emulator.
c                                           1/5/90       jcd
c
c  Modified 11/5/90 to allow arbitrary order of expansion for tau
c  in s/r f4mhd (modification necessary because of unavailability
c  of real*16).                                          jcd
c
c  Modified 5/6/90 to reset constants with JC-D values for consistency.
c  Note: this requires that s/r setcnsc be called before calling setf4.
c  In JC-D usage, setf4 would normally be called from setcnsc.
c
c======================================================================
c
      implicit real*8 (a-h,o-z)
	implicit integer(i-n)
c
      parameter ( mchem =  3, mspes = 18, mz = 9, mion = mz + 1 )
      parameter ( mlam  = mspes - mchem -1 )
      parameter ( mfe   =   4, mfd    =  5)
c
      common /cstmhd/
     .                                camu   ,               cc     ,
     .                                ce     ,               ch     ,
     .                                ck     ,               cme    ,
     .                                cpi    ,               cevw   ,
     .                                cx     ,               cf4    ,
     .                                carad
c
      common /degmhd/
     .                                eta    ,               detdn  ,
     .                                detdt  ,               detdv  ,
     .                                d2etdn2,               d2etdnt,
     .                                d2etdnv,               d2etdt2,
     .                                d2etdtv,               d2etdv2,
     .                                exeta  ,               thet   ,
     .                                dthet  ,               d2thet ,
     .                                fd    (mfd   )
      common /freecl/
     .                         e     (mfe   ),        p     (mfe   ),
     .                         f     (mfe   ),               d2fdt2 ,
     .                                d2fdtv ,               d2fdv2 ,
     .                         fscr  (mspes ),        dfdn  (mspes ),
     .                         d2fdnt(mspes ),        d2fdnv(mspes ),
     .                         d2fdn2( mspes ,mspes )
      common /mapcl /
     .                                nchem  ,               nspes  ,
     .                                ise
      common /spescl/
     .                         zs    (mspes ),        zsq   (mspes ),
     .                         zmask (mspes ),        sn    (mspes )
c
c  commons of fundamental constants set by JC-D routine setcnsc.
c  for consistency, replace WD values by these.
c
      common/constsc/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev1,exhm
c
c  common giving fundamental constants
c
      common /fconst/ pi, ebase, amu, ame, clight, planck, boltzm,
     *  cgrav, amsun, echar, ergev, syear, iver
c
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
c
      if(idgeos.ge.1) write(6,*) 'Calling setf4'
c
      camu  =  1.6605655d-24
      cc    =  2.9979246d+10
      ce    =  4.8032426d-10
      ch    =  6.6261764d-27
c
      ck    =  1.3806624d-16
      cme   =  9.1095345d-28
      cpi   =  3.141592654d0
      cevw  =  8.0654653d+03
c
c  resetting with JC-D values
c
      camu  =  amu
      cc    =  clight
      ce    =  echar
      ch    =  planck
      ck    =  boltzm
      cme   =  ame
      cpi   =  pi
c
c  Note: the meaning of cevw is currently unclear. On the other
c  hand, it appears not to be used. Should be fixed up.
 
c
      cf4   = 2.0d0 * dsqrt( cpi ) * ce**3 / ( 3.0d0 * dsqrt( ck ) )
      cx    = 3.0d0 * cf4 / ck
c
      carad =  7.56567  d-15
c
c  JC-D value
c
      carad =  car
c
      nspes =  nspe
      ise   =  nspe
c
      zmask(1) = 0.d0
      zmask(2) = 1.d0
      zmask(3) = 0.d0
      zmask(4) = 1.d0
      zmask(5) = 1.d0
c
      zs(1) = 0.d0
      zs(2) = 1.d0
      zs(3) = 0.d0
      zs(4) = 1.d0
      zs(5) = 2.d0
c
      zsq(1) = 0.d0
      zsq(2) = 1.d0
      zsq(3) = 0.d0
      zsq(4) = 1.d0
      zsq(5) = 4.d0
c
c................... following values should never be used; if none the less,
c................... they ought to provoke a crash of the program!
      zmask(6) = -1.d200
      zs(6)    = -1.d200
      zsq(6)   = -1.d200
c
      return
      end
