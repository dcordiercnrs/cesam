c
c********************************************************************************
c
	subroutine f4der(rhol,tl,sn,nspe,f4,e4,p4,h4,d2f4r2,d2f4t2,
     .           d2f4rt,df4,d2f4i,d2f4f,d2f4x,d2f4r,d2f4t,
     .           dp4i,dp4dr,dp4dt,dp4dx,dh4i,dh4dr,dh4dt,dh4dx,npar)
c
c     nspe  : number of particles actually used
c     npar  : first dimension of d2f4i,d2f4f (=number of ionization
c             degrees = number of number fractions)
c
c     f4    : Coulomb configurational free energy (Debye-Huckel)
c     e4    : internal energy
c     p4    : pressure
c     h4    : enthalpy
c
c     df4   : derivatives with respect to reaction parameters ("Saha equations")
c     d2f4i : derivatives of df4 with respect to ionization degrees
c     d2f4f : derivatives of df4 with respect to number fractions
c     d2f4x : derivatives of df4 with respect to X (Y varying accordingly,
c             and at fixed ionization degrees)
c
c     d2f4r : derivatives of df4 with respect to log10 rho
c     d2f4t : derivatives of df4 with respect to log10 t
c     d2f4r2: derivative  of  f4 with respect to log10 rho ** 2
c     d2f4t2: derivative  of  f4 with respect to log10 t   ** 2
c     d2f4rt: derivative  of  f4 with respect to log10 rho and log10 t
c
c     in both d2f4i and d2f4f the element (i,j) denotes the derivative
c     of df4(i) with respect to parameter j (ionization degree or number
c     fraction)
c
c     dp4i  : derivatives of p4 with respect to ionization degrees
c     dp4dr : derivative  of p4 with respect to log10 rho
c     dp4dt : derivative  of p4 with respect to log10 t
c     dp4dx : derivative  of p4 with respect to X (Y varying accordingly
c             and at fixed ionization degrees)
c
c     dh4i  : derivatives of h4 with respect to ionization degrees
c     dh4dr : derivative  of h4 with respect to log10 rho
c     dh4dt : derivative  of h4 with respect to log10 t
c     dh4dx : derivative  of h4 with respect to X (Y varying accordingly
c             and at fixed ionization degrees)
c
c=========================================================================
c
      implicit real*8 (a-h,o-z)
	implicit integer(i-n)
c
      parameter ( mchem =  3, mspes = 18, mz = 9, mion = mz + 1 )
      parameter ( mlam  = mspes - mchem -1 )
      parameter ( mfe   =   4, mfd    =  5)
c
      common /freecl/
     .                         e     (mfe   ),        p     (mfe   ),
     .                         f     (mfe   ),               d2fdt2 ,
     .                                d2fdtv ,               d2fdv2 ,
     .                         fscr  (mspes ),        dfdn  (mspes ),
     .                         d2fdnt(mspes ),        d2fdnv(mspes ),
     .                         d2fdn2( mspes ,mspes )
      common/constsc/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev1,exhm
c
      dimension sn(nspe),df4(npar),d2f4i(npar,npar),d2f4f(npar,npar),
     .          d2f4x(npar),d2f4r(npar),d2f4t(npar),
     .          dp4i(npar),dh4i(npar)
c
      data umod /2.302585092994046 d0/
c
c...assume ah = 1.008, ahe = 4.0026
      data hovhe,heovh /0.251836d0,3.97083d0/
c
c  reset with JC-D values
c
      hovhe = ah/ahe
      heovh = ahe/ah
c
      t   = 10.d0**tl
      vol = 10.d0**(-rhol)
c
      call f4n(t,vol,sn,nspe)
c
c     f4    : Coulomb configurational free energy (Debye-Huckel)
c     e4    : internal energy
c     p4    : pressure
c     p4    : enthalpy
c
      f4     =  f(4)
      e4     =  e(4)
      p4     =  p(4)
      h4     =  e(4) + p(4)*vol
c
c     d2f4r2: derivative of  f4 with respect to log10 rho ** 2
c     d2f4t2: derivative of  f4 with respect to log10 t   ** 2
c     d2f4rt: derivative of  f4 with respect to log10 rho and log10 t
c
      d2f4r2 =  d2fdv2*vol*vol*umod*umod
      d2f4rt =  d2fdtv*  t*vol*umod*umod
      d2f4t2 =  d2fdt2*  t*  t*umod*umod
c
c     df4   : derivatives with respect to reaction parameters ("Saha equations")
c
      df4(1) =  dfdn(2) + dfdn(6)
      df4(2) =  dfdn(4) + dfdn(6)
      df4(3) = -dfdn(4) + dfdn(5) + dfdn(6)
c
c     d2f4r : derivatives of df4 with respect to log10 rho
c
      d2f4r(1) = -(             d2fdnv(2) + d2fdnv(6))*vol*umod
      d2f4r(2) = -(             d2fdnv(4) + d2fdnv(6))*vol*umod
      d2f4r(3) = -(-d2fdnv(4) + d2fdnv(5) + d2fdnv(6))*vol*umod
c
c     d2f4t : derivatives of df4 with respect to log10 t
c
      d2f4t(1) =  (             d2fdnt(2) + d2fdnt(6))*t  *umod
      d2f4t(2) =  (             d2fdnt(4) + d2fdnt(6))*t  *umod
      d2f4t(3) =  (-d2fdnt(4) + d2fdnt(5) + d2fdnt(6))*t  *umod
c
      toth   =  sn(1) + sn(2)
      tothe  =  sn(3) + sn(4) + sn(5)
c
c     d2f4i : derivatives of df4 with respect to ionization degrees
c
      d2f4i(1,1) = toth *(-d2fdn2(2,1)+d2fdn2(2,2)+d2fdn2(2,6)
     .                    -d2fdn2(6,1)+d2fdn2(6,2)+d2fdn2(6,6))
      d2f4i(1,2) = tothe*(-d2fdn2(2,3)+d2fdn2(2,4)+d2fdn2(2,6)
     .                    -d2fdn2(6,3)+d2fdn2(6,4)+d2fdn2(6,6))
      d2f4i(1,3) = tothe*(-d2fdn2(2,3)+d2fdn2(2,5)+d2fdn2(2,6)*2.d0
     .                    -d2fdn2(6,3)+d2fdn2(6,5)+d2fdn2(6,6)*2.d0)
c
      d2f4i(2,1) = toth *(-d2fdn2(4,1)+d2fdn2(4,2)+d2fdn2(4,6)
     .                    -d2fdn2(6,1)+d2fdn2(6,2)+d2fdn2(6,6))
      d2f4i(2,2) = tothe*(-d2fdn2(4,3)+d2fdn2(4,4)+d2fdn2(4,6)
     .                    -d2fdn2(6,3)+d2fdn2(6,4)+d2fdn2(6,6))
      d2f4i(2,3) = tothe*(-d2fdn2(4,3)+d2fdn2(4,5)+d2fdn2(4,6)*2.d0
     .                    -d2fdn2(6,3)+d2fdn2(6,5)+d2fdn2(6,6)*2.d0)
c
      d2f4i(3,1) = toth *( d2fdn2(4,1)-d2fdn2(4,2)-d2fdn2(4,6)
     .                    -d2fdn2(5,1)+d2fdn2(5,2)+d2fdn2(5,6)
     .                    -d2fdn2(6,1)+d2fdn2(6,2)+d2fdn2(6,6))
      d2f4i(3,2) = tothe*( d2fdn2(4,3)-d2fdn2(4,4)-d2fdn2(4,6)
     .                    -d2fdn2(5,3)+d2fdn2(5,4)+d2fdn2(5,6)
     .                    -d2fdn2(6,3)+d2fdn2(6,4)+d2fdn2(6,6))
      d2f4i(3,3) = tothe*( d2fdn2(4,3)-d2fdn2(4,5)-d2fdn2(4,6)*2.d0
     .                    -d2fdn2(5,3)+d2fdn2(5,5)+d2fdn2(5,6)*2.d0
     .                    -d2fdn2(6,3)+d2fdn2(6,5)+d2fdn2(6,6)*2.d0)
c
      fach       =  sn(1)*sn(1)/toth
c
c  Note: temporary fudge to avoid problems when sn(3) or sn(4)
c  are zero. Need correction from WD           10/5/90
c
      if(sn(3).ne.0) then
        ff2      =  sn(4)/sn(3)
      else
        ff2      =  1
      end if
      if(sn(4).ne.0) then
        ff3      =  sn(5)/sn(4)
      else
        ff3      =  1
      end if
      ff22       =  ff2*ff2
      fac2       =  1.d0 + ff2
      fac3       =  1.d0 + ff3
      fache      =  sn(3)*sn(3)/tothe
      dhed2      = -fac3*fache
      dhed3      = -ff2*fache
      dhepd2     =  fache
      dhepd3     = -ff22*fache
      dhe2d2     =  ff3*fache
      dhe2d3     =  ff2*fac2*fache
      ded2       =  (1.d0 + ff3*2.d0)*fache
      ded3        = ff2*(2.d0 + ff2)*fache
c
c     d2f4f : derivatives of df4 with respect to number fractions
c
      d2f4f(1,1) = fach*(-d2fdn2(2,1)+d2fdn2(2,2)+d2fdn2(2,6)
     .                   -d2fdn2(6,1)+d2fdn2(6,2)+d2fdn2(6,6))
      d2f4f(2,1) = fach*(-d2fdn2(4,1)+d2fdn2(4,2)+d2fdn2(4,6)
     .                   -d2fdn2(6,1)+d2fdn2(6,2)+d2fdn2(6,6))
      d2f4f(3,1) = fach*( d2fdn2(4,1)-d2fdn2(4,2)-d2fdn2(4,6)
     .                   -d2fdn2(5,1)+d2fdn2(5,2)+d2fdn2(5,6)
     .                   -d2fdn2(6,1)+d2fdn2(6,2)+d2fdn2(6,6))
c
      d2f4f(1,2) =    ( d2fdn2(2,3)+d2fdn2(6,3))*dhed2
     .              + ( d2fdn2(2,4)+d2fdn2(6,4))*dhepd2
     .              + ( d2fdn2(2,5)+d2fdn2(6,5))*dhe2d2
     .              + ( d2fdn2(2,6)+d2fdn2(6,6))*ded2
      d2f4f(2,2) =    ( d2fdn2(4,3)+d2fdn2(6,3))*dhed2
     .              + ( d2fdn2(4,4)+d2fdn2(6,4))*dhepd2
     .              + ( d2fdn2(4,5)+d2fdn2(6,5))*dhe2d2
     .              + ( d2fdn2(4,6)+d2fdn2(6,6))*ded2
      d2f4f(3,2) =    (-d2fdn2(4,3)+d2fdn2(5,3)+d2fdn2(6,3))*dhed2
     .              + (-d2fdn2(4,4)+d2fdn2(5,4)+d2fdn2(6,4))*dhepd2
     .              + (-d2fdn2(4,5)+d2fdn2(5,5)+d2fdn2(6,5))*dhe2d2
     .              + (-d2fdn2(4,6)+d2fdn2(5,6)+d2fdn2(6,6))*ded2
c
      d2f4f(1,3) =    ( d2fdn2(2,3)+d2fdn2(6,3))*dhed3
     .              + ( d2fdn2(2,4)+d2fdn2(6,4))*dhepd3
     .              + ( d2fdn2(2,5)+d2fdn2(6,5))*dhe2d3
     .              + ( d2fdn2(2,6)+d2fdn2(6,6))*ded3
      d2f4f(2,3) =    ( d2fdn2(4,3)+d2fdn2(6,3))*dhed3
     .              + ( d2fdn2(4,4)+d2fdn2(6,4))*dhepd3
     .              + ( d2fdn2(4,5)+d2fdn2(6,5))*dhe2d3
     .              + ( d2fdn2(4,6)+d2fdn2(6,6))*ded3
      d2f4f(3,3) =    (-d2fdn2(4,3)+d2fdn2(5,3)+d2fdn2(6,3))*dhed3
     .              + (-d2fdn2(4,4)+d2fdn2(5,4)+d2fdn2(6,4))*dhepd3
     .              + (-d2fdn2(4,5)+d2fdn2(5,5)+d2fdn2(6,5))*dhe2d3
     .              + (-d2fdn2(4,6)+d2fdn2(5,6)+d2fdn2(6,6))*ded3
c
c     d2f4x : derivatives of df4 with respect to X (Y varying accordingly)
c
      dnhdx    =  (toth + tothe*heovh)/toth
      dnhedx   = -(toth*hovhe + tothe)/tothe
c
      d2f4x(1) = dnhdx *((d2fdn2(2,1)+d2fdn2(6,1))* sn(1) +
     .                   (d2fdn2(2,2)+d2fdn2(6,2))* sn(2) +
     .                   (d2fdn2(2,6)+d2fdn2(6,6))* sn(2))+
     .           dnhedx*((d2fdn2(2,3)+d2fdn2(6,3))* sn(3) +
     .                   (d2fdn2(2,4)+d2fdn2(6,4))* sn(4) +
     .                   (d2fdn2(2,5)+d2fdn2(6,5))* sn(5) +
     .                   (d2fdn2(2,6)+d2fdn2(6,6))*(sn(4)+2.d0*sn(5)))
c
      d2f4x(2) = dnhdx *((d2fdn2(4,1)+d2fdn2(6,1))* sn(1) +
     .                   (d2fdn2(4,2)+d2fdn2(6,2))* sn(2) +
     .                   (d2fdn2(4,6)+d2fdn2(6,6))* sn(2))+
     .           dnhedx*((d2fdn2(4,3)+d2fdn2(6,3))* sn(3) +
     .                   (d2fdn2(4,4)+d2fdn2(6,4))* sn(4) +
     .                   (d2fdn2(4,5)+d2fdn2(6,5))* sn(5) +
     .                   (d2fdn2(4,6)+d2fdn2(6,6))*(sn(4)+2.d0*sn(5)))
c
      d2f4x(3) = dnhdx *((-d2fdn2(4,1)+d2fdn2(5,1)+d2fdn2(6,1))*
     .                     sn(1) +
     .                   (-d2fdn2(4,2)+d2fdn2(5,2)+d2fdn2(6,2))*
     .                     sn(2) +
     .                   (-d2fdn2(4,6)+d2fdn2(5,6)+d2fdn2(6,6))*
     .                     sn(2))+
     .           dnhedx*((-d2fdn2(4,3)+d2fdn2(5,3)+d2fdn2(6,3))*
     .                     sn(3) +
     .                   (-d2fdn2(4,4)+d2fdn2(5,4)+d2fdn2(6,4))*
     .                     sn(4) +
     .                   (-d2fdn2(4,5)+d2fdn2(5,5)+d2fdn2(6,5))*
     .                     sn(5) +
     .                   (-d2fdn2(4,6)+d2fdn2(5,6)+d2fdn2(6,6))*
     .                    (sn(4)+2.d0*sn(5)))
c
c     dp4i  : derivatives of p4 with respect to ionization degrees
c     dp4dr : derivative  of p4 with respect to log10 rho
c     dp4dt : derivative  of p4 with respect to log10 t
c     dp4dx : derivative  of p4 with respect to X (Y varying accordingly
c             and at fixed ionization degrees)
c
      dp4i(1) = - toth *( -d2fdnv(1)+d2fdnv(2)+d2fdnv(6))
      dp4i(2) = - tothe*( -d2fdnv(3)+d2fdnv(4)+d2fdnv(6))
      dp4i(3) = - tothe*( -d2fdnv(3)+d2fdnv(5)+d2fdnv(6)*2.d0)
c
      dp4dr   =    umod  * vol * d2fdv2
      dp4dt   =  - umod  *   t * d2fdtv
c
      dp4dx   = -dnhdx * (d2fdnv(1)* sn(1) +
     .                    d2fdnv(2)* sn(2) +
     .                    d2fdnv(6)* sn(2))
     .          -dnhedx* (d2fdnv(3)* sn(3) +
     .                    d2fdnv(4)* sn(4) +
     .                    d2fdnv(5)* sn(5) +
     .                    d2fdnv(6)*(sn(4)+2.d0*sn(5)))
c
c     ... analogous for h4 ...
c
      dh4i(1) =   toth *((-dfdn  (1)+dfdn  (2)+dfdn  (6))
     .                  -(-d2fdnt(1)+d2fdnt(2)+d2fdnt(6))*t
     .                  -(-d2fdnv(1)+d2fdnv(2)+d2fdnv(6))*vol  )
c
      dh4i(2) =   tothe*((-dfdn  (3)+dfdn  (4)+dfdn  (6))
     .                  -(-d2fdnt(3)+d2fdnt(4)+d2fdnt(6))*t
     .                  -(-d2fdnv(3)+d2fdnv(4)+d2fdnv(6))*vol  )
c
      dh4i(3) =   tothe*((-dfdn  (3)+dfdn  (5)+dfdn  (6)*2.d0)
     .                  -(-d2fdnt(3)+d2fdnt(5)+d2fdnt(6)*2.d0)*t
     .                  -(-d2fdnv(3)+d2fdnv(5)+d2fdnv(6)*2.d0)*vol  )
c
      dh4dr   =    umod * vol * t * d2fdtv  +  vol * dp4dr
      dh4dt   =  - umod *  t  * t * d2fdt2  +  vol * dp4dt
c
      dh4dx   = dnhdx *((dfdn(1)-t*d2fdnt(1))* sn(1) +
     .                  (dfdn(2)-t*d2fdnt(2))* sn(2) +
     .                  (dfdn(6)-t*d2fdnt(6))* sn(2))+
     .          dnhedx*((dfdn(3)-t*d2fdnt(3))* sn(3) +
     .                  (dfdn(4)-t*d2fdnt(4))* sn(4) +
     .                  (dfdn(5)-t*d2fdnt(5))* sn(5) +
     .                  (dfdn(6)-t*d2fdnt(6))*(sn(4) + 2.d0*sn(5)))+
     .          vol * dp4dx
c
      return
      end
