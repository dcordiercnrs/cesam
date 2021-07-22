c
c*****************************************************************************
c
	subroutine f4mhd(t,vol)
c     ................................
c     (modified from MHD package)
c     ................................
      implicit real*8 (a-h,o-z)
	implicit integer(i-n)
c
c     free energy of coulomb interactions and derivatives              c
c
      parameter ( mchem =  3, mspes = 18, mz = 9, mion = mz + 1 )
      parameter ( mlam  = mspes - mchem -1 )
      parameter ( mfe   =   4, mfd    =  5)
c
c  parameters for tau expansion. xtautr gives transition between
c  expansion and direct expression.
c  With nmax = 15, xtautr = 0.1, the error at the transition point
c  is around 1.e-13 in tau, 1.e-10 in dtau and 1.e-9 in d2tau.
c
      parameter ( nmax = 15, xtautr = 0.1)
c
c-------------------------------------------------------------------
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
c
c  for running on Suns, take out quadruple precision.
c
c  *** Note: we may need to fix this up in expansion later, by
c      including more terms
c
c..      real*16 dpx, dp1, dp2, dp3
c
      dimension   dxdn(mspes), d2xdnt(mspes), d2xdnv(mspes),
     .            d2xdn2(mspes, mspes)
      dimension   ctau(nmax), cdtau(nmax), cd2tau(nmax)
c
      equivalence (dxdn  , dzdn  ), (d2xdnt, d2zdnt), (d2xdnv, d2zdnv),
     .            (d2xdn2, d2zdn2)
c
      data initcf /0/
c
      save initcf, ctau, cdtau, cd2tau
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     sums over charges                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      sumn = dmax1( sdot(nspes - 1, zmask, 1, sn, 1), 1.d-70 )
      zn   = dmax1( sdot(nspes - 1, zs   , 1, sn, 1), 1.d-70 )
      znt  = dmax1( sdot(nspes - 1, zsq  , 1, sn, 1) + sn(ise)*thet,
     .                                              1.d-70 )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     parameter x and its derivatives                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      frat = fd(2) / fd(3)
      x    = cx * frat * (zn/sumn) * dsqrt( znt ) / dsqrt( vol*t**3 )
c
c-----------------------------------------------------------------------
c     zero everything
c-----------------------------------------------------------------------
c
      do 2 is = 1, nspes
      fscr  (is) = 0.0d0
      dxdn  (is) = 0.0d0
      d2xdnt(is) = 0.0d0
      d2xdnv(is) = 0.0d0
      do 1 js = 1, nspes
      d2xdn2(is, js) = 0.0d0
    1 continue
    2 continue
c
c-----------------------------------------------------------------------
c     first derivatives                                                c
c-----------------------------------------------------------------------
c
      do 3 is = 1, nspes - 1
      dxdn(is) = x * ( zs(is)/zn - zmask(is)/sumn+0.5d0*zsq(is)/znt )
    3 continue
c
      dxdn(ise) = x * ( 0.5d0* (thet + sn(ise)*dthet*detdn)/znt
     .                - 1.5d0* frat * detdn
     .                + 1.0d0 / sn(ise) )
c
      dxdt = x * ( detdt*(thet - 1.5d0*frat + 0.5d0*sn(ise)*dthet/znt)
     .           - 1.5d0/t   )
      dxdv = x * ( detdv*(thet - 1.5d0*frat + 0.5d0*sn(ise)*dthet/znt)
     .           - 0.5d0/vol )
c
c-----------------------------------------------------------------------
c     second derivatives
c-----------------------------------------------------------------------
c
      do 5 js = 1, nspes - 1
c
      if ( zs (js) .ne. 0.0d0 ) then
            do 4 is = 1, js
            d2xdn2(is, js) = ( dxdn(is) * dxdn(js) / x
     .                   + x * (       zmask(is) * zmask(js) / sumn**2
     .                         -       zs   (is) * zs   (js) / zn**2
     .                       - 0.5d0 * zsq  (is) * zsq  (js) / znt**2) )
    4       continue
c
            if( js .gt. 1 ) call scopy( js - 1, d2xdn2( 1,js), 1,
     .                                          d2xdn2(js, 1), mspes )
      d2xdnt(js) = dxdt * dxdn(js) / x
     .               - 0.5d0 * x * zsq(js)*sn(ise) * dthet*detdt/znt**2
      d2xdnv(js) = dxdv * dxdn(js) / x
     .               - 0.5d0 * x * zsq(js)*sn(ise) * dthet*detdv/znt**2
      end if
    5 continue
c
      do 6 is = 1, nspes - 1
      d2xdn2(is, ise) = ( dxdn(is) * dxdn(ise) / x
     .   - 0.5d0 * x * zsq(is) * (thet + sn(ise)*dthet*detdn)/znt**2 )
    6 continue
      call scopy( nspes - 1, d2xdn2(1,ise), 1, d2xdn2(ise,1), mspes )
c
      d2xdn2(ise, ise) =
     .    dxdn(ise)**2/x
     .    - x * ( 1.0d0/sn(ise)**2
     .       + 1.5d0 * frat * (d2etdn2 + (thet - 1.5d0*frat)*detdn**2)
     .         - 0.5d0 * ( d2thet * sn(ise) * detdn**2
     .                  + dthet * (2.d0*detdn + sn(ise)*d2etdn2)
     .                - (thet + sn(ise)*dthet*detdn)**2/znt ) / znt  )
c
      d2xdnt(ise) = dxdt*dxdn(ise)/x + x * (
     .  -1.5d0*frat* (d2etdnt + (thet - 1.5d0*frat)*detdn*detdt)
     .  +0.5d0*( d2thet * sn(ise)*detdn*detdt
     .       + dthet * (detdt + sn(ise)*d2etdnt)
     .     - (thet+sn(ise)*dthet*detdn)*sn(ise)*dthet*detdt/znt)/znt )
c
      d2xdnv(ise) = dxdv*dxdn(ise)/x + x * (
     .  -1.5d0*frat* (d2etdnv + (thet - 1.5d0*frat)*detdn*detdv)
     .  +0.5d0*( d2thet * sn(ise)*detdn*detdv
     .       + dthet * (detdv + sn(ise)*d2etdnv)
     .     - (thet+sn(ise)*dthet*detdn)*sn(ise)*dthet*detdv/znt)/znt )
c
      d2xdt2 = dxdt**2/x + x * ( 1.5d0/t**2
     .     + d2etdt2 * ( thet - 1.5d0*frat + 0.5d0*sn(ise)*dthet/znt )
     .       + detdt**2 * ( dthet - 1.5d0*frat*(thet - 1.5d0 * frat)
     .         + 0.5d0*sn(ise)*(d2thet - sn(ise)*dthet**2/znt)/znt ) )
c
      d2xdtv = dxdt*dxdv/x + x * (
     .       d2etdtv * ( thet - 1.5d0*frat + 0.5d0*sn(ise)*dthet/znt )
     .    + detdt * detdv * ( dthet - 1.5d0*frat*(thet - 1.5d0 * frat)
     .         + 0.5d0*sn(ise)*(d2thet - sn(ise)*dthet**2/znt)/znt ) )
c
      d2xdv2 = dxdv**2/x + x * ( 0.5d0/vol**2
     .     + d2etdv2 * ( thet - 1.5d0*frat + 0.5d0*sn(ise)*dthet/znt )
     .     + detdv**2 * ( dthet - 1.5d0*frat*(thet - 1.5d0 * frat)
     .         + 0.5d0*sn(ise)*(d2thet - sn(ise)*dthet**2/znt)/znt ) )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     tau and its derivatives                                          c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c..      if ( x .gt. 1.0d-1 ) then
c
c
      if ( x .gt. xtautr ) then
          dpx =   x
          dp1 =   3.d0 * ( log(1.d0 + dpx)
     .                   - dpx*(1.d0 - 0.5d0*dpx) ) / dpx**3
          dp2 = - 3.d0 * ( dp1 - 1.d0/(1.d0 + dpx)    ) / dpx
          dp3 =   - ( 4.d0*dp2 + 3.d0/(1.d0 + dpx)**2 ) / dpx
          tau   = dp1
          dtau  = dp2
          d2tau = dp3
      else
c..          tau  = ((((((((x/11. - 1.d0/10.)*x +  1.d0/9.)*x - 1.d0/8.)*x
c..     .                         + 1.d0/7. )*x -  1.d0/6.)*x + 1.d0/5.)*x
c..     .                         - 1.d0/4. )*x +  1.d0/3.)*3.d0
c..          dtau = (((((((8.d0*x/11. - 7.d0/10.)*x +2.d0/3.)*x-5.d0/8.)*x
c..     .                              + 4.d0/7.)*x -1.d0/2.)*x+2.d0/5.)*x
c..     .                              - 1.d0/4. )*3.d0
c..          d2tau= ((((((56.d0*x/11. - 21.d0/5.)*x+10.d0/3.)*x-5.d0/2.)*x
c..     .                             + 12.d0/7.)*x-1.d0)*x+2.d0/5.)*3.d0
c
c  test for setting coefficients
c
        if(initcf.eq.0) then
c
c  set coefficients for expansion
c
          isg=6*mod(nmax,2)-3
          do n=nmax,3,-1
            ctau(n)  =dfloat(isg)/dfloat(n)
            cdtau(n) =dfloat(isg*(n-3))/dfloat(n)
            cd2tau(n)=dfloat(isg*(n-3)*(n-4))/dfloat(n)
            isg=-isg
          end do
          initcf = 1
        end if
c
c  do expansion as do loop, to allow arbitrary high order
c
        tau=ctau(nmax)
        do n=nmax-1,3,-1
          tau=tau*x+ctau(n)
        end do
c
        dtau=cdtau(nmax)
        do n=nmax-1,4,-1
          dtau=dtau*x+cdtau(n)
        end do
c
        d2tau=cd2tau(nmax)
        do n=nmax-1,5,-1
          d2tau=d2tau*x+cd2tau(n)
        end do
      end if
c
c..      write(6,*) 'x, tau, etc =',x,tau, dtau, d2tau
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     free energy, pressure, internal energy                           c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c..      write(6,*) 'x, tau, dtau',x, tau,dtau
c     free energy
      f(4) = - cf4 * tau * dsqrt( znt**3 / (t * vol) )
c
c     pressure
      p(4) = f(4) * ( 0.5d0/vol-1.5d0* sn(ise) * dthet * detdv / znt
     .                          - dxdv * dtau / tau )
c
c     internal energy
      e(4) = f(4) * ( 1.5d0/t -1.5d0* sn(ise) * dthet * detdt / znt
     .                          - dxdt * dtau / tau ) * t
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     derivatives of f4                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c-----------------------------------------------------------------------
c     first derivatives                                                c
c-----------------------------------------------------------------------
c
      df4dt =   f(4) * ( -0.5d0/ t + 1.5d0*sn(ise) * dthet * detdt / znt
     .                            + dxdt * dtau / tau )
      df4dv = - p(4)
c
      do 7 is = 1, nspes - 1
      fscr(is) = f(4)*(1.5d0*zsq(is)/znt + dxdn(is)*dtau/tau)
      dfdn(is) = fscr(is)
    7 continue
      fscr(ise) = f(4)*( 1.5d0* ( thet + sn(ise)*dthet*detdn )/znt
     .                   + dxdn(ise) * dtau/tau )
      dfdn(ise) = fscr(ise)
c
c-----------------------------------------------------------------------
c     second derivatives                                               c
c-----------------------------------------------------------------------
c
      do 9 js = 1, nspes - 1
c
      if ( zs(js)  .ne. 0.d0 ) then
            do 8 is = 1, js
            d2fdn2(is, js) = ( fscr(is) * fscr(js) / f(4)
     .            + ( (d2tau/tau - (dtau/tau)**2) * dxdn(is) * dxdn(js)
     .            + d2xdn2(is, js) * dtau / tau
     .            - 1.5d0* zsq(is) * zsq(js) / znt**2 ) * f(4) )
    8       continue
            if( js .gt. 1 ) call scopy( js - 1, d2fdn2(1 ,js), 1,
     .                                          d2fdn2(js, 1), mspes )
c
            d2fdnt(js) = fscr(js) * df4dt / f(4)
     .              + ( (d2tau/tau - (dtau/tau)**2) * dxdn(js) * dxdt
     .              + d2xdnt(js) * dtau/tau
     .              - 1.5d0*zsq(js)*sn(ise)*dthet*detdt/znt**2)*f(4)
c
            d2fdnv(js) = fscr(js) * df4dv / f(4)
     .              + ( (d2tau/tau - (dtau/tau)**2) * dxdn(js) * dxdv
     .              + d2xdnv(js) * dtau/tau
     .              - 1.5d0*zsq(js)*sn(ise)*dthet*detdv/znt**2)*f(4)
      end if
    9 continue
c
      do 10 is = 1, nspes - 1
      d2fdn2(is, ise) = ( fscr(is) * fscr(ise) / f(4)
     .    + f(4)*( (d2tau/tau - (dtau/tau)**2) * dxdn(is) * dxdn(ise)
     .    + d2xdn2(is, ise) * dtau / tau
     .    -1.5d0*zsq(is)*(thet+sn(ise)*dthet*detdn)/znt**2))
   10 continue
      call scopy( nspes - 1, d2fdn2(1,ise), 1, d2fdn2(ise,1), mspes )
c
      d2fdn2(ise, ise) = fscr(ise)**2 / f(4)
     .    + f(4) * ( (d2tau/tau - (dtau/tau)**2) * dxdn(ise)**2
     .    + d2xdn2(ise, ise) * dtau / tau
     .    + 1.5d0*( d2thet * sn(ise) * detdn**2
     .    + dthet * ( 2.d0*detdn + sn(ise) * d2etdn2 )
     .    - (thet + sn(ise)*dthet*detdn)**2/znt) / znt )
c
      d2fdnt(ise) = fscr(ise)*df4dt/f(4) + f(4) *
     . ( (d2tau/tau - (dtau/tau)**2)*dxdt*dxdn(ise)
     . + d2xdnt(ise) * dtau/tau
     . + 1.5d0 *( d2thet*sn(ise)*detdt*detdn
     . +dthet*(detdt + sn(ise)*d2etdnt)
     . -(thet+sn(ise)*dthet*detdn)*sn(ise)*dthet*detdt/znt)/znt )
c
      d2fdnv(ise) = fscr(ise)*df4dv/f(4) + f(4) *
     . ( (d2tau/tau - (dtau/tau)**2)*dxdv*dxdn(ise)
     . + d2xdnv(ise) * dtau/tau
     . + 1.5d0 *( d2thet*sn(ise)*detdv*detdn
     . +dthet*(detdv + sn(ise)*d2etdnv)
     . -(thet+sn(ise)*dthet*detdn)*sn(ise)*dthet*detdv/znt)/znt )
c
      d2fdtv = df4dt*df4dv/f(4) + f(4) *
     . ( (d2tau/tau - (dtau/tau)**2)*dxdt*dxdv
     . + d2xdtv * dtau/tau
     . + 1.5d0 * sn(ise) * ( d2thet*detdt*detdv
     . + dthet*(d2etdtv-sn(ise)*dthet*detdt*detdv/znt))/znt )
c
      d2fdt2 = df4dt**2/f(4) + f(4) *
     . ( 0.5d0 / t**2
     . + (d2tau/tau - (dtau/tau)**2)*dxdt**2
     . + d2xdt2 * dtau/tau
     . + 1.5d0 * sn(ise) * ( d2thet*detdt**2
     . + dthet*(d2etdt2-sn(ise)*dthet*detdt**2/znt))/znt )
c
      d2fdv2 = df4dv**2/f(4) + f(4) *
     . ( 0.5d0 / vol**2
     . + (d2tau/tau - (dtau/tau)**2)*dxdv**2
     . + d2xdv2 * dtau/tau
     . + 1.5d0 * sn(ise) * ( d2thet*detdv**2
     . + dthet*(d2etdv2-sn(ise)*dthet*detdv**2/znt))/znt )
c
      return
      end
