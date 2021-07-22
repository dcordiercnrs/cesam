c
c******************************************************************************
c
	subroutine neweta (t,vol,sne,ier)
c
c     calculate degeneracy parameter and its                           c
c     derivatives. evaluate theta and its derivatives                  c
c     ................................
c     (modified from MHD package)
c     ................................
      implicit real*8 (a-h,o-z)
	implicit integer(i-n)
      logical ok
      common /marche/ok
c
      parameter (mfd  =  5)
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
      common /cstmhd/
     .                                camu   ,               cc     ,
     .                                ce     ,               ch     ,
     .                                ck     ,               cme    ,
     .                                cpi    ,               cevw   ,
     .                                cx     ,               cf4    ,
     .                                carad
c
      data   niter,   conv
     .    /    20, 1.d-10 /
c	logical ok
c	common /marche/ok
c
c........... the following factorizing in gse is made to avoid
c........... overflows on simple machines with ranges of about 1.e37
c
      gse  = 2.*cpi*(cme/ch)*(ck/ch)
      gse  = gse**1.5
c
      ceta = sqrt( cpi ) / ( 4.0 * gse )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialization                                                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      ier = 0
      rhs = ceta * sne   / ( t**1.5 * vol )
ccc      write(6,*) ' rhs,ceta = ',rhs,ceta
c
c........... initialization to a few per cent accuracy
c........... (see ref. in dappen, astron. astrphys., 1980)
c
      if(rhs.lt.1.d-6) then
            eta = dlog(rhs)
      else if(rhs.lt.100.d0) then
            g   =rhs
            g2  =rhs*rhs
            et0 =g+(g2/(4.45+2.257*g+g2))*dexp((1.32934*g)**0.66666667)
            eta = dlog(et0)
      else
            eta = (1.32934*rhs)**0.66666667
      end if
      eta = eta + 0.120782
c
ccc      write(6,*) 'eta,rhs ini: ',eta,rhs
      call ferdir( eta, fd )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     newton-raphson iteration
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      iter  = 0
    1 deta  = 2.0d0 * ( rhs - fd(2) ) / fd(1)
      eta = eta + deta
      call ferdir( eta, fd )
      if( dabs(deta) .le. conv ) go to 3
      iter = iter + 1
      if( iter .le. niter ) go to 1
c
c     failure to converge
      write  ( 6, 2 ) t, vol, sne, eta
    2 format ( ' nonconvergence of degeneracy parameter ' ,/,
     .         ' t =', 1pg10.2 ,' vol =', g10.2, ' ne =', g10.3,
     .         ' eta   =,', g12.4 )
	ok=.false.
      ier = 1
      return
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     convergence                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
    3 exeta = dexp( - eta )
c
      thet   =  0.5d0 *   fd(1) / fd(2)
      dthet  = thet * ( fd(4) / fd(1) - thet )
      d2thet = thet * ( fd(5) / fd(1) - 3.0d0 * dthet - thet**2 )
c
      detdv   = - 1.0d0 / ( thet * vol )
      detdt   = - 1.5d0 / ( thet * t )
      detdn   =   1.0d0 / ( thet * sne )
      d2etdn2 = - detdn**2 * fd(4) / fd(1)
      d2etdnt = - dthet * detdn * detdt / thet
      d2etdnv = - dthet * detdn * detdv / thet
      d2etdtv = - dthet * detdt * detdv / thet
      d2etdt2 = - ( 1.0d0 / t   + dthet * detdt / thet ) * detdt
      d2etdv2 = - ( 1.0d0 / vol + dthet * detdv / thet ) * detdv
c
      if(iter.gt.5) write (6,*)' slow convergence in neweta: iter,eta',
     . ',t,vol = ',iter,eta,t,vol
c
ccc      write(6,*) 'eta,fd ',eta,fd
      return
      end
