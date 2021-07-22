

      subroutine opa_ff_test( xh, t, ro, kappa, dkapdt, dkapdr, dkapdx)

c VERSION TEST : on a supprime des securites pour pouvoir faire des
c                profils


c Opacite de type Kramers (absorption free-free des electrons et
c diffusion par ces memes electrons (cf. Kippenhahn p. 137 et 138)
c pour les formules)

c But : prolonger les tables d'Yveline au-dela de Z=0.1
c       (cas du coeur d'une cepheide, ou on a combustion de He)

c 29 Mai 1997

c Domaine de validite : X=0. !!!!

      implicit none

      include 'cesam_3.parametres'
      include 'ctephy.common'
      include 'modele_3.common'
      include 'evol_chim_3.common'

      integer i

      real*8 xh(1), t, ro, kappa1, B, kappa, dkapdt, dkapdr, dkapdx

      real*8 num_atom(18), kappa_ff, kappa_sc

      data num_atom/1., 2., 2., 6., 6., 7., 7., 7.,
     +              8., 8., 8., 9., 10., 10., 12., 12.,
     +              12., 1./


      B=0.

c      if ( nchim .ne. 18 ) then   ! Mesure de securite
c         print*, 'Erreur 1 dans opa_ff !'
c         stop
c      end if
c      if ( ihe4 .ne. 3 ) then     ! Mesure de securite
c         print*, 'Erreur 2 dans opa_ff !'
c         stop
c      end if

      ihe4=3   ! version TEST!!!!!

      do i= ihe4+1, nchim, 1
         B=B+xh(i)*(num_atom(i))**2./nucleo(i)
      end do

      kappa1=3.8E22*(1.+xh(1))*(xh(1)+xh(ihe4)+B)

      kappa_ff=kappa1*ro*t**(-7./2.)

      kappa_sc=0.2

      kappa=kappa_ff+kappa_sc

      dkapdt=-7./2.*kappa1*ro*t**(-9./2.)

      dkapdr=kappa1*t**(-7./2.)

      dkapdx=3.8E22*(1.+2.*xh(1)+xh(ihe4)+B)*ro*t**(-7./2.)

      return


      end
