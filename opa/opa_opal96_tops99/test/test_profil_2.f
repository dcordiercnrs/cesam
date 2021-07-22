c*********************************************************************

c Programme de comparaison et de test d'opacite sur un <<profil>>
c d'etoile, cad sur un ensemble de valeurs m/Mstar, ro, T, X et Z
c correspondant a un modele reel

c             VERSION  2       (utilisation de <<opa_ff>> au lieu
c                      <<TOPS>> pour le programme <<test_profil.f>>)


c ATTENTION : ce programme utilise en entree le fichier de type
c             <<profil_2.dat>> au lieu de <<profil.dat>> (il faut
c             la compo. chim. complete!!

c*********************************************************************


c D. Cordier, Fevrier 1998

      implicit none

      integer i, nl_max, n_c, jx

      real*8 lr, t6, bid1, bid2, bid3, kap1, kap2

      character opa1*31, opa2*31

      parameter ( nl_max=5000, opa1='../../cesam3.2/opa/opacite.dat',
     +            opa2='opa_tops.dat')

      real*8 ro(nl_max), t(nl_max), z(nl_max), x(nl_max), 
     +       m_mstar(nl_max), xh(12), xh_c(nl_max,12)

      open(1,status='old',file='profil_2.dat')

      n_c=1

      do i= 1, nl_max

         read(1,1000,end=100) m_mstar(i), ro(i), t(i),
     +                       (xh_c(i,jx),jx=37,48)
         n_c=n_c+1

      end do

100   print*, 'Fichier de profil lu!'


c Ecriture du fichier de sortie des opacites TOPS
      
      open(2,status='old',file='k_opaff.dat')

      do i= 1, nl_max

         do jx= 1, 12
            xh(jx)=xh_c(i,jx)
         end do

      print*, 'T= ', t(i), ' ro= ', ro(i)

      call opa_ff_test(xh,t(i),ro(i),kap1,bid1,bid2,bid3) 

      print*, 'PB 1'

      write(2,2000) m_mstar(i), kap1
2000  format(2D10.3)

      end do    

      close(2)  

c Ecriture du fichier de sortie des opacites OPAL

      open(3,status='old',file='k_opal.dat')

      do i= 1, nl_max

      x(i)=xh_c(i,1)
      z(i)=1.d0-x(i)-xh_c(i,2)-xh_c(i,3)

      if ( z(i) .le. 0.1d0 ) then

      t6=t(i)/10.**6
      lr=log10(ro(i)/t6**3)

      call kappa_opal_3(lr,t6,ro(i),x(i),z(i),kap2,bid1,bid2,
     +                  bid3,opa1)

      write(3,2000) m_mstar(i), kap2

      end if

      end do

      close(3)

1000  format(15D10.3)

      end
