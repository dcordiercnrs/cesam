c*********************************************************************

c Programme de comparaison et de test d'opacite sur un <<profil>>
c d'etoile, cad sur un ensemble de valeurs m/Mstar, T, ro, Xi
c correspondant a un modele reel

c             VERSION  4  utilisation de <<TOPS>> version 4.0 !!!!

c ATTENTION : ce programme utilise en entree le fichier de type
c             <<profil_2.dat>> au lieu de <<profil.dat>> (il faut
c             la compo. chim. complete!!

c*********************************************************************


c D. Cordier, Mai 1998

      implicit none

      include 'cesam_3.parametres'
      include 'modele_3.common'
      include 'evol_chim_3.common'

      integer i, nl_max, n_c, jx

      real*8 kappa, dkapdt, dkapdr, dkapdx

      character filename1*50, filename2*50

      parameter(nl_max=5000)

      real*8 ro(nl_max), t(nl_max), z(nl_max), x(nl_max), 
     +       m_mstar(nl_max), xh(18), xh_c(nl_max,18)

      print*, 'Fichier de profil :'
      read(5,'(a)') filename1
      print*, ' '
      print*, 'Fichier de sortie :'
      read(5,'(a)') filename2

      print*, ' '
      print*, 'Valeur de z0 = ?'
      read*, z0
      print*, 'Z0= ', z0

      ihe4=3

      print*, ' '
      opa1='../../../cesam3.2/opa/opacite.dat'
      print*, 'opa1= ', opa1

      print*, ' '
      opa2='../data/opa_tops4_0.02.dat'
      print*, 'opa2= ', opa2
      opa3='../data/opa_tops4_0.008.dat'
      print*, 'opa3= ', opa3
      opa4='../data/opa_tops4_0.004.dat'
      print*, 'opa4= ', opa4

      open(1,status='old',file=filename1)

      n_c=0

      do i= 1, nl_max

         read(1,1000,end=100) m_mstar(i), t(i), ro(i),
     +                       (xh_c(i,jx),jx=1,18)
         n_c=n_c+1

      end do

100   print*, 'Fichier de profil lu!'

c Ecriture du fichier de sortie des opacites TOPS
      
      open(2,status='unknown',form='formatted',file=filename2)

      do i= 1, n_c

         do jx= 1, 18
            xh(jx)=xh_c(i,jx)
         end do
      
c         print*, 'xh(1)= ', xh(1)
c         print*, 'z= ', 1.d0-xh(1)-xh(2)-xh(3)
 
c         print*, t(i), ro(i)
         call opal96_tops41(xh,t(i),ro(i),kappa,dkapdt,dkapdr,
     +                      dkapdx)

c      print*, 'm_mstar= ', m_mstar(i),' kappa : ', kappa

      write(2,2000) m_mstar(i), kappa

1000  format(1p21E19.11)
2000  format(2D10.3)

      end do    

      close(2)  

      end
