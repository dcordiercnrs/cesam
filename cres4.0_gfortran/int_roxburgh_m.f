c**************************************************************************

      subroutine int_roxburgh_m( jlim, r_zc1, m_zc1, Hp1, n, bp, q, qt, 
     +              knot, chim, mc, mct, nc, knotc, r_rox, m_rox, etat, 
     +              opa, conv, test, calib1, calib2 )

c ATTENTION : toutes les routines specifiques appelees par 
c 'int_roxburgh_m' sont dans 'int_roxburgh'.

c Programme de calcul des integrales de Roxburgh ---> estimation de la
c taille du coeur convectif.

c auteur : D. Cordier, Fevrier 1999, Octobre 1999, Novembre 1999

c DANS CETTE VERSION on melange toute la zone ou on calcule les integrales
c de Roxburgh.

c Entrees :
c     * jlim : indice de la limite classique (det. dans 'lim_zc_3')
c     * r_zc1: rayon a la limite classique (en 'rsol')
c     * Hp1  : echelle de hauteur de pression a la limite classique (en cm)
c     * n : nbre de couches
c     * n, bp, q, qt, knot : donnent p, t, r, l, m
c       ('ne' est transporte par COMMON)
c     * chim, mc, mct, nc, knotc : donnent la composition chimique.
c     * etat, opa, conv : routines d'equation d'etat, d'opacite et de
c       convection.

c Sortie(s) :
c     * r_rox : rayon du coeur convectif estime par les integrales de Roxburgh.
c       ATTENTION : 'r_rox' est en Rsol (pas besoin de convertir dans 'lim_zc_3'

c Options : suivant les valeurs de 'para_rox1'(actif) et para_rox2'(inactif pour le
c           moment) on tient compte de la dissipation visqueuse ou non.

c           En pratique a la date du 23 Fev. 99 : si 'para_rox1' .gt. 0. alors
c           il y a dissipation visqueuse (cf. le code pour la methode exacte)

c Deroulement du calcul:

c (1) on alcule differentes grandeurs : delad, delrad, Hp, delconv, Vconv 
c     on en deduite Lrad et Lnuc en tout point
c     ATTENTION : on detecte s'il y a un coeur convectif, s'il n'y en a pas on
c                 fait 'r_rox'=-100. et RETURN

c (2) le critere de Roxburgh est applique en calculant les integrales

c -----------------------------------------------------------------------------

c                         AVERTISSEMENTS DIVERS :

c  (i)   l'integrande 1 (left hand side) est dans une unite exotique !!!
c  (ii)  l'integrande 2 (right hand side)est dans une unite exotique !!!
c  (iii) dans cette version on ne se sert pas de Vconv

c

      implicit none

      include 'cesam_3.parametres'

      include 'atmosphere_3.common'
      include 'ctephy.common'
      include 'modele_3.common'
      include 'evol_chim_3.common'

      integer n, ic, knot, lq, nc, knotc, lc, ie, i_coeurconv, i_limzc,
     +        jlim(1), irox, imaxYc

      parameter ( i_coeurconv= 5 ) ! Il faut un coeur conv. de au moins 3 couches

      logical init, first, test, calib1, calib2, mix

c      parameter ( test= .true. ) 
c, calib1=.false., calib2=.false. )
c     'test = .true.'     : differentes grandeurs sont listees dans des
c                           fichiers, la routine est stoppee ensuite.
c     'calib1'            : option permettant de calibrer le parametre 'para_rox1'
c                           afin d'obtenir un overshoot egale a 'ov_sht'
c     'calib2'            : option permettant de calibrer le parametre 'para_rox2'
c                           afin d'obtenir un overshoot egale a 'ov_sht'
      real*8 ov_sht

      real*8 bp(1),q(1),qt(1),f(pne),dfdq(pne),chim(1),mc(1),mct(1),
     +   integ_rox1(pn),compx(pnelem),dcompx(pnelem),integ_rox2(pn),
     +   nuc_m(pnelem),petit_g,Fconv(pn),masse(pn),rayon(pn),ro,rho(pn),
     +   drop,drot,drox,drott,drotp,drotx,u,dup,dut,dux,dutt,dutp,dutx,
     +   nh1,nhe1,nhe2,lamb,delta,cpp,cp,kappa(pn),dkapdt,dkapdr,dkapdx,
     +   delad(pn),delrad(pn),delconv(pn),p(pn),t(pn),l(pn),m(pn),r(pn),
     +   vconv(pn),krad,gravite,taur,s_rox1,s_rox2,Lrad(pn),Lnuc(pn),
     +   xchim(pnelem),dxchim(pnelem),pt(pn),vc_moy, Hp(pn),s_vconv, 
     +   mass, x(pn,pnelem), xh(pn), y, yc

      real*8 r_zc1, r_ov1, m_zc1, Hp1, pourcent, largeur, A1, A2, l1, l2

      real*8 s1, s2, x1, x2, y11, y21, y12, y22

      real*8 bid1,bid2,bid3,bid4,bid5,bid6,bid7,bid8,bid9

      external etat, opa, conv

      real*8 r_rox, m_rox, dm_rox, diff, age, mroxSmsch

      common/roxburgh/diff, age, mroxSmsch
 
      data init /.true./

c-------------------------------------------------------------------------

      ov_sht=OVSHTS
c (0) Message de bienvenue ...

      if ( init ) then
         write(6,*) ' '
         write(6,*) '     *****************************   '
         write(6,*) ' '
         write(6,*) '         Integrales de Roxburgh'
         write(6,*) ' '

         write(6,*) '       ---------------------------'
         write(6,*) '       Version avec MELANGE ! Ok ?'
         write(6,*) '       ---------------------------'
         write(6,*) ' '
         pause

         if ( calib1 ) then
            write(6,*) 'ATT. :  on calibre le parametre 1 '
            write(6,*) 'a la valeur : ', ov_sht
            write(6,*) 'Ok ?'
            pause
         end if
         if ( calib2 ) then
            write(6,*) 'ATT. :  on calibre le parametre 2 '
            write(6,*) 'a la valeur : ', ov_sht
            write(6,*) 'Ok ?'
            pause
         end if

         write(6,*) '     *****************************   '
         write(6,*) ' '
         init=.false.
      end if

      mix=.true.

c (1) Calcul des differentes grandeurs : p, t, ..., une fois pour toutes 
c     (on garde ces grandeurs meme lorsqu'on melange)

      do ic= 1, n ! Iterations sur les couches du modele

         call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,q(ic),lq,f,dfdq) !'mpr' est transportee par COMMON

         p(ic)=exp(f(1)) ! en cgs
         t(ic)=exp(f(2)) ! en K
         l(ic)=sqrt(abs(f(4)))**3 ! en lsol
         m(ic)=sqrt(abs(f(5)))**3 ! en msol
         r(ic)=max(sqrt(abs(f(3))),1.d-30) ! Pour eviter les divisions par zero
                                       ! r est en rsol

c        La chimie :

         call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,knotc,.true.,min(f(5),
     +                mc(nc)),lc,compx,dcompx) !'nbelem' transporte per COMMON

         do ie= 1, nbelem
            xchim(ie)=abs(compx(ie))
            dxchim(ie)=dcompx(ie)
         end do

         call chim_gram_3(xchim,dxchim,nuc_m)

         do ie= 1, nbelem
            x(ic,ie)=xchim(ie)
         end do
c         print*, 'ic= ', ic, ' y= ', xchim(ihe4)
         xh(ic)=xchim(1)

      end do

c (2) Détermination de l'indice maximum 'imaxYc' où on a encore Y=Yc

      Yc=x(3,ihe4-1)+x(3,ihe4)
      Y=Yc
      ic=1
      do while ( abs(y-yc)/yc .le. 1.d-3 ) 
         ic=ic+1
         y=x(ic,ihe4-1)+x(ic,ihe4)
      end do
      imaxYc=ic-1
c      print*, 'imaxYc= ', imaxYc
c      print*, 'y=', x(imaxYc,ihe4-1)+x(imaxYc,ihe4)
c      pause

c (3) Calcul de Delad, Delrad, Hp, rho, Delreel, Vconv, Lrad, Lnuc et
c     des integrandes de Roxburgh pours toutes les couches telles que
c     ic .le. imaxYc

 100  do ic= 1, imaxYc

         do ie= 1, nbelem
            xchim(ie)=x(ic,ie)
         end do

         call etat(p(ic),t(ic),xchim,.false.,ro,drop,drot,drox,drott,
     +             drotp,drotx,u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,
     +             nhe2,lamb)

         delta=-t(ic)/ro*drot
         cpp=p(ic)/ro/t(ic)*delta
         cp=dut+cpp

         delad(ic)=p(ic)/ro*delta/cp/t(ic)

         call opa(xchim,t(ic),ro,kappa(ic),dkapdt,dkapdr,dkapdx)

         delrad(ic)=3./16./pi/aradia/clight/g*kappa(ic)*l(ic)*lsol*p(ic)
     +              /m(ic)/msol/t(ic)**4

         Hp(ic)=p(ic)*rsol**2*r(ic)**2/g/msol/m(ic)/ro ! Hp est ici en cm

         rho(ic)=ro

c        Detection du coeur convectif :

         if ( ic .eq. i_coeurconv ) then
            if ( delrad(ic) .lt. delad(ic) ) then
               r_rox=-100.d0
c	       pause
               return
            end if
         end if

c        Calcul du gradient donne par la MLT

         krad=4./3.*aradia*clight/kappa(ic)/ro*t(ic)**3 !Conductivite radiative

         gravite=G*m(ic)*msol/rsol**2/r(ic)**2

         taur=kappa(ic)*ro*alpha*Hp(ic)  !Epaisseur optique de la bulle


         call conv(krad,gravite,delta,cp,ro,Hp(ic),taur,delrad(ic),
     +             delad(ic),.false.,delconv(ic),bid1,bid2,bid3,bid4,
     +             bid5,bid6,bid7,bid8,bid9)

c        Flux convectif :

         if ( delrad(ic) .ge. delconv(ic) ) then
            Fconv(ic)=4.*aradia*clight*G/3.*t(ic)**4*m(ic)*msol/P(ic)
     +                /kappa(ic)/r(ic)**2/rsol**2
     +                *(delrad(ic)-delconv(ic))

            petit_g=G*m(ic)*msol/r(ic)**2/rsol**2

            Vconv(ic)=1./sqrt(8.)*(petit_g*delta*alpha*Hp(ic)*
     +                4.*sqrt(2.)/ro/cp/t(ic)*Fconv(ic))**(1./3.)
         else
            vconv(ic)=0.d0
         end if

c        Calcul de Lrad :

         Lrad(ic)=16.*pi/3.*aradia*clight*G*m(ic)*msol*t(ic)**4
     +            /kappa(ic)/p(ic)*delad(ic)

c        Calcul de Lnuc :

         Lnuc(ic)=l(ic)*lsol

c        La masse :

         masse(ic)=m(ic)*msol

c        Le rayon :

         rayon(ic)=r(ic)*rsol

c        Le produir pression x temperature :

         pt(ic)=p(ic)*t(ic)

c        Dans le cas ou on met de la dissipation visqueuse :

         if ( para_rox1 .ne. 0.d0 .AND. para_rox2 .ne. 0.d0 ) then
            write(6,*) 'Pb. dans ''int_roxburgh'' !'
            write(6,*) ' '
            write(6,*) 'para_rox1 = ', para_rox1
            write(6,*) 'para_rox2 = ', para_rox2
            write(6,*) ' '
            stop
         else
            if ( para_rox1 .eq. 0.d0 .AND. para_rox2 .eq. 0.d0 ) then
               integ_rox2(ic)=0.d0
            else
               if ( para_rox1 .ne. 0.d0 ) then
                  if ( .NOT. calib1 ) then
                  integ_rox2(ic)=para_rox1/t(ic)*ro*4*pi*rayon(ic)**2 ! cf. cahier XII p. 38v
     +                           /rsol**2 ! en unite solaire            et cah. XIII p. 26v
                  else
c                 CALIBRATION, a ce stade 'para_rox=1' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  integ_rox2(ic)=1./t(ic)*ro*4*pi*rayon(ic)**2 ! cf. cahier XII p. 38 verso
     +                /rsol**2 ! en unite solaire                et cah. XIII p. 26v
                  end if
               end if
               if ( para_rox2 .ne. 0.d0 ) then
                  integ_rox2(ic)=1./2./Hp(ic)/t(ic)*rho(ic)*4.*pi*
     +                           rayon(ic)**2/rsol**2
               end if
            end if
         end if

      end do ! ic 


c ------------------------------------------------------------------------------
c ------------------------------------------------------------------------------

c         *  det. de l'indice du bord du coeur convectif
c         *  det. de la vitesse convective moyenne
c         *  det. extrapolation eventuelle de log10(Lrad) :
c            - calcul des coeff. de la droite d'extrapolation
c            - calcul des valeurs extrapolee de Lrad
c         *  l'integrande de l'integrale de Roxburgh est calcule

c     det. de l'indice de la limite ZC/ZR

      first=.true.

      do ic= 1, n
         if ( delrad(ic) .le. delad(ic) .AND. first ) then
            first=.false.
            i_limzc=ic-1
         end if
      end do

      if ( i_limzc .ne. jlim(1) ) then
          print*, ' '
          print*, 'jlim(1) .ne. i_limzc !'
          print*, 'jlim(1)= ', jlim(1)
          print*, 'i_limzc= ', i_limzc
c          print*, 'Que fait-on ?'
c          pause
          print*, ' '
          print*, 'On fait : i_limzc=jlim(1)'
          i_limzc=jlim(1)
      end if

c     det. de la vitesse convective moyenne

      if ( para_rox2 .ne. 0.d0 ) then
         call integrale(Vconv,rayon,n,i_limzc,s_vconv)
         vc_moy=s_vconv/rayon(i_limzc)
      end if

c      print*, 'Vitesse convective moyenne : ', vc_moy
c      pause

c     On calcule l'integrande de l'integrale de Roxburgh

      integ_rox1(1)=0.d0

      do ic= 1, n
         integ_rox1(ic)=-(Lrad(ic)-Lnuc(ic))*G*delad(ic)*masse(ic)
     +    /rayon(ic)**2/pt(ic)*rho(ic) ! cf. formule p.38-39, cahier XII, pt: produit P*T
     +   /lsol/G/msol*rsol**2 ! en unite de G et solaires !!! (cf. cahier XII p. 38 verso)
      end do

c     Deuxieme integrale de Roxburgh (cas ou on utilise la vitesse moyenne)

      integ_rox2(1)=0.d0

      if ( para_rox2 .ne. 0.d0 ) then
         do ic= 1, n
c           ATTENTION !!! Le parametre n'est pas au cube !!!
            if ( .NOT. calib2 ) then
               integ_rox2(ic)=integ_rox2(ic)*vc_moy**3*para_rox2
            else
               integ_rox2(ic)=integ_rox2(ic)*vc_moy**3
            end if
         end do
      end if

c     ecriture dans des fichiers en cas de test

      if ( test ) then
         call ouverture_fichiers
         print*, ' '
         print*, 'Ecriture des fichiers Roxburgh ...'
         print*, 'ATTENTION on utilise ''mtot'' et pas ''mstar'' !'
         print*, ' '
         do ic= 1, n
         mass=masse(ic)/msol/mtot
         call ecriture( mass,rayon(ic)/rsol,log10(lrad(ic)/lsol),
     +       log10(lnuc(ic)/lsol), delad(ic),delrad(ic),kappa(ic),
     +       xh(ic),Vconv(ic), Fconv(ic),Hp(ic), Integ_rox1(ic), 
     +       Integ_rox2(ic) )
         end do
      end if
c --------------------------------------------------------------------------
c --------------------------------------------------------------------------

c Calcul des integrales de Roxburgh et recherche de la limite
c              du coeur convectif.

      do ic= 1, n
         rayon(ic)=rayon(ic)/rsol
      end do

      if ( (.NOT. calib1) .AND. (.NOT. calib2) ) then

      first=.true.

      do ic= 1, imaxYc
         call integrale(integ_rox1,rayon,n,ic,s_rox1)
c         print*, ' '
c         print*, 's_rox1= ', s_rox1/1.d-32
         call integrale(integ_rox2,rayon,n,ic,s_rox2)
c        print*, 's_rox2= ', s_rox2/1.d-32

         if ( s_rox1 .lt. s_rox2 ) then
            if ( ic .ge. i_limzc ) then
               if ( first ) then
            call integrale(integ_rox1,rayon,n,ic-1,s1)
            call integrale(integ_rox2,rayon,n,ic-1,s2)
                  if ( s1 .gt. s2 ) then
                     first=.false.
                     irox=ic
                     print*, 'irox= ', irox
                  end if
               end if
            end if
         end if

         if ( test ) then
            s_rox1=s_rox1*1.d+22
            s_rox2=s_rox2*1.d+22
            write(54,1000) rayon(ic), s_rox1
            write(55,1000) rayon(ic), s_rox2
            s_rox1=s_rox1/1.d+22
            s_rox2=s_rox2/1.d+22
         end if

      end do ! ic

 1000 format(1p2d13.6)

      if ( first ) then
         if ( imaxYc .eq. n ) then
            print*, ' '
            print*, 'Pb. dans ''int_roxburgh'' :'
            print*, 'La limite n''a pas ete trouvee !'
            print*, ' '
            print*, 'On fait r_rox=-100., ok ?'
            pause
            return
c           stop
         else
            imaxYc=imaxYc+1
c            print*, 'Avant melange :'
c            print*, 'x(5,3)= ', x(5,3)
            call melange( x, masse, imaxYc )
c            print*, 'Apres melange :'
c            print*, 'x(5,3)= ', x(5,3)
c            pause
            go to 100
         end if
      end if

c Affinement de la limite de Roxburgh :

c Par dichotomie (pour tests) :
c      call rech_dich_lim_rox( integ_rox1, integ_rox2, n, irox, rayon, 
c     +                              r_rox )

c Methode analytique :

      call integrale(integ_rox1,rayon,n,irox-1,s1)

      call integrale(integ_rox2,rayon,n,irox-1,s2)

      x1=rayon(irox-1) ! en rsol
      x2=rayon(irox)   ! en rsol

      y11=integ_rox1(irox-1) ! INTEGRANDE, le deuxieme indice est celui de
                             ! l'integrande de Roxburgh
      y21=integ_rox1(irox)   ! idem

      y12=integ_rox2(irox-1) ! idem
      y22=integ_rox2(irox)   ! idem

      print*, ' '
      print*, 'Appel de <<affine>> ...'

      call affine( s1, s2, x1, x2, y11, y21, y12, y22, r_rox )
      pourcent=int((r_rox-x1)/(x2-x1)*100.)
      print*, 'Limite situee a  : ', pourcent, ' % entre les couches ',
     + irox-1, ' et ', irox
      print*, '... apres <<affine>> , r_rox= ', r_rox
      largeur=(rayon(irox)-rayon(irox-1))/Hp(irox-1)*rsol
      print*, 'Largeur de la couche en Hp : ', largeur

c en cas de test on ferme tous les fichiers

      if ( test ) then
         call fermeture
      end if

c     On interpole pour obtenir la masse de la limite de Roxburgh
      call interpol_lin(masse, rayon, n, r_rox, m_rox, dm_rox)
      m_rox=m_rox/msol
c      print*, 'm_rox= ', m_rox
c      print*, 'm_rox/m_zc1= ', m_rox/m_zc1
      if ( m_zc1 .gt. 0.d0 ) then
         mroxSmsch=m_rox/m_zc1
      else
         print*, 'Pb. dans ''int_roxburgh'' : m_zc1=0'
         stop
      end if

      return

      end if ! if .NOT. calib1 .AND. .NOT. calib2

c Dans le cas d'une calibration du parametre 'para_rox2' : --------------------
      if ( calib1 ) then
c        recherche de l'indice de la limite overshootée
         r_ov1=r_zc1+ov_sht*Hp1/rsol ! Hp1 entre dans la routine en 'cm'
         first=.true.
         do ic= 1, n-1
            if ( (r_ov1 .ge. rayon(ic)) .AND. ! rayon est ici en rsol
     +           (r_ov1 .lt. rayon(ic+1)) ) then
               if ( first ) then
                  first=.false.
                  irox=ic+1
               end if
            end if
         end do ! ic

         if ( first ) then
            print*, 'Pb. dans ''int_roxburgh'' : on ne trouve pas la couc
     +he correspondant a r_ov1'
            stop
         end if
         if ( mix ) then
            mix=.false.
            call melange( x, masse, irox )
            go to 100
         end if
c        Calcule de la premiere integrale
         call integrale(integ_rox1,rayon,n,irox-1,s1)
         call integrale2(rayon(irox-1),rayon(irox),integ_rox1(irox-1),
     +                   integ_rox1(irox),r_ov1,A1)
c        l'integrale vaut alors 's1+A1'
c        Calcule de la deuxieme integrale
         call integrale(integ_rox2,rayon,n,irox-1,s2)
         call integrale2(rayon(irox-1),rayon(irox),integ_rox2(irox-1),
     +                   integ_rox2(irox),r_ov1,A2)
c        l'integrale vaut alors 's2+A2'
         if ( s1+A1 .ne. 0.d0 ) then
	open(unit=7,status='unknown',access='append',form='formatted',
     +       file='para_calib1.dat')
            l1=(s1+A1)/(s2+A2)
            r_rox=r_ov1
            write(6,'(A31,D13.6)') '>>>>>>>>>>>>>>>>>>> lambda 1 = ', l1
	    write(7,'(d13.6)') l1, age
         close(unit=7)

c     On interpole pour obtenir la masse de la limite de Roxburgh
         call interpol_lin(masse, rayon, n, r_rox, m_rox, dm_rox)
         m_rox=m_rox/msol
         if ( m_zc1 .gt. 0.d0 ) then
            mroxSmsch=m_rox/m_zc1
         else
            print*, 'Pb. dans ''int_roxburgh'' : m_zc1=0'
            stop
         end if
            return
         else
            print*, 'Pb. dans ''int_roxburgh'' : ''s1+A1''=0. !'
            stop
         end if ! ( s1+A1 .ne. 0.d0 )
      end if ! if calib1 = .true.
c ++++++++++++++++++++++++++++++++++
      if ( calib2 ) then
c        recherche de l'indice de la limite overshootée
         r_ov1=r_zc1+ov_sht*Hp1/rsol ! Hp1 entre dans la routine en 'cm'
         first=.true.
         do ic= 1, n-1
            if ( (r_ov1 .ge. rayon(ic)) .AND. ! rayon est ici en rsol
     +           (r_ov1 .lt. rayon(ic+1)) ) then
               if ( first ) then
                  first=.false.
                  irox=ic+1
               end if
            end if
         end do ! ic

         if ( first ) then
            print*, 'Pb. dans ''int_roxburgh'' : on ne trouve pas la couc
     +he correspondant a r_ov1'
            stop
         end if
         if ( mix ) then
            mix=.false.
            call melange( x, masse, irox )
            go to 100
         end if
c        Calcule de la premiere integrale
         call integrale(integ_rox1,rayon,n,irox-1,s1)
         call integrale2(rayon(irox-1),rayon(irox),integ_rox1(irox-1),
     +                   integ_rox1(irox),r_ov1,A1)
c        l'integrale vaut alors 's1+A1'
c        Calcule de la deuxieme integrale
         call integrale(integ_rox2,rayon,n,irox-1,s2)
         call integrale2(rayon(irox-1),rayon(irox),integ_rox2(irox-1),
     +                   integ_rox2(irox),r_ov1,A2)
c        l'integrale vaut alors 's2+A2'
         if ( s1+A1 .ne. 0.d0 ) then
	open(unit=7,status='unknown',access='append',form='formatted',
     +       file='para_calib2.dat')
            l1=(s1+A1)/(s2+A2)
            r_rox=r_ov1
            write(6,'(A31,D13.6)') '>>>>>>>>>>>>>>>>>>> lambda 2 = ', l2
	    write(7,'(2d13.6)') l2, age
         close(unit=7)
c     On interpole pour obtenir la masse de la limite de Roxburgh
         call interpol_lin(masse, rayon, n, r_rox, m_rox, dm_rox)
         m_rox=m_rox/msol
         if ( m_zc1 .gt. 0.d0 ) then
            mroxSmsch=m_rox/m_zc1
         else
            print*, 'Pb. dans ''int_roxburgh'' : m_zc1=0'
            stop
         end if
            return
         else
            print*, 'Pb. dans ''int_roxburgh'' : ''s1+A1''=0. !'
            stop
         end if ! ( s1+A1 .ne. 0.d0 )
      end if ! if calib2 = .true.

c -----------------------------------------------------------------------------

      end

c**************************************************************************

      subroutine melange( x, m, imax )

c Routine de melange de la composition chimique jusqu'a l'indice 'imax'
      implicit none

      include 'cesam_3.parametres'

      include 'evol_chim_3.common'

      integer imax, i, ie

      real*8 x(pn,pnelem), m(pn), m_tot, mtotx

c m_tot : masse de l'étoile jusqu'à 'imax'

      m_tot=0.d0

      do i= 1, imax
         m_tot=m_tot+m(i)
      end do

      do ie= 1, nbelem
         mtotx=0.d0
         do i= 1, imax
            mtotx=mtotx+x(i,ie)*m(i)
         end do
         do i= 1, imax
            x(i,ie)=mtotx/m_tot
         end do
      end do

      return

      end

c**************************************************************************












