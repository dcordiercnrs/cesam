c**************************************************************************

      subroutine int_roxburgh( jlim, r_zc1, m_zc1, Hp1, n, bp, q, qt, 
     +       knot, chim, mc, mct, nc, knotc, r_rox, m_rox, etat, opa, 
     +       conv, test, calib1, calib2 )

c Programme de calcul des integrales de Roxburgh ---> estimation de la
c taille du coeur convectif.

c auteur : D. Cordier, Fevrier 1999, Octobre 1999

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

c (2) le critere de Roxburg est applique en calculant les integrales

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

      integer n, ic, lq, nc, knotc, lc, ie, i_coeurconv, i_limzc,
     +        i_recule, jlim(1), irox, knot

c      character rep*1

      parameter ( i_coeurconv= 5 ) ! Il faut un coeur conv. de au moins 3 couches

      logical init, first, extrapol, extrapol_bis, test, calib1, calib2

      parameter ( extrapol = .false. )
c, test= .true. )
c calib1=.false., 
c     +            calib2=.false. )
c     'extrapol = .true.' : on extrapole lineairement log Lrad
c ATTENTION : cette option d'extrapolation existe uniquement pour des
c tests !!!
c     'test = .true.'     : differentes grandeurs sont listees dans des
c                           fichiers, la routine est stoppee ensuite.
c     'calib1'            : option permettant de calibrer le parametre 'para_rox1'
c                           afin d'obtenir un overshoot egale a 'ov_sht'
c     'calib2'            : option permettant de calibrer le parametre 'para_rox2'
c                           afin d'obtenir un overshoot egale a 'ov_sht'

      real*8 ov_sht
c      parameter ( ov_sht = 0.2 )

      real*8 bp(1),q(1),qt(1),f(pne),dfdq(pne),chim(1),mc(1),mct(1),
     +   integ_rox1(pn),compx(pnelem),dcompx(pnelem),integ_rox2(pn),
     +   nuc_m(pnelem),masse(pn),rayon(pn),ro,rho(pn),fconv(pn),
     +   drop,drot,drox,drott,drotp,drotx,u,dup,dut,dux,dutt,dutp,dutx,
     +   nh1,nhe1,nhe2,lamb,delta,cpp,cp,kappa(pn),dkapdt,dkapdr,dkapdx,
     +   delad(pn),delrad(pn),delconv(pn),p,t,l,m,r,Hp(pn),vconv(pn),
     +   krad,gravite,taur,s_rox1,s_rox2,Lrad(pn),
     +   Lnuc(pn),xchim(pnelem),dxchim(pnelem),pt(pn),
     +   a,b, x(pn), y(pn), Amoy, n_e, Zmoy, lambda_c, 
     +   visco_cine(pn), epsi_sur, Phi

      real*8 r_zc1, r_ov1, m_zc1, Hp1, pourcent, largeur, A1, A2, l1, l2

      real*8 s1, s2, x1, x2, y11, y21, y12, y22

      real*8 bid1,bid2,bid3,bid4,bid5,bid6,bid7,bid8,bid9

      external etat, opa, conv

      real*8 r_rox, m_rox, dm_rox, diff, age, mroxSmsch

      common/roxburgh/diff, age, mroxSmsch
 
      data init /.true./, extrapol_bis /.false./

c-------------------------------------------------------------------------

      ov_sht=OVSHTS
c (0) Message de bienvenue ...

      if ( init ) then
         write(6,*) ' '
         write(6,*) '     *****************************   '
         write(6,*) ' '
         write(6,*) '         Integrales de Roxburgh'
         write(6,*) ' '

         write(6,*) ' '
         write(6,*) '       -------------------------'
         write(6,*) '       Version avec SANS melange'
         write(6,*) '       -------------------------'
         write(6,*) ' '

c         if ( test ) then
c            write(6,*) 'Extrapole-t-on ? (o/n)'
c            read(5,'(A)') rep
c            if ( rep .eq. 'o' ) then
c               extrapol_bis = .true.
c            else
c               extrapol_bis = .false.
c            end if
c            write(6,*) ' '
c         end if
         if ( extrapol ) then
            write(6,*) 'ATT. : on extrapole log10(Lrad) !'
            write(6,*) ' ...et seulement cette grandeur !'
            write(6,*) ' ---> OK ?'
            pause
         end if
         if ( calib1 ) then
            write(6,*) 'ATT. :  on calibre le parametre 1 '
            write(6,*) 'a la valeur : ', OVSHTS
            write(6,*) 'Ok ?'
            pause
         end if
         if ( calib2 ) then
            write(6,*) 'ATT. :  on calibre le parametre 2 '
            write(6,*) 'a la valeur : ', OVSHTS
            write(6,*) 'Ok ?'
            pause
         end if

         write(6,*) '     *****************************   '
         write(6,*) ' '
c         if ( extrapol .OR. extrapol_bis ) then
c            pause
c         end if
         init=.false.
      end if

c (1) Calcul des differentes grandeurs : p, t, delad, ...

      do ic= 1, n ! Iterations sur les couches du modele

         call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,q(ic),lq,f,dfdq) !'mpr' est transportee par COMMON

         p=exp(f(1)) ! en cgs
         t=exp(f(2)) ! en K
         l=sqrt(abs(f(4)))**3 ! en lsol
         m=sqrt(abs(f(5)))**3 ! en msol
         r=max(sqrt(abs(f(3))),1.d-30) ! Pour eviter les divisions par zero
                                       ! r est en rsol

c         if ( ic .eq. n ) then
c            print*, ' '
c            print*, 'p= ', p
c            print*, 't= ', t
c            print*, 'l= ', l
c            print*, 'm= ', m
c            print*, 'r= ', r
c            pause
c         end if

c         if ( ic .eq. 2 ) then
c            print*, 'On verifie que m .ne. 0. en i= 2'
c            print*, 'm(1)= ', m
c            print*, 'On verifie que r .ne. 0. en i= 2'
c            print*, 'r(1)= ', r
c            pause
c         end if

c         if ( ic .eq. n ) then
c            print*, 'Verification de l''unite de m, en i=n'
c            print*, 'm(n)= ', m
c            pause
c         end if

c La chimie :

         call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,knotc,.true.,min(f(5),
     +                mc(nc)),lc,compx,dcompx) !'nbelem' transporte per COMMON

         do ie= 1, nbelem
            xchim(ie)=abs(compx(ie))
            dxchim(ie)=dcompx(ie)
         end do

         call chim_gram_3(xchim,dxchim,nuc_m)

         x(ic)=xchim(1)
         y(ic)=xchim(2)+xchim(3)

c Verification de la compo. chim.

c         if ( (ic .eq. 2) .OR. (ic .eq. n) ) then
c         print*, 'Fraction massique pour la compo. chim. ?'
c         print*, 'ic= ', ic
c         do ie= 1, nbelem
c            print*, 'xchim(',ie,')= ',xchim(ie)
c         end do
c         s_rox1=0.
c         do ie= 1, nbelem
c            s_rox1=s_rox1+xchim(ie)
c         end do
c         print*, 'La somme des Xi= ', s_rox1
c         s_rox1=0.
c         pause
c         end if

c Calcul de ro, ...

         call etat(p,t,xchim,.false.,ro,drop,drot,drox,drott,drotp,drotx
     +             ,u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)

         delta=-t/ro*drot
         cpp=p/ro/t*delta
         cp=dut+cpp

         delad(ic)=p/ro*delta/cp/t

         call opa(xchim,t,ro,kappa(ic),dkapdt,dkapdr,dkapdx)

         delrad(ic)=3./16./pi/aradia/clight/g*kappa(ic)*l*lsol*p
     +              /m/msol/t**4

         Hp(ic)=p*rsol**2*r**2/g/msol/m/ro ! Hp est ici en cm

         rho(ic)=ro

c Detection du coeur convectif :

         if ( ic .eq. i_coeurconv ) then
            if ( delrad(ic) .lt. delad(ic) ) then
               r_rox=-100.d0
c	       pause
               return
            end if
         end if

c Calcul du gradient donne par la MLT

         krad=4./3.*aradia*clight/kappa(ic)/ro*t**3 !Conductivite radiative

         gravite=G*m*msol/rsol**2/r**2

         taur=kappa(ic)*ro*alpha*Hp(ic)  !Epaisseur optique de la bulle


         call conv(krad,gravite,delta,cp,ro,Hp(ic),taur,delrad(ic),
     +             delad(ic),.false.,delconv(ic),bid1,bid2,bid3,bid4,
     +             bid5,bid6,bid7,bid8,bid9)

c Calcul de Lrad :

         Lrad(ic)=16.*pi/3.*aradia*clight*G*m*msol*t**4/kappa(ic)
     +            /p*delad(ic)

c Calcul de Lnuc :

         Lnuc(ic)=l*lsol

c La masse :

         masse(ic)=m*msol

c Le rayon :

         rayon(ic)=r*rsol

c Le produit pression x temperature :

         pt(ic)=p*t

c Calcul de la viscosite cinetique, formule de Spitzer (cf. Hansen & Kawaler p.184
c et cahier XV p. 2, 3, 4) :

         Amoy=(x(ic)*nucleo(1)+xchim(2)*nucleo(2)+xchim(3)*nucleo(3)
     +     +(1.-x(ic)-y(ic))*20.)/1.

         n_e=(6.02d23*(x(ic)/nucleo(1)+xchim(2)/nucleo(2)+
     +        xchim(3)/nucleo(3))
     +        +2.3d23*(1.-x(ic)-y(ic))
     +       )*rho(ic)

         lambda_c=10.**4*T**(3./2.)/sqrt(n_e)

         Zmoy=(x(ic)*1+y(ic)*2+(1.-x(ic)-y(ic))*9.9)/1.

         visco_cine(ic)=2.d-15/rho(ic)*T**(5./2.)*sqrt(Amoy)
     +                  /Zmoy**4/log(lambda_c)

c         if ( ic .eq. 1 ) then
c            print*, 'masse(1)      = ', masse(1)
c            print*, 'visco_cine(1) = ', visco_cine(1)
c            pause
c         end if

c Dans le cas ou on met de la dissipation visqueuse :

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
                  integ_rox2(ic)=para_rox1/t*ro*4*pi*rayon(ic)**2 ! cf. cahier XII p. 38v
     +                           /rsol**2 ! en unite solaire        et cah. XIII p.26v
                  else
c                 CALIBRATION, a ce stade 'para_rox=1' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  integ_rox2(ic)=1./t*ro*4*pi*rayon(ic)**2 ! cf. cahier XII p. 38 verso
     +                /rsol**2 ! en unite solaire            et cah. XIII p.26v
                  end if
               end if
               if ( para_rox2 .ne. 0.d0 ) then
                  if ( .NOT. calib2 )  then
c                    Cf. formule p. 4 cahier XV
                     epsi_sur=para_rox2*(2.*visco_cine(ic)*3.87E+10)**3
     +                        *(2.*pi)**4/Hp(ic)**4
                     Phi=rho(ic)*epsi_sur
                     integ_rox2(ic)=Phi/T*4.*pi*rayon(ic)**2
                  else
                     epsi_sur=(2.*visco_cine(ic)*3.87E+10)**3
     +                        *(2.*pi)**4/Hp(ic)**4
                     Phi=rho(ic)*epsi_sur
                     integ_rox2(ic)=Phi/T*4.*pi*rayon(ic)**2
                  end if
               end if
            end if
         end if

      end do ! Iteration sur les couches du modele

c ------------------------------------------------------------------------------
c ----------------- Fin de la partie (1) ---------------------------------------
c ------------------------------------------------------------------------------

c Partie (2) : 
c         *  det. de l'indice du bord du coeur convectif
c         *  det. de la vitesse convective moyenne
c         *  det. extrapolation eventuelle de log10(Lrad) :
c            - calcul des coeff. de la droite d'extrapolation
c            - calcul des valeurs extrapolee de Lrad
c         *  l'integrande de l'integrale de Roxburgh est calcule

c det. de l'indice de la limite ZC/ZR

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

c det. des coeff. de la droite d'extrapolation
      
      if ( extrapol .OR. extrapol_bis) then
         print*, 'extrapol= ', extrapol
         print*, 'extrapol_bis= ', extrapol_bis
         print*, 'On va extrapoler !'
         i_recule=0

 10      ic=i_limzc-i_recule

         a=(log10(Lrad(ic))-log10(Lrad(ic-1)))/(masse(ic)-masse(ic-1))

         b=log10(Lrad(ic))-a*masse(ic)

         if ( a .le. 0 ) then
            i_recule=i_recule+1
            go to 10
         end if

c         print*, 'a= ', a
c         print*, 'b= ', b
c         print*, 'ic= ', ic
c         print*, 'masse= ', masse(ic)
c         print*, 'log10Lard= ', log10(Lrad(ic))

c On change les valeurs de Lrad pour ic > i_limzc (extrapolation)

         do ic= i_limzc, n
            Lrad(ic)=10.**(a*masse(ic)+b)
         end do

      end if

c On calcule l'integrande de l'integrale de Roxburgh

      integ_rox1(1)=0.d0

      do ic= 1, n
         integ_rox1(ic)=-(Lrad(ic)-Lnuc(ic))*G*delad(ic)*masse(ic)
     +    /rayon(ic)**2/pt(ic)*rho(ic) ! cf. formule p.38-39, cahier XII, pt: produit P*T
     +   /lsol/G/msol*rsol**2 ! en unite de G et solaires !!! (cf. cahier XII p. 38 verso)
      end do

      if ( test ) then
         call ouverture_fichiers
      end if

c --------------------------------------------------------------------------
c ----------------- Fin de la partie (2) -----------------------------------
c --------------------------------------------------------------------------

c Partie (3) : calcul des integrales de Roxburgh et recherche de la limite
c              du coeur convectif.

      do ic= 1, n
         rayon(ic)=rayon(ic)/rsol
      end do

      if ( (.NOT. calib1) .AND. (.NOT. calib2) ) then

      first=.true.

      do ic= 1, n
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
            if ( (ic .ge. 2) .AND. (ic .le. irox+10) ) then
               s_rox1=s_rox1*1.d+17
               s_rox2=s_rox2*1.d+17
               write(54,1000) rayon(ic), s_rox1
               write(55,1000) rayon(ic), s_rox2
               s_rox1=s_rox1/1.d+17
               s_rox2=s_rox2/1.d+17
            end if
         end if

      end do ! ic

      if ( test ) then
         print*, ' '
         print*, 'Ecriture des fichiers Roxburgh ...'
         print*, 'ATTENTION on utilise ''mtot'' et pas ''mstar'' !'
         print*, ' '
         do ic= 2, irox+5
         call ecriture( masse(ic)/msol/mtot, rayon(ic)/rsol, 
     +       log10(lrad(ic)/lsol),
     +       log10(lnuc(ic)/lsol), delad(ic),delrad(ic),kappa(ic),x(ic)
     +      ,Vconv(ic), Fconv(ic), Hp(ic), Integ_rox1(ic), 
     +       Integ_rox2(ic) )
         end do
      end if



 1000 format(1p2d13.6)

      if ( first ) then
         print*, ' '
         print*, 'Pb. dans ''int_roxburgh'' :'
         print*, 'La limite n''a pas ete trouvee !'
         print*, ' '
         print*, 'On fait r_rox=-100., ok ?'
c         pause
         return
c         stop
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

      end if ! if (.NOT. calib1) .AND. (.NOT. calib2)

c Dans le cas d'une calibration du parametre 'para_rox1' : --------------------
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
	    write(7,'(2d13.6)') l1, age
         close(unit=7)
c        On interpole pour obtenir la masse de la limite de Roxburgh
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
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
            l2=(s1+A1)/(s2+A2)
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

      subroutine affine( s1, s2, x1, x2, y11, y21, y12, y22, x0 )

c But : on cherche la limite de Roxburgh entre deux abscisses x1 et x2 (x1 < x2)
c       de facon analytique.

c DC, Octobre 1999

c cf. page 36 cahier n. XII

c Entrees :

c s1 et s2 sont les integrales de Roxburgh (respectivement membre de gauche et membre
c de droite jusqu'a indice 'irox-1' (cad valeur d'abscisse 'x1')

c y11 : valeur de l'INTEGRANDE de l'integrale 1 (mbre de gauche) en x1
c y21 : valeur de l'INTEGRANDE de l'integrale 1 (mbre de gauche) en x2

c y12 : valeur de l'INTEGRANDE de l'integrale 2 (mbre de droite) en x1
c y22 : valeur de l'INTEGRANDE de l'integrale 2 (mbre de droite) en x2

c Sorties :

c x0  : valeur de la limite de Roxburgh.

      implicit none

      logical in1, in2

      real*8 x1, x2, y11, y21, y12, y22, a1, b1, a2, b2, alpha, beta, 
     +       gamma, delta, s1, s2, x01, x02, x0

      if ( x1 .eq. x2 ) then
         print*, ' '
         print*, 'Pb. dans affine : x1=x2 !!'
         stop
      end if

      if ( x1 .gt. x2 ) then
         print*, ' '
         print*, 'Pb. dans affine : x1 > x2 !!'
         stop
      end if

      in1=.true. ! si in1=T alors x01 est dans l'intervalle x1, x2

      in2=.true. ! si in2=T alors x02 est dans l'intervalle x1, x2

      a1=(y21-y11)/(x2-x1)
      b1=y11-a1*x1

      a2=(y22-y12)/(x2-x1)
      b2=y12-a2*x1

      alpha=(a1-a2)/2.

      beta=b1-b2

      gamma=s1-s2+(a2-a1)*x1**2/2.+(b2-b1)*x1

      delta=beta**2-4.*alpha*gamma

      if ( delta .lt. 0.d0 ) then
         print*, ' '
         print*, 'Pb. dans affine : delta < 0 !!'
         print*, 'delta= ', delta
         stop
      end if

      if ( alpha .eq. 0.d0 ) then
         if ( beta .ne. 0.d0 ) then
            print*, ' '
            print*, 'alpha=0 dans affine !!'
            pause
            x0=-gamma/beta
            return
         else
            print*, ' '
            print*, 'Pas de solution dans affine !!'
            print*, 'alpha=0 et beta=0 !!!'
            stop
         end if
      end if

c les solutions du polynome du deuxieme degre :

      x01=(-beta-sqrt(delta))/2./alpha

      x02=(-beta-sqrt(delta))/2./alpha

c on teste l'appartenance de x01 et x02 a l'intervalle x1, x2

      if ( ( x01 .gt. x2 ) .OR. ( x01 .lt. x1 ) ) then
         in1=.false.
      end if

      if ( ( x02 .gt. x2 ) .OR. ( x02 .lt. x2 ) ) then
         in2=.false.
      end if

c Cas ou aucune solution n'est dans l'intervalle :

      if ( (.NOT. in1) .AND. (.NOT. in2) ) then
         print*, ' '
         print*, 'Pb. dans affine : aucune sol. dans l''intervalle x1, x
     +2!'
         print*, 'x1= ', x1
         print*, 'x2= ', x2
         print*, 'x01= ', x01
         print*, 'x02= ', x02

         stop
       end if

c Cas ou les deux solutions sont ensemble dans x1, x2 :

       if ( (in1 .AND. in2) .AND. (x01 .ne. x02 ) ) then
         print*, ' '
         print*, 'Pb. dans affine : les deux  sol. sont dans l''interval
     +le x1, x2!'
         print*, 'x1= ', x1
         print*, 'x2= ', x2
         print*, 'x01= ', x01
         print*, 'x02= ', x02
         stop
       end if

c Si pas de probleme :

       if ( in1 ) then
          x0=x01
          return
       end if

       if ( in2 ) then
          x0=x02
          return
       end if

       end
      


      

c*********************************************************************

      subroutine rech_dich_lim_rox( introx1, introx2, n, irox, r, 
     +                              r_rox )

c Affinement par dichotomie de la valeur de la masse pour laquelle
c on a egalite des integrales de Roxburgh.

c Entrees :
c           * introx1 : integrande du membre de gauche
c           * introx2 :  "          "  "      " droite
c           * irox    : premier indice ou Integrale 1 < Integrale 2
c           * r       : valeurs du rayon (en rsol)
c           * n       : nombre de couches dans le modele

c Sortie  :
c           * r_roc   : le rayon de Roxburgh en "rsol"
 
      implicit none

      include 'cesam_3.parametres'

      integer irox, n, ntour, ntour_max

      real*8 introx1(pn), introx2(pn), r(pn), srox1, srox2

      real*8 x10,x20, x1, x2, yrox1_1, yrox1_2,yrox2_1,
     +  yrox2_2, s1, s2, x0, r_rox

      real*8 precis, precis_mini

      parameter ( precis_mini= 1.d-3, ntour_max= 100 )

      call integrale(introx1,r,n,irox-1,srox1) 

      call integrale(introx2,r,n,irox-1,srox2)

      srox1=srox1/1.d-32
      srox2=srox2/1.d-32

      x10=r(irox-1)
      x20=r(irox)

      print*, ' '
      print*, 'r(irox-1)= ', r(irox-1)
      print*, 'r(irox  )= ',r(irox  )
      print*, ' '

      x1=x10
      x2=x20

      yrox1_1=introx1(irox-1)/1.d-32
      yrox1_2=introx1(irox)/1.d-32

      yrox2_1=introx2(irox-1)/1.d-32
      yrox2_2=introx2(irox)/1.d-32

      ntour=0

      precis=100.

      do while ( precis .gt. precis_mini )

         x0=(x1+x2)/2.

         call integrale2(x10,x20,yrox1_1,yrox1_2,x0,s1)
         s1=s1+srox1
         call integrale2(x10,x20,yrox2_1,yrox2_2,x0,s2)
         s2=s2+srox2

         if ( s1 .ge. s2 ) then
            x1=x0
         else
            x2=x0
         end if

         print*, 'x0= ', x0

         precis=abs(s1-s2/s1)
         print*, 'precis= ', precis

         ntour=ntour+1
         print*, 'ntour= ', ntour

         if ( ntour .gt. ntour_max ) then
            print*, ' '
            print*, 'Non convergence dans ''rech_dich_lim_rox'''
            stop
         end if
      end do

      r_rox=x0

      return

      end

c*********************************************************************

      subroutine integrale( f, x, nmax, ns, sum )

c Integration de 'f' jusqu'a l'abscisse d'indice 'ns'

c Auteur : D. Cordier, Fevrier 1999

      implicit none

      include 'cesam_3.parametres'

      integer nmax, ns, i

      real*8 f(pn), x(pn), sum, dx, fact

      sum=0.d0

      do i= 2, ns
         dx=x(i)-x(i-1)
         fact=(f(i)+f(i-1))/2.
         sum=sum+fact*dx
      end do

      return

      end

c*********************************************************************

      subroutine integrale2( x1, x2, y1, y2, x0, s )

c Calcul de "l'integrale" entre x1 et x0 pour une droite passant
c par les points (x1,y1) et (x2,y2).

      implicit none

      real*8 x1, x2, y1, y2, x0, s, a, b, y0, dx,  moy

      if ( x1 .eq. x2 ) then
         print*, ' '
         print*, 'Pb. dans ''integrale2'' (int_roxburgh) :'
         print*, 'x1=x2 !'
         print*, ' '
         stop
      end if

      a=(y2-y1)/(x2-x1)
      b=y1-a*x1

      y0=a*x0+b

      dx=x0-x1

      moy=(y1+y0)/2.

      s=moy*dx

      return

      end

c**********************************************************************

      subroutine ouverture_fichiers

      implicit none

c Les unites logiques deja utilisees par CESAM 3.2 (en Mai 1999) :
c 1, 2, 3, 4, 11 (opa), 24, 25, 26, 30, 31, 48, 53, 60, 66, 71, 76

         open(unit=41,status='unknown',form='formatted',
     +        file='loglrad.roxdat')

         open(unit=42,status='unknown',form='formatted',
     +        file='loglnuc.roxdat')

         open(unit=43,status='unknown',form='formatted',
     +        file='delad.roxdat')

         open(unit=44,status='unknown',form='formatted',
     +        file='delrad.roxdat')

         open(unit=45,status='unknown',form='formatted',
     +        file='kappa.roxdat')

         open(unit=46,status='unknown',form='formatted',
     +        file='x.roxdat')

         open(unit=47,status='unknown',form='formatted',
     +        file='Vconv.roxdat')

c 48 est deja pris
         open(unit=49,status='unknown',form='formatted',
     +        file='Fconv.roxdat')

         open(unit=50,status='unknown',form='formatted',
     +        file='Hp.roxdat')

         open(unit=51,status='unknown',form='formatted',
     +        file='Int_Rox1.roxdat')

         open(unit=52,status='unknown',form='formatted',
     +        file='Int_Rox2.roxdat')

         open(unit=54,status='unknown',form='formatted',
     +        file='srox1.roxdat')

         open(unit=55,status='unknown',form='formatted',
     +        file='srox2.roxdat')

         print*, ' '
         print*, 'Les fichiers du test Roxburgh ont l''extension :'
         print*, ' '
         print*, '            ''.roxdat''                        '
         print*, ' '

         print*, 'les grandeurs listees sont :'
         print*, ' '
         print*, 'Log Lrad, Log Lnuc, delad, delrad, kappa, X'
         print*, 'Fconv, Vconv, Int_Rox1 et Int_Rox2'

         return

         end

c*************************************************************************

      subroutine ecriture( mass, renrsol, loglrad, loglnuc, delad, 
     +       delrad, kappa, x, Vconv, Fconv, Hp, Int_rox1, Int_rox2 )

c Ecriture de certaines grandeurs lors du test de la routine

c Avril 99

      implicit none

      logical enmasse

      parameter ( enmasse= .false. )

      real*8 mass, loglrad, loglnuc, delad, delrad, kappa, x, Vconv,
     +       Fconv, Hp, Int_rox1, Int_rox2, renrsol

      if ( enmasse ) then
         write(41,1000) mass, loglrad
         write(42,1000) mass, loglnuc
         write(43,1000) mass, delad
         write(44,1000) mass, delrad
         write(45,1000) mass, kappa
         write(46,1000) mass, x
         write(47,1000) mass, Vconv

         write(49,1000) mass, Fconv
         write(50,1000) mass, Hp
         write(51,1000) mass, Int_rox1
         write(52,1000) mass, Int_rox1
      else
         write(41,1000) renrsol, loglrad
         write(42,1000) renrsol, loglnuc
         write(43,1000) renrsol, delad
         write(44,1000) renrsol, delrad
         write(45,1000) renrsol, kappa
         write(46,1000) renrsol, x
         write(47,1000) renrsol, Vconv

         write(49,1000) renrsol, Fconv
         write(50,1000) renrsol, Hp
         write(51,1000) renrsol, Int_rox1
         write(52,1000) renrsol, Int_rox1
      end if

 1000 format(1p2d13.6)

      return

      end

c**********************************************************************

      subroutine fermeture

c fermeture des fichiers ouverts pour le test de la routine.

      implicit none

      close(41)
      close(42)
      close(43)
      close(44)
      close(45)
      close(46)
      close(47)

      close(49)
      close(50)
      close(51)
      close(52)

      close(54)
      close(55)

      return

      end



















