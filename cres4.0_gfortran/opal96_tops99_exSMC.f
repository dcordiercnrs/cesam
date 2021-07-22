c--------------------------------------------------------------------------
 
      subroutine opal96_tops99_exSMC(xh,t,ro,kappa,dkapdt,dkapdr,dkapdx)
 

c Routine d'opacite avec opacite TOPS version 4.0, uniquement pour
c etoiles avec Z0 = 0.02, 0.008 et 0.004

c ====================================================================
c ====================================================================
c ATTENTION : cette version permet d'extrapoler pour des valeurs de Z0
c             proches de celle du SMC !
c ====================================================================
c ====================================================================
 
c Modif. Decembre 1998 : on stoppe lorsque on depasse les domaines de
c validite des tables

c Janvier 1999 : utilisation de nouvelles tables (etendues par rapport
c aux anciennes)

c Fevrier 99 : on rebatise 'quad' pour eviter des conflits avec 
c 'quad' de la source CESAM
 
c Daniel, Mai 1998, Mai 2000.
 
c
c      calcul de l'opacite repris du code de Geneve
c      polynomes de Lagrange pour log T6 et log R
c      quadratique pour X et Z
c       Yveline (aout 1991-code de Geneve--> 10/95 CESAM)
c
 
c Remarques : * l'opacite conductive est incluse!!!
c             * les routines supplementaires ont ete mises a la
c               fin du fichier.
 
c-------------------------------------------------------------------
c      entree :
c      xh(1)=X : comp. chim. en H
c      t : temperature K
c      ro : densite cgs
 
c      sortie :
c      kappa : opacite gr / cm2)
c      dkapdt : kappa / d t
c      dkapdr : kappa / d densite
c      dkapdx : kappa / d xchim(1)
 
c      Z est obtenu par 1-X-Y
c-------------------------------------------------------------------
 
      implicit none
 
      include 'cesam_3.parametres'
      include 'modele_3.common'
      include 'evol_chim_3.common'
 
      logical tops, init, MW, LMC, SMC, extrapol_tops_smc
 
      real*8 t6, lr, xh(1), t, ro, kappa, dkapdt, dkapdr,
     +       dkapdx, xi, zi, y, dc12, do16, xc12_0, xo16_0,
     +       log_do16sdc12, ymax, delta, bid1, bid2, bid3,
     +       kap_opal, kap_tops, limite, limite_smc
 
      real*8 a, b, xc12_0_mw, xc12_0_lmc, xc12_0_smc,
     +             xo16_0_mw, xo16_0_lmc, xo16_0_smc,
     +       ktlmc, ktsmc

      character*50 file_opa_tops, file_opa_tops1, file_opa_tops2
 
      parameter( limite= 1.d-10, limite_smc= 0.003001d0  )
 
      data init/.true./

      data MW/.false./
      data LMC/.false./
      data SMC/.false./

      data xc12_0_mw, xc12_0_lmc, xc12_0_smc/0.34664719d-02, 
     +                        0.13865888d-02,0.69329439d-03/

      data xo16_0_mw, xo16_0_lmc, xo16_0_smc/0.96438571d-02,
     +                        0.38575428d-02,0.19287714d-02/

c     ------ calcul des valeurs pour lesquelles on veut l'opacite ---
 
      t6=t*1.d-6			!T6
      lr=log10(ro/t6**3)		!log R
      xi=xh(1)		!determination de X et de Z
 
      y=xh(ihe4)+xh(ihe4-1)
 
      zi=1.d0-xi-xh(ihe4)-xh(ihe4-1)

c Quelques initialisations suplementaires :

      extrapol_tops_smc=.false.

      if( abs(z0-0.004d0).le. limite_smc ) then
         if ( init .AND. abs(z0-0.004d0) .gt. limite ) then
         print*, ' '
         print*, '*************************************************'
         print*, 'ATTENTION : '
         print*, '   * Valeur de Z0= ', z0
         print*, '     non-compatible avec les tables d''opacite!'
         print*, '     <<TOPS>> pour le SMC'
         print*, ' '
         print*, '   ===>  on va donc donc extrapoler !'
         print*, ' '
         print*, '*************************************************'
         print*, ' '
         print*, 'Ok ?'
         pause
         end if
      else
         if( abs(z0-0.02d0) .ge. limite .AND.
     +        abs(z0-0.008d0).ge. limite ) then
            print*, ' '
            print*, '*************************************************'
            print*, 'Valeur de Z0= ', z0
            print*, 'non compatible avec les tables d''opacite <<TOPS>>!
     +'
            print*, ' '
            print*, 'Il faut construire une routine qui extrapole'
            print*, 'Comme cela a ete fait pour les metallicites'
            print*, 'proche de celle du SMC avec ''opal96_tops99_exSMC''
     +'
            print*, '*************************************************'
            print*, ' '
            stop
         end if
      end if

      if( abs(z0-0.02d0) .le. limite ) then
        xc12_0=xc12_0_mw
        xo16_0=xo16_0_mw
        dc12=xh(4)-xc12_0
        do16=xh(9)-xo16_0
        file_opa_tops=opa2
        if(init)then
          print*, 'OK zO= 0.02'
          print*, 'file_opa= ', file_opa_tops
          MW=.true.
          print*, 'MW= ', MW, ' , LMC= ', LMC, ' , SMC= ', SMC
          if ( rep_pause .ne. 'n') then
          pause
          end if
          init=.false.
        end if
      end if
 
      if( abs(z0-0.008d0) .le. limite ) then
        xc12_0=0.13865888d-02
        xo16_0=0.38575428d-02
        dc12=xh(4)-xc12_0
        do16=xh(9)-xo16_0
        file_opa_tops=opa3
        if(init)then
          print*, 'OK zO= 0.008'
          print*, 'file_opa= ', file_opa_tops
          LMC=.true.
          print*, 'MW= ', MW, ' , LMC= ', LMC, ' , SMC= ', SMC
          if ( rep_pause .ne. 'n' ) then
          pause
          end if
          init=.false.
        end if
      end if
 
      if( abs(z0-0.004d0) .le. limite ) then
        xc12_0=0.69329439d-03
        xo16_0=0.19287714d-02
        dc12=xh(4)-xc12_0
        do16=xh(9)-xo16_0
        file_opa_tops=opa4
        if(init)then
          print*, 'OK z0= 0.004'
          print*, 'file_opa= ', file_opa_tops
          SMC=.true.
          print*, 'MW= ', MW, ' , LMC= ', LMC, ' , SMC= ', SMC
          if ( rep_pause .ne. 'n' ) then
          pause
          end if
          init=.false.
        end if
      end if

      if ( abs(z0-0.004d0) .gt. limite .AND.
     +     abs(z0-0.004d0) .le. limite_smc ) then ! Cas ou on interpole pour le SMC
         extrapol_tops_smc=.true.

c         a=(xc12_0_lmc-xc12_0_smc)/(0.008-0.004)
c         b=xc12_0_smc-a*0.004
c         xc12_0=a*z0+b
         xc12_0=z0/0.004d0*xc12_0_smc
c         a=(xo16_0_lmc-xo16_0_smc)/(0.008-0.004)
c         b=xo16_0_smc-a*0.004
c         xo16_0=a*z0+b
         xo16_0=z0/0.004d0*xo16_0_smc

         dc12=xh(4)-xc12_0
         do16=xh(9)-xo16_0

         file_opa_tops1=opa3 ! LMC
         file_opa_tops2=opa4 ! SMC

c        if ( dc12 .lt. 0. ) then
c           print*, 'dc12 .lt. 0.'
c           print*, 'dc12= ', dc12
c           pause
c        end if
c        if ( do16 .lt. 0. ) then
c           print*, 'do16 .lt. 0.'
c           print*, 'do16= ', do16
c           pause
c        end if

        if(init)then
          print*, 'OK z0= ', z0
          print*, 'les fichiers "file_opa" : '
          print*, file_opa_tops1
          print*, file_opa_tops2
          print*, 'dc12= ', dc12
          print*, 'do16= ', do16
          print*, 'xh(4)= ', xh(4)
          print*, 'xh(9)= ', xh(9)
          if ( rep_pause .ne. 'n' ) then
             pause
          end if
          init=.false.
        end if

      end if

         
      if( xi .le. 1.d-3 .AND. zi .gt. z0 ) then
c        if ( zi .gt. 0.1 ) then
c           print*, 'y   = ', y
c           print*, 'zi  = ', zi
c           print*, 'xi   = ', xi
c           print*, 'T    = ', t
c           print*, 'ro   = ', ro
c           print*, 'dc12 = ', dc12
c           print*, 'do16 = ', do16
c           pause
c       end if
        if( dc12 .gt. 0.d0 .AND. do16 .gt. 0.d0 ) then
            log_do16sdc12=log10(do16/dc12)
        else
            if ( dc12 .gt. 0.d0 .AND. do16 .le. 0.d0 ) then
               if ( SMC ) then
            log_do16sdc12=-10.d0
               else
            log_do16sdc12=-4.d0
               end if
c            else
c            print*, 'Pb. dans ''opal96_tops99'' : '
c            print*, 'y    = ', y
c            print*, 'zi   = ', zi
c            print*, 'xi   = ', xi
c            print*, 'T    = ', t
c            print*, 'ro   = ', ro
c            print*, 'dc12 = ', dc12
c            print*, 'do16 = ', do16
c           stop
            end if
            if ( dc12 .le. 0.d0 .AND. do16 .le. 0.d0 ) then
               if ( SMC ) then
            log_do16sdc12=-10.d0 ! peu importe la valeur ici !
               else
            log_do16sdc12=-4.d0
               end if
            end if
            if ( dc12 .le. 0.d0 .AND. do16 .gt. 0.d0 ) then
               if ( SMC ) then
            log_do16sdc12=+5.d0
               else
            log_do16sdc12=+1.d0
               end if
            end if
        end if
        tops=.true.
      else
        tops=.false.
      end if

c Cas ou on N'extrapole PAS TOPS pour une metallicite voisine de celle
c du SMC
      if ( .NOT. extrapol_tops_smc ) then
      if( tops ) then
c        On calcule d'abord l'opacite avec Y=Ymax=1.-Z0 avec OPAL96
c        et TOPS
         xi=0.d0
         call kappa_opal_3(lr,t6,ro,xi,z0,kap_opal,bid1,bid2,
     +                     bid3,opa1)
         ymax=1.d0-z0
         call kappa_tops_99(ymax,log_do16sdc12,t,ro,file_opa_tops,
     +           kap_tops,dkapdt,dkapdr)
c        On calcule la difference "delta" (normalement ces deux opacites
c        devrait coincider mais il reste une tres petite difference qu'on
c        va retirer a l'opacite finale par soucis de coherence entre TOPS
c        et OPAL
         delta=kap_tops-kap_opal

         call kappa_tops_99(y,log_do16sdc12,t,ro,file_opa_tops,
     +           kappa,dkapdt,dkapdr)

         kappa=kappa-delta

         if( kappa .le. 0.d0 ) then
             print*, 'Pb. avec kappa .le. 0. !'
             stop
         end if
         dkapdx=0.d0
         call kappa_cond_3(xh,t,ro,kappa,dkapdt,dkapdr,dkapdx)
         return
      else
         if( zi .ge. 0.025d0 ) then
           print*, 'Pb. OPAL avec zi .ge. 0.025! Ici : zi= ', zi
           stop
         end if
         call kappa_opal_3(lr,t6,ro,xi,zi,kappa,dkapdr,dkapdt,
     +                     dkapdx,opa1)
         call kappa_cond_3(xh,t,ro,kappa,dkapdt,dkapdr,dkapdx)
c         print*, 'kappa= ', kappa
c         return
      end if
      end if ! .NOT. extrapol_tops_smc

c Cas ou on EXTRAPOLE <<TOPS>> pour une metallicite voisine de celle
c du SMC
      if ( extrapol_tops_smc ) then
      if( tops ) then
c        On calcule d'abord l'opacite avec Y=Ymax=1.-Z0 avec OPAL96
c        et TOPS
         xi=0.d0
         call kappa_opal_3(lr,t6,ro,xi,z0,kap_opal,bid1,bid2,
     +                     bid3,opa1)
         ymax=1.d0-z0
         call kappa_tops_99(ymax,log_do16sdc12,t,ro,file_opa_tops1,
     +           ktlmc,dkapdt,dkapdr) ! LMC
         call kappa_tops_99(ymax,log_do16sdc12,t,ro,file_opa_tops2,
     +           ktsmc,dkapdt,dkapdr) ! SMC
         a=(ktlmc-ktsmc)/(0.008-0.004)
         b=ktsmc-a*0.004
         kap_tops=a*z0+b
c        On calcule la difference "delta" (normalement ces deux opacites
c        devrait coincider mais il reste une tres petite difference qu'on
c        va retirer a l'opacite finale par soucis de coherence entre TOPS
c        et OPAL
         delta=kap_tops-kap_opal

         call kappa_tops_99(y,log_do16sdc12,t,ro,file_opa_tops1,
     +           ktlmc,dkapdt,dkapdr) ! LMC
         call kappa_tops_99(y,log_do16sdc12,t,ro,file_opa_tops2,
     +           ktsmc,dkapdt,dkapdr) ! SMC

         a=(ktlmc-ktsmc)/(0.008-0.004)
         b=ktsmc-a*0.004

         kappa=a*z0+b

         kappa=kappa-delta

         if( kappa .le. 0.d0 ) then
             print*, 'Pb. avec kappa .le. 0. !'
             stop
         end if
         dkapdx=0.d0
         call kappa_cond_3(xh,t,ro,kappa,dkapdt,dkapdr,dkapdx)
         return
      else
         if( zi .ge. 0.025d0 ) then
           print*, 'Pb. OPAL avec zi .ge. 0.025! Ici : zi= ', zi
           stop
         end if
         call kappa_opal_3(lr,t6,ro,xi,zi,kappa,dkapdr,dkapdt,
     +                     dkapdx,opa1)
         call kappa_cond_3(xh,t,ro,kappa,dkapdt,dkapdr,dkapdx)
c         print*, 'kappaOPAL96= ', kappa
c         pause
         return
      end if
      end if ! extrapol_tops_smc

      end
