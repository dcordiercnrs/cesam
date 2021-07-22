c--------------------------------------------------------------------------
 
      subroutine opal96_tops99_2000(xh,t,ro,kappa,dkapdt,dkapdr,dkapdx)
 

c Routine d'opacite avec opacites OPAL96 (Yveline) et opacites de Los Alamos
c pour Z élevée (serveur <<TOPS>>, Daniel) et interpolation/extrapolation (de TOPS) 
c lorsque Z0 ne vaut pas 0.02, 0.008 ou 0.004

c Remarque IMPORTANTE : les noms des fichiers contenant les tables doivent etre
c                       rentres dans l'ordre suivant dans le fichier '*.don' :
c    OPA1 = table OPAL
c    OPA2 = table TOPS pour Z0=0.02  (MW)
c    OPA3 = table TOPS pour Z0=0.008 (SMC)
c    OPA4 = table TOPS pour Z0=0.004 (LMC)

c Daniel, Mai 2000.
 
c
c      calcul de l'opacite repris du code de Geneve
c      polynomes de Lagrange pour log T6 et log R
c      quadratique pour X et Z
c       Yveline (aout 1991-code de Geneve--> 10/95 CESAM)
c
 
c Remarques : * l'opacite conductive est incluse!!!

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
 
      logical tops, init, MW, LMC, SMC
 
      real*8 t6, lr, xh(1), t, ro, kappa, dkapdt, dkapdr,
     +       dkapdx, xi, zi, y, dc12, do16, xc12_0, xo16_0,
     +       log_do16sdc12, ymax, delta, bid1, bid2, bid3,
     +       kap_opal, kap_tops, lim_smc1, lim_smc2,
     +       lim_lmc1, lim_lmc2, lim_mw1, lim_mw2 
 
      real*8 a, b, xc12_0_mw, xc12_0_lmc, xc12_0_smc,
     +             xo16_0_mw, xo16_0_lmc, xo16_0_smc,
     +       kt1, kt2, z01, z02

      character*50 file_tops1, file_tops2

c On définit ici les intervalles de métallicité dans lesquels
c on va utiliser les différentes tables de Los Alamos :
      parameter( lim_smc1= 0.0009d0, lim_smc2= 0.006d0,
     +           lim_lmc1= 0.006d0 , lim_lmc2= 0.010d0,
     +           lim_mw1 = 0.010d0 , lim_mw2 = 0.035d0  )
 
      data init /.true./

      data MW   /.false./
      data LMC  /.false./
      data SMC  /.false./

c Estimations des fractions massiques initiales en C12 et O16 pour MW,
c LMC et SMC :
      data xc12_0_mw, xc12_0_lmc, xc12_0_smc/0.34664719d-02, 
     +                        0.13865888d-02,0.69329439d-03/

      data xo16_0_mw, xo16_0_lmc, xo16_0_smc/0.96438571d-02,
     +                        0.38575428d-02,0.19287714d-02/

c----------------------------------------------------------------------------
c On calcule des valeurs pour lesquelles on veut l'opacite
 
      t6=t*1.d-6			!T6
      lr=log10(ro/t6**3)		!log R
      xi=xh(1)		!determination de X et de Z
 
      y=xh(ihe4)+xh(ihe4-1)
 
      zi=1.d0-xi-xh(ihe4)-xh(ihe4-1)

c----------------------------------------------------------------------------
c On cherche les tables a utiliser pour <<TOPS>> :

      if( init ) then
         print*, 'Routine avec extrapol/interpol en Z0'
         print*, 'il vaut mieux la tester avant utilisation !'
         print*, 'Ok ?'
c         pause
         if ( ( lim_smc1 .le. z0 ) .AND.
     +        ( z0 .lt. lim_smc2 )      ) then ! Cas du SMC
            file_tops1=opa4 ! SMC
            z01=0.004d0
            file_tops2=opa3 ! LMC
            z02=0.008d0

            xc12_0=z0/z01*xc12_0_smc
            xo16_0=z0/z01*xo16_0_smc

            dc12=xh(4)-xc12_0
            do16=xh(9)-xo16_0

            SMC=.true. ! On a une metallicite proche de celle du SMC
         end if

         if ( ( lim_lmc1 .le. z0 ) .AND.
     +        ( z0 .lt. lim_lmc2 )      ) then ! Cas du LMC
            if ( z0 .ge. 0.008d0 ) then
               file_tops1=opa3   ! LMC
               z01=0.008d0
               file_tops2=opa2  ! MW
               z02=0.02d0

               xc12_0=z0/z01*xc12_0_lmc
               xo16_0=z0/z01*xo16_0_lmc
            else
               file_tops1=opa4  ! SMC
               z01=0.004d0
               file_tops2=opa3  ! LMC
               z02=0.008d0

               xc12_0=z0/z02*xc12_0_lmc
               xo16_0=z0/z02*xo16_0_lmc
            end if
               dc12=xh(4)-xc12_0
               do16=xh(9)-xo16_0
               LMC=.true.       ! On a une metallicite proche de celle du SMC
               if ( SMC ) then
                  print*, 'Pb. de chevauchement des domaines de Z0 entre 
     +SMC'
                  print*, 'et LMC : ajuster les valeurs de lim_smc et'
                  print*, 'lim_lmc'
               end if
         end if

         if ( ( lim_mw1 .le. z0 ) .AND.
     +        ( z0 .le. lim_mw2 )      ) then ! Cas de la Galaxie
            file_tops1=opa3 ! LMC
            z01=0.008d0
            file_tops2=opa2 ! MW
            z02=0.02d0

            xc12_0=z0/z02*xc12_0_mw
            xo16_0=z0/z02*xo16_0_mw

            dc12=xh(4)-xc12_0
            do16=xh(9)-xo16_0

            MW=.true. ! On a une metallicite proche de celle du SMC
            if ( SMC ) then
                  print*, 'Pb. de chevauchement des domaines de Z0 entre 
     +SMC'
                  print*, 'et MW : ajuster les valeurs de lim_smc et'
                  print*, 'lim_mw'
            end if
            if ( LMC ) then
                  print*, 'Pb. de chevauchement des domaines de Z0 entre 
     +LMC'
                  print*, 'et MW : ajuster les valeurs de lim_lmc et'
                  print*, 'lim_mw'
            end if
         end if

          print*, 'OK z0= ', z0
          print*, 'les tables <<TOPS>> utilisées sont : '
          print*, file_tops1
          print*, file_tops2
          print*, 'dc12= ', dc12
          print*, 'do16= ', do16
          print*, 'xh(4)= ', xh(4)
          print*, 'xh(9)= ', xh(9)
          print*, 'MW= ', MW, ' , LMC= ', LMC, ' , SMC= ', SMC
          if ( rep_pause .ne. 'n' ) then
             pause
          end if
          init=.false.

      end if

c-----------------------------------------------------------------------
c On cherche si on doit utiliser OPAL (tops=.false.) ou TOPS (tops=.true.)
         
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

c-----------------------------------------------------------------------
c On calcule la valeur de l'opacite :

      if( .NOT. tops ) then ! Cas ou on utilise OPAL
         if( zi .ge. 0.025d0 ) then
           print*, 'Pb. OPAL avec zi .ge. 0.025! Ici : zi= ', zi
           stop
         end if
         call kappa_opal_3(lr,t6,ro,xi,zi,kappa,dkapdr,dkapdt,
     +                     dkapdx,opa1)
         call kappa_cond_3(xh,t,ro,kappa,dkapdt,dkapdr,dkapdx)
         return
      else ! Cas ou on utilise TOPS
c        On calcule d'abord l'opacite avec Y=Ymax=1.-Z0 avec OPAL96
c        et TOPS
         xi=0.d0
         call kappa_opal_3(lr,t6,ro,xi,z0,kap_opal,bid1,bid2,
     +                     bid3,opa1)
         ymax=1.d0-z01
         call kappa_tops_99(ymax,log_do16sdc12,t,ro,file_tops1,
     +           kt1,dkapdt,dkapdr) ! Z01
         ymax=1.d0-z02
         call kappa_tops_99(ymax,log_do16sdc12,t,ro,file_tops2,
     +           kt2,dkapdt,dkapdr) ! Z02
c        On interpole/extrapole
         a=(kt2-kt1)/(z02-z01)
         b=kt1-a*z01
         kap_tops=a*z0+b
c        On calcule la difference "delta" (normalement ces deux opacites
c        devraient coincider mais il reste une tres petite difference qu'on
c        va retirer a l'opacite finale par soucis de coherence entre TOPS
c        et OPAL
         delta=kap_tops-kap_opal

c        On calcule Kappa en interpolant/extrapolant si necessaire :
         call kappa_tops_99(y,log_do16sdc12,t,ro,file_tops1,
     +           kt1,dkapdt,dkapdr) ! Z01
         call kappa_tops_99(y,log_do16sdc12,t,ro,file_tops2,
     +           kt2,dkapdt,dkapdr) ! Z02

         a=(kt2-kt1)/(z02-z01)
         b=kt1-a*z01

         kappa=a*z0+b

c         print*, 'Opcaite TOPS = ', kappa
c         pause

c        On retire le petit decalage entre TOPS et OPAL
         kappa=kappa-delta

c        Au cas ou ... ;-)
         if( kappa .le. 0.d0 ) then
             print*, 'Pb. avec kappa .le. 0. !'
             stop
         end if
         dkapdx=0.d0
         call kappa_cond_3(xh,t,ro,kappa,dkapdt,dkapdr,dkapdx)
         return
      end if

      end
