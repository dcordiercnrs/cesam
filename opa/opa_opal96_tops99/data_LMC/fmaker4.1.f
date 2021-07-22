c****************************************************************************

c                             fmaker4.1.f


c   Programme de conversion des tables d'opacites "TOPS" en fichier binaire

c   VERSION 4.0 : * nouveau melange (GN'93), nouveaux domaines de temperature
c                   et de densite ; nouvelles variables : ro, T, X_He4,
c                   log(Delta_O16/Delta_C12), log(Delta_Ne22/Delta_C12), Z0


c****************************************************************************


c Daniel Cordier, Mai 1998, avril 2001

c Janvier 1999 : petites modifications, correction de petits bugs

c ATTENTION : ce programme puise des arguments dans un fichier
c             de type <<filename002.dat>>

c ATTENTION BIS : les data doivent etre ecrites avec un <<d>> devant
c                 l'exposant!!!

      implicit none

      integer jt, jr, jf, j_y, j_o16

      integer nt, nr, nf, n_y, n_o16, nf_max, ny_max, no16_max


      parameter( nt=26, nr=30, nf_max=50, ny_max=10, no16_max=20 )

c nf_max = nombre maximum de fichiers de sortie de TOPS
c ny_max et no16_max= nombre maximum de valeurs pour :
c  Y=X_he4 et log(DO16/DC12), ici on a log(DNe22/DC12)=-1.4 

      real*8 tkev(nt), t6(nt), ro(nr), vlro(nr), eve, kbol, z0,
     +       k(nf_max,nt,nr), vlt6(nt), bid, val_y(ny_max),
     +       val_o16(no16_max)

      character filename1*26, filename2*26, filename3*26, ligne*75


      data eve / 1.60217733d-19 /, kbol/ 1.380658d-23 /

      data tkev /1.0000d-01,1.2500d-01,1.5000d-01,2.0000d-01,
     +           2.5000d-01,3.0000d-01,4.0000d-01,5.0000d-01,
     +           6.0000d-01,
     +           8.0000d-01,1.0000d+00,1.2500d+00,1.5000d+00,
     +           2.0000d+00,2.5000d+00,3.0000d+00,4.0000d+00,
     +           5.0000d+00,6.0000d+00,8.0000d+00,1.0000d+01,
     +           1.5000d+01,2.5000d+01,4.0000d+01,6.0000d+01,
     +           1.0000d+02/

      data ro /1.0000d+02,1.3180d+02,1.7370d+02,2.2893d+02,3.0172d+02,
     +         3.9765d+02,
     +         5.2409d+02,6.9072d+02,9.1034d+02,1.1998d+03,1.5813d+03,
     +         2.0840d+03,
     +         2.7467d+03,3.6200d+03,4.7710d+03,6.2880d+03,8.2873d+03,
     +         1.0922d+04,
     +         1.4395d+04,1.8972d+04,2.5004d+04,3.2955d+04,4.3433d+04,
     +         5.7242d+04,
     +         7.5443d+04,9.9430d+04,1.3104d+05,1.7271d+05,2.2763d+05,
     +         3.0000d+05/

c Initialisation des parametres de composition chimique

      print*, ' '
      print*, '***************************************'
      print*, ' '
      print*, ' Tabulation des opacités pour stades'
      print*, '       avancés (Los Alamos)'
      print*, ' '
      print*, '***************************************'
 10   print*, ' '
      print*, ' - Nom du fichier de composition chimique : '
      print*, '   (Du type "compo_LMC.dat")'
      read(5,'(a)',err=10) filename1

      open(1,status='old',form='formatted',file=filename1)

           read(1,1000) nf
           read(1,1000) n_y
           read(1,1000) n_o16

      if ( nf .gt. nf_max ) then
         print*, 'ATTENTION : nf .gt. nf_max !'
         print*, 'Redimensionner le parametre.'
         stop
      end if

      if ( n_y .gt. ny_max ) then
         print*, 'ATTENTION : n_y .gt. ny_max !'
         print*, 'Redimensionner le parametre.'
         stop
      end if

      if ( n_o16 .gt. no16_max ) then
         print*, 'ATTENTION : n_o16 .gt. no16_max !'
         print*, 'Redimensionner le parametre.'
         stop
      end if


           read(1,2000) ligne

           read(1,3000) z0

           read(1,2000) ligne

      
           do j_y= 1, n_y
              read(1,3000) val_y(j_y)
           end do

           read(1,2000) ligne

           do j_o16= 1, n_o16
              read(1,3000) val_o16(j_o16)
           end do

1000  format(I4)
2000  format(a)
3000  format(D11.3)

      close(1)

c Conversion des valeurs des temperatures keV ---> MK

      do jt= 1, nt
         t6(jt)=tkev(jt)*1.d+03*eve/kbol/1.d+06
      end do

c Test
c      print*, 't6= ', t6
c      stop

c Conversion de t6 ---> log10 t6

      do jt= 1, nt
         vlt6(jt)=log10(t6(jt))
      end do
c Test
c      print*, 'vlt6= ', vlt6
c      stop

c Conversion de ro ---> log10 ro

      do jr= 1, nr
         vlro(jr)=log10(ro(jr))
      end do

c-------------------------------------------

 20   print*, ' '
      print*, ' - Nom du fichier contenant les noms de fichiers :'
      print*, '   (du type "filename_LMC.dat")'
      read(5,'(a)',err=20) filename2

 30   print*, ' '
      print*, ' - Nom du fichier de sortie :'
      print*, '   (fichier qui sera lu par CESAM)'
      read(5,'(a)',err=30) filename3


      open( 1, status='old', form='formatted', file=filename2)
      print*, 'ATTENTION verifier que les fichiers indiques dans'
      print*, filename2, ' y sont inscrits dans le bon'
      print*, 'ordre!'
      print*, ' '
      open( 2, status='unknown', form='unformatted', file=filename3 )

      write(2) nf, n_y, n_o16, nt, nr
      write(2) z0
      write(2) (val_y(j_y), j_y= 1, n_y)
      write(2) (val_o16(j_o16), j_o16= 1, n_o16)
      write(2) (vlt6(jt), jt= 1, nt)
      write(2) (vlro(jr), jr= 1, nr)

      jf=1

      do j_y= 1, n_y

         do j_o16= 1, n_o16

         read(1,'(a)') filename1

         print*, ' '
         print*, 'Nom du fichier : ', filename1
         print*, 'Composition : '
         print*, 'He4 : ', val_y(j_y), ' log(DO16/DC12) : ',
     +           val_o16(j_o16), ' Z0 : ', z0

         open(3, status='old', file=filename1)

         do jt= 1, nt
            read(3,2000) ligne
            do jr= 1, nr
               read(3,4000) bid, k(jf,jt,jr)
            end do
         end do

         close(3)

        jf=jf+1

        end do

      end do

      write(2) (((k(jf,jt,jr),jr=1,nr),jt=1,nt),jf=1,nf)

      close(2)
      close(1)

      print*, ' '
      print*, 'Fichier binaire : ', filename3, ' cree!'
      print*, ' '

4000  format(2E12.4)

      end
