c**************************************************************************

c                           profil_tops4.f

c**************************************************************************

c Programme de construction de profiles avec les opacites TOSP 4.0

c Daniel Cordier
c 26 Mai 1998

      implicit none

      integer np, nl_max, n_c, j

      real*8 ro,logt_min,logt_max,y,logt,t,pas_logt,log_do16sdc12,
     +       kappa, bid1, bid2

      character*50 filename1, filename2, file_opa

      parameter( logt_min=6.07E0, logt_max=8.06E0, np=60,
     +     file_opa='../data/opa_tops_0.02.dat',nl_max=5000)

      real*8 m_mstar1(nl_max), t1(nl_max), ro1(nl_max)

c      print*, 'Valeur de ro = '
c      read(5,*) ro

      print*, 'Valeur de Y= '
      read(5,*) y

      print*, 'Valeur de log_do16sdc12= '
      read(5,*) log_do16sdc12

      print*, 'Fichier de profile :'
      read(5,'(a)') filename1

      print*, 'Nom du fichier de sortie :'
      read(5,'(a)') filename2


      open(1,status='old',form='formatted',file=filename1)

          n_c=0
          do j= 1, nl_max
             read(1,1000,end=100) m_mstar1(j), t1(j), ro1(j)
c             print*, 'm_mstar= ', m_mstar1(j)
             n_c=n_c+1
c             print*, 'n_c= ', n_c
          end do
1000      format(1p3E19.11)
100   print*, 'Fichier lu!'
      close(1)
        pas_logt=(logt_max-logt_min)/(np-1)
        logt=logt_min

      open(1,status='unkown',form='formatted',file=filename2)

        do j= 1, n_c
           t=10.**logt
           t=t1(j)
           ro=ro1(j)
           if( ro .ge. 100.d0 .AND. t .ge. 10.**logt_min) then
        call kappa_tops_4(y,log_do16sdc12,t,ro,file_opa,kappa,
     +                   bid1,bid2)
           write(1,2000) m_mstar1(j), kappa
2000       format(2E15.4)
           end if
           logt=logt+pas_logt
        end do

       close(1)

       end
