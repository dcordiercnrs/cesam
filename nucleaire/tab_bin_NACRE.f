c****************************************************************************
c
c  
c              Programme de tabulation des taux de reaction 
c                      dans un fichier binaire
c
c                Taux de réaction de la compilation NACRE
c                     (sauf quelques exceptions)
c
c****************************************************************************

c Daniel Cordier Septembre 1997, Aout 2000.
c Version avec NACRE : Avril 2003, Varsovie.

      implicit none

      integer i, i_t, i_reac, i_cycle, pnt, nr_pp, nr_cno, nr_3a,
     &        num_reac, i_reac_max, nt, long, annee_c12ago16 

      real*8 t_tmin, t_tmax, past_t, Q, rt1, rt2, coefftest,
     &       coefftest3a

      logical existe, debut_pp, debut_cno, debut_3a

      character filename*100, cycle*3, le_cycle*3, rep*1

      parameter(pnt=300, t_tmin=1.d-3, t_tmax=1.d+1 )

c Rq. : t_t = T / 1.d+9
c       pnt= nombre de valeurs de la temperature

      parameter (nr_pp=7, nr_cno=13, nr_3a=10)

c nombre de reactions dans PP, CNO et 3 alpha

      real*8 t_t(pnt), r_t(pnt*(nr_pp+nr_cno+nr_3a)), tt_t(pnt),
     &                tr_t(pnt*(nr_pp+nr_cno+nr_3a)),
     &     coeff_c12ago16

      dimension cycle(3), num_reac(nr_pp+nr_cno+nr_3a)

      common /c12ago16/coeff_c12ago16,annee_c12ago16

      data cycle/'pp ', 'cno', '3a '/

      data num_reac/10,20,30,40,60,50,70,   ! reactions chaines PP
     &              130,140,110,160,120,170,190,150,180,200, ! cycle CNO
     &              210,220,230,
     &              1010,1020,1030,1040,1050,1060,1070,1080,1090,
     &              1100/

      external long

      open(unit=88,status='unknown',form='formatted',access='append',
     &     file='rapp_8588.dat')
      close(88,status='delete')

      print*, ' '
      print*, ' **********************************'
      print*, ' '
      print*, '      Fabrication des tables' 
      print*, '               de'
      print*, '    taux de reactions nucleaires'
      print*, ' '
      print*, ' **********************************'
      print*, ' '
      print*, 'Nombre de reactions pour PP  : ', nr_pp
      print*, 'Nombre de reactions pour CNO : ', nr_cno
      print*, 'Nombre de reactions pour 3a  : ', nr_3a
      print*, ' '
 5    print*, ' '
      print*, '==>> Nom du fichier de sortie : '
      read(5,'(A)',err=5) filename

c Initialisations

      debut_pp = .true.
      debut_cno= .true.
      debut_3a = .true.

c Ouverture du fichier de sortie

      open(unit=1,form='unformatted',status='unknown',file=
     &     filename(:long(filename)))

c Formation du tableau des T/10.**9

      print*, 'Formation du tableau t_t '

      past_t=(t_tmax-t_tmin)/(pnt-1)

      t_t(1)=t_tmin

      do i_t= 2, pnt 
         t_t(i_t)=t_t(i_t-1)+past_t
      end do

c Formation du tableau 'r_t'

        do i_reac= 1, nr_pp+nr_cno+nr_3a

           if(i_reac.le.nr_pp)then
             le_cycle=cycle(1)
             if ( debut_pp ) then
                print*, ' '
                print*, 'Debut cycle PP ---------------------------'
                print*, ' '
                debut_pp=.false.
             end if
           end if
           if(i_reac.gt.nr_pp .and. 
     &        i_reac.le.nr_pp+nr_cno) then
             le_cycle=cycle(2)
             if ( debut_cno ) then
                print*, ' '
                print*, 'Debut cycle CNO --------------------------'
                print*, ' '
                debut_cno=.false.
             end if
           end if
           if(i_reac.gt.nr_pp+nr_cno)then
             le_cycle=cycle(3)             
             if ( debut_3a ) then
                print*, ' '
                print*, 'Debut cycle 3a ---------------------------'
                print*, ' '
                debut_3a=.false.
             end if

           end if

            write(6,*) 'num. reac. :', num_reac(i_reac),' cycle= ',
     &                 le_cycle

            do i_t= 1, pnt

      call taux(le_cycle, num_reac(i_reac),t_t(i_t)*1.d+9,
     &          r_t((i_reac-1)*pnt+i_t),Q,existe)


            if(.not.existe)then
              print*,'Probleme : ',' cycle= ',le_cycle,
     &        ' num_reac = ', num_reac(i_reac)
              stop
            end if

            if(r_t((i_reac-1)*pnt+i_t) .le. 0.d+0) then
              r_t((i_reac-1)*pnt+i_t)=1.d-300
            end if

          end do

       end do


       write(1) pnt,t_t,r_t

      close(1)

      print*, ' '
      print*, 'Test de lecture du fichier binaire ----------------'
      print*, '(A titre de vérification)'
      print*, ' '

      open(unit=20,status='old',form='unformatted',
     &     file=filename(:long(filename)))
       read(20) nt,tt_t,tr_t
c      read(20) nt
c      read(20) (tt_t(i_t), i_t= 1,pnt)
c      read(20) (tr_t(i), i= 1,30*pnt)
      close(20)

      print*, 'nt lu = ', nt
      write(6,*) (tt_t(i), i=1, 10)
      write(6,*) (tt_t(i), i=nt-10,nt)

c      open(unit=30,status='unknown',form='formatted',
c     +     file='test_lecture_nuc.dat')
c      do i_reac= 1, nr_pp+nr_cno+nr_3a
c         do i_t= 1, pnt
c         rt1=tr_t((i_reac-1)*pnt+i_t)
c         rt2= r_t((i_reac-1)*pnt+i_t)
c         if(rt1 .le. 0.d0)then
c      write(30,*) 'i_reac= ', i_reac,' i_t= ',i_t,' itr_t= ', rt1,
c     +            'r_t= ', rt2
c         end if 
c         print*,'tr_t= ', tr_t(i)
c1000  format(2(a,i4),2(a,D20.10))
c         end do
c      end do

      print*, ' '

      close(unit=30)

      print*, '---------------------------------------------------'
      print*, ' '
      print*, ' C''est fini ! On a écrit : ', filename(:long(filename))
      print*, ' (Fichier utilisable par CRES)'
      print*, ' '

      end

c----------------------------------------------------------------------------
  
      subroutine taux( cycle, num_reac, T, NsigmaV, Q, existe)
  
c Donne le ''taux'' NsigmaV pour la reaction numero num_reac du ''cycle''
c (Q: qte d energie degagee, 'existe': variable logique qui porte sur
c l'existence de la reaction)
  
      implicit none
  
      integer num_reac
  
      logical existe
  
      real*8 T, NsigmaV, Q
  
      character*3 cycle
  
  
      if ( cycle .eq. 'pp ') then
	 goto 100
      end if
  
      if ( cycle .eq. 'cno') then
	 goto 200
      end if
  
      if ( cycle .eq. '3a ') then
	 goto 300
      else
	 print*, 'Pour appeler ''taux'' la variable ''cycle'' doit etre'
	 print*, 'egale a ''pp '', ''cno'' ou ''3a '''
	 stop
      end if
  
100   call reac_pp_NACRE( T, num_reac, NsigmaV, Q, existe)
  
      return
  
200   call reac_cno_NACRE( T, num_reac, NsigmaV, Q, existe)
  
      return
  
300   call reac_3a_NACRE( T, num_reac, NsigmaV, Q, existe)
c      if ( num_reac .gt. 1010 ) then
c         print*, 'NsigmaV= ', NsigmaV
c         print*, 'Q      = ', Q
c      end if

      return
  
      end

c***********************************************************************
  
      subroutine reac_pp_NACRE( T, n, NsigmaV, Q, existe)
  
c Taux des reactions nucleaires des chaines PP
c----------------------------------------------
  
c entrees :
c  T   : la temperature en K.
c domaine de validite des expressions analytiques selon Caughlan et
c Fowler :
c                   10**6  <  T  < 10**10  K
  
c  n   : le numero de la reaction.
  
c sorties :
c  NsigmaV : Na*<sigma*v>
c  Q       : quantite de chaleur en MeV
c  existe  :   existe=.true. si la reaction numero n est presente.
c              existe=.false. si la reaction numero n est absente.
  
c Daniel Cordier Decembre 1996
  
      implicit none
  
      integer n, table_n, nbr_reac
  
      logical existe
  
      real*8 T, t9, t9a, NsigmaV, Q, rev_ratio

      real*8 NsigmaV_NACRE, Q_NACRE

      parameter ( nbr_reac= 8 )
  
      dimension table_n(nbr_reac)
  
      data table_n/ 10, 20, 30, 40, 50, 60, 70, 80/
  
      existe=.false.
  
c test  de l'appartenance de T au domaine de validite
  
c      if ( T .gt. 1.d+10 ) then
c         print*, 'Temperature trop elevee!'
c         print*, 'T= ', T
c         stop
c      endif
c      if ( T .lt. 1.d+6 ) then
c         print*, 'Temperature trop basse!'
c         stop
c      endif
  
c On teste si ''n'' appartient a table_n l'ensemble des numeros des
c reactions
  
      if ( n .eq. 10 ) then
         t9=T/10.**9
         existe=.true.
         goto 10
      end if
  
  
      if ( n .eq. 20 ) then
         t9=T/10.**9
         existe=.true.
         goto 20
      end if
  
      if ( n .eq. 30 ) then
         t9=T/10.**9
         existe=.true.
         goto 30
      end if
  
      if ( n .eq. 40 ) then
         t9=T/10.**9
         existe=.true.
         goto 40
      end if
  
      if ( n .eq. 50 ) then
         t9=T/10.**9
         existe=.true.
         goto 50
      end if
  
      if ( n .eq. 60 ) then
         t9=T/10.**9
         existe=.true.
         goto 60
      end if
  
      if ( n .eq. 70 ) then
         t9=T/10.**9
         existe=.true.
         goto 70
      end if
  
      if ( n .eq. 80 ) then
         t9=T/10.**9
         existe=.true.
         goto 80
      end if
  
      if ( .NOT. existe ) then
         return
      end if
  
c------------------------------------------
c reaction 10 : H1(p,e+nu)H2
  
10    NsigmaV=4.01E-15/t9**(2./3.)*exp(-3.380/t9**(1./3.))
     &        *(1.+0.123*t9**(1./3.)+1.09*t9**(2./3.)+0.938*t9)
  
      Q= 1.442d0

      NsigmaV_NACRE=4.08d-15/t9**(2./3.)*exp(-3.381/t9**(1./3.))
     &        *(1.d0+3.82*t9+1.51d0*t9**2+0.144d0*t9**3
     &          -1.14d-2*t9**4)
      Q_NACRE= 1.442d0
 
      call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,10,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return

c------------------------------------------
c reaction 20 : H2(p,g)He3
  
20    NsigmaV=2.24E+03/t9**(2./3.)*exp(-3.720/t9**(1./3.))
     &        *(1.+0.112*t9**(1./3.)+3.38*t9**(2./3.)+2.65*t9)
  
      Q= 5.494d0

      if ( t9 .le. 0.11d0 ) then
         NsigmaV_NACRE=1.81d3/t9**(2./3.)*exp(-3.721d0/t9**(1./3.))
     &        *(1.d0+14.3d0*t9-90.5d0*t9**2+395.d0*t9**3)
      else
         NsigmaV_NACRE=2.58d3/t9**(2./3.)*exp(-3.721d0/t9**(1./3.))
     &        *(1.d0+3.96d0*t9+0.116d0*t9**2)
      end if

      Q_NACRE= 5.493d0
 
      call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,20,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
  
c------------------------------------------
c reaction 30 : He3(He3,2p)He4
  
30    NsigmaV=6.04E+10/t9**(2./3.)*exp(-12.276/t9**(1./3.))
     &        *(1.+0.034*t9**(1./3.)-0.522*t9**(2./3.)-0.124*t9+0.353*
     &        t9**(4./3.)+0.213*t9**(5./3.))
     
      Q= 12.860d0

      NsigmaV_NACRE= 5.59d10/t9**(2./3.)*exp(-12.277/t9**(1./3.))
     &        *(1.d0-0.135d0*t9+2.54d-2*t9**2-1.29d-3*t9**3)

      Q_NACRE= 12.859d0
 
      !call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,30,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
  
c------------------------------------------
c reaction 40 : He3(a,g)Be7 (ecrite He4(He3,g)Be7 dans Caughlan et Fowler)
  
40    t9a=t9/(1.+4.95E-02*t9)
  
      NsigmaV=5.61E+06*t9a**(5./6.)/t9**(3./2.)*
     &        exp(-12.826/t9a**(1./3.))
  
      Q= 1.588d0

      NsigmaV_NACRE=5.46d6/t9**(2./3.)*exp(-12.827/t9**(1./3.))
     &       * (1.d0-0.307d0*t9+8.81d-2*t9**2-1.06d-2*t9**3
     &         +4.46d-4*t9**4)

      Q_NACRE= 1.587d0
 
      call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,40,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
  
c------------------------------------------
c reaction 50 : Be7(e-,nu+g)Li7
  
50    NsigmaV=1.34E-10/t9**(1./2.)*(1.-0.537*t9**(1./3.)+
     &        3.86*t9**(2./3.)
     &        +0.0027/t9*exp(2.515E-03/t9))
  
      Q= 0.862d0

      NsigmaV_NACRE= NsigmaV ! Pas trouvé dans NACRE
      Q_NACRE= Q
c      print*, 'T9 less than or equal to 3'
c      print*, 'Q=0.049 exclusive of nu energy'
c      print*, 'rate must not exceed 1.51E-07/(rho*(1.+X)/2.)'
c      print*, 'for t9 less than 0.001'
   
      call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,50,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
  
c------------------------------------------
c reaction 60 : Li7(p,a)He4
  
60    t9a=t9/(1.+0.759*t9)
  
      NsigmaV=1.096E+09/t9**(2./3.)*exp(-8.472/t9**(1./3.))
     &        -4.830E+08*t9a**(5./6.)/t9**(3./2.)*
     &        exp(-8.472/t9a**(1./3.))
     &        +1.06E+10/t9**(3./2.)*exp(-30.442/t9)
  
      Q= 17.346d0

      NsigmaV_NACRE= 7.20d8/t9**(2./3.)*exp(-8.473d0/t9**(1./3.)
     &      -(t9/6.5)**2)
     &      *(1.d0+1.05d0*t9-0.653d0*t9**2+0.185*t9**3
     &      -2.12d-2*t9**4+9.30d-4*t9**5)
     &      +9.85d6*t9**(0.576)*exp(-10.415d0/t9)

      Q_NACRE= 17.347d0
 
      call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,60,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
  
c------------------------------------------
c reaction 70 : Be7(p,g)B8
  
70    NsigmaV=3.11E+05/t9**(2./3.)*exp(-10.262/t9**(1./3.))
     &        +2.53E+03/t9**(3./2.)*exp(-7.306/t9)

      Q= 0.137

      NsigmaV_NACRE= 2.61d5/t9**(2./3.)*exp(-10.264/t9**(1./3.))
     &        *(1.d0-5.11d-2*t9+4.68d-2*t9**2-6.60d-3*t9**3
     &        +3.12d-4*t9**4)+2.05d3/t9**(3./2.)*exp(-7.345d0/t9)

      Q_NACRE= 0.137d0
 
      call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,70,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
  
c------------------------------------------
c reaction 80 : Be8(,a)He4 (ecrite comme He4(a)Be8 par Caughlan et Fowler)
  
80    print*, '-----------------------------------------------------'
      print*, 'Il y a un probleme avec le REV_RATIO de cette reac'
      print*, 'on a exp(1.0663/t9) ce qui depasse les possibilites'
      print*, 'de la machine pour t9=0.001 par exemple!'
      print*, '-----------------------------------------------------'
      print*, ' '
   
      stop
  
      rev_ratio=(1.40E+10)*t9**(3./2.)*exp( 1.0663/t9)
  
      NsigmaV=  rev_ratio*
     &        7.40E+05/t9**(3./2.)*exp(-1.0663/t9)
     &        +4.164E+09/t9**(2./3.)*exp(-13.490/t9**(1./3.)-(t9/0.098)
     &        **2.)
     &        *(1.+0.031*t9**(1./3.)+8.009*t9**(2./3.)+1.732*t9+49.883
     &        *t9**(4./3.)+27.426*t9**(5./3.))
  
      Q= 0.092
  
c      print*, 'Hypothese : pour une reac inverse on utilise -Q !'

      NsigmaV_NACRE= NsigmaV ! Pas trouvée dans NACRE

      Q_NACRE= Q
 
      call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,80,t9)

      return
  
      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      end
  
c***********************************************************************
  
      subroutine reac_cno_NACRE( T, n, NsigmaV, Q, existe)
  
c Taux des reactions nucleaires du tri-cycle CNO
c------------------------------------------------
  
c entrees :
c  T   : la temperature en K.
c domaine de validite des expressions analytiques selon Caughlan et
c Fowler :
c                   10**6  <  T  < 10**10  K
  
c  n   : le numero de la reaction.
  
c sorties :
c  NsigmaV : Na*<sigma*v>
c  Q       : quantite de chaleur en MeV
c  existe  :   existe=.true. si la reaction numero n est presente.
c              existe=.false. si la reaction numero n est absente.
  
c Daniel Cordier Decembre 1996
  
      implicit none
  
      integer n, table_n, nbr_reac
  
      logical existe
  
      real*8 T, t9, t9a, gt9, NsigmaV, Q, rev_ratio, f1, f2

      real*8 NsigmaV_NACRE, Q_NACRE
  
      parameter ( nbr_reac= 13 )
  
      dimension table_n(nbr_reac)
  
      data table_n/ 110, 120, 130, 140, 150, 160, 170, 180, 190, 200,
     &              210, 220, 230/
  
      existe=.false.
  
c test  de l'appartenance de T au domaine de validite
  
c      if ( T .gt. 1.d+10 ) then
c         print*, 'Temperature trop elevee!'
c         write(6,*) 'T= ', T
c         stop
c      endif
      if ( T .lt. 1.d+6 ) then
         print*, 'Temperature trop basse!'
         write(6,*) 'T= ', T
         stop
      endif
  
c On teste si ''n'' appartient a table_n l'ensemble des numeros des
c reactions
  
      if ( n .eq. 110 ) then
         t9=T/10.**9
         existe=.true.
         goto 110
      end if
  
  
      if ( n .eq. 120 ) then
         t9=T/10.**9
         existe=.true.
         goto 120
      end if
  
      if ( n .eq. 130 ) then
         t9=T/10.**9
         existe=.true.
         goto 130
      end if
  
      if ( n .eq. 140 ) then
         t9=T/10.**9
         existe=.true.
         goto 140
      end if
  
      if ( n .eq. 150 ) then
         t9=T/10.**9
         existe=.true.
         goto 150
      end if
  
      if ( n .eq. 160 ) then
         t9=T/10.**9
         existe=.true.
         goto 160
      end if
  
      if ( n .eq. 170 ) then
         t9=T/10.**9
         existe=.true.
         goto 170
      end if
  
      if ( n .eq. 180 ) then
         t9=T/10.**9
         existe=.true.
         goto 180
      end if
  
      if ( n .eq. 190 ) then
         t9=T/10.**9
         existe=.true.
         goto 190
      end if
  
      if ( n .eq. 200 ) then
         t9=T/10.**9
         existe=.true.
         goto 200
      end if
  
      if ( n .eq. 210 ) then
         t9=T/10.**9
         existe=.true.
         goto 210
      end if
  
      if ( n .eq. 220 ) then
         t9=T/10.**9
         existe=.true.
         goto 220
      end if
  
      if ( n .eq. 230 ) then
         t9=T/10.**9
         existe=.true.
         goto 230
      end if
  
      if ( .NOT. existe ) then
         return
      end if
  
c------------------------------------------
c reaction 110 : N14(p,g)O15
  
110   NsigmaV=4.90E+07/t9**(2./3.)*exp(-15.228/t9**(1./3.)-
     &        (t9/3.294)**2.)
     &        *(1.+0.027*t9**(1./3.)-0.778*t9**(2./3.)-0.149*t9+0.261*
     &        t9**(4./3.)+0.127*t9**(5./3.))
     &        +2.37E+03/t9**(3./2.)*exp(-3.011/t9)+2.19E+04*
     &        exp(-12.530/t9)
  
      Q= 7.297d0

      NsigmaV_NACRE= 4.83d7/t9**(2./3.)*exp(-15.231/t9**(1./3.)
     &        -(t9/0.8d0)**2)
     &        *(1.d0-2.00d0*t9+3.41d0*t9**2-2.43d0*t9**3)
     &        +2.36d3/t9**(3./2.)*exp(-3.010d0/t9)
     &        +6.72d3*t9**(0.380d0)*exp(-9.530d0/t9)

      Q_NACRE= 7.297d0
 
      call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,110,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
  
c------------------------------------------
c reaction 120 : N15(p,a)C12
  
120   NsigmaV=1.08E+12/t9**(2./3.)*exp(-15.251/t9**(1./3.)-
     &        (t9/0.522)**2.)
     &        *(1.+0.027*t9**(1./3.)+2.62*t9**(2./3.)+
     &        0.501*5.36*t9**(4./3.)
     &        +2.60*t9**(5./3.))
     &        +1.19E+08/t9**(3./2.)*exp(-3.676/t9)+5.41E+08/t9**(1./2.)
     &        *exp(-8.926/t9)
     &        +               0.0*
     &        4.72E+08/t9**(3./2.)*exp(-7.721/t9)+2.20E+09/t9**(3./2.)
     &        *exp(-11.418/t9)
  
      Q= 4.966d0

      if ( t9 .le. 2.5d0 ) then
         NsigmaV_NACRE= 1.12d12/t9**(2./3.)*exp(-15.253d0/t9**(1./3.)-
     &        (t9/0.28d0)**2)*(1.d0+4.95d0*t9+143.d0*t9**2)
     &        +1.01d8/t9**(3./2.)*exp(-3.643d0/t9)
     &        +1.19d9/t9**(3./2.)*exp(-7.406d0/t9)
      else
         NsigmaV_NACRE= 4.17d7*t9*(0.917d0)*exp(-3.292d0/t9)
      end if

      Q_NACRE= 4.966d0
 
      call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,120,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
  
c------------------------------------------
c reaction 130 : C12(p,g)N13
  
130   NsigmaV=2.04E+07/t9**(2./3.)*exp(-13.690/t9**(1./3.)-
     &        (t9/1.500)**2.)
     &        *(1.+0.030*t9**(1./3.)+1.19*t9**(2./3.)+0.254*t9+2.06*
     &        t9**(4./3.)+1.12*t9**(5./3.))
     &        +1.08E+05/t9**(3./2.)*exp(-4.925/t9)+2.15E+05/t9**(3./2.)
     &        *exp(-18.179/t9)
  
      Q= 1.944d0

      NsigmaV_NACRE= 2.00d7/t9**(2./3.)*exp(-13.692/t9**(1./3.)
     &               -(t9/0.46d0)**2)
     &             * (1.d0+9.89d0*t9-59.8d0*t9**2+266.d0*t9**3)
     &             + 1.00d5/t9**(3./2.)*exp(-4.913d0/t9)+
     &             4.24d5/t9**(3./2.)*exp(-21.62d0/t9) 

      Q_NACRE= 1.943d0
 
      call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,130,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
  
c------------------------------------------
c reaction 140 : C13(p,g)N14
  
140   NsigmaV=8.01E+07/t9**(2./3.)*exp(-13.717/t9**(1./3.)-
     &        (t9/2.000)**2.)
     &        *(1.+0.030*t9**(1./3.)+0.958*t9**(2./3.)+0.204*t9+1.39*
     &        t9*(4./3.)+0.753*t9**(5./3.))
     &        +1.21E+06/t9**(6./5.)*exp(-5.701d0/t9)
  
      Q= 7.551d0

      NsigmaV_NACRE= 9.57d7/t9**(2./3.)*(1.d0+3.56d0*t9)*
     &        exp(-13.720d0/t9**(1./3.)-t9**2)
     &        +1.50d6/t9**(3./2.)*exp(-5.930d0/t9) + 
     &        6.83d5*t9**(-0.864d0) * exp(-12.057d0/t9)

      Q_NACRE= 7.551d0
 
      call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,140,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
  
c------------------------------------------
c reaction 150 : N13(p,g)O14
  
150   NsigmaV=4.04E+07/t9**(2./3.)*exp(-15.202/t9**(1./3.)-
     &        (t9/1.191)**2.)
     &        *(1.+0.027*t9**(1./3.)-0.803*t9**(2./3.)-0.154*t9+5.00
     &        *t9**(4./3.)+2.44*t9**(5./3.))
     &        +2.43E+05/t9**(3./2.)*exp(-6.348/t9)
  
      Q= 4.628d0

      NsigmaV_NACRE= 4.02d7/t9**(2./3.) * exp(-15.205d0/t9**(1./3.) -
     &        -(t9/0.54d0)**2) 
     &        * (1.d0+3.81d0*t9+18.6d0*t9**2+32.3*t9**3)
     &        +3.25d5 * t9**(-1.35d0) * exp(-5.926d0/t9)

      Q_NACRE= 4.628d0
 
      call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,150,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
  
c------------------------------------------
c reaction 160 : N15(p,g)O16
  
160   NsigmaV=9.78E+08/t9**(2./3.)*exp(-15.251/t9**(1./3.)-
     &        (t9/0.450)**2.)
     &        *(1.+0.027*t9**(1./3.)+0.219*t9**(2./3.)+0.042*t9+6.83
     &        *t9**(4./3.)+3.32*t9**(5./3.))
     &        +1.11E+04/t9**(3./2.)*exp(-3.328/t9)+1.49E+04/t9**(3./2.)
     &        *exp(-4.665/t9)
     &        +3.80E+06/t9**(3./2.)*exp(-11.048/t9)
  
      Q= 12.128d0

      if ( t9 .le. 3.5d0 ) then
         NsigmaV_NACRE= 1.08d9/t9**(2./3.)*exp(-15.254d0/t9**(1./3.)-
     &        (t9/0.34d0)**2) * (1.d0+6.15d0 * t9 + 16.4d0 * t9**2)
     &        +9.23d3/t9**(3./2.)*exp(-3.597/t9) + 3.27d6/t9**(3./2.)
     &        *exp(-11.024d0/t9)
      else
         NsigmaV_NACRE= 3.54d4 * t9**(0.095d0) * exp(-2.306d0/t9)
      end if

      Q_NACRE= 12.127d0
 
      call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,160,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
  
c------------------------------------------
c reaction 170 : O16(p,g)F17
  
170   NsigmaV=1.50E+08/(t9**(2./3.)*(1.+2.13*
     &        (1.-exp(-0.728*t9**(2./3.)) ) ))
     &        *exp(-16.692/t9**(1./3.))
  
      Q= 0.600

      NsigmaV_NACRE= 7.37d7 * exp(-16.696d0/t9**(1./3.)) * t9**(-0.82d0)

      Q_NACRE= 0.600d0
 
      call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,170,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
  
c------------------------------------------
c reaction 180 : O17(p,g)F18
  
180   t9a=t9/(1.+2.69*t9)
  
c f1 et f2 sont des facteurs d'incertitude
      f1= 0.6
      f2= 0.5
  
c      print*, 'Taux de Landre 1990!'
  
      NsigmaV=7.97E+07*t9a**(5./6.)/t9**(3./2.)
     &        *exp(-16.712/t9a**(1./3.))
     &        +1.51E+08/t9**(2./3.)*exp(-16.712/t9**(1./3.))
     &        *(1.+0.025*t9**(1./3.)-0.051*t9**(2./3.)
     &        -8.82E-03*t9)
     &        +1.56E+05/t9*exp(-6.272/t9)
     &        +f1*3.16E-05*t9**(-3./2.)*exp(-0.767/t9)
     &        +f2*98./t9**(-3./2.)*exp(-2.077/t9)
  
c Taux de Caughlan et Fowler 1988
c      NsigmaV=7.97E+07*t9a**(5./6.)/t9**(3./2.)*
c     &        exp(-16.712/t9a**(1./3.))
c     &        +1.51E+08/t9**(2./3.)*exp(-16.712/t9**(1./3.))
c     &        *(1.+0.025*t9**(1./3.)-0.051*t9**(2./3.)-8.82E-03*t9)
c     &        +1.56E+05/t9*exp(-6.272/t9)
c     &        +             0.5*
c     &        1.31E+01/t9**(3./2.)*exp(-1.961/t9)
  
      Q= 5.607d0

      if ( t9 .le. 3.d0 ) then
         NsigmaV_NACRE= 1.50d8/t9**(2./3.) * exp(-16.710d0/t9**(1./3.)
     &        -(t9/0.2d0)**2) + 9.79d-6/t9**(3./2.)
     &        *exp(-0.7659d0/t9)
     &        + 4.15d0/t9**(3./2.) * exp(-2.083d0/t9)
     &        + 7.74d4 * t9**(1.16d0) * exp(-6.342d0/t9)
      else
         NsigmaV_NACRE= 1.74d3 * t9**(0.700d0) * exp(-1.072d0/t9)
      end if

      Q_NACRE= 5.606d0
 
      !call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,180,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
    
c------------------------------------------
c reaction 190 : O17(p,a)N14
  
c Taux donne dans Landre 1990  (f1 et f2 sont deux facteurs d'incertitude)
190   f1= 0.6
      f2= 0.5
  
      NsigmaV=1.53E+07/t9**(2./3.)*exp(-16.712/t9**(1./3.)
     &        -(t9/0.565)**2.)
     &        *(1.+0.025*t9**(1./3.)+5.39*t9**(2./3.)+0.940*t9
     &        +13.5*t9**(4./3.)+5.98*t9**(5./3.))
     &        +2.92E+06*t9*exp(-4.247/t9)
     &        +1.78E+05/t9**(2./3.)*exp(-16.67/t9**(1./3.))/
     &        (0.479*t9**(2./3.)+0.00312)**2.
     &        +f1*2.8E+11*t9*exp(-16.67/t9**(1./3.)-(t9/0.040)**2.)
     &        +f1*2.94E-03/t9**(3./2.)*exp(-0.767/t9)
     &        +f2*98./t9**(3./2.)*exp(-2.077/t9)
  
c       print*, 'Valeur de Landre 1990'
  
c      NsigmaV=1.53E+07/t9**(2./3.)*exp(-16.712/t9**(1./3.)-
c     &        (t9/0.565)**2.)
c     &        *(1.+0.025*t9**(1./3.)+5.39*t9**(2./3.)+
c     &        0.940*t9+13.5*(4./3.)
c     &        +5.98*t9**(5./3.))
c     &        +2.92E+06*t9*exp(-4.247/t9)
c     &        +    0.5*
c     &        (4.81E+10*t9*exp(-16.712/t9**(1./3.)-(t9/0.040)**2.)
c     &        +5.05E-05/t9**(3./2.)*exp(-0.723/t9))
c     &        +       0.5*
c     &        1.31E+01/t9**(3./2.)*exp(-1.961/t9)

      Q= 1.191d0

      if ( t9 .le. 6.d0 ) then
         NsigmaV_NACRE= 9.20d8/t9**(2./3.)*exp(-16.715/t9**(1./3.)
     &        -(t9/0.06d0)**2) 
     &        * (1.d0-80.31d0*t9+2211.*t9**2)
     &        + 9.13d-4/t9**(3./2.) * exp(-0.7667d0/t9)
     &        + 9.68d0/t9**(3./2.) * exp(-2.083d0/t9)
     &        + 8.13d6/t9**(3./2.) * exp(-5.685d0/t9)
     &        + 1.85d6 * t9**(1.591d0) * exp(-4.848d0/t9)
      else
         NsigmaV_NACRE= 8.73d6 * t9**(0.950d0) * exp(-7.508d0/t9)
      end if

      Q_NACRE= 1.192d0
 
      !call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,190,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
  
c------------------------------------------
c reaction 200 : O18(p,g)F19
  
200   NsigmaV=3.45E+08/t9**(2./3.)*exp(-16.729/t9**(1./3.)-
     &        (t9/0.139)**2.)
     &        *(1.+0.025*t9**(1./3.)+2.26*t9**(2./3.)+0.394*t9+30.56*
     &        t9**(4./3.)+13.55*t9**(5./3.))
     &        +1.25E-15/t9**(3./2.)*exp(-0.231/t9)
     &        +1.64E+02/t9**(3./2.)*exp(-1.670/t9)
     &        +1.28E+04*t9**(1./2.)*exp(-5.098/t9)
  
      Q= 7.994d0

      if ( t9 .le. 2.d0 ) then
         NsigmaV_NACRE= 4.59d8/t9**(2./3.) * exp(-16.732/t9**(1./3.)
     &        -(t9/0.15d0)**2) * (1.d0-9.02d0*t9+506.d0*t9**2
     &        2400.d0*t9**3)
     &        +9.91d-17/t9**(3./2.)*exp(-0.232d0/t9) 
     &        + 3.30d-3/t9**(3./2.) * exp(-1.033d0/t9)
     &        + 1.61d2/t9**(3./2.) * exp(-1.665d0/t9) 
     &        + 1.25d4 * t9**(0.458d0) * exp(-5.297d0/t9)
      else
         NsigmaV_NACRE= 1.38d4 * t9**(0.829d0) * exp(-5.919d0/t9)
      end if

      Q_NACRE= 7.994d0

      call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,200,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
  
c------------------------------------------
c reaction 210 : F19(p,a)O16
  
210   NsigmaV=3.55E+11/t9**(2./3.)*exp(-18.113/t9**(1./3.)-
     &        (t9/0.845)**2.)
     &        *(1.+0.023*t9**(1./3.)+1.96*t9**(2./3.)+0.316*t9+2.86*t9
     &        **(4./3.)+1.17*t9**(5./3.))
     &        +3.67E+06/t9**(3./2.)*exp(-3.752/t9)+3.07E+08*exp(-6.019
     &        /t9)
     
      gt9=1.+4.*exp(-2.090/t9)+7.*exp(-16.440/t9)
  
      NsigmaV=NsigmaV/gt9
  
      Q= 8.114d0

      NsigmaV_NACRE= 2.62d11/t9**(2./3.) * exp(-18.116d0*t9**(1./3.)
     &        -(t9/0.185d0)**2)
     &        * (1.d0+6.26d-2*t9+0.285d0*t9**2+4.94d-3*t9**3
     &           + 11.5d0 * t9**4 + 7.40d4 * t9**5)
     &        + 3.80d6/t9**(3./2.) * exp(-3.752d0/t9) 
     &        + 3.27d7 * t9**(-0.193d0) * exp(-6.587d0/t9)
     &        + 7.30d8 * t9**(-0.201d0) * exp(-16.249d0/t9)

      Q_NACRE= 8.114d0

      !call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,210,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
  
c------------------------------------------
c reaction 220 : F19(p,g)Ne20
  
220   NsigmaV=6.04E+07/t9**(2./3.)*exp(-18.113/t9**(1./3.)-
     &        (t9/0.416)**2.)
     &        *(1.+0.023*t9**(1./3.)+2.06*t9**(2./3.)+0.332*t9+3.16
     &        *t9**(4./3.)+1.30*t9**(5./3.))
     &        +6.32E+02/t9**(3./2.)*exp(-3.752/t9)+7.56E+04/t9**(2./7.)*
     &        exp(-5.722/t9)
  
      gt9=1.+4.*exp(-2.090/t9)+7.*exp(-16.440/t9)
  
      NsigmaV=NsigmaV/gt9

      Q= 12.848d0

      if ( t9 .le. 1.5d0 ) then
         NsigmaV_NACRE= 6.37d7/t9**(2./3.) * exp(-18.116/t9**(1./3.))
     &       * (1.d0 + 0.775d0 * t9 + 36.1d0 * t9**2)
     &       +8.27d2/t9**(3./2.) * exp(-3.752d0/t9) 
     &       + 1.28d6 * t9**(-3.667d0) * exp(-9.120d0/t9)
      else
         NsigmaV_NACRE= 3.66d3 * t9**(0.947d0) * exp(-2.245d0/t9)
      end if

      Q_NACRE= 12.843d0

      call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,220,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
  
c------------------------------------------
c reaction 230 : O18(p,a)N15
  
230   NsigmaV=3.63E+11/t9**(2./3.)*exp(-16.729/t9**(1./3.)-
     &        (t9/1.361)**2.)
     &        *(1.+0.025*t9**(1./3.)+1.88*t9**(2./3.)+0.327*t9+4.66*
     &        t9**(4./3.)+2.06*t9**(5./3.))
     &        +9.90E-14/t9**(3./2.)*exp(-0.231/t9)+2.66E+04/t9**(3./2.)
     &        *exp(-1.670/t9)
     &        +2.41E+09/t9**(3./2.)*exp(-7.638/t9)+
     &        1.46E+09/t9*exp(-8.310/t9)
  
      Q= 3.980d0

      NsigmaV_NACRE= 5.58d11/t9**(2./3.) * exp(-16.732d0/t9**(1./3.)
     &        - (t9/0.51d0)**2)
     &        * (1.d0 + 3.2d0 * t9 + 21.8d0 * t9**2)
     &        + 9.91d-14/t9**(3./2.) * exp(-0.232/t9)
     &        + 2.58d4 * exp(-1.665/t9)
     &        + 3.24d8/t9**(0.378d0) * exp(-6.395d0/t9)

      Q_NACRE= 3.981d0

      !call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,230,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
  
      end
  
c***********************************************************************
  
      subroutine reac_3a_NACRE( T, n, NsigmaV, Q, existe)
  
c Taux des reactions nucleaires des chaines 3 alpha
c---------------------------------------------------
  
c entrees :
c  T   : la temperature en K.
c domaine de validite des expressions analytiques selon Caughlan et
c Fowler :
c                   10**6  <  T  < 10**10  K
  
c  n   : le numero de la reaction.
  
c sorties :
c  NsigmaV : Na*<sigma*v>
c  Q       : quantite de chaleur en MeV
c  existe  :   existe=.true. si la reaction numero n est presente.
c              existe=.false. si la reaction numero n est absente.
  
c Daniel Cordier Decembre 1996
  
      implicit none
  
      integer n, table_n, nbr_reac, annee_c12ago16 
  
      logical existe, first
  
      real*8 T, t9, t9a, ft9a, gt9, fpt9a, NsigmaV, Q, 
     &       rev_ratio, He4aBe8, Be8agC12,
     &       NsigmaV_85, NsigmaV_88

      real*8 NsigmaV_NACRE, Q_NACRE
      real*8 aa, aan, a8Be, E1, E2, res

      parameter ( nbr_reac= 10 )
  
      dimension table_n(nbr_reac)

      data table_n/ 1010, 1020, 1030, 1040, 1050, 1060, 1070, 1080,
     &              1090, 1100/

      data first /.true. /

      existe=.false.

c test  de l'appartenance de T au domaine de validite
  
c      if ( T .gt. 1.d+10 ) then
c         print*, 'Temperature trop elevee!'
c         print*, 'T= ', T
c         stop
c      endif
      if ( T .lt. 1.d+6 ) then
         print*, 'Temperature trop basse!'
         print*, 'T= ', T
         stop
      endif
  
c On teste si ''n'' appartient a table_n l'ensemble des numeros des
c reactions
  
      if ( n .eq. 1010 ) then
         t9=T/10.**9
         existe=.true.
         goto 1010
      end if
  
      if ( n .eq. 1020 ) then
         t9=T/10.**9
         existe=.true.
         goto 1020
      end if
  
      if ( n .eq. 1030 ) then
         t9=T/10.**9
         existe=.true.
         goto 1030
      end if
  
      if ( n .eq. 1040 ) then
         t9=T/10.**9
         existe=.true.
         goto 1040
      end if
  
      if ( n .eq. 1050 ) then
         t9=T/10.**9
         existe=.true.
         goto 1050
      end if
  
      if ( n .eq. 1060 ) then
         t9=T/10.**9
         existe=.true.
         goto 1060
      end if
  
      if ( n .eq. 1070 ) then
         t9=T/10.**9
         existe=.true.
         goto 1070
      end if
  
      if ( n .eq. 1080 ) then
         t9=T/10.**9
         existe=.true.
         goto 1080
      end if
  
      if ( n .eq. 1090 ) then
         t9=T/10.**9
         existe=.true.
         goto 1090
      end if
  
      if ( n .eq. 1100 ) then
         t9=T/10.**9
         existe=.true.
         goto 1100
      end if
  
      if ( .NOT. existe ) then
         return
      end if
  
c------------------------------------------
c reaction 1010 : He4(2a,g)C12
  
c He4(a)Be8 
1010  He4aBe8=7.40E+05/t9**(3./2.)*exp(-1.0663/t9)
     &        +4.164E+09/t9**(2./3.)*exp(-13.490/t9**(1./3.)-(t9/0.098)
     &        **2)
     &        *(1.+0.031*t9**(1./3.)+8.009*t9**(2./3.)+1.732*t9+49.883
     &        *t9**(4./3.)+27.426*t9**(5./3.))
  
      Be8agC12=1.30E+02/t9**(3./2.)*exp(-3.3364/t9)
     &         +2.510E+07/t9**(2./3.)*exp(-23.570/t9**(1./3.)-
     &         (t9/0.235)**2.)
     &         *(1.+0.018*t9**(1./3.)+5.249*t9**(2./3.)+0.650*t9+19.176
     &         *t9**(4./3.)+6.034*t9**(5./3.))
     
      NsigmaV=2.90E-16*He4aBe8*Be8agC12
     &        *(0.01+0.2*(1.+4.*exp(-(0.025/t9)**3.263))/
     &        (1.+4.*exp(-(t9/0.025)**9.227)))
     &        +           0.0*
     &        1.35E-07/t9**(3./2.)*exp(-24.811/t9)
  
      if ( t9 .gt. 0.08 ) then
         NsigmaV=2.79E-08/t9**3*exp(-4.4027/t9)
     &           +      0.1*
     &           1.35E-07/t9**(3./2.)*exp(-24.811/t9)
      end if
  
      Q= 14.437d0
c      print*, 'Q=14.437 if C12(a,g)O16 always follows'

c     NACRE :
      aa = 2.43d9/t9**(2./3.) * exp(-13.490/t9**(1./3.)-(t9/0.15d0)**2)
     &     * (1.d0+74.5d0 * t9) + 6.09d5/t9**(3./2.) * exp(-1.054d0/t9)

      if ( t9 .le. 0.03d0) then
         aan= aa * 6.69d-12 * (1.-192.*t9+2.48d4*t9**2-1.50d6*t9**3
     &        + 4.13d7 * t9**4 - 3.90d8 * t9**5)
      else
         aan= aa * 2.42d-12 * (1.d0-1.52*log10(t9) 
     &        + 0.448d0*(log10(t9))**2 + 0.435d0 * (log10(t9))**3)
      end if

      a8Be = 2.76d7/t9**(2./3.) * exp(-23.570d0/t9**(1./3.)
     &       -(t9/0.4d0)**2) * (1.d0+5.47d0 * t9 + 326.d0 * t9**2)
     &       +130.7d0/t9**(3./2.) * exp(-3.338d0/t9)+2.51d4/t9**(3./2.)
     &       *exp(-20.307d0/t9)
      if ( t9 .le. 0.03d0 ) then
         NsigmaV_NACRE=aa * a8Be *3.07d-16*(1.d0-29.1d0*t9+1308.*t9**2)
      else
         NsigmaV_NACRE=aa * a8Be *3.44d-16*(1.+0.0158d0*t9**(-0.65d0))
      end if

      Q_NACRE= 7.274 ! donné par NACRE

      !call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,1010,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
  
c------------------------------------------
c reaction 1020 : C12(a,g)O16
  
c Taux tel qu'il est dans Caughlan et Fowler 1985 :

 1020 NsigmaV_85=2.93E+08/t9**2./(1.+0.0489/t9**(2./3.))**2.
     &        *exp(-32.120/t9**(1./3.)-(t9/3.496)**2.)
     &        +3.14E+08/t9**2./(1.+0.2654/t9**(2./3.))**2.
     &        *exp(-32.120/t9**(1./3.))
     &        +1.25E+03/t9**(3./2.)*exp(-27.499/t9)+1.43E-02
     &        *t9**5.*exp(-15.541/t9)


c Taux tel qu'il est dans Caughlan et Fowler 1988 :
  
      NsigmaV_88=1.04E+08/t9**2./(1.+0.0489/t9**(2./3.))**2.
     &        *exp(-32.120/t9**(1./3.)-(t9/3.496)**2.   )
     &        +1.76E+08/t9**2/(1.+0.2654/t9**(2./3.))**2*exp(-32.120/t9
     &        **(1./3.))
     &        +1.25E+03/t9**(3./2.)*exp(-27.499/t9)+1.43E-02*t9**5.
     &        *exp(-15.541/t9)

      Q= 7.162d0

c NACRE :
      E1 = 6.66d7/t9**2 * exp(-32.123/t9**(1./3.) - (t9/4.6d0)**2)
     &     * (1.d0+2.54d0*t9+1.04d0*t9**2-0.226d0*t9**3)
     &     + 1.39d3/t9**(3./2.)*exp(-28.930/t9)

      E2 = 6.56d7/t9**2 * exp(-32.123/t9**(1./3.) -(t9/1.3)**2)
     &     *(1.d0+9.23d0*t9-13.7d0*t9**2+7.4d0*t9**3)

      res= 19.2d0 * t9**2 * exp(-26.9d0/t9)

      NsigmaV_NACRE= E1 + E2 + res

      Q_NACRE= 7.162d0

      call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,1020,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return

c------------------------------------------
c reaction 1030 : O16(a,g)Ne20
  
1030  NsigmaV=9.37E+09/t9**(2./3.)*exp(-39.757/t9**(1./3.)-
     &        (t9/1.586)**2.)
     &        +6.21E+01/t9**(3./2.)*exp(-10.297/t9)+5.38E+02/t9**(3./2.)
     &        *exp(-12.226/t9)
     &        +1.30E+01*t9**2.*exp(-20.093/t9)
      
      Q= 4.734d0

      NsigmaV_NACRE= 2.68d10/t9**(2./3.) * exp(-39.760d0/t9**(1./3.)
     &        -(t9/1.6d0)**2)
     &        +51.1d0/t9**(3./2.) * exp(-10.32d0/t9)
     &        + 616.1d0/t9**(3./2.) * exp(-12.200d0/t9) + 0.41d0
     &        * t9**(2.966d0) * exp(-11.900/t9)

      Q_NACRE= 4.730d0

      call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,1030,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return

c------------------------------------------
c reaction 1040 : Ne20(a,g)Mg24
  
1040  NsigmaV=4.11E+11/t9**(2./3.)*exp(-46.766/t9**(1./3.)-
     &        (t9/2.219)**2.)
     &        *(1.+0.009*t9**(1./3.)+0.882*t9**(2./3.)+0.055*t9+0.749
     &        *t9**(4./3.)+0.119*t9**(5./3.))
     &        +5.27E+03/t9**(3./2.)*exp(-15.869/t9)+6.51E+03*t9**(1./2.)
     &        *exp(-16.223/t9)
     &        +                       0.5*
     &        4.21E+01/t9**(3./2.)*exp(-9.115/t9)+3.20E+01/t9**(2./3.)
     &        *exp(-9.383/t9)
  
      gt9=1.+5.*exp(-18.960/t9)
  
      NsigmaV=NsigmaV/gt9
  
      Q= 9.312d0

c NACRE :
      if ( t9 .le. 1.d0 ) then
         NsigmaV_NACRE= 8.72/t9**(0.532d0) * exp(-8.995d0/t9)
      else
         NsigmaV_NACRE= 3.74d2 * t9**(2.229d0) * exp(-12.681/t9)
      end if

      Q_NACRE= 9.316d0

      !call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,1040,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
  
c------------------------------------------
c reaction 1050 : C13(a,n)O16
  
1050  NsigmaV=6.77E+15/t9**(2./3.)*exp(-32.329/t9**(1./3.)-
     &        (t9/1.284)**2.)
     &        *(1.+0.013*t9**(1./3.)+2.04*t9**(2./3.)+0.184*t9)
     &        +3.82E+05/t9**(3./2.)*exp(-9.373/t9)+1.41E+06/t9**(3./2.)
     &        *exp(-11.873/t9)
     &        +2.00E+09/t9**(3./2.)*exp(-20.409/t9)+2.92E+09/t9**(3./2.)
     &        *exp(-29.283/t9)
  
      Q= 2.216

c NACRE :
      if ( t9 .le. 4.0d0 ) then
         NsigmaV_NACRE= 3.78d14/t9**2 * exp(-32.33d0/t9**(1./3.) 
     &        -(t9/0.71d0)**2)
     &        *(1.+46.8d0*t9-292.d0*t9**2+738.d0*t9**3)
     &        +2.30d7 * t9**(0.45d0) * exp(-13.03d0/t9)
      else
         NsigmaV_NACRE= 7.59d6 * t9**(1.078) * exp(-12.056/t9)
      end if

      Q_NACRE= 2.216d0
 
      call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,1050,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
  
c------------------------------------------
c reaction 1060 : O17(a,n)Ne20
  
1060  t9a=t9/(1.+0.0268*t9+0.0232*t9**(5./3.)/(1.+0.0268*t9)**(2./3.))
  
      gt9=1.+exp(-10.106/t9)/3.
  
      NsigmaV=1.03E+18/gt9*t9a**(5./6.)/t9**(3./2.)*
     &        exp(-39.914/t9a**(1./3.))
  
      Q= 0.590d0

c NACRE :
      NsigmaV_NACRE= 4.38d17/t9**(2./3.) * exp(-39.918d0/t9**(1./3.)
     &  -(t9/1.1d0)**2)
     &  +1.73d3/t9**(3./2.) * exp(-8.55d0/t9) 
     &  +7.50d5 * t9**(1.83d0) * exp(-13.8d0/t9)

      Q_NACRE= 0.586d0

      call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,1060,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
  
c------------------------------------------
c reaction 1070 : N14(a,g)F18
  
1070  NsigmaV=7.78E+09/t9**(2./3.)*exp(-36.031/t9**(1./3.)-
     &        (t9/0.881)**2.)
     &        *(1.+0.012*t9**(1./3.)+1.45*t9**(2./3.)+0.117*t9+1.97*
     &        t9**(4./3.)+0.406*t9**(5./3.))
     &        +2.36E-10/t9**(3./2.)*exp(-2.798/t9)+2.03E+00/t9**(3./2.)
     &        *exp(-5.054/t9)
     &        +1.15E+04/t9**(2./3.)*exp(-12.310/t9)
  
      Q= 4.415d0

c NACRE :
      if ( t9 .le. 2.d0 ) then
         NsigmaV_NACRE= 7.93d11/t9**(2./3.) * exp(-36.035d0/t9**(1./3.)
     &        -(t9/0.07d0)**2)
     &        +1.85d-10/t9**(3./2.) * exp(-2.750d0/t9)
     &        +2.62d0/t9**(3./2.) * exp(-5.045d0/t9)
     &        +2.93d3*t9**(0.344d0) * exp(-10.561d0/t9)
      else
         NsigmaV_NACRE= 1.52d2 * t9**(1.567d0) * exp(-6.315/t9)
      end if

      Q_NACRE= 4.415d0

      call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,1070,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
  
c------------------------------------------
c reaction 1080 : O18(a,g)Ne22
  
1080  NsigmaV=1.82E+12/t9**(2./3.)*exp(-40.057/t9**(1./3.)-
     &        (t9/0.343)**2.)
     &        *(1.+0.010*t9**(1./3.)+0.988*t9**(2./3.)+0.072*t9
     &        +3.17*t9**(4./3.)+0.586*t9**(5./3.))
     &        +7.54E+00/t9**(3./2.)*exp(-6.228/t9)+3.48E+01/t9**(3./2.)
     &        *exp(-7.301/t9)
     &        +6.23E+03*t9*exp(-16.987/t9)
     &        +                0.5*
     &        1.00E-11/t9**(3./2.)*exp(-1.994/t9)
  
      Q= 9.669d0

c NACRE :
      if ( t9 .le. 6.d0 ) then
         NsigmaV_NACRE= 1.95d-13/t9**(3./2.) * exp(-2.069d0/t9)
     &      + 1.56d-2/t9**(3./2.) * exp(-4.462/t9)
     &      + 10.1d0/t9**(3./2.) * exp(-6.391d0/t9) + 44.1d0/t9**(3./2.)
     &      * exp(-7.389d0/t9)
     &      +3.44d5 * t9**(-0.5d0) * exp(-22.103d0/t9)
      else
         NsigmaV_NACRE= 3.31d5 * t9**(-0.221d0) * exp(-24.990d0/t9)
      end if

      Q_NACRE= 9.667d0

      !call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,1080,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
  
c------------------------------------------
c reaction 1090 : Ne22(a,n)Mg25
  
1090  t9a=t9/(1.+0.0548*t9)
      gt9=1.+5.0*exp(-14.781/t9)
      ft9a=exp(-(0.197/t9a)**4.82)
  
      NsigmaV=4.16E+19*ft9a/gt9*t9a**(5./6.)/t9**(3./2.)*exp(-47.004/
     &        t9a**(1./3.))
     &        +1.44E-04/gt9*exp(-5.577/t9)
  
      Q= -0.481d0

c NACRE :
      if ( t9 .le. 2.d0 ) then
         NsigmaV_NACRE= 7.40d0 * exp(-7.79d0/t9) + 1.30d-4*t9**(0.83d0)
     &      * exp(-5.52/t9)
     &      +9.41d3 *t9**(2.78d0) * exp(-11.7d0/t9)
     &      +8.59d6 *t9**(0.892d0) * exp(-24.4d0/t9)
      else
         NsigmaV_NACRE= 1.51d5 *t9**(2.879d0) * exp(-16.717d0/t9)
      end if

      Q_NACRE= -0.478d0

      !call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,1090,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
     
c------------------------------------------
c reaction 1100 : Ne22(a,g)Mg26
  
1100  t9a=t9/(1.+0.0548*t9)
      gt9=1.+5.0*exp(-14791/t9)
      ft9a=exp(-(0.197/t9a)**4.82)
      fpt9a=exp(-(t9a/0.249)**2.31)
  
      NsigmaV=4.16E+19*fpt9a/gt9*t9a**(5./6.)/t9**(3./2.)*exp(-47.004/
     &        t9a**(1./3.))
  
      Q= 10.612

c NACRE :
      if ( t9 .le. 1.25d0 ) then
         NsigmaV_NACRE= 3.55d-9/t9**(3./2.) * exp(-3.927d0/t9)
     &     +7.07d-1 *t9**(-1.064d0)* exp(-7.759d0/t9)
     &     +1.27d-3 *t9**(-2.556d0)* exp(-6.555d0/t9)
      else
         NsigmaV_NACRE= 1.76d0 *t9**(3.322d0) * exp(-12.412d0/t9)
      end if

      Q_NACRE= 10.615d0

      !call check(NsigmaV,NsigmaV_NACRE,Q,Q_NACRE,1100,t9)

      NsigmaV = NsigmaV_NACRE
      Q       = Q_NACRE

      return
  
      end

c***********************************************************************

      subroutine check(Nsig1,Nsig2,Q1,Q2,nreac,t9)

c Vérification : on cherche d'éventuelles différences supérieures
c                à 5% entre Caughlan & Fowler et NACRE.
c D. Cordier, 2 avril 2003.

      implicit none
      integer nreac
      real*8 Nsig1,Nsig2,Q1,Q2,t9

      if ( (Nsig2 .ne. 0.d0) .AND. 
     &     (abs(Nsig1-Nsig2)/Nsig2 .gt. 1.50d0) ) then
         print*, 'Pour la réaction num.  = ', nreac
         print*, 'abs(Nsig1-Nsig2)/Nsig2 = ', abs(Nsig1-Nsig2)/Nsig2
         print*, 't9                     = ', t9
         print*, ' '
         print*, ' Nsig1   = ', Nsig1
         print*, ' Nsig2   = ', Nsig2
         print*, ' '
         print*, ' Q1      = ', Q1
         print*, ' Q2      = ', Q2
         print*, ' '
         stop
      end if

      if ( (Q2 .ne. 0.d0) .AND. 
     &     (abs(Q1-Q2)/Q2 .gt. 0.05d0) ) then
         print*, 'Pour la réaction num.  = ', nreac
         print*, 'abs(Q1-Q2)/Q2          = ', abs(Q1-Q2)/Q2
         print*, 't9                     = ', t9
         print*, ' '
         print*, ' Nsig1   = ', Nsig1
         print*, ' Nsig2   = ', Nsig2
         print*, ' '
         print*, ' Q1      = ', Q1
         print*, ' Q2      = ', Q2
         print*, ' '
         stop
      end if

      return

      end

c***********************************************************************
 
	function long(a)
 
c	Determination de la longueur effective d'une chaine
 
c	Auteur: G.Gonczi, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91
 
	implicit none
 
	integer*4 long
 
	character*(*)a
 
	long=len(a)
 
	do while(a(long:long) .eq. ' ' .and. long .gt. 1)
	 long=long-1
	enddo
 
	return
 
	end

  
  
  
