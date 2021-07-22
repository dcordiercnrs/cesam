c****************************************************************************

  
c           Programme de tabulation des taux de reaction 
c                   dans un fichier binaire


c Possibilites de faire des tests :
c      (1) choix entre des taux pour C12(a,g)O16
c      (2) multiplication de C12(a,g)O16 par un coefficient arbitraire.
c      (3) suppression des reactions autre que 3 alpha dans la combustion
c          de He

c****************************************************************************

c Daniel Cordier Septembre 1997, Aout 2000.

      implicit none

      integer i, i_t, i_reac, i_cycle, pnt, nr_pp, nr_cno, nr_3a,
     +        num_reac, i_reac_max, nt, long, annee_c12ago16 

      real*8 t_tmin, t_tmax, past_t, Q, rt1, rt2, coefftest,
     +       coefftest3a

      logical existe, debut_pp, debut_cno, debut_3a, test, only_3a,
     +        test2

      character filename*40, cycle*3, le_cycle*3, rep*1

      parameter(pnt=300, t_tmin=1.d-3, t_tmax=1.d+1 )

c Rq. : t_t = T / 1.d+9
c       pnt= nombre de valeurs de la temperature

      parameter (nr_pp=7, nr_cno=13, nr_3a=10)

c nombre de reactions dans PP, CNO et 3 alpha

      real*8 t_t(pnt), r_t(pnt*(nr_pp+nr_cno+nr_3a)), tt_t(pnt),
     +                tr_t(pnt*(nr_pp+nr_cno+nr_3a)),
     +     coeff_c12ago16

      dimension cycle(3), num_reac(nr_pp+nr_cno+nr_3a)

      common /trois_alpha/ only_3a

      common /c12ago16/coeff_c12ago16,annee_c12ago16

      data cycle/'pp ', 'cno', '3a '/

      data num_reac/10,20,30,40,60,50,70,   ! reactions chaines PP
     +              130,140,110,160,120,170,190,150,180,200, ! cycle CNO
     +              210,220,230,
     +              1010,1020,1030,1040,1050,1060,1070,1080,1090,
     +              1100/

      external long

      open(unit=88,status='unknown',form='formatted',access='append',
     +     file='rapp_8588.dat')
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

 10   print*, ' '
      print*, '==>> Annee pour C12(a,g)O16 : (1985 ou 1988)'
      print*, '     (la valeur "standard" est 1985)'
      read(5,*,err=10) annee_c12ago16
c     -----------------------------------------------------------------
      print*, ' '
      print*, '==>> Fait-on un test ? (o/n)'
      print*, '     (multiplication de C12(a,g)O16 par un coeff arbitrai
     +re, défaut : n)'
      rep='n'
      read(5,'(A)') rep
      if ( rep .eq. 'o' ) then
         test=.true.
 20      print*, 'Coeff. multiplicateur pour C12(a,g)O16 : '
         read(5,*) coeff_c12ago16
         print*, 'Coeff. multiplicateur pour les reactions du cycle CNO 
     +:'
         read(5,*,err=20) coefftest
         print*, ' '
         print*, 'On va faire un test avec :'
         print*, '--------------------------'
         print*, 'coeff_c12ago16= ', coeff_c12ago16
         print*, 'coefftest (pour CNO)= ',  coefftest
         print*, ' '
         print*, 'Ok ?'
         pause
      else
         test=.false.
         coeff_c12ago16=1.d0
         coefftest=1.d0
      end if
c     -----------------------------------------------------------------
      print*, ' '
      print*, '==>> Multiple-t-on  3alpha par un coefficient arbitraire 
     +? (Défaut : n)'
      rep='n'
      test2=.false.
      read(5,'(A)') rep
      if ( rep .eq. 'o' ) then
         print*, ' '
         print*, '     ******************************************'
         print*, '     ******************************************'
         print*, '==>> ATTENTION on multiple 3 alpha par un coeff'
         print*, '     ******************************************'
         print*, '     ******************************************'
         print*, ' '
         print*, '==>> OK ?'
         pause
         test2=.true.
 30      print*, 'Valeur du coefficient : ?'
         read(5,*,err=30) coefftest3a
      else
         test2=.false.
         coefftest3a=1.d0
      end if
c     -----------------------------------------------------------------
      print*, ' '
      print*, '==>> Supprime-t-on les reactions autre que 3alpha dans'
      print*, '     la combustion de He ? (o/n) défaut : n'
      rep='n'
      only_3a=.false.
      read(5,'(A)') rep
      if ( rep .eq. 'o' ) then
         print*, ' '
         print*, '     ******************************************'
         print*, '     ******************************************'
         print*, '==>> ATTENTION on met les taux de reactions des'
         print*, '      reactions autre que 3 alpha a ZERO !!!'
         print*, '     ******************************************'
         print*, '     ******************************************'
         print*, ' '
         print*, '==>> OK ?'
         pause
         only_3a=.true.
      else
         only_3a=.false.
      end if

c Initialisations

      debut_pp = .true.
      debut_cno= .true.
      debut_3a = .true.

c Ouverture du fichier de sortie

      open(unit=1,form='unformatted',status='unknown',file=
     +     filename(:long(filename)))

c Formation du tableau des T/10.**9

      print*, 'Formation du tableau t_t '

      past_t=(t_tmax-t_tmin)/(pnt-1)

      t_t(1)=t_tmin

      do i_t= 2, pnt 
         t_t(i_t)=t_t(i_t-1)+past_t
      end do

c      print*, 't_t= ', t_t
c      stop

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
     +        i_reac.le.nr_pp+nr_cno) then
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
     +                 le_cycle

            do i_t= 1, pnt

      call taux(le_cycle, num_reac(i_reac),t_t(i_t)*1.d+9,
     +          r_t((i_reac-1)*pnt+i_t),Q,existe)

            if ( test .AND. le_cycle .eq. 'cno') then
c               print*, 'TEST1'
               r_t((i_reac-1)*pnt+i_t)=r_t((i_reac-1)*pnt+i_t)*coefftest
            end if

            if ( test2 .AND. le_cycle .eq. '3a') then
c               print*, 'TEST2       :'
c               print*, 'coefftest3a = ', coefftest3a
               r_t((i_reac-1)*pnt+i_t)=r_t((i_reac-1)*pnt+i_t)*
     +coefftest3a
            end if

            if(.not.existe)then
              print*,'Probleme : ',' cycle= ',le_cycle,
     +        ' num_reac = ', num_reac(i_reac)
              stop
            end if
            
c            if ( num_reac(decal+i_reac) .eq. 140 .and.
c     +           i_t .eq. 1)then
c               print*, 'reaction 140, i_t=1, r_t=', 
c     +                 r_t((i_reac-1)*pnt+i_t)
c            end if

            if(r_t((i_reac-1)*pnt+i_t) .le. 0.d+0) then
c              print*, 'Pb. avec le taux suivant :'
c              print*, 'Cycle= ',le_cycle
c              print*, 'Reaction num.: ', num_reac(decal+i_reac)
c              print*, 'r_t= ',r_t((i_reac-1)*pnt+i_t)
              r_t((i_reac-1)*pnt+i_t)=1.d-300
c              stop
            end if

          end do

       end do


       write(1) pnt,t_t,r_t
c      write(1) pnt
c      write(1) (t_t(i_t), i_t= 1, pnt)
c      write(1) (r_t(i), i= 1, 30*pnt)

      close(1)

      print*, ' '
      print*, 'Test de lecture du fichier binaire ----------------'
      print*, '(A titre de vérification)'
      print*, ' '

      open(unit=20,status='old',form='unformatted',
     +     file=filename(:long(filename)))
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
      print*, ' (Fichier utilisable par CESAM 3.2, routine D.C.)'
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
  
100   call reac_pp( T, num_reac, NsigmaV, Q, existe)
  
      return
  
200   call reac_cno( T, num_reac, NsigmaV, Q, existe)
  
      return
  
300   call reac_3a( T, num_reac, NsigmaV, Q, existe)
c      if ( num_reac .gt. 1010 ) then
c         print*, 'NsigmaV= ', NsigmaV
c         print*, 'Q      = ', Q
c      end if

      return
  
      end


c***********************************************************************
  
      subroutine reac_pp( T, n, NsigmaV, Q, existe)
  
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
     +        *(1.+0.123*t9**(1./3.)+1.09*t9**(2./3.)+0.938*t9)
  
      Q= 1.442
  
c      print*, 'T9 less than or equal to 3'
  
      return
  
c------------------------------------------
c reaction 20 : H2(p,g)He3
  
20    NsigmaV=2.24E+03/t9**(2./3.)*exp(-3.720/t9**(1./3.))
     +        *(1.+0.112*t9**(1./3.)+3.38*t9**(2./3.)+2.65*t9)
  
      Q= 5.494
  
      return
  
c------------------------------------------
c reaction 30 : He3(He3,2p)He4
  
30    NsigmaV=6.04E+10/t9**(2./3.)*exp(-12.276/t9**(1./3.))
     +        *(1.+0.034*t9**(1./3.)-0.522*t9**(2./3.)-0.124*t9+0.353*
     +        t9**(4./3.)+0.213*t9**(5./3.))
     
      Q= 12.860
  
      return
  
c------------------------------------------
c reaction 40 : He3(a,g)Be7 (ecrite He4(He3,g)Be7 dans Caughlan et Fowler)
  
40    t9a=t9/(1.+4.95E-02*t9)
  
      NsigmaV=5.61E+06*t9a**(5./6.)/t9**(3./2.)*
     +        exp(-12.826/t9a**(1./3.))
  
      Q= 1.588
  
      return
  
c------------------------------------------
c reaction 50 : Be7(e-,nu+g)Li7
  
50    NsigmaV=1.34E-10/t9**(1./2.)*(1.-0.537*t9**(1./3.)+
     +        3.86*t9**(2./3.)
     +        +0.0027/t9*exp(2.515E-03/t9))
  
      Q= 0.862
  
c      print*, 'T9 less than or equal to 3'
c      print*, 'Q=0.049 exclusive of nu energy'
c      print*, 'rate must not exceed 1.51E-07/(rho*(1.+X)/2.)'
c      print*, 'for t9 less than 0.001'
  
      return
  
c------------------------------------------
c reaction 60 : Li7(p,a)He4
  
60    t9a=t9/(1.+0.759*t9)
  
      NsigmaV=1.096E+09/t9**(2./3.)*exp(-8.472/t9**(1./3.))
     +        -4.830E+08*t9a**(5./6.)/t9**(3./2.)*
     +        exp(-8.472/t9a**(1./3.))
     +        +1.06E+10/t9**(3./2.)*exp(-30.442/t9)
  
      Q= 17.346
  
      return
  
c------------------------------------------
c reaction 70 : Be7(p,g)B8
  
70    NsigmaV=3.11E+05/t9**(2./3.)*exp(-10.262/t9**(1./3.))
     +        +2.53E+03/t9**(3./2.)*exp(-7.306/t9)
  
      Q= 0.137
  
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
     +        7.40E+05/t9**(3./2.)*exp(-1.0663/t9)
     +        +4.164E+09/t9**(2./3.)*exp(-13.490/t9**(1./3.)-(t9/0.098)
     +        **2.)
     +        *(1.+0.031*t9**(1./3.)+8.009*t9**(2./3.)+1.732*t9+49.883
     +        *t9**(4./3.)+27.426*t9**(5./3.))
  
      Q= 0.092
  
c      print*, 'Hypothese : pour une reac inverse on utilise -Q !'
  
      return
  
      end
  
c***********************************************************************
  
      subroutine reac_cno( T, n, NsigmaV, Q, existe)
  
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
  
      parameter ( nbr_reac= 13 )
  
      dimension table_n(nbr_reac)
  
      data table_n/ 110, 120, 130, 140, 150, 160, 170, 180, 190, 200,
     +              210, 220, 230/
  
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
     +        (t9/3.294)**2.)
     +        *(1.+0.027*t9**(1./3.)-0.778*t9**(2./3.)-0.149*t9+0.261*
     +        t9**(4./3.)+0.127*t9**(5./3.))
     +        +2.37E+03/t9**(3./2.)*exp(-3.011/t9)+2.19E+04*
     +        exp(-12.530/t9)
  
      Q= 7.297
  
      return
  
c------------------------------------------
c reaction 120 : N15(p,a)C12
  
120   NsigmaV=1.08E+12/t9**(2./3.)*exp(-15.251/t9**(1./3.)-
     +        (t9/0.522)**2.)
     +        *(1.+0.027*t9**(1./3.)+2.62*t9**(2./3.)+
     +        0.501*5.36*t9**(4./3.)
     +        +2.60*t9**(5./3.))
     +        +1.19E+08/t9**(3./2.)*exp(-3.676/t9)+5.41E+08/t9**(1./2.)
     +        *exp(-8.926/t9)
     +        +               0.0*
     +        4.72E+08/t9**(3./2.)*exp(-7.721/t9)+2.20E+09/t9**(3./2.)
     +        *exp(-11.418/t9)
  
      Q= 4.966
  
      return
  
c------------------------------------------
c reaction 130 : C12(p,g)N13
  
130   NsigmaV=2.04E+07/t9**(2./3.)*exp(-13.690/t9**(1./3.)-
     +        (t9/1.500)**2.)
     +        *(1.+0.030*t9**(1./3.)+1.19*t9**(2./3.)+0.254*t9+2.06*
     +        t9**(4./3.)+1.12*t9**(5./3.))
     +        +1.08E+05/t9**(3./2.)*exp(-4.925/t9)+2.15E+05/t9**(3./2.)
     +        *exp(-18.179/t9)
  
      Q= 1.944
  
      return
  
c------------------------------------------
c reaction 140 : C13(p,g)N14
  
140   NsigmaV=8.01E+07/t9**(2./3.)*exp(-13.717/t9**(1./3.)-
     +        (t9/2.000)**2.)
     +        *(1.+0.030*t9**(1./3.)+0.958*t9**(2./3.)+0.204*t9+1.39*
     +        t9*(4./3.)+0.753*t9**(5./3.))
     +        +1.21E+06/t9**(6./5.)*exp(-5.701/t9)
  
      Q= 7.551
  
      return
  
c------------------------------------------
c reaction 150 : N13(p,g)O14
  
150   NsigmaV=4.04E+07/t9**(2./3.)*exp(-15.202/t9**(1./3.)-
     +        (t9/1.191)**2.)
     +        *(1.+0.027*t9**(1./3.)-0.803*t9**(2./3.)-0.154*t9+5.00
     +        *t9**(4./3.)+2.44*t9**(5./3.))
     +        +2.43E+05/t9**(3./2.)*exp(-6.348/t9)
  
      Q= 4.628
  
      return
  
c------------------------------------------
c reaction 160 : N15(p,g)O16
  
160   NsigmaV=9.78E+08/t9**(2./3.)*exp(-15.251/t9**(1./3.)-
     +        (t9/0.450)**2.)
     +        *(1.+0.027*t9**(1./3.)+0.219*t9**(2./3.)+0.042*t9+6.83
     +        *t9**(4./3.)+3.32*t9**(5./3.))
     +        +1.11E+04/t9**(3./2.)*exp(-3.328/t9)+1.49E+04/t9**(3./2.)
     +        *exp(-4.665/t9)
     +        +3.80E+06/t9**(3./2.)*exp(-11.048/t9)
  
      Q= 12.128
  
      return
  
c------------------------------------------
c reaction 170 : O16(p,g)F17
  
170   NsigmaV=1.50E+08/(t9**(2./3.)*(1.+2.13*
     +        (1.-exp(-0.728*t9**(2./3.)) ) ))
     +        *exp(-16.692/t9**(1./3.))
  
      Q= 0.600
  
      return
  
c------------------------------------------
c reaction 180 : O17(p,g)F18
  
180   t9a=t9/(1.+2.69*t9)
  
c f1 et f2 sont des facteurs d'incertitude
      f1= 0.6
      f2= 0.5
  
c      print*, 'Taux de Landre 1990!'
  
      NsigmaV=7.97E+07*t9a**(5./6.)/t9**(3./2.)
     +        *exp(-16.712/t9a**(1./3.))
     +        +1.51E+08/t9**(2./3.)*exp(-16.712/t9**(1./3.))
     +        *(1.+0.025*t9**(1./3.)-0.051*t9**(2./3.)
     +        -8.82E-03*t9)
     +        +1.56E+05/t9*exp(-6.272/t9)
     +        +f1*3.16E-05*t9**(-3./2.)*exp(-0.767/t9)
     +        +f2*98./t9**(-3./2.)*exp(-2.077/t9)
  
c Taux de Caughlan et Fowler 1988
c      NsigmaV=7.97E+07*t9a**(5./6.)/t9**(3./2.)*
c     +        exp(-16.712/t9a**(1./3.))
c     +        +1.51E+08/t9**(2./3.)*exp(-16.712/t9**(1./3.))
c     +        *(1.+0.025*t9**(1./3.)-0.051*t9**(2./3.)-8.82E-03*t9)
c     +        +1.56E+05/t9*exp(-6.272/t9)
c     +        +             0.5*
c     +        1.31E+01/t9**(3./2.)*exp(-1.961/t9)
  
      Q= 5.607
  
      return
    
c------------------------------------------
c reaction 190 : O17(p,a)N14
  
c Taux donne dans Landre 1990  (f1 et f2 sont deux facteurs d'incertitude)
190   f1= 0.6
      f2= 0.5
  
      NsigmaV=1.53E+07/t9**(2./3.)*exp(-16.712/t9**(1./3.)
     +        -(t9/0.565)**2.)
     +        *(1.+0.025*t9**(1./3.)+5.39*t9**(2./3.)+0.940*t9
     +        +13.5*t9**(4./3.)+5.98*t9**(5./3.))
     +        +2.92E+06*t9*exp(-4.247/t9)
     +        +1.78E+05/t9**(2./3.)*exp(-16.67/t9**(1./3.))/
     +        (0.479*t9**(2./3.)+0.00312)**2.
     +        +f1*2.8E+11*t9*exp(-16.67/t9**(1./3.)-(t9/0.040)**2.)
     +        +f1*2.94E-03/t9**(3./2.)*exp(-0.767/t9)
     +        +f2*98./t9**(3./2.)*exp(-2.077/t9)
  
c       print*, 'Valeur de Landre 1990'
  
c      NsigmaV=1.53E+07/t9**(2./3.)*exp(-16.712/t9**(1./3.)-
c     +        (t9/0.565)**2.)
c     +        *(1.+0.025*t9**(1./3.)+5.39*t9**(2./3.)+
c     +        0.940*t9+13.5*(4./3.)
c     +        +5.98*t9**(5./3.))
c     +        +2.92E+06*t9*exp(-4.247/t9)
c     +        +    0.5*
c     +        (4.81E+10*t9*exp(-16.712/t9**(1./3.)-(t9/0.040)**2.)
c     +        +5.05E-05/t9**(3./2.)*exp(-0.723/t9))
c     +        +       0.5*
c     +        1.31E+01/t9**(3./2.)*exp(-1.961/t9)
  
      Q= 1.191
  
      return
  
c------------------------------------------
c reaction 200 : O18(p,g)F19
  
200   NsigmaV=3.45E+08/t9**(2./3.)*exp(-16.729/t9**(1./3.)-
     +        (t9/0.139)**2.)
     +        *(1.+0.025*t9**(1./3.)+2.26*t9**(2./3.)+0.394*t9+30.56*
     +        t9**(4./3.)+13.55*t9**(5./3.))
     +        +1.25E-15/t9**(3./2.)*exp(-0.231/t9)
     +        +1.64E+02/t9**(3./2.)*exp(-1.670/t9)
     +        +1.28E+04*t9**(1./2.)*exp(-5.098/t9)
  
      Q= 7.994
  
      return
  
c------------------------------------------
c reaction 210 : F19(p,a)O16
  
210   NsigmaV=3.55E+11/t9**(2./3.)*exp(-18.113/t9**(1./3.)-
     +        (t9/0.845)**2.)
     +        *(1.+0.023*t9**(1./3.)+1.96*t9**(2./3.)+0.316*t9+2.86*t9
     +        **(4./3.)+1.17*t9**(5./3.))
     +        +3.67E+06/t9**(3./2.)*exp(-3.752/t9)+3.07E+08*exp(-6.019
     +        /t9)
     
       gt9=1.+4.*exp(-2.090/t9)+7.*exp(-16.440/t9)
  
       NsigmaV=NsigmaV/gt9
  
       Q= 8.114
  
       return
  
c------------------------------------------
c reaction 220 : F19(p,g)Ne20
  
220   NsigmaV=6.04E+07/t9**(2./3.)*exp(-18.113/t9**(1./3.)-
     +        (t9/0.416)**2.)
     +        *(1.+0.023*t9**(1./3.)+2.06*t9**(2./3.)+0.332*t9+3.16
     +        *t9**(4./3.)+1.30*t9**(5./3.))
     +        +6.32E+02/t9**(3./2.)*exp(-3.752/t9)+7.56E+04/t9**(2./7.)*
     +        exp(-5.722/t9)
  
      gt9=1.+4.*exp(-2.090/t9)+7.*exp(-16.440/t9)
  
      NsigmaV=NsigmaV/gt9
  
      Q= 12.848
  
      return
  
c------------------------------------------
c reaction 230 : O18(p,a)N15
  
230   NsigmaV=3.63E+11/t9**(2./3.)*exp(-16.729/t9**(1./3.)-
     +        (t9/1.361)**2.)
     +        *(1.+0.025*t9**(1./3.)+1.88*t9**(2./3.)+0.327*t9+4.66*
     +        t9**(4./3.)+2.06*t9**(5./3.))
     +        +9.90E-14/t9**(3./2.)*exp(-0.231/t9)+2.66E+04/t9**(3./2.)
     +        *exp(-1.670/t9)
     +        +2.41E+09/t9**(3./2.)*exp(-7.638/t9)+
     +        1.46E+09/t9*exp(-8.310/t9)
  
      Q= 3.980
  
      return
  
      end
  
c***********************************************************************
  
      subroutine reac_3a( T, n, NsigmaV, Q, existe)
  
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
  
      logical existe, first, only_3a
  
      real*8 T, t9, t9a, ft9a, gt9, fpt9a, NsigmaV, Q, 
     +       rev_ratio, He4aBe8, Be8agC12, coeff_c12ago16,
     +       NsigmaV_85, NsigmaV_88, rapp_8588
  
      parameter ( nbr_reac= 10 )
  
      dimension table_n(nbr_reac)

      common /trois_alpha/ only_3a

      common /c12ago16/coeff_c12ago16, annee_c12ago16

      data table_n/ 1010, 1020, 1030, 1040, 1050, 1060, 1070, 1080,
     +              1090, 1100/

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
         if ( .NOT. only_3a ) then
            t9=T/10.**9
            existe=.true.
            goto 1020
         else
            NsigmaV=1.d-33
            Q=0.d0
            existe=.true.
            return
         end if
      end if
  
      if ( n .eq. 1030 ) then
         if ( .NOT. only_3a ) then
            t9=T/10.**9
            existe=.true.
            goto 1030
         else
            NsigmaV=1.d-33
            Q=0.d0
            existe=.true.
            return
         end if
      end if
  
      if ( n .eq. 1040 ) then
         if ( .NOT. only_3a ) then
            t9=T/10.**9
            existe=.true.
            goto 1040
         else
            NsigmaV=1.d-33
            Q=0.d0
            existe=.true.
            return
         end if
      end if
  
      if ( n .eq. 1050 ) then
         if ( .NOT. only_3a ) then
            t9=T/10.**9
            existe=.true.
            goto 1050
         else
            NsigmaV=1.d-33
            Q=0.d0
            existe=.true.
            return
         end if
      end if
  
      if ( n .eq. 1060 ) then
         if ( .NOT. only_3a ) then
            t9=T/10.**9
            existe=.true.
            goto 1060
         else
            NsigmaV=1.d-33
            Q=0.d0
            existe=.true.
            return
         end if
      end if
  
      if ( n .eq. 1070 ) then
         if ( .NOT. only_3a ) then
            t9=T/10.**9
            existe=.true.
            goto 1070
         else
            NsigmaV=1.d-33
            Q=0.d0
            existe=.true.
            return
         end if
      end if
  
      if ( n .eq. 1080 ) then
         if ( .NOT. only_3a ) then
            t9=T/10.**9
            existe=.true.
            goto 1080
         else
            NsigmaV=1.d-33
            Q=0.d0
            existe=.true.
            return
         end if
      end if
  
      if ( n .eq. 1090 ) then
         if ( .NOT. only_3a ) then
            t9=T/10.**9
            existe=.true.
            goto 1090
         else
            NsigmaV=1.d-33
            Q=0.d0
            existe=.true.
            return
         end if
      end if
  
      if ( n .eq. 1100 ) then
         if ( .NOT. only_3a ) then
            t9=T/10.**9
            existe=.true.
            goto 1100
         else
            NsigmaV=1.d-33
            Q=0.d0
            existe=.true.
            return
         end if
      end if
  
      if ( .NOT. existe ) then
         return
      end if
  
c------------------------------------------
c reaction 1010 : He4(2a,g)C12
  
c He4(a)Be8 
1010  He4aBe8=7.40E+05/t9**(3./2.)*exp(-1.0663/t9)
     +        +4.164E+09/t9**(2./3.)*exp(-13.490/t9**(1./3.)-(t9/0.098)
     +        **2)
     +        *(1.+0.031*t9**(1./3.)+8.009*t9**(2./3.)+1.732*t9+49.883
     +        *t9**(4./3.)+27.426*t9**(5./3.))
  
      Be8agC12=1.30E+02/t9**(3./2.)*exp(-3.3364/t9)
     +         +2.510E+07/t9**(2./3.)*exp(-23.570/t9**(1./3.)-
     +         (t9/0.235)**2.)
     +         *(1.+0.018*t9**(1./3.)+5.249*t9**(2./3.)+0.650*t9+19.176
     +         *t9**(4./3.)+6.034*t9**(5./3.))
     
      NsigmaV=2.90E-16*He4aBe8*Be8agC12
     +        *(0.01+0.2*(1.+4.*exp(-(0.025/t9)**3.263))/
     +        (1.+4.*exp(-(t9/0.025)**9.227)))
     +        +           0.0*
     +        1.35E-07/t9**(3./2.)*exp(-24.811/t9)
  
      if ( t9 .gt. 0.08 ) then
         NsigmaV=2.79E-08/t9**3*exp(-4.4027/t9)
     +           +      0.1*
     +           1.35E-07/t9**(3./2.)*exp(-24.811/t9)
      end if
  
      Q= 14.437
  
c      print*, 'Q=14.437 if C12(a,g)O16 always follows'
  
      return
  
c------------------------------------------
c reaction 1020 : C12(a,g)O16
  
c Taux tel qu'il est dans Caughlan et Fowler 1985 :

 1020 NsigmaV_85=2.93E+08/t9**2./(1.+0.0489/t9**(2./3.))**2.
     +        *exp(-32.120/t9**(1./3.)-(t9/3.496)**2.)
     +        +3.14E+08/t9**2./(1.+0.2654/t9**(2./3.))**2.
     +        *exp(-32.120/t9**(1./3.))
     +        +1.25E+03/t9**(3./2.)*exp(-27.499/t9)+1.43E-02
     +        *t9**5.*exp(-15.541/t9)


c Taux tel qu'il est dans Caughlan et Fowler 1988 :
  
      NsigmaV_88=1.04E+08/t9**2./(1.+0.0489/t9**(2./3.))**2.
     +        *exp(-32.120/t9**(1./3.)-(t9/3.496)**2.   )
     +        +1.76E+08/t9**2/(1.+0.2654/t9**(2./3.))**2*exp(-32.120/t9
     +        **(1./3.))
     +        +1.25E+03/t9**(3./2.)*exp(-27.499/t9)+1.43E-02*t9**5.
     +        *exp(-15.541/t9)

      Q= 7.162

      rapp_8588=NsigmaV_85/NsigmaV_88

      open(unit=88,status='unknown',form='formatted',access='append',
     +     file='rapp_8588.dat')
      write(88,88) log10(t), rapp_8588
 88   format(2d13.4)
      close(88)

      if ( annee_c12ago16 .eq. 1985 ) then

         if ( first ) then
            print*, 'Valeur de 1985 pour C12(a,g)O16 !'
            print*, 'Coefficient multiplicateur : ', coeff_c12ago16
            first=.false.
         end if

         NsigmaV=NsigmaV_85*coeff_c12ago16
  
         return

      else
         if ( annee_c12ago16 .eq. 1988 ) then

            if ( first ) then
               print*, 'Valeur de 1988 pour C12(a,g)O16 !'
               first=.false.
               print*, 'Coefficient multiplicateur : ', coeff_c12ago16
               print*, 'See comment in CFHZ85'
            end if
            
            NsigmaV=NsigmaV_88*coeff_c12ago16

            return

         else
            print*, 'Il faut choisir 1985 ou 1988'
            print*, 'pour la reaction C12(a,g)O16'
            stop
         end if

      end if

c------------------------------------------
c reaction 1030 : O16(a,g)Ne20
  
1030  NsigmaV=9.37E+09/t9**(2./3.)*exp(-39.757/t9**(1./3.)-
     +        (t9/1.586)**2.)
     +        +6.21E+01/t9**(3./2.)*exp(-10.297/t9)+5.38E+02/t9**(3./2.)
     +        *exp(-12.226/t9)
     +        +1.30E+01*t9**2.*exp(-20.093/t9)
      
      Q= 4.734
  
      return
  
  
c------------------------------------------
c reaction 1040 : Ne20(a,g)Mg24
  
1040  NsigmaV=4.11E+11/t9**(2./3.)*exp(-46.766/t9**(1./3.)-
     +        (t9/2.219)**2.)
     +        *(1.+0.009*t9**(1./3.)+0.882*t9**(2./3.)+0.055*t9+0.749
     +        *t9**(4./3.)+0.119*t9**(5./3.))
     +        +5.27E+03/t9**(3./2.)*exp(-15.869/t9)+6.51E+03*t9**(1./2.)
     +        *exp(-16.223/t9)
     +        +                       0.5*
     +        4.21E+01/t9**(3./2.)*exp(-9.115/t9)+3.20E+01/t9**(2./3.)
     +        *exp(-9.383/t9)
  
      gt9=1.+5.*exp(-18.960/t9)
  
      NsigmaV=NsigmaV/gt9
  
      Q= 9.312
  
      return
  
c------------------------------------------
c reaction 1050 : C13(a,n)O16
  
1050  NsigmaV=6.77E+15/t9**(2./3.)*exp(-32.329/t9**(1./3.)-
     +        (t9/1.284)**2.)
     +        *(1.+0.013*t9**(1./3.)+2.04*t9**(2./3.)+0.184*t9)
     +        +3.82E+05/t9**(3./2.)*exp(-9.373/t9)+1.41E+06/t9**(3./2.)
     +        *exp(-11.873/t9)
     +        +2.00E+09/t9**(3./2.)*exp(-20.409/t9)+2.92E+09/t9**(3./2.)
     +        *exp(-29.283/t9)
  
      Q= 2.216
  
      return
  
c------------------------------------------
c reaction 1060 : O17(a,n)Ne20
  
1060  t9a=t9/(1.+0.0268*t9+0.0232*t9**(5./3.)/(1.+0.0268*t9)**(2./3.))
  
      gt9=1.+exp(-10.106/t9)/3.
  
      NsigmaV=1.03E+18/gt9*t9a**(5./6.)/t9**(3./2.)*
     +        exp(-39.914/t9a**(1./3.))
  
      Q= 0.590
  
      return
  
c------------------------------------------
c reaction 1070 : N14(a,g)F18
  
1070  NsigmaV=7.78E+09/t9**(2./3.)*exp(-36.031/t9**(1./3.)-
     +        (t9/0.881)**2.)
     +        *(1.+0.012*t9**(1./3.)+1.45*t9**(2./3.)+0.117*t9+1.97*
     +        t9**(4./3.)+0.406*t9**(5./3.))
     +        +2.36E-10/t9**(3./2.)*exp(-2.798/t9)+2.03E+00/t9**(3./2.)
     +        *exp(-5.054/t9)
     +        +1.15E+04/t9**(2./3.)*exp(-12.310/t9)
  
      Q= 4.415
  
      return
  
c------------------------------------------
c reaction 1080 : O18(a,g)Ne22
  
1080  NsigmaV=1.82E+12/t9**(2./3.)*exp(-40.057/t9**(1./3.)-
     +        (t9/0.343)**2.)
     +        *(1.+0.010*t9**(1./3.)+0.988*t9**(2./3.)+0.072*t9
     +        +3.17*t9**(4./3.)+0.586*t9**(5./3.))
     +        +7.54E+00/t9**(3./2.)*exp(-6.228/t9)+3.48E+01/t9**(3./2.)
     +        *exp(-7.301/t9)
     +        +6.23E+03*t9*exp(-16.987/t9)
     +        +                0.5*
     +        1.00E-11/t9**(3./2.)*exp(-1.994/t9)
  
      Q= 9.669
  
      return
  
c------------------------------------------
c reaction 1090 : Ne22(a,n)Mg25
  
1090  t9a=t9/(1.+0.0548*t9)
      gt9=1.+5.0*exp(-14.781/t9)
      ft9a=exp(-(0.197/t9a)**4.82)
  
      NsigmaV=4.16E+19*ft9a/gt9*t9a**(5./6.)/t9**(3./2.)*exp(-47.004/
     +        t9a**(1./3.))
     +        +1.44E-04/gt9*exp(-5.577/t9)
  
      Q= -0.481
  
      return
     
c------------------------------------------
c reaction 1100 : Ne22(a,g)Mg26
  
1100  t9a=t9/(1.+0.0548*t9)
      gt9=1.+5.0*exp(-14791/t9)
      ft9a=exp(-(0.197/t9a)**4.82)
      fpt9a=exp(-(t9a/0.249)**2.31)
  
      NsigmaV=4.16E+19*fpt9a/gt9*t9a**(5./6.)/t9**(3./2.)*exp(-47.004/
     +        t9a**(1./3.))
  
      Q= 10.612
  
      return
  
  
      end
 
c**********************************************************************
 
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

  
  
  
