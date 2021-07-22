c************************************************************************
 
	subroutine taux_nacre(x,t,ro,r,rt,rro,rx,rxx,q,qt,qro,qx,deriv)
 
c	cycles PP, CNO et 3 alpha
 
c	Taux de reaction thermonucleaires de la database NACRE.

 
c                   ********************************
 
c ATTENTION : cette routine utilise le fichier <<nuc1>> (specifie dans *.don)
c             ou sont tabules les taux de reaction.
 
c                   ********************************
 
c	Auteur:
c              P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice i
c              D. Cordier, ENSCR
 
 
c references :

c Toutes les réactions :
c ----------------------
c @ARTICLE{angulo_etal_99,
c          author  = {C. Angulo and M. Arnould and M. Rayet and et al}, 
c          year    = {1999}, 
c          journal = {Nucl. Phys. A} ,
c          volume  = {656} , 
c          pages   = {3} }

c Exceptions :
c ------------
c  Be7(e-,nu+g)Li7 : Caughlan et Fowler Atomic Data and Nuclear Data Tables, 40, 283, 1988
c  Be8(,a)He4      : Caughlan et Fowler Atomic Data and Nuclear Data Tables, 40, 283, 1988

c entree :
c	x : proportion d'hydrogene en mole /g
c	t : temperature
c	ro : densite
 
c sortie :
c	r(i) : taux de reaction g /s /mole/mole (/mole)
c	rt, rro, rx(i), rxx(i) : derivee /t, /ro, /x, /x2  (due a X dans ecran)
c	q(i): energie ( erg /sec /gr )
c	qt, qro, qx : derivee /t, /ro, /x
 
c initialisation de COMMON
 
c	/reac_nuc/
c	nreac : nombre de reactions
c	t_inf : au dessous pas reactions nulles
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'modele_3.common'
	include 'ctephy.common'
	include 'evol_chim_3.common'

	integer nc_max
	parameter ( nc_max = 100 )

	integer pnt,pmrd
	parameter (pnt=300,pmrd=4)
 
	integer j,i,nrt,knot,l,nt,mrd,pecran(pnreac),nreac_tab,long
 
	logical init,deriv
	real*8 t_t(pnt),r_t(pnt*pnreac),f(pnt*pnreac),lnt9,t_infi,
     1	tt(pnt+pmrd),tr(pnt),drt(pnreac),st(pmrd*pnt*pnreac)
	real*8 t,ro,r(1),rt(1),rro(1),rx(1),q(1),qt(1),qro(1),qx(1),
     1	x,q0(pnreac),ecran,decrx,decr,zeta,decrt,rxx(1),decrxx,
     2	decrro,coef(pnreac)
 
	character table_name*(nc_max)

	external long

	data init/.false./
 
	save
 
2000	format((1x,1p8d10.3))
 
	if(.not. init)then	!en cas de modifications assurer la coherence
	 init=.true.		!avec reac_nuc
 
c Lecture des tables
 
      open(1,form='unformatted',status='old',file=nuc1)
	read(1,err=200) table_name
	if ( table_name(:long(table_name)) .ne. 'NACRE_01' ) then
	   print*, ' '
	   print*, ' ==>> Dans "taux_nacre" les tables de réactions nucléaires ne sont'
           print*, '      pas des "NACRE_01" !'
	   print*, ' '
	   stop
	end if
        read(1,err=100) nt, t_t, r_t
        if( nt .ne. pnt ) then
          print*, 'Probleme avec taux_tab : pnt=',pnt,' et nt=',nt
          stop
        end if
      close(1)

c	 ATTENTION nreac: nb. de reac. utilisees, nreac_tab: nb. reac. tabulees
 
	 nreac_tab=30	!nombre de reactions tabulees
	 nreac=30	!nombre de reactions utilisees
	
	 lnt9=log(1.d9)
	 do i=1,nt
	  t_t(i)=log(t_t(i))+lnt9
c      write(6,'(31A10)') ' ',
c     &      '1', '2', '3', '4', '5', '6', '7', '8', '9','10',
c     &     '11','12','13','14','15','16','17','18','19','20',
c     &     '21','22','23','24','25','26','27','28','29','30'
c	  write(6,'(31d10.3)') t_t(i), (r_t(nt*(j-1)+i),j=1,nreac_tab)
	  do j=1,nreac_tab
	   if(r_t(nt*(j-1)+i) .le. 0.d0)then
	    write(6,*)'r .le. 0; i temp., j reac., r',i,j,r(nt*(j-1)+i)
	    stop
	   endif
	   f(nreac_tab*(i-1)+j)=log(r_t(nt*(j-1)+i))	!r(t,reac), f(reac,t)
	  enddo !i
	 enddo	!j
c	 write(6,2000)(t_t(i),i=1,nt)
c	 write(6,2000)(exp(t_t(i)),i=1,nt)
 
c	 calcul des coefficients d'interpolation des reactions
 
	 mrd=4
	 call pp1dn(nreac_tab,t_t,tt,tr,nt,nrt,mrd,t_t(1),knot,l,drt,f,st,rt,
     1	.false.)
c	 write(6,2000)(tr(i),i=1,nrt)
c	 write(6,2000)(tt(i),i=1,knot)
 
c	 initialisation des energies, des puissances d'ecran, des coefficients
c	 les indices correspondent aux numeros des reactions
 
	 q0(1)=1.442d0		!PP
	 pecran(1)=1
	 coef(1)=1./2.
 
	 q0(2)=5.494d0
	 pecran(2)=1
	 coef(2)=1.
 
	 q0(3)=12.860d0
	 pecran(3)=4
	 coef(3)=1./2.
 
	 q0(4)=1.588d0
	 pecran(4)=4
	 coef(4)=1.
 
	 q0(5)=17.346d0
	 pecran(5)=3
	 coef(5)=1./2.
 
	 q0(6)=0.0497
	 pecran(6)=0
	 coef(6)=1.
 
	 q0(7)=0.137+11.261+0.092
	 pecran(7)=4
	 coef(7)=1.
 
	 q0(8)=1.9435+1.508
	 pecran(8)=6
	 coef(8)=1.
 
	 q0(9)=7.55063
	 pecran(9)=6
	 coef(9)=1.
 
	 q0(10)=7.2971+1.757
	 pecran(10)=7
	 coef(10)=1.
 
	 q0(11)=12.128
	 pecran(11)=7
	 coef(11)=1.
 
	 q0(12)=4.96561
	 pecran(12)=7
	 coef(12)=1.
 
	 q0(13)=0.60035+1.761
	 pecran(13)=8
	 coef(13)=1.
 
	 q0(14)=1.192d0
	 pecran(14)=8
	 coef(14)=1.
 
         q0(15)=0.		!CNO chaud non incorpore
         pecran(15)=7
         coef(15)=1.
 
         q0(16)=5.607+1.76
         pecran(16)=8
         coef(16)=1.
 
         q0(17)=7.994d0
         pecran(17)=8
         coef(17)=1.
 
         q0(18)=8.114d0
         pecran(18)=9
         coef(18)=1.
 
         q0(19)=12.843d0
         pecran(19)=9
         coef(19)=1.
 
         q0(20)=3.980d0
         pecran(20)=8
         coef(20)=1.
 
	 q0(21)=7.274d0	!3 alpha
	 pecran(21)=8	!facteur ro**2
	 coef(21)=1./6.
 
	 q0(22)=7.162d0
	 pecran(22)=12
	 coef(22)=1.
 
	 q0(23)=4.730d0
	 pecran(23)=16
	 coef(23)=1.
 
         q0(24)=9.316d0
         pecran(24)=20
         coef(24)=1.
 
         q0(25)=2.216d0
         pecran(25)=12
         coef(25)=1.
 
         q0(26)=0.586d0
         pecran(26)=16
         coef(26)=1.
 
         q0(27)=4.415d0+1.76d0		!1.76 pour beta+
         pecran(27)=14
         coef(27)=1.
 
         q0(28)=9.667d0
         pecran(28)=16
         coef(28)=1.
 
         q0(29)=-0.481d0
         pecran(29)=20
         coef(29)=1.
 
         q0(30)=10.615d0
         pecran(30)=20
         coef(30)=1.
 
 
	 do i=1,nreac
	  q0(i)=q0(i)*eve*1.d6/amu	!energie des reactions en erg/reaction
	 enddo
c	 write(6,2000)q0
 
	 t_inf=exp(t_t(1))		!temperature minimale pour les reactions
	 t_infi=t_inf			!tinfi: t_inf interne a la routine
c	 write(6,2000)t,t_inf
c	 write(6,*)'fin initialisation taux_reac'
c	 pause'taux_reac'
 
	endif
 
c	write(6,2000)t,t_inf,t_infi
c	pause'taux_reac'
	if(t .lt. t_infi)then
	 do i=1,nreac
	  r(i)=1.d-100
	  rt(i)=0.
	  rro(i)=0.
	  rx(i)=0.
	  q(i)=0.
	  qt(i)=0.
	  qro(i)=0.
	  qx(i)=0.
	 enddo
	 return
	endif
 
c--------------------------------------------------------------------------
c	effet d'ecran
 
	zeta=(x+3.)/2.
	decr=.188*sqrt(ro*zeta/(t*1.d-6)**3)
	ecran=exp(decr)
	if(deriv)then
	 decrt=-3./2.*decr/t			!a ecran pres
	 decrro=decr/ro/2.			!a ecran pres
	 decrx=decr/zeta/2./2.			!a ecran pres
	 decrxx=decrx/4./zeta*(decr-1./zeta)
c	 write(6,*)'t, ro, x, ecran, decrx'
c	 write(6,2000)t, ro, x, ecran, decrx
	endif
 
c--------------------------------------------------------------------------
c	interpolation en ln t
 
	call pp1dn(nreac_tab,t_t,tt,tr,nt,nrt,mrd,log(t),knot,l,drt,f,st,r,
     1	.true.)
 
c        do i=1, nreac
c        print*, 'r(', i, ')', r(i)
c        end do
 
c	reactions et derivees
 
	do i=1,nreac
	 r(i)=exp(r(i))*ro*ecran**pecran(i)*coef(i)
	 if(deriv)then
	  rt(i)=r(i)*(drt(i)/t+pecran(i)*decrt)
	  rro(i)=r(i)*(pecran(i)*decrro+1./ro)
	  rx(i)=r(i)*pecran(i)*decrx
	  rxx(i)=r(i)*pecran(i)*decrxx
	 endif
c	 write(6,2000)r(i)
	enddo	!i
 
c	write(6,*)'t,ro,r(15)'
c	write(6,2000)t,ro,r(15)
 
	r(21)=r(21)*ro
	if(deriv)then
	 rt(21)=rt(21)*ro
	 rro(21)=r(21)*(pecran(21)*decrro+2./ro)	!ro**2 pour reaction 15
	 rx(21)=rx(21)*ro
	 rxx(21)=rxx(21)*ro
	endif
 
	do i=1,nreac
	 q(i)=r(i)*q0(i)
	 if(deriv)then
	  qt(i)=rt(i)*q0(i)
	  qro(i)=rro(i)*q0(i)
	  qx(i)=rx(i)*q0(i)
	 endif
	enddo
 
	return
 
100     print*, ' '
        print*, 'Erreur de lecture de le fichier <<nuc1>> !'
        stop

200	print*, ' '
        print*, 'Dans le fichier <<nuc1>> il n''y a probablement des'
        print*, 'taux NACRE_01 !'
        stop

	end
 
