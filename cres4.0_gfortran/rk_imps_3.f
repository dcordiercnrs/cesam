 
c**************************************************************************
 
	subroutine rk_imps_3(t__,ro__,dt__,t_,ro_,dt_,
     1	t_t,ro_t,compx,t,ro,compy,dt,esti,ok,nuc,kk,dm)
 
c	integration de la composition chimique par divers schemas RKi simplifiee
c	implementation selon Hairer et Wanner p.128
c	contrairement a rk_imp il n'y a pas de test de precision
c	par comparaison avec un second schema
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	CESAM: version 3
 
c entrees
c	t__, t_, t_t, t : temperatures aux temps t-(dt__+dt_), t_dt_, t et t+dt
c	ro__, ro_, ro_t, ro : densites aux temps t-(dt__+dt_), t_dt_, t et t+dt
c	compx : composition chimique au temps t
c	dt__, dt_, dt: pas temporel
c	dm : epaisseur en masse des couches des ZC
c	kk : nombre de couches a melanger
c	rk_n: indice de la formule de Runge Kutta implicite
 
c entrees / sorties
c	compy : composition chimique au temps t+dt
c	ok = .true. : nouveau modele / la precision est atteinte
 
c sorties
c	esti : precision estimee par element
 
c externe
c	nuc : calcul des reactions thermonucleaires
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'modele_3.common'
	include 'evol_chim_3.common'
 
	integer ps,pns,pkk
	parameter (ps=6,pkk=pn, pns=ps*pnelem)
 
	integer i,tourmax,tour,j,kk,s,ns,l,indice,m,rk_n,k,nt
 
	real*8	esti(1),dt,t(1),ro(1),t_t(1),ro_t(1),corr(pns),
     1	compy(1),dcomp(pnelem),ti(ps*pkk),roi(ps*pkk),zi(pns),
     2	e(5),et,ero,hh,be7,b8,n13,o15,f17,epsi,ex(pnelem),compx(1),
     3	jac(pnelem*pnelem),dm(1),jacob((pns)**2),corr_max,
     4	scale(pnelem),dt_,dcomps(pnelem),jacs(pnelem*pnelem),
     5	compx0(pnelem),a(ps*ps),c(ps),
     6	dt__,t_(1),ro_(1),t__(1),ro__(1)
 
	real*8	neuville,x(0:ps),f(0:ps)
 
	logical init,ok
c	data init/.true./
 
	external nuc
 
	data init/.true./
 
	save init,a,c,s,ns
 
2000	format((1x,1p8d10.3))
2001	format((1x,1p9d8.1))
2003	format(1x,1p12d8.1)
 
	if(init)then
	 init=.false.
	 rk_n=max(1,ordre)
	 if(rk_n .eq. 1)then
	  print*,'integration par Euler implicite d''ordre 1'
	  s=1		!p pour precision
	  a(1)=1.d0
	  c(1)=1.d0
 
	 elseif(rk_n .eq. 2)then
	  print*,'integration par Lobatto IIIc d''ordre 2'
	  s=2
	  a(1)=0.5d0
	  a(2)=0.5d0
	  a(3)=-0.5d0
	  a(4)=0.5d0
	  c(1)=0.d0
	  c(2)=1.d0
 
	 elseif(rk_n .eq. 3)then
	  print*,'integration par Radau IIa d''ordre 3'
	  s=2
	  a(1)=5.d0/12.d0
	  a(2)=0.75d0
	  a(3)=-1.d0/12.d0
	  a(4)=0.25d0
	  c(1)=1.d0/3.d0
	  c(2)=1.d0
 
	 elseif(rk_n .eq. 4)then
	  print*,'integration par Lobatto IIIc d''ordre 4'
	  s=3
	  a(1)=1.d0/6.d0
	  a(2)=1.d0/6.d0
	  a(3)=1.d0/6.d0
	  a(4)=-1.d0/3.d0
	  a(5)=5.d0/12.d0
	  a(6)=2.d0/3.d0
	  a(7)=1.d0/6.d0
	  a(8)=-1.d0/12.d0
	  a(9)=1.d0/6.d0
	  c(1)=0.d0
	  c(2)=0.5d0
	  c(3)=1.d0
 
	 elseif(rk_n .eq. 5)then
	  print*,'integration par Radau IIa d''ordre 5'
	  s=3
	  a(1)=(88.d0-7.d0*sqrt(6.d0))/360.d0
	  a(2)=(296.d0+169.d0*sqrt(6.d0))/1800.d0
	  a(3)=(16.d0-sqrt(6.d0))/36.d0
	  a(4)=(296.d0-169.d0*sqrt(6.d0))/1800.d0
	  a(5)=(88.d0+7.d0*sqrt(6.d0))/360.d0
	  a(6)=(16.d0+sqrt(6.d0))/36.d0
	  a(7)=(-2.d0+3.d0*sqrt(6.d0))/225.d0
	  a(8)=(-2.d0-3.d0*sqrt(6.d0))/225.d0
	  a(9)=1.d0/9.d0
	  c(1)=(4.d0-sqrt(6.d0))/10.d0
	  c(2)=(4.d0+sqrt(6.d0))/10.d0
	  c(3)=1.d0
 
	 elseif(rk_n .eq. 6)then
	  print*,'integration par Lobatto IIIc d''ordre 6'
	  s=4
	  a(1)=1.d0/12.d0
	  a(2)=1.d0/12.d0
	  a(3)=1.d0/12.d0
	  a(4)=1.d0/12.d0
	  a(5)=-sqrt(5.d0)/12.d0
	  a(6)=1.d0/4.d0
	  a(7)=(10.d0+7.d0*sqrt(5.d0))/60.d0
	  a(8)=5.d0/12.d0
	  a(9)=sqrt(5.d0)/12.d0
	  a(10)=(10.d0-7.d0*sqrt(5.d0))/60.d0
	  a(11)=1.d0/4.d0
	  a(12)=5.d0/12.d0
	  a(13)=-1.d0/12.d0
	  a(14)=sqrt(5.d0)/60.d0
	  a(15)=-sqrt(5.d0)/60.d0
	  a(16)=1.d0/12.d0
	  c(1)=0.d0
	  c(2)=(5.d0-sqrt(5.d0))/10.d0
	  c(3)=(5.d0+sqrt(5.d0))/10.d0
	  c(4)=1.d0
	 else
	  print*,'il n''y a pas de formule pour rk_n=',rk_n
	  stop
	 endif
 
	 ns=nbelem*s
 
	 tourmax=30	!nombre de tours max pour N-R
	 epsi=1.d-6
	 do i=1,nbelem	!abondances de reference
	  scale(i)=ab_min(i)*100.
	 enddo
 
	endif
 
c	valeur moyenne en t
 
c	print*,'kk,dm',kk
c	write(6,2000)(dm(k),k=1,kk)
 
	do l=1,nbelem		!initialisation de somme Xi sur la ZM
	 compx0(l)=0.d0
	 do k=1,kk		!comp(element,ZM)
	  compx0(l)=compx0(l)+compx(nbelem*(k-1)+l)*dm(k)
	 enddo
c	 print*,'compx0',compx0(l)
	enddo
 
c	interpolation des temperatures et densites
 
	x(0)=1.d0		!les abscisses
	x(1)=0.d0
	nt=1			!ordre d'interpolation
	if(rk_n .gt. 2)then	!les points intermediaire ne sont
	 if(dt_ .gt. 0.d0)then	!necessaires que si l'ordre > 2
	  x(2)=-dt_/dt
	  nt=2
	 endif
	 if(dt__ .gt. 0.d0)then
	  x(3)=x(2)-dt__/dt
	  nt=3
	 endif
	endif
	nt=min(nt,s)
	if(t(1) .eq. t_t(1))nt=1		!interpolation lineaire
 
c	les valeurs de T et de ro aux points intermediaires par interpolation
c	ti(kk,s)=ti(point de ZM,etape RK)
 
c	print*,'ordre superieur',sp,kk,nt
	do i=1,s
	 do k=1,kk
	  f(0)=t(k)
	  f(1)=t_t(k)
	  if(nt .ge. 2)f(2)=t_(k)
	  if(nt .eq. 3)f(3)=t__(k)
	  ti(kk*(i-1)+k) =neuville(c(i),x,f,nt)
	  f(0)=ro(k)
	  f(1)=ro_t(k)
	  if(nt .ge. 2)f(2)=ro_(k)
	  if(nt .eq. 3)f(3)=ro__(k)
	  roi(kk*(i-1)+k) =neuville(c(i),x,f,nt)
c	  print*,'i, c, t, ro',i,c(i),k,ti(kk*(i-1)+k),roi(kk*(i-1)+k)
	 enddo
c	 write(6,2000)(ti(kk*(i-1)+k),k=1,kk)
c	 write(6,2000)(roi(kk*(i-1)+k),k=1,kk)
	enddo
 
c	initialisation du vecteur Zi: z(l,s), z(variable,etape RK)
 
	do i=1,ns
	 zi(i)=0
	enddo
c	write(6,2000)(zi(i),i=1,ns)
c	pause'Z initial'
	
c	iterations NR
 
	tour=0
	corr_max=1.d5
	do while (tour .lt. tourmax .and. corr_max. gt. epsi)
	 if(tour .gt. tourmax)then
	  write(6,115)rk_n,t(1),ro(1),compx(1),dt
115	  format(1x,'pas de conv. de N_R dans RK_IMPS_3 numero:',i3,/,
     1	1x,'t=',1pd10.3,' ro=',1pd10.3,' X=',1pd10.3,' dt=',
     2	1pd10.3)
	  write(6,*)'x / y / esti'
	  write(6,2001)(compx(i),i=1,nbelem)
	  write(6,2001)(compy(i),i=1,nbelem)
	  write(6,2001)(esti(i),i=1,nbelem)
	  write(6,*)' '
	  ok=.false.
c	  pause'ne marche pas'
	  if(dt/2. .gt. dtmin)return
	  write(6,*)'le dt < dtmin, abandon, pour ordre sup. avec RK_IMPS_3 ',
     1	rk_n
	  stop
	 endif
 
	 do i=1,(ns)**2	!mise a 0 du JACOBIEN
	  jacob(i)=0.d0
	 enddo
 
	 do i=1,ns			!initialisations
	  corr(i)=zi(i)			!second membre
	  jacob(ns*(i-1)+i)=1.d0	!diagonale JACOBIEN
	 enddo
	 do j=1,s		!pour chaque colonne
	  do l=1,nbelem		!pour chaque inconnue
	   compy(l)=compx0(l)+zi(nbelem*(j-1)+l)	!compy = y0+zi VT
	  enddo
c	  print*,'j,zi',j,(zi(nbelem*(j-1)+l),l=1,nbelem)
 
	  do l=1,nbelem		!initialisation
	   dcomps(l)=0.d0
	  enddo
	  do l=1,nbelem*nbelem
	   jacs(l)=0.d0
	  enddo
 
	  do k=1,kk	!pour chaque couche les reactions thermonucleaires
c	   write(6,2000)ti(kk*(j-1)+k),roi(kk*(j-1)+k),(compy(l),l=1,nbelem)
	   call nuc(ti(kk*(j-1)+k),roi(kk*(j-1)+k),compy,dcomp,jac,.true.,2,
     1	e,et,ero,ex,hh,be7,b8,n13,o15,f17)	!F ( x+Ci h,y0+zi )
c	   print*,'j,(z,dcomp)',j,(compy(l),dcomp(l),l=1,nbelem)
	   do l=1,nbelem
	    dcomps(l)=dcomps(l)+dcomp(l)*dm(k)*dt
	   enddo	!l
	   do l=1,nbelem**2
	    jacs(l)=jacs(l)+jac(l)*dm(k)*dt
	   enddo	!l
	  enddo		!k
 
c	  jacob(variable,etape,variable,etape)
c	  2d corr(variable,etape)
	
	  do l=1,nbelem
	   do i=1,s	!zi - A*fi*h	pour chaque etape RK
	    corr(nbelem*(i-1)+l)=corr(nbelem*(i-1)+l)-a(s*(j-1)+i)*dcomps(l)
c	    print*,'l,i,j,corr',l,i,j,corr(nbelem*(i-1)+l)
	    do m=1,nbelem	!contribution au jacobien
	     indice=nbelem*(s*(nbelem*(j-1)+m-1)+i-1)+l
	     jacob(indice)=jacob(indice)-jacs(nbelem*(m-1)+l)*a(s*(j-1)+i)
	    enddo		!m
	   enddo		!i
	  enddo			!l
	 enddo			!j des etapes RK
 
c	 do i=1,ns
c	  write(6,2000)(jacob(ns*(j-1)+i),j=1,ns),corr(i)
c	 enddo
c	 pause'avant simq'
 
	 call simq(jacob,corr,ns)
c	 print*,'correction',tour
c	 do i=1,s
c	  write(6,2000)(corr(nbelem*(i-1)+l),l=1,nbelem)
c	 enddo
c	 pause'iteration'
	 corr_max=0
	 do l=1,nbelem
	  do i=1,s
	   if(abs(compx0(l)) .gt. ab_min(l))corr_max=max(corr_max,
     1	abs(corr(nbelem*(i-1)+l)/scale(l)))
	   zi(nbelem*(i-1)+l)=zi(nbelem*(i-1)+l)-corr(nbelem*(i-1)+l)
	  enddo
	 enddo
c	 print*,'corr_max',corr_max
c	 write(6,2002)(corr(i),i=1,ns)
2002	 format(1x,'corr',1p12d8.1)
c	 write(6,2001)(zi(i),i=1,ns)
 
	 tour=tour+1
c	 pause'les zi'
	enddo		!second while
 
c	print*,'tour',tour,rk_n,ordre
c	write(6,2000),t(1),ro(1)
c	pause
 
c	solution et estimation de la variation relative
 
	do l=1,nbelem
	 compy(l)=compx0(l)+zi(nbelem*(s-1)+l)
	 esti(l)=abs(zi(nbelem*(s-1)+l))/scale(l)
	enddo
 
c	print*,'solution'
c	write(6,2000)(compy(l),l=1,nbelem)
c	pause'solution'
 
	ok=.true.
	
c	write(6,*)'sortie RK_IMPS_3, kk, ok/compx/compy',kk,ok
c	write(6,2003)(compx0(l),l=1,nbelem)
c	write(6,2003)(compy(l),l=1,nbelem)
c	write(6,2003)(zi(nbelem*(s-1)+l),l=1,nbelem)
c	write(6,2003)(esti(l),l=1,nbelem)
c	pause
 
	return
 
	end
 
