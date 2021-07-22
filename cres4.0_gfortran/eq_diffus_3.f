 
c***********************************************************************
 
	 subroutine eq_diffus_3(y,cx,aes,aed,bes,bed,fait,mstar,
     1	bp,q,qt,knot,n,r2,m23,mc_t,mct_t,nc_t,knotc_t,chim_t,dt,
     2	a_mix,nmix,mix,a_mixt,knot_mix,
     3	old_m23,new_m23,nm,new_m23t,knotm,	
     4	etat,opa,conv,nuc,coeff_diff)
 
c	formation des equations pour diffus.for
c	methode de Galerkin ou de Petrov-Galerkin
c	les derivees/xchim sont calculees dans la routine de calcul
c	du coefficient de diffusion
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	CESAM, version 3
 
c entrees
c	y: solution
c	cx: point de calcul en m**2/3
c	dt: pas temporel
c	bp,r2,m23,q,qt,knot,n: modele au temps t+dt
c	mc_t,mct_t,nc_t,knotc_t,chim_t(1): int. comp. chim. temps t
c	a_mix,nmix,mix,a_mixt,knot_mix: interpolation de la nature de la zone
c	fait=1: point interieur, fait=2: point exterieur (perte de masse)
c	mstar: masse temps t+dt, avec perte de masse
c	old_m23,new_m23,nm,new_m23t,knotm : interpolation m(t+dt)-->m(t)
 
c sorties	
c	indice s/d: coefficients de la spline/derivee
c	aes, aed: coefficients de la spline/derivee
c	bes, bed: seconds membres
 
c external
c	nuc: reactions thermonucleaires
c	coeff_diff: calcul du coefficient de diffusion
c	etat: equation d'etat
c	opa: opacite
c	conv: routine de convection
 
c_______________________________________________________________________
 
c	pour chaque spline Nj
c	le systeme est constitue de nbelem produits scalaires
 
c	les equations:
 
c	< 3/2 nu^1/2 [ Xi(t+1)-Xi(t) - Psi dt ] , Fi_j > +
c	+ < 32/(3 nu^1/2)pi^2r^4ro^2 dt[D1i dX1/dnu + D2i dXi/dnu],dFi_j/dnu >-
c	- < 4 pi r^2 ro dt Vi Xi , dFi_j/dnu >=0
 
c	d(1,i): coefficient de dX1/dnu pour l'element i, d(1,1)=0
c	d(2,i): coefficient de dXi/dnu pour l'element i, ie diffusion turbulente
c	+diffusion atomique element i + diffusion de melange (Dt+Di+Dm)
c	V: vitesse de diffusion
c	nu: m**2/3 = cx
 
c	Dt, D1i, D2i, Dm et Vi sont calcules par coeff_diff,
c	dependent de l'element Xi, Xdot est calcule par nuc
 
c	les coefficients de diffusion microscopiques D et V dependent de
c	Xi et de X1. Au lieu d'avoir une matrice pour D on utilise deux vecteurs
c	d(1,i): coefficient de dX1/dnu pour l'element i,
c	d(2,i): coefficient de dXi/dnu pour l'element i
c	afin ne ne pas compter 2 fois d(1,1)=d(2,1) on impose d(1,1)=0
 
c_______________________________________________________________________
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'ctephy.common'
	include 'modele_3.common'
	include 'evol_chim_3.common'
	
	integer	i,knotc_t,lx,n,knot,nc_t,nmix,knot_mix,fait,j,k,kd,nm,knotm
 
	real*8	aed(1),aes(1),bed(1),bes(1),y(1),d(2,pnelem),cx,mn,
     1	dcomp(pnelem),p,t,r,l,dtp,mstar,old_m23(1),new_m23(1),
     2	jac(pnelem**2),epsilon(5),a_mix(1),mix(1),a_mixt(1),
     4	dt,chim_t(1),mc_t(1),mct_t(1),cte11,cte21,cte31,new_m23t(1),
     5	comp_t(pnelem),f(pne),dfdx(pne),dcomp_t(pnelem),
     6	r2(1),m23(1),bp(1),q(1),qt(1),ro,drop,drot,drox,drott,
     7	drotp,drotx,xchim(pnelem),v(pnelem),kap,dkapt,dkapro,dkapx,
     8	cte1,cte2,cte3,ddx(3,pnelem),dvdx(pnelem),nuc_m(pnelem),
     9	cte4,u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb,mstarp
	
	real*8	deltax,dcpx,dgradadx,dgradradx,stor,unpdx,dstor,
     1	bess(pnelem),beds(pnelem)
c	data unpdx/1.000001d0/,dtp/-100.d0/,cte1/1.5d0/,mstarp/-1.d15/
 
	real*8	depst,depsro,depsx(pnelem),hhe,be7e,
     1	b8e,n13e,o15e,f17e,delta,cpp,cp,gradad,gradrad,
     2	cte80,cte9,cte8
 
	external coeff_diff,etat,opa,conv,nuc
 
	logical melange,init
	data init/.true./
	data unpdx/1.000001d0/,dtp/-100.d0/,cte1/1.5d0/,mstarp/-1.d15/
 
	save cte1,cte2,cte3,cte4,dtp,mstarp,init
 
c________________________________________________________________________
 
c	le nombre d'inconnues est nbelem
 
c	y(var,der):
c	var=1... nbelem pour les Xi
c	var=1+nbelem, 2+nbelem... 2*nbelem pour les Fi
 
c	der=1 pour la fonction
c	der=2 pour derivee premiere, etc..
 
c  	y(nbelem*(der-1)+var) := derivee d'ordre der(1 pour 0), (2 pour ')
c	de la variable var pour X: 1, .. nbelem, pour F: 1+nbelem, .. 2*nbelem
 
2000	format(1x,1p8d10.3)
2001	format(1x,1p10d8.1)
 
c	write(6,*)fait,cx
c	write(6,20000)dt,dtp
c	write(6,2001)(y(nbelem*(1-1)+i),i=1,nbelem)
c	pause'entree eq_diffus'
 
	if(init)then
	 init=.false.
	 cte9=3./16./pi/aradia/clight/g
	 cte80=cte9*lsol/msol		!calcul direct de gradrad
	endif
	
	goto(100,200),fait
 
100	if(dt .ne. dtp .or. mstar .ne. mstarp)then
	 dtp=dt
	 mstarp=mstar
	 cte2=2.d0/3.d0*(4.d0*pi*rsol**2/msol)**2*secon6*dt
	 cte3=4.d0*pi*rsol**2/msol*secon6*dt
	 cte4=mdot*1.d6*dt		!perte de masse
	 cte8=cte80/mstar		!calcul direct de gradrad	
c	 write(6,2000)cte1,cte2,cte3,cte4
	endif
 
c	M, P, T, R, L au point mu=m**2/3  (m**2/3=cx)
 
	mn=sqrt(cx)**3
c	prnt*,'mu',cx
 
	call inter_3('mu',bp,q,qt,n,knot,min(cx,m23(n)),f,dfdx,r2,m23)
 
	p=exp(f(1))
	t=exp(f(2))
	r=sqrt(abs(f(3)))
c	print*,'r',r	
	l=sqrt(abs(f(4)))**3
c	write(6,2000)p,t,r,l,mn
 
c	est-on dans une zone melangee ?
 
	u=min(cx,a_mix(nmix))	!u: vT
	call sbsp1dn(1,mix,a_mix,a_mixt,nmix,2,knot_mix,.true.,u,lx,f,dfdx)
	melange=f(1) .lt. 0.d0
c	print*,cx,f(1),melange
 
c	la composition chimique au temps t
 
	u=min(cx,new_m23(nm))	!u:vT
	call sbsp1dn(1,old_m23,new_m23,new_m23t,nm,2,knotm,.true.,
     1	u,lx,f(1),dfdx)	!masse au temps t f(1)
	f(1)=min(f(1),mc_t(nc_t))
	call sbsp1dn(nbelem,chim_t,mc_t,mct_t,nc_t,m_ch,
     1	knotc_t,.true.,f(1),lx,comp_t,dcomp_t)
	if(iw .gt. 1)comp_t(iw)=comp_t(iw)*cx/f(1)	!perte de Mw
c	write(6,2000)(comp_t(i),i=1,nbelem),(y(i),i=1,nbelem)
	
c	la composition chimique par gramme
 
	do i=1,nbelem
	 xchim(i)=y(i)
	enddo
	call chim_gram_3(xchim,dcomp,nuc_m)	!dcomp VT
 
c	la densite, le gradient adiabatique et leurs derivees/X
 
	call etat(p,t,xchim,.false.,
     1       ro,drop,drot,drox,drott,drotp,drotx,
     2       u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
 
	delta=-t/ro*drot
	deltax=delta*(-drox/ro+drotx/drot)
			
	cpp=p/ro/t*delta
	dcpx=dutx+cpp*(-drox/ro+deltax/delta)
	cp=dut+cpp
	
	gradad=p/ro*delta/cp/t			!gradient adiabatique
	dgradadx=gradad*(-drox/ro+deltax/delta-dcpx/cp)	
	 		
c	write(6,2000)ro,gradad
	
c	l'opacite et le gradient radiatif
 
	call opa(xchim(1),t,ro,kap,dkapt,dkapro,dkapx)
	
	gradrad=cte8*l*p*kap/mn/t**4
	dgradradx=gradrad/kap*dkapx
	
c	write(6,2000)ro,gradad,gradrad
	
c	le coefficient de diffusion
 
	call coeff_diff(melange,p,t,r,l,mn,ro,drox,kap,dkapx,
     1	gradrad,dgradradx,gradad,dgradadx,
     2	xchim,mstar,d,ddx,v,dvdx)
c	print*,'d1i,d2i,vi'
c	do i=1,nbelem	
c	 write(6,2000)d(1,i),d(2,i),v(i),ddx(1,i),ddx(2,i),ddx(3,i),dvdx(i)
c	enddo
	
c	les reactions thermonucleaires Xdot
 
	call nuc(t,ro,y,dcomp,jac,.true.,2,
     1	epsilon,depst,depsro,depsx,hhe,be7e,b8e,n13e,o15e,f17e)
 
c	do i=1,nbelem	
c	 write(6,2000)d(1,i),d(2,i),v(i),ddx(1,i),ddx(2,i),ddx(3,i),dvdx(i)
c	enddo
c	write(6,2000)ro,(dcomp(i),i=1,nbelem)
 
c	les equations:
 
c	< 3/2 nu^1/2 [ Xi(t+1)-Xi(t) - Psi dt ] , Fi_j > +
c	+ < 32/(3 nu^1/2) pi^2 r^4 ro^2 dt
c	[D1i dX1/dnu+D2i dXi/dnu],dFi_j/dnu >
c	- < 4 pi r^2 ro dt Vi Xi , dFi_j/dnu >=0
 
	cte11=cte1*sqrt(cx)
	cte21=cte2*r**4*ro**2/sqrt(cx)
	cte31=-cte3*r**2*ro
	drox=drox/ro
 
c	pour bes(spline):	3/2 nu^1/2 [ Xi(t+1)-Xi(t) - Psi dt ]
c	pour bed(derivee de la spline): 32/(3 nu^1/2) pi^2 r^4 ro^2 dt *
c				*  (Dt+Di+Dm) dXi/dnu  -4 pi r^2 ro dt Vi Xi
 
c	bes : coefficient de la spline du produit scalaire
c	bed : coefficient de la derivee de la spline du produit scalaire
	
	do i=1,nbelem	!d(1,1)=0 pour ne pas etre compte 2 fois
	 bes(i)=cte11*(y(i)-comp_t(i)-dt*dcomp(i))
	 bed(i)=cte31*v(i)*y(nbelem*(1-1)+i)	!coeff. de Xi: vitesse
     1	+cte21*(d(1,i)*y(nbelem*(2-1)+1)+d(2,i)*y(nbelem*(2-1)+i))
	
c	 derivees
 
c  	 y(nbelem*(der-1)+var) := derivee d'ordre der(1 pour 0), (2 pour ')
 
c	 ae(nbelem,nbelem,2)=
c	 ae(equation,variable,derivee+1)=ae(nbelem*(nbelem*(k-1)+j-1)+i)
c	 coefficient de la (k-1)-ieme derivee de la j-ieme variable de la i-ieme
c	 equation
 
c	 aes : derivees/ Xi, Xi' des coeff. de la spline du produit scalaire
c	 aed : derivees/ Xi, Xi' des coeff. de la der. du produit scalaire
 
c	 d(1,i) : coeff de d X1 pour element i
c	 d(2,i) : coeff de d Xi pour element i
c	 ddx(1,i) : derivee de d(1,i) / X1
c	 ddx(2,i) : derivee de d(2,i) / X1
c	 ddx(3,i) : derivee de d(2,i) / Xi (nul, sauf pour He4)
c	 v(i) : coeff de Xi
c	 dvdx(i) : derivee de v / X1
 
c	 derivee dans bed / Xi et derivee dans bed de d et de v / Xi
	
	 aed(nbelem*(nbelem*(1-1)+i-1)+i)=cte31*v(i)
     1	+cte21*ddx(3,i)*y(nbelem*(2-1)+i)
	
c	 derivee dans bed / X1' (d(1,1)=0 pour ne pas compter 2fois)
	
	 aed(nbelem*(nbelem*(2-1)+1-1)+i)=cte21*d(1,i)
	
c	 derivee dans bed / Xi'
	
	 aed(nbelem*(nbelem*(2-1)+i-1)+i)=cte21*d(2,i)
	
c	 derivee dans bed de d et de v / X1 et derivee dans bed / ro puis /X1
 
	 aed(nbelem*(nbelem*(1-1)+1-1)+i)=aed(nbelem*(nbelem*(1-1)+1-1)+i)
     1	+cte31*y(nbelem*(1-1)+i)*(dvdx(i)+v(i)*drox)
     2	+cte21*(y(nbelem*(2-1)+1)*(ddx(1,i)+d(1,i)*2.d0*drox)
     3	+y(nbelem*(2-1)+i)*(ddx(2,i)+d(2,i)*drox))
	
c	 derivee dans bes de y(i) / Xi (la diagonale)
 
	 aes(nbelem*(nbelem*(1-1)+i-1)+i)=cte11
	 	
	 do j=1,nbelem
	
c	  derivee dans bes de dcomp(i) / Xj (le jacobien)
 
	  aes(nbelem*(nbelem*(1-1)+j-1)+i)=aes(nbelem*(nbelem*(1-1)+j-1)+i)
     1	-cte11*dt*jac(nbelem*(j-1)+i)
	 enddo	
 
	enddo
	
c	write(6,2001)(bes(i),i=1,nbelem)
c	write(6,2001)(bed(i),i=1,nbelem)	
c	do kd=1,2
c	 print*,kd
c	 do k=1,nbelem
c	  write(6,2001)(aes(nbelem*(nbelem*(kd-1)+k-1)+i),i=1,nbelem)
c	  print*
c	  write(6,2001)(aed(nbelem*(nbelem*(kd-1)+k-1)+i),i=1,nbelem)
c	 enddo
c	enddo
c	pause	
 
c	verification des derivees
 
	if(.false.)then
c	if(.true.)then
c	if(.not.melange)then
 
	 print*,'nbelem=',nbelem
	 print*,'p,t,r,l,mn,cx,ro,cte11'	
	 write(6,2000)p,t,r,l,mn,cx,ro,cte11
	 print*,'d1i,d2i,vi,Xi,dXi'
	 do i=1,nbelem
	  write(6,2000)d(1,i),d(2,i),v(i),y(nbelem*(1-1)+i),y(nbelem*(2-1)+i),
     1	comp_t(i)
	 enddo
	 print*,'ddx(1,i),ddx(2,i),ddx(3,i),dvdx(i),dvdx(i)'
	 do i=1,nbelem
	  write(6,2000)ddx(1,i),ddx(2,i),ddx(3,i),dvdx(i)
	 enddo
	
c	 les equations et les derivees numeriques en Xi
 
	 print*,'bes/bed initiaux'
	 write(6,2000)(bes(i),i=1,nbelem)
	 write(6,2000)(bed(i),i=1,nbelem)
	
	 do kd=1,2		!derivees 0 et '
	  do k=0,nbelem		!0 pour la fonction
	
	   if(k .ne. 0)then
	    stor=y(nbelem*(kd-1)+k)
	    if(abs(stor) .lt. 1.d-10)then
	     y(nbelem*(kd-1)+k)=1.d-10
	    else
	     y(nbelem*(kd-1)+k)=stor*unpdx
	    endif
	    dstor=y(nbelem*(kd-1)+k)-stor
	   endif
 
c	   la composition chimique par gramme
 
	   do i=1,nbelem
	    xchim(i)=y(i)
	   enddo
	   call chim_gram_3(xchim,dcomp,nuc_m)	!dcomp VT
 
c	   la densite et le gradient adiabatique
 
	   call etat(p,t,xchim,.false.,
     1       ro,drop,drot,drox,drott,drotp,drotx,
     2       u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
 
	   delta=-t/ro*drot
	   cpp=p/ro/t*delta
	   cp=dut+cpp
	   gradad=p/ro*delta/cp/t			!gradient adiabatique
	 	
c	   l'opacite et le gradient radiatif
 
	   call opa(xchim(1),t,ro,kap,dkapt,dkapro,dkapx)
	   gradrad=cte8*l*p*kap/mn/t**4
	
c	   le coefficient de diffusion
 
	   call coeff_diff(melange,p,t,r,l,mn,ro,drox,kap,dkapx,
     1	gradrad,dgradradx,gradad,dgradadx,
     2	xchim,mstar,d,ddx,v,dvdx)
 
c	   les reactions thermonucleaires Xdot
 
	   call nuc(t,ro,y,dcomp,jac,.false.,2,
     1	epsilon,depst,depsro,depsx,hhe,be7e,b8e,n13e,o15e,f17e)
 
c	   write(6,2000)ro,(d(i),i=1,nbelem),(dcomp(i),i=1,nbelem)
 
	   cte21=cte2*r**4*ro**2/sqrt(cx)
	   cte31=-cte3*r**2*ro
	
	   do i=1,nbelem	!d(1,1)=0
	    bess(i)=cte11*(y(i)-comp_t(i)-dt*dcomp(i))
	    beds(i)=cte31*v(i)*y(nbelem*(1-1)+i)	!coeff. de Xi: vitesse
     1	+cte21*(d(1,i)*y(nbelem*(2-1)+1)+d(2,i)*y(nbelem*(2-1)+i))
	   enddo
 
	   if(k .eq. 0)then
	    if(kd .eq. 1)then
	     do i=1,nbelem	
	      bes(i)=bess(i)
	      bed(i)=beds(i)
	     enddo
	     print*,'nouvelles valeurs (k=0)'
	     write(6,2000)(bes(i),i=1,nbelem)
	     write(6,2000)(bed(i),i=1,nbelem)
	     print*
	    endif
	
	   else		!derivees numerique / Xk
	    print*,'derivee aes, k, kd',k,kd
	    write(6,2000)(aes(nbelem*(nbelem*(kd-1)+k-1)+i),
     1	(bess(i)-bes(i))/dstor,i=1,nbelem)
	    print*,'derivee aed'
	    write(6,2000)(aed(nbelem*(nbelem*(kd-1)+k-1)+i),
     1    (beds(i)-bed(i))/dstor,i=1,nbelem)
	    print*	
	    y(nbelem*(kd-1)+k)=stor
	   endif
	  enddo	!kd
	 enddo	!k
	 pause'eq_diffus'
	endif		!test
 
	return
 
c	perte de masse on ajoute X_i Mdot dt aux derniers produits scalaires
 
200	do i=1,nbelem
	 bes(i)=y(i)*cte4
	 aes(nbelem*(nbelem*(1-1)+i-1)+i)=cte4	!derivee/Xi
	enddo
c	write(6,2000)mstar,cte4,bes(1),y(1),mdot
 
	return
 
	end
 
 
