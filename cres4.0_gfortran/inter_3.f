c********************************************************************
 
	subroutine inter_3(mu_ou_r2,bp,q,qt,n,knot,x,f,dfdx,r2,mu)
	
c	en x=m**23 ou x=r2
c	on fait une interpolation inverse mu-->q-->f ou r2-->q-->f
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version 3
 
c	ATTENTION : il faut initialiser et utiliser les bons bp, mu et r2 !!!
 
c entree
c	mu_ou_r2 = 'mu' : en m23, = 'r2' : en r2
c	ne, bp,  : nb de fonctions (6)
c	bp,q,qt,n,m,knot : spline a interpoler
c	x : point d'interpolation
 
c entree/sortie
c	r2, mu : abscisses
 
c sorties
c	f, dfdx : valeurs et derivees/r2 ou mu
c	ATTENTION f(6)=q_int : point d'interpoletion
 
	implicit none
	
	include 'cesam_3.parametres'
	include 'modele_3.common'
 
	integer knot,i,ntour,l,n,inc,l0
 
	real*8 bp(1),q(1),qt(1),x,mu(1),r2(1),q_int,q0,corr,f(1),dfdx(1)
 
	character*2 mu_ou_r2
 
2000	format((1x,1p8d10.3))
 
c	cas ou on donne le rayon ou la masse
 
	if(mu_ou_r2 .eq. 'mu')then
	 inc=5	!indice de l'inconnue pour l'interpolation inverse
	 if(x .le. mu(1))then
	  q_int=1.d0
	 elseif(x .ge. mu(n))then
	  q_int=dfloat(n)
	 else
	  call slinf(x,mu,n,l)		!entre l et l+1
	  q_int=(x-mu(l))/(mu(l+1)-mu(l))+q(l)	!approx. lin.
	 endif
	 q_int=max(1.d0,min(q_int,dfloat(n)))
	elseif(mu_ou_r2 .eq. 'r2')then
	 inc=3
	 if(x .le. r2(1))then
	  q_int=1.d0
	 elseif(x .ge. r2(n))then
	  q_int=dfloat(n)
	 else
	  call slinf(x,r2,n,l)		!entre l et l+1
	  q_int=(x-r2(l))/(r2(l+1)-r2(l))+q(l)	!approx. lin.
	 endif
	 q_int=max(1.d0,min(q_int,dfloat(n)))	
	else
	 print*,'erreur d''identification dans inter_3, mu_ou_r2=',mu_ou_r2
	 stop
	endif
	
c	determination de l'indice pour interpolation inverse
 
c	print*,mu_ou_r2,inc
c	write(6,2000)q_int,x	
	
	if(q_int .gt. 1.d0 .and. q_int .lt. dfloat(n))then
	 q0=q_int
	 l0=l
	 corr=1.d30
	 ntour=0
	 do while(ntour .lt. 30 .and. abs(corr) .gt. 1.d-6)
	  call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,q_int,l,f,dfdx)
	  corr=f(inc)-x
	  q_int=q_int-corr/dfdx(inc)	!NR: Xn+1=Xn-f(Xn)/f'(Xn)
	  q_int=max(1.d0,min(q_int,dfloat(n)))
	  ntour=ntour+1
	 enddo
c	 write(6,*)ntour,corr
c	 write(6,100)ntour,inc,x,f(inc),f(inc)-x,dfdx(inc),corr,q_int,q0
100	 format(1x,2i4,1p8d10.3)
	 if(ntour .ge. 30)then
	  print*,'pas de conv. newton dans inter_3, approx. lin. ',mu_ou_r2,inc
	  print*,x,f(inc),q_int,q0
c	  if(mu_ou_r2 .eq. 'mu')then
c	   write(6,2000)mu(l0),x,mu(l0+1),q_int,q0,corr,dfdx(inc)
c	  else
c	   write(6,2000)r2(l0),x,r2(l0+1),q_int,q0,corr,dfdx(inc)
c	  endif	
c	  write(6,2000)(mu(l),l=1,n)
c	  write(6,2000)(r2(l),l=1,n)
c	  write(6,2000)(q(l),l=1,n)
c	  call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,q0,l,f,dfdx)
c	  print*,l,knot,mpr,ne,n
c	  write(6,2000)f(inc),q0,x	
	
c	  q_int=q0	!on ecrit la suite des iterations
c	  corr=1.d30
c	  ntour=0
c	  do while(ntour .lt. 30 .and. abs(corr) .gt. 1.d-6)
c	   call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,q_int,l,f,dfdx)
c	   corr=f(inc)-x
c	   q_int=q_int-corr/dfdx(inc)		!NR: Xn+1=Xn-f(Xn)/f'(Xn)
c	   write(6,2000)q_int,corr,f(inc),x,dfdx(inc)
c	   q_int=max(1.d0,min(q_int,dfloat(n)))
c	   ntour=ntour+1
c	  enddo
c	  pause'inter_3'
	  q_int=q0
	 endif
	endif
 
c	on retrouve la solution en q_int et derivees
 
c	print*,ntour
c	write(6,2000)corr,q_int
 
	call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,q_int,l,f,dfdx)
	
	do i=1,ne
	 if(i .ne. inc)dfdx(i)=dfdx(i)/dfdx(inc)
	enddo
	f(6)=q_int
c	write(6,2000)(f(i),i=1,ne)
c	write(6,2000)(dfdx(i),i=1,ne)	
 
	return
 
	end
