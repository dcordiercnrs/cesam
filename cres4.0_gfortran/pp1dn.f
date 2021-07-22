 
c***************************************************************
 
	subroutine pp1dn(n,x,xt,xr,nx,nrx,mx,xx,knotx,lx,dfdx,f,s,fx,init)
 
c	interpolation PP1D pour un ensemble de n tables
 
c	gain de temps par rapport a n fois PP1D car :
c	1-- calcul simultane des coefficients
c	2-- une seule recherche des indices
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91
 
c entree
c	n : nombre de tables
c	x : abscisses
c	mx : ordre des B-splines
c	xx : abscisse d'interpolation
c	init=.false. : il faut calculer les coefficients d'interpolation
 
c entrees / sorties
c	lx : encadrement d'indice approximatif de xx
c	f(n,nx) : tables a interpoler/coeff des splines
c	xt : points de table			(entree si init=.true.)
c	xr : limites des PP			(entree si init=.true.)
c	nrx : nombre de limites xr		(entree si init=.true.)
c	knotx : nombre de points de table	(entree si init=.true.)
c	s(mx,nrx-1,n) : table des coefficients	(entree si init=.true.)
 
c sortie
c	dfdx(n) : derivees premieres
c	fx(n) : interpoles
 
	implicit none
 
	integer pn,pm
	parameter (pn=800, pm=4)
 
	integer nx,mx,lx,knotx,indpc(pn),ix,sx,i,j,nrx,n,k
 
	real*8 x(1),xt(1),xx,dfdx(1),f(1),s(1),ax(pn*pm),qx(pm*pm),
     1	deriv,fac(pm),xr(1),fx(1),xx1
 
	logical init
 
	if(.not. init)then		!les coefficient sont a calculer
 
	 if(nx .gt. pn .or. mx .gt. pm)then	!test de dimensions
	  write(6,*)'dans pp1dn ajustement de parametres'
	  write(6,*)'mettre le parametre pm a ',mx
	  write(6,*)'mettre le parametre pn a ',nx
	  stop
	 endif
 
c	 calcul des B-splines
 
	 call snoein(x,xt,nx,mx,knotx)
	 lx=mx
	 do i=1,nx
	  call slinf(x(i),xt,knotx,lx)
	  call bval(x(i),xt,mx,lx,qx)
	  do j=1,mx
	   ax(nx*(j-1)+i)=qx(j)
	  enddo	!j
	  indpc(i)=lx-mx+1
c	  write(6,*)indpc(i)
c	  write(6,2000)(ax(nx*(j-1)+i),j=1,mx)
2000	  format(1x,1p8d10.3)
	 enddo	!i
	 call gausdn(ax,f,indpc,nx,mx,n)
 
c	 initialisation des PP
 
	 fac(1)=1
	 fac(2)=1
	 do i=3,mx
	  fac(i)=fac(i-1)*float(i-1)
	 enddo	!i
 
	 nrx=nx-mx+2
	 do ix=1,nrx	!les xr sont identiques aux xt sans multiplicite
	  xr(ix)=xt(ix+mx-1)
	 enddo
c	 write(6,*)'les xr'
c	 write(6,2000)(xr(i),i=1,nrx)
c	 write(6,*)'les x'
c	 write(6,2000)(x(i),i=1,nx)
c	 write(6,*)' s '
 
	 do k=1,n	!pour chaque table
	  lx=mx
	  do i=1,nrx-1
	   call slinf(xr(i),xt,knotx,lx)
	   call bvald(xr(i),xt,mx,lx,mx,qx)
	   do ix=1,mx		!ix-1=ordre de la derivee
	    deriv=0.
	    do sx=1,mx		!numero de la spline
	     deriv=deriv+f(n*(lx-mx+sx-1)+k)*qx(mx*(sx-1)+ix)
	    enddo	!sx
 
c	s(mx,nrx-1,n) --> s(mx-ix+1,i,k)=s(mx*((nrx-1)*(k-1)+i-1)+mx-ix+1)
c	pour la fonction k, coefficient de degre ix pour la i-ieme portion
 
	    s(mx*((nrx-1)*(k-1)+i-1)+mx-ix+1)=deriv/fac(ix)
	   enddo	!ix
	  enddo	!i
	 enddo	!k	pour toutes le n variables
 
	endif	!pour l'initialisation
 
c	localisation et interpolation
 
	if(xx .lt. xt(1))then
	 write(6,*)'extrapolation dans PP1DN xx < xt(1),n=',n
	 write(6,*)xx,xt(1)
c	 write(6,*)'n,nrx',n,nrx
c	 write(6,2000)(xr(i),i=1,nrx)
	 lx=1
	 xx1=xr(1)
	elseif(xx .gt. xr(nrx))then
	 write(6,*)'extrapolation dans PP1DN xr(n) < xx,n=',n
	 write(6,*)xr(nrx),xx
	 lx=nrx-1
	 xx1=xr(nrx)
	else
	 call slinf(xx,xr,nrx,lx)
	 xx1=xx
	endif
 
c	interpolation
 
c	write(6,*)'xr,xx'
c	write(6,2000)xr(lx),xx
c	write(6,*)'lx,mx*(lx-1)+1',lx,mx*(lx-1)+1
 
	do k=1,n	!pour chaque fonction
	 call horder(mx,mx*((k-1)*(nrx-1)+(lx-1))+1,s,xx1-xr(lx),
     1	fx(k),dfdx(k))
	enddo	!k
 
	return
 
	end
