 
c**********************************************************************
 
	subroutine colloc(mx,rx,px,x,xt,knox,xx,dx,derxx,nlx,lderxx,xl)
 
c	initialisation des valeurs des B-splines et de leurs derivees aux
c	points de collocation en x et aux limites pour vecteur nodal de de Boor
 
c entrees:
c	mx : ordre des splines
c	rx : ordre des equ.diff
c	px : nb. de points
c	x : abscisses
c	xt : points de table
c	knox : nb. de points de tables
c	dx(mx+rx**2) : vecteur de travail
c	nlx : nombre de limites
c	xl : limites
 
c sorties
c	xx : points de collocation
c	derxx : derivees
c	lderxx : derivees aux limites
 
c	initialisation des valeurs des B-splines et de leurs derivees aux
c	points de collocation en x, xx(px-1,mx)
c	indices : a(i,j,k,l)=a(ni{nj[nk(l-1)+k-1]+j-1}+i)
c	derxx(der,B-spl.,pt.coll.,pt.rac.) : derxx(rx+1,mx+rx,mx,px-1)
 
c	initialisation des valeurs des B-splines et de leurs derivees aux
c	points limite en x
c	lderxx(rx+1,mx+rx,rx*eq) : lderxx(der,B-spl.,pt.lim)
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91
 
	implicit none
 
	integer mx,rx,lx,ipx,knox,ix,cx,derx,indice,nlx,px,splx
 
	real*8 x(1),xt(1),xx(1),colpnt,rho(10),dx(1),derxx(1),lderxx(1),xl(1)
 
	logical init
 
	lx=mx+rx
	init=.false.	!pour appel a colpnt
	do ipx=1,px-1
	 call slinf(x(ipx),xt,knox,lx)
	 do ix=1,mx
	  cx=(px-1)*(ix-1)+ipx	!indice du point de collocation
	  xx(cx)=colpnt(x,ix,ipx,mx,rho,init)
	  call bvald(xx(cx),xt,mx+rx,lx,rx+1,dx)
	  do splx=1,mx+rx
	   do derx=1,rx+1
	    indice=(rx+1)*((mx+rx)*(mx*(ipx-1)+ix-1)+splx-1)+derx
	    derxx(indice)=dx((mx+rx)*(splx-1)+derx)
	   enddo	!derx
	  enddo		!splx
	 enddo		!ix
	enddo		!ipx
 
	do ix=1,nlx
	 call slinf(xl(ix),xt,knox,lx)
	 call bvald(xl(ix),xt,mx+rx,lx,rx+1,dx)
	 do splx=1,mx+rx
	  do derx=1,rx+1
	   indice=(rx+1)*((mx+rx)*(ix-1)+splx-1)+derx
	   lderxx(indice)=dx((mx+rx)*(splx-1)+derx)
	  enddo		!derx
	 enddo		!splx
	enddo		!ix
 
c	do ipx=1,px-1
c	 do ix=1,mx
c	  cx=(px-1)*(ix-1)+ipx
c	write(6,*)'cx,xx(cx)',cx,xx(cx)
c	write(6,1000)((derxx((rx+1)*((mx+rx)*(mx*(ipx-1)+ix-1)+splx-1)+derx),
c	1	splx=1,mx+rx),derx=1,rx+1)
c1000	format((1x,1p10e10.3))
c	 enddo	!ix
c	enddo	!ipx
 
c	do ix=1,nlx
c	 write(6,*)'ix,xl(ix)',ix,xl(ix)
c	 write(6,1000)((lderxx((rx+1)*((mx+rx)*(ix-1)+splx-1)+derx),
c	1	splx=1,mx+rx),derx=1,rx+1)
c	enddo 	!ix
 
	return
 
	end
