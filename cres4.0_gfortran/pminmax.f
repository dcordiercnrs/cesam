	
c***********************************************************************
 
	subroutine pminmax(x,n,xmax,xmin)
 
c	calcul de min/max pour le plot
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91
 
	implicit none
 
	integer i,n
 
	real*4 x(1),xmax,xmin,dx
 
	xmax=-1.d37
	xmin=1.d37
	do i=1,n
	 xmax=max(xmax,x(i))
	 xmin=min(xmin,x(i))
	enddo
	if(xmax .eq. xmin)then
	 if(x(1) .eq. 0.)then
	  dx=.051
	 else
	  dx=.05*x(1)
	 endif
	elseif(xmax+xmin .eq. 0.)then
	 dx=(xmax-xmin)/20.
	elseif((xmax-xmin)/(xmax+xmin) .lt. 1.d-2)then
	 dx=.05*(xmax+xmin)			
	else
	 dx=(xmax-xmin)/20.
	endif
	xmax=xmax+abs(dx)
	xmin=xmin-abs(dx)
 
	return
 
	end
