 
c********************************************************************
 
	subroutine slinf(x,y,n,l)
 
c	y(n) : suite croissante (non strictement) des points de table
c	recherche de l'indice l tel que :
c		       y(l) .le. x .lt. y(l+1)
 
c	s'il y a debordement a gauche ou a droite on prend:
c	x < y(1) <= y(l) < y(l+1) ou y(l) < y(l+1) <= y(n) < x
 
c	on commence la recherche au voisinage de l
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91
 
	implicit none
 
	integer n,l
 
	real*8 y(1),x
 
 
	if(x .le. y(1))then
	 l=1
	 if(x .lt. y(1))write(6,*)'slinf x .lt. y(1)',x,' < ',y(1)
	 do while(y(l+1) .le. y(1))
	  l=l+1
	 enddo
	elseif(x .ge. y(n))then
	 l=n-1
	 if(x .gt. y(n))write(6,*)'slinf x .gt. y(n)',x,' > ',y(n)
	 do while(y(l) .ge. y(n))
	  l=l-1
	 enddo
	else
	 l=max(1,min(l,n-1))		!initialisation de l
	 if(x .ge. y(l))then
	  do while(x .ge. y(l+1))
	   l=l+1
	  enddo
	 else
	  do while(x .lt. y(l))
	   l=l-1
	  enddo
	 endif
	endif
 
	return
 
	end
