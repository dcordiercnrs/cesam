 
c*************************************************************************
 
	subroutine smatpr(a,la,calb,b,cb,ab)
 
c	produit de matrices ab=a*b
c	la : nombre de lignes de a
c	calb : nombre de colonnes de a et de lignes de b
c	cb : nombre de colonnes de b
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91
 
	implicit none
 
	integer*4 la,calb,cb
	integer*4 i,j,k
 
	real*8 a(1),b(1),ab(1)
 
	do i=1,la
		do j=1,cb
		ab(la*(j-1)+i)=0.
			do k=1,calb
	ab(la*(j-1)+i)=ab(la*(j-1)+i)+a(la*(k-1)+i)*b(calb*(j-1)+k)
			enddo
		enddo
	enddo
 
	return
 
	end
