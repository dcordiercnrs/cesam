 
c************************************************************
 
	subroutine sbval(x,y,m,l,q,d)
 
c	calcul les B-splines normalisees d'ordre m>1
c	et leurs derivees premieres au point x tel que y(l).le.x.lt.y(l+1)
c	si x est l'extremite droite : y(l).lt.x.eq.y(l+1) on obtient
c	la limite a droite
c	les valeurs des m b-splines non nulles : l-m+1, ..., l sont
c	dans q(1), ..., q(m), les derivees dans d(1), ... , d(m)
c	!! ATTENTION : q doit etre dimensionne au moins a m+1 !!
 
c	d'apres l'algorithme 5-5 p.192 de schumaker
c	et addition d'une instruction pour les derivees
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91
 
	implicit none
 
	integer*4 m,l
	integer*4 j,i
	
	real*8 y(1),q(1),d(1),x
	real*8 mm1,denom,a1,a2
 
	do j=1,m-1
	q(j)=0.
	end do
	
	q(m)=1./(y(l+1)-y(l))
	q(m+1)=0.
 
	do j=2,m-1
		do i=m-j+1,m
		denom=y(i+l-m+j)-y(i+l-m)
		a1=(x-y(i+l-m))/denom
		a2=1.d0-a1
		q(i)=a1*q(i)+a2*q(i+1)
		enddo
	enddo
 
	mm1=float(m-1)
	do i=1,m
	d(i)=mm1*(q(i)-q(i+1))
	q(i)=(x-y(i+l-m))*q(i)+(y(i+l)-x)*q(i+1)
	enddo
 
	return
 
	end
