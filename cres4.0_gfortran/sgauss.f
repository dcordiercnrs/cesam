 
c*************************************************************
	
	subroutine sgauss(a,b,n,m)
 
c	resolution par la methode de Gauss avec pivot partiel, du
c	systeme lineaire bande obtenu avec les B-splines, d'apres De Boor
 
c	a : matrice des coefficients
c	b : second membre
c	n : ordre du systeme
c	m : ordre des splines
 
c	le rangement dans a est effectue par ligne de 2m-1 elements,
c	l'element diagonal etant le m-ieme, chaque ligne i correspondant
c	a la i-ieme equation. Exemple m=3 n=6 :
 
c		.    .  a11 a12 a13		a11 a12 a13  .   .   .
c		.   a21 a22 a23 a24		a21 a22 a23 a24  .   .
c		a31 a32 a33 a34 a35		a31 a32 a33 a34 a35  .
c		a42 a43 a44 a45 a46		 .  a42 a43 a44 a45 a46
c		a53 a54 a55 a56  .		 .   .  a53 a54 a55 a56
c		a64 a65 a66  .   .		 .   .   .  a64 a65 a66
 
c	Suivant la repartition des points de table, certains aij des bords
c	sont nuls car il y a, au plus, m elements par colonne et 2m-1 elements
c	par ligne, qui peuvent etre non nuls.
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91
 
	implicit none
 
	integer*4 n,m
	integer*4 i,j,k,ipivot
 
	real*8 a(1),b(1)
	real*8 pivot,sto
 
c	equilibrage
 
	do i=1,n
	 sto=0.
	 do j=1,2*m-1
	  sto=max(sto,abs(a(n*(j-1)+1)))
	 enddo	!j
	 if(sto .eq. 0.)then
	  write(6,*)'dans sgauss la ligne',i,'est nulle'
	  stop
	 else
	  do j=1,2*m-1
	   a(n*(j-1)+i)=a(n*(j-1)+i)/sto
	  enddo	!j
	 endif
	enddo	!i
 
c	decalage a gauche pour les m-1 iere lignes
 
	do i=1,m-1
	 do k=1,i
	  do j=m-i+1,2*m-1
	   a(n*(j-2)+k)=a(n*(j-1)+k)
	   a(n*(j-1)+k)=0.
	  enddo
	 enddo
	enddo
 
c	triangulation	ipivot : indice de la ligne du pivot
 
	do i=1,n
	 pivot=0.
	 do k=i,min(n,i+m-1)
	  if(pivot .lt. abs(a(k)))then
	   pivot=abs(a(k))
	   ipivot=k
	  endif
	 enddo
 
c	 cas du pivot nul
 
	 if(pivot .eq. 0.)then
	  write(6,*)'pivot nul dans sgauss; ligne i, n, m',i,n,m
	  stop
	 endif
 
c	 permutation et division par pivot
 
	 if(ipivot .eq. i)then
	  b(i)=b(i)/a(ipivot)
	  do j=2*m-1,1,-1
	   a(n*(j-1)+i)=a(n*(j-1)+ipivot)/a(ipivot)
	  enddo
	 else
	  sto=b(i)
	  b(i)=b(ipivot)/a(ipivot)
	  b(ipivot)=sto
	  do j=2*m-1,1,-1
	   sto=a(n*(j-1)+i)
	   a(n*(j-1)+i)=a(n*(j-1)+ipivot)/a(ipivot)
	   a(n*(j-1)+ipivot)=sto
	  enddo
	 endif
	
c	 elimination de gauss
 
	 if(i .eq. n)goto 10
	 do k=i+1,min(n,i+m-1)
	  b(k)=b(k)-b(i)*a(k)
	  do j=2*m-1,2,-1
	   a(n*(j-1)+k)=a(n*(j-1)+k)-a(n*(j-1)+i)*a(k)
	  enddo
 
c	  decalage
 
	  do j=2,2*m-1
	   a(n*(j-2)+k)=a(n*(j-1)+k)
	   a(n*(j-1)+k)=0.
	  enddo
	 enddo
10	enddo
 
c	on remonte
 
	do i=n-1,1,-1
	 do j=2,2*m-1
	  if(i+j-1 .le. n)b(i)=b(i)-a(n*(j-1)+i)*b(i+j-1)
	 enddo
	enddo
 
	return
 
	end
