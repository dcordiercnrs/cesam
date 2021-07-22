 
c**************************************************************
 
	subroutine gausdn(a,b,indpc,nl,nemr,n)
 
c	resolution par la methode du pivot d'un systeme lineaire bande
c	avec n second membres
c	les indpc etant reordonnes
c	on ne garde que les coefficients non identiquement nuls
 
c	methode: pivot partiel avec equilibrage
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91
 
c entree
c	a(nl,nemr) : matrice des coefficients non identiquement nuls
c	indpc(nl) : indice de la premiere colonne de chaque ligne
c	nl : nombre de lignes
c	nemr : longueur fixe d'une ligne
c	n : nombre de colonnes au second membre
 
c entree/sortie
c	b(n,nl) : second membre/solution
 
	implicit none
 
	integer indpc(1),nl,nemr,ligne,ipivot,i,j,k,n,isto
 
	real*8 a(1),b(1),pivot,sto,ai
 
c	equilibrage de la matrice
 
	do ligne=1,nl
	 sto=0.
	 do i=1,nemr
	  sto=max(sto,abs(a(nl*(i-1)+ligne)))
	 enddo	!i
	 if(sto .eq. 0.)then
	  write(6,*)'ligne identiquement nulle dans GAUSDN ligne=',ligne
	  stop
	 endif
	 do i=1,nemr
	  a(nl*(i-1)+ligne)=a(nl*(i-1)+ligne)/sto
	 enddo	!i
	 do k=1,n
	  b(n*(ligne-1)+k)=b(n*(ligne-1)+k)/sto
	 enddo		!k
	enddo		!ligne
 
c	mise en ordre des lignes pour que indpc soit croissant
 
	do ligne=1,nl-1
	 if(indpc(ligne) .gt. indpc(ligne+1))then
	  i=ligne
	  do while(i .ge. 1)
	   if(indpc(i) .gt. indpc(i+1))then
	    do k=1,n
	     sto=b(n*(i-1)+k)
	     b(n*(i-1)+k)=b(n*i+k)
	     b(n*i+k)=sto
	    enddo
	    isto=indpc(i)
	    indpc(i)=indpc(i+1)
	    indpc(i+1)=isto
	    do j=1,nemr
	     sto=a(nl*(j-1)+i)
	     a(nl*(j-1)+i)=a(nl*(j-1)+i+1)
	     a(nl*(j-1)+i+1)=sto
	    enddo	!j
	   endif
	   i=i-1
	  enddo		!while
	 endif
	enddo		!ligne
 
c	elimination de gauss avec pivots
 
	do ligne=1,nl
c	 write(6,*)' '
c	 write(6,*)ligne
c	 write(6,*)' '
	 pivot=0.
	 do i=ligne,nl
	  if(indpc(i) .gt. ligne)goto 10
c	  write(6,*)i,a(i)
	  if(pivot .lt. abs(a(i)))then
	   pivot=abs(a(i))
	   ipivot=i
c	   write(6,*)'pivot',pivot,ipivot
	  endif
	 enddo	!i
	 i=nl		!sortie de boucle avec i+1
10	 if(pivot .eq. 0.d0)then
	  write(6,*)'pivot nul dans GAUSDN ligne : ',i
	  stop
	 endif
c	 write(6,*)'pivot pour ligne',ligne,pivot,ipivot
 
c	 permutation de ligne et de ipivot, division par pivot
 
	 do k=1,n
	  sto=b(n*(ligne-1)+k)
	  b(n*(ligne-1)+k)=b(n*(ipivot-1)+k)/a(ipivot)
	  if(ligne .ne. ipivot)b(n*(ipivot-1)+k)=sto
	 enddo
	 do j=nemr,1,-1
	  sto=a(nl*(j-1)+ligne)
	  a(nl*(j-1)+ligne)=a(nl*(j-1)+ipivot)/a(ipivot)
	  if(ligne .ne. ipivot)a(nl*(j-1)+ipivot)=sto
	 enddo	!j
 
c	 triangulation et decalage pour avoir la diagonale en colonne 1
 
	 do i=ligne+1,nl
	  if(indpc(i) .gt. ligne)goto 20
c	  write(6,*)'i de perm',i
	  ai=a(i)
	  do k=1,n
	   b(n*(i-1)+k)=b(n*(i-1)+k)-b(n*(ligne-1)+k)*ai
	  enddo		!k
	  do j=2,nemr
	   a(nl*(j-2)+i)=a(nl*(j-1)+i)-a(nl*(j-1)+ligne)*ai
	  enddo	!j
	  a(nl*(nemr-1)+i)=0.
	 enddo	!i
 
c	do i=1,nl
c	 write(6,1000)(a(nl*(j-1)+i),j=1,nemr),(b(n*(i-1)+k),k=1,n)
c1000	 format((1x,1p16e8.1))
c	enddo	!i
 
20	enddo	!ligne
 
c	on remonte
 
	do i=nl-1,1,-1
	 do j=2,nemr
	  if(i+j-1 .le. nl)then
	   do k=1,n
	    b(n*(i-1)+k)=b(n*(i-1)+k)-a(nl*(j-1)+i)*b(n*(i+j-2)+k)
	   enddo		!k
	  endif
	 enddo	!j
	enddo	!i
 
	return
 
	end
