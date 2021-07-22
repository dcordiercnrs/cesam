 
c**************************************************************
 
	subroutine gausdp_g(a,b,indpc,nl,n,nemr,inversible)
 
c	resolution par la methode du pivot d'un systeme lineaire bande de
c	en ne gardant que les coefficients non identiquement nuls,
c	avec un nombre d'equations superieur au rang comme dans la methode de
c	Galerkin (convient egalement et sans surcout, si nl=n)
 
c	methode: pivot partiel avec equilibrage
 
c	Les indpc ne sont pas necessairement croissants.
 
c	A l'issue du pivotage, le rang etant n, les nl-n dernieres equations
c	sont 0=0 aux arrondis pres
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
 
c entree
c	n: nombre de colonnes = rang
c	nl : nombre de lignes .ge. au rang = n
c	nemr : longueur max. d'une ligne
 
c entree/sortie
c	a(nl,nemr) : matrice des coefficients non identiquement nuls
c	indpc(nl) : indice de la premiere colonne de chaque ligne
c	b(nl) : second membre de 1 a nl/solution de 1 a n
 
c sortie
c	inversible=.true. : solution obtenue
 
	implicit none
 
	integer indpc(1),nl,nemr,ligne,ipivot,i,j,isto,n
 
	real*8 a(1),b(1),pivot,sto,ai
 
	logical inversible
 
c	equilibrage de la matrice
 
	do ligne=1,nl
	 sto=0.
	 do i=1,nemr
	  sto=max(sto,abs(a(nl*(i-1)+ligne)))
	 enddo	!i
	 if(sto .eq. 0.)then
	  write(6,*)'ligne identiquement nulle dans gausdp_g, ligne=',ligne
	  sto=1.
	 endif
	 do i=1,nemr
	  a(nl*(i-1)+ligne)=a(nl*(i-1)+ligne)/sto
	 enddo	!i
	 b(ligne)=b(ligne)/sto
	enddo	!ligne
 
c	mise en ordre des lignes pour que indpc soit croissant
 
	do ligne=1,nl-1
	 if(indpc(ligne) .gt. indpc(ligne+1))then
	  i=ligne
	  do while(i .ge. 1)
	   if(indpc(i) .gt. indpc(i+1))then
	    sto=b(i)
	    b(i)=b(i+1)
	    b(i+1)=sto
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
 
	do ligne=1,n			!elimination jusqu'a ligne n=rang
c	 write(6,*)' '
c	 write(6,*)ligne
c	 write(6,*)' '
	 pivot=0.
	 do i=ligne,nl			!recherche du pivot jusqu'a nl
	  if(indpc(i) .gt. ligne)goto 10
c	  write(6,*)i,a(i)
	  if(pivot .lt. abs(a(i)))then
	   pivot=abs(a(i))
	   ipivot=i
c	   write(6,*)'pivot',pivot,ipivot
	  endif
	 enddo	!i
	 i=nl		!sortie de boucle avec i+1
10	 inversible=pivot .gt. 0.
	 if(.not. inversible)then
	  write(6,*)'pivot nul dans gausdp_g, ligne',ligne,' pivot',pivot
	  return
	 endif
c	 write(6,*)'pivot pour ligne',ligne,pivot,ipivot
 
c	 permutation de ligne et de ipivot, division par pivot
 
	 sto=b(ligne)
	 b(ligne)=b(ipivot)/a(ipivot)
	 if(ligne .ne. ipivot)b(ipivot)=sto
	 do j=nemr,1,-1
	  sto=a(nl*(j-1)+ligne)
	  a(nl*(j-1)+ligne)=a(nl*(j-1)+ipivot)/a(ipivot)
	  if(ligne .ne. ipivot)a(nl*(j-1)+ipivot)=sto
	 enddo	!j
 
c	 decalage pour avoir la diagonale en colonne 1
 
	 do i=ligne+1,nl
	  if(indpc(i) .gt. ligne)goto 20
c	  write(6,*)'i de perm',i
	  ai=a(i)
	  b(i)=b(i)-b(ligne)*ai
	  do j=2,nemr
	   a(nl*(j-2)+i)=a(nl*(j-1)+i)-a(nl*(j-1)+ligne)*ai
	  enddo	!j
	  a(nl*(nemr-1)+i)=0.
	 enddo	!i
 
c	do i=1,nl
c	 write(6,1000)(a(nl*(j-1)+i),j=1,nemr),b(i)
c1000	 format((1x,1p16e8.1))
c	enddo	!i
 
20	enddo	!ligne
 
c	on remonte de n (seulement) a 1
 
	do i=n-1,1,-1
	 do j=2,nemr
	  if(i+j-1 .le. n)b(i)=b(i)-a(nl*(j-1)+i)*b(i+j-1)
	 enddo	!j
	enddo	!i
 
	return
 
	end
