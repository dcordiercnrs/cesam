 
c**********************************************************
 
	subroutine sn1dn_dis(nf,f,x,xt,n,knot,dfdx1,dfdxn,
     1	nd,id,x_g,df_g,x_d,df_d)
 
c	calcul des coefficients de la
c	spline naturelle (ordre m=4) d'interpolation dans le tableau f(nf,n)
c	avec comme conditions limites f'(nf,x0) et f'(nf,xn) connus
c	et avec discontinuites des derivees premieres en certains points
 
c	s'exploite avec sbsp1dn
 
c	knot=n+6+2*nd, formation du tableau xt(knot) des points de table
c	!!le tableau des donnees f est donc modifie!!
 
c	!!!!!!!!ATTTENTION: il faut dimensionner f(nf,n+4+2*nd)!!!!!!!!!!
c	et xt a n+6+2*nd
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 22 07 95
 
c entrees
c	nf: nombre de fonctions
c	x: abscisses
c	n: nombre  de points
c	dfdx1(nf): derivees 1-ieres en x0
c	dfdxn(nf): derivees 1-ieres en xn
c	nd: nombre de discontinuites de la derivee premiere
c	id: indices des discontinuites
c	df_g: derivee premiere a gauche en x_g
c	df_d: derivee premiere a droite en x_d
 
c entrees/sorties
c	f(nf,n+2*(nd+1)): valeurs a droite a interpoler/coefficients des splines
 
c sorties
c	knot: nb de noeuds
c	xt: points de raccord
 
	implicit none
 
	integer pn,pnd
	parameter (pn=850, pnd=15)
 
	integer n,l,knot,i,j,nf,indpc(pn+2*(pnd+1)),ligne,mi(pn),k,
     1	n_ligne,nd,id(1)
 
	real*8 f(1),xt(1),x(1),dfdx1(1),d(16),x_d(1),x_g(1),
     1	a((pn+2*(pnd+1))*4),dfdxn(1),df_g(1),df_d(1)
 
 
2000	format((1x,1p8d10.3))
 
c	write(6,*)'entree dans s...'
 
	if(n .gt. pn .or. nd .gt. pnd)then	!test de dimensions
	 write(6,*)'dans sn1dn_dis ajustement de parametres'
	 write(6,*)'mettre le parametre pn a ',n	
	 write(6,*)'mettre le parametre pnd a ',nd
	 stop
	endif
 
c	calcul des B-splines
 
	mi(1)=4 	!vecteur nodal: continuite jusqu'a l'ordre 2
	do i=2,n
	 mi(i)=1
	enddo
	mi(n)=4
	do i=1,nd
	 mi(id(i))=3		!discontinuite de la derivee premiere
	enddo
	call noeud(x,xt,mi,n,knot)
c	pause'noeud'
	ligne=knot-4
c	print*,knot,ligne
c	write(6,2000)(xt(i),i=1,knot)
c	pause'1'
 
	l=4
	do i=1,n
	 call slinf(x(i),xt,knot,l)
	 call bval(x(i),xt,4,l,d)
	 do j=1,4
	  a(ligne*(j-1)+i)=d(j)
	 enddo	!j
	 indpc(i)=l-4+1
c	 write(6,*)indpc(i),i
c	 write(6,2000)(a(ligne*(j-1)+i),j=1,4)
	enddo	!i
c	pause'2'
 
c	les derivees aux extremites
 
	l=4
	call bvald(x(1),xt,4,l,2,d)		!a gauche
	do j=1,4
	 a(ligne*(j-1)+n+1)=d(4*(j-1)+2)		!derivee 1-iere
	enddo
	indpc(n+1)=l-4+1	
c	write(6,*)'premier point',indpc(n+1),l
c	write(6,2000)(a(ligne*(j-1)+n+1),j=1,4)
	do i=1,nf
	 f(nf*n+i)=dfdx1(i)
	enddo
c	pause'3'
	
	l=knot-4
	call bvald(x(n),xt,4,l,2,d)		!a droite
	do j=1,4
	 a(ligne*(j-1)+n+2)=d(4*(j-1)+2)		!derivee 1-iere
	enddo
	indpc(n+2)=l-4+1
c	write(6,*)'dernier point',indpc(n+2),l
c	write(6,2000)(a(ligne*(j-1)+n+2),j=1,4)	
	do i=1,nf
	 f(nf*(n+1)+i)=dfdxn(i)
	enddo
c	pause'derivees'
 
c	derivee a gauche et a droite discontinues
 
	n_ligne=n+2		!indice de ligne
	do k=1,nd
	 call slinf(x_d(k),xt,knot,l)	!a droite
	 call bvald(x_d(k),xt,4,l,2,d)
	 n_ligne=n_ligne+1
	 do j=1,4
	  a(ligne*(j-1)+n_ligne)=d(4*(j-1)+2)		!derivee 1-iere
	 enddo
	 indpc(n_ligne)=l-4+1	
	 do i=1,nf
	  f(nf*(n_ligne-1)+i)=df_d(nf*(k-1)+i)
	 enddo	
c	 write(6,*)'a droite',indpc(n_ligne),n_ligne,l,id(k)
c	 write(6,2000)(a(ligne*(j-1)+n_ligne),j=1,4),x_d(k)	
	
	 n_ligne=n_ligne+1
	 call slinf(x_g(k),xt,knot,l)	!indice a gauche
	 call bvald(x_g(k),xt,4,l,2,d)	!limite a gauche
	 do j=1,4
	  a(ligne*(j-1)+n_ligne)=d(4*(j-1)+2)		!derivee 1-iere
	 enddo
	 indpc(n_ligne)=l-4+1	
c	 write(6,*)'a gauche',indpc(n_ligne),n_ligne,l,id(k)
c	 write(6,2000)(a(ligne*(j-1)+n_ligne),j=1,4),x_g(k)	
	 do i=1,nf
	  f(nf*(n_ligne-1)+i)=df_g(nf*(k-1)+i)
	 enddo
	enddo
c	pause'5'	
 
c	write(6,*)'ligne,nf',ligne,nf
c	write(6,*)(indpc(i),i=1,ligne)
c	do i=1,ligne
c	 write(6,2000)(a(ligne*(j-1)+i),j=1,4),(f(nf*(i-1)+j),j=1,nf)
c	enddo	
 
	call gausdn(a,f,indpc,ligne,4,nf)
 
	return
 
	end
