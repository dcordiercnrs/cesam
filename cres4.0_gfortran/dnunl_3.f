	
c**************************************************************************
 
	subroutine dnunl_3(r,c,n,rtot,nu0,dnu02,dnu13,a)
	
c	calcul de nu0 et delta nu  n l, formule 100, 101, 102 p.389
c	du Schatzman et Praderie,
c	Formule 102: delta nu 02 ~ 6A/20, delta nu 13 ~ 10A/20, J. Provost
 
c	doit etre utilise avec CESAM_3
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	CESAM version 3
 
c entrees:
c	r(n): rayons du centre a l'exterieur (Rsol)
c	c(n): vitese du son (cgs)
c	n: nombre de points
c	rtot: rayon (peut differer de r(n)
 
c sorties
c	nu0: fondamental
c	dnu02, dnu13: ecarts approximatifs de frequence delta nu 02 et 13
c	a: formule 100, ~ - somme 1/r dc/dr = 2 somme dc /dr**2
 
	implicit none
	
	include 'cesam_3.parametres'	
	include 'ctephy.common'
	
	integer n,ig,i,m,knotr,l,mg
	
	real*8 r(1),c(1),rtot,nu0,dnu02,dnu13,cr(1),dcdr(1),
     1	rg(6),wg(6),r2t(pn+4),a,r2(pn)
	
 
2000	format(1x,1p8d10.3)
 
c	print*,n
c	write(6,2000)rtot
c	write(6,2000)(c(i),i=1,n)
c	write(6,2000)(r(i),i=1,n)
		
c	tabulation de c(r**2)
 
	do i=1,n
	 r2(i)=r(i)**2
	enddo
 
	m=4	!ordre de l'interpolation et m/2 pour l'integrale de Gauss
	call sbsp1dn(1,c,r2,r2t,n,m,knotr,.false.,r2(1),l,cr,dcdr)
 
c	recherche de la limite et calcul des integrales
c	calcul approximatif: on ne tient pas compte de l'atmosphere
 
	a=0
	nu0=0
	mg=2
	do i=1,n-1
	 call intgauss(r(i),r(i+1),rg,wg,mg) !int.Gauss
	  do ig=1,mg
	   call sbsp1dn(1,c,r2,r2t,n,m,knotr,.true.,rg(ig)**2,l,cr,dcdr)
	   nu0=nu0+wg(ig)/cr(1)
	   a=a+wg(ig)*dcdr(1)
	  enddo
	enddo
	a=2.*a
c	pause
	
	nu0=2.*nu0
	nu0=1.d6/nu0
	
	call sbsp1dn(1,c,r2,r2t,n,m,knotr,.true.,r(n)**2,l,cr,dcdr)
c	write(6,2000)a,cr(1),r(n)
	a=(cr(1)/r(n)-a)/(2.*pi)**2*1.d6	!formule 100, p.389, modifiee
c	write(6,2000)a
c	pause
	dnu02=6.*a/20
	dnu13=10.*a/20
	
c	write(6,2000)cr(1),a,nu0,dnu02
 
	return
	
	end
	
	
	
	
	
		
