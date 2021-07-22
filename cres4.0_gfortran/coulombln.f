 
c*********************************************************************
	
	subroutine coulombln(zi,zj,thetae,ro,x,t,lnlambij,cij)
 
c	calcul du logarithme de Coulomb
c	d'apres Michaux-Proffit, IAU 137, p.250
 
c	Auteur: J. Mathias, DASGAL, pour H et He
c	adaptation: P. Morel, Departement J.D. Cassini, O.C.A.,
c	Observatoire de Nice
 
c	CESAM, Version 3
 
	real*8	zi,zj,thetae,ro,x,t,lnlambij,cij
 
 
	lnlambij = - 19.26 - log(zi*zj) -0.5*log(ro)
     1		   - 0.5*log(1.+((x+1.)/2.)*thetae) + 1.5*log(t)
 
	cij = log(exp(1.2*lnlambij)+1.)/1.2
 
c	print*,cij,lnlambij
c	write(6,2000)cij,lnlambij,zi,zj,thetae,ro,t,x
2000	format(1x,1p8d10.3)
 
	return
	
	end
 
