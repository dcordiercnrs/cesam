c**********************************************************************
 
	function long(a)
 
c	Determination de la longueur effective d'une chaine
 
c	Auteur: G.Gonczi, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91
 
	implicit none
 
	integer*4 long
 
	character*(*)a
 
	long=len(a)
 
	do while(a(long:long) .eq. ' ' .and. long .gt. 1)
	 long=long-1
	enddo
 
	return
 
	end
