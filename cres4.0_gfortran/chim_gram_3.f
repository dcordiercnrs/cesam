 
c****************************************************************
 
	subroutine chim_gram_3(xchim,dxchim,nuc_m)
 
c	transforme les abondances par mole ---> par gramme
	
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
 
c	CESAM, version 3
 
c entrees/sorties
c	xchim,dxchim : comp. chim et derivee par mole ---> par gramme
 
c sortie
c	nuc_m : nucleo(i)/mass
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'evol_chim_3.common'
	include 'ctephy.common'
	
	real*8	mass_z,mass,xchim(1),dxchim(1),nuc_m(1)
c	data mass_z/0.d0/
 
	integer i
	
	logical init
 
	data init/.true./
	data mass_z/0.d0/
 
	save
	
	if(init)then
	 init=.false.
	 mass_z=1.d0
	 do i=1,nchim
	  mass_z=mass_z-ab_ini(i)
	 enddo
	 mass_z=max(mass_z,0.d0)
	 write(6,1)mass_z	
	 write(2,1)mass_z
1	 format(1x,/,1x,'-------------------------------------------',//,
     1	1x,'dans Z, masse elem. chim. invariants :',1pd10.3)
	endif
 
	if(nchim .gt. 1)then	
	 mass=mass_z
	 do i=1,nchim
	  mass=mass+xchim(i)*nucleo(i)
	 enddo
	else
	 mass=1.d0
	endif
	
	do i=1,nchim
	 nuc_m(i)=nucleo(i)/mass
	 xchim(i)=xchim(i)*nuc_m(i)
	 dxchim(i)=dxchim(i)*nuc_m(i)
	enddo
	
	
	return
	
	end	
	
