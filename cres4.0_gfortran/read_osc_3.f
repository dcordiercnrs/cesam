c*****************************************************************
 
	subroutine read_osc_3(nom_fich,itot,iconst,ivar,nbelem,nchim,iw,iz,
     1	modele,nom_elem,glob,var)
 
c	lecture du fichier d'oscillation de cesam_3
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	CESAM, version 3
 
c entrees:
c	nom_fich: nom du fichier d'oscillation sans l'extension .osc
 
c sorties:
c	var: variables
c	itot: nombre de points
c	iconst: nombre de constantes global
c	ivar: nombre de variables du tableau var sans les elements chimiques
c	modele: nom du modele
c	nom_elem: nom des elements chimiques utilises
c	nbelem: nombre d'elements chimiques utilises
c	glob: variables globales
c		glob(1)=mstar*msol
c		glob(2)=rtot*rsol
c		glob(3)=ltot*lsol
c		glob(4)=z0
c		glob(5)=x0
c		glob(6)=alpha
c		glob(7)=9./4.
c		glob(8)=1./162.
c		glob(9)=X dans ZC
c		glob(10)=Y dans ZC
c		glob(11)=d2p
c		glob(12)=d2ro
c		glob(13)=age
c		glob(14)=vsal
c		glob(15)=w_rot initial
c	var: variables
c		var(1,i)=r*rsol
c	 	var(2,i)=log(m) -1.d38 au centre
c		var(3,i)=t
c		var(4,i)=p
c		var(5,i)=ro
c		var(6,i)=gradient reel d ln T / d ln P
c		var(7,i)=l
c		var(8,i)=kap
c		var(9,i)=energie thermo+gravifique
c		var(10,i)=grand Gamma1
c		var(11,i)=gradient adiabatique
c		var(12,i)=delta
c		var(13,i)=cp
c		var(14,i)=1 / (mu elec.)
c		var(15,i)=vaissala, 0 au centre
c	 	var(16,i)=vitesse angulaire, radian/sec
c	 	var(17,i)=d ln kappa / d ln T
c	 	var(18,i)=d ln kappa / d ln ro
c	 	var(19,i)=d epsilon(nuc) / d ln T
c	 	var(20,i)=d epsilon(nuc) / d ln ro
c	  	var(20+j,i)=xchim(j)*nucleo(j), j=1,nbelem
 
	implicit none
 
	include 'cesam_3.parametres'
 
	integer iconst,ivar,ivers,i,itot,j,long,nbelem,nchim,iw,iz

	external long

	real*8 glob(15),var(20+pnelem,1)
 
	character*3	nom_elem(pnelem)
	character*80 abid,modele
	character*40 nom_fich
 
	print*,'nom du fichier d''oscillations a lire (sans l''extension .osc)'
	read(5,'(a)')nom_fich
	open(unit=30,form='formatted',status='old',
     1	file=nom_fich(:long(nom_fich))//'.osc')
 
	read(30,'(a)')abid	!m 100X707a173_$.dat du 18/ 8/1993 etc...
c	print*,abid
	read(30,'(a)')abid	!fichier pour le calcul des oscillations etc..
c	print*,abid
	read(30,'(a)')modele	! fichier: diff0.osc
c	print*,modele
	read(30,'(a)')abid	!methode: CESAM 3.0.0.0 colloc. etc...
c	print*,abid
	read(30,90)nbelem,(nom_elem(i),i=1,nbelem)	!  9  H  He3 He4 etc...
90	format(i3,14(1x,a3))
c	print*,nbelem,(nom_elem(i),i=1,nbelem)
	read(30,137)itot,iconst,ivar,ivers,nchim,iw,iz	!349  14  20  3  (9)
137	format(7i10)
	read(30,138)(glob(i),i=1,iconst)
138	format(1p5d19.12)
	do j=1,itot
	 read(30,138)(var(i,j),i=1,ivar+nbelem)
	enddo
	close(unit=30)
 
	return
 
	end
