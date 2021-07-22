 
c***************************************************************
 
	subroutine ctes_ceph
 
c	Version destinee a etre utilisee avec 'energie_nuc' et
c	'taux-reac'
 
c	constantes physiques de Christensen-Dalsgaard
 
c	"call ctes" initialise le common /ctephy/
 
c	les constantes fondamentales utilisees sont ecrites dans le listing
c	de sortie au premier appel de cette routine
 
c	ce listing est dans le fichier FORTRAN "2"
c	qui, eventuellement, doit etre ouvert precedemment
 
c	Auteurs:
c               P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c               D. Cordier, ENSCR
c	version: 21 01 97
 
c COMMON /ctephy/
c	initialisation de toutes les variables
 
c APPELS externes:
c	neant
 
	implicit none
 
	include 'ctephy.common'
 
	logical init
	data init/.true./
 
	pi=acos(-1.d0)
	echarg=4.8032068d-10
	kbol=1.380658d-16	!boltzman, R=k/mh
	hpl=6.6260755d-27		!planck
	clight=2.99792458d10
	eve=1.60217733d-12	!electron volt
	kih=13.595d0		!potentiels d'ionisation
	kihe=24.580d0
	kihe1=54.403d0
	amu=1.6605402d-24	!masse atom. unite, Avogadro=1/amu
	me=9.1093897d-28	!masse electron
	g=6.67232d-8		!gravite
	aradia=(2.*pi*kbol/clight/hpl)**3*pi*pi*kbol/15.	!sigma=a c /4
	ah=1.00813
	ah2=2.01471
	ahe3=3.01693
	ahe4=4.00388
	ali7=7.01818
	abe7=7.
	ac12=12.00386
	ac13=13.00756
        an13=13.00991
	an14=14.00754
	an15=15.00489
	ao16=16.00000
	ao17=17.00450
        ao18=18.0049
        af19=19.00454
        ane20=19.99872
        ane22=21.99844
        amg24=23.9924
        amg25=24.9938
        amg26=25.9898
        an=1.00896
 
	msol=1.989d33		!le soleil
	lsol=3.846d33
	rsol=6.9599d10
 
	secon6=365.24219878*1.d6*24.*3600.	!nombre de seconde en 10**6 ans
 
	if(init)then	!ecritures au premier passage
	 init=.false.
	 write(2,1)
1	 format(//,t3,'------------------------------------'
     1	//,t5,'Constantes fondamentales utilisees',//)
 
	 write(2,2)msol,lsol,rsol,echarg,kbol,hpl,g,aradia,clight,eve,amu,me
2	 format(t4,'masse solaire=',1pd15.8,t45,'luminosite solaire=',
     1	1pd15.8,t85,'rayon solaire=',1pd15.8,/,
     2	t4,'charge de l''electron=',1pd15.8,t45,
     3	'constante de Boltzmann=',1pd15.8,t85,'constante de Planck=',
     4	1pd15.8,/,t4,'constante de la gravite=',1pd15.8,t45,
     5	'constante radiation=',1pd15.8,t85,'celerite de la lumiere=',
     6	1pd15.8,/,t4,'electron volt=',1pd15.8,t45,
     7	'masse atomique unite=',1pd15.8,t85,'masse de l''electron=',
     8	1pd15.8,//)
	 write(2,3)kih,kihe,kihe1
3	 format(t4,'potentiel d''ionisation de l''hydrogene=',1pd15.8,/,
     1	t4,'premier potentiel d''ionisation de l''helium=',1pd15.8,
     2	t65,'second potentiel d''ionisation de l''helium=',1pd15.8,//)
	 write(2,4)ah,ah2,ahe3,ahe4,ali7,abe7,ac12,ac13,an14
4	 format(t4,'masse atome d''H',1pd15.8,
     1	t45,'masse atome d''H2',1pd15.8,
     2	t85,'masse atome d''He3',1pd15.8,/,
     3	t4,'masse atome d''He4',1pd15.8,
     4	t45,'masse atome d''Li7',1pd15.8,
     5	t85,'masse atome d''Be7',1pd15.8,/,
     7	t4,'masse atome d''C12',1pd15.8,
     8	t45,'masse atome d''C13',1pd15.8,
     9	t85,'masse atome d''N14',1pd15.8)
 
	 write(2,5)an15,ao16,ao17,secon6*1.d-6
5	format(	t4,'masse atome d''N15',1pd15.8,
     1	t45,'masse atome d''O16',1pd15.8,
     2	t85,'masse atome d''O17',1pd15.8,//,
     3	t4,'nombre de seconde par annee tropique :',1pd15.8)
	 write(2,6)
6	 format(//,t40,'ATTENTION CONSTANTES ''CEPHE''',//)
	endif
 
	return
 
	end
