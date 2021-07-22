c**************************************************************************
 
	subroutine nrj_nuc_cn(t,ro,comp,dcomp,jac,deriv,fait,
     1	epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)
 

c Routine avec rapport Nc/Nn (rapport des nombres de particules) variable, 
c la fraction massique totale de C et N etant gardée constante 
c (cad egale a sa valeur solaire).

c Daniel, Juin 2000.

c	cycles PP, CNO et 3 alpha
 
c               Clayton p. 380, 392 et 430,
c               B. Pichon : cours d'Aussois p. 354 (1989)
c               A. Maeder A&A 120, 113-129 (1983)
c               G. Schaller et al A&A Suppl. Ser. 96, 269-231 (1992)
 
c	les tabulations de taux_reac
c	tables de Caughlan et Fowler 1988
 
c	Composition Chimique : Grevesse & Noels 1993
c	Rapports isotopiques : Maeder A&A, 120, 113 (1993)
 
c       Auteurs:
c          P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c          D. Cordier, ENSCR Campus de Beaulieu, Rennes
c		       DASGAL, Observatoire de Paris - Meudon
 
c	CESAM, Version 3
 
c entree :
c	t : temperature cgs
c	ro : densite cgs
c	comp : abondances
c	deriv=.true. : on calcule le jacobien
c	fait=1 : initialisation de la composition chimique
c	    =2 : calcul de dcomp et jacobien si deriv
c	    =3 : energie nucleaire et derivees / t et ro
c	    =4 : production de neutrinos  ( FICTIF ICI ! )
c	    =5 : calcul de dcomp et jacobien si deriv energie et derivee/ t, ro
 
c sorties
c	dcomp : derivee temporelle (unite de temps : 10**6 ans)
c	jac : jacobien (unite de temps : 10**6 ans)
c	e, et, ero, ex : energie thermonucleaire (unite de temps : s)
c			   : et derivees /t, ro ,X
c	hhe, be7e, b8e, n13e, o15e, f17e : nombre de neutrinos g/s
c	pour les reactions designees par le symbole (e pour electron +)
 
c initialisation de COMMON
 
c	/evol_chim_3/
c	nucleo : masse des noyaux
c	ab_min : abondance negligeable
c	nchim : nombre d'elements chimiques hors equilibre
c	nreac : nombre de reactions
c	nom_elem : symboles des elements chimiques utilises
 
c	r(1) : reaction H(H,e+ nu)H2			   PP
c	r(2) : reaction H2(H,g)H3
c	r(3) : reaction He3(He3,2H)He4
c	r(4) : reaction He4(He3,g)Be7
c	r(5) : reaction Li7(H,He4)He4
c	r(6) : reaction Be7(e-,nu g)Li7
c	r(7) : reaction Be7(H,g)B8(,e+ nu)Be8(,He4)He4
 
c	r(8) : reaction C12(H,g)N13(,e+ nu)C13	           CNO
c	r(9) : reaction C13(H,g)N14
c	r(10) : reaction N14(H,g)O15(e+,nu)N15
c	r(11) : reaction N15(H,g)O16
c	r(12) : reaction N15(H,He4)C12
c	r(13) : reaction O16(H,g)F17(,e+ nu)O17
c	r(14) : reaction O17(H,He4)N14
c	r(15) : reaction N13(H,g)O14(e+nu)N14   CNO <<chaud>>
c	r(16) : reaction O17(H,g)F18(e+nu)O18
c	r(17) : reaction O18(H,g)F19
c	r(18) : reaction F19(H,a)O16
c	r(19) : reaction F19(H,g)Ne20
c	r(20) : reaction O18(H,a)N15
 
c	r(21) : reaction He4(2He4,g)C12		           3 alpha
c	r(22) : reaction C12(He4,g)O16
c	r(23) : reaction O16(He4,g)Ne20
c	r(24) : reaction Ne20(a,g)Mg24
c	r(25) : reaction C13(a,n)O16
c	r(26) : reaction O17(a,n)Ne20
c	r(27) : reaction N14(a,g)F18(e+nu)O18
c	r(28) : reaction O18(a,g)Ne22
c	r(29) : reaction Ne22(a,n)Mg25
c	r(30) : reaction Ne22(a,g)Mg26
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'ctephy.common'
	include 'modele_3.common'
	include 'evol_chim_3.common'
 
	integer fait,i
 
        integer pel
 
	real*8 t,ro,dcomp(18),jac(18*18),dbe7y,dli7y,
     +  dbe7he3,dli7he3,
     1	r(pnreac),rt(pnreac),rro(pnreac),rx(pnreac),ddenx,
     2	q(pnreac),qt(pnreac),qro(pnreac),qx(pnreac),
     3	nel,epsilon(5),et,ero,ex(1),den,
     4	hhe,be7e,b8e,n13e,o15e,f17e,rxx(pnreac),dnelx,
     5	h2,dh2t,dh2ro,dh2x,be7,dbe7t,dbe7ro,dbe7x,
     6	li7,dli7t,dli7ro,dli7x,zs2
 
	real*8 abon_nbre(21), abon_nbre_suppl(18), comp(18),
     +         n_4He
 
	real*8 sss, n_pel, rapp_abon(36), log_A(21), log_A_suppl(18),
     +         masse_mol(18)
 
	real*8 Xc, Xn, Xcn, Mc, Mn, McsMn, NcsNn, Xc12, Xc13, Xn13, Xn14,
     +         Xn15

c	pour He3 on a inclus le D protosolaire 8.40d5 + 3.00d5
c	Le N13 est mis a 0.
c	Les neutrons sont mis a 0.
	
	logical deriv
 
	data log_A/12.000,1.160,1.150,2.600,7.301,10.990,8.5427,
     1         6.7718,-30.00,7.9686,5.4914,8.8689,5.6659,6.1710,
     2         4.5600,8.0289,7.1150,7.4760,6.5856,
     3         6.6281,-30.00/
 
        data log_A_suppl/5.5114,6.3300,6.4700,7.5500,5.4500,7.2100,
     1                   5.5000,6.5200,5.1200,6.3600,3.1700,5.0200,
     2                   4.0000,5.6700,5.3900,7.5000,4.9200,6.2500/
 
 
        data masse_mol/20.0000,22.989768,26.981539,28.0855,30.9738,
     1             32.0666,35.4528,39.9481,39.0983,40.0784,44.9559,
     2             47.8830,50.9415,51.9961,54.9381,55.8473,58.9332,
     3             58.6934/
 
	save
 
2000	format((1x,1p8d10.3))
 
 
	goto(100,200,300,400,200),fait
 
c---------------------------------------------------------------------------
c	initialisations
 
100	nchim=18	!nombre d'elements chimiques
	nreac=30	!nombre de reactions utilisees
 
	call taux_tab(.7d0,1.d7,1.d0,r,rt,rro,rx,rxx,q,qt,qro,qx,.false.)
 
	nucleo(1)=ah		!initialisation de nucleo: masses des noyaux
	nucleo(2)=ahe3
	nucleo(3)=ahe4
	nucleo(4)=ac12
	nucleo(5)=ac13
        nucleo(6)=an13
	nucleo(7)=an14
	nucleo(8)=an15
	nucleo(9)=ao16
	nucleo(10)=ao17
        nucleo(11)=ao18
        nucleo(12)=af19
        nucleo(13)=ane20
        nucleo(14)=ane22
        nucleo(15)=amg24
        nucleo(16)=amg25
        nucleo(17)=amg26
        nucleo(18)=an
	
	ihe4=3	!indice de He4
 
	nom_elem(1)= ' H '
	nom_elem(2)= 'He3'
	nom_elem(3)= 'He4'
	nom_elem(4)= 'C12'
	nom_elem(5)= 'C13'
        nom_elem(6)= 'N13'
	nom_elem(7)= 'N14'
	nom_elem(8)= 'N15'
	nom_elem(9)= 'O16'
	nom_elem(10)='O17'
        nom_elem(11)='O18'
        nom_elem(12)='F19'
        nom_elem(13)='N20'		!Ne
        nom_elem(14)='N22'
        nom_elem(15)='M24'
        nom_elem(16)='M25'
        nom_elem(17)='M26'		!Mg
        nom_elem(18)=' n '
c----------------------------------------------------------
 
c Initialisation de la composition chimique
 
c Attention : * le tableau 'comp' contient des fractions massiques divisees
c               par les masses molaires = nombre de moles !!
 
c             * H2, Li7 et Be7 ne sont pas dans 'comp'
c
 
c Initialisation de Z0 (se fait aussi dans ''resout_c2'') :
 
        z0=1.d0-x0-y0
 
c	Conversion des <<log_A=log Nelem/NH+12.00>> en <<Nelem>>
c	(abondance en nombre)
 
        do i = 1, 21
           abon_nbre(i)=10.**(log_A(i)-12.d0)
        end do
 
        do i= 1, 18
           abon_nbre_suppl(i)=10.**(log_A_suppl(i)-12.d0)
        end do
 
c        print*, 'abon_nbre= ', abon_nbre
 
c	Composition en hydrogene :
 
        comp(1)=x0/nucleo(1)
 
c	Composition en helium :
 
        rapp_abon(3)=abon_nbre(2+3)/abon_nbre(3+3)
 
        n_4He=y0/(nucleo(3)+rapp_abon(3)*nucleo(2))
 
        comp(3)=n_4He
 
        comp(2)=(y0-nucleo(3)*n_4He)/nucleo(2)
 
c	Composition des elements lourds :
c	---------------------------------
 
c determination des rapports d'abondances en nombre :
 
      pel=4 ! Premier element lourd, ici : C12
 
      do i= pel, 18, 1
         rapp_abon(i)=abon_nbre(i+3)/abon_nbre(pel+3)
      end do
 
      do i= 1, 18
         rapp_abon(i+18)=abon_nbre_suppl(i)/abon_nbre(pel+3)
      end do
 
c      print*, 'rapp_abon= ', rapp_abon
 
 
c Calcul de la somme <<sss>> des Mi (masses molaires) * rapp_abon
c (Nelem/N12C)
 
      sss=0.d0
 
      do i= pel, 18
         sss=sss+rapp_abon(i)*nucleo(i)
      end do
 
      do i= 1, 18
         sss=sss+rapp_abon(i+18)*masse_mol(i)
      end do
 
c Nombre de moles <<n_pel>> du premier element lourd
 
      n_pel=z0/sss
 
c Nombre de moles des differents elements lourds :
 
      do i= pel, 18
         comp(i)= n_pel*rapp_abon(i)
      end do
 
 
c Initialisation de 'ab_ini'
 
      do i= 1, nchim
         ab_ini(i)=comp(i)*nucleo(i)
      end do
 
c#######################################################################
c Calcul des nouvelles fractions massiques du Carbone et de l'Azote
c en tenant compte d'un rapport en nombre Nc/Nn non-solaire :

c fraction massique du carbone et de l'azote :

	Xc=comp(4)*nucleo(4)+comp(5)*nucleo(5)
	Xn=comp(6)*nucleo(6)+comp(7)*nucleo(7)+comp(8)*nucleo(8)

c la somme :

	Xcn=Xc+Xn

c Les masses molaires moyennes :

	Mc=Xc/(comp(4)+comp(5))
	Mn=Xn/(comp(6)+comp(7)+comp(8))

c leur rapport :

	McsMn=Mc/Mn

c le rapport des nombres de particules :

	print*, ' '
	print*, '#############################################'
	print*, ' '
	print*, '  Rapport des nbres Nc/Nn (carbone/azote) :'
        print*, ' '
        print*, '  (pour le Soleil on a 3.27, pour VX Pup 0.44)'
        print*, ' '
	read(5,*) NcsNn

c Fraction massique totale du carbone et de l'azote :

	Xc=Xcn*NcsNn*McsMn/(1.+NcsNn*McsMn)
	Xn=Xcn/(1.+NcsNn*McsMn)

c Les nouvelles fractions massiques des isotopes sont :

	Xc12=Xc*comp(4)/(comp(4)+comp(5))
	Xc13=Xc*comp(5)/(comp(4)+comp(5))

	Xn13=Xn*comp(6)/(comp(6)+comp(7)+comp(8))
	xn14=Xn*comp(7)/(comp(6)+comp(7)+comp(8))
	xn15=Xn*comp(8)/(comp(6)+comp(7)+comp(8))

c Vérification

	print*, 'xc12= ', xc12, ', comp(4)*Mc12= ', comp(4)*nucleo(4)
	print*, 'Xc13= ', Xc13, ', comp(5)*Mc13= ', comp(5)*nucleo(5)

        print*, 'Xn13= ', Xn13, ', comp(6)*Mn13= ', comp(6)*nucleo(6)
        print*, 'Xn14= ', Xn14, ', comp(7)*Mn14= ', comp(7)*nucleo(7)
        print*, 'Xn15= ', Xn15, ', comp(8)*Mn15= ', comp(8)*nucleo(8)
        
        print*, ' '
	print*, 'Ok ?'

	pause

c Affectation :
c pour le carbone :

	comp(4)=Xc12/nucleo(4)
	comp(5)=Xc13/nucleo(5)

c pour l'azote :

	comp(6)=Xn13/nucleo(6)
	comp(7)=Xn14/nucleo(7)
	comp(8)=xn15/nucleo(8)

c#######################################################################
c-----------------------------------------------------------------------
 
	ab_min(1)=1.d-10	!abondances negligeables
	ab_min(2)=5.d-10
	ab_min(3)=1.d-10
	ab_min(4)=5.d-10
	ab_min(5)=1.d-10
        ab_min(6)=3.d-23	!N13 non pris en compte pour l'instant
	ab_min(7)=1.d-10
	ab_min(8)=5.d-20
	ab_min(9)=1.d-10
	ab_min(10)=4.d-20
        ab_min(11)=2.d-20
        ab_min(12)=4.d-20
        ab_min(13)=2.d-20
        ab_min(14)=1.d-20
        ab_min(15)=5.d-20
        ab_min(16)=7.d-20
        ab_min(17)=8.d-20
        ab_min(18)=1.d-25
 
	write(2,*)' '
	write(2,*)'Reactions thermonucleaires des cycles PP, CNO, 3 Alpha'
	write(2,*)'par interpolation des tables de Caughlan et Fowler 1988'
        write(2,*)'et des taux de V. Landre et al'
	write(2,*)' '
	write(2,*)'nombre de reactions : ',nreac
	write(2,*)'nombre d''elements chimiques : ',nchim
	write(2,*)' '
	write(2,20)
20	format(1x,'les abondances initiales sont deduites de X suivant',/,
     1	1x,'Grevesse & Noels 1993,',/,
     2	1x,'mutipliees par les rapports isotopiques pris par Maeder 1993'
     3  )
	write(2,*)		
	write(2,1)(comp(i)*nucleo(i),i=1,nchim)
1	format(3x,' H : ',1pd10.3,3x,' He3 : ',1pd10.3,3x,' He4 : ',1pd10.3,3x,
     1	' C12 : ',1pd10.3,3x,' C13 : ',1pd10.3,
     2	/,1x,' N13 : ',1pd10.3,3x,' N14 : ',1pd10.3,3x,
     3	' N15 : ',1pd10.3,3x,' O16 : ',1pd10.3,3x,' O17 : ',1pd10.3,
     3  /,1x,
     4  ' O18 : ',1pd10.3,3x,' F19 : ',1pd10.3,3x,' Ne20: ',1pd10.3,3x,
     5  ' Ne22: ',1pd10.3,3x,' Mg24: ',1pd10.3,
     6  /,1x,' Mg25: ',1pd10.3,3x,
     6  ' Mg26: ',1pd10.3,3x,' n   : ',1pd10.3)
 
	write(2,*)' '
	write(2,*)'abondances negligeables:'
	write(2,1)(ab_min(i),i=1,nchim)
	write(2,*)' '
	write(2,*)'H2, Li7, Be7 sont pris a l''equilibre'
	write(2,*)'on utilise une table'
	write(2,*)' '
	write(2,*)'pour l''evolution temporelle, test de precision sur H et He4'
	write(2,*)' '
	
	write(6,*)' '
	write(6,*)'Reactions thermonucleaires des cycles PP, CNO, 3 Alpha'
	write(6,*)'par interpolation des tables de Caughlan et Fowler 1988'
	write(6,*)' '
	write(6,*)'nombre de reactions : ',nreac
	write(6,*)'nombre d''elements chimiques : ',nchim
	write(6,*)' '	
	write(6,20)
	write(6,*)	
	write(6,1)(comp(i)*nucleo(i),i=1,nchim)
	write(6,*)' '
	write(6,*)'abondances negligeables:'
	write(6,1)(ab_min(i),i=1,nchim)
	write(6,*)' '
	write(6,*)'H2, Li7, Be7 sont pris a l''equilibre'
	write(6,*)'on utilise une table'
	write(6,*)' '
	write(6,*)'pour l''evolution temporelle, test de precision sur H et He4'
	write(6,*)' '
	write(2,*)' '
 
	do i=1,nchim
	 ab_min(i)=ab_min(i)/nucleo(i)
	enddo
 
	zs2=8./16.*z0		!on ne suppose pas Z diffuse
	dnelx=1.-2.*ah/ahe4	!d nel / dX en supposant Y=1-X-Z
 
	return
 
c-------------------------------------------------------------------------
c	reactions
 
200	if(t .lt. t_inf)then	!si t<t_inf
	 do i=1,nbelem
	  dcomp(i)=0.
	  ex(i)=0
	 enddo
	 do i=1,nbelem*nbelem
	  jac(i)=0
	 enddo
	 do i=1,4
	  epsilon(i)=0
	 enddo
	 et=0
	 ero=0
	 return
	endif
	call taux_tab(comp(1),t,ro,r,rt,rro,rx,rxx,q,qt,qro,qx,deriv)
 
c	write(6,*)'comp'
c	write(6,2000)(comp(i),i=1,nchim)
c	write(6,*)'reactions'
c	write(6,2000)(r(i),i=1,nreac)
 
c************************************************************************
 
c	equations d'evolution
 
c
c equation de H
 
      dcomp(1)=-(3.*r(1)*comp(1)+r(8)*comp(4)+r(9)*comp(5)+r(10)*comp(7)
     +   +(r(11)+r(12))*comp(8)+r(13)*comp(9)
     +   +(r(14)+r(16))*comp(10)+(r(17)+r(20))*comp(11)
     +   +(r(18)+r(19))*comp(12))*comp(1)+(2.*r(3)*comp(2)
     +   -r(4)*comp(3))*comp(2)
 
 
c equation de He3 (elem. no 2)
 
      dcomp(2)=r(1)*comp(1)**2.-(2.*r(3)*comp(2)+r(4)*comp(3))*comp(2)
 
 
c equation de He4 (elem. no 3)
 
      dcomp(3)=(r(3)*comp(2)+r(4)*comp(3))*comp(2)
     +   +(r(12)*comp(8)+r(14)*comp(10)+r(18)*comp(12)+r(20)*comp(11))
     +   *comp(1)
     +   -(3.*r(21)*comp(3)**2.+r(22)*comp(4)+r(23)*comp(9)+r(24)
     +   *comp(13)
     +   +r(25)*comp(5)
     +   +r(26)*comp(10)+r(27)*comp(7)+r(28)*comp(11)+(r(29)+r(30))
     +   *comp(14))*comp(3)
 
 
c equation de C12 (elem. no 4)
 
      dcomp(4)=r(12)*comp(1)*comp(8)+r(21)*comp(3)**3.
     +   -(r(8)*comp(4)*comp(1)+r(22)*comp(4)*comp(3))
 
 
c equation de C13 (elem. no 5)
 
      dcomp(5)=(r(8)*comp(4)-r(9)*comp(5))*comp(1)
     +   -r(25)*comp(5)*comp(3)
 
 
c equation de N13 (elem. no 6)
 
      dcomp(6)=0.
c     le CNO chaud n'est pas encore pris en compte!
 
 
c equation de N14 (elem. no 7)
 
      dcomp(7)=(r(9)*comp(5)+r(14)*comp(10))*comp(1)
     +   -(r(10)*comp(1)+r(27)*comp(3))*comp(7)
 
 
c equation de N15 (elem. no 8)
 
      dcomp(8)=(r(10)*comp(7)+r(20)*comp(11)-(r(11)+r(12))*comp(8))
     +   *comp(1)
 
 
c equation de O16 (elem. no 9)
 
      dcomp(9)=(r(11)*comp(8)+r(18)*comp(12))*comp(1)+(r(22)*comp(4)
     +   +r(25)*comp(5))*comp(3)
     +   -(r(23)*comp(3)+r(13)*comp(1))*comp(9)
 
 
c equation de O17 (elem. no 10)
 
      dcomp(10)=r(13)*comp(9)*comp(1)-((r(16)+r(14))*comp(1)+
     +   r(26)*comp(3))*comp(10)
 
 
c equation de O18 (elem. no 11)
 
      dcomp(11)=r(16)*comp(1)*comp(10)+r(27)*comp(3)*comp(7)
     +   -((r(17)+r(20))*comp(1)*comp(11)+r(28)*comp(3)*comp(11))
 
 
c equation de F19 (elem. no 12)
 
      dcomp(12)=(r(17)*comp(11)-(r(18)+r(19))*comp(12))*comp(1)
 
 
c equation de Ne20 (elem. no 13)
 
      dcomp(13)=r(19)*comp(1)*comp(12)+(r(23)*comp(9)+r(26)*comp(10))
     +   *comp(3)
     +   -r(24)*comp(3)*comp(13)
 
 
c equation de Ne22 (elem. no 14)
 
      dcomp(14)=r(28)*comp(3)*comp(11)-(r(29)+r(30))*comp(3)*comp(14)
 
 
c equation de Mg24 (elem. no 15)
 
      dcomp(15)=r(24)*comp(3)*comp(13)
 
 
c equation de Mg25 (elem. no 16)
 
      dcomp(16)=r(29)*comp(3)*comp(14)
 
 
c equation de Mg26 (elem. no 17)
 
      dcomp(17)=r(30)*comp(3)*comp(14)
 
 
c equation des neutrons (''elem.'' no 18)
 
      dcomp(18)=(r(25)*comp(5)+r(26)*comp(10)+r(29)*comp(14))*comp(3)
 
 
c--------------------------------------------------------------------------	
c	pour le moment angulaire ou Z
	
	if(iw .gt. 0)dcomp(iw)=0.	!pour MA
	if(iz .gt. 0)dcomp(iz)=0.	!pour Z		
 
c	write(6,*)'somme << dcomp / dcomp'
c	sum=0.
c	do i=1,nchim
c	 sum=sum+nucleo(i)*dcomp(i)
c	enddo
c	write(6,2000)sum
c	write(6,2000)(dcomp(i),i=1,nchim)
c	write(6,*)'reac,rx'
c	write(6,2000)(r(i),i=1,nreac)
c	write(6,2000)(rx(i),i=1,nreac)
 
 
c-----------------------------------------------------------------------------
c	calcul du jacobien
 
	if(deriv .or. fait .eq. 6)then	!jacp(i,j) : equation, j : element i
 
c	a travers l'effet d'ecran, les reactions dependent de H par rx
 
	 do i=1,nbelem*nbelem
	  jac(i)=0.
	 enddo		!i
 
 
c       print*, 'valeur de nbelem dans energie_nuc : ', nbelem
 
c**************
c Jacobien de H
 
c H
 
      jac(nbelem*(1-1)+1)=-(3.*rx(1)*comp(1)
     +   +rx(8)*comp(4)+rx(9)*comp(5)+rx(10)*comp(7)+(rx(11)+rx(12))
     +   *comp(8)+rx(13)*comp(9)
     +   +(rx(14)+rx(16))*comp(10)+
     +   +(rx(17)+rx(20))*comp(11)
     +   +(rx(18)+rx(19))*comp(12))*comp(1)
     +   -(6.*r(1)*comp(1)+r(8)*comp(4)+r(9)*comp(5)+r(10)*comp(7)+
     +   (r(11)+r(12))*comp(8)+r(13)*comp(9)+(r(14)+r(16))*comp(10)
     +   +(r(17)+r(20))*comp(11)+(r(18)+r(19))*comp(12)
     +   )+(2.*rx(3)*comp(2)-rx(4)*comp(3))*comp(2)
 
c He3
 
      jac(nbelem*(2-1)+1)=4.*r(3)*comp(2)-r(4)*comp(3)
 
c He4
 
      jac(nbelem*(3-1)+1)=-r(4)*comp(2)
 
c C12
 
      jac(nbelem*(4-1)+1)=-r(8)*comp(1)
 
c C13
 
      jac(nbelem*(5-1)+1)=-r(9)*comp(1)
 
c N13 ( ATTENTION UNIQUEMENT POUR LE CNO ''CHAUD'' )
 
c      jac(nbelem*(6-1)+1)=-r(15)*comp(1)
      jac(nbelem*(6-1)+1)=0.
 
c N14
 
      jac(nbelem*(7-1)+1)=-r(10)*comp(1)
 
c N15
 
      jac(nbelem*(8-1)+1)=-(r(11)+r(12))*comp(1)
 
c O16
 
      jac(nbelem*(9-1)+1)=-r(13)*comp(1)
 
c O17
 
      jac(nbelem*(10-1)+1)=-(r(14)+r(16))*comp(1)
 
c O18
 
      jac(nbelem*(11-1)+1)=-(r(17)+r(20))*comp(1)
 
c F19
 
      jac(nbelem*(12-1)+1)=-(r(18)+r(19))*comp(1)
 
 
c*************
c Jacobien He3
 
c H
 
      jac(nbelem*(1-1)+2)=rx(1)*comp(1)**2.+r(1)*2.*comp(1)-(2.*rx(3)
     +   *comp(2)+rx(4)*comp(3))*comp(2)
 
c He3
 
      jac(nbelem*(2-1)+2)=-(4.*r(3)*comp(2)+r(4)*comp(3))
 
c He4
 
      jac(nbelem*(3-1)+2)=-r(4)*comp(2)
 
 
c****************
c Jacobien de He4
 
c H
 
      jac(nbelem*(1-1)+3)=(rx(3)*comp(2)+rx(4)*comp(3))*comp(2)
     +   +(rx(12)*comp(8)+rx(14)*comp(10)+rx(18)*comp(12)+rx(20)
     +   *comp(11))*comp(1)
     +   +(r(12)*comp(8)+r(14)*comp(10)+r(18)*comp(12)+r(20)*comp(11))
     +   -(3.*rx(21)*comp(3)**2.+rx(22)*comp(4)+rx(23)*comp(9)
     +   +rx(24)*comp(13)+rx(25)*comp(5)
     +   +rx(26)*comp(10)+rx(27)*comp(7)+rx(28)*comp(11)
     +   +(rx(29)+rx(30))*comp(14))*comp(3)
 
c He3
 
      jac(nbelem*(2-1)+3)=2.*r(3)*comp(2)+r(4)*comp(3)
 
c He4
 
      jac(nbelem*(3-1)+3)=r(4)*comp(2)
     +   -(9.*r(21)*comp(3)**2.+r(22)*comp(4)+r(23)*comp(9)+r(24)
     +   *comp(13)+r(25)*comp(5)
     +   +r(26)*comp(10)+r(27)*comp(7)+r(28)*comp(11)+(r(29)+r(30))
     +   *comp(14))
 
c C12
 
      jac(nbelem*(4-1)+3)=-r(22)*comp(3)
 
c C13
 
      jac(nbelem*(5-1)+3)=-r(25)*comp(3)
 
c N14
 
      jac(nbelem*(7-1)+3)=-r(27)*comp(3)
 
c N15
 
      jac(nbelem*(8-1)+3)=r(12)*comp(1)
 
c O16
 
      jac(nbelem*(9-1)+3)=-r(23)*comp(3)
 
c O17
 
      jac(nbelem*(10-1)+3)=r(14)*comp(1)-r(26)*comp(3)
 
c O18
 
      jac(nbelem*(11-1)+3)=r(20)*comp(1)-r(28)*comp(3)
 
c F19
 
      jac(nbelem*(12-1)+3)=r(18)*comp(1)
 
c Ne20
 
      jac(nbelem*(13-1)+3)=-r(24)*comp(3)
 
c Ne22
 
      jac(nbelem*(14-1)+3)=-(r(29)+r(30))*comp(3)
 
 
c****************
c Jacobien de C12
 
c H
 
      jac(nbelem*(1-1)+4)=rx(12)*comp(1)*comp(8)+rx(21)*comp(3)**3.
     +   -(rx(8)*comp(4)*comp(1)+rx(22)*comp(4)*comp(3))
     +   +r(12)*comp(8)-r(8)*comp(4)
 
c He4
 
      jac(nbelem*(3-1)+4)=3.*r(21)*comp(3)**2.-r(22)*comp(4)
 
c C12
 
      jac(nbelem*(4-1)+4)=-(r(8)*comp(1)+r(22)*comp(3))
 
c N15
 
      jac(nbelem*(8-1)+4)=r(12)*comp(1)
 
 
c****************
c Jacobien de C13
 
c H
 
      jac(nbelem*(1-1)+5)=(rx(8)*comp(4)-rx(9)*comp(5))*comp(1)-rx(25)
     +   *comp(5)*comp(3)
     +   +r(8)*comp(4)-r(9)*comp(5)
 
c He4
 
      jac(nbelem*(3-1)+5)=-r(25)*comp(5)
 
c C12
 
      jac(nbelem*(4-1)+5)=r(8)*comp(1)
 
c C13
 
      jac(nbelem*(5-1)+5)=-r(9)*comp(1)-r(25)*comp(3)
 
 
c****************
c Jacobien de N13   ( FACTICE TANT QUE LE CNO CHAUD N'EST PAS PRIS EN
c                    COMPTE! )
 
c H
 
      jac(nbelem*(1-1)+6)=0.
 
 
c****************
c Jacobien de N14
 
c H
      jac(nbelem*(1-1)+7)=(rx(9)*comp(5)+rx(14)*comp(10))*comp(1)
     +   -(rx(10)*comp(1)+rx(27)*comp(3))*comp(7)
     +   +r(9)*comp(5)+r(14)*comp(10)-r(10)*comp(7)
 
c He4
 
      jac(nbelem*(3-1)+7)=-r(27)*comp(7)
 
c C13
 
      jac(nbelem*(5-1)+7)=r(9)*comp(1)
 
c N14
 
      jac(nbelem*(7-1)+7)=-(r(10)*comp(1)+r(27)*comp(3))
 
c O17
 
      jac(nbelem*(10-1)+7)=r(14)*comp(1)
 
 
c****************
c Jacobien de N15
 
c H
 
      jac(nbelem*(1-1)+8)=(rx(10)*comp(7)+rx(20)*comp(11)-(rx(11)
     +   +rx(12))*comp(8))*comp(1)
     +   +r(10)*comp(7)+r(20)*comp(11)-(r(11)+r(12))*comp(8)
 
c N14
 
      jac(nbelem*(7-1)+8)=r(10)*comp(1)
 
c N15
 
      jac(nbelem*(8-1)+8)=-(r(11)+r(12))*comp(1)
 
c O18
 
      jac(nbelem*(11-1)+8)=r(20)*comp(1)
 
 
c****************
c Jacobien de O16
 
c H
 
      jac(nbelem*(1-1)+9)=(rx(11)*comp(8)+rx(18)*comp(12))*comp(1)
     +   +(rx(22)*comp(4)+rx(25)*comp(5))*comp(3)
     +   -(rx(23)*comp(3)+rx(13)*comp(1))*comp(9)
     +   +(r(11)*comp(8)+r(18)*comp(12))-r(13)*comp(9)
 
c He4
 
      jac(nbelem*(3-1)+9)=r(22)*comp(4)+r(25)*comp(5)-r(23)*comp(9)
 
c C12
 
      jac(nbelem*(4-1)+9)=r(22)*comp(3)
 
c C13
 
      jac(nbelem*(5-1)+9)=r(25)*comp(3)
 
c N15
 
      jac(nbelem*(8-1)+9)=r(11)*comp(1)
 
c O16
 
      jac(nbelem*(9-1)+9)=-(r(23)*comp(3)+r(13)*comp(1))
 
c F19
 
      jac(nbelem*(12-1)+9)=r(18)*comp(1)
 
 
c****************
c Jacobien de O17
 
c H
 
      jac(nbelem*(1-1)+10)=rx(13)*comp(9)*comp(1)-((rx(16)+rx(14))
     +   *comp(1)+
     +   +rx(26)*comp(3))*comp(10)
     +   +r(13)*comp(9)-r(16)*comp(10)
 
c He4
 
      jac(nbelem*(3-1)+10)=-r(26)*comp(10)
 
c 016
 
      jac(nbelem*(9-1)+10)=r(13)*comp(1)
 
c O17
 
      jac(nbelem*(10-1)+10)=-((r(16)+r(14))*comp(1)+r(26)*comp(3))
 
 
c****************
c Jacobien de O18
 
c H
 
      jac(nbelem*(1-1)+11)=rx(16)*comp(1)*comp(10)+rx(27)*comp(3)
     +   *comp(7)
     +   -((rx(17)+rx(20))*comp(1)*comp(11)+rx(28)*comp(3)*comp(11))
     +   +r(16)*comp(10)-(r(17)+r(20))*comp(11)
 
c He4
 
      jac(nbelem*(3-1)+11)=r(27)*comp(7)-r(28)*comp(11)
 
c N14
 
      jac(nbelem*(7-1)+11)=r(27)*comp(3)
 
c O17
 
      jac(nbelem*(10-1)+11)=r(16)*comp(1)
 
c O18
 
      jac(nbelem*(11-1)+11)=-((r(17)+r(20))*comp(1)+r(28)*comp(3))
 
 
c****************
c Jacobien de F19
 
c H
 
      jac(nbelem*(1-1)+12)=(rx(17)*comp(11)-(rx(18)+rx(19))*comp(12))
     +   *comp(1)
     +   +r(17)*comp(11)-(r(18)+r(19))*comp(12)
 
c O18
 
      jac(nbelem*(11-1)+12)=r(17)*comp(1)
 
c F19
 
      jac(nbelem*(12-1)+12)=-(r(18)+r(19))*comp(1)
 
 
c*****************
c Jacobien de Ne20
 
c H
 
      jac(nbelem*(1-1)+13)=rx(19)*comp(1)*comp(12)+(rx(23)*comp(9)
     +   +rx(26)*comp(10))*comp(3)
     +   -rx(24)*comp(3)*comp(13)+r(19)*comp(12)
 
c He4
 
      jac(nbelem*(3-1)+13)=r(23)*comp(9)+r(26)*comp(10)-r(24)*comp(13)
 
c O16
 
      jac(nbelem*(9-1)+13)=r(23)*comp(3)
 
c O17
 
      jac(nbelem*(10-1)+13)=r(26)*comp(3)
 
c F19
 
      jac(nbelem*(12-1)+13)=r(19)*comp(1)
 
c Ne20
 
      jac(nbelem*(13-1)+13)=-r(24)*comp(3)
 
 
c*****************
c Jacobien de Ne22
 
c H
 
      jac(nbelem*(1-1)+14)=rx(28)*comp(3)*comp(11)-(rx(29)+rx(30))
     +   *comp(3)*comp(14)
 
c He4
 
      jac(nbelem*(3-1)+14)=r(28)*comp(11)-(r(29)+r(30))*comp(14)
 
c O18
 
      jac(nbelem*(11-1)+14)=r(28)*comp(3)
 
c Ne22
 
      jac(nbelem*(14-1)+14)=-(r(29)+r(30))*comp(3)
 
 
c*****************
c Jacobien de Mg24
 
c H
 
      jac(nbelem*(1-1)+15)=rx(24)*comp(3)*comp(13)
 
c He4
 
      jac(nbelem*(3-1)+15)=r(24)*comp(13)
 
c Ne20
 
      jac(nbelem*(13-1)+15)=r(24)*comp(13)
 
c*****************
c Jacobien de Mg25
 
c H
 
      jac(nbelem*(1-1)+16)=rx(29)*comp(3)*comp(14)
 
c He4
 
      jac(nbelem*(3-1)+16)=r(29)*comp(14)
 
c Ne22
 
      jac(nbelem*(14-1)+16)=r(29)*comp(3)
 
 
c*****************
c Jacobien de Mg26
 
c H
 
      jac(nbelem*(1-1)+17)=rx(30)*comp(3)*comp(14)
 
c He4
 
      jac(nbelem*(3-1)+17)=r(30)*comp(14)
 
c Ne22
 
      jac(nbelem*(14-1)+17)=r(30)*comp(3)
 
 
c**********************
c Jacobien des neutrons
 
c H
 
      jac(nbelem*(1-1)+18)=(rx(25)*comp(5)+rx(26)*comp(10)+rx(29)
     +   *comp(14))*comp(3)
 
c He4
 
      jac(nbelem*(3-1)+18)=r(25)*comp(5)+r(26)*comp(10)+r(29)*comp(14)
 
c C13
 
      jac(nbelem*(5-1)+18)=r(25)*comp(3)
 
c 017
 
      jac(nbelem*(10-1)+18)=r(26)*comp(3)
 
c Ne22
 
      jac(nbelem*(14-1)+18)=r(29)*comp(3)
 
 
c----------------------------------------------------------------------------
c	 pour le moment angulaire ou Z jac=0
 
 
 
c----------------------------------------------------------------------------
c	 unites de temps pour integration temporelle
 
	 do i=1,nbelem*nbelem
	  jac(i)=jac(i)*secon6
	 enddo
	endif		!deriv
 
	do i=1,nbelem
	 dcomp(i)=dcomp(i)*secon6
	enddo
 
	if(fait .ne. 5)return
 
c-----------------------------------------------------------------------------
c	calcul de la production d'energie nucleaire et derivees
c	pour H2(H,g)He3, q(2)H**2=q(2)*r(1)/r(2)
 
300	if(t .ge. t_inf)then
	 if(fait .ne. 5)call taux_tab(comp(1),t,ro,r,rt,
     1	rro,rx,rxx,q,qt,qro,qx,deriv)
 
c------------------------------------------------------------------------------
c	 nombre d'electrons / mole/g Y est compte comme 1-X-Z
 
	 nel=comp(1)+2.*(1.-comp(1)*ah-z0)/ahe4+zs2
 
	 if(comp(1) .ne. 0.d0)then
	  h2=r(1)/r(2)*comp(1)
	  den=r(6)*nel+r(7)*comp(1)
	  be7=r(4)*comp(2)*comp(3)/den
	  li7=r(6)*be7*nel/r(5)/comp(1)
	 else
	  h2=0.d0
	  be7=0.d0
	  li7=0.d0
	 endif
 
	 epsilon(1)=0.
 
	 epsilon(2)=(q(1)*comp(1)+q(2)*h2+q(5)*li7+q(7)*be7)*comp(1)
     1	+(q(3)*comp(2)+q(4)*comp(3))*comp(2)+q(6)*nel*be7
 
	 epsilon(3)=(q(8)*comp(4)+q(9)*comp(5)+q(10)*comp(7)
     +              +(q(11)+q(12))*comp(8)+q(13)*comp(9)+q(14)*comp(10)
     +              +q(16)*comp(10)+q(17)*comp(11)
     +              +(q(18)+q(19))*comp(12)+q(20)*comp(11))*comp(1)
 
	 epsilon(4)=(q(21)*comp(3)**2+q(22)*comp(4)+q(23)*comp(9)
     +              +q(24)*comp(13)+q(25)*comp(5)+q(26)*comp(10)
     +              +q(27)*comp(7)+q(28)*comp(11)+(q(29)+q(30))
     +              *comp(14))*comp(3)
 
 
	 do i=2,4
	  epsilon(1)=epsilon(1)+epsilon(i)
	 enddo
 
c	 write(6,2000)(epsilon(i),i=1,4)
c	 write(6,2000)(comp(i),i=1,nchim)
c	 write(6,2000)(q(i),i=1,nreac)
c	 write(6,2000)nel,h2,li7,be7
c	 write(6,2000)q(1)*comp(1)*comp(1)
c	 write(6,2000)q(2)*h2*comp(1)
c	 write(6,2000)q(5)*li7*comp(1)
c	 write(6,2000)q(7)*be7*comp(1)
c	 write(6,2000)q(3)*comp(2)*comp(2)
c	 write(6,2000)q(4)*comp(3)*comp(2)
c	 write(6,2000)q(6)*nel*be7*comp(1)
c	 write(6,2000)q(8)*comp(4)*comp(1)
c	 write(6,2000)q(9)*comp(5)*comp(1)
c	 write(6,2000)q(10)*comp(6)*comp(1)
c	 write(6,2000)q(11)*comp(7)*comp(1)
c	 write(6,2000)q(12)*comp(7)*comp(1)
c	 write(6,2000)q(13)*comp(8)*comp(1)
c	 write(6,2000)q(14)*comp(9)*comp(1)
c	 write(6,2000)q(15)*comp(3)**2*comp(3)
c	 write(6,2000)q(16)*comp(4)*comp(3)
c	 write(6,2000)q(17)*comp(8)*comp(3)
c	 pause
 
	 if(deriv)then
	  if(h2 .ne. 0.d0)then
	   dh2t=h2*(rt(1)/r(1)-rt(2)/r(2))
	   dh2ro=h2*(rro(1)/r(1)-rro(2)/r(2))
	   dh2x=h2*(rx(1)/r(1)-rx(2)/r(2)+1./comp(1))
	  else
	   dh2t=0.
	   dh2ro=0.
	   dh2x=0.
	  endif
 
	  if(be7 .ne. 0.d0)then
	   ddenx=r(6)*dnelx+r(7)
	   dbe7t=be7*(rt(4)/r(4)-(rt(6)*nel+rt(7)*comp(1))/den)
	   dbe7ro=be7*(rro(4)/r(4)-(rro(6)*nel+rro(7)*comp(1))/den)
	   dbe7x=be7*(rx(4)/r(4)-(rx(6)*nel+rx(7)*comp(1)+ddenx)/den)
	   dbe7y=be7/comp(3)
	   dbe7he3=be7/comp(2)
	  else
	   dbe7t=0.
	   dbe7ro=0.
	   dbe7x=0.
	   dbe7y=0.
	   dbe7he3=0.
	  endif
 
	  if(li7 .ne. 0.d0)then
	   dli7t=li7*(rt(6)/r(6)+dbe7t/be7-rt(5)/r(5))
	   dli7ro=li7*(rro(6)/r(6)+dbe7ro/be7-rro(5)/r(5))
	   dli7x=li7*(rx(6)/r(6)+dbe7x/be7+dnelx/nel-rx(5)/r(5)-1./comp(1))
	   dli7y=li7*dbe7y/be7
	   dli7he3=li7*dbe7he3/be7
	  else
	   dli7t=0.
	   dli7ro=0.
	   dli7x=0.
	   dli7y=0.
	   dli7he3=0.
	  endif
 
	  et=(qt(1)*comp(1)+qt(2)*h2+q(2)*dh2t+qt(5)*li7+q(5)*dli7t
     1	     +qt(7)*be7+q(7)*dbe7t)*comp(1)
     2	     +(qt(3)*comp(2)+qt(4)*comp(3))*comp(2)
     3	     +nel*(qt(6)*be7+q(6)*dbe7t)+
     +       +(qt(8)*comp(4)+qt(9)*comp(5)+qt(10)*comp(7)
     +       +(qt(11)+qt(12))*comp(8)+qt(13)*comp(9)+qt(14)*comp(10)
     +       +qt(16)*comp(10)+qt(17)*comp(11)
     +       +(qt(18)+qt(19))*comp(12)+qt(20)*comp(11))*comp(1)
     +       +(qt(21)*comp(3)**2.+qt(22)*comp(4)+qt(23)*comp(9)
     +       +qt(24)*comp(13)+qt(25)*comp(5)+qt(26)*comp(10)
     +       +qt(27)*comp(7)+qt(28)*comp(11)+qt(29)*comp(14)
     +       +qt(30)*comp(14))*comp(3)
 
 
	  ero=(qro(1)*comp(1)+qro(2)*h2+q(2)*dh2ro+qro(5)*li7+q(5)*dli7ro
     1	     +qro(7)*be7+q(7)*dbe7ro)*comp(1)
     2	     +(qro(3)*comp(2)+qro(4)*comp(3))*comp(2)
     3	     +nel*(qro(6)*be7+q(6)*dbe7ro)
     +       +(qro(8)*comp(4)+qro(9)*comp(5)+qro(10)*comp(7)
     +       +(qro(11)+qro(12))*comp(8)+qro(13)*comp(9)+qro(14)*comp(10)
     +       +qro(16)*comp(10)+qro(17)*comp(11)
     +       +(qro(18)+qro(19))*comp(12)+qro(20)*comp(11))*comp(1)
     +       +(qro(21)*comp(3)**2.+qro(22)*comp(4)+qro(23)*comp(9)
     +       +qro(24)*comp(13)+qro(25)*comp(5)+qro(26)*comp(10)
     +       +qro(27)*comp(7)+qro(28)*comp(11)+qro(29)*comp(14)
     +       +qro(30)*comp(14))*comp(3)
 
 
	  ex(1)=(qx(1)*comp(1)+qx(2)*h2+q(2)*dh2x+qx(5)*li7+q(5)*dli7x
     1	+qx(7)*be7+q(7)*dbe7x)*comp(1)
     2	+2.*q(1)*comp(1)+q(2)*h2+q(5)*li7+q(7)*be7
     3	+(qx(3)*comp(2)+qx(4)*comp(3))*comp(2)+q(6)*dnelx*be7
     4	+nel*(qx(6)*be7+q(6)*dbe7x)
     +          +(qx(8)*comp(4)+qx(9)*comp(5)+qx(10)*comp(7)
     +          +(qx(11)+qx(12))*comp(8)+qx(13)*comp(9)+qx(14)*comp(10)
     +          +qx(16)*comp(10)+qx(17)*comp(11)
     +          +(qx(18)+qx(19))*comp(12)+qx(20)*comp(11))*comp(1)
     +          +q(8)*comp(4)+q(9)*comp(5)+q(10)*comp(7)
     +          +(q(11)+q(12))*comp(8)+q(13)*comp(9)+q(14)*comp(10)
     +          +q(16)*comp(10)+q(17)*comp(11)
     +          +(q(18)+q(19))*comp(12)+q(20)*comp(11)
     +          +(qx(21)*comp(3)**2.+qx(22)*comp(4)+qx(23)*comp(9)
     +          +qx(24)*comp(13)+qx(25)*comp(5)+qx(26)*comp(10)
     +          +qx(27)*comp(7)+qx(28)*comp(11)+qx(29)*comp(14)
     +          +qx(30)*comp(14))*comp(3)
 
 
	  ex(2)=(q(5)*dli7he3+q(7)*dbe7he3)*comp(1)+2.*q(3)*comp(2)
     1	+q(4)*comp(3)+q(6)*nel*dbe7he3
 
	  ex(3)=(q(5)*dli7y+q(7)*dbe7y)*comp(1)+q(4)*comp(2)+q(6)*nel*dbe7y
     +          +3.*q(21)*comp(3)**2.+q(22)*comp(4)+q(23)*comp(9)
     +          +q(24)*comp(13)+q(25)*comp(5)+q(26)*comp(10)
     +          +q(27)*comp(7)+q(28)*comp(11)+q(29)*comp(14)
     +          +q(30)*comp(14)
 
	  ex(4)=q(8)*comp(1)+q(22)*comp(3)
 
	  ex(5)=q(9)*comp(1)+q(25)*comp(3)
 
	  ex(6)=0. ! r(15) : CNO chaud, q0(15)=0.
 
	  ex(7)=q(10)*comp(1)+q(27)*comp(3)
 
	  ex(8)=(q(11)+q(12))*comp(1)
 
	  ex(9)=q(13)*comp(1)+q(23)*comp(3)
 
          ex(10)=(q(14)+q(16))*comp(1)+q(26)*comp(3)
 
          ex(11)=(q(17)+q(20))*comp(1)+q(28)*comp(3)
 
          ex(12)=(q(18)+q(19))*comp(1)
 
          ex(13)=q(24)*comp(3)
 
          ex(14)=(q(29)+q(30))*comp(3)
 
          ex(15)=0.
 
          ex(16)=0.
 
          ex(17)=0.
 
          ex(18)=0.
	
	 endif	!deriv
 
	else	!t<t_inf
	 do i=1,4
	  epsilon(i)=0.
	 enddo
	 et=0.
	 ero=0.
	 do i=1,nbelem
	  ex(i)=0.
	 enddo
	endif
 
	return
 
c--------------------------------------------------------------------------
c	production de neutrinos     ( Ici aucune amelioration par rapport
c                                     a 'ppcno3a_9_3' )
 
400	if(t .ge. t_inf)then
	 call taux_tab(comp(1),t,ro,r,rt,rro,rx,rxx,q,qt,qro,qx,.false.)
 
	 be7=r(4)*comp(2)*comp(3)/(r(6)*nel+r(7)*comp(1))
 
	 hhe=r(1)*comp(1)**2/amu
	 be7e=r(6)*nel*be7/amu
	 b8e=r(7)*comp(1)*be7/amu
	 n13e=r(8)*comp(1)*comp(4)/amu
	 o15e=r(10)*comp(1)*comp(6)/amu
	 f17e=r(13)*comp(1)*comp(8)/amu
	else
	 hhe=0.
	 be7e=0.
	 b8e=0.
	 n13e=0.
	 o15e=0.
	 f17e=0.
	endif
 
	return
 
	end
 
 
