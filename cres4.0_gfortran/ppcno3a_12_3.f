 
c**************************************************************************
 
	subroutine ppcno3a_12_3(t,ro,comp,dcomp,jac,deriv,fait,
     1	epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)
 
c	cycles PP, CNO et 3 alpha
c	cf. Clayton p. 380, 392 et 430,
 
c	on utilise, suivant que table=.true./.false. (initialise ci avant):
c	soit les tables de Caughlan et Fowler 1988
c	soit les tables crees par tab_reac
 
c	les abondances initiales sont deduites de celles de H et de He4
c	suivant les rapports d'abondance d'Anders & Grevesse
c	Geochimica et Cosmochimica Acta, 53, 197, 1989
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
 
c	CESAM, Version 3
 
c entree :
c	t : temperature cgs
c	ro : densite cgs
c	comp : abondances
c	deriv=.true. : on calcule le jacobien
c	fait=1 : initialisation de la composition chimique
c	    =2 : calcul de dcomp et jacobien si deriv
c	    =3 : energie nucleaire et derivees / t et ro
c	    =4 : production de neutrinos
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
 
c	r(1) : reaction H(H,e+ nu)H2			PP
c	r(2) : reaction H2(H,g)H3
c	r(3) : reaction He3(He3,2H)He4
c	r(4) : reaction He4(He3,g)Be7
c	r(5) : reaction Li7(H,He4)He4
c	r(6) : reaction Be7(e-,nu g)Li7
c	r(7) : reaction Be7(H,g)B8(,e+ nu)Be8(,He4)He4
 
c	r(8) : reaction C12(H,g)N13(,e+ nu)C13	CNO
c	r(9) : reaction C13(H,g)N14
c	r(10) : reaction N14(H,g)O15(e+,nu)N15
c	r(11) : reaction N15(H,g)O16
c	r(12) : reaction N15(H,He4)C12
c	r(13) : reaction O16(H,g)F17(,e+ nu)O17
c	r(14) : reaction O17(H,He4)N14
 
c	r(15) : reaction He4(2He4,g)C12		3 alpha
c	r(16) : reaction C12(He4,g)N13
c	r(17) : reaction O16(He4,g)Ne20
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'modele_3.common'
	include 'evol_chim_3.common'
	include 'ctephy.common'
 
	integer fait,i
 
	real*8 t,ro,dcomp(1),jac(1),comp(1),dnelx,zs2,
     1	r(pnreac),rt(pnreac),rro(pnreac),rx(pnreac),
     2	q(pnreac),qt(pnreac),qro(pnreac),qx(pnreac),
     3	nel,epsilon(5),et,ero,ex(1),abon_meteor(pnchim),
     4	hhe,be7e,b8e,n13e,o15e,f17e,rxx(pnreac)
	
c	abondances des nuclides d'Anders & Grevesse,
c	Geochimica et Cosmochimica Acta, 53, 197, 1989
c	donnees de la table 1 page 198, mutipliees par les rapports
c	isotopiques de la table 3 p. 200
c	pour Be7: valeur arbitraire non nulle : 0.73d-5
	
	logical deriv,table
 
	data abon_meteor/2.79d10,9.49d5,3.86d5,2.72d9,52.82d0,0.73d-5,
     1	 9.99d6,1.11d5,3.12d6,1.15d4,2.37d7,9.04d3,9*0.d0/	
 
c	logical deriv,table
 
	save table,zs2,dnelx
 
2000	format((1x,1p8d10.3))
 
c	pause'mettre a jour ppcno3a_12_3'
 
	goto(100,200,300,400,200),fait
 
c	initialisations
 
100	nchim=12			!nombre total d'elements
	nreac=17			!nombre de reactions utilisees
 
	table=.true.	!on utilise les tables: reac_gf_3
c	table=.false.	!on utilise la tabulation des formules: reac_t_3
 
	if(table)then		 !fictif
	 call reac_gf_3(.7d0,1.d7,1.d0,r,rt,rro,rx,rxx,q,qt,qro,qx,.false.)
	else
	 call reac_t_3(.7d0,1.d7,1.d0,r,rt,rro,rx,rxx,q,qt,qro,qx,.false.)
	endif
 
	nucleo(1)=ah		!initialisation de nucleo: masses des noyaux
	nucleo(2)=ah2
	nucleo(3)=ahe3
	nucleo(4)=ahe4
	nucleo(5)=ali7
	nucleo(6)=abe7
	nucleo(7)=ac12
	nucleo(8)=ac13
	nucleo(9)=an14
	nucleo(10)=an15
	nucleo(11)=ao16
	nucleo(12)=ao17
	
	ihe4=4	!indice de He4
 
	nom_elem(1)=' H '
	nom_elem(2)=' H2'
	nom_elem(3)='He3'
	nom_elem(4)='He4'
	nom_elem(5)='Li7'
	nom_elem(6)='Be7'
	nom_elem(7)='C12'
	nom_elem(8)='C13'
	nom_elem(9)='N14'
	nom_elem(10)='N15'
	nom_elem(11)='O16'
	nom_elem(12)='O17'
	
c	abondances initiales en masse a partir de X: comp(1), et Y=comp(4)
 
	comp(1)=x0/nucleo(1)
	comp(2)=x0/nucleo(1)*abon_meteor(2)/abon_meteor(1)
	comp(3)=y0/nucleo(4)*abon_meteor(3)/abon_meteor(4)
	comp(4)=y0/nucleo(4)		
	do i=5,nchim
	 comp(i)=x0/nucleo(1)*abon_meteor(i)/abon_meteor(1)
	enddo
	do i=1,nchim
	 ab_ini(i)=comp(i)*nucleo(i)
	enddo
	
c	la somme des elements lourds est Z0=1-X0-Y0	
	
	nel=0.d0		!nel: vt
	do i=ihe4+1,nchim		!elements lourds
	 nel=nel+ab_ini(i)
	enddo
c	write(6,2000)1.d0-x0-y0,nel
	nel=(1.d0-x0-y0)/nel
	write(6,21)nel
	write(2,21)nel		
21	format(1x,'multiplication abondances elements lourds par',1pd10.3)	
	do i=ihe4+1,nchim
	 ab_ini(i)=ab_ini(i)*nel
	 comp(i)=ab_ini(i)/nucleo(i)
	enddo
	
c	nel=0.d0		!verif
c	do i=ihe4+1,nchim
c	 nel=nel+ab_ini(i)
c	enddo
c	write(6,2000)(comp(i),i=1,nchim)
c	write(6,2000)(ab_ini(i),i=1,nchim)
c	write(6,2000)nel
c	pause
 
c	abondances minimales
 
	ab_min(1)=1.d-3		!H
	ab_min(2)=1.d-20	!D
	ab_min(3)=5.d-7		!He3
	ab_min(4)=1.d-3		!He4
	ab_min(5)=1.d-14	!Li7
	ab_min(6)=1.d-17	!Be7
	ab_min(7)=5.d-6		!C12
	ab_min(8)=1.d-7		!C13
	ab_min(9)=1.d-6		!N14
	ab_min(10)=5.d-9	!N15
	ab_min(11)=1.d-5	!O16
	ab_min(12)=5.d-9	!O17
 
	write(2,*)' '
	write(2,*)'Reactions thermonucleaires des cycles PP, CNO, 3 Alpha'
	write(2,*)'par interpolation des tables de Caughlan et Fowler 1988'
	write(2,*)' '
	write(2,*)'nombre de reactions : ',nreac
	write(2,*)'nombre d''elements chimiques : ',nchim
	write(2,*)' '
	write(2,20)
20	format(1x,'les abondances initiales sont deduites de X suivant',/,
     1	1x,'Anders & Grevesse, Geochimica et Cosmochimica Acta,',/,
     2	1x,'53, 197, 1989. Donnees de la table 1 page 198,',/,
     3	1x,'mutipliees par les rapports isotopiques de la',/,
     4	1x,'table 3 p. 200; pour Be7 on utilise une valeur arbitraire',/,
     5	1x,'non nulle : 0.73d-5')
	write(2,*)		
	write(2,1)(comp(i)*nucleo(i),i=1,nchim)
1	format(3x,' H : ',1pd10.3,3x,' H2 : ',1pd10.3,3x,' He3 : ',
     1	1pd10.3,3x,' He4 : ',1pd10.3,3x,' Li7 : ',1pd10.3,
     2	3x,' Be7 : ',1pd10.3,/,1x,' C12 : ',1pd10.3,3x,' C13 : ',
     3	1pd10.3,3x,' N14 : ',1pd10.3,3x,' N15 : ',1pd10.3,
     4	3x,' O16 : ',1pd10.3,3x,' O17 : ',1pd10.3)
	write(2,*)' '
	write(2,*)'abondances negligeables:'
	write(2,1)(ab_min(i),i=1,nchim)
	write(2,*)' '
	write(2,*)'aucun element a l''equilibre'
	write(2,*)' '
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
	write(6,*)'aucun element a l''equilibre'
	write(6,*)' '
	write(6,*)'on utilise une table'
	write(6,*)' '
	write(6,*)'pour l''evolution temporelle, test de precision sur H et He4'
	write(6,*)' '
	write(6,*)' '
 
	do i=1,nchim
	 ab_min(i)=ab_min(i)/nucleo(i)
	enddo
 
	zs2=8./16.*z0		!on ne suppose pas Z diffuse
	dnelx=1.-2.*ah/ahe4	!d nel / dX en supposant Y=1-X-Z
 
	return
 
c	reactions
 
200	if(t .lt. t_inf)then	!si t<t_inf
	 do i=1,nbelem
	  dcomp(i)=0.
	  ex(i)=0
	 enddo
	 do i=1,nbelem*nbelem
	  jac(i)=0.
	 enddo
	
	 do i=1,4
	  epsilon(i)=0
	 enddo
	 et=0
	 ero=0
	 return
	endif
	if(table)then
	 call reac_gf_3(comp(1),t,ro,r,rt,rro,rx,rxx,q,qt,qro,qx,deriv)
	else
	 call reac_t_3(comp(1),t,ro,r,rt,rro,rx,rxx,q,qt,qro,qx,deriv)
	endif
 
c	nombre d'electrons / mole/g Y est compte comme 1-X-Z0
 
	nel=comp(1)+2.*(1.-comp(1)*ah-z0)/ahe4+zs2
 
c	write(6,*)'comp'
c	write(6,2000)(comp(i),i=1,nchim)
c	write(6,*)'reactions'
c	write(6,2000)(r(i),i=1,nreac)
 
c	equations d'evolution
 
	dcomp(1)=-(2.*r(1)*comp(1)+r(2)*comp(2)+r(5)*comp(5)
     1	+r(7)*comp(6)+r(8)*comp(7)+r(9)*comp(8)+r(10)*comp(9)
     2	+(r(11)+r(12))*comp(10)+r(13)*comp(11)
     3	+r(14)*comp(12))*comp(1)+2.*r(3)*comp(3)**2		!H
	dcomp(2)=(r(1)*comp(1)-r(2)*comp(2))*comp(1)			!H2
	dcomp(3)=r(2)*comp(1)*comp(2)-(2.*r(3)*comp(3)
     1	+r(4)*comp(4))*comp(3)					!He3
	dcomp(4)=(r(3)*comp(3)-r(4)*comp(4))*comp(3)
     1	+(2.*(r(5)*comp(5)+r(7)*comp(6))+r(12)*comp(10)
     2	+r(14)*comp(12))*comp(1)-(3.*r(15)*comp(4)**2
     3	+r(16)*comp(7)+r(17)*comp(11))*comp(4)			!He4
	dcomp(5)=-r(5)*comp(1)*comp(5)+r(6)*comp(6)*nel			!Li7
	dcomp(6)=r(4)*comp(3)*comp(4)-(r(6)*nel+r(7)*comp(1))*comp(6)	!Be7
	dcomp(7)=(-r(8)*comp(7)+r(12)*comp(10))*comp(1)
     1	+(r(15)*comp(4)**2-r(16)*comp(7))*comp(4)		!C12
	dcomp(8)=(r(8)*comp(7)-r(9)*comp(8))*comp(1)			!C13
	dcomp(9)=(r(9)*comp(8)-r(10)*comp(9)+r(14)*comp(12))*comp(1)	!N14
	dcomp(10)=(r(10)*comp(9)-(r(11)+r(12))*comp(10))*comp(1)	!N15
	dcomp(11)=(r(11)*comp(10)-r(13)*comp(11))*comp(1)
     1	+(r(16)*comp(7)-r(17)*comp(11))*comp(4)			!O16
	dcomp(12)=(r(13)*comp(11)-r(14)*comp(12))*comp(1)		!O17
 
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
 
c	calcul du jacobien
 
	if(deriv .or. fait .eq. 6)then	!jacp(i,j) : equation, j : element i
 
c	a travers l'effet d'ecran, les reactions dependent de H par rx
 
	 do i=1,nbelem*nbelem
	  jac(i)=0.
	 enddo		!i
	 jac(nbelem*(1-1)+1)=-r(2)*comp(2)-r(5)*comp(5)-r(7)*comp(6)
     1	-r(8)*comp(7)-r(9)*comp(8)-r(10)*comp(9)
     2	-(r(11)+r(12))*comp(10)-r(13)*comp(11)-r(14)*comp(12)
     3	-(4.*r(1)+2.*rx(1)*comp(1)+rx(2)*comp(2)+rx(5)*comp(5)
     4	+rx(7)*comp(6)+rx(8)*comp(7)+rx(9)*comp(8)
     5	+rx(10)*comp(9)+(rx(11)+rx(12))*comp(10)
     6	+rx(13)*comp(11)+rx(14)*comp(12))*comp(1)
     7	+2.*rx(3)*comp(3)**2					!d /H
	 jac(nbelem*(2-1)+1)=-r(2)*comp(1)				!d/dH2
	 jac(nbelem*(3-1)+1)=4.*r(3)*comp(3)				!d /He3
	 jac(nbelem*(5-1)+1)=-r(5)*comp(1)				!d /Li7
	 jac(nbelem*(6-1)+1)=-r(7)*comp(1)				!d /Be7
	 jac(nbelem*(7-1)+1)=-r(8)*comp(1)				!d /C12
	 jac(nbelem*(8-1)+1)=-r(9)*comp(1)				!d /C13
	 jac(nbelem*(9-1)+1)=-r(10)*comp(1)				!d /N14
	 jac(nbelem*(10-1)+1)=-(r(11)+r(12))*comp(1)			!d /N15
	 jac(nbelem*(11-1)+1)=-r(13)*comp(1)				!d /O16
	 jac(nbelem*(12-1)+1)=-r(14)*comp(1)				!d /O17
 
	 jac(nbelem*(1-1)+2)=-r(2)*comp(2)+(2.*r(1)+rx(1)*comp(1)
     1	-rx(2)*comp(2))*comp(1)					!d /H
	 jac(nbelem*(2-1)+2)=-r(2)*comp(1)				!d /H2
 
	 jac(nbelem*(1-1)+3)=(r(2)+rx(2)*comp(1))*comp(2)
     1	-(2.*rx(3)*comp(3)+rx(4)*comp(4))*comp(3)		!d /H
	 jac(nbelem*(2-1)+3)=r(2)*comp(1)				!d /H2
	 jac(nbelem*(3-1)+3)=-4.*r(3)*comp(3)-r(4)*comp(4)		!d /He3
	 jac(nbelem*(4-1)+3)=-r(4)*comp(3)				!d /He4
 
	 jac(nbelem*(1-1)+4)=(rx(3)*comp(3)-rx(4)*comp(4))*comp(3)
     1	+2.*((r(5)+rx(5)*comp(1))*comp(5)
     2	+(r(7)+rx(7)*comp(1))*comp(6))
     3	+r(12)*comp(10)+r(14)*comp(12)
     4	+(rx(12)*comp(10)+rx(14)*comp(12))*comp(1)
     5	-(3.*rx(15)*comp(4)**2+rx(16)*comp(7)
     6	+rx(17)*comp(11))*comp(4)				!d /H
	 jac(nbelem*(3-1)+4)=2.*r(3)*comp(3)-r(4)*comp(4)		!d /He3
	 jac(nbelem*(4-1)+4)=-r(4)*comp(3)-9.*r(15)*comp(4)**2
     1	-r(16)*comp(7)-r(17)*comp(11)				!d /He4
	 jac(nbelem*(5-1)+4)=r(5)*comp(1)*2.				!d /Li7
	 jac(nbelem*(6-1)+4)=r(7)*comp(1)*2.				!d /Be7
	 jac(nbelem*(7-1)+4)=-r(16)*comp(4)				!d /C12
	 jac(nbelem*(10-1)+4)=r(12)*comp(1)				!d /N15
	 jac(nbelem*(11-1)+4)=-r(17)*comp(4)				!d /O16
	 jac(nbelem*(12-1)+4)=r(14)*comp(1)				!d /O17
 
	 jac(nbelem*(1-1)+5)=-(r(5)+rx(5)*comp(1))*comp(5)
     1	+comp(6)*(rx(6)*nel+r(6)*dnelx)				!d /H
	 jac(nbelem*(5-1)+5)=-r(5)*comp(1)				!d /Li7
	 jac(nbelem*(6-1)+5)=r(6)*nel					!d /Be7
 
	 jac(nbelem*(1-1)+6)=rx(4)*comp(3)*comp(4)-(rx(6)*nel+r(6)*dnelx
     1	+r(7)+rx(7)*comp(1))*comp(6)				!d /H
	 jac(nbelem*(3-1)+6)=r(4)*comp(4)				!d /He3
	 jac(nbelem*(4-1)+6)=r(4)*comp(3)				!d /He4
	 jac(nbelem*(6-1)+6)=-r(6)*nel-r(7)*comp(1)			!d /Be7
 
	 jac(nbelem*(1-1)+7)=-r(8)*comp(7)+r(12)*comp(10)
     1	-(rx(8)*comp(7)-rx(12)*comp(10))*comp(1)
     2	+(rx(15)*comp(4)**2-rx(16)*comp(7))*comp(4)		!d /H
	 jac(nbelem*(4-1)+7)=3.*r(15)*comp(4)**2-r(16)*comp(7)		!d /He4
	 jac(nbelem*(7-1)+7)=-r(8)*comp(1)-r(16)*comp(4)		!d /C12
	 jac(nbelem*(10-1)+7)=r(12)*comp(1)				!d /N15
 
	 jac(nbelem*(1-1)+8)=r(8)*comp(7)-r(9)*comp(8)
     1	+(rx(8)*comp(7)-rx(9)*comp(8))*comp(1)			!d /H
	 jac(nbelem*(7-1)+8)=r(8)*comp(1)				!d /C12
	 jac(nbelem*(8-1)+8)=-r(9)*comp(1)				!d /C13
 
	 jac(nbelem*(1-1)+9)=r(9)*comp(8)-r(10)*comp(9)+r(14)*comp(12)
     1	+(rx(9)*comp(8)-rx(10)*comp(9)+rx(14)*comp(12))*comp(1)	!d /H
	 jac(nbelem*(8-1)+9)=r(9)*comp(1)				!d /C13
	 jac(nbelem*(9-1)+9)=-r(10)*comp(1)				!d /N14
	 jac(nbelem*(12-1)+9)=r(14)*comp(1)				!d /O17
 
	 jac(nbelem*(1-1)+10)=r(10)*comp(9)-(r(11)+r(12))*comp(10)
     1	+(rx(10)*comp(9)-(rx(11)+rx(12))*comp(10))*comp(1)	!d /H
	 jac(nbelem*(9-1)+10)=r(10)*comp(1)				!d /N14
	 jac(nbelem*(10-1)+10)=-(r(11)+r(12))*comp(1)			!d /N15
 
	 jac(nbelem*(1-1)+11)=r(11)*comp(10)-r(13)*comp(11)+(rx(11)*comp(10)
     1	-rx(13)*comp(11))*comp(1)+(rx(16)*comp(7)
     2	-rx(17)*comp(11))*comp(4)				!d /H
	 jac(nbelem*(4-1)+11)=r(16)*comp(7)-r(17)*comp(11)		!d /He4
	 jac(nbelem*(7-1)+11)=r(16)*comp(4)				!d /C12
	 jac(nbelem*(10-1)+11)=r(11)*comp(1)				!d /N15
	 jac(nbelem*(11-1)+11)=-r(13)*comp(1)-r(17)*comp(4)		!d /O16
 
	 jac(nbelem*(1-1)+12)=r(13)*comp(11)-r(14)*comp(12)
     1	+(rx(13)*comp(11)-rx(14)*comp(12))*comp(1)		!d /H
	 jac(nbelem*(11-1)+12)=r(13)*comp(1)				!d /O16
	 jac(nbelem*(12-1)+12)=-r(14)*comp(1)				!d /O17
	
c	 pour le moment angulaire ou Z jac=0
 
c	 unites de temps pour integration temporelle
 
	 do i=1,nbelem*nbelem
	  jac(i)=jac(i)*secon6
	 enddo
 
	endif
 
	do i=1,nbelem
	 dcomp(i)=dcomp(i)*secon6
	enddo
 
	if(fait .ne. 5)return
 
c	calcul de la production d'energie nucleaire et derivees
c	pour H2(H,g)He3, q(2)H**2=q(2)*r(1)/r(2)
 
300	if(t .ge. t_inf)then
	 if(fait .ne. 5)then
	  if(table)then
	   call reac_gf_3(comp(1),t,ro,r,rt,rro,rx,rxx,q,qt,qro,qx,deriv)
	  else
	   call reac_t_3(comp(1),t,ro,r,rt,rro,rx,rxx,q,qt,qro,qx,deriv)
	  endif
	endif
 
c	 nombre d'electrons / mole/g Y est compte comme 1-X-Z
 
	 nel=comp(1)+2.*(1.-comp(1)*ah-z0)/ahe4+zs2
 
	 epsilon(1)=0.
	 epsilon(2)=(q(1)*comp(1)+q(2)*comp(2)+q(5)*comp(5)+q(7)*comp(6))
     1	*comp(1)+(q(3)*comp(3)+q(4)*comp(4))*comp(3)+q(6)*nel*comp(6)
	 epsilon(3)=(q(8)*comp(7)+q(9)*comp(8)+q(10)*comp(9)+
     1	(q(11)+q(12))*comp(10)+q(13)*comp(11)+q(14)*comp(12))*comp(1)
	 epsilon(4)=(q(15)*comp(4)**2+q(16)*comp(7)+q(17)*comp(11))*comp(4)
	 do i=2,4
	  epsilon(1)=epsilon(1)+epsilon(i)
	 enddo
 
c	 write(6,2000)(epsilon(i),i=1,4)
c	 write(6,2000)(comp(i),i=1,nchim)
c	 write(6,2000)(q(i),i=1,nreac)
c	 write(6,2000)nel
c	 write(6,2000)q(1)*comp(1)*comp(1)
c	 write(6,2000)q(2)*comp(2)*comp(1)
c	 write(6,2000)q(5)*comp(5)*comp(1)
c	 write(6,2000)q(7)*comp(6)*comp(1)
c	 write(6,2000)q(3)*comp(3)*comp(3)
c	 write(6,2000)q(4)*comp(4)*comp(3)
c	 write(6,2000)q(6)*nel*comp(6)
c	 write(6,2000)q(8)*comp(7)*comp(1)
c	 write(6,2000)q(9)*comp(8)*comp(1)
c	 write(6,2000)q(10)*comp(9)*comp(1)
c	 write(6,2000)q(11)*comp(10)*comp(1)
c	 write(6,2000)q(12)*comp(10)*comp(1)
c	 write(6,2000)q(13)*comp(11)*comp(1)
c	 write(6,2000)q(14)*comp(12)*comp(1)
c	 write(6,2000)q(15)*comp(4)**2*comp(4)
c	 write(6,2000)q(16)*comp(7)*comp(4)
c	 write(6,2000)q(17)*comp(11)*comp(4)
c	 pause
 
	 if(deriv)then
	  et=(qt(1)*comp(1)+qt(2)*comp(2)+qt(5)*comp(5)+qt(7)*comp(6))
     1	*comp(1)+(qt(3)*comp(3)+qt(4)*comp(4))*comp(3)
     2	+qt(6)*nel*comp(6)+(qt(8)*comp(7)+qt(9)*comp(8)
     3	+qt(10)*comp(9)+(qt(11)+qt(12))*comp(10)+qt(13)*comp(11)
     4	+qt(14)*comp(12))*comp(1)+(qt(15)*comp(4)**2
     5	+qt(16)*comp(7)+qt(17)*comp(11))*comp(4)
	  ero=(qro(1)*comp(1)+qro(2)*comp(2)+qro(5)*comp(5)+qro(7)*comp(6))
     1	*comp(1)+(qro(3)*comp(3)+qro(4)*comp(4))*comp(3)
     2	+qro(6)*nel*comp(6)+(qro(8)*comp(7)+qro(9)*comp(8)
     3	+qro(10)*comp(9)+(qro(11)+qro(12))*comp(10)+qro(13)*comp(11)
     4	+qro(14)*comp(12))*comp(1)+(qro(15)*comp(4)**2
     5	+qro(16)*comp(7)+qro(17)*comp(11))*comp(4)
	  ex(1)=(qx(1)*comp(1)+qx(2)*comp(2)+qx(5)*comp(5)+qx(7)*comp(6))
     1	*comp(1)+(qx(3)*comp(3)+qx(4)*comp(4))*comp(3)
     2	+(qx(6)*nel+q(6)*dnelx)*comp(6)+(qx(8)*comp(7)+qx(9)*comp(8)
     3	+qx(10)*comp(9)+(qx(11)+qx(12))*comp(10)+qx(13)*comp(11)
     4	+qx(14)*comp(12))*comp(1)+(qx(15)*comp(4)**2
     5	+qx(16)*comp(7)+qx(17)*comp(11))*comp(4)
     6	+2.*q(1)*comp(1)+q(2)*comp(2)+q(5)*comp(5)+q(7)*comp(6)
     7	+q(8)*comp(7)+q(9)*comp(8)+q(10)*comp(9)
     8	+(q(11)+q(12))*comp(10)+q(13)*comp(11)+q(14)*comp(12)
	  ex(2)=q(2)*comp(1)
	  ex(3)=2.*q(3)*comp(3)+q(4)*comp(4)
	  ex(4)=q(4)*comp(3)+3.*q(15)*comp(4)**2+q(16)*comp(7)+q(17)*comp(11)
	  ex(5)=q(5)*comp(1)
	  ex(6)=q(6)*nel+q(7)*comp(1)
	  ex(7)=q(8)*comp(1)+q(16)*comp(4)
	  ex(8)=q(9)*comp(1)
	  ex(9)=q(10)*comp(1)
	  ex(10)=(q(11)+q(12))*comp(1)
	  ex(11)=q(13)*comp(1)+q(17)*comp(4)
	  ex(12)=q(14)*comp(1)
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
 
c	production de neutrinos
 
400	if(t .ge. t_inf)then
	 if(table)then
	  call reac_gf_3(comp(1),t,ro,r,rt,rro,rx,rxx,q,qt,qro,qx,.false.)
	 else
	  call reac_t_3(comp(1),t,ro,r,rt,rro,rx,rxx,q,qt,qro,qx,.false.)
	 endif
 
	 hhe=r(1)*comp(1)**2/amu
	 be7e=r(6)*nel*comp(6)/amu
	 b8e=r(7)*comp(1)*comp(6)/amu
	 n13e=r(8)*comp(1)*comp(7)/amu
	 o15e=r(10)*comp(1)*comp(9)/amu
	 f17e=r(13)*comp(1)*comp(10)/amu
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
