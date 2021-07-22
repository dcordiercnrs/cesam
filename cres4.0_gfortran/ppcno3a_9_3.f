c**************************************************************************
 
	subroutine ppcno3a_9_3(t,ro,comp,dcomp,jac,deriv,fait,
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
	include 'ctephy.common'
	include 'modele_3.common'
	include 'evol_chim_3.common'
 
	integer fait,i
 
	real*8 t,ro,dcomp(1),jac(1),comp(1),dbe7y,dli7y,dbe7he3,dli7he3,
     1	r(pnreac),rt(pnreac),rro(pnreac),rx(pnreac),ddenx,
     2	q(pnreac),qt(pnreac),qro(pnreac),qx(pnreac),
     3	nel,epsilon(5),et,ero,ex(1),den,abon_meteor(pnchim),
     4	hhe,be7e,b8e,n13e,o15e,f17e,rxx(pnreac),dnelx,
     5	h2,dh2t,dh2ro,dh2x,be7,dbe7t,dbe7ro,dbe7x,
     6	li7,dli7t,dli7ro,dli7x,zs2
 
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
 
c	pause'modifier ppcno3a_9_3'
 
	goto(100,200,300,400,200),fait
 
c	initialisations
 
100	nchim=9			!nombre d'elements chimiques
	nreac=17		!nombre de reactions utilisees
 
	table=.true.	!on utilise les tables: reac_gf3
c	table=.false.	!on utilise la tabulation des formules: reac_t_3
 
	if(table)then		 !fictif
	 call reac_gf_3(.7d0,1.d7,1.d0,r,rt,rro,rx,rxx,q,qt,qro,qx,.false.)
	else
	 call reac_t_3(.7d0,1.d7,1.d0,r,rt,rro,rx,rxx,q,qt,qro,qx,.false.)
	endif
 
	nucleo(1)=ah		!initialisation de nucleo: masses des noyaux
	nucleo(2)=ahe3
	nucleo(3)=ahe4
	nucleo(4)=ac12
	nucleo(5)=ac13
	nucleo(6)=an14
	nucleo(7)=an15
	nucleo(8)=ao16
	nucleo(9)=ao17
	
	ihe4=3	!indice de He4
 
	nom_elem(1)=' H '
	nom_elem(2)='He3'
	nom_elem(3)='He4'
	nom_elem(4)='C12'
	nom_elem(5)='C13'
	nom_elem(6)='N14'
	nom_elem(7)='N15'
	nom_elem(8)='O16'
	nom_elem(9)='O17'
	
c	abondances initiales en masse a partir de X: comp(1), et Y=comp(3)
 
	comp(1)=x0/nucleo(1)		!H
	comp(2)=y0/nucleo(3)*abon_meteor(3)/abon_meteor(4)	!He3
	comp(3)=y0/nucleo(3)		
	do i=4,nchim	!Li7 et Be7 a l'equilibre
	 comp(i)=x0/nucleo(1)*abon_meteor(i+3)/abon_meteor(1)
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
 
	ab_min(1)=1.d-3		!abondances negligeables
	ab_min(2)=5.d-7
	ab_min(3)=1.d-3
	ab_min(4)=5.d-6
	ab_min(5)=1.d-7
	ab_min(6)=1.d-6
	ab_min(7)=5.d-9
	ab_min(8)=1.d-5
	ab_min(9)=5.d-9
 
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
     4	1x,'table 3 p. 200.')
	write(2,*)		
	write(2,1)(comp(i)*nucleo(i),i=1,9)
1	format(3x,' H : ',1pd10.3,3x,' He3 : ',1pd10.3,3x,' He4 : ',1pd10.3,3x,
     1	' C12 : ',1pd10.3,3x,' C13 : ',1pd10.3,3x,
     2	/,1x,' N14 : ',1pd10.3,3x,
     3	' N15 : ',1pd10.3,3x,' O16 : ',1pd10.3,3x,' O17 : ',1pd10.3)
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
	write(6,1)(comp(i)*nucleo(i),i=1,9)
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
	if(table)then
	 call reac_gf_3(comp(1),t,ro,r,rt,rro,rx,rxx,q,qt,qro,qx,deriv)
	else
	 call reac_t_3(comp(1),t,ro,r,rt,rro,rx,rxx,q,qt,qro,qx,deriv)
	endif
 
c	write(6,*)'comp'
c	write(6,2000)(comp(i),i=1,nchim)
c	write(6,*)'reactions'
c	write(6,2000)(r(i),i=1,nreac)
 
c	equations d'evolution
 
	dcomp(1)=-(3.*r(1)*comp(1)+r(8)*comp(4)+r(9)*comp(5)+r(10)*comp(6)
     1	+(r(11)+r(12))*comp(7)+r(13)*comp(8)+r(14)*comp(9))*comp(1)
     2	+(2.*r(3)*comp(2)-r(4)*comp(3))*comp(2)			!H
	dcomp(2)=r(1)*comp(1)**2-(2.*r(3)*comp(2)+r(4)*comp(3))*comp(2)	!He3
	dcomp(3)=(r(3)*comp(2)+r(4)*comp(3))*comp(2)
     1	+(r(12)*comp(7)+r(14)*comp(9))*comp(1)-(3.*r(15)*comp(3)**2
     3	+r(16)*comp(4)+r(17)*comp(8))*comp(3)			!He4
	dcomp(4)=(-r(8)*comp(4)+r(12)*comp(7))*comp(1)
     1	+(r(15)*comp(3)**2-r(16)*comp(4))*comp(3)		!C12
	dcomp(5)=(r(8)*comp(4)-r(9)*comp(5))*comp(1)			!C13
	dcomp(6)=(r(9)*comp(5)-r(10)*comp(6)+r(14)*comp(9))*comp(1)	!N14
	dcomp(7)=(r(10)*comp(6)-(r(11)+r(12))*comp(7))*comp(1)		!N15
	dcomp(8)=(r(11)*comp(7)-r(13)*comp(8))*comp(1)
     1	+(r(16)*comp(4)-r(17)*comp(8))*comp(3)			!O16
	dcomp(9)=(r(13)*comp(8)-r(14)*comp(9))*comp(1)			!O17
	
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
	 jac(nbelem*(1-1)+1)=-r(8)*comp(4)-r(9)*comp(5)-r(10)*comp(6)
     2	-(r(11)+r(12))*comp(7)-r(13)*comp(8)-r(14)*comp(9)
     3	-(6.*r(1)+3.*rx(1)*comp(1)+rx(8)*comp(4)+rx(9)*comp(5)
     5	+rx(10)*comp(6)+(rx(11)+rx(12))*comp(7)
     6	+rx(13)*comp(8)+rx(14)*comp(9))*comp(1)
     7	+(2.*rx(3)*comp(2)-rx(4)*comp(3))*comp(2)		!d /H
	 jac(nbelem*(2-1)+1)=4.*r(3)*comp(2)-r(4)*comp(3)		!d /He3
	 jac(nbelem*(3-1)+1)=-r(4)*comp(2)				!d /He4
	 jac(nbelem*(4-1)+1)=-r(8)*comp(1)				!d /C12
	 jac(nbelem*(5-1)+1)=-r(9)*comp(1)				!d /C13
	 jac(nbelem*(6-1)+1)=-r(10)*comp(1)				!d /N14
	 jac(nbelem*(7-1)+1)=-(r(11)+r(12))*comp(1)			!d /N15
	 jac(nbelem*(8-1)+1)=-r(13)*comp(1)				!d /O16
	 jac(nbelem*(9-1)+1)=-r(14)*comp(1)				!d /O17
 
	 jac(nbelem*(1-1)+2)=(2.*r(1)+rx(1)*comp(1))*comp(1)
     1	-(2.*rx(3)*comp(2)+rx(4)*comp(3))*comp(2)		!d /H
	 jac(nbelem*(2-1)+2)=-4.*r(3)*comp(2)-r(4)*comp(3)		!d /He3
	 jac(nbelem*(3-1)+2)=-r(4)*comp(2)				!d /He4
 
	 jac(nbelem*(1-1)+3)=(rx(3)*comp(2)+rx(4)*comp(3))*comp(2)
     1	+r(12)*comp(7)+r(14)*comp(9)
     2	+(rx(12)*comp(7)+rx(14)*comp(9))*comp(1)
     3	-(3.*rx(15)*comp(3)**2+rx(16)*comp(4)
     4	+rx(17)*comp(8))*comp(3)				!d /H
	 jac(nbelem*(2-1)+3)=2.*r(3)*comp(2)+r(4)*comp(3)		!d /He3
	 jac(nbelem*(3-1)+3)=r(4)*comp(2)-9.*r(15)*comp(3)**2
     1	-r(16)*comp(4)-r(17)*comp(8)				!d /He4
	 jac(nbelem*(4-1)+3)=-r(16)*comp(3)				!d /C12
	 jac(nbelem*(7-1)+3)=r(12)*comp(1)				!d /N15
	 jac(nbelem*(8-1)+3)=-r(17)*comp(3)				!d /O16
	 jac(nbelem*(9-1)+3)=r(14)*comp(1)				!d /O17
 
	 jac(nbelem*(1-1)+4)=-r(8)*comp(4)+r(12)*comp(7)
     1	-(rx(8)*comp(4)-rx(12)*comp(7))*comp(1)
     2	+(rx(15)*comp(3)**2-rx(16)*comp(4))*comp(3)		!d /H
	 jac(nbelem*(3-1)+4)=3.*r(15)*comp(3)**2-r(16)*comp(4)		!d /He4
	 jac(nbelem*(4-1)+4)=-r(8)*comp(1)-r(16)*comp(3)		!d /C12
	 jac(nbelem*(7-1)+4)=r(12)*comp(1)				!d /N15
 
	 jac(nbelem*(1-1)+5)=r(8)*comp(4)-r(9)*comp(5)
     1	+(rx(8)*comp(4)-rx(9)*comp(5))*comp(1)			!d /H
	 jac(nbelem*(4-1)+5)=r(8)*comp(1)				!d /C12
	 jac(nbelem*(5-1)+5)=-r(9)*comp(1)				!d /C13
 
	 jac(nbelem*(1-1)+6)=r(9)*comp(5)-r(10)*comp(6)+r(14)*comp(9)
     1	+(rx(9)*comp(5)-rx(10)*comp(6)+rx(14)*comp(9))*comp(1)	!d /H
	 jac(nbelem*(5-1)+6)=r(9)*comp(1)				!d /C13
	 jac(nbelem*(6-1)+6)=-r(10)*comp(1)				!d /N14
	 jac(nbelem*(9-1)+6)=r(14)*comp(1)				!d /O17
 
	 jac(nbelem*(1-1)+7)=r(10)*comp(6)-(r(11)+r(12))*comp(7)
     1	+(rx(10)*comp(6)-(rx(11)+rx(12))*comp(7))*comp(1)	!d /H
	 jac(nbelem*(6-1)+7)=r(10)*comp(1)				!d /N14
	 jac(nbelem*(7-1)+7)=-(r(11)+r(12))*comp(1)			!d /N15
 
	 jac(nbelem*(1-1)+8)=r(11)*comp(7)-r(13)*comp(8)+(rx(11)*comp(7)
     1	-rx(13)*comp(8))*comp(1)+(rx(16)*comp(4)
     2	-rx(17)*comp(8))*comp(3)				!d /H
	 jac(nbelem*(3-1)+8)=r(16)*comp(4)-r(17)*comp(8)		!d /He4
	 jac(nbelem*(4-1)+8)=r(16)*comp(3)				!d /C12
	 jac(nbelem*(7-1)+8)=r(11)*comp(1)				!d /N15
	 jac(nbelem*(8-1)+8)=-r(13)*comp(1)-r(17)*comp(3)		!d /O16
 
	 jac(nbelem*(1-1)+9)=r(13)*comp(8)-r(14)*comp(9)
     1	+(rx(13)*comp(8)-rx(14)*comp(9))*comp(1)		!d /H
	 jac(nbelem*(8-1)+9)=r(13)*comp(1)				!d /O16
	 jac(nbelem*(9-1)+9)=-r(14)*comp(1)				!d /O17
 
c	 pour le moment angulaire ou Z jac=0
 
 
c	 unites de temps pour integration temporelle
 
	 do i=1,nbelem*nbelem
	  jac(i)=jac(i)*secon6
	 enddo
	endif		!deriv
 
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
	 epsilon(3)=(q(8)*comp(4)+q(9)*comp(5)+q(10)*comp(6)+
     1	(q(11)+q(12))*comp(7)+q(13)*comp(8)+q(14)*comp(9))*comp(1)
	 epsilon(4)=(q(15)*comp(3)**2+q(16)*comp(4)+q(17)*comp(8))*comp(3)
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
     1	+qt(7)*be7+q(7)*dbe7t)*comp(1)
     2	+(qt(3)*comp(2)+qt(4)*comp(3))*comp(2)
     3	+nel*(qt(6)*be7+q(6)*dbe7t)+(qt(8)*comp(4)+qt(9)*comp(5)
     4	+qt(10)*comp(6)+(qt(11)+qt(12))*comp(7)+qt(13)*comp(8)
     5	+qt(14)*comp(9))*comp(1)
     1	+(qt(15)*comp(3)**2+qt(16)*comp(4)+qt(17)*comp(8))*comp(3)
 
	  ero=(qro(1)*comp(1)+qro(2)*h2+q(2)*dh2ro+qro(5)*li7+q(5)*dli7ro
     1	+qro(7)*be7+q(7)*dbe7ro)*comp(1)
     2	+(qro(3)*comp(2)+qro(4)*comp(3))*comp(2)
     3	+nel*(qro(6)*be7+q(6)*dbe7ro)+(qro(8)*comp(4)+qro(9)*comp(5)
     4	+qro(10)*comp(6)+(qro(11)+qro(12))*comp(7)+qro(13)*comp(8)
     5	+qro(14)*comp(9))*comp(1)
     1	+(qro(15)*comp(3)**2+qro(16)*comp(4)+qro(17)*comp(8))*comp(3)
 
	  ex(1)=(qx(1)*comp(1)+qx(2)*h2+q(2)*dh2x+qx(5)*li7+q(5)*dli7x
     1	+qx(7)*be7+q(7)*dbe7x)*comp(1)
     2	+2.*q(1)*comp(1)+q(2)*h2+q(5)*li7+q(7)*be7
     3	+(qx(3)*comp(2)+qx(4)*comp(3))*comp(2)+q(6)*dnelx*be7
     4	+nel*(qx(6)*be7+q(6)*dbe7x)+(qx(8)*comp(4)+qx(9)*comp(5)
     5	+qx(10)*comp(6)+(qx(11)+qx(12))*comp(7)+qx(13)*comp(8)
     6	+qx(14)*comp(9))*comp(1)
     7	+q(8)*comp(4)+q(9)*comp(5)+q(10)*comp(6)+(q(11)+q(12))*comp(7)
     8	+q(13)*comp(8)+q(14)*comp(9)
     9	+(qx(15)*comp(3)**2+qx(16)*comp(4)+qx(17)*comp(8))*comp(3)
	  ex(2)=(q(5)*dli7he3+q(7)*dbe7he3)*comp(1)+2.*q(3)*comp(2)
     1	+q(4)*comp(3)+q(6)*nel*dbe7he3
	  ex(3)=(q(5)*dli7y+q(7)*dbe7y)*comp(1)+q(4)*comp(2)+q(6)*nel*dbe7y
     1	+3.*q(15)*comp(3)**2+q(16)*comp(4)+q(17)*comp(8)
	  ex(4)=q(8)*comp(1)+q(16)*comp(3)
	  ex(5)=q(9)*comp(1)
	  ex(6)=q(10)*comp(1)
	  ex(7)=(q(11)+q(12))*comp(1)
	  ex(8)=q(13)*comp(1)+q(17)*comp(3)
	  ex(9)=q(14)*comp(1)
	
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
