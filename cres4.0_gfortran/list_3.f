c******************************************************************
 
	subroutine list_3(n,m,p,ro,t,r,l,compg,age,teff,nom,rstar,hp,mstar,
     1	methode,preci,ecritout,kap,epsilon,eint,xdot,ydot,rdot,
     2	nh1,nhe1,nhe2,alfa,beta,delta,cp,vaissala,conv,dlpdxt,dlpdxr,
     3	gradconv,gradrad,gradad,dcapdt,dcapdr,depsdr,depsdt,d2p,d2ro,
     4	anupp,anupep,anub8,anube7,anun13,anu015,degene,w,z,chaine,
     5	p_atm,t_atm,m_atm,tau,r_atm,ro_atm,k_atm,gradr_atm,grada_atm,
     6	gradc_atm,n_atm,i_pp,i_cno,i_3a,i_gr)
 
c	LISTING des resultats
 
c	P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
 
c	CESAM, Version 3
 
c	07 03 97 : ecriture de "Z_opa"
 
c	n : nombre de couches, n : le raccord a l'enveloppe, 1 : le centre
c	m(.) : masse/mtot
c	p(.) : pression
c	ro(.) : densite
c	t(.) : temperature
c	r(.) : rayon/rsoleil
c	l(.) : luminosite/lsoleil
c	compg(nbelem,.) : composition chimique par gramme
c	teff : temperature effective
c	mtot : masse totale
c	alpha : l/Hp
c	nbelem : nombre d'elements diffuses
c	nchim : nombre d'elements chimiques
c	nom : identification du modele
c	ecritout : .true. / .false. pour listing complet / partiel
c	p_atm : pression dans l'atmsphere
c	t_atm : temperature
c	m_atm : masse
c	tau : epaisseur optique
c	r_atm : rayon
c	n_atm : nombre de couches
c	epsilon : energie thermonucleaire et gravitationnelle
c	eint : energie interne specifique
c	degene: degenerescence
c	mstar: masse au temps t+dt, avec perte de masse
c	w(.): rotation
c	degene(.): degenerescence
c	z(.) : le Z
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'ctephy.common'
	include 'modele_3.common'
	include 'evol_chim_3.common'
 
	integer i,n,iln,j,n_atm,i_pp,i_cno,i_3a,i_gr,
     1	j2,j3,j4,j7,j7b,j12,j13,j14,j15,j16,j17
c	data j2,j3,j4,j7,j7b,j12,j13,j14,j15,j16,j17/11*0/
 
	real*8 m(1),r(1),p(1),ro(1),t(1),l(1),age,rstar,preci,mstar,
     1	teff,dlpdxr(1),dlpdxt(1),xdot(1),ydot(1),rdot(1),
     2	kap(1),epsilon(1),cte2,cte3,depsdr(1),depsdt(1),w(1),
     3	nh1(1),nhe1(1),nhe2(1),vaissala(1),alfa(1),beta(1),delta(1),
     4	gradad(1),gradrad(1),gradconv(1),d2p,d2ro,mu_mol,logg,
     5	hp(1),cp(1),dcapdt(1),dcapdr(1),eint(1),cv,dm,degene(1),
     6	anupp(1),anupep(1),anub8(1),anube7(1),anun13(1),anu015(1),
     7	compg(1),gamma1,p_atm(1),t_atm(1),m_atm(1),tau(1),r_atm(1),
     8	anufp,anufq,z(1),
     9	ro_atm(1),k_atm(1),gradr_atm(1),grada_atm(1),gradc_atm(1),
     1	anupps,anupeps,anub8s,anube7s,anun13s,anu015s,anufs,r_son(pn),
     2	v_son(pn),nu0,dnu02,dnu13,a,cte4,zsx
 
	logical conv(1),convp,init,ecritout,sautp,abon
c	data init/.false./
 
	character*30 chaine
	character*50 nom
	character*60 methode
 
	data j2,j3,j4,j7,j7b,j12,j13,j14,j15,j16,j17/11*0/
	data init/.false./
 
	save
 
2000	format((1p8d10.3))
 
	if(.not.init)then
	 init=.true.
	 sautp=.false.
	 cte2=4.*pi*(1.5d13)**2 !4 pi (R orb. terre)**2, pour neutrinos
	 cte3=mstar*msol/2.	!pour neutrinos
	 cte4=0.107d-6/cte2	!pour neutrino du B8 (G. Berthomieu)
 
c	 caracteristiques du modele
 
	 WRITE(2,4)
4	 FORMAT(//,1x,'---------------------------------------------',//,
     1	1x,'Constantes du modele',//)
	 WRITE(2,5)mtot,alpha,x0,1.d0-x0-z0,z0
5	 format(t2,'masse totale/Msol=',f6.3,t30,'l/Hp=',f4.2,
     1	t45,'Composition chimique initiale X=',f6.3,2x,'Y=',f6.3,
     2	2x,'Z=',f6.3)
	 write(2,6)methode,preci,m(n)*.1	!*.1 a cause du format 1P
6	 format(///,t6,'Parametres numeriques',//,t2,'methode:',a60,/,
     1	t2,'precision=',1pd10.3,', masse de la couche externe ',
     2	f5.2,'Msol',//,
     3	1x,'---------------------------------------------------',///)
	
	 do j=2,nchim
	  if(nom_elem(j) .eq. ' H2')then
	   j2=j	
	  elseif(nom_elem(j) .eq. 'He3')then
	   j3=j
	  elseif(nom_elem(j) .eq. 'He4')then
	   j4=j
	  elseif(nom_elem(j) .eq. 'Li7')then
	   j7=j	
	  elseif(nom_elem(j) .eq. 'Be7')then
	   j7b=j	
	  elseif(nom_elem(j) .eq. 'C12')then
	   j12=j
	  elseif(nom_elem(j) .eq. 'C13')then
	   j13=j
	  elseif(nom_elem(j) .eq. 'N14')then
	   j14=j
	  elseif(nom_elem(j) .eq. 'N15')then
	   j15=j
	  elseif(nom_elem(j) .eq. 'O16')then
	   j16=j
	  elseif(nom_elem(j) .eq. 'O17')then
	   j17=j
	  endif
	 enddo
	 abon=j3*j4*j12*j13*j14*j15*j16*j17 .gt. 0
	endif
	logg=log10(msol*mstar/rsol/rsol*g)
	
	if(ecritout)then	!listing ou en tete, suivant le cas
	 if(sautp)then	!on vient de faire un saut de page
	  WRITE(2,13)age*1.d-3,log10(teff),log10(l(n)),rstar,
     1	logg+log10(m(n)/rstar**2),p(1),t(1),
     2	compg(nbelem*(n-1)+1),compg(nbelem*(n-1)+ihe4),d2p,d2ro,
     3	i_pp,i_cno,i_3a,i_gr,chaine
13	  format(1x,'age(milliard annees)=',1pd10.3,1x,'Log(temp. eff.)=',
     1	1pd11.4,1x,'Log(luminosite/Lsol)=',1pd11.4,3x,'Rstar/Rsol=',
     2	1pd10.3,/,1x,'Log10g=',1pd10.3,1x,'Pc=',1pd10.3,1x,
     3	'Tc=',1pd10.3,1x,'Xc=',1pd10.3,1x,'Yc=',1pd10.3,/,
     4	1x,'(Rstar**2/p d2p/dr2)c=',1pd10.3,1x,
     5	'(Rstar**2/ro d2ro/dr2)c=',1pd10.3,/,
     6	1x,'en. PP=',i3,'%, en. CNO=',i3,'%, en. 3 alpha=',i3,
     7	'%, en. grav=',i3,'%, ',a30,/)
	 else
	  WRITE(2,7)age*1.d-3,log10(teff),log10(l(n)),rstar,
     1	logg+log10(m(n)/rstar**2),p(1),t(1),
     2	compg(nbelem*(1-1)+1),compg(nbelem*(1-1)+ihe4),d2p,d2ro,
     3	i_pp,i_cno,i_3a,i_gr,chaine
7	  format(1x,'age(milliard annees)=',1pd10.3,1x,'Log(temp. eff.)=',
     1	1pd11.4,1x,'Log(luminosite/lsol)=',1pd11.4,3x,'Rstar/Rsol=',
     2	1pd10.3,/,1x,'Log10g=',1pd10.3,1x,'Pc=',1pd10.3,1x,
     3	'Tc=',1pd10.3,1x,'Xc=',1pd10.3,1x,'Yc=',1pd10.3,/,
     4	1x,'(Rstar**2/p d2p/dr2)c=',1pd10.3,1x,
     5	'(Rstar**2/ro d2ro/dr2)c=',1pd10.3,/,
     6	1x,'en. PP=',i3,'%, en. CNO=',i3,'%, en. 3 alpha=',i3,
     7	'%, en. grav=',i3,'%, ',a30,/)
	 endif
	 write(2,26)1.d0-mstar/mtot,mstar,mtot
	 write(2,34)(nom_elem(i),i=1,min(nchim,12))
34	 format(t20,'abondances/volume des elements en surface',/,
     1	1x,12(4x,a3,3x))
	 write(2,31)(compg(nbelem*(n-1)+i)*ro(n)/amu/nucleo(i),
     1	i=1,min(nchim,12))
 
	 if (abon)then	
	  write(2,35)compg(nbelem*(n-1)+j3)/compg(nbelem*(n-1)+j4)*
     1	nucleo(j4)/nucleo(j3),
     2	compg(nbelem*(n-1)+j13)/compg(nbelem*(n-1)+j12)*
     3	nucleo(j12)/nucleo(j13),
     4	compg(nbelem*(n-1)+j15)/compg(nbelem*(n-1)+j14)*
     5	nucleo(j14)/nucleo(j13),	
     6	compg(nbelem*(n-1)+j17)/compg(nbelem*(n-1)+j16)*
     7	nucleo(j16)/nucleo(j17)	
35	  format(1x,/,t20,'Rapports isotopiques en nombre',/,
     1	'He3/He4=',1pd10.3,', C13/C12=',1pd10.3,
     2	', N15/N14=',1pd10.3,', O17/O16=',1pd10.3,/)
	 endif	
	else
	 WRITE(2,14)age*1.d-3,log10(teff),log10(l(n)),rstar,
     1	logg+log10(m(n)/rstar**2),p(1),t(1),compg(nbelem*(1-1)+1),
     2	compg(nbelem*(1-1)+ihe4),i_pp,i_cno,i_3a,i_gr,chaine
14	 format(1x,'age(milliard annees)=',1pd10.3,1x,'Log(temp. eff.)=',
     1	1pd11.4,1x,'Log(luminosite/lsol)=',1pd11.4,3x,'Rstar/Rsol=',
     2	1pd10.3,/,1x,'Log10g=',1pd10.3,1x,'Pc=',1pd10.3,2x,
     3	'Tc=',1pd10.3,2x,'Xc=',1pd10.3,2x,'Yc=',1pd10.3,/,
     4	1x,'en. PP=',i3,'%, en. CNO=',i3,'%, en. 3 alpha=',i3,
     5	'%, en. grav=',i3,'%, ',a30,/)
	 sautp=.false.	!il n'y a pas eu saut de page
	 write(2,26)mstar/mtot-1.d0,mstar,mtot
26	 format(1x,'Var. rel. de masse=',1pd10.3,' mstar=',1pd10.3,
     1	' mtot=',1pd010.3,/)
	
	 write(2,34)(nom_elem(i),i=1,min(nchim,12))
	 write(2,31)(compg(nbelem*(n-1)+i)*ro(n)/amu/nucleo(i),
     1	i=1,min(nchim,12))
	 if (abon)then	
	  write(2,35)compg(nbelem*(n-1)+j3)/compg(nbelem*(n-1)+j4)*
     1	nucleo(j4)/nucleo(j3),
     2	compg(nbelem*(n-1)+j13)/compg(nbelem*(n-1)+j12)*
     3	nucleo(j12)/nucleo(j13),
     4	compg(nbelem*(n-1)+j15)/compg(nbelem*(n-1)+j14)*
     5	nucleo(j14)/nucleo(j13),	
     6	compg(nbelem*(n-1)+j17)/compg(nbelem*(n-1)+j16)*
     7	nucleo(j16)/nucleo(j17)	
	 endif	
	 return
	endif
	write(2,8)
8	format(1x,'couche',t10,'m/Mstar',t21,'R/Rstar',t30,'pression',
     1	t42,'temp.',t50,'densite',t60,'L/Lsol',t69,'en. int.',
     2	t80,'gradrad',t91,'kappa',t99,'epsilon',t112,
     3	'TdS',t121,'eps3AL',/,
     4	t9,'1-m/Mstar',t21,'mu mol.',t33,'h+',t43,'he+',t52,
     5	'he++',t60,'pgas/P',t70,'gamma1',t80,'gradconv',
     6	t89,'(dK/dT)ro',t99,'desp/dT',t111,'epsPP',
     7	t121,'Z_opa',/,
     8	t10,'alpha',t21,'delta',t33,'cp',t40,'ech.ht.p',t51,'degene',
     9	t60,'v. son',t69,'vaissala',t81,'gradad',t89,
     1	'(dK/dro)T',t99,'deps/dro',t111,'epsCNO',t121,'Omega')
	
c	pour faire un saut de page, remplacer 1x par '1' dans le format 81	
	
81	format(1x,'couche',t10,'m/Mstar',t21,'R/Rstar',t30,'pression',
     1	t42,'temp.',t50,'densite',t60,'L/Lsol',t69,'en. int.',
     2	t80,'gradrad',t91,'kappa',t99,'epsilon',t112,
     3	'TdS',t121,'eps3AL',/,
     4	t9,'1-m/Mstar',t21,'mu mol.',t33,'h+',t43,'he+',t52,
     5	'he++',t60,'pgas/P',t70,'gamma1',t80,'gradconv',
     6	t89,'(dK/dT)ro',t99,'desp/dT',t111,'epsPP',
     7	t121,'Z_opa',/,
     8	t10,'alpha',t21,'delta',t33,'cp',t40,'ech.ht.p',t51,'degene',
     9	t60,'v. son',t69,'vaissala',t81,'gradad',t89,
     1	'(dK/dro)T',t99,'deps/dro',t111,'epsCNO',t121,'Omega')
	write(2,9)(nom_elem(i),i=1,min(nchim,12))
9	format(7x,12(4x,a3,3x))
	iln=2+1+3+4 !on a ecrit 10 lignes (2 lignes identification machine)
	convp=.false.	!pour reperer les zones convectives
	anupps=0.	!pour integrer le flux de neutrinos
	anupeps=0.
	anub8s=0.
	anube7s=0.
	anun13s=0.
	anu015s=0.
	anufs=0.           !modif GB
 
	do i=n,1,-1	!pour chaque couche
	 if(i .ne. 1)then	!integration de la production de neutrinos
	  dm=(m(i)-m(i-1))*mstar	!formule des trapezes
	  anupps=anupps+(anupp(i-1)+anupp(i))*dm
	  anupeps=anupeps+(anupep(i-1)+anupep(i))*dm
	  anub8s=anub8s+(anub8(i-1)+anub8(i))*dm
	  anube7s=anube7s+(anube7(i-1)+anube7(i))*dm
 
	  anun13s=anun13s+(anun13(i-1)+anun13(i))*dm
	  anu015s=anu015s+(anu015(i-1)+anu015(i))*dm
	  anufq=anupp(i)*1.102d-4*(1.+0.02*t(i)*1.d-6)*ro(i)*	!modif GB
     1	(1.+compg(nbelem*(i-1)+1))/2./dsqrt(t(i)*1.d-6)
	  anufp=anupp(i-1)*1.102d-4*(1.+0.02*t(i-1)*1.d-6)*ro(i-1)* !modif GB
     1	(1.+compg(nbelem*(i-2)+1))/2./dsqrt(t(i-1)*1.d-6)
	  anufs=anufs+(anufp+anufq)*dm
	 endif
 
	 iln=iln+5		!on va ajouter 5 lignes
	 if(iln .gt. 60)then	!peux-t-on ?
 
	  write(2,81)		!en tete
	  write(2,9)(nom_elem(j),j=1,min(nchim,12))
 
	  iln=2+4+5		!2,+ en-tete, +6 lignes
	 endif
	 if(conv(i) .NEQV. convp)then
	  iln=iln+2
	  if(iln .gt. 60)then	!peux-t-on ?
	   write(2,81)		!en tete 	
	   write(2,9)(nom_elem(j),j=1,min(nchim,12))
	   iln=2+4+4+1		!2,+ en-tete, +4 lignes
	  endif
	  write(2,*)' '
	  if(conv(i))then
	   write(2,*)'---------------debut de zone convective--------------'
		else
	   write(2,*)'---------------fin de zone convective----------------'
	  endif
	 endif
	 convp=conv(i)
	 cv=cp(i)-p(i)*delta(i)**2/ro(i)/t(i)/alfa(i)
	 gamma1=1./(alfa(i)-delta(i)*gradad(i))
 
	 if(max(nh1(i),nhe1(i),nhe2(i)) .lt. 0.)then !taux d'ionis. ???
	  mu_mol=1./(compg(nbelem*(i-1)+1)/nucleo(1)*(2./ah-3./ahe4)+
     1	3./ahe4*(1.-z(i))+.5*z(i))
	 else
	  mu_mol=compg(nbelem*(i-1)+1)/nucleo(1)*(1./ah-1./ahe4)+(1-z(i))/ahe4+
     1	z0/16.+compg(nbelem*(i-1)+1)/nucleo(1)*
     2	(nh1(i)/ah-(nhe1(i)+2.*nhe2(i))/ahe4)
     3	+(1.-z(i))*(nhe1(i)+2.*nhe2(i))/ahe4+z(i)*.5
	  mu_mol=1./mu_mol
	 endif
	 v_son(i)=sqrt(gamma1*p(i)/ro(i))
	 r_son(i)=r(i)*rsol
	 write(2,12)i,m(i)/mstar,r(i)/rstar,p(i),t(i),ro(i),
     1	l(i),eint(i),gradrad(i),kap(i),epsilon(5*(i-1)+1),
     2	epsilon(5*(i-1)+5),epsilon(5*(i-1)+4),
     3	1.d0-m(i)/mstar,mu_mol,nh1(i),nhe1(i),nhe2(i),beta(i),gamma1,
     4	gradconv(i),dcapdt(i),depsdt(i),epsilon(5*(i-1)+2),z(i),
     5	alfa(i),delta(i),cp(i),hp(i),degene(i),v_son(i),
     6	vaissala(i),gradad(i),dcapdr(i),depsdr(i),epsilon(5*(i-1)+3),
     7	w(i)
	 write(2,16)(compg(nbelem*(i-1)+j),j=1,min(nchim,12))
 
12	 format(/,1x,i6,1p12d10.3,/,(7x,1p12d10.3))
16	 format(7x,1p12d10.3)
	enddo	!i
 
	write(2,*)' '
	write(2,*)'------------fin du modele-----------------'
 
	if(n_atm .gt. 0)then
	 write(2,23)
23	 format(//,t40,'ATMOSPHERE',//,t6,'ep.opt.',3x,'pression',1x,
     1	'temperature',2x,'R/R* - 1',3x,'M/M* - 1',5x,'ro',6x,
     2	'log10(k)',2x,'grad.rad.',1x,'grad.conv',3x,'grad.ad',
     3	2x,'convec',/)
	 do i=1,n_atm
	  write(2,24)tau(i),p_atm(i),t_atm(i),r_atm(i)/rstar-1.d0,
     1	m_atm(i)/mstar-1.d0,ro_atm(i),log10(k_atm(i)),gradr_atm(i),
     2	gradc_atm(i),grada_atm(i),(gradr_atm(i) .gt. grada_atm(i))
 
24	  format(1x,1p10d11.4,3X,l1)
	  if(mod(i,5) .eq. 0)write(2,*)' '
	 enddo
	endif
 
c	depletion des elements
	
	write(2,30)(nom_elem(i),i=1,min(nchim,12))
30	format(//,t20,'DEPLETION DES ELEMENTS EN SURFACE',//,1x,12(4x,a3,3x))
	write(2,31)(12.-log10(compg(nbelem*(n-1)+1)/nucleo(1))
     1	+log10(compg(nbelem*(n-1)+i)/nucleo(i)),i=1,min(nchim,12))
31	format(1x,1p12d10.3)
 
c	le Z/X final
 
	if(ihe4 .le. 1)then
	 zsx=z0
	else
	 zsx=1.d0 -compg(nbelem*(n-1)+ihe4)-compg(nbelem*(n-1)+1)
	endif
	zsx=zsx/compg(nbelem*(n-1)+1)
	if(w(n) .gt. 0.d0)then
	 write(2,32)zsx,2.*pi/24./3600./w(n),rstar*w(n)*rsol*1.d-5
32	 format(//,t4,' Z / X du modele:',1pd10.3,//,
     1	t4,'Periode de rotation a la surface',1pd10.3,' jours,',
     2	' vitesse de rotation',1pd10.3,' km/s',//)
	else
	 write(2,33)zsx
33	 format(//,t4,' Z / X du modele:',1pd10.3,//)
	endif	
 
c	calcul des flux de neutinos a une distance de 1UA
 
	anupps=anupps*cte3
	anupeps=anupeps*cte3
	anub8s=anub8s*cte3
	anube7s=anube7s*cte3
	anun13s=anun13s*cte3
	anu015s=anu015s*cte3
        anufs=anufs*cte3	!modif GB
 
c	write(2,15)anupps,anupeps,anub8s,anube7s,anun13s,anu015s,
c	1	anupps/cte2,anupeps/cte2,anub8s/cte2,anube7s/cte2,
c	2	anun13s/cte2,anu015s/cte2
	write(2,15)anupps,anufs,anupeps,anub8s,anube7s,anun13s, !modif GB
     1    anu015s,anupps/cte2,anufs/cte2,anupeps/cte2,anub8s/cte2,
     2	anube7s/cte2,anun13s/cte2,anu015s/cte2,
     3	anupps*11.8d-10/cte2,anufs*215.d-10/cte2,
     4	anupeps*72.7d-10/cte2,anub8s*2.43d-6/cte2,
     5	anube7s*61.8d-10/cte2,anun13s*116.d-10/cte2,
     6	anu015s*117.d-10/cte2,
     7	anupps*0.d0/cte2,anufs*16.d-10/cte2,
     8	anupeps*2.38d-10/cte2,anub8s*1.06d-6/cte2,
     9	anube7s*1.66d-10/cte2,anun13s*6.61d-10/cte2,
     1	anu015s*6.67d-10/cte2
15	 format(//,t30,'N E U T R I N O S',//,			!modif GB
     1	t24,'anupp',t34,'anupep',t44,'anube7',t54,'anub8',
     2	t64,'anun13',
     3	t74,'anuo15',t84,'anuf17',/,1x,'production',9x,1p7d10.3,/
     4	1x,'flux sur terre',5x,1p7d10.3,/1x,'pour le gallium',4x,
     5	1p7d10.3,/1x,'pour le chlore',5x,1p7d10.3)
	
c	flux en SNU pour le chlore et le gallium, sections efficaces selon
c	G. Berthomieu
 
	write(2,17)(anupps*11.8d-10+anupeps*72.7d-10+anub8s*2.43d-6+ !GB
     1	anube7s*61.8d-10+anun13s*116.d-10+anufs*215.d-10
     2	+anu015s*117.d-10)/cte2,
     3	(anupps*0.d0+anupeps*2.38d-10+anub8s*1.06d-6+
     4	anube7s*1.66d-10+anun13s*6.61d-10 +anufs*16.d-10
     5	+anu015s*6.67d-10)/cte2,anub8s*cte4
17	format(/,1x,'Pour le Gallium:',1pd10.3,' SNU,',' pour le Chlore: ',
     1	1pd10.3,' SNU,', ' pour B8/modele:',1pd10.3)
	
c	estimation de Nu0 et dnu02 et dnu13
 
	call dnunl_3(r_son,v_son,n,rstar,nu0,dnu02,dnu13,a)	
	write(2,18)nu0,dnu02,dnu13,a
18	format(///,'Valeurs approx. de Nu0, delta Nu02, delta Nu13 et A',
     1	//,1x,'Nu0=',1pd10.3,' dnu02=',1pd10.3,' dnu13=',1pd10.3,
     2	' A=',1pd10.3)
 
c	gestion des pages	
	
	write(2,19)
19	format(//,1x,'---------------------------------------------------',//)	
	write(2,11)	!page suivante
11	format('1')
	sautp=.true.	!on se souvient qu'il y a eu saut de page
 
	return
 
 
	end
