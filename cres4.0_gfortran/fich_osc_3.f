c******************************************************************
 
	subroutine fich_osc_3(etat,opa,conv,nuc,cte)
 
c	formation d'un fichier pour calcul des oscillations
c	est appele par le programme for037
 
c	12 09 96, modif dans la recherche des limites ZR/ZC
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	CESAM3
 
c	itot: nombre de points
c	iconst: nombre de constantes global
c	ivar: nombre de variables du tableau var sans les elements chimiques
c	modele: nom du modele
c	nom_elem: nom des elements utilises
c	nbelem: nombre d'elements utilises
c	nchim: nombre d'elements chimiques utilises	
c	glob: variables globales
c		glob(1)=mstar*msol
c		glob(2)=rstar*rsol
c		glob(3)=ltot*lsol
c		glob(4)=z
c		glob(5)=x0
c		glob(6)=alpha
c		glob(7)=9./4.
c		glob(8)=1./162.
c		glob(9)=1.
c		glob(10)=1.
c		glob(11)=d2p
c		glob(12)=d2ro
c		glob(13)=age
c		glob(14)=vsal
c		glob(15)=0.
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
c	 	var(16,i)=w(i)
c	 	var(17,i)=d ln kappa / d ln T
c	 	var(18,i)=d ln kappa / d ln ro
c	 	var(19,i)=d epsilon(nuc) / d ln T
c	 	var(20,i)=d epsilon(nuc) / d ln ro
c	  	var(20+j,i)=xchim(j), j=1,nbelem
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'modele_3.common'
	include 'evol_chim_3.common'
	include 'atmosphere_3.common'
	include 'ctephy.common'
 
	integer pnosc
	parameter (pnosc=3500)		!nb de points maximum
 
	integer i,n_add,itot,j,n,long,ns,iconst,ivar,ivers,lim,
     1	nd,jlim(2*pnzc),k,compt,compt_max,n_lim,nadl

	external long

	real*8	p,t,r,l,m,ro,drop,drot,drox,u,dup,dut,dux,pint(pnosc),
     2	gradient,dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,m_zc(2*pnzc),
     3	epsilon(5),depsp,depst,depsx(pnelem),kap,alfa,w,z,
     4	delta,cp,hp,gradad,gradrad,vaissala,beta,d2p,d2ro,mstar,
     5	nh1,nhe1,nhe2,xchim(pnelem),dxchim(pnelem),q_int,cte5,
     6	rstar,ltot,dkapp,dkapt,dkapx,age,r_dble,q_lim(2*pnzc),
     7	fonc(0:5),absc(0:5),poly(0:5),dcapdt,dcapdr,m_dble,depsdt,
     8	depsdr,y,grad_mj,glob(15),var(20+pnelem,pnosc),
     9	qint(pnosc),nuc_m(pnelem)
	
	real*8	v1,v2,v3,qzc1,qzc2,qzc3,dqz,exit,adl,add
 
	character*1 oui
	character*2 c2
	character*3 c3,cpt
	character*4 c4
	character*5 m_ou_r
	character*9 today
	character*50 nom_fich,fich_int
	character*80 titre
 
	logical atm,pt_dble,pt_dbl
	
c	pt_dble=.true. : le point est point double
c	pt_dbl=.true. : on met des points doubles aux limites ZR/ZC
 
	external etat,opa,conv,nuc,cte
 
2000	format((1x,1p8d10.3))
 
	write(6,*)'entrer le nom generique des fichiers du modele'
	write(6,*)'Exemple: soleil'
	read(5,'(a)')nom_fich2
	write(6,*)'le nom du fichier binaire de la structure interne est donc:'
	write(6,*)nom_fich2(:long(nom_fich2))//'_B.dat'
	write(6,*)'o/n ?'
	read(5,'(a)')oui
	if(oui .ne. 'o')then
	 write (6,*)'renommer le fichier de la structure interne'
	 stop
	endif
	
c	lecture du fichier de donnees	
	
	call lit_nl_3
	
c	avec ou sans atmosphere		!
 
	write(6,*)'le modele a-t-il une atmosphere reconstituee ? o/n'
	read(5,'(a)')oui
	atm=oui .eq. 'o'
	if(atm)then
	 write(6,*)'le nom du fichier binaire de l''atmosphere est donc:'
	 write(6,*)nom_fich2(:long(nom_fich2))//'_B.atm'
	 write(6,*)'o/n ?'
	 read(5,'(a)')oui
	 if(oui .ne. 'o')then
	  write (6,*)'renommer le fichier de l''atmosphere'
	  stop
	 endif
	 write(6,*)'faut-il tenir compte de l''atmosphere dans le fichier'
	 write(6,*)'d''oscillations a creer ? o/n'
	 read(5,'(a)')oui
	 atm=oui .eq. 'o'
	endif
 
	q_int=2.d0		!initialisation, calcul de v1 (d_grad) et qzc1
	m_ou_r='  q  '
	call intertot_3(atm,m_ou_r,dxchim,q_int,nd,jlim,lim,m_zc,
     1	p,t,xchim,m,l,r,ro,drop,drot,drox,u,dup,dut,dux,rstar,ltot,
     2	gradient,dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,mstar,w,
     3	epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,alfa,grad_mj,
     4	delta,cp,hp,gradad,gradrad,vaissala,beta,age,n,ns,
     5	nh1,nhe1,nhe2,etat,opa,conv,nuc,cte)
	v1=gradrad-gradad
	qzc1=q_int
 
	write(6,*)'faut-il multiplier/interpoler le modele?'
	write(6,*)'multiplier ? o/n'
	read(5,'(a)')oui
	if(oui .eq. 'o')then
	 m_ou_r='  q  '
	 write(6,*)'nombre de points a inserer entre deux couches ?'
	 read(5,*)n_add
 
c	 on construit le fichier des qint des points d'interpolation
 
	 if(pnosc .lt. (ns+1)*n_add+4*pnzc)then
	  write(6,*)'dans fich_osc_3 mettre le parametre pnosc a ',
     1	(ns+1)*n_add+4*pnzc
	  stop
	 endif
 
	 add=1./float(n_add+1)
	 nadl=9		!ajout autour des limites ZR/ZC
	 adl=1./float(nadl+1)
	 k=1	
	 itot=1
	 qint(itot)=1.01d0
	 i=1
	 do while(i .le. n-1)
	  if(jlim(k)-1 .eq. i)then
	   k=min(lim,k+1)
	   itot=itot+1
	   qint(itot)=i		!le point de grille avant la limite	
	   do j=1,nadl		!nadl points intermediaires avant
	    itot=itot+1
	    qint(itot)=i+adl*j
	   enddo
	   i=i+1
	   itot=itot+1
	   qint(itot)=i		!le point de grille de la limite
	   do j=1,nadl		!nadl points intermediaires apres
	    itot=itot+1
	    qint(itot)=i+adl*j
	   enddo
 
	  else
	   itot=itot+1
	   qint(itot)=i		!les points de la grille
	   do j=1,n_add		!points intermediaires
	    itot=itot+1
	    qint(itot)=i+add*j
	   enddo
	  endif
	  i=i+1
	 enddo
	
c	 atmosphere	
 
	 do i=n,ns	
	  itot=itot+1
	  qint(itot)=i
	 enddo
	
c	 print*,n,ns,itot,lim,(jlim(i),i=1,lim)
c	 write(6,2000)(qint(i),i=1,itot)
c	 pause
	
	 call shell(itot,qint)
	
	else
	 write(6,*)'on interpole donc, soit en masse/Mtot soit en rayon/Rsol'
	 write(6,*)'en masse ? o/n'
	 read(5,'(a)')oui
	 if(oui .eq. 'o')then
	  m_ou_r='masse'
	 else
	  m_ou_r='rayon'
	 endif
	 write(6,*)'nom du fichier ASCII ou se trouvent les points'
	 write(6,*)'ou interpoler ?'
	 read(5,'(a)')fich_int
	 open(unit=25,form='formatted',status='old',file=fich_int)
	 do itot=1,pnosc
	  read(25,*,end=10)i,pint(itot)
c	  print*,itot,pint(itot)
	 enddo
	 write(6,*)'dans fich_osc_3 mettre le parametre pnosc au moins egal'
	 write(6,*)'au nombre+1 de points d''interpolation desires'
	 close(unit=25)
	 stop
10	 itot=itot-1
c	 print*,itot
	 close(unit=25)
 
c	 on retourne les donnees
 
	 do i=1,itot		!p: vt
	  p=pint(i)
	  pint(i)=pint(itot-i+1)
	  pint(itot-i+1)=p
	 enddo
	endif
 
c	points doubles	
c	pt_dble=.true. : le point est point double
c	pt_dbl=.true. : on met des points doubles aux limites ZR/ZC
 
	if(m_ou_r .eq. '  q  ')then
	  	
	 write(6,*)'met-on un point double a chaque limite ZR/ZC? o/n'
	 read(5,'(a)')oui
	 pt_dbl=oui .eq. 'o'
	 if(pt_dbl)then
	  cpt='_db'
	 else		!sans point double
	  cpt='_sp'
	 endif
 
c	 on affine les limites ZC/ZR pour y mettre un point
c	 v1, v2, v3 : d_grad aux points qzc1 < qzc2 < qzc3
c	 quand on avance v1 est le d_grad du point precedent
 
	 if(pt_dbl)then
	  dqz=1.d-4	!precision pour la dichotomie
	  compt_max=30	!nombre maximal de bissections
	  exit=adl/3.	!pour eliminer les points trop pres des limites	
	  n_lim=0
	  do i=3,ns
	   qzc2=i
	   call intertot_3(atm,'  q  ',dxchim,qzc2,nd,jlim,lim,m_zc,
     1	p,t,xchim,m,l,r,ro,drop,drot,drox,u,dup,dut,dux,rstar,ltot,
     2	gradient,dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,mstar,w,
     3	epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,alfa,grad_mj,
     4	delta,cp,hp,gradad,gradrad,vaissala,beta,age,n,ns,
     5	nh1,nhe1,nhe2,etat,opa,conv,nuc,cte)
	   v2=gradrad-gradad
c	   write(6,2000)v1,v2,v1*v2,float(i)
	   if(v2*v1 .lt. 0.)then
	    compt=1		!dichotomie
	    do while(abs(qzc1-qzc2) .gt. dqz .and. compt .le. compt_max)
	     qzc3=(qzc1+qzc2)/2.
	     call intertot_3(atm,'  q  ',dxchim,qzc3,nd,jlim,lim,m_zc,
     1	p,t,xchim,m,l,r,ro,drop,drot,drox,u,dup,dut,dux,rstar,ltot,
     2	gradient,dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,mstar,w,
     3	epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,alfa,grad_mj,
     4	delta,cp,hp,gradad,gradrad,vaissala,beta,age,n,ns,
     5	nh1,nhe1,nhe2,etat,opa,conv,nuc,cte)
	     v3=gradrad-gradad	
c	     print*,'dichotomie',qzc1,qzc2,qzc3,compt
c	     write(6,2000)v1,v2,v3
	     if(v1*v3 .gt. 0.)then
	      qzc1=qzc3
	      v1=v3
	     else
	      qzc2=qzc3
	      v2=v3
	     endif
	     compt=compt+1
	    enddo		!while de la dichotomie
	
c	    print*,v1,qzc1,v2,qzc2
c	    pause
	
	    if(compt .gt. compt_max)then
	     print*,'pas de conv dans la dicho, localis. approx. pour limite',i
 
	    else
 
c	     on met en 3 (v3, qz3) le cote positif (ZC) de d_grad
	
	     if(v2 .gt. 0.)then
	      v3=v2
	      qzc3=qzc2
	     else
	      v3=v1
	      qzc3=qzc1
	     endif
	
	     itot=itot+1	
	     qint(itot)=qzc3	!on garde le point cote d_grad negatif
	     n_lim=n_lim+1
	     q_lim(n_lim)=qzc3
	     print*,'qzc3',qzc3
	
	    endif		!convergence de la dichotomie
	   endif		!inversion de signe  v1*v2 .lt. 0.
	   v1=v2		!2--->1
	   qzc1=qzc2
	  enddo		!boucle sur ns
	
c	  on sort les couches trop pres des q_lim, on trie on elimine les extra	
 
	  do i=1,itot
	   do j=1,n_lim
	    if(abs(qint(i)-q_lim(j)) .le. exit .and. qint(i) .ne. q_lim(j))
     1	qint(i)=1.d5
	   enddo
	  enddo
	
	  call shell(itot,qint)	!tri
	
	  do while (qint(itot) .gt. ns)
	   itot=itot-1
	  enddo
	
c	  installation des points doubles
	
	  do i=1,n_lim	
	   itot=itot+1
	   qint(itot)=q_lim(i)+.01
	  enddo
	  call shell(itot,qint)	!tri
	 endif
	
	else		!pas de point doubles si m_ou_r .ne. ' q '
	 cpt='_sp'
	 pt_dbl=.false.
	endif		!si m_ou_r .eq. ' q '
	
	if(itot .ge. pnosc)then
	 print*,'trop de points demandes, arret'
	 print*,'dans fich_osc_3, mettre le parametre pnosc=',itot+10
	 stop
	endif
 
	write (6,*)'nb. de pts. du fichier d''oscillations a creer:',itot
	if(.true.)then	!suppression du nb. de couches dans le nom du fichier
c	if(.false.)then
	 if(itot .le. 99)then
c	  encode(2,2,c2)itot
          write(c2,2) (itot)
2	  format(i2)
	  nom_fich=nom_fich2(:long(nom_fich2))//'_'//c2//cpt//'.osc'
	 elseif(itot .le. 999)then
c	  encode(3,3,c3)itot
          write(c3,3) (itot)
3	  format(i3)
	  nom_fich=nom_fich2(:long(nom_fich2))//'_'//c3//cpt//'.osc'
	 elseif(itot .le. 9999)then
c	  encode(4,4,c4)itot
          write(c4,4) (itot)
4	  format(i4)
	  nom_fich=nom_fich2(:long(nom_fich2))//'_'//c4//cpt//'.osc'
	 endif
	else	
	 nom_fich=nom_fich2(:long(nom_fich2))//cpt//'.osc'
	endif	!suppression
	write(6,*)'on cree le fichier d''oscillations : ',
     1	nom_fich(:long(nom_fich))
	open(unit=30,form='formatted',status='unknown',file=nom_fich)
 
	titre='extension du modele: '//nom_fich2(:long(nom_fich2))
	if(atm)titre=titre(:long(titre))//' avec atm.'
	if(pt_dbl)then
	 titre=titre(:long(titre))//', avec pt. dbles'
	else
	 titre=titre(:long(titre))//', sans pt. dbles'	
	endif
	call date_and_time(today)
	titre=titre(:long(titre))//', du '//today
	
c	au centre pour d2p, d2ro
 
	cte5=-4.*pi*g/3.*rsol**2
 
	do i=0,5		!d2ro par newton au centre
	 q_int=1.d0+.05d0*i	!on utilise des points intermediaires
	 call intertot_3(atm,'  q  ',dxchim,q_int,nd,jlim,lim,m_zc,
     1	p,t,xchim,m,l,absc(i),fonc(i),drop,drot,drox,u,dup,dut,dux,
     2	rstar,ltot,
     3	gradient,dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,mstar,w,
     4	epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,alfa,grad_mj,
     5	delta,cp,hp,gradad,gradrad,vaissala,beta,age,n,ns,
     6	nh1,nhe1,nhe2,etat,opa,conv,nuc,cte)
	 if(i .eq. 0)then
	  ro=fonc(0)
	  d2p=cte5*ro**2*rstar**2/p
	 endif
c	 write(6,*)'i/p,t,xchim(1),m,l,q_int',i
c	 write(6,2000)p,t,xchim(1),m,l,q_int
	enddo
 
c	write(6,*)'ro/r/deriv'
c	write(6,2000)(fonc(i),i=0,5)
c	write(6,2000)(absc(i),i=0,5)
 
	call newton(0,m_qs+1,fonc,absc,0.d0,poly,2)
	d2ro=rstar**2/ro*poly(2)
 
c	write(6,2000)(poly(i),i=0,2)
c	write(6,*)'rstar,d2ro,ro,poly(2)'
c	write(6,2000)rstar,d2ro,ro,poly(2)
 
c	la composition chimique externe pour glob(9) et glob(10)
 
	q_int=ns	!a l'exterieur
	call intertot_3(atm,'  q  ',dxchim,q_int,nd,jlim,lim,m_zc,
     1	p,t,xchim,m,l,absc(i),fonc(i),drop,drot,drox,u,dup,dut,dux,
     2	rstar,ltot,
     3	gradient,dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,mstar,w,
     4	epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,alfa,grad_mj,
     5	delta,cp,hp,gradad,gradrad,vaissala,beta,age,n,ns,
     6	nh1,nhe1,nhe2,etat,opa,conv,nuc,cte)
	call chim_gram_3(xchim,dxchim,nuc_m)
	if(ihe4 .le. 1)then
	 y=1.-z-xchim(1)
	else
	 y=xchim(ihe4)
	endif
 
	write(30,*)titre(:long(titre))
	if(m_ou_r .eq. '  q  ')then
	 write(30,20)n_add
20 	 format(1x,'nb. de pts. entre 2 couches:',i3)
	elseif(m_ou_r .eq. 'masse')then
	 write(30,*)'inter. aux masses de '//fich_int(:long(fich_int))
	else
	 write(30,*)'inter. aux rayons de '//fich_int(:long(fich_int))
	endif
	write(30,*)'fichier d''oscillations: '//nom_fich(:long(nom_fich))
	write(30,*)'modele calcule par CESAM3 Splines/Collocation'
	write(30,90)nbelem,(nom_elem(i),i=1,nbelem)
90	format(i3,14(1x,a3))
 
	glob(1)=mstar*msol
	glob(2)=rstar*rsol
	glob(3)=ltot*lsol
	glob(4)=z0
	glob(5)=x0
	glob(6)=alpha
	glob(7)=9./4.
	glob(8)=1./162.
	glob(9)=xchim(1)
	glob(10)=y
	glob(11)=d2p
	glob(12)=d2ro
	glob(13)=age
	glob(14)=vsal
	glob(15)=w_rot
 
c	pause
 
c	le fichier d'oscillations
 
c	write(6,2000)(qint(k),k=1,itot)
 
	pt_dble=.false.
	do k=itot,1,-1
	 if(m_ou_r .eq. '  q  ')then
	  q_int=qint(k)
c	  write(6,2000)q_int
c	  write(6,*)'k,q_int,m,r,p,t,vaissala',k
	  call intertot_3(atm,m_ou_r,dxchim,qint(k),nd,jlim,lim,m_zc,
     1	p,t,xchim,m,l,r,ro,drop,drot,drox,u,dup,dut,dux,rstar,ltot,
     2	gradient,dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,mstar,w,
     3	epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,alfa,grad_mj,
     4	delta,cp,hp,gradad,gradrad,vaissala,beta,age,n,ns,
     5	nh1,nhe1,nhe2,etat,opa,conv,nuc,cte)
	
c	  if(q_int .lt. 3. .and. m .gt. 0.)write(6,2000)q_int,l,p,kap,m,t,l*p*kap/m/t**4
	
	  if(pt_dble)then	!on identifie masse et rayon
	   m=m_dble
	   r=r_dble
	  endif
c	  write(6,2000)q_int,m,r,p,t,vaissala
	 elseif(m_ou_r .eq. 'masse')then
	  m=pint(itot-k+1)
	  call intertot_3(atm,m_ou_r,dxchim,q_int,nd,jlim,lim,m_zc,
     1	p,t,xchim,m,l,r,ro,drop,drot,drox,u,dup,dut,dux,rstar,ltot,
     2	gradient,dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,mstar,w,
     3	epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,alfa,grad_mj,
     4	delta,cp,hp,gradad,gradrad,vaissala,beta,age,n,ns,
     5	nh1,nhe1,nhe2,etat,opa,conv,nuc,cte)
c	  write(6,*)'en m int,m,r',m,r
	 else
	  r=pint(itot-k+1)
c	  write(6,*)'k,r',k,r
	  call intertot_3(atm,m_ou_r,dxchim,q_int,nd,jlim,lim,m_zc,
     1	p,t,xchim,m,l,r,ro,drop,drot,drox,u,dup,dut,dux,rstar,ltot,
     2	gradient,dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,mstar,w,
     3	epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,alfa,grad_mj,
     4	delta,cp,hp,gradad,gradrad,vaissala,beta,age,n,ns,
     5	nh1,nhe1,nhe2,etat,opa,conv,nuc,cte)
c	  write(6,*)'en r int,m,r',m,r
	 endif
	
	 call chim_gram_3(xchim,dxchim,nuc_m)
	 if(ihe4 .le. 1)then
	  y=1.-z0-xchim(1)
	 else
	  y=xchim(ihe4)
	 endif
	 if(iz .gt. 1)then
	  z=xchim(iz)
	 else
	  z=1.-xchim(1)-y
	 endif
	
c	 print*,q_int
	
	 i=itot-k+1
	 var(1,i)=r*rsol
	 if(m .ne. 0.)then
	  var(2,i)=log(m/mstar)
	 else
	  var(2,i)=-1.d38
	 endif
	 var(3,i)=t
	 var(4,i)=p
	 var(5,i)=ro
	 var(6,i)=grad_mj
	 var(7,i)=l*lsol
	 var(8,i)=kap
	 var(9,i)=epsilon(1)
	 var(10,i)=1./(alfa-delta*gradad)
	 var(11,i)=gradad
	 var(12,i)=delta
	 var(13,i)=cp
	 var(14,i)=nh1*xchim(1)+(nhe1+2.*nhe2)*y+z*.5
	 if(m*r .ne. 0.)then
	  var(15,i)=vaissala
	 else
	  var(15,i)=0.
	 endif
	 var(16,i)=w
	 dcapdr=dkapp/drop
	 dcapdt=dkapt-dcapdr*drot
	 var(17,i)=dcapdt*t/kap
	 var(18,i)=dcapdr*ro/kap
	 depsdr=depsp/drop
	 depsdt=depst-depsdr*drot
	 var(19,i)=depsdt*t
	 var(20,i)=depsdr*ro
	 do j=1,nbelem
	  var(20+j,i)=xchim(j)
	 enddo
	
	 if(pt_dbl)then
	  pt_dble=.false.
	  do j=1,n_lim		!s'il y a un point double
	   if(q_int .eq. q_lim(j))pt_dble=.true.
	  enddo
	  if(pt_dble)then
	   r_dble=r	
	   m_dble=m
	  endif
	 endif
	
	enddo
	
	iconst=15		!nb de global
	ivar=20			!nb de varaiables
	ivers=3			!numero de version
	write(30,137)itot,iconst,ivar,ivers,nchim,iw,iz
137	format(7i10)
	write(30,138)(glob(i),i=1,iconst)
138	format(1p5d19.12)
	do j=1,itot
	 write(30,138)(var(i,j),i=1,ivar+nbelem)
	enddo
	close(unit=30)
 
	stop
 
	end
