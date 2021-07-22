c***********************************************************************
 
	subroutine desfull_3(etat,opa,conv,nuc,cte)
 
c	dessin d'un modele
 
 
c	interpolation du modele par splines de la base de de Boor
c	l'algorithme d'interpolation est repris
 
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
 
 
c	CESAM, version 3
 
 
	implicit none
 
 
	include 'cesam_3.parametres'
 
	include 'modele_3.common'
	include 'evol_chim_3.common'
 
	include 'atmosphere_3.common'
	include 'ctephy.common'
 
 
	integer pndes
 
	parameter (pndes=1000)
 
	integer i,itot,j,k,n,ns,long,j12,j13,j14,j15,j16,j17,nd,
     1	lim,jlim(2*pnzc)

	external long

c	data j12,j13,j14,j15,j16,j17/6*0/
 
	real*8	p,t,r,l,m,ro,drop,drot,drox,u,dup,dut,dux,gamma1,age,sum,
     2	gradient,dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,p8(pndes),
     3	epsilon(5),depsp,depst,depsx(pnelem),kap,alfa,m_zc(2*pnzc),
     4	delta,cp,hp,gradad,gradrad,vaissala,beta,n_nu(pnelem),grad_mj,
     5	nh1,nhe1,nhe2,xchim(pnelem),dxchim(pnelem),q_int,nuc_m(pnelem),
     6	rtot,ltot,dkapp,dkapt,dkapx,compchim(pnelem*pndes),mstar,
     7	q_inf,q_sup,pas,m_inf,m_sup,compg(pnelem*pndes)
 
	real*4	ymax,ymin,point(pndes),p4(pndes),xmax,xmin,c4(pndes),
     1	ro4(pndes),epsilo4(pndes),epsilo4pp(pndes),m4(pndes),
     2	epsilo4cno(pndes),epsilo43a(pndes),epsilo4tds(pndes),
     3	vaissala4(pndes),gradad4(pndes),gradrad4(pndes),t4(pndes),
     4	r4(pndes),l4(pndes),gradient4(pndes),kap4(pndes),w4(pndes)
 
	logical cno,atm
 
	character*1 oui
	character*5 m_ou_r
	character*10 elem(pnelem)
	character*50 text,en_x
	character*100 model,titre
	character*128 device
 
	external etat,opa,conv,nuc,cte
	data j12,j13,j14,j15,j16,j17/6*0/
 
2000	format((1x,1p10d10.3))
2	format(i2)
 
c	write(6,*)'pndes',pndes
 
	write(6,*)'le fichier des_full3.f correspond au modele a dessiner ?'
	write(6,*)'o/n ?'
	read(5,'(a)')oui
	if(oui .ne. 'o')then
	 write (6,*)'effectuer les rectifications necessaires'
	 stop
	endif
 
	write(6,*)'entrer le nom du fichier en binaire du modele a dessiner'
	write(6,*)'exemple: soleil pour soleil_B.dat'
	read(5,'(a)')nom_fich2
	write(6,'(a)')nom_fich2
	write(6,*)'le nom du fichier binaire du modele a dessiner est donc:'
	write(6,*)nom_fich2(:long(nom_fich2))//'_B.dat'
	write(6,*)'o/n ?'
	read(5,'(a)')oui
	if(oui .ne. 'o')then
	 write (6,*)'renommer le fichier du modele a dessiner'
	 stop
	endif
	
	call lit_nl_3		!lecture des namelist
	
	write(6,*)'le modele a-t-il une atmosphere reconstituee ? o/n'
	read(5,'(a)')oui
	atm=oui .eq. 'o'
	if(atm)then
	 write(6,*)'le nom du fichier binaire de l''atmosphere est donc:'
	 write(6,*)nom_fich2(:long(nom_fich2))//'_$.atm'
	 write(6,*)'o/n ?'
	 read(5,'(a)')oui
	 if(oui .ne. 'o')then
	  write (6,*)'renommer le fichier de l''atmosphere'
	  stop
	 endif
c	 write(6,'(a)')model_atm
	 write(6,*)'faut-il tenir compte de l''atmosphere pour le dessin, o/n'
	 read(5,'(a)')oui
	 atm=oui .eq. 'o'
	endif
 
	q_int=1.d0		!initialisation
	call intertot_3(atm,'  q  ',dxchim,q_int,nd,jlim,lim,m_zc,
     1	p,t,xchim,m,l,r,ro,drop,drot,drox,u,dup,dut,dux,rtot,ltot,
     2	gradient,dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,mstar,
     3	epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,alfa,grad_mj,
     4	delta,cp,hp,gradad,gradrad,vaissala,beta,age,n,ns,
     5	nh1,nhe1,nhe2,etat,opa,conv,nuc,cte)
 
	write(6,*)'Zoom ? o/n'
	read(5,'(a)')oui
	if(oui .eq. 'o')then	!cas du zoom determination de q_sup et q_inf
	 write(6,*)'entrer m inf, m sup'
	 read(5,*)m_inf,m_sup
c	 write(6,*)m_inf,m_sup
	 write(6,21)m_inf,m_sup
21	 format(1x,'Zoom entre m=',1pd10.3,' et m=',1pd10.3)
	 call intertot_3(atm,'masse',dxchim,q_inf,nd,jlim,lim,m_zc,
     1	p,t,xchim,m_inf,l,r,ro,drop,drot,drox,u,dup,dut,dux,rtot,ltot,
     2	gradient,dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,mstar,
     3	epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,alfa,grad_mj,
     4	delta,cp,hp,gradad,gradrad,vaissala,beta,age,n,ns,
     5	nh1,nhe1,nhe2,etat,opa,conv,nuc,cte)
c	 write(6,2000)q_inf,r,m_inf
	 call intertot_3(atm,'masse',dxchim,q_sup,nd,jlim,lim,m_zc,
     1	p,t,xchim,m_sup,l,r,ro,drop,drot,drox,u,dup,dut,dux,rtot,ltot,
     2	gradient,dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,mstar,
     3	epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,alfa,grad_mj,
     4	delta,cp,hp,gradad,gradrad,vaissala,beta,age,n,ns,
     5	nh1,nhe1,nhe2,etat,opa,conv,nuc,cte)
c	 write(6,2000)q_sup,r,m_sup
	else
	 q_inf=1
	 q_sup=ns
c	 write(6,*)n,ns
c	 write(6,2000)rtot
	endif
	print*,'nombre de points pour le dessin?'
	read*,itot
	itot=min(itot,pndes-1)
	pas=(q_sup-q_inf)/float(itot-1)	!pas en q
c	print*,q_inf,q_sup,pas,itot
 
c	le type de trace
 
	write(6,*)'trace en masse? entrer o/n'
	read(5,'(a)')oui
	if(oui .eq. 'o')then
	 m_ou_r='masse'
	else
	 write(6,*)'trace en rayon? entrer o/n'
	 read(5,'(a)')oui
	 if(oui .eq. 'o')then
	  m_ou_r='rayon'
	 else
	  write(6,*)'trace en q'
	  m_ou_r='  q  '
	 endif
	endif
 
	i=0
	q_int=q_inf
	do while (q_int .le. q_sup)
	 call intertot_3(atm,'  q  ',dxchim,q_int,nd,jlim,lim,m_zc,
     1	p,t,xchim,m,l,r,ro,drop,drot,drox,u,dup,dut,dux,rtot,ltot,
     2	gradient,dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,mstar,
     3	epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,alfa,grad_mj,
     4	delta,cp,hp,gradad,gradrad,vaissala,beta,age,n,ns,
     5	nh1,nhe1,nhe2,etat,opa,conv,nuc,cte)
 
c	 write(6,2000)q_int,(xchim(j),j=1,min(7,nbelem))
 
	 i=i+1
 
	 p4(i)=log10(p)
	 t4(i)=log10(t)
	 r4(i)=r/rtot
	 l4(i)=l
	 m4(i)=m
	 if(m_ou_r .eq. 'masse')then
	  point(i)=m4(i)
	 elseif(m_ou_r .eq. 'rayon')then
	  point(i)=r4(i)
	 else
	  point(i)=q_int
	 endif
	 gradient4(i)=gradient
	 gradad4(i)=gradad
	 gradrad4(i)=gradrad
	 vaissala4(i)=vaissala
	 ro4(i)=log10(ro)
	 gamma1=1./(alfa-delta*gradad)
	 c4(i)=sqrt(gamma1*p/ro)
	 epsilo4(i)=epsilon(1)
	 epsilo4pp(i)=epsilon(2)
	 epsilo4cno(i)=epsilon(3)
	 epsilo43a(i)=epsilon(4)
	 epsilo4tds(i)=epsilon(5)
	 kap4(i)=kap
	 if(iw .gt. 1 .and. r .gt. 0.)then
	  w4(i)=xchim(iw)/r**2
	 else
	  w4(i)=w_rot
	 endif
	 do j=1,nbelem			!par mole
	  compchim(nbelem*(i-1)+j)=xchim(j)
	 enddo	!j
	 call chim_gram_3(xchim,dxchim,nuc_m)
	 do j=1,nbelem			!par gramme
	  compg(nbelem*(i-1)+j)=xchim(j)
	 enddo	!j
	 	
	 q_int=q_int+pas
	enddo
	itot=i
	if(itot .gt. pndes)then
	 write(6,*)'trop de points: mettre le parametre pndes >',itot
	 stop
	endif
 
c	write(6,*)'sortie de la boucle'
c	if(.true.)stop
 
	write(text,10)age*1.d-3
10	format('   age=',f8.4,'10**9 ans')
	model=nom_fich2(:long(nom_fich2))//text
 
c	device='/XSERVE'
	write(6,*)'device ? /XSERVE, /SUN, /PS'
	read(5,'(a)')device
	call pgbegin(0,device,2,2)
 
	if(m_ou_r .eq. 'masse')then
	 en_x='Masse/Mtot'
	 xmin=point(1)
	 xmax=point(itot)
	elseif(m_ou_r .eq. 'rayon')then
	 en_x='Rayon/Rtot'
	 xmin=point(1)
	 xmax=point(itot)
	else
	 en_x='couche'
	 xmin=point(1)
	 xmax=point(itot)
	endif
 
	titre='Log Pression cgs'
	call pminmax(p4,itot,ymax,ymin)
	call pgenv(xmin,xmax,ymin,ymax,0,0)
	call pglabel(en_x,titre,model)
	call pgline(itot,point,p4)
 
	titre='Log Temperature cgs'
	call pminmax(t4,itot,ymax,ymin)
	call pgenv(xmin,xmax,ymin,ymax,0,0)
	call pglabel(en_x,titre,model)
	call pgline(itot,point,t4)
 
	titre='Log Densite cgs'
	call pminmax(ro4,itot,ymax,ymin)
	call pgenv(xmin,xmax,ymin,ymax,0,0)
	call pglabel(en_x,titre,model)
	call pgline(itot,point,ro4)
 
	titre='Rayon/Rtot'
	call pminmax(r4,itot,ymax,ymin)
	call pgenv(xmin,xmax,ymin,ymax,0,0)
	call pglabel(en_x,titre,model)
	call pgline(itot,point,r4)
 
	titre='Masse/Mtot'
	call pminmax(m4,itot,ymax,ymin)
	call pgenv(xmin,xmax,ymin,ymax,0,0)
	call pglabel(en_x,titre,model)
	call pgline(itot,point,m4)
 
	titre='Luminosite/Lsol'
	call pminmax(l4,itot,ymax,ymin)
	call pgenv(xmin,xmax,ymin,ymax,0,0)
	call pglabel(en_x,titre,model)
	call pgline(itot,point,l4)
 
	ymin=0.
	ymax=.6
	call pgsls(1)
	titre='Gradient ___ , ad _ _ , rad _._'
	call pgenv(xmin,xmax,ymin,ymax,0,0)
	call pglabel(en_x,titre,model)
	call pgline(itot,point,gradient4)
	call pgsls(2)
	call pgline(itot,point,gradad4)
	call pgsls(3)
	call pgline(itot,point,gradrad4)
 
	titre='V son cgs'
	call pgsls(1)
	call pminmax(c4,itot,ymax,ymin)
	call pgenv(xmin,xmax,ymin,ymax,0,0)
	call pglabel(en_x,titre,model)
	call pgline(itot,point,c4)
 
	titre='Vaissala = N'
	call pminmax(vaissala4,itot,ymax,ymin)
	ymax=min(ymax,2.5)
	call pgenv(xmin,xmax,ymin,ymax,0,0)
	call pglabel(en_x,titre,model)
	call pgline(itot,point,vaissala4)
 
	titre='opacite Rosseland'
	call pminmax(kap4,itot,ymax,ymin)
	call pgenv(xmin,xmax,ymin,ymax,0,0)
	call pglabel(en_x,titre,model)
	call pgline(itot,point,kap4)
	
	titre='rotation'
	call pminmax(w4,itot,ymax,ymin)
	if(max(abs(ymin),abs(ymax)) .gt. 0.)then
	 call pgenv(xmin,xmax,ymin,ymax,0,0)
	 call pglabel(en_x,titre,model)
	 call pgline(itot,point,w4)
	endif
 
	titre='Epsilon ___, PP _ _, CNO _._, 3 al .., grav _..._'
	call pminmax(epsilo4,itot,ymax,ymin)
	call pgenv(xmin,xmax,ymin,ymax,0,0)
	call pglabel(en_x,titre,model)
	call pgline(itot,point,epsilo4)
	call pgsls(2)
	call pgline(itot,point,epsilo4pp)
	call pgsls(3)
	call pgline(itot,point,epsilo4cno)
	call pgsls(4)
	call pgline(itot,point,epsilo43a)
	call pgsls(5)
	call pgline(itot,point,epsilo4tds)
	call pgsls(1)
 
	write(6,*)'abondances ? entrer o/n'
	read(5,'(a)')oui
	if(oui .ne. 'o')goto100
 
	do i=1,nbelem
	 elem(i)='Log 10 '//nom_elem(i)
	 if(nom_elem(i) .eq. ' H ')then
	  n_nu(i)=1
c	  ymax=0.5
c	  ymin=-5.1
	 elseif(nom_elem(i) .eq. ' H2')then
	  n_nu(i)=2
c	  ymax=-10.9
c	  ymin=-20.1
	 elseif(nom_elem(i) .eq. 'He3')then
	  n_nu(i)=3
c	  ymax=-1.9
c	  ymin=-9.1
	 elseif(nom_elem(i) .eq. 'He4')then
	  n_nu(i)=4
c	  ymax=0.1
c	  ymin=-2.1
	 elseif(nom_elem(i) .eq. 'Li7')then
	  n_nu(i)=7
c	  ymax=-6.9
c	  ymin=-16.1
	 elseif(nom_elem(i) .eq. 'Be7')then
	  n_nu(i)=7
c	  ymax=-7.9
c	  ymin=-17.1
	 elseif(nom_elem(i) .eq. 'C12')then
	  n_nu(i)=12
c	  ymax=0.5
c	  ymin=-6.1
	 elseif(nom_elem(i) .eq. 'C13')then
	  n_nu(i)=13
c	  ymax=-2.9
c	  ymin=-7.1
	 elseif(nom_elem(i) .eq. 'N14')then
	  n_nu(i)=14
c	  ymax=0.1
c	  ymin=-6.1
	 elseif(nom_elem(i) .eq. 'N15')then
	  n_nu(i)=15
c	  ymax=-3.9
c	  ymin=-10.1
	 elseif(nom_elem(i) .eq. 'O16')then
	  n_nu(i)=16
c	  ymax=0.1
c	  ymin=-6.1
	 elseif(nom_elem(i) .eq. 'O17')then
	  n_nu(i)=17
c	  ymax=-0.9
c	  ymin=-7.1
	 endif
c	 xxmax=-1.d30
c	 do j=1,itot
c	  xxmax=max(xxmax,p4(j))
c	 enddo
	 do j=1,itot
c	  p4(j)=log10(abs(compchim(nbelem*(j-1)+i)))
	  p4(j)=compg(nbelem*(j-1)+i)
	 enddo	!j
	 elem(i)=nom_elem(i)
	 call pminmax(p4,itot,ymax,ymin)
	 call pgenv(xmin,xmax,ymin,ymax,0,0)
	 call pglabel(en_x,elem(i),model)
	 call pgline(itot,point,p4)
	enddo
 
	write(6,*)'conservation et rapports ? entrer o/n'
	read(5,'(a)')oui
	if(oui .ne. 'o')goto100
 
c	write(6,2000)(n_nu(i),i=1,nbelem)
 
c	conservation totale
 
	do i=1,itot	
	 sum=0
	 do k=1,nbelem
	  sum=sum+compchim(nbelem*(i-1)+k)*n_nu(k)
	 enddo
	 p8(i)=sum	!somme des nucleons couche i
	enddo
	sum=0
	do i=1,itot	!somme totale des nucleons
	 sum=sum+p8(i)
	enddo
	sum=sum/float(itot)	!valeur moyenne
	do i=1,itot
	 p4(i)=p8(i)-sum	!ecart a la valeur moyenne couche i
c	 t4(i)=p4(i)
	enddo
	call pminmax(p4,itot,ymax,ymin)
	call pgenv(xmin,xmax,ymin,ymax,0,0)
	call pglabel(en_x,'conservation totale',model)
	call pgline(itot,point,p4)
 
 
c	le CNO
 
	do j=1,nbelem
	 if(nom_elem(j) .eq. 'C12')then
	  j12=j
	  cno=.true.
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
 
	if(.not. cno)goto100
 
c	conservation CNO
 
	do i=1,itot
	 do j=1,nbelem
	  xchim(j)=compchim(nbelem*(i-1)+j)
	 enddo
	 sum=xchim(j12)+xchim(j13)+xchim(j14)+xchim(j16)+xchim(j17)
	 if(j15 .gt. 0)sum=xchim(j15)+sum
	 p8(i)=sum
	enddo
	sum=0
	do i=1,itot
	 sum=sum+p8(i)
	enddo
	sum=sum/float(itot)
	do i=1,itot
	 p4(i)=p8(i)-sum
c	 write(6,2000)t4(i),p4(i),point(i)
	enddo
	call pminmax(p4,itot,ymax,ymin)
	call pgenv(xmin,xmax,ymin,ymax,0,0)
	call pglabel(en_x,'conservation CNO',model)
	call pgline(itot,point,p4)
 
 
c	C12/C13
 
	do i=1,itot	
	 do j=1,nbelem
	  xchim(j)=compg(nbelem*(i-1)+j)
	 enddo
	 p4(i)=log10(xchim(j12)/xchim(j13))
	enddo
	call pminmax(p4,itot,ymax,ymin)
	call pgenv(xmin,xmax,ymin,ymax,0,0)
	call pglabel(en_x,
     1	'Log10 '//nom_elem(j12)//'/'//nom_elem(j13),model)
	call pgline(itot,point,p4)
 
 
c	C13/N14
 
	do i=1,itot	
	 do j=1,nbelem
	  xchim(j)=compg(nbelem*(i-1)+j)
	 enddo
	 p4(i)=log10(xchim(j13)/xchim(j14))
	enddo
	call pminmax(p4,itot,ymax,ymin)
	call pgenv(xmin,xmax,ymin,ymax,0,0)
	call pglabel(en_x,
     1	'Log10 '//nom_elem(j13)//'/'//nom_elem(j14),model)
	call pgline(itot,point,p4)
 
 
c	N14/N15 puis N15/C12  ou N14/C12
 
	if(j15 .gt. 0)then
	 do i=1,itot	
	  do j=1,nbelem
	   xchim(j)=compg(nbelem*(i-1)+j)
	  enddo
	  p4(i)=log10(xchim(j14)/xchim(j15))
	 enddo
	 call pminmax(p4,itot,ymax,ymin)
	 call pgenv(xmin,xmax,ymin,ymax,0,0)
	 call pglabel(en_x,
     1	'Log10 '//nom_elem(j14)//'/'//nom_elem(j15),model)
	 call pgline(itot,point,p4)
	 do i=1,itot		!N15/C12
	  do j=1,nbelem
	   xchim(j)=compg(nbelem*(i-1)+j)
	  enddo
	  p4(i)=log10(xchim(j15)/xchim(j12))
	 enddo
	 call pminmax(p4,itot,ymax,ymin)
	 call pgenv(xmin,xmax,ymin,ymax,0,0)
	 call pglabel(en_x,
     1	'Log10 '//nom_elem(j15)//'/'//nom_elem(j12),model)
	 call pgline(itot,point,p4)
	 do i=1,itot		!N15/O16
	  do j=1,nbelem
	   xchim(j)=compg(nbelem*(i-1)+j)
	  enddo
	  p4(i)=log10(xchim(j15)/xchim(j16))
	 enddo
	 call pminmax(p4,itot,ymax,ymin)
	 call pgenv(xmin,xmax,ymin,ymax,0,0)
	 call pglabel(en_x,
     1	'Log10 '//nom_elem(j15)//'/'//nom_elem(j16),model)
	 call pgline(itot,point,p4)
	else
	 do i=1,itot		!N14/C12
	  do j=1,nbelem
	   xchim(j)=compg(nbelem*(i-1)+j)
	  enddo
	  p4(i)=log10(xchim(j14)/xchim(j12))
	 enddo
	 call pminmax(p4,itot,ymax,ymin)
	 call pgenv(xmin,xmax,ymin,ymax,0,0)
	 call pglabel(en_x,
     1	'Log10 '//nom_elem(j14)//'/'//nom_elem(j12),model)
	 call pgline(itot,point,p4)
	 do i=1,itot		!N14/O16
	  do j=1,nbelem
	   xchim(j)=compg(nbelem*(i-1)+j)
	  enddo
	  p4(i)=log10(xchim(j14)/xchim(j16))
	 enddo
	 call pminmax(p4,itot,ymax,ymin)
	 call pgenv(xmin,xmax,ymin,ymax,0,0)
	 call pglabel(en_x,
     1	'Log10 '//nom_elem(j14)//'/'//nom_elem(j16),model)
	 call pgline(itot,point,p4)
	endif
 
c	O16/O17
 
	do i=1,itot	
	 do j=1,nbelem
	  xchim(j)=compg(nbelem*(i-1)+j)
	 enddo
	 p4(i)=log10(xchim(j16)/xchim(j17))
	enddo
	call pminmax(p4,itot,ymax,ymin)
	call pgenv(xmin,xmax,ymin,ymax,0,0)
	call pglabel(en_x,
     1	'Log10 '//nom_elem(j16)//'/'//nom_elem(j17),model)
	call pgline(itot,point,p4)
 
 
c	O17/N14
 
	do i=1,itot	
	 do j=1,nbelem
	  xchim(j)=compg(nbelem*(i-1)+j)
	 enddo
	 p4(i)=log10(xchim(j17)/xchim(j14))
	enddo
	call pminmax(p4,itot,ymax,ymin)
	call pgenv(xmin,xmax,ymin,ymax,0,0)
	call pglabel(en_x,
     1	'Log10 '//nom_elem(j17)//'/'//nom_elem(j14),model)
	call pgline(itot,point,p4)
 
100	continue
	call pgend
 
	stop
 
	end
 
