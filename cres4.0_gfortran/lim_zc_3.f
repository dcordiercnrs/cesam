c********************************************************************
 
	subroutine lim_zc_3(bp,q,qt,knot,n,jlim,lim,lconv,derxx,lderxx,
     1	xx,xl,new_n,ni,nl,fac,
     2	mc,mct,nc,knotc,chim,m_zc,mstar,r2,m23,dim,
     3	r_zc,r_ov,bloc,etat,opa,conv,nuc)

c--------------------------------------------------------------------
c Modifications Daniel Cordier :

c       * Initialisation correcte des distance 'r_ov(i)', etc ...
c         ceci afin d'eviter des overshoot et undershoot delirant
c         apres une reprise.

c       * On met de l'OVERshoot uniquement si "i=1", pour eviter
c         d'en avoir sur les couches convectives qui apparaissent
c         quelques fois pres du coeur lors de la combustion de He

c       * estimation de la taille du coeur convectif (lorsqu'il y en
c         a un!) avec le critere de Roxburgh (cf. routine 'int_roxburgh')

c       * releves des temperatures en bord de zone convective : avec
c         et sans over/undershoot
c         Mot cle : 'tzc' (Fev. 1999)

c       * possibilite de "bloquer" (ou du moins de donner une borne sup.)
c         a la Zone Convective centrale. (Septembre 1999)

c-------------------------------------------------------------------

c	repartition des couches et determination des limites ZR/ZC.
c	Le nombre est modifie si new_n .ne. n
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
 
c	17 09 96 : correction bug dans declaration de vt; signale par D. Cordier
c	08 10 96 modif de la gestion de somme dm/somme r2 dm dans ZC
 
c	CESAM3
 
c entree
c	new_n : nouveau nombre de couches
c	mc,mct,nc,knotc,chim : pour interpolation
c	mstar: masse au temps t+dt, avec perte de masse
 
c entree/sortie
c	bp : solution en B-spline
c	q, qt, knot : abscisses, points de table, nombre (1, 2..n ou masse)
c	n : nombre de couches
c	fac : constantes de repartition
c	bloc: longueur du bloc du Jacobien
c	nl: nombre de lignes du Jacobien
c	dim : dimension de la base
c	r2, m23 : r**2, m**2/3
 
c sorties
c	jlim(i) : plus proche indice de la i-ieme limite ZR/ZC
c	lim : nombre de limites ZR/ZC
c	lconv(i)=.true. : dedut d'une ZC a la i-ieme limite ZR/ZC
c	derxx, lderxx : coefficients pour les splines
c	xx, xl : points de collocation et limites
c	ni(i) : nombre de limite en i
c	nl : nombre de lignes du Jacobien
c	m_zc : masses des limites ZR/ZC eventuellement overshootees
c	r_zc : rayons des limites ZR/ZC sans overshoot
c	r_ov : rayons des limites ZR/ZC eventuellement overshootees
c	fac(i) : poids de la couche [i, i+1]
 
c external
c	etat,opa,conv,nuc
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'ctephy.common'
	include 'modele_3.common'
	include 'evol_chim_3.common'
 
	integer knot,n,jlim(1),lim,nc,knotc,new_knot,indpc(pqt),
     1	new_n,i,j,ni(1),bloc,nl,lq,dim,lc,compt, iii

	integer idebut
c tzc
	real*8 logtzc(2*pnzc), logrozc(2*pnzc), logtov(2*pnzc)

	real*8 bp(1),q(1),qt(1),derxx(1),lderxx(1),xx(1),xl(1),
     1	mc(1),mct(1),m_zc(1),j1,j2,j3,dgrad,
     2	new(pn),chim(1),newt(pqt),f(pne),dfdq(pne),r2(1),m23(1),
     3	s1(pbp),vt(pbp*pm_qs),esp(pn),fac(1),r_zc(1),r_ov(1),
     4	p,t,m,l,r,d_grad(pn),a,j_new,gradconv,w,
     5	grp,gr,rap(2*pnzc),psi,rn,
     6  r_zc1, m_zc1, Hp1
c	data a/.45d0/				!limites ZR/ZC a 5%
 
	real*8 compx(pnelem),ro,drop,drot,drox,u,dup,dut,dux,mstar,
     1	gradient,dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,m13,
     2	epsilon(5),depsp,depst,depsx(pnelem),kap,dkapp,dkapt,dkapx,
     3	delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,bid,
     4	hp,gradad,gradrad,nh1,nhe1,nhe2,kip,lamb,
     5	dcompx(pnelem)

	real*8 r_rox, m_rox, diff, age, r_ov_max1, mroxSmsch

	real*8 mfix, mmstar(pn)


      common/roxburgh/diff, age, mroxSmsch
c tzc
	common /temp_limZCZR/ logtzc, logrozc, logtov

	external etat,opa,conv,nuc
	
	logical lconv(1),initr,ovsht,tot_conv,init,tot_rad, melange,
     +          rox_test

	logical blocage_core, limmfix, calib1, calib2

	data init/.true./
	data a/.45d0/				!limites ZR/ZC a 5%
 
	save
 
2000	format((1x,1p8d10.3))
 
	if(init)then
	 init=.false.

c        Causerie lors du premier appel de 'lmi_zc_3'

	call  causerie1_limzc( r2, n, lim, r_zc, m_zc, r_ov, 
     +              blocage_core, r_ov_max1, melange, limmfix,
     +              mfix, rox_test, calib1, calib2 )

c--------definition des points de collocation-------------	 	
 
c	 pause'lim_zc_3'
 
	 xl(1)=q(1)		!r au centre
	 xl(2)=q(1)		!l au centre
	 xl(3)=q(1)		!m au centre
	 xl(4)=q(n)		!P a l'exterieur
	 xl(5)=q(n)		!T a l'exterieur
	 xl(6)=q(n)		!m a l'exterieur
 
	 call colloc(m_qs,1,n,q,qt,knot,xx,vt,derxx,ne,lderxx,xl)
c	 pause'appel a colloc'
 
c	 localisation des limites dans la grille de points de raccord
c	 ni(i) : nombre de limite en i
 
	 do i=1,n
	  ni(i)=0
	 enddo	!i
	 do i=1,n-1
	  do j=1,1*ne
	   if(q(i) .le. xl(j) .and. q(i+1) .gt. xl(j))ni(i)=ni(i)+1
	  enddo	!j
	 enddo		!i
	 do j=1,1*ne
	  if(q(n) .eq. xl(j))ni(n-1)=ni(n-1)+1
	 enddo	!j
c	 write(6,*)'ni',ni(1),ni(2),ni(n-1),ni(n)
 
	 do i=1,n		!extraction de r2, m23 pour inter_3
	  r2(i) =bp(ne*(i-1)*m_qs+3)
	  m23(i)=bp(ne*(i-1)*m_qs+5)
	 enddo		!i
	endif		!init
	
c-----------modification du nombre de couches--------------------------

c	if ( .false. ) then ! TEST TEST TEST TEST
	
	if(n .ne. new_n)then

c	 on determine les new_n points de masse new pour
c	 assurer une repartition approximativement uniforme en
c	 lnP/(delta lnP)+lnT/(delta lnT)+(r/(delta r)**2+l/(delta l)+m**2/3))
 
	 do i=1,n		!fonction de repartition
	  call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,q(i),lq,f,dfdq)
c	  esp(i)=ctep*f(1)+ctem*f(5)
	  esp(i)=ctep*f(1)+cter*f(3)+ctel*f(4)+ctem*f(5) ! TEST
c	  print*, ' '
c	  print*, 'f(1)= ', f(1) ! ln P
c	  print*, 'f(5)= ', f(5) ! m^2/3
c	  print*, 'f(3)= ', f(3) ! r**2
c	  write(6,2000)f,q(i),esp(i)
	  if(i .gt. 1)then
	   do while(esp(i) .lt. esp(i-1))
	    esp(i)=esp(i)+0.01
	   enddo
	  endif
	 enddo	!i
c	 pause
c	 print*,'redistribution'
c	 pause'repartit1'
	 if(new_n .lt. 10)then
	  write(6,*)'on utilise 10 couches bien que ncouche=',new_n
	  new_n=10
	 endif
 
c	 on disposee les ncouche pour assurer une repartition
c	 approximativement uniforme de la fonction de repartition
 
	 call zone_3(n,q,esp,new_n,new)	!choix des nouveaux q
	
c	 write(6,*)'anciens q',n
c	 write(6,2000)(q(i),i=1,n)
c	 write(6,*)'esp'
c	 write(6,2000)(esp(i),i=1,n)
c	 write(6,*)'nouveaux q'
c	 write(6,2000)(new(i),i=1,new_n)
c	 pause'apres zone'
	
c	 on se place aux new_n points new, on met la spline dans s1
 
	 do i=1,new_n
	  call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,new(i),lq,f,dfdq)
	  do j=1,ne-1		!sauf pour espacement
	   s1(ne*(i-1)+j)=f(j)
	  enddo
c	  esp(i)=ctep*f(1)+ctem*f(5)
	  esp(i)=ctep*f(1)+cter*f(3)+ctel*f(4)+ctem*f(5) ! TEST	
	  if(i .gt. 1)then
	   do while(esp(i) .lt. esp(i-1))
	    esp(i)=esp(i)+0.01
	   enddo
	  endif
c	  write(6,2000)(f(j),j=1,ne-1),esp(i)
	 enddo
c	 pause'nouvel esp'
 
c	 s1 dans la base de snoein d'ordre mpr sur new_n couches
 
	 do i=1,new_n
	  new(i)=i
	  q(i)=i		!de facon previsionnelle
	 enddo
	
c	 pour la fonction d'espacement aux nouveaux points
 
	 call sbsp1dn(1,esp,new,newt,new_n,mpr,new_knot,.false.,
     1	new(1),lq,f,dfdq)
	
c	 derivee de la fonction d'espacement
 
	 do i=1,new_n
	  call sbsp1dn(1,esp,new,newt,new_n,mpr,new_knot,.true.,
     1	new(i),lq,f,dfdq)
	  s1(ne*(i-1)+ne)=dfdq(1)
	 enddo
 
c	 la nouvelle spline dans la base de snoein	
 
	 call sbsp1dn(ne,s1,new,newt,new_n,mpr,new_knot,.false.,
     1	new(1),lq,f,dfdq)
	
c	 do i=1,new_n
c	  call sbsp1dn(ne,s1,new,newt,new_n,mpr,new_knot,.true.,
c	1	new(i),lq,f,dfdq)
c	  write(6,2000)(f(j),j=1,ne)
c	 enddo
c	 pause'les s1 verification'
 
c	 base pour la collocation aux nouveaux points
 
	 call noedif(q,qt,new_n,m_qs,1,knot)
	
c	 changement de base
	
	 call newspn(ne,new,newt,new_knot,mpr,qt,knot,m_qs+1,s1,bp,vt,indpc)
	
c	 nouvelle distribution
 
	 n=new_n
	 dim=(n-1)*m_qs+1	!dimension
	 nl=ne*dim		!nombre de lignes
	 bloc=ne*mpr		!longueur d'un bloc
	 do i=3,5		!evite des erreurs d'arrondi
	  bp(i)=0.d0
	 enddo
c	 do i=1,n
c	  call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,q(i),lq,f,dfdq)
c	  write(6,2000)(f(j),j=1,ne)
c	 enddo
c	 print*,ne,knot,mpr
c	 pause'nouvelle base'
	
c--------redefinition des points de collocation-------------	 	
 
c	 pause'lim_zc_3, 2'
 
	 xl(1)=q(1)		!r au centre
	 xl(2)=q(1)		!l au centre
	 xl(3)=q(1)		!m au centre
	 xl(4)=q(n)		!P a l'exterieur
	 xl(5)=q(n)		!T a l'exterieur
	 xl(6)=q(n)		!m a l'exterieur
 
	 call colloc(m_qs,1,n,q,qt,knot,xx,vt,derxx,ne,lderxx,xl)
c	 pause'appel a colloc'
 
c	 localisation des limites dans la grille de points de raccord
c	 ni(i) : nombre de limite en i
 
	 do i=1,n
	  ni(i)=0
	 enddo	!i
	 do i=1,n-1
	  do j=1,1*ne
	   if(q(i) .le. xl(j) .and. q(i+1) .gt. xl(j))ni(i)=ni(i)+1
	  enddo	!j
	 enddo		!i
	 do j=1,1*ne
	  if(q(n) .eq. xl(j))ni(n-1)=ni(n-1)+1
	 enddo	!j
c	 write(6,*)'ni',ni(1),ni(2),ni(n-1),ni(n)
 
	 do i=1,n		!extraction de r2, m23 pour inter_3
	  r2(i) =bp(ne*(i-1)*m_qs+3)
	  m23(i)=bp(ne*(i-1)*m_qs+5)
	 enddo		!i
	endif		!modification du nombre de couches


c	end if ! if .false. TEST TEST TEST
 
c------------------------limites ZR/ZC----------------------------
 
c	recherche les limites zones convectives / zones radiatives
c	en fonction de q, determination des facteurs de repartition
c	overshooting inferieur et superieur : extension de la zone melangee
c	ou overshooting inferieur PC de JPZ avec saut du gradient
 
c	modele totalement convectif: lim=1,jlim(1)=n,lconv(1)=.false.
c	modele totalement radiatif: lim=0,jlim(i)=-100,lconv(i)=.false.
 
c	pause'debut recherche limites ZR/ZC'
 
	write(6,*)' '
	write(6,*)'------------- La Convection (debut) -----------------'
	write(6,*)' '
	
	tot_conv=.false.
	tot_rad =.false.
	ovsht=max(ovshts,ovshti) .gt. 0.d0
 
c	calcul de grad-grad
 
c	write(6,*)nc,m_ch,knotc,n,nbelem
c	write(6,2000)(mc(j),j=1,nc)
c	write(6,2000)(mct(j),j=1,knotc)
c	print*,mstar
c	pause'debut lim_zc'
 
	lim=0		!initialisation	
	do i=1,n		!calcul de grad-grad
	 call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,q(i),lq,f,dfdq)
c	 write(6,2000)f
	 p=exp(f(1))
	 t=exp(f(2))
	 r=sqrt(abs(f(3)))
	 l=sqrt(abs(f(4)))**3
	 m=sqrt(abs(f(5)))**3
	 mmstar(i)=m
	 m13=max(m**(1.d0/3.d0),1.d-5)	!en 0 decalage pour calcul derivee
	 call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,min(f(5),mc(nc)),lc,compx,dcompx)
c	 write(6,2000)p,t,r,l,m,compx(1),q(i),mc(nc),f(5)
c	 write(6,*)mc(1),mc(nc),f(5),compx(1)
c	 print*
c	 pause'boucle lim_zc'
	 do j=1,nbelem		!dX/dm
	  dcompx(j)=dcompx(j)*2.d0/3.d0/m13
	 enddo
	
	 if(iw .gt. 1 .and. f(3) .ne. 0.)then		!rotation
	  w=compx(iw)/f(3)
	 else
	  w=w_rot
	 endif
	
	 d_grad(i)=dgrad(p,t,compx,m,l,r,dcompx,mstar,w,etat,opa,nuc)
	enddo	!i
	
	rn=r		!r a l'exterieur
c	pause'lim_zc'

c----------------------------------------------------------------------------
        write(6,*)' ' 
	write(6,*)'********************************************************* '
	write(6,*)' '
	write(6,*)'* Rayon de l''etoile ''rn''= ', rn            
	write(6,*)' '
	write(6,*)'********************************************************* '
	write(6,*)' '
c----------------------------------------------------------------------------

	grp=d_grad(1)
	if(grp .eq. 0. .and. d_grad(2) .ge. 0.)grp=-1  !ZC avec zero exact en 1
 
	do j=2,n			!recherche des limites ZR/ZC
	 gr=d_grad(j)
	 if(grp .ne. 0. .and. gr*grp .le. 0.)then
	  lim=lim+1
	  if(lim .gt. 2*pnzc)then
	   write(6,*)'nombre de limites ZR/ZC superieur a',pnzc*2
	   write(6,*)'limites ZR/ZC imprecises, floues ?'
	   write(6,*)'difficulte d''initialisation ? erreur ? '
	   write(6,*)'sinon augmenter le parametre pnzc et recompiler'
	   write(6,*)'poursuite considerant le modele totalement radiatif:'
	   write(6,*)'pour la composition chimique i.e. en ne melangeant'
	   write(6,*)'pas les ZC pour ce modele'
	   tot_rad=.true.
	   lim=0
	   do i=1,pnzc
	    jlim(i)=-100
	    lconv(i)=.false.
	   enddo
	   goto10
	  endif
	  lconv(lim)=grp .lt. 0.	!.true. passage ZR-->ZC
	  rap(lim)=abs(grp/gr)
	  rap(lim)=rap(lim)/(1.d0+rap(lim))
	  jlim(lim)=j				!limite ZR/ZC entre j et j-1
	  if(lconv(lim))then
	   write(6,3)int(rap(lim)*100.),j-1,j,lim
3	   format(1x,'debut de ZC localisee a',i3,'% entre les couches #',
     1	i3,' et #',i3,' limite ZR/ZC #',i3)
	  else
	   write(6,4)int(rap(lim)*100.),j-1,j,lim
4	   format(1x,'fin de ZC localisee a',i3,'% entre les couches #',
     1	i3,' et #',i3,' limite ZR/ZC #',i3)
	  endif
	 endif
	 grp=gr
	enddo			!j

c	write(6,*)'sortie du tri, lim',(lconv(i),i=1,lim),lim
c	pause
 
10	write(6,*)' '
	do i=1,n
	 new(i)=q(i)		!indices des points de grille
	 fac(i)=1.d0		!facteurs de repartittion
	enddo
c========================================================================================
c Option de blocage de la limite (en MASSE) du coeur convectif)
	if ( limmfix ) then
c	   if ( d_grad(3) .gt. 0. ) then ! Y-a-t-il un coeur convectif
	   jlim(1)=0
	      do i= 1, n-1
		 if ( ( mfix .ge. mmstar(i)) .AND. 
     +                ( mfix .lt. mmstar(i+1)) ) then
		    if ( abs(mfix-mmstar(i)) .gt. abs(mfix-mmstar(i-1)) ) then
		       jlim(1)=i+1
		    else
		       jlim(1)=i
		    end if
		 end if
	      end do
	      m_zc(1)=mfix
	      print*, ' '
	      print*, 'jlim(1)= ', jlim(1)
c	      print*, 'lconv(1)= ', (lconv(i),i=1,lim)
	      pause
c	   end if
	end if
c=======================================================================================
c	il y a une limite ZR/ZC entre jlim(i)-1 et jlim(i), i=1,lim
c	on determine l'incide new(jlim(i)) a la limite ZR/ZC
c	on fixe la limite soit a doite soit a gauche
c	jlim(i) : indice le plus proche de la limite
 
	do i=1,lim
	 j_new=new(jlim(i)-1)+
     1	rap(i)*(new(jlim(i))-new(jlim(i)-1))	!indice de la lim. ZR/ZC
c	 print*,i,j_new,new(jlim(i)-1),rap(i),new(jlim(i))
	 jlim(i)=nint(j_new)	
	 new(jlim(i))=j_new		!nouvel indice pour les limites ZR/ZC
c	 print*,i,j_new,jlim(i)
	enddo
 
c	print*
c	print*,lim
c	print*,(new(jlim(i)),i=1,lim)
c	print*,(jlim(i),i=1,lim)	
 
c	suppression des limites trop proches
 
	i=1
	do while (i .le. lim)
	 if(jlim(i) .le. 3)then		!limites sur premieres couches
	  write(6,*)'suppression d''une limite ZR/ZC trop centrale #',i
	  lim=lim-1				!supprimees
	  do j=1,lim
	   jlim(j)=jlim(j+1)
	   new(jlim(j))=new(jlim(j+1))
	   lconv(j)=lconv(j+1)
	  enddo
	 else
	  i=i+1
	 endif
	enddo
 
c	print*
c	print*,lim
c	print*,(new(jlim(i)),i=1,lim)
c	print*,(jlim(i),i=1,lim)	
 
	i=2
	do while (i .le. lim)
	 if(jlim(i)-jlim(i-1) .le. 4)then	!limites trop proches
	  write(6,*)'limites trop proches, suppression des limites #',i-1,i
	  lim=lim-2			!suppression de i-1 et i
	  do j=i-1,lim
	   jlim(j)=jlim(j+2)
	   new(jlim(j))=new(jlim(j+2))
	   lconv(j)=lconv(j+2)
	  enddo
	 else
	  i=i+1
	 endif			!pour la suppression de limites trop rapprochees
	enddo
 
c	print*
c	print*,lim
c	print*,(new(jlim(i)),i=1,lim)
c	print*,(jlim(i),i=1,lim)	
 
	i=1
	do while (i .le. lim)
 	 if(jlim(i) .ge. n-4)then
	  write(6,*)'suppression d''une limite ZR/ZC trop externe #',lim
	  new(jlim(i))=jlim(i)
	  lim=lim-1
	 else
	  i=i+1
	 endif
	enddo
 
c	print*
c	print*,lim
c	print*,(new(jlim(i)),i=1,lim)
c	print*,(jlim(i),i=1,lim)	
 
	if(lim .eq. 0 .and. jlim(1) .ne. -100)then
	 if(d_grad(3) .ge. 0.)then
	  write(6,*)'modele completement convectif'
	  tot_conv=.true.
	  lim=1
	  jlim(1)=n
	  lconv(1)=.false.
	  m_zc(1)=mstar
	  r_zc(1)=rn
	 else
	  write(6,*)'modele completement radiatif'
	  tot_rad=.true.
	  lim=0
	  do i=1,pnzc
	   jlim(i)=-100
	   lconv(i)=.false.
	  enddo
	  m_zc(1)=-100
	 endif
	 do i=1,n
	  fac(i)=1.d0
	 enddo
	
	else			!pas de limites ZR/ZC
	
c-----------on affine, par dichotomie, la position des limites retenues	
c==========================================================
c Quelques petites modifs concernant l'option de blocage
c (en MASSE) du coeur convectif)
c	   print*, ' '
c	   print*, 'limmfix= ', limmfix
c	   pause
	   if ( limmfix ) then
	      idebut=lim+1
	   else
	      idebut=1
	   end if
	   print*, 'idebut= ', idebut
	   print*, 'lim= ', lim

c	 do i= idebut, lim
	 do while ( .false. )
c        =================
	  j1=new(jlim(i))
	  call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,j1,lq,f,dfdq)
c	  write(6,2000)f
	  p=exp(f(1))
	  t=exp(f(2))
	  r=sqrt(abs(f(3)))
	  l=sqrt(abs(f(4)))**3
	  m=sqrt(abs(f(5)))**3
	  m13=max(m**(1.d0/3.d0),1.d-5)	!en 0 decalage pour calcul derivee
	  call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,min(f(5),mc(nc)),lc,compx,dcompx)
	  do j=1,nbelem		!dX/dm
	   dcompx(j)=dcompx(j)*2.d0/3.d0/m13
	  enddo
	
	  if(iw .gt. 1 .and. f(3) .ne. 0.)then		!rotation
	   w=compx(iw)/f(3)
	  else
	   w=w_rot
	  endif
	
	  d_grad(1)=dgrad(p,t,compx,m,l,r,dcompx,mstar,w,etat,opa,nuc)
	
	  j2=jlim(i)
	  call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,j2,lq,f,dfdq)
c	  write(6,2000)f
	  p=exp(f(1))
	  t=exp(f(2))
	  r=sqrt(abs(f(3)))
	  l=sqrt(abs(f(4)))**3
	  m=sqrt(abs(f(5)))**3
	  m13=max(m**(1.d0/3.d0),1.d-5)	!en 0 decalage pour calcul derivee
	  call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,min(f(5),mc(nc)),lc,compx,dcompx)
	  do j=1,nbelem		!dX/dm
	   dcompx(j)=dcompx(j)*2.d0/3.d0/m13
	  enddo
	
	  if(iw .gt. 1 .and. f(3) .ne. 0.)then		!rotation
	   w=compx(iw)/f(3)
	  else
	   w=w_rot
	  endif
 
	  d_grad(2)=dgrad(p,t,compx,m,l,r,dcompx,mstar,w,etat,opa,nuc)	
	
	  if(d_grad(1)*d_grad(2) .gt. 0)then
	   j2=jlim(i)-1
	   call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,j2,lq,f,dfdq)
c	   write(6,2000)f
	   p=exp(f(1))
	   t=exp(f(2))
	   r=sqrt(abs(f(3)))
	   l=sqrt(abs(f(4)))**3
	   m=sqrt(abs(f(5)))**3
	   m13=max(m**(1.d0/3.d0),1.d-5) !en 0 decalage pour calcul derivee
	   call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,min(f(5),mc(nc)),lc,compx,dcompx)
	   do j=1,nbelem		!dX/dm
	    dcompx(j)=dcompx(j)*2.d0/3.d0/m13
	   enddo
	
	   if(iw .gt. 1 .and. f(3) .ne. 0.)then		!rotation
	    w=compx(iw)/f(3)
	   else
	    w=w_rot
	   endif
 
	   d_grad(2)=dgrad(p,t,compx,m,l,r,dcompx,mstar,w,etat,opa,nuc)
	
	   if(d_grad(1)*d_grad(2) .gt. 0)then	
	    j2=jlim(i)+1
	    call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,j2,lq,f,dfdq)
c	    write(6,2000)f
	    p=exp(f(1))
	    t=exp(f(2))
	    r=sqrt(abs(f(3)))
	    l=sqrt(abs(f(4)))**3
	    m=sqrt(abs(f(5)))**3
	    m13=max(m**(1.d0/3.d0),1.d-5) !en 0 decalage pour calcul derivee
	    call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,min(f(5),mc(nc)),lc,compx,dcompx)
	    do j=1,nbelem		!dX/dm
	     dcompx(j)=dcompx(j)*2.d0/3.d0/m13
	    enddo
	
	    if(iw .gt. 1 .and. f(3) .ne. 0.)then		!rotation
	     w=compx(iw)/f(3)
	    else
	     w=w_rot
	    endif
	
	    d_grad(2)=dgrad(p,t,compx,m,l,r,dcompx,mstar,w,etat,opa,nuc)
	    	
	    if(d_grad(1)*d_grad(2) .gt. 0)then
	     print*,'lim_zc_3 : erreur dans la dichotomie, limite ZR/ZC',i
	     print*,'j1, jlim(i)',j1,jlim(i)
	     !pause'on arrete'
	     print*,'on arrete'
	     stop
	    endif	!jlim(i)+1
	   endif	!jlim(i)-1
	  endif	!on a encadre la limite entre j1 et j2
	  compt=0	
	  j3=(j1+j2)/2.d0	
	  do while(abs(j1-j2) .gt. 1.d-4)
	   call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,j3,lq,f,dfdq)
c	   write(6,2000)f
	   p=exp(f(1))
	   t=exp(f(2))
	   r=sqrt(abs(f(3)))
	   l=sqrt(abs(f(4)))**3
	   m=sqrt(abs(f(5)))**3
	   m13=max(m**(1.d0/3.d0),1.d-5) !en 0 decalage pour calcul derivee
	   call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,min(f(5),mc(nc)),lc,compx,dcompx)
	   do j=1,nbelem		!dX/dm
	    dcompx(j)=dcompx(j)*2.d0/3.d0/m13
	   enddo
 
	   if(iw .gt. 1 .and. f(3) .ne. 0.)then		!rotation
	    w=compx(iw)/f(3)
	   else
	    w=w_rot
	   endif
	
	   d_grad(3)=dgrad(p,t,compx,m,l,r,dcompx,mstar,w,etat,opa,nuc)	
	
c	   print*,j1,j2,j3,compt
c	   write(6,2000)d_grad(1),d_grad(2),d_grad(3)
	
	   if(d_grad(2)*d_grad(3) .lt. 0.)then
	    j1=j3
	    d_grad(1)=d_grad(3)
	   elseif(d_grad(1)*d_grad(3) .le. 0.)then		
	    j2=j3
	    d_grad(2)=d_grad(3)
	   endif
	   j3=(j1+j2)/2.d0
	   compt=compt+1
	   if(compt .gt. 30)then
	    print*,'dichotomie dans lim_zc_3, compt>30',compt
	    print*,j1,d_grad(1),j2,d_grad(2)
	    !pause'abandon'
	    print*,'abandon'
	    stop
	   endif
	  enddo		!while	
	 enddo

c================================
c On met les deux lignes suivantes en commentaires DC 12/12/99
c	 jlim(i)=nint(j3)	
c	 new(jlim(i))=j3
	
c-------------------ecritures-------------------------------	
	
	 write(6,*)' '
	 write(6,*)'Limites ZR/ZC retenues apres affinement de la localisation:'
	 do i=1,lim		!masses des limites
	  call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,new(jlim(i)),lq,f,dfdq)
	  r_zc(i)=sqrt(abs(f(3)))
	  m_zc(i)=sqrt(abs(f(5)))**3
c tzc
	  logtzc(i)=log10(exp(f(2)))
	  logrozc(i)=log10( exp(f(1))/exp(f(2))*
     +               (nucleo(1)/2.+(nucleo(2)+nucleo(3))/2.)*amu/kbol ) 

c	  write(6,*)'masse r_zc,m_zc',i,jlim(i),new(jlim(i))
 
	  if(lconv(i))then	!debut de ZC
	   write(6,24)i,int((new(jlim(i))-int(new(jlim(i))))*100.),
     1	int(new(jlim(i))),int(new(jlim(i)))+1,r_zc(i)/rn,
     2	m_zc(i)/mstar,1.d0-m_zc(i)/mstar
24	   format(2x,'pour limite ZR/ZC #',i4,' a ',i3,'% entre couches',i4,
     1	' et',i4,', debut de ZC',/,3x,
     2	'R_zc/Rtot=',1pd10.3,', M_zc/Mstar=',1pd10.3,
     3	' 1-M_zc/Mstar=',1pd10.3)
	  else
	   write(6,25)i,int((new(jlim(i))-int(new(jlim(i))))*100.),
     1	int(new(jlim(i))),int(new(jlim(i)))+1,r_zc(i)/rn,
     2	m_zc(i)/mstar,1.d0-m_zc(i)/mstar	
25	   format(2x,'pour limite ZR/ZC #',i4,' a ',i3,'% entre couches #',i4,
     1	' et #',i4,', fin de ZC',/,3x,
     2	'R_zc/Rtot=',1pd10.3,', M_zc/Mstar=',1pd10.3,
     3	' 1-M_zc/Mstar=',1pd10.3)
	
	  endif
 	 enddo		!i

c	 write(6,*) 'jlim(i)= ', (jlim(i),i=1,lim)
c	 write(6,*) 'lconv(i)= ', (lconv(i),i=1,lim)
c	 write(6,*) 'm_zc(i)= ', (m_zc(i),i=1,lim)
c	 pause

c--------------------------------------------------------------------
c Critere de Roxburgh

	if ( Rox ) then
c On calcule r_zc1 et Hp1 (respectivement rayon classique du coeur et
c echelle de hauteur de pression a cet endroit.
	   call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,new(jlim(1)),lq,f,dfdq)
	   p=exp(f(1))
	   t=exp(f(2))
	   l=sqrt(abs(f(4)))**3
	   r_zc1=sqrt(abs(f(3)))
	   if ( r_zc1 .ne. r_zc(1) ) then
	      print*, 'Pb. r_zc(1) .ne. r_zc(1) !!'
	      print*, 'r_zc(1)= ', r_zc(1)
	      print*, 'r_zc1  = ', r_zc1
	   end if
	   m_zc1=sqrt(abs(f(5)))**3
	   m13=max(m_zc1**(1.d0/3.d0),1.d-5)
	   call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,min(f(5),mc(nc)),lc,compx,dcompx)
	   do j=1,nbelem		!dX/dm
	    dcompx(j)=dcompx(j)*2.d0/3.d0/m13
	   enddo
	
	   call thermo_3(p,t,compx,m_zc1,l,r_zc1,.false.,dcompx,mstar,
     1	ro,drop,drot,drox,u,dup,dut,dux,lamb,r_zc,r_ov,0,
     2	gradient,dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,
     3	epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,
     4	delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
     5	hp1,gradad,gradrad,gradconv,bid,w,0.d0,nh1,nhe1,nhe2,
     6	etat,opa,conv,nuc)

	   print*, 'Appel de ''int_roxburgh'''
	   if ( .NOT. melange ) then
	      call int_roxburgh( jlim, r_zc1, m_zc1, Hp1, n, bp, q, qt, knot, 
     +                   chim, mc, mct, nc, knotc, r_rox, m_rox, etat, 
     +                   opa, conv, rox_test, calib1, calib2 )
	   else
	      call int_roxburgh_m( jlim, r_zc1, m_zc1, Hp1, n, bp, q, qt, 
     +                knot, chim, mc, mct, nc, knotc, r_rox, etat, 
     +                opa, conv, rox_test, calib1, calib2 )
	   end if

	   print*, 'l''appel de Roxburgh s''est bien passe voici ''r_rox'':'
           print*, 'r_rox (en rsol)= ', r_rox
	   print*, 'r_zc(1)= ', r_zc(1)
	   print*, 'm_rox (en msol)= ', m_rox
	   if ( m_zc1 .ne. 0.d0 ) then
	      print*, 'm_rox/m_zc1= ', m_rox/m_zc1
	   end if
	      
	   if ( r_zc(1) .gt. 0.) then
	      if ( (r_rox-r_zc(1))/r_zc(1)*100. .le. -10. ) then
		 if ( r_rox .ne. -100.d0 ) then
		    print*, 'Rayon de Roxburgh tres petit !'
c		    pause
		 end if
	      end if
	   end if
c           print*, 'rn= ', rn
c           print*, 'r_ov(1)= ', r_ov(1)
c	   pause
	end if

	      
c___________  calcul de r_zc et de r_ov si overshooting ou PC _____________
 
	 if(ovsht)then		!s'il y a overshooting ou PC
	  initr=.false.		!pour interpolation de R en Q
	  do i=1,lim
c	   print*,'limite ',i,' ',lconv(i)
	   r_zc(i)=-100
	   r_ov(i)=-100
	   call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,new(jlim(i)),lq,f,dfdq)
	   p=exp(f(1))
	   t=exp(f(2))
	   l=sqrt(abs(f(4)))**3
	   r_zc(i)=sqrt(abs(f(3)))
	   m_zc(i)=sqrt(abs(f(5)))**3
	   m13=max(m_zc(i)**(1.d0/3.d0),1.d-5)	!en 0 decalage
	   call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,min(f(5),mc(nc)),lc,compx,dcompx)
	   do j=1,nbelem		!dX/dm
	    dcompx(j)=dcompx(j)*2.d0/3.d0/m13
	   enddo
c	   write(6,*)'m_zc,r_zc,new(jlim(i))',lconv(i),i
c	   write(6,2000)m_zc(i),r_zc(i),new(jlim(i))
c	   write(6,2000)p,t,r_zc(i),l,m_zc(i),compx(1)
 
	   if(iw .gt. 1 .and. f(3) .ne. 0.)then		!rotation
	    w=compx(iw)/f(3)
	   else
	    w=w_rot
	   endif
	
	   call thermo_3(p,t,compx,m_zc(i),l,r_zc(i),.false.,dcompx,mstar,
     1	ro,drop,drot,drox,u,dup,dut,dux,lamb,r_zc,r_ov,0,
     2	gradient,dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,
     3	epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,
     4	delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
     5	hp,gradad,gradrad,gradconv,bid,w,0.d0,nh1,nhe1,nhe2,
     6	etat,opa,conv,nuc) !fin de ZC, ext. vers le haut:rext>r_ov>r_zc
 
c	   overshoot: ovshti(s)*Hp ou ovshti*KIp de JPZ si T>5d5 a la limite
 
	   if(t .gt. 5.d5)then
	    if(lconv(i))then !debut de ZC, ext. vers le bas
	     if(ovshti .gt. 0.d0)then
	      if(jpz)then
	       kip=(3.-t*(drot/ro+dkapt/kap))*gradad-p*(drop/ro+dkapp/kap)
	       r_ov(i)=max(r_zc(i)-hp*ovshti/rsol/kip,0.d0)
c	       write(6,*)'i/r_zc(i),hp*ovshti/rsol/kip,hp,ovshti,rsol,kip',i
c	       write(6,2000)r_zc(i),hp*ovshti/rsol/kip,hp,ovshti,rsol,kip
c	       write(6,2000)r_ov(i)
	      else		!ext. vers le bas: 0.<r_ov<r_zc
	       r_ov(i)=max(r_zc(i)-hp*ovshti/rsol,0.d0)
	      endif	!jpz
c	      write(6,*)'overshoot inferieur'
c	      write(6,2000)r_ov(i),r_zc(i),hp
	     endif		!ovshti
	    else	!fin de ZC overshoot superieur: rext>r_ov>r_zc
	     if(ovshts .gt. 0.d0)then		!il y a ovsht sup.
	      if(jpz)then
	       r_ov(i)=r_zc(i)*(1.d0+ovshts) !overshoot sup.= R_noyau X ovshts
	      else
	       if ( i .eq. 1 ) then ! On met de l'overshoot uniquement sur le coeur !!
		  if ( Rox ) then ! On utilise le critere de Roxburgh
		     if ( r_rox .gt. 0.d0 ) then
              do iii= 1, lim
                 print*, 'r_zc(',iii,')= ', r_zc(iii)
              end do
			if ( r_rox .ge. r_zc(i) ) then
                  print*,'dif. en Hp= ',(r_rox-r_zc(1))/Hp*rsol
                  diff=(r_rox-r_zc(1))/Hp*rsol
c                  pause
			   r_ov(i)=r_rox
			else
			   print*, 'Pb. dans lim_zc_3 :'
			   print*, 'r_rox < r_zc(1) !'
			   print*, 'On fait : r_ov=r_zc !'
			   diff=0.d0
			   r_ov(i)=r_zc(i)
c			   stop
			end if ! if ( r_rox .ge. r_zc(i) )
                     else ! if ( r_rox .gt. 0.d0 )
                       print*, ' '
                       print*, 'TEST TEST TEST ******************'
                       print*, ' '
                       print*, 'r_rox est .le. 0. voici sa valeur :'
		       print*, '(pas de coeur convectif detecte dans ''int_roxburgh'')'
                       print*, 'r_rox=  ', r_rox
                       print*, 'On va faire r_ov(1)=-100. ; Ok ?'
c                       pause
                       r_ov(i)=-100.
		     end if
		  else ! if ( Rox )
		     r_ov(i)=min(r_zc(i)+hp*ovshts/rsol,r_zc(i)*(1.d0+ovshts))
		  end if ! if ( Rox )
	       end if ! if ( i .eq. 1 )
	      endif !jpz
	      r_ov(i)=min(r_ov(i),rn)
c	      write(6,*)'overshoot superieur'
c	      write(6,2000)r_ov(i),r_zc(i),hp,rn
	     endif		!ovshts
	    endif		!lconv
c	    pause'apres if sur lconv'
 
c	    masse d'overshooting ou de penetration

c------------------------------------------------------
c Blocage du Convective Core :
      if ( blocage_core )  then
	  call bloc_cc( r_ov, r_ov_max1 )
      end if
c------------------------------------------------------

	    if(r_ov(i) .ge. 0.)then	!r_ov=-100: no overshoot en r_zc(i)
	     call inter_3('r2',bp,q,qt,n,knot,r_ov(i)**2,f,dfdq,r2,m23)
	     m_zc(i)=sqrt(f(5))**3
c tzc
	     logtov(i)=log10(exp(f(2)))
	     j_new=f(ne)
c	     write(6,2000)(f(j),j=1,ne)
 
c	     l'ancien new(jlim(i)) est desormais un point ordinaire
c	     qui redevient l'indice entier jlim(i)
c	     la limite est desormais en jlim(i)= + proche entier de j_new
 
c	     write(6,*)'jlim(i),new(jlim(i)),j_new',jlim(i),new(jlim(i)),j_new
	     jlim(i)=nint(j_new)	
	     new(jlim(i))=j_new		!nouvel indice pour les limites ZR/ZC
c	     write(6,*)'jlim(i),new(jlim(i)),j_new',jlim(i),new(jlim(i)),j_new
c	     write(6,*)'j_new/M_ov,R_ov',j_new
c	     write(6,2000)m_zc(i),r_ov(i)
c	     pause'apres j_new'
	    endif		!r_ov .ne. 0
	   endif		!t > 1.d5
	  enddo		!i sur lim
 
c	  il ne doit pas y avoir de chevauchement apres les extensions
 
c	  write(6,*)lim,ovsht,(jlim(i),i=1,lim),(lconv(i),i=1,lim)
c	  write(6,2000)(r_ov(j),j=1,lim)
c	  write(6,2000)(r_zc(j),j=1,lim)
c	  write(6,2000)(m_zc(j),j=1,lim)
 
	  if(lconv(1))then	!debut de ZC a la premiere limite
	   i=2
	  else
	   i=1
	  endif
	  do while (i .lt. lim)
	   if(jlim(i) .ge. jlim(i+1))then
	    do j=i,lim-2	!chevauchement: suppression d'une ZR decalage
	     jlim(j)=jlim(j+2)
	     m_zc(j)=m_zc(j+2)
	     lconv(j)=lconv(j+2)
	     r_zc(j)=r_zc(j+2)
	     r_ov(j)=r_ov(j+2)
	    enddo
	    lim=lim-2
	   else
	    i=i+2
	   endif
	  enddo		!while
	
c	  a cause d'un overshoot, le modele devient totalement convectif
 
	  if(lim .eq. 1 .and. jlim(1) .eq. 1 .and. lconv(1))then
	   jlim(1)=n
	   lconv(1)=.false.
	   tot_conv=.true.
	   r_zc(1)=rn
	  else
	   tot_conv=.false.
	
c	   pas de limite ZR/ZC aux extremites
 
	   if(jlim(1) .eq. 1)then		!en 1
	    do i=1,lim-1
	     jlim(i)=jlim(i+1)
	     lconv(i)=lconv(i+1)
	     m_zc(i)=m_zc(i+1)
	     r_zc(i)=r_zc(i+1)
	     r_ov(i)=r_ov(i+1)
	    enddo
	    lim=lim-1
	   endif
	   if(jlim(lim) .eq. n)lim=lim-1			!en n
	  endif
	
	  write(6,*)' '
	  write(6,*)'limites ZR/ZC apres examen de l''overshooting:'
	  write(6,*)' '
	
c	  print*,lim,(jlim(i),i=1,lim)
	
	  if(tot_conv)then
	   write(6,33)
33	   format(1x,'modele totalement convectif',/)
	  else
	   do j=1,lim
	    if(lconv(j))then		!debut de ZC
	     if(ovshti .gt. 0.d0)then
	      if(jpz)then
	       write(6,31)j,jlim(j),r_zc(j)/rn,r_ov(j)/rn,m_zc(j)/mstar
31	       format(2x,'Penetration Convective JPZ limite ZC/ZR #',i4,
     1	', couche #',i4,/,3x,'rayon ZC/Rtot=',1pd10.3,
     2	', rayon reduit/Rtot=',1pd10.3,/,3x,
     3	'm/Mstar a la limite de la penetration convective=',1pd10.3,/)
	      else		!jpz
	       write(6,2)j,jlim(j),r_zc(j)/rn,r_ov(j)/rn,m_zc(j)/mstar
2	       format(2x,'Overshooting inferieur limite ZC/ZR #',i4,
     1	', couche #',i4,/,3x,'rayon ZC/Rtot=',1pd10.3,
     2	', rayon reduit/Rtot=',1pd10.3,/,3x,
     3	'm/Mstar a la limite de l''overshooting=',1pd10.3,/)
	      endif		!jpz
	     else		!sans overshoot
	      write(6,11)j,jlim(j),r_zc(j)/rn,m_zc(j)/mstar
11	      format(2x,'Limite sans overshooting ni PC, limite ZR/ZC #:',i4,
     1	', couche #',i4,/,3x,'rayon ZC/Rtot=',1pd10.3,
     3	' m/Mstar a la fin de la ZC=',1pd10.3,/)
	     endif		!ovshti
	    else		!fin de ZC
	     if(ovshts .gt. 0.d0)then
	      if(jpz)then
	       write(6,12)j,jlim(j),r_zc(j)/rn,r_ov(j)/rn,m_zc(j)/mstar
12	       format(2x,'Overshoot superieur * R noyau, limite ZC/ZR #',i4,
     1	', couche #',i4,/,3x,'rayon ZC/Rtot=',1pd10.3,
     2	', rayon etendu/Rtot=',1pd10.3,/,3x,
     3	'm/Mstar a la limite de l''overshoot=',1pd10.3,/)
	      else		!jpz
	       write(6,1)j,jlim(j),r_zc(j)/rn,r_ov(j)/rn,m_zc(j)/mstar
1	       format(2x,'Overshooting superieur * Hp, limite ZC/ZR #',i4,
     1	', couche #',i4,/,3x,'rayon ZC/Rtot=',1pd10.3,
     2	', rayon etendu/Rtot=',1pd10.3,/,3x,
     3	'm/Mstar a la limite de l''overshooting=',1pd10.3,/)
	      endif		!jpz
	     else		!pas d'overshoot
	      write(6,21)j,jlim(j),r_zc(j)/rn,m_zc(j)/mstar
21	      format(2x,'Limite sans overshooting ni PC, limite ZR/ZC #',i4,
     1	', couche #',i4,/,3x,'rayon ZC/Rtot=',1pd10.3,
     3	' m/Mstar au debut de la ZC=',1pd10.3,/)
	     endif		!ovshts
	    endif		!lconv
	   enddo	!j limites ZR/ZC
	  endif		!tot_conv	
	 endif		!ovsht
	
c_______________determination de la fonction de poids fac __________________
 
c	 interpolation aux nouveaux q de P, T, R, L, M pour le calcul de la
c	 fonction d'espacement aux new
c	 write(6,*)'ctep,ctem,lnp,lnt,r2,l23,m23,esp,new'
c	 write(6,2000)ctep,ctem
 
c	 print*
c	 print*,lim
c	 print*,(new(jlim(i)),i=1,lim)
c	 print*,(jlim(i),i=1,lim)	
 
	 do i=1,n
	  call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,new(i),lq,f,dfdq)
c	  esp(i)=ctep*f(1)+ctem*f(5)
	  esp(i)=ctep*f(1)+cter*f(3)+ctel*f(4)+ctem*f(5)
c	  write(6,2000)(f(j),j=1,5),esp(i),new(i)
	  if(i .gt. 1)then
	   do while(esp(i) .lt. esp(i-1))
	    esp(i)=esp(i)+0.01
	   enddo
	  endif
	 enddo
	 psi=bp(6)
 
c	 esp(i) est la valeur de la fonction d'espacement aux nouveaux q
c	 determinee avec les anciens facteurs fac(i)
c	 calcul des nouveaux fac(i) de facon a amener un noeud sur chaque
c	 nouveau q
 
	 do i=1,n-1
	  fac(i)=psi/(esp(i+1)-esp(i))
	 enddo
c	 write(6,2000)(fac(i),new(i),i=1,n)
c	 print*,(lim,jlim(i),i=1,lim)
c	 print*,psi,(fac(jlim(i)),new(jlim(i)),i=1,lim)
c	 print*,(esp(jlim(i)-1),esp(jlim(i)),esp(jlim(i)+1),i=1,lim)
c	 pause
 
	endif		!lim=0
	
c	write(6,2000)(fac(i),i=1,n)
c	pause'sortie lim_zc'
 
	write(6,*)' '
	write(6,*)'------------- La Convection (fin) -----------------'
c	write(6,*)' '
 
	return
 
	end
 
c***********************************************************************

	subroutine causerie1_limzc( r2, n, lim, r_zc, m_zc, r_ov, 
     +              blocage_core, r_ov_max1, melange, limmfix,
     +              mfix, rox_test, calib1, calib2 )

c Routine d'interactivite en debut (premier appel) de la routine
c de repartition des couches et de determination des Zones Convectives
c 'lim_zc_3'.

c n : nombre de couches du modele
c r2 : les valeurs du carre du rayon (exprime en rsol !)
c 'blocage_core' : .true. : la taille du coeur convectif ne peut depasser 
c                           la valeur 'r_ov_max1'


	implicit none
 
	include 'cesam_3.parametres'
	include 'modele_3.common'

	integer i, lim, n, ipara

	logical ovsht, blocage_core, melange, limmfix, calib1, calib2,
     +          rox_test

	real*8 r2(1), r_zc(1), r_ov(1), r_ov_max1, mfix,
     +        m_zc(1)

	character rep*1

c n : nombre de couches

         write(6,*)' ' 
	 write(6,*)'********************************************************* '
	 write(6,*)'*                                                       *' 
	 write(6,*)'*              REPARTITION et OVERSHOOTING              *'
	 write(6,*)'*                                                       *' 
	 print*,   '*               "lim_zc_3" version D.C. !               *'
	 write(6,*)'*                                                       *'
	 write(6,*)'********************************************************* '
	 write(6,*)' '
c         write(6,*)'                    AVERTISSEMENTS :'
c         write(6,*)' ' 
c	 print*, '        (1) on met du "Core OVERshoot" seulement si i=1 !'
c        write(6,*)' '
c       print*, '        (2) on a supprime les ecritures dans ''*.lis'''
c       write(6,*)' '
c       write(6,*)'       ****************************************** '
c       write(6,*)' '
c       write(6,*)'       >>  contenu du fichier de reprise :'
c       write(6,*)' '
c        print*, '             ovshts= ', ovshts
c        print*, '             ovshti= ', ovshti
c       write(6,*)' '
	 print*, ' '
	 print*, ' (1) Affiche-t-on les infos ? (o/n)'
	 read(5,*) rep
	 if ( rep .eq. 'o' ) then
	    do i= 1, lim
	       print*, '          r_zc(',i,')= ', r_zc(i)
	    end do
	    print*, ' '
	    do i= 1, lim
	       print*, '          r_ov(',i,')= ', r_ov(i)
	    end do
	    do i= 1, lim
	       print*, '          m_zc(',i,')= ', m_zc(i)
	    end do
	    print*, ' '
	    write(6,*)'       >>  rayon total :'
	    print*, '             Rstar= ', sqrt(r2(n))
	    write(6,*)' '
	 else
	    rep='n'
	 end if
c      write(6,*) '       >> Pour eviter des rayons d''overshoot deliran
c     +ts '
c      write(6,*)'          et de l''undershoot artificiel,on initialise'
c       write(6,*)  '          ''r_zc'' et ''r_ov'' a -100., ...'
c       write(6,*)  ' ' 
         do i=1,2*pnzc
            r_ov(i)=-100.d0
            r_zc(i)=-100.d0
         end do
c         print*, 'nouveau      r_zc= ', r_zc
c         print*, 'nouveau      r_ov= ', r_ov
         print*, ' '
	 print*, ' (2) Bloque-t-on la taille du coeur convectif ? (o/n)'
	 print*, '     (en fait on impose une borne sup. en rayon)'
	 read(5,*) rep
	 if ( rep .eq. 'o' ) then
	    blocage_core=.true.
	    print*, ' '
	    print*, 'Valeur maximale (en rsol) du coeur conv. :'
	    read(5,*) r_ov_max1
	 else
	    rep='n'
	    blocage_core=.false.
	    r_ov_max1=0.d0
	 end if
c        Option de Blocage en MASSE !
         print*, ' '
	 print*, ' (3) Fixe-t-on la taille du coeur convectif en MASSE ? (o/n)'
	 print*, '     (le coeur convectif a la meme taille en masse a tout instant)'
	 read(5,*) rep
	 if ( rep .eq. 'o' ) then
	    limmfix=.true.
	    print*, 'ATTENTION pour utiliser cette option il faut :'
	    print*, '    (i)  ne pas mettre d''overshooting.'
	    print*, '    (ii) ne pas utiliser Roxburgh.'
	    print*, 'Ok ?'
	    pause
	    print*, 'Valeur de Mzc a utiliser :'
	    read(5,*) mfix
	 else
	    rep='n'
	    limmfix=.false.
	    mfix=-100.
	 end if
c       Option de prise en compte du RAYON dans la distribution des points
         print*, ' '
	 print*, ' (4) Prend-t-on en compte le RAYON dans la disribution des points ? (o/n)'	 
	 read(5,*) rep
	 if ( rep .eq. 'o' ) then
	    print*, '  Si on desire tenir compte du rayon pour'
	    print*, '  la distribution des points, indiquer une'
	    print*, '  valeur non-nulle pour ''cter'' :'
	    print*, ' '
	    print*, '  Valeur ce cter ='
	    read(5,*) cter
	 else
	    cter=0.0d0
	    rep='n'
	 end if
c       Options concernant la determiantion de la taille du coeur par integrales
c       de Roxburgh
	if ( Rox ) then
	   print*, ' '
	   print*, ' (5) Dans la routine de Roxburgh melange-t-on ? (o/n)'	 
	   read(5,*) rep
	   if ( rep .eq. 'o' ) then
	      melange=.true.
	   else
	      melange=.false.
	      rep='n'
	   end if
	   if ( melange ) then
	      print*, 'ON MELANGE DANS ROXBURGH OK ?'
	      pause
	   else
	      print*, 'ON ne MELANGE pas DANS ROXBURGH OK ?'
	      pause
	   end if
	   print*, ' '
	   print*, ' (6) Dans Roxburgh : calibre-t-on ? (o/n)'	 
	   read(5,*) rep
c          ---------------------------------------
	   if ( rep .eq. 'o' ) then
 10	      print*, '  * le parametre 1 ? (1) '
	      print*, '  * le parametre 2 ? (2) '
	      read(5,*) ipara
	      if ( ipara .eq. 1 ) then
		 calib1=.true.
		 calib2=.false.
	      end if
	      if ( ipara .eq. 2 ) then
		 calib1=.false.
		 calib2=.true.
	      end if
	      if ( (ipara .ne. 1) .AND. (ipara .ne. 2) ) then
		 print*, 'Pb. reponses permises : 1 ou 2 !'
		 go to 10
	      end if
	   else
	      calib1=.false.
	      calib2=.false.
	      rep='n'
	   end if
c          ---------------------------------------
	   print*, ' '
	   print*, ' (7) Active-t-on l''option de "test" ? (o/n)'	 
	   read(5,*) rep
	   if ( rep .eq. 'o' ) then
	      rox_test=.true.
	   else
	      rox_test=.false.
	      rep='n'
	   end if
	end if

         print*, '       ****************************************** '
	 print*, ' '

c        Eventuellement une pause ...
c	 if ( rep_pause .ne. 'n' ) then
c	    pause
c	 end if

c------------------------------------------------------------------------

	 write(6,201)ctep,cter,ctel,ctem
201	 format(1x,/,1x,'Constantes de repartition utilisees',/,
     &	1x,' ctep=',1pd9.2,'cter=',1pd9.2,'ctel=',1pd9.2,' ctem=',
     &  1pd9.2,/)
	
	 write(6,*)' '
	 write(6,*)'limites ZR/ZC fixees sur un noeud'
	 write(6,*)' '

c-------------------------------------------------------------------------
	 ovsht=ovshti .gt. 0.d0 .or. ovshts .gt. 0.d0
	 if(ovsht)then
	  if(jpz)then
	   if(ovshti .gt. 0.)then
	    write(6,*)'overshooting inferieur de JPZ des ZC: ',
     1	int(ovshti*100),'% Hp/Kip'
	   endif
	   if(ovshts .gt. 0.d0)then
	    write(6,*)'overshooting superieur de JPZ des ZC: ',
     1	int(ovshts*100),'% R noyau'
	   endif
	  else
	   if(ovshti .gt. 0.)then
	    write(6,*)'overshooting inferieur des ZC: ',
     1	int(ovshti*100),'% Hp'
	   endif
	   if(ovshts .gt. 0.)then
	    write(6,*)'overshooting superieur des ZC: ',
     1	int(ovshts*100),'% Hp'
	   endif
	  endif
	 else
	  write(6,*)'modele sans overshooting'
	 endif
	 write(6,*)' '


	 if ( blocage_core ) then
	    print*, ' '
	    print*, '                    ATTENTION '
	    print*, 'La taille du coeur conv. va etre bornee a :'
	    print*, 'r_ov_max1= ', r_ov_max1, ' rsol'
	    print*, 'Cette option est a utiliser avec les plus'
	    print*, 'grandes precautions ! On peut avoir apparition'
	    print*, 'D''une convection centrale en couche !'
	    print*, 'De plus lorsque r_ov= -100. ...'
	    print*, ' '
	    print*, 'OK ?'
	    pause
	 end if

	 return

	 end

c********************************************************************

	subroutine bloc_cc( r_ov, r_ov_max1 )

c Blocage de la taille du rayon du 'Convective Core'.

c r_ov : rayons des limites des Zones Convectives (inclus l'overshoot), en unite solaire.
c r_ov_max1 : taille max. du Convective Core (en rsol) admise.

	implicit none

	include 'modele_3.common'
	include 'ctephy.common'

	logical init

	real*8 r_ov(1), r_ov_max1

	data init / .true. /

	print*, ' '
	print*, '****************************************************'
	print*, 'ATTENTION : blocage du Convective Core ---ACTIF--- !'
	print*, '****************************************************'
	print*, ' '

	if ( r_ov(1) .gt. r_ov_max1 ) then
	   print*, 'On bloque le Convective Core !'
	   print*, 'Ancienne valeur de r_ov(1)= ', r_ov(1), ' rsol'
	   print*, ' '
	   r_ov(1)=r_ov_max1
	   print*, 'Nouvelle valeur de r_ov(1)= ', r_ov(1), ' rsol'
	   if ( init .AND. ( rep_pause .ne. 'n' ) ) then
	      init=.false.
	      pause
	   end if
	end if

	return

	end 
