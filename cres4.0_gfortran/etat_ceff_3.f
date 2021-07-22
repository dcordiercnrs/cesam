 
c***************************************************************************
 
	subroutine etat_ceff_3(p,t,xchim,deriv,
     1       ro,drop,drot,drox,drott,drotp,drotx,
     2       u,dup,dut,dux,dutt,dutp,dutx,dh1,dhe1,dhe2,degene)
c
c       adaptation du logiciel ceff ecrit par jorgen christensen_dalsgaard
c       a partir de p et t, il calcule ro et u et toutes les derivees
c       jusqu'au deuxieme ordre
c       annie baglin le 1 09 1991
c
c	MODIF:
c	ajout d'un appel a EFF puis ETAT_GONG2 en cas d'echec de convergence
c	P. Morel octobre 91
c	appel systematique a EFF si T <3000 ou si X <.1
 
c	version 3
 
	implicit none
 
	integer nnmax,npar
	parameter(nnmax=1301, npar = 3)
 
	include 'cesam_3.parametres'
	include 'modele_3.common'
	include 'evol_chim_3.common'
 
	integer icmlit,nit
 
	real*8 p,t,xchim(1),ro,drop,drot,drox,drott,drotp,drotx,
     1       u,dup,dut,dux,dutt,dutp,dutx,dh1,dhe1,dhe2,degene,
     2	dropp,dupp,pl,tl,fl,f,drofn,dropl,dropn,drotn,droxn,
     3	dpffn,droffn,drottn,droxxn,droftn,drotxn,dptxn,drofxn,
     4	dpfxn,dpftn,dfpl,dftl,dfxl,dfpn,dpttn,dftn,dpxn,dfxn,
     5	dptn,droppn,drotpn,betbet,betgam,unsro,h,ro2,psro2,
     6	dhfn,dhpn,dupn,dhtn,dutn,dhxn,duxn,dhttn,dhtxn,dhffn,
     7	dhftn,dhfxn,dhtpn,dhppn,ro3,psro3,dutpn,duttn,dutxn,duppn
	
	logical nosd,notd,deriv,pass,init
c	data init/.true./,pass/.true./
 
	real*8 datmod(31),idatmd(3)
	equivalence(datmod(29),idatmd(1))
 
	real*8 amm,amm2,amm3
	common/ln10c/amm,amm2,amm3
 
	integer ihvz,iprrad,ihmin
	real*8 anh0,anhe0
	common/eqscntc/ anh0,anhe0,ihvz,iprrad,ihmin
 
	real*8  av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     1	ct,crho,cpe,che,ca03,caa,ckh,car,ergev,exhm
	common/constsc/av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     1	ct,crho,cpe,che,ca03,caa,ckh,car,ergev,exhm
 
	real*8 xii1(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     1	dlt(4),gm1,tprh,trhp,rhxp,gmm1(4)
	common/eqstdc/xii1,ane,rho,ht,pt,cp,dad,dlt,gm1,tprh,trhp,rhxp,gmm1
 
 
	real*8 east(30),xii(30),dne(10),dph(20),he(20),pe(10),
     1	hi(20),pi(10),hh(20),ph(10),hr(20),pr(10),pcoul(10),hcoul(10)
	common/eqsoutc/east,xii,dne,dph,he,pe,hi,pi,hh,ph,hr,pr,pcoul,hcoul
 
	integer iomfll
	common/hvomcl/iomfll
 
	integer iab
	real*8 ab(10)
	common/hvabndc/ab,iab
 
	integer icnthv,iwrthv
	real*8  dptst0,dptst1
	common /hvcntl/icnthv,iwrthv,dptst0,dptst1
 
	integer icoulm,iclmit,iclmsd,nitdmu
	real*8 epsdmu,epssdr
	common/ccoulm/epsdmu,icoulm,iclmit,iclmsd,nitdmu,epssdr
 
	integer idgbcs,idgrhs,idgeos,idgopc,idgeng
	common/diagns/idgbcs,idgrhs,idgeos,idgopc,idgeng
 
	logical ok
	common /marche/ok
 
c	common pour la tabulation uniquement
	common /dersec/dropp,dupp
 
	data init/.true./,pass/.true./
 
	save init,pass
 
	if(init)then
	 call ctes_3
	 call setcnsc		!initialisation des ctes
	 anh0=.5		!initialisation ?
	 init=.false.
	 write(2,1)
	 write(6,1)
1	 format(//,5x,'equation d''etat CEFF',//)
	endif
 
	if(t .lt. 3000. .or. xchim(1) .lt. .2 .or. p .lt. 1.d3)then
	 if(pass)write(6,20)t,xchim(1),p
	 pass=.false.
20	 format(1x,'divers appels a EFF car 3000 > t=',1pd10.3,
     1	' ou .2 > X(1)=',1pd10.3, ' ou 1.d3 > P=',1pd10.3)
	 call etat_eff_ps_3(p,t,xchim,deriv,
     1       ro,drop,drot,drox,drott,drotp,drotx,
     2       u,dup,dut,dux,dutt,dutp,dutx,dh1,dhe1,dhe2,degene)
 
	 if(.not. ok)then	!probleme : appel etat_gong2
	  call etat_gong2_3(p,t,xchim,deriv,
     1       ro,drop,drot,drox,drott,drotp,drotx,
     2       u,dup,dut,dux,dutt,dutp,dutx,dh1,dhe1,dhe2,degene)
	  write(6,200)p,t,xchim(1),z0,ro,u
	 endif
	 return
	endif
 
c	index divers
 
	idgeos=3		!lie a l'appel de setf4?
	icmlit=0		!trouve nulle part!
	ihmin=0			!associe a h-
	iomfll=1		! utilise dans le calcul de omegac
	amm2=amm*amm
c
	ok=.true.
        pl=log10(p)
        tl=log10(t)
	nosd=.false.
	notd=.true.
	if(ihe4 .gt. 1)then
	 if(iz .gt. 1)then		!ihe4-1: indice de He3
          call eqstpc(pl,tl,xchim(1),xchim(ihe4)+xchim(ihe4-1),
     1	xchim(iz),nosd,notd,fl,nit)	
	 elseif(z_cte)then
          call eqstpc(pl,tl,xchim(1),xchim(ihe4)+xchim(ihe4-1),
     1	z0,nosd,notd,fl,nit)
         else
          call eqstpc(pl,tl,xchim(1),xchim(ihe4)+xchim(ihe4-1),
     1	1.d0-xchim(1)-xchim(ihe4)-xchim(ihe4-1),nosd,notd,fl,nit)
c         print*,'on passe par la?'
         endif
        else		!X uniquement
         call eqstpc(pl,tl,xchim(1),1.d0-xchim(1)-z0,
     1	z0,nosd,notd,fl,nit)
        endif
  	if(ok)then
 
c        sorties
 
c        ro et ses derivees
 
         ro=rho(1)
	 f=10**fl
 
c        derivees premieres de rho
 
	 drofn=rho(2)*rho(1)/f
         dropl=rho(2)/pt(2)
         dropn=(rho(1)/pt(1))*dropl
         drotn=(rho(1)/t)*(rho(3)-dropl*pt(3))
         droxn=rho(1)*(rho(4)-dropl*pt(4))*amm
         drop=dropn
         drot=drotn
         drox=droxn
 
c        derivees secondes de rho et p
 
c        1. derivees secondes naturelles dans le systeme ftx
 
	 dpffn=(pt(1)/(f*f))*((pt(5)/amm)+pt(2)*(pt(2)-1.d0))
         droffn=(rho(1)/(f*f))*((rho(5)/amm)+rho(2)*(rho(2)-1.d0))
         drottn=(rho(1)/(t*t))*((rho(8)/amm)+rho(3)*(rho(3)-1.d0))
	 dpttn=(pt(1)/(t*t))*((pt(8)/amm+pt(3)*(pt(3)-1.d0)))
         droxxn=rho(1)*(rho(10)+rho(4)*rho(4)*amm)*amm
         droftn=(rho(1)/(f*t))*(rho(2)*rho(3)+rho(6)/amm)
         drotxn=(rho(1)/t)*(rho(9)+rho(4)*rho(3)*amm)
	 dptxn=(pt(1)/t)*(pt(9)+pt(4)*pt(3)*amm)
         drofxn=(rho(1)/f)*(rho(7)+rho(4)*rho(2)*amm)
	 dpfxn=(pt(1)/f)*(pt(7)+pt(4)*pt(2)*amm)
         dpftn=(pt(1)/(f*t))*(pt(2)*pt(3)+pt(6)/amm)
 
c        2. derivees premieres naturelles de f dans le systeme ptx
 
         dfpl=1.d0/pt(2)
         dftl=-pt(3)/pt(2)
         dfxl=-pt(4)/pt(2)
         dfpn=(f/pt(1))*dfpl
         dftn=(f/t)*dftl
	 dpxn=pt(4)*pt(1)*amm
	 dfxn=-dfpn*dpxn
         dfxn=f*dfxl
	 dptn=pt(3)*pt(1)/t
 
c        3. derivees secondes naturelles dans le systeme ptx
 
         droppn=dfpn*dfpn*(droffn-dpffn*drofn*dfpn)
	 dropp=droppn
         drotpn=dfpn*(droftn-dfpn*dptn*droffn+drofn*dfpn*
     1	(dpffn*dptn*dfpn-dpftn))
	 betbet=-(dpttn+2.d0*dftn*dpftn+dftn*dftn*dpffn)
	 betgam=-(dptxn+dfxn*dpftn+dftn*dpfxn+dftn*dfxn*dpffn)
	 drotxn=drotxn+dftn*drofxn+dfxn*droftn+dftn*dfxn*droffn+betgam*dropn
	 drottn=drottn+2.d0*dftn*droftn+dftn*dftn*droffn+betbet*dropn
         drott=drottn
         drotp=drotpn
         drotx=drotxn
 
c        energie interne  u   ! attention jorgen calcule l'enthalpie
 
         unsro=1.d0/ro
         h=ht(1)
	 u=ht(1)-pt(1)*unsro
 
c        derivees premieres de u
 
         ro2=ro*ro
         psro2=pt(1)/ro2
	 dhfn=ht(2)/(f*amm)
         dhpn=dhfn*dfpn
	 dupn=dhpn-unsro+psro2*dropn
	 dhtn=(ht(3)-ht(2)*pt(3)/pt(2))/(t*amm)
	 dutn=dhtn+psro2*drotn
	 dhxn=ht(4)-ht(2)*pt(4)/pt(2)
	 duxn=dhxn+psro2*droxn
	 dux=duxn
	 dut=dutn
	 dup=dupn	
 
c        derivees secondes de u
 
c	 1- derivees naturelles de h dans le systeme ftx
 
	 dhttn=((ht(8)/amm)-ht(3))/(t*t*amm)
	 dhtxn=ht(9)/(t*amm)
	 dhffn=((ht(5)/amm)-ht(2))/(f*f*amm)
	 dhftn=ht(6)/(f*t*amm2)
	 dhfxn=ht(7)/(f*amm)
 
c	 2- derivees naturelles de h dans le systeme ptx
 
	 dhttn=dhttn+2.d0*dftn*dhftn+dftn*dftn*dhffn+betbet*dhpn
	 dhtpn=dfpn*(dhftn-dfpn*dptn*dhffn+dhfn*dfpn*
     1	(dpffn*dptn*dfpn-dpftn))
	 dhtxn=dhtxn+dftn*dhfxn+dfxn*dhftn+dftn*dfxn*dhffn+betgam*dhpn
	 dhppn=dfpn*dfpn*(dhffn-dpffn*dhfn*dfpn)
 
c	 3- derivees de u dans le systeme ptx
 
         ro3=rho(1)*ro2
	 psro3=2.d0*psro2/rho(1)
	 dutpn=dhtpn+drotn/ro2-psro3*dropn*drotn+psro2*drotpn
	 duttn=dhttn-psro3*drotn*drotn+psro2*drottn
	 dutxn=dhtxn-psro3*drotn*droxn+psro2*drotxn+drotn*dptn/ro2
	 duppn=dhppn+(2.d0/ro2)*drop+psro2*dropp-psro3*drop*drop
	 dupp=duppn
	 dutt=duttn
	 dutp=dutpn
	 dutx=dutxn
 
c        degres d'ionisation
 
	 dh1=xii1(1)
	 dhe1=xii1(2)
	 dhe2=xii1(3)
	 degene=5.
 
	else	!appel a EFF car il y a pb avec CEFF pour les cas extremes
 
	 call etat_eff_ps_3(p,t,xchim,deriv,
     1       ro,drop,drot,drox,drott,drotp,drotx,
     2       u,dup,dut,dux,dutt,dutp,dutx,dh1,dhe1,dhe2,degene)
 
	 if(.not. ok)then	!probleme : appel etat_gong2
	  call etat_gong2_3(p,t,xchim,deriv,
     1       ro,drop,drot,drox,drott,drotp,drotx,
     2       u,dup,dut,dux,dutt,dutp,dutx,dh1,dhe1,dhe2,degene)
	 endif
 
	 write(6,200)p,t,xchim(1),z0,ro,u
200	format(1x,' p=',1pd10.3,' t=',1pd10.3,' x=',1pd10.3,
     1	' z0=',1pd10.3,' ro=',1pd10.3,' u=',1pd10.3)
 
	endif
 
        return
 
        end
