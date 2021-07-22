c******************************************************************
	
c	etat_eff : equation d'atat
c	opa_int_zsx : opacites
c	conv_jmj : convection
c	ppcno9 : reac. thermo. PP+CNO
c	lim_atm : reconstitution d'une atmosphere
c	hopf : loi T(tau)
c	diffm_mp : diffusion microscopique
c	difft_cte0 : diffusion turbulente
c	perte_ext : perte de masse
c	ctes : constantes de physique
c	des_r : dessin en cours d'execution	
	
	implicit none
	
	include 'modele.common'

	external etat_eff_dc,opa_daniel,conv_jmj,ppcno3a_dc,
	1	lim_atm,edding,diffm_mp,difft_cte,perte_masse,ctes,
	2	dess_ol_dig01
		
	sub_phys(1)='etat_eff_dc,opa_daniel,conv_jmj,ppcno3a_dc'
	sub_phys(2)='lim_atm,edding,diffm_mp,difft_cte,perte_masse,ctes,
     1 dess_ol_dig01'	 

	call cesam(etat_eff_dc,opa_daniel,conv_jmj,ppcno3a_dc,
	1	lim_atm,edding,diffm_mp,difft_cte,perte_masse,ctes,
	2	dess_ol_dig01)

	stop
	
	end
