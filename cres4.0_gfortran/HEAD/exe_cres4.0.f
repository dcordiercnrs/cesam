c--------------------------------------------------------------------------
c
c                   Programme principal de CRES 4.0
c
c--------------------------------------------------------------------------
c
		
	character*70 sub_phys(2)
	common/subphys/sub_phys

	external etat_eff_dc,opal96_tops99_2003,conv_jmj_3,nrj_nuc,
     1	lim_teff_3,edding_3,diffus_3,perte_masse,ctes_ceph,dess_ol_dig01

	sub_phys(1)='etat_eff_dc,opal96_tops99_2003,conv_jmj_3,nrj_nacre'
	sub_phys(2)='lim_teff_3,edding_3,diffus_3,perte_masse,ctes_ceph,dess_ol_dig01'

	call cres(etat_eff_dc,opal96_tops99_2003,conv_jmj_3,nrj_nacre,
     1	lim_teff_3,edding_3,diffus_3,perte_masse,ctes_ceph,dess_ol_dig01)
	
	stop
	
	end


