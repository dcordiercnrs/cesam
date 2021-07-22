c*******************************************************************
c Standard Main Program for Yellow CESAM
c
		
	character*70 sub_phys(2)
	common/subphys/sub_phys

	external etat_eff_dc,opal96_tops99_2000,conv_jmj_3,nrj_nuc,
     1	lim_teff_3,edding_3,diffus_3,perte_masse,ctes_ceph,dess_ol_dig01

	sub_phys(1)='etat_eff_dc,opal96_tops99_2000,conv_jmj_3,nrj_nuc'
	sub_phys(2)='lim_teff_3,edding_3,diffus_3,perte_masse,ctes_ceph,dess_ol_dig01'

	call ycesam(etat_eff_dc,opal96_tops99_2000,conv_jmj_3,nrj_nuc,
     1	lim_teff_3,edding_3,diffus_3,perte_masse,ctes_ceph,dess_ol_dig01)
	
	stop
	
	end


