		
	character*70 sub_phys(2)
	common/subphys/sub_phys

	external etat_opal_dc,opal96_tops99_2001,conv_jmj_3,nrj_nacre,
     1	lim_teff_3,edding_3,diffus_3,gain_masse,ctes_ceph,no_des_3

	sub_phys(1)='etat_opal_dc,opal96_tops99_2001,conv_jmj_3,nrj_nacre'
	sub_phys(2)='lim_teff_3,edding_3,diffus_3,gain_masse,ctes_ceph,no_des_3'

	call cres(etat_opal_dc,opal96_tops99_2001,conv_jmj_3,nrj_nacre,
     1	lim_teff_3,edding_3,diffus_3,gain_masse,ctes_ceph,no_des_3)
	
	stop
	
	end


