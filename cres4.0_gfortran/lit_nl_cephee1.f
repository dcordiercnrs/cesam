c*********************************************************************
       
      subroutine lit_nl_cephee1
 
 
c	lecture des NAMELISTs
 
c	Auteurs : P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	     	D. Cordier, DASGAL, Observatoire de Paris - Meudon
c               Version du 2 Novembre 1998
c	CESAM, version 3


	implicit none
 
	include 'cesam_3.parametres'
	include 'modele_3.common'
	include 'atmosphere_3.common'
	include 'evol_chim_3.common'
	
	integer long

        external long

        logical estilla
 
	character*100 fich(8)
	common/fich_etat/fich

        logical enre_bin, lcroissant, tdstest
 
	real*8 log_teff, log_lslsol, y_stop, x_stop,t_stop !pour test d'arret d'evolution
 
        real*8 dtrep

	character arret*4, key_word*5
 
        common/enrep/dtrep,enre_bin,key_word
 
        common/tdstest/tdstest
 
	common/arret/log_teff,log_lslsol,y_stop,
     +               x_stop,t_stop,arret,lcroissant                                    	
	namelist/nl_espace/mtot,mdot,der_num,tau_max
        namelist/nl_num/precision,dtmax,dtmin,psi0_c,precix,precit,
     +                  d_grav,repri,tdstest,dty
 
	namelist/nl_temps/agemax,dtlist,x_stop,y_stop,log_teff,
     +                    lcroissant, log_lslsol,t_stop,arret
c       y_stop et loglslsol : nouveaux parametres. (DC 2/11/98)
 
        namelist/nl_enreg_mod/enre_bin,key_word,dtrep
 
	namelist/nl_chim/x0,y0,z_cte,diffusion
	namelist/nl_rot/w_rot,rot_solid	
	namelist/nl_conv/alpha,ovshts,ovshti,para_rox1,para_rox2,jpz,rox,ledoux
	namelist/nl_etat/fich
        namelist/nl_nuc/nuc1,nuc2,nuc3
	namelist/nl_opa/opa1,opa2,opa3,opa4,opa5,opa6,opa7,opa8
	
        inquire(file=nom_fich2(:long(nom_fich2))//'.don3',exist=estilla)
        if( .NOT. estilla ) then
           print *, ' '
           print *, '*************************************************'
           print *, ' '
           print *, 'pb. : le fichier : ',
     +               nom_fich2(:long(nom_fich2))//'.don3'
           print * ,'est absent !'
           print *, ' '
           print *, '**************************************************'
           print *, ' '
           stop
        end if
 
	open(unit=3,form='formatted',status='old',
     1	file=nom_fich2(:long(nom_fich2))//'.don3')
	
	write(6,*)'	donnees luees dans :'
	write(6,*)
	write(6,*)'     ',nom_fich2(:long(nom_fich2))//'.don3'
	write(6,*)	
	read(3,nl_espace)
	write(6,nl_espace)
        read(3,nl_num)
        write(6,nl_num)	
	read(3,nl_temps,err=100)
	write(6,nl_temps)
 
        read (3,nl_enreg_mod,err=200)
        if ( long(key_word) .gt. 5 ) then
           print *, ' '
           print *, 'Pb. <<key_word>> est trop long (sup. a 5 caract.)'
           print *, ' '
           stop
        end if
 
        write(6,nl_enreg_mod)
 
	read(3,nl_chim)
	write(6,nl_chim)	
	read(3,nl_rot)
	write(6,nl_rot)	
	read(3,nl_conv)

c------------------------------------------------------------------------
        if ( jpz .AND. rox ) then
           write(6,*) 'Pb. : on ne peut pas faire du Roxburg et du Zahn'
           write(6,*) 'en meme temps ! Corrigez *.don3'
           stop
        end if

        if ( (ovshts .le. 0.d0) .AND. Rox ) then
           print*, 'ATTENTION pour utiliser le critere de Roxburg'
           print*, 'On doit avoir ovshts > 0. !'
           stop
        end if
c-----------------------------------------------------------------------

	write(6,nl_conv)	
	read(3,nl_etat)
	write(6,nl_etat)
        read(3,nl_nuc)
        write(6,nl_nuc)	
	read(3,nl_opa)
	write(6,nl_opa)	
	close(unit=3)
	write(6,*)
	
	write(2,*)'	donnees luees dans :'
	write(2,*)
	write(2,*)'     ',nom_fich2(:long(nom_fich2))//'.don3'
	write(2,*)	
	
	write(2,nl_espace)
	write(2,nl_num)
	write(2,nl_temps)	
	write(2,nl_chim)	
	write(2,nl_rot)	
	write(2,nl_conv)	
	write(2,nl_etat)	
	write(2,nl_opa)	

c Resume

        write(6,*) ' '
        write(6,*) '*** Petit resume ***'
        write(6,*) ' '
        write(6,*) 'M= ', mtot
        write(6,*) 'X0=', x0, ' , Y0=', y0, ', Z0=', 1.-x0-y0
        write(6,*) 'dtmax= , ', dtmax
        write(6,*) 'd_grav= ', d_grav
        if ( (.NOT. rox) .AND. (.NOT. jpz) ) then
           write(6,*) 'Critere de Sch., OVERshoot= ', ovshts, ' Hp'
        end if
        if ( rox ) then
           write(6,*) 'Limite ZC/ZR (coeur) avec integrale de Roxburg'
           write(6,*) 'Valeurs des parametres utilises :'
           write(6,*) 'para_rox1= ', para_rox1
           write(6,*) 'para_rox2= ', para_rox2
        end if
        if ( jpz ) then
           write(6,*) 'Critere de J.P. Zahn'
           write(6,*) 'OVERshoot= ', ovshts, ' Hp'
        end if

        if ( rep_pause .ne. 'n' ) then
           pause
        end if

	return

 100    print *, ' '
        print *, 'Pb. : le fichier : ',
     +           nom_fich2(:long(nom_fich2))//'.don3'
        print *, 'n''est pas au bon format pour le namelist'
        print *, ' '
        print *, '         << nl_temps >>'
        print *, ' '
        stop
 
 200    print *, ' '
        print *, 'Pb. : le fichier : ',
     +           nom_fich2(:long(nom_fich2))//'.don3'
        print *, 'n''est pas au bon format pour le namelist'
        print *, ' '
        print *, '         << nl_enreg_mod >>'
        print *, ' '
        stop

	end
