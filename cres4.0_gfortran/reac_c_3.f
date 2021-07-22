 
c******************************************************************
 
	subroutine reac_c_3(t,r,nreac)
 
c	cycles PP, CNO et 3 alpha
c	taux de reaction thermonucleaires
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version 3
 
c references :
 
c	Caughlan et Al. Atomic Data and Nuclear Data Tables, 40, 283, 1988
 
c entree :
c	t : temperature
c	nreac : nombre de reactions
 
c sortie :
c	r(i) : ln( taux de reaction)
 
	implicit none
 
	include 'cesam_3.parametres'
 
	integer nreac,i
 
	real*8 t,r(1),a(0:5),p(0:5),t913p
 
	real*8 t9,lnt9a,t923,lnt9,lnt912,lnt932,lnt913,t913,lnt923,
     1	v0,v1,v2,v3,v4,he4abe8,be8agc12,ln38,r16(pnreac)
 
	logical init
	data init/.false./
 
	if(.not. init)then
	 init=.true.
 
	 ln38=log(1.d-100)	!log(0)
 
	 if(nreac .ne. pnreac)then
	  write(6,*)'dans reac_c le nombre de reaction differe de nreac',nreac
	  stop
	 endif
	endif
 
c	les facteurs t9
 
	t9=t*1.d-9
	lnt9=log(t9)
	lnt912=lnt9/2.
	lnt932=lnt912*3.
	lnt913=lnt9/3.
	t913=exp(lnt913)
	t913p=t913
	lnt923=lnt913*2.
	t923=exp(lnt923)
 
c=======================================================================
 
c	cycle PP
 
c========================================================================
 
c	reaction 1 : H(H,e+ nu)H2	z0=1, z1=1
 
	if(t9 .gt. 3.)then
	 r16(1)=ln38
 
	else
	 a(0)=1.		!coefficients du polynome en t913
	 a(1)=.123
	 a(2)=1.09
	 a(3)=.938
	 call polyder(a,3,0,t913P,p)	!algorithme de Horner
 
	 r16(1)=log(p(0))+log(4.01d-15)-lnt923-3.38d0/t913
	endif
c	write(6,2000)r16(1)
 
c	reaction 2 : H2(H,g)He3	z0=1, z1=1
 
	a(0)=1.		!coefficients du polynome en t913
	a(1)=.112
	a(2)=3.38
	a(3)=2.65
	call polyder(a,3,0,t913P,p)	!algorithme de Horner
 
	r16(2)=log(p(0))+log(2.24d+03)-lnt923-3.720d0/t913
c	write(6,2000)r16(2)
 
c	reaction 3 : He3(He3,H H)He4	z0=2, z1=2
 
	a(0)=1.		!coefficients du polynome en t913
	a(1)=.034
	a(2)=-0.522
	a(3)=-0.124
	a(4)=.353
	a(5)=.213
	call polyder(a,5,0,t913P,p)	!algorithme de Horner
 
	r16(3)=log(p(0))+log(6.04d10)-lnt923-12.276d0/t913
c	write(6,2000)r16(3)
 
c	reaction 4 : He4(He3,g)Be7	z0=2, z1=2
 
	lnt9a=lnt9-log(1.+4.95d-02*t9)
	r16(4)=log(5.61d6)+5./6.*lnt9a-lnt932-12.826d0/exp(lnt9a/3.)
c	write(6,2000)r16(4)
 
c	reaction 5 : Li7(H,He4)He4	z0=3, z1=1
 
	lnt9a=lnt9-log(1.+.759*t9)
 
	v0=exp(log(1.096d9)-lnt923-8.472d0/t913)
	v1=log(4.830d8)+5./6.*lnt9a-lnt932-8.472/exp(lnt9a/3.)
	v1=-exp(v1)
	v2=exp(log(1.06d10)-lnt932-30.442/t9)
 
	r16(5)=v0+v1+v2
	if(r16(5) .gt. 0.)then
	 r16(5)=log(r16(5))
	else
	 r16(5)=ln38
	endif
c	write(6,2000)r16(5)
 
c	reaction 6 : Be7(e-,nu g)Li7     pas de saturation t9 .ge. .001
 
	a(0)=1.		!coefficients du polynome en t913
	a(1)=-0.537
	a(2)=3.86
	call polyder(a,2,0,t913P,p)	!algorithme de Horner
	p(0)=log(p(0)+.0027/t9*exp(2.515d-3/t9))
 
	r16(6)=p(0)+log(1.34d-10)-lnt912
	if(t9 .gt. 3.)r16(6)=ln38
c	write(6,2000)r16(6)
 
c	reaction 7 : Be7(H,g)B8(,e+ nu)Be8(,He4)He4	z0=4, z1=1
 
	v0=exp(log(3.11d5)-lnt923-10.262/t913)
	v1=exp(log(2.53d3)-lnt932-7.306/t9)
 
	r16(7)=v0+v1		!Be7(H,g)B8
	if(r16(7) .gt. 0.)then
	 r16(7)=log(r16(7))
	else
	 r16(7)=ln38
	endif
c	write(6,2000)r16(7)
 
c======================================================================
 
c	cycle CNO
 
c=========================================================================
 
c	reaction 8 : C12(H,g)N13(,e+ nu)C13	z0=6, z1=1
 
	a(0)=1.
	a(1)=.03
	a(2)=1.19
	a(3)=.254
	a(4)=2.06
	a(5)=1.12
	call polyder(a,5,0,t913P,p)	!algorithme de Horner
 
	v0=exp(log(p(0))+log(2.04d7)-lnt923-13.690/t913-(t9/1.5)**2)
	v1=exp(log(1.08d5)-lnt932-4.925/t9)
	v2=exp(log(2.15d5)-lnt932-18.179/t9)
 
	r16(8)=v0+v1+v2
	if(r16(8) .gt. 0.)then
	 r16(8)=log(r16(8))
	else
	 r16(8)=ln38
	endif
c	write(6,2000)r16(8)
 
c	reaction 9 : C13(H,g)N14	z0=6, z1=1
 
	a(0)=1.
	a(1)=.03
	a(2)=0.958
	a(3)=0.204
	a(4)=1.39
	a(5)=0.753
	call polyder(a,5,0,t913P,p)	!algorithme de Horner
 
	v0=exp(log(p(0))+log(8.01d7)-lnt923-13.717/t913-(t9/2.)**2)
	v1=exp(log(1.21d6)-6./5.*lnt9-5.701/t9)
 
	r16(9)=v0+v1
	if(r16(9) .gt. 0.)then
	 r16(9)=log(r16(9))
	else
	 r16(9)=ln38
	endif
c	write(6,2000)r16(9)
 
c	reaction 10 : N14(H,g)O15(e+,nu)N15	z0=7, z1=1
 
	a(0)=1.
	a(1)=0.027
	a(2)=-0.778
	a(3)=-.149
	a(4)=0.261
	a(5)=0.127
	call polyder(a,5,0,t913P,p)	!algorithme de Horner
 
	v0=exp(log(p(0))+log(4.9d7)-lnt923-15.228/t913-(t9/3.294)**2)
	v1=exp(log(2.37d3)-lnt932-3.011/t9)
	v2=exp(log(2.19d4)-12.53/t9)
 
	r16(10)=v0+v1+v2
	if(r16(10) .gt. 0.)then
	 r16(10)=log(r16(10))
	else
	 r16(10)=ln38
	endif
 
c	write(6,2000)r16(10)
 
c	reaction 11 : N15(H,g)O16	z0=7, z1=1
 
	a(0)=1.
	a(1)=.027
	a(2)=0.219
	a(3)=0.042
	a(4)=6.83
	a(5)=3.32
	call polyder(a,5,0,t913P,p)	!algorithme de Horner
 
	v0=exp(log(p(0))+log(9.78d8)-lnt923-15.251/t913-(t9/0.450)**2)
	v1=exp(log(1.11d4)-lnt932-3.328/t9)
	v2=exp(log(1.49d4)-lnt932-4.665/t9)
	v3=exp(log(3.80d6)-lnt932-11.048/t9)
 
	r16(11)=v0+v1+v2+v3
	if(r16(11) .gt. 0.)then
	 r16(11)=log(r16(11))
	else
	 r16(11)=ln38
	endif
c	write(6,2000)r16(11)
 
c	reaction 12 : N15(H,He4)C12	z0=7, z1=1
 
	a(0)=1.
	a(1)=.027
	a(2)=2.62
	a(3)=0.501
	a(4)=5.36
	a(5)=2.60
	call polyder(a,5,0,t913P,p)	!algorithme de Horner
 
	v0=exp(log(p(0))+log(1.08d12)-lnt923-15.251/t913-(t9/0.522)**2)
	v1=exp(log(1.19d8)-lnt932-3.676/t9)
	v2=exp(log(5.41d8)-lnt912-8.926/t9)
	v3=exp(log(4.72d8)-lnt932-7.721/t9)
	v4=exp(log(2.2d9)-lnt932-11.418/t9)
 
	r16(12)=v0+v1+v2+.1*(v3+v4)	!0 to 1 = .1
	if(r16(12) .gt. 0.)then
	 r16(12)=log(r16(12))
	else
	 r16(12)=ln38
	endif
 
c	write(6,2000)r16(12)
 
c	reaction 13 : O16(H,g)F17(,e+ nu)O17	z0=8, z1=1
 
	r16(13)=log(1.5d8)-log(t923*(1.+2.13*(1.-exp(-.728*t923))))-16.692/t913
c	write(6,2000)r16(13)
 
c	reaction 14 : O17(H,He4)N14	z0=8, z1=1
 
	a(0)=1.
	a(1)=.025
	a(2)=5.39
	a(3)=0.940
	a(4)=13.5
	a(5)=5.98
	call polyder(a,5,0,t913P,p)	!algorithme de Horner
 
	v0=exp(log(p(0))+log(1.53d7)-lnt923-16.712/t913-(t9/0.565)**2)
	v1=exp(log(2.92d6)+lnt9-4.247/t9)
	v2=exp(log(4.81d10)+lnt9-16.712/t913-(t9/0.04)**2)
	v3=exp(log(5.05d-5)-lnt932-.723/t9)
	v4=exp(log(1.31d1)-lnt932-1.961/t9)
 
	r16(14)=v0+v1+(v3+v2+v4)*.1	!0 to 1 = .1
	if(r16(14) .gt. 0.)then
	 r16(14)=log(r16(14))
	else
	 r16(14)=ln38
	endif
c	write(6,2000)r16(14)
 
c====================================================================
 
c	cycle 3 alpha
 
c=====================================================================
 
c	reaction 15 : He4(He4 He4,g)C12	z0=2, z1=2, z2=2
 
	if(t9 .le. .08)then
	 a(0)=1.
	 a(1)=.031
	 a(2)=8.009
	 a(3)=1.732
	 a(4)=49.883
	 a(5)=27.426
	 call polyder(a,5,0,t913P,p)	!algorithme de Horner
 
	 v0=exp(log(7.4d5)-lnt932-1.0663/t9)
	 v1=exp(log(p(0))+log(4.164d9)-lnt923-13.49/t913-(t9/.098)**2)
 
	 he4abe8=log(v0+v1)		!ecran?
c	 write(6,*)'he4abe8'
c	 write(6,2000)he4abe8
 
	 a(0)=1
	 a(1)=.018
	 a(2)=5.249
	 a(3)=.65
	 a(4)=19.176
	 a(5)=6.034
	 call polyder(a,5,0,t913P,p)	!algorithme de Horner
 
	 v0=exp(log(1.3d2)-lnt932-3.3364/t9)
	 v1=exp(log(p(0))+log(2.51d7)-lnt923-23.57/t913-(t9/.0235)**2)
 
	 be8agc12=log(v0+v1)
c	 write(6,*)'be8agc12'
c	 write(6,2000)be8agc12
 
	 v0=log(2.9d-16)+he4abe8+be8agc12+
     1	log(.01+.2*(1.+4.*exp(-(.025/t9)**3.263)))-
     2	log(1.+4.*exp(-(t9/0.025)**9.227))
	 v1=log(1.35d-7)-lnt932-24.811/t9
 
	 r16(15)=exp(v0)+.1*exp(v1)		!0 to 1 = .1
	 if(r16(15) .gt. 0.)then
	  r16(15)=log(r16(15))
	 else
	  r16(15)=ln38
	 endif
c	 write(6,2000)r16(15)
 
	else
	 v0=exp(log(2.79d-8)-lnt9*3.-4.4027/t9)
	 v1=exp(log(1.35d-7)-lnt932-24.811/t9)
 
	 r16(15)=v0+.1*v1			!0 to 1 = .1
	 if(r16(15) .gt. 0.)then
	  r16(15)=log(r16(15))
	 else
	  r16(15)=ln38
	 endif
c	 write(6,2000)r16(15)
	endif
 
c	reaction 16 : C12(He4,g)O16	z0=6, z1=2
 
	v0=exp(log(1.04d8)-2.*lnt9-2.*log(1.+0.0489*t923)
     1	-32.120/t913-(t9/3.496)**2)
	v1=exp(log(1.76d8)-2.*t9-2.*log(1.+.2654/t923)-32.120/t913)
	v2=exp(log(1.25d3)-lnt932-27.499/t9)
	v3=exp(log(1.43d-2)+5.*lnt9-15.541/t9)
 
	r16(16)=v0+v1+v2+v3
	if(r16(16) .gt. 0.)then
	 r16(16)=log(r16(16))
	else
	 r16(16)=ln38
	endif
c	write(6,2000)r16(16)
 
c	reaction 17 : O16(He4,g)Ne20	z0=8, z1=2
 
	v0=exp(log(9.37d9)-lnt923-39.757/t913-(t9/1.586)**2)
	v1=exp(log(6.21d1)-lnt932-10.297/t9)
	v2=exp(log(5.38d2)-lnt932-12.226/t9)
	v3=exp(log(1.3d1)+2.*t9-20.093/t9)
 
	r16(17)=v0+v1+v2+v3
	if(r16(17) .gt. 0.)then
	 r16(17)=log(r16(17))
	else
	 r16(17)=ln38
	endif
c	write(6,2000)r16(17)
 
2000	format(1x,1p8d10.3)
 
	do i=1,nreac
	 r(i)=r16(i)
	enddo
 
	return
 
	end
