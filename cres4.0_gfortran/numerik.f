c*********************************************************************
 
c  Block contenant les routines numeriques:
 
c  interpol_lin    : interpolation lineaire
c  interpol_Tcheby : interpolation par polynomes de Tchebychev
c  resol_syst      : resolution d un systeme lineaire par Gauss
c  int_ext_parab   : interpolation et extrapolation parabolique
c  locate          : localisation d'une valeur dans un tableau
 
c**********************************************************************
 
      subroutine interpol_lin(f, x, nx, x0, fx0, dfsdx)
 
      implicit none
 
      integer i, nx
 
      real*8 f, x, x0, fx0, dfsdx,
     +       a, b, xi, xi1, yi, yi1
 
      dimension f(nx), x(nx)
 
      do 100 i=2, nx
 
         if (x0.le.x(i)) then
         xi=x(i-1)
         xi1=x(i)
         yi=f(i-1)
         yi1=f(i)
         goto 200
         endif
 
100   continue
 
200   if(xi.eq.xi1) then
        print*,'Division par 0 dans interpol_lin !'
        stop
      endif
 
      a=(yi-yi1)/(xi-xi1)
      b=yi-a*xi
 
      fx0=a*x0+b
 
      dfsdx=a
 
      return
 
      end
 
c***********************************************************************
 
      subroutine  interpol_Tcheby ( f, x, nx, x0, fx0, dfsdx )
 
c-----------------------------------------------------------------------
 
c Programme d interpolation par polynomes de Tchebychev
 
c Ce programme fournit f(x0) et df/dx](x0)
c f(x0) est obtenue par la methode de Tchebychev
c df/dx](x0) obtenue par un taux de variation
 
c Domaine d utilisation: eviter les fonctions presentant une asymptote
c verticale ( meme en travaillant sur les log comme ici, avec par ex.
c x/(x^2-1)^(.5) avec x proche de 1 on a des problemes.
 
c subroutine: resol_syst
 
c Reference : Methodes de calcul numerique
c             ----------------------------
 
c             J. P. Nougier 3ieme edition Masson 1991
c                     pages 84 - 85
 
c Remarque : _ les x(j) doivent tous etre differents !
c            _ le nombre max de points d interpolation est fixe par nmax
 
 
      implicit none
 
      integer nx, i1, i2, nmax
 
      logical zerr
 
      real*8 f, x, x0, fx0, dfsdx, x_bis, x0_bis, x0_g, x0_d, x0_g_bis,
     +       x0_d_bis, Tl, Tl_yg, Tl_y, Tl_yd, sum, a_l, fx0_g, fx0_d,
     +       Tl_t, decal, lnf
 
      parameter ( nmax=50 )
 
      dimension  f(nx), x(nx), x_bis(nmax), a_l(nmax), 
     +           Tl(nmax,nmax), Tl_yg(nmax), Tl_y(nmax), Tl_yd(nmax),
     +           Tl_t(nmax,nmax), lnf(nmax)
 
c Test nx < nmax
 
      if ( nx .gt. nmax ) then
         print*,'Dans interpol_Tcheby nx est superieur a nmax'
         stop
      end if
 
c On verifie que: x(1) .lt. x0  et x0.lt.x(nx)
 
      if ( x0.le.x(1) .or. x0.ge.x(nx) ) then
         print*,'Dans interpol_Tcheby x0 sort du domaine des x(j)'
         stop
      end if
 
c Determination de decal
 
      decal=0.D+00
 
      do i1= 1, nx, 1
         if ( f(i1).le.decal ) then
            decal=-f(i1)
         end if
      end do
 
c Decalage et passage aux log
 
      do i1= 1, nx
         lnf(i1)=log(1.D+00+decal+f(i1))
      end do
 
 
c Passage des x(j) aux x_bis(j) ( cf. formule 3.57 p. 83 )
 
      x0_bis=(2.D+00*x0-x(nx)-x(1))/(x(nx)-x(1))
 
      do i1= 1, nx, 1
         x_bis(i1)=(2.D+00*x(i1)-x(nx)-x(1))/(x(nx)-x(1))
      end do
 
c Recherche des x(j) tels que x(j) < x0 < x(j+1)
 
      do i1= 1, nx-1, 1
 
         if ( x(i1) .le. x0  .and.  x0 .le. x(i1+1) ) then
 
            x0_g=x0-(x0-x(i1))/1000.D+00
            x0_d=x0+(x(i1+1)-x0)/1000.D+00
 
            goto 100
 
         end if
 
      end do
 
c Conversion de x0_g et x0_d
 
100   x0_g_bis=(2.D+00*x0_g-x(nx)-x(1))/(x(nx)-x(1))
      x0_d_bis=(2.D+00*x0_d-x(nx)-x(1))/(x(nx)-x(1))
 
c calcul des Tl(xj_bis)
 
      do i1= 1, nx, 1
 
         Tl(i1,1)=1.
         Tl(i1,2)=x_bis(i1)
 
         do i2= 3, nx, 1
 
            Tl(i1,i2)=2.D+00*x_bis(i1)*Tl(i1,i2-1)-Tl(i1,i2-2)
 
         end do
 
      end do
 
c calcul des Tl(x0_g_bis)
 
      Tl_yg(1)=1.
      Tl_yg(2)=x0_g_bis
 
      do i1= 3, nx, 1
 
         Tl_yg(i1)=2.D+00*x0_g_bis*Tl_yg(i1-1)-Tl_yg(i1-2)
 
      end do
 
c calcul des Tl(x0_bis)
 
      Tl_y(1)=1.
      Tl_y(2)=x0_bis
 
      do i1= 3, nx, 1
 
         Tl_y(i1)=2.D+00*x0_bis*Tl_y(i1-1)-Tl_y(i1-2)
 
      end do
 
c calcul des Tl(x0_d_bis)
 
      Tl_yd(1)=1.
      Tl_yd(2)=x0_d_bis
 
      do i1= 3, nx, 1
 
         Tl_yd(i1)=2.D+00*x0_d_bis*Tl_yd(i1-1)-Tl_yd(i1-2)
 
      end do
 
c Calcul des a_l ( cf. systeme 3.60 p. 84 )
c     on transpose Tl en Tl_t car resol_syst modifie la matrice d entree
c     de plus dans resol_syst les indices de ligne et de colonne sont
c     permutes
 
      do i1= 1, nx
         do i2= 1, nx
            Tl_t(i1,i2)=Tl(i2,i1)
         end do
      end do
 
      call resol_syst ( Tl_t, nx, lnf, nx, zerr )
 
c attention f est modifie
 
      if ( zerr ) then
         print*,'le systeme lineaire a resoudre dans interpol_Tcheby est
     + singulier'
      stop
      end if
 
      do i1= 1, nx, 1
         a_l(i1)=lnf(i1)
      end do
 
c Calcul de fx0 ( formule 3.65 p. 85 )
 
      sum=0.D+00
 
      do i1= 1, nx, 1
         sum=sum+a_l(i1)*Tl_y(i1)
      end do
 
      fx0=exp(sum)-1.D+00-decal
 
c Calcul de df/dx ( elucubration personnelle )
 
      sum=0.D+00
 
      do i1= 1, nx, 1
         sum=sum+a_l(i1)*Tl_yg(i1)
      end do
 
      fx0_g=exp(sum)-1.D+00-decal
 
      sum=0.D+00
 
      do i1= 1, nx, 1
         sum=sum+a_l(i1)*Tl_yd(i1)
      end do
 
      fx0_d=exp(sum)-1.D+00-decal
      
      dfsdx=(fx0_d-fx0_g)/(x0_d-x0_g)
 
      return
 
      end
 
c***********************************************************************
 
      subroutine resol_syst(a, maxcol, bx, n, zerr)
 
c-----------------------------------------------------------------------
 
c Programme de resolution de systeme lineaire A.X = B
 
c entrees / sorties : 
c           _ a : matrice input ( modifiee au cours du calcul )
c                                 ----------------------------
c           Remarque: a( colonne , ligne )
c
c           _ maxcol : nbre de colonnes par ligne
c           _ bx : matrice colonne ( modifiee pendant le calcul )
c           _ n : ordre du systeme
c           _ zerr : variable logique zerr=.true. si le syst. est singulier
c                    zerr=.false. s il est soluble.
 
c Reference: 'La pratique du FORTRAN 77' P. Lignelet 2ieme edition Masson 1991
c            ( pages 126 - 127 )
 
      implicit none
 
      integer maxcol, n, k, l, lmax, kol, k1
 
      real*8 a, bx, epsi, biga, aux, bk, xl
 
      logical zerr
 
      dimension a(maxcol,*), bx(*)
 
      parameter ( epsi= 1.00D-10 )
 
 
      if ( n.lt.1  .or.  n.gt.maxcol ) goto 99
 
      do 10 k= 1, n
 
c Recherche dans la colonne k de la ligne lmax du + grand element
 
      biga=0.
 
      do 2 l= k, n
 
         aux=abs( a(k,l) )
         if ( aux.gt.biga ) then
            biga=aux
            lmax=l
         end if
 
2     continue
 
      if ( biga.lt.epsi ) goto 99
 
c Echange des lignes k et lmax, puis normalisation de la ligne k
 
      biga=1.D+00/a(k,lmax)
 
      do 3 kol= k, n
           aux=a(kol,lmax)
           a(kol,lmax)=a(kol,k)
           a(kol,k)=aux*biga
3     continue
 
      aux=bx(lmax)
      bx(lmax)=bx(k)
      bk=aux*biga
      bx(k)=bk
 
c Mise a zero de la zone (k+1:n) de la colonne k
 
      k1=k+1
 
      do 10 l= k1, n
         aux=a(k,l)
         do 5 kol= k1, n
5           a(kol,l)=a(kol,l)-aux*a(kol,k)
         bx(l)=bx(l)-aux*bk
10    continue
 
c Final: calcul 'en remontant' du vecteur des solutions
 
      do 20 l= n-1, 1, -1
         xl=bx(l)
         do 25 k= l+1, n
25          xl=xl-a(k,l)*bx(k)
20    bx(l)=xl
 
      zerr=.false.
      return
 
c lorsque le systeme est singulier zerr= .true.
 
99    zerr=.true.
 
      end
 
c***********************************************************************
 
      subroutine int_ext_parab( f, x, nx, x0, fx0, dfsdx )
 
c-----------------------------------------------------------------------
 
c Ce programme fournit f(x0) et df/dx](x0), f(x0) est obtenu par 
c interpolation parabolique.
 
c Rq.: dans le choix des 3 plus proches voisins on n envisage pas le cas
c      ou x0 n est pas dans l intervalle couvert par les 3 plus proches
c      voisins ( et c est peut etre mieux comme cela!)
 
c ATTENTION: les x(j) doivent etre tous differents avec des valeurs
c            croissantes ou decroissantes.
 
c subroutines appelees: * resol_syst33
c                       * locate
 
c entrees: * f : tableau des valeurs de la fonction a interpolee
c          * x : tableau des nx valeurs de l abscisse
c          * x0: points ou on interpole ou extrapole
 
c sortie: * fx0: valeur approchee de la fonction en x0
c         * dfsdx: valeur approchee de la derivee en x0
 
 
      implicit none
 
      integer nx, j, i1
 
      real*8 f, x, x0, fx0, dfsdx, u, v, tab_u, a, dist_1, dist_2,
     +       alpha, beta, x1_d, xn_g, f1_d, fn_g
 
      logical ext_inf, ext_sup
 
      dimension f(nx), x(nx), u(3), v(3), tab_u(3,3), a(3)
 
         ext_inf=.false.
         ext_sup=.false.
 
c Cas ou x0 est inferieur ou egal a x(1)
 
      if ( x0 .le. x(1) ) then
         u(1)=x(1)
         u(2)=x(2)
         u(3)=x(3)
 
         v(1)=f(1)
         v(2)=f(2)
         v(3)=f(3)
 
         ext_inf=.true.
 
         goto 100
      end if
 
c Cas ou x0 est superieur ou egal a x(nx)
 
      if ( x0 .ge. x(nx) ) then
         u(1)=x(nx-2)
         u(2)=x(nx-1)
         u(3)=x(nx)
 
         v(1)=f(nx-2)
         v(2)=f(nx-1)
         v(3)=f(nx)
 
         ext_sup=.true.
 
         goto 100
      end if
 
c Cas ou x(1) < x0 < x(nx)
 
      call locate(x, nx, x0, j)
      
      if ( j .eq. 1 ) then
         u(1)=x(1)
         u(2)=x(2)
         u(3)=x(3)
 
         v(1)=f(1)
         v(2)=f(2)
         v(3)=f(3)
 
         goto 100
      end if
 
      if ( j .eq. nx-1 ) then
         u(1)=x(nx-2)
         u(2)=x(nx-1)
         u(3)=x(nx)
 
         v(1)=f(nx-2)
         v(2)=f(nx-1)
         v(3)=f(nx)
 
         goto 100
      end if
 
      dist_1=abs(x(j-1)-x(j))
      dist_2=abs(x(j)-x(j+2))
 
      if ( dist_1 .le. dist_2 ) then
         u(1)=x(j-1)
         u(2)=x(j)
         u(3)=x(j+1)
 
         v(1)=f(j-1)
         v(2)=f(j)
         v(3)=f(j+1)
 
         goto 100
      else
         u(1)=x(j)
         u(2)=x(j+1)
         u(3)=x(j+2)
 
         v(1)=f(j)
         v(2)=f(j+1)
         v(3)=f(j+2)
 
         goto 100
      end if
 
c Construction du tableau tab_u
 
100   do i1= 1, 3, 1
         tab_u(i1,1)=(u(i1))**2
         tab_u(i1,2)=u(i1)
         tab_u(i1,3)=1.
      end do
       
c       print*,'u= ', u
c       print*,'tab_u= ', tab_u
c       print*,'v= ', v
 
c Resolution du systeme
 
      call resol_syst33( tab_u, a, v)
 
c traitement de l extrapolation inferieure
 
c     extrapolation par une droite d equation y=alpha*x+beta
c     avec alpha=dfsdx(x(1)) et beta= y1-alpha*x1
 
      if ( ext_inf ) then
         x1_d=x(1)+(x(2)-x(1))/1000.D+00
         f1_d=a(1)*x1_d**2+a(2)*x1_d+a(3)
         alpha=(f1_d-f(1))/(x1_d-x(1))
         beta=f(1)-alpha*x(1)
 
         fx0=alpha*x0+beta
         dfsdx=alpha
 
         return
      end if
 
c Traitement de l extrapolation superieure
 
      if ( ext_sup ) then
         xn_g=x(nx)-(x(nx)-x(nx-1))/1000.D+00
         fn_g=a(1)*xn_g**2+a(2)*xn_g+a(3)
         alpha=(fn_g-f(1))/(xn_g-x(1))
         beta=f(nx)-alpha*x(nx)
 
         fx0=alpha*x0+beta
         dfsdx=alpha
 
         return
      end if
 
c Cas de l extrapolation
 
      fx0=a(1)*x0**2+a(2)*x0+a(3)
      dfsdx=(2.D+00)*a(1)*x0+a(2)
 
      return
 
      end
 
 
c***********************************************************************
 
      subroutine resol_syst33( a, x, b )
 
c-----------------------------------------------------------------------
 
c sous-programme de resolution de systeme 3*3
 
c ATTENTION: le premier indice est l indice de COLONNE et le  deuxieme
c            l indice de LIGNE ! (ainsi les donnees sont rangees en
c            memoire ligne par ligne.
 
c a: matrice 3*3
c x: matrice colonne (3,1) solution du systeme
c b: matrice colonne (3,1)
 
 
      implicit none
 
      real*8 a, x, b, Da, D1, D2, D3
 
      dimension a(3,3), x(3,1), b(3,1)
 
c Calcul du determinant de A:
 
      Da=a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(2,1)*a(1,2)*a(3,3)+
     +   a(2,1)*a(1,3)*a(3,2)+a(3,1)*a(1,2)*a(2,3)-a(3,1)*a(1,3)*a(2,2)
 
c      print*,'Da= ', Da
 
      if ( Da .eq. 0. ) then
         print*,'Dans resol_syst33 le systeme a resoudre est singulier'
         stop
      end if
 
c Calcul des 3 determinants:
 
      D1=b(1,1)*a(2,2)*a(3,3)-b(1,1)*a(2,3)*a(3,2)-
     +   b(2,1)*a(1,2)*a(3,3)+b(2,1)*a(1,3)*a(3,2)+b(3,1)   
     +   *a(1,2)*a(2,3)-b(3,1)*a(1,3)*a(2,2)
 
      D2=a(1,1)*b(2,1)*a(3,3)-a(1,1)*a(2,3)*b(3,1)-a(2,1)*b(1,1)*a(3,3)+
     +   a(2,1)*a(1,3)*b(3,1)+a(3,1)
     +   *b(1,1)*a(2,3)-a(3,1)*a(1,3)*b(2,1)
 
      D3=a(1,1)*a(2,2)*b(3,1)-a(1,1)*b(2,1)*a(3,2)-a(2,1)*a(1,2)*b(3,1)+
     +   a(2,1)*b(1,1)*a(3,2)+a(3,1)
     +   *a(1,2)*b(2,1)-a(3,1)*b(1,1)*a(2,2)
 
c Ecriture de la solution:
 
      x(1,1)=D1/Da
      x(2,1)=D2/Da
      x(3,1)=D3/Da
 
      return
 
      end
 
c***********************************************************************
 
      subroutine locate(xx,n,x,j)
 
 
      integer j, n, j1, jm, ju
 
      real*8 x, xx(n)
 
c subroutine prise dans page 111 de 'Numerical Recipes in FORTRAN'
 
c Given an array xx(1:n), and given a value x, returns a value l such that
c x is between xx(j) and xx(j+1). xx(1:n) must be monotonic, either increasing
c or decreasing. j=0 or j=n is returned to indicate that x is out of range.
 
      j1=0
      ju=n+1
 
10    if ( ju-j1 .gt. 1 ) then
         jm=(ju+j1)/2
         if ( (xx(n).gt.xx(1)).eqv.(x.gt.xx(jm))) then
            j1=jm
         else
            ju=jm
         end if
      goto 10
 
      end if
 
      j=j1
 
      return
 
      end
         
 
