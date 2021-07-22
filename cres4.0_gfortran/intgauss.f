 
c*************************************************************************
 
	subroutine intgauss(a,b,x,w,n)
 
c	initialisation des poids et des abscisses pour l'integration de gauss
 
	implicit none
 
	integer n,i
 
	real*8 a,b,x(1),w(1)
 
	if(n .gt. 6)then
	 write(6,*)'pas assez de cas dans intgauss n=',n
	 stop
	endif
 
	goto (1,2,3,4,5,6),n
1	x(1)=0.
	w(1)=1.
	goto 10
2	x(1)= .57735 02691 89626
	x(2)=-.57735 02691 89626
	w(1)=1.
	w(2)=1.
	goto 10
3	x(1)=0.
	x(2)= .77459 66692 41483
	x(3)=-.77459 66692 41483
	w(1)= .88888 88888 88889
	w(2)= .55555 55555 55556
	w(3)= .55555 55555 55556
	goto10
4	x(1)= .33998 10435 84856
	x(2)= .86113 63115 94053
	x(3)=-.33998 10435 84856
	x(4)=-.86113 63115 94053
	w(1)= .65214 51548 62546
	w(2)= .34785 48451 37454
	w(3)= .65214 51548 62546
	w(4)= .34785 48451 37454
	goto10
5	x(1)=0.
	x(2)= .53846 93101 05683
	x(3)= .90617 98459 38664
	x(4)=-.53846 93101 05683
	x(5)=-.90617 98459 38664
	w(1)= .56888 88888 88889
	w(2)= .47862 86704 99366
	w(3)= .23692 68850 56189
	w(4)= .47862 86704 99366
	w(5)= .23692 68850 56189
	goto10
6	x(1)= .23861 91860 82197
	x(2)= .66120 93864 66265
	x(3)= .93246 95142 03152
	x(4)=-.23861 91860 82197
	x(5)=-.66120 93864 66265
	x(6)=-.93246 95142 03152
	w(1)= .46791 39345 72691
	w(2)= .36076 15730 48139
	w(3)= .17132 44923 79190
	w(4)= .46791 39345 72691
	w(5)= .36076 15730 48139
	w(6)= .17132 44923 79190
10	do i=1,n
	 x(i)=(b-a)/2.*x(i)+(b+a)/2.
	 w(i)=w(i)/2.*(b-a)
	enddo
 
	return
 
	end
