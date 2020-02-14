
	function integra(alpha,nu,taud,tk,nn,n)
	dimension taud(nn),tk(nn)
	real*4 nu,alpha,dummy,delta,integra
	character*24 d_and_t
	
c	call fdate(d_and_t)
c	write(*,*)'TIME begin integra:',d_and_t(12:19)
	xsuma=0.
	do i=1,n-1
	   delta=taud(i+1)-taud(i)
	   
	   dummy=plank(nu,tk(i))*exp(-alpha*taud(i))*delta
c	   write(*,*)'plank(nu,tk(i))',plank(nu,tk(i))
c	   write(*,*)'exp(-alpha*taud(i))',exp(-alpha*taud(i))
c	   write(*,*)'plank(nu,tk(i))*exp(-alpha*taud(i))',plank(nu,tk(i))*exp(-alpha*taud(i))
	   xsuma=xsuma+dummy
	end do

	integra=alpha*xsuma
c	call fdate(d_and_t)
c	write(*,*)'TIME end integra:',d_and_t(12:19)
c	stop
c	write(*,*)'taud(i),tk(i)'
c	write(*,*)taud(i),tk(i)
c	write(*,*)'integra'
c	write(*,*)integra
	
	return
	end 

