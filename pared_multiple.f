c----   viene de pared_ramiro pero calcula para varios valores
c----   de tdp, con polvo (amin,amax,p,filesil) fijo

	program pared
	 
	integer ncdsx,i,j,kk,escribir,nwlirs,nwl,est00,jpaola
	parameter (ncdsx=100)
	parameter(ntdpt=15)

	character*100 wlfile,fileimagen,filespectro,filestar
	character*24 d_and_t

	real*4 step,nu,tauext,tau00,gg,mui,integra
	real*4 dummy(10)
	real*4 cs,Rd,mu,seni,sent,altura,delta,alt
	real*4 Rstar,Mstar,Tstar,mmu,xlbb,xnu
	real*4 X(0:ncdsx)
	real*4 Y(-ncdsx:ncdsx)
	real*4 Z
	real*4 iYup,sYup,sYdown
	real*4 Iterm(0:ncdsx,-ncdsx:ncdsx),Itotp(0:ncdsx,-ncdsx:ncdsx)
	real*4 Iscat(0:ncdsx,-ncdsx:ncdsx),Itot(0:ncdsx,-ncdsx:ncdsx)  
	real*4 W,Bdel,Bscat
	
 
	real*4 Msol,Rsol,G,au,mh
		parameter (Msol=1.989e33,Rsol=6.95508e10, 
     1             G=6.67e-8,au=1.495979e13,mh=1.67e-24)		
	real*4 h,k,c,Pi
	parameter (h=6.6260755e-27,  k=1.380658e-16, c=2.99792458e10,
     1         Pi=3.14159)      
	parameter(nwl0=4000,ntaumaxx=5000,nflll=nwl0,nmax=nwl0)

	character*120 filesil,filetroi,fileice1,fileice2
	common/dustfiles/filesil,filetroi,fileice1,fileice2

c       almacena los datos necesarios para calcular la I y el
c       flujo en el principal
        common/tdpvarios/tdpt(ntdpt),tkt(ntaumaxx,ntdpt),rdt(ntdpt),
     .  alturat(ntdpt),xkdt(ntdpt),ntdp
c       almacena las opacidades individuales para sumarlas despues
c       con las tsub de acuerdo con las tk
        common/opacidadesvarias/polsil(nmax),police(nmax),
     .  poltroi(nmax),polorg(nmax),albesil(nmax),albeice(nmax),
     .  albetroi(nmax),albeorg(nmax)

c.	sale de pared_arbitraria:
	dimension taud(ntaumaxx),tk(ntaumaxx)
	dimension wl(nwl0),alb(nwl0),xk(nwl0),xls(nwl0),dusto(nwl0)
	dimension xlbb(nwl0)
c	rd en R_*, altura desde el plano medio en R_*
	common/salida/rd,altura,xkd,taud,tk,xk,xls,ntaumax
c.	entra a pared_arbitraria:
c.	amin y amax en micras (tipico: 0.005 y 0.25), p=3.5
	common/sizes/amin,amax,p
c.	rstar y mstar en unidades de radio estelares en pared_arbitraria_iso
	common/estrella/Rstar,Mstar
	common/otros/xmu,alt,mmu,wlfile,filestar
	common/masotros/tdp,xmp,alfa
	common/alt1/fz0nada,scalenada,rdcmnada,deltar
	common/salva/dusto,alb
	common/ww/wl,nwl
	common/lambdas/nwlirs
	common/tstar1/tstar


	real*4 fluu(nflll,ntdpt)
	real*8 flujo(nflll),faprox(nflll),area,angsol

	external fscat, plank, integra
	data torg,ttroi,tice/425.,680,1./

	read(5,*)ntdp
	read(*,*)(tdpt(i),i=1,ntdp)
		write(6,*)(tdpt(i),i=1,ntdp)
	read(*,*)xmp,alfa,ntaumax
c       en unidades solares
	read(*,*)Mstar, Tstar, Rstar
	rstarori=rstar
	read(*,*)distancia
c	definimos cosas para la subtrutina pared arbitraria
	read(*,*)mui
	xmu=mui
	mmu=2.36
c.	en unidades de la escala de altura a teff. Poner = 1, el flujo escala
	read(*,*)alt
	alt0=alt
c	para que lambdas queremos
	read(*,'(a)')wlfile
c----   propiedades del polvo
        read(*,*)amin,amax,p
	write(17,*)amin,amax,p
	read(*,'(a)')filesil
	read(*,'(a)')filetroi
	read(*,'(a)')fileice1
	read(*,'(a)')fileice2
c       estrella
c	read(*,'(a)')filestar
c       Se escribe el archivo del mapa ? 1=si
c	'Estrella en la imagen (si=1, extingida=2)?'
	read(*,*)est00
	write(6,*)'est00',est00

	if(est00.eq.2) then
	   write(*,*)'tau?'
	   read(*,*)tauext
	endif

c call and print time to screen	
c	call fdate(d_and_t)
c	write(*,*)'TIME calling tyopas: ', d_and_t(12:19)
c	write(*,*)'TIME calling tyopas: '
	
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c 	write(*,*)'begin with tyopas'
	call tyopas
c	write(*,*)'done with tyopas'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c call and print time to screen	
c	call fdate(d_and_t)
c	write(*,*)'TIME done with tyopas:',d_and_t(12:19)
c    write(*,*)'TIME done with tyopas'
c	stop
	
	
c	pasamos a cm lo que sale de pared_arbitraria

c	for input to propio3
	write(20,*)(rdt(itdp)*rstar*rsol/au,itdp=1,ntdp)

	do itdp=1,ntdp
	rstar=rstarori
c-------------------------
c	llena los arreglos correspondientes
	tdp=tdpt(itdp)
	rd=rdt(itdp)
	altura=alturat(itdp)
	xkd=xkdt(itdp)
		do itau=1,ntaumax
		tk(itau)=tkt(itau,itdp)
		enddo
c-------opacidades, kappa, alb monocromaticos
	do il=1,nwlirs
		dustsil1=polsil(il)
        	xksil1=(1.-albesil(il))*dustsil1

                if (tdp.lt.torg) then
        	dustorg1=polorg(il)
        	xkorg1=(1.-albeorg(il))*dustorg1
                else
        	dustorg1=0.
         	xkorg1=0.
                endif

        	if (tdp.lt.tice) then
        	dustice1=police(il)
        	xkice1=(1.-albeice(il))*dustice1
        	else
        	dustice1=0.
        	xkice1=0.
        	endif

		if (tdp.lt.ttroi) then
        	dusttroi1=poltroi(il)
        	xktroi1=(1.-albetroi(il))*dusttroi1
        	else
        	dusttroi1=0.
        	xktroi1=0.
		endif
        dusto(il)=dustsil1+dustice1+dustorg1+dusttroi1
        xk(il)=xksil1+xkice1+xkorg1+xktroi1
	alb(il)=(dusto(il)-xk(il))/dusto(il)
c	de il
	enddo
c-------------------------

	rd0=rd
	write(*,*)'rd=',rd0
	rdau=rd0*rstar*rsol/au

	write(6,*)'rstar,rd',rstar,rd
	Rstar=rstar*rsol
	Rd=rd*Rstar
	write(6,*)'rd,rstar,altura',rd,rstar,altura
	altura=altura*Rstar
	write(6,*)'altura',altura
	
	seni=sqrt(1-mui**2)
	delta=(altura/Rd)*(seni/mui)
	step=Rd/float(ncdsx)
	area=step**2

	write(6,*)'seni,delta,step,area,altura'
	write(6,*)seni,delta,step,area,altura

c	calculamos el factor de dilucion de una vez para 'Rd' dado
	dummy(1)=Rstar/Rd
	W=0.5*(1-sqrt(1-(dummy(1))**2))

	angsol=area/(distancia*206265.*au)/(distancia*206265.*au)

	fac=(rd/(distancia*206265.*au))**2

	if(delta.gt.1.) then
	   angsolp=pi*fac*mui
	   else
	   angsolp=2.*fac*mui*
     * (delta*sqrt(1.-delta**2)+asin(delta))
	   endif

c	write(6,*)'w,angsol,fac,angsolp',w,angsol,fac,angsolp

c.	OJO: esta version no calcula el factor de asimetria consistentemente !!
	gg=0.24

c	tau del pixel central
	tau00=0.

cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c       Aqui empieza el bucle para las diferentes lambdas
c call and print time to screen	
c	call fdate(d_and_t)
c	write(*,*)'TIME begin pared do loop:',d_and_t(12:19)
c	write(*,*)'TIME begin pared loop:'
	
	do kk=1,nwlirs

c	lambda en micras
	nu=c/(wl(kk)*1e-4)

c	Bscat=plank(nu,Tstar)
c	if (kk.eq.1) then 
c	call fdate(d_and_t)
c	write(*,*)'TIME begin Bscat:',d_and_t(12:19)
c	end if
	Bscat=bnu(tstar,wl(kk))
c	if (kk.eq.1) then 
c	call fdate(d_and_t)
c	write(*,*)'TIME end Bscat:',d_and_t(12:19)
c	end if
	
c       ahora calculamos para todos los puntos y escribimos
	
	flujo(kk)=0.
	do i=0,ncdsx
	   X(i)=step*float(i)
	   do j=-ncdsx,ncdsx
c	if((kk/100)*100.eq.kk)write(6,*)itdp,kk,i,j
	      Y(j)=step*float(j)
              Z=Y(j)*mui/seni - sqrt(Rd**2-X(i)**2)/seni
c       el seno del angulo theta  (azimutal) en las coordenadas del
c       disco
	      sent=(Y(j)*mui-Z*seni)/Rd
	      mu=seni*sent

c	write(6,*)'sent',sent

	      dummy(2)=dusto(kk)/xkd/mu	 
c		  write(6,*)'dusto(kk),xkd,mu',dusto(kk),xkd,mu   
	  
c       primero por los casos posibles en delta
	      sYup= mui*sqrt(1-(X(i)/Rd)**2)*Rd+altura*seni
	      iYup=-mui*sqrt(1-(X(i)/Rd)**2)*Rd+altura*seni
	      sYdown= mui*sqrt(1-(X(i)/Rd)**2)*Rd-altura*seni
	      
c       en este bloque preguntamos si el punto esta dentro de la zona 
c       visible y evaluamos si no I=0

	      if (delta.gt.1.0) then 
c       toda la elipse de arriba
		 if ((Y(j).gt.iYup).and.(Y(j).lt.sYup)) then
c call and print time to screen	
c	call fdate(d_and_t)
c	write(*,*)'TIME begin integra:',d_and_t(12:19)
		    Iterm(i,j)=
     .			integra(dummy(2),nu,taud,tk,ntaumaxx,ntaumax)
c call and print time to screen	
c	call fdate(d_and_t)
c	write(*,*)'TIME end integra:',d_and_t(12:19)
		    Iscat(i,j)=W*Bscat*alb(kk)*fscat(gg,mu)/(1.+mu)
		    Itot(i,j)=Iterm(i,j)+Iscat(i,j)
		 else
		    Itot(i,j)=0.
		 end if
		 
	      else			
c       ya se ve el hueco y preguntamos
		 if(sent.lt.delta) then 
		    if ((Y(j).gt.iYup).and.(Y(j).lt.sYup)) then
c call and print time to screen	
c	call fdate(d_and_t)
c	write(*,*)'TIME begin integra:',d_and_t(12:19)
		       Iterm(i,j)=
     .			integra(dummy(2),nu,taud,tk,ntaumaxx,ntaumax)
c call and print time to screen	
c	call fdate(d_and_t)
c	write(*,*)'TIME end integra:',d_and_t(12:19)
		       Iscat(i,j)=W*Bscat*alb(kk)*fscat(gg,mu)/(1.+mu)
		       Itot(i,j)=Iterm(i,j)+Iscat(i,j)
c.	prueba
c	Itot(i,j)=1.
		    else
		       Itot(i,j)=0.
		    end if
		 else
		    if ((Y(j).gt.sYdown).and.(Y(j).lt.sYup)) then
c call and print time to screen	
c	call fdate(d_and_t)
c	write(*,*)'TIME begin integra:',d_and_t(12:19)
		       Iterm(i,j)=
     .			integra(dummy(2),nu,taud,tk,ntaumaxx,ntaumax)
c call and print time to screen	
c	call fdate(d_and_t)
c	write(*,*)'TIME end integra:',d_and_t(12:19)
		       Iscat(i,j)=W*Bscat*alb(kk)*fscat(gg,mu)/(1.+mu)
		       Itot(i,j)=Iterm(i,j)+Iscat(i,j)
c.	prueba
c	Itot(i,j)=1.
		    else
		       Itot(i,j)=0.
		    end if
c	if(i.lt.65)then
c	write(6,*)'i,j,kk,wl(kk),Bscat,Iterm(i,j),Iscat(i,j),Itot(i,j)'
c	write(6,*)i,j,kk,wl(kk),Bscat,Iterm(i,j),Iscat(i,j),Itot(i,j)
c	endif
c	if(i.lt.65)then
c	write(6,*)'i,Iterm(i,j)'
c	write(6,*)i,Iterm(i,j)
c	endif

		 end if
		 
	      end if

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    aqui extinguimos la intesidad de la PARED y la ESTRELLA
c    si es que se pide.   1.085.. es  2.5log(e)

	      if(est00.eq.2)  then
		 Itot(i,j)=Itot(i,j)*exp(-1.08573*tauext)
		 Bscat=Bscat*exp(-1.08573*tauext)
	      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


	      if((i.eq.0).and.(j.eq.0)) then
		 xfactor=Pi*rstar**2/(step**2)
		 Itotp(i,j)=Itot(i,j)+Bscat*xfactor
	      else
		 Itotp(i,j)=Itot(i,j)
	      end if
	      flujo(kk)=flujo(kk)+2*Itot(i,j)*angsol
	   end do

	end do
	
	fluu(kk,itdp)=flujo(kk)
c	if(itdp.eq.2)write(6,*)wl(kk),flujo(kk),fluu(kk,itdp)
	end do
c 	termino el ciclo de las lambdas 
c	call fdate(d_and_t)
c	write(*,*)'TIME end pared do loop:',d_and_t(12:19)
c    write(*,*)'TIME finish pared loop:'
	write(6,*)'itdp',itdp

c	de itdp
	enddo

	
c	write(*,*)'writing output'
	write(17,*)ntdp
	write(17,*)'Tdp'
	write(17,200)(tdpt(itdp),itdp=1,ntdp)
200	format(10(1pe10.2))
c	write(17,201)(rdt(itdp)*rstar*rsol/au,itdp=1,ntdp)
	write(17,*)'Rd(AU)'
	write(17,200)(rdt(itdp)*rstarori*rsol/au,itdp=1,ntdp)
	write(17,*)'zwall(AU)(1 H)'
	write(17,200)(alturat(itdp)*rstarori*rsol/au,itdp=1,ntdp)
	write(17,*)'wl(mic),log(nuFnu)'

	do i=1,nwlirs
	 xnu=c/(wl(i)*1e-4);
	 write(17,202)wl(i),(alog10(xnu*fluu(i,itdp)),itdp=1,ntdp)
202	 format(15(0pf8.3))
	end do
	close(17)

	
	end program pared
	
