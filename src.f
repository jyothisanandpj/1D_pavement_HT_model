program pavement_model

	integer, parameter :: dp = kind(1.d0), grid = 32
	integer :: i, dt = 1, period = 60*60*24*5, n = 1, m = 0, length, inter2, inter3, por, j
	real(kind=dp), dimension(4) :: row, cp, k, l, dx, alpha
	real(kind=dp), dimension(144) :: Tdb, RH, Spd, Sol, Tdp, TS
	real(kind=dp), dimension(432000) :: TAtm, WSpd, QSol, TDew, TSky, QTherm, HCoeff, QConv, QFlux
	real(kind=dp), dimension(126) :: a, b, c, d
	real(kind=dp), dimension(128,432001) :: T = 0
	real(kind=dp) :: albedo, emissivity, TGround = 38.5, psurf, prop, frac_step
	real(kind=dp), dimension(48) :: Temp_hourly, Flux_hourly
	
	open(unit=1, file="mat_prop.csv")
	read(unit=1, fmt=*) !avoid headers
	do i = 1,4
		read(unit=1, fmt=*) row(i),cp(i),k(i),l(i)
	enddo
	close(unit=1)
	
	open(unit=2, file="rad_prop.csv")
	read(unit=2, fmt=*) !avoid headers
	read(unit=2, fmt=*) albedo, emissivity
	close(unit=2)
	
	open(unit=3, file="weatherfile.data")
	read(unit=3, fmt=*) !avoid headers
	do i = 1, 144
		read(unit=3, fmt=*) Tdb(i), RH(i), Spd(i), Sol(i), Tdp(i) 
	enddo
	close(unit=3)
		
	TS = ((Tdb+273.15)*(0.004*Tdp + 0.8)**0.25) - 273.15

	dx = l/grid
	alpha = k/(row*cp)
	
	
	T(:,1) = 40 ! initial conditions
	inter2 = grid*2
	inter3 = grid*3
	length = grid*4
	
	print *, NEW_LINE('a'), NEW_LINE('a'), "The program is executing, please wait..",  NEW_LINE('a')
	
	!computation starts from here
	timestep: do i = 1, period
		frac_step = real(m,dp)/3600
		TAtm(i) = Tdb(n+1)*frac_step + Tdb(n)*(1-frac_step)
		QSol(i) = Sol(n+1)*frac_step + Sol(n)*(1-frac_step)
		WSpd(i) = Spd(n+1)*frac_step + Spd(n)*(1-frac_step)
		TSky(i) = TS(n+1)*frac_step  + TS(n)*(1-frac_step)
		m = m+1
		if (mod(i+1,3600) .EQ. 0.0) then
			m = 0
			n = n+1
		endif
		
		QTherm(i) = LW_rad(T(1,i), TSky(i), emissivity)
		HCoeff(i) = HCombined(TAtm(i), T(1,i), WSpd(i))
		QConv(i)  = HCoeff(i)*(TAtm(i) - T(1,i))
		QFlux(i) = QConv(i) + (1-albedo)*QSol(i) + QTherm(i)
		
		
	
	!CN method starts here for each timestep
		por = 1
		psurf = alpha(1)*dt/(2*dx(1)**2)
		
		gridspacing: do j = 1, 126
			
			if (j.EQ.grid .or. j.EQ.inter2 .or. j.EQ.inter3) then
				por = por + 1
			endif
			prop = alpha(por)*dt/(2*dx(por)**2)
			if (j.EQ.1) then
				a(j)=0
				b(j)=1+prop
				c(j)=-prop
			elseif (j.EQ.126) then
				a(j)=-prop
				b(j)=1+3*prop
				c(j)=0
			else
				a(j)=-prop
				b(j)=1+2*prop
				c(j)=-prop
			endif
			
			d(j) = prop*T(j,i) + (1-2*prop)*T(j+1,i) + prop*T(j+2,i)
	
		enddo gridspacing
		
		d(1) = d(1) + psurf*QFlux(i)*(dx(1)/k(1))
		d(126) = d(125) + 2*prop*TGround
		
		call tridiagonal(a,b,c,d)
	!t = 128 j values
		T(2:127,i+1) = d 
		
		!BC
		T(1,i+1) = T(2,i+1) + QFlux(i)*(dx(1)/k(1))
		T(128,i+1) = 2*TGround - T(127,i+1)
		
		!Flux balance for each layer speration
		T(grid+1,i+1) = T(grid+2,i+1) + (dx(2)/k(2))*(k(1)/dx(1))*(T(grid-1,i+1)-T(grid,i+1))
		T(inter2+1, i+1) = T(inter2+2,i+1) + (dx(3)/k(3))*(k(2)/dx(2))*(T(inter2-1,i+1)-T(inter2,i+1))
		T(inter3+1, i+1) = T(inter3+2,i+1) + (dx(4)/k(4))*(k(3)/dx(3))*(T(inter3-1,i+1)-T(inter3,i+1))
		
		!print *, T(1,i), QFlux(i), QTherm(i), QConv(i)
		
	enddo timestep
	
	Temp_hourly = Average(T(1,432001-172800:432001))
	Flux_hourly = Average(QConv(432000-172800:432000))	
	
	call write_to_csv(Temp_hourly,Flux_hourly)
	
	print *, "Done!!"
	print *, "Please close this window and look for OUTPUT.CSV file in the current folder", NEW_LINE('a'), NEW_LINE('a')
	print *, "To modify the pavement material properties (conductivity, density and heat capacity), please edit MAT_PROP.CSV file."
	print *, "To modify radiative properties (albedo and emissivity), please edit RAD_PROP.CSV file."
	print *, NEW_LINE('a'), NEW_LINE('a')
	print *, "                                ==============================================================="
	print *, "                                    Thanks for using UCRC 1-D pavement heat transfer model     "
	print *, "                                ==============================================================="
	
	read (*,*)

	contains
	
	subroutine write_to_csv(var1,var2)
		
		implicit none
		integer :: i, n , hour = 0
		real(kind=dp), dimension(:), intent(in) :: var1, var2
		character(56) :: header
		character(9) :: dat 
		n = size(var1)
		dat = "6/25/2004"
		open(unit=1, file="output.csv")
		15 format(A9,",",i2,",",F7.2,",",F7.2)
		5 format(A56)
		header = "Date,Time(h),Surface Temperature(C),Surface flux(W/sq.m)"
		write(unit=1, fmt=5) header
		
		do i = 1, n 
			write(unit=1, fmt=15) dat, hour, var1(i), var2(i)
			hour = hour + 1
			if (hour .GT. 23) then
				hour = 0
				dat = "6/26/2004"
			endif
		enddo 
		close(unit=1)
		
	end subroutine write_to_csv
	
	subroutine tridiagonal (a,b,c,d)
	
		implicit none
		integer :: i, p = 126 
		real(kind=dp), intent(in), dimension(126) :: a, c
		real(kind=dp), intent(inout), dimension(126) :: b, d 
		
		forward_loop: do i = 2, p
			b(i) = b(i) - c(i-1)*a(i)/b(i-1)
			d(i) = d(i) - d(i-1)*a(i)/b(i-1)
		enddo forward_loop
		
		d(p) = d(p)/b(p)
		
		backward_loop: do i = p-1, 1, -1
			d(i) = (d(i) - c(i)*d(i+1))/b(i)
		enddo backward_loop
		
		
	end subroutine tridiagonal
	
	function Average(Propsurf) result(AvgProp)
		
		implicit none
		integer :: i, n, m 
		real(kind=dp), intent(in), dimension(172800) :: Propsurf
		real(kind=dp), dimension(48) :: AvgProp
		
		do i = 1, 48
			n = ((i-1)*60*60)+1
			m = i*60*60
			AvgProp(i) = real(sum(Propsurf(n:m)),dp)/3600
		enddo
		
		
	end function Average
	
	real(kind=dp) function HCombined(T1, T2, Spd)
	
		implicit none
		real(kind=dp), intent(in) :: T1, T2, Spd
		real, parameter :: a = 3.26, b = 0.89, beta = 0.0, rf =1.52
		real(kind=dp) :: natural, h_cg
		
		if (T1 .ge. T2) then
		natural = 1.810*(abs(T1-T2)**(0.3333333333334))/(1.382 + abs(cos(beta)))
		else
		natural = 9.482*(abs(T1-T2)**(0.3333333333334))/(7.283 - abs(cos(beta)))
		endif
		
		h_cg = sqrt(natural**2 + (a*(Spd**b))**2)
		HCombined = natural + rf*(h_cg-natural)
	
	
	end function HCombined
	
	
	real(kind=dp) function LW_rad(TSurf, TSky, E)
	
		implicit none
		real(kind=dp), intent(in) :: TSurf, TSky, E
		real(kind=dp), parameter :: sigma = 5.67/10**(8)
		
		LW_rad = E * sigma * ((TSky+273.15)**(4) - (TSurf+273.15)**(4))
		
	end function LW_rad


end program pavement_model