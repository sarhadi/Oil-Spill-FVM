	
	PROGRAM Hydrodynamic

!	Model.f90 
!
!*****************************************************************************************************
!
!	PROGRAM: Depth Average Finite Volume Model with Alternating Direction Implicit Method
!
!	PURPOSE: 2D Modeling of Flow and Water Quality Processes due to Transport and Fate of Oil Spills
!
!*****************************************************************************************************

!	***   January 2009
!
!	***	Ehsan Sarhadi zadeh
!
!	***   M.Sc. Project, KNTU University, Tehran, Iran

!*********************************
!	Definition of Variables
!*********************************

!	|Z0			=Initial Water Surface Elevation
!	|Cd			=Air/Fluid Drag Coefficient
!	|VolOIL		=Total Volume of the Oil Slick
!	|HWave		=Offshore Deep Water Wave Height
!	|Period		=Wave Period
!	|Slope		=Slope of the Plane Beach
!	|Rip		=Rip Current Spacing
!	|Dx			=Grid Spacing in X-Direction
!	|Dy			=Grid Spacing in Y-Direction
!	|Dt			=Time Step
!	|U			=X-Velocity Component for each Velocity Control Volume in X-Direction
!	|V			=Y-Velocity Component for each Velocity Control Volume in Y-Direction
!	|FluxU		=Amount of U Crossing the Unit of Surface per Unit of Time in X-Direction
!	|FluxV		=Amount of V Crossing the Unit of Surface per Unit of Time in Y-Direction
!	|FUStar		=Temporary Values of FluxU for Correction Algorithm
!	|FVStar		=Temporary Values of FluxV for Correction Algorithm
!	|Time		=Simulation Time
!	|Grav		=Earth's Gravity Acceleration
!	|Cori		=Coriolis Parameter due to the Earth's Rotation
!	|AngLat		=Geographical Angel of Latitude
!	|Omega		=Angular Rotation Speed of the Earth
!	|Pw			=Density of Fluid
!	|Pa			=Density of Air
!	|Po			=Density of OIL
!	|Pi			=Pi Number 3.1415926535897932
!	|Thick		=OIL Slick Thickness
!	|Cfilm		=OIL Film/Water Surface Friction Coefficient
!	|MINoil		=Minimum OIL Slick Thickness Distribution
!	|CENTERx	=OIL Slick Coordinate in X-Direction
!	|CENTERy	=OIL Slick Coordinate in Y-Direction
!	|Diameter	=OIL Slick Diameter
!	|Radius		=OIL Slick Radius
!	|SIGMA		=OIL Slick Distribution Parameter
!	|OILSource	=OIL Spill Source (OIL Injection)
!	|OILSink	=OIL Spill Sink (Evaporation, Emulsification, ...)
!	|Beta		=Momentum Correction Factor for a non-Uniform Vertical Velocity Profile
!	|Karman		=von-Karman Coefficient
!	|Ke			=Roughness Lenght
!	|Ce			=Eddy Viscosity Coefficient
!	|A1,...,H1	=Temporary Array for Matrix Coefficients
!	|A0,...,H0	=Temporary Array for Matrix Coefficients for First Cell
!	|AX,BX,CX,ZX=Matrix Coefficients
!	|CX1		=Continuity Components for Matrix Coefficients
!	|MX0		=Momentum Components for First Cell
!	|MX1		=Momentum Components in X-Direction for Matrix Coefficients
!	|MX2		=Momentum Components in Y-Direction for Matrix Coefficients
!	|Vw			=Wind Velocity at 10 m above the water surface
!	|Twx		=Shear Stresses due to Wind in X-Direction
!	|Twy		=Shear Stresses due to Wind in Y-Direction
!	|WINDteta	=Direction to which the Wind blows measured Counter Clockwise from East
!	|LANDA		=Slick Spreading/Diffusion Function
!	|EVAP		=Evaporation due to Wind 
!	|Z			=Water Surface Elevation
!	|H			=Total Water Depth
!	|D			=Water Depth Below Datum
!	|WINDx		=Temporary Array for Wind Terms in X-Direction
!	|WINDy		=Temporary Array for Wind Terms in Y-Direction
!	|Friction	=Temporary Array for Friction Terms
!	|Vegetation	=Temporary Array for Vegetation Terms
!	|Cvg		=Vegetation Drag Coefficient
!	|Mvg		=Vegetation Density
!	|Dvg		=Vegetation Diameter
!	|ADVv		=Temporary Array for Advection Terms 
!	|DIFx		=Temporary Array for Diffusion Terms in X-Direction
!	|DIFy		=Temporary Array for Diffusion Terms in Y-Direction
!	|FUp		=Temporary Array for Averaged Flux Velocity in X-Direction
!	|FVp		=Temporary Array for Averaged Flux Velocity in Y-Direction
!	|Coriolis	=Temporary Array for Coriolis Terms
!	|AAA,...,DDD=Temporary Array for Matrix Coefficients
!	|C,CXX		=Temporary Array for Matrix Coefficients
!	|UU,VV,ZZ	=Temporary Values for Variables at the end of Time Step
!	|HH,HHX		=Temporary Values for Variables at the end of Time Step
!	|S			=OIL Slick Thickness in Grid Cells
!	|SX			=Temporary Values for OIL Slick Thickness at the end of Time Step
!	|WET		=Wet Cells
!	|Ah,Ahh		=Depth Averaged Turbulent Eddy Viscosity in Grid Cells & Its Temporary Values
!	|Cchezy		=Chezy Roughness Coefficient
!	|X			=Cell Number in X-Direction (Domain Length)
!	|Y			=Cell Number in Y-Direction (Domain Width)
!	|M			=Numbered Length
!	|N			=Numbered Width
!	|Bath		=Bathymetry Data
!	|HPl		=Half Periodic Sinusoidal Bed Length
!	|HPc		=Half Periodic Sinusoidal Lower Beach Coefficients
!	|CL			=Coastline Distance from Datum
!	|UB			=Upper Beach Elevation
!	|FBd		=Flat Bed Depth
!	|Count		=Counter Values
!	|Direction	=Direction of Velocity or Flux for Solving Shape Function
!	|Correct	=Counter Values for Correction Algorithm
!	|CorrectTest=Minimum Values for Correction Algorithm Test
!	|TECz,TECs	=Temporary Values for TECPLOT

!*********************************
!	Declarations
!*********************************

	IMPLICIT NONE 

	REAL  Z0,Cd,VolOIL,HWave,Period,Slope,Rip
	REAL  Dx,Dy,Dt,Time,Grav,Cori
   	REAL  AngLat,Omega,Pw,Pa,Po,Pi
	REAL  Thick,Cfilm,MINoil,CENTERx,CENTERy,Diameter,Radius,SIGMA
	REAL  OILSource,OILSink,Karman,Ke,Ce,ZX,HPl,HPc,CL,UB,FBd
	REAL,ALLOCATABLE:: A1(:,:),B1(:,:),C1(:,:),D1(:,:),E1(:,:),F1(:,:),G1(:,:),H1(:,:)
	REAL,ALLOCATABLE:: F0(:),G0(:),H0(:),XX(:),BX(:),AX(:,:),CX1(:,:),MX1(:,:),MX2(:,:),MX0(:)
	REAL,ALLOCATABLE:: Vw(:,:),Twx(:,:),Twy(:,:),WINDteta(:,:),LANDA(:,:),EVAP(:,:)
	REAL,ALLOCATABLE:: Z(:,:),H(:,:),D(:,:),U(:,:),V(:,:),FluxU(:,:),FluxV(:,:),WINDx(:,:),WINDy(:,:)
	REAL,ALLOCATABLE:: Friction(:,:),ADVv(:,:),DIFx(:,:),DIFy(:,:),FVp(:,:),FUp(:,:),Coriolis(:,:)
	REAL,ALLOCATABLE:: AAA(:),BBB(:),CCC(:),DDD(:),CXX(:),Bath(:,:),Beta(:,:),TECz(:,:),TECs(:,:)
	REAL,ALLOCATABLE:: Vegetation(:,:),Cvg(:,:),Mvg(:,:),Dvg(:,:),FUStar(:,:),FVStar(:,:),CorrectTest(:,:)
	REAL,ALLOCATABLE:: UU(:,:),VV(:,:),ZZ(:,:),Ahh(:,:),S(:,:),SX(:,:),WET(:,:),Ah(:,:),Cchezy(:,:)
	INTEGER,ALLOCATABLE:: UDirection(:,:),VDirection(:,:)
	INTEGER X,Y,Q,M,N,I,J,P,R,Count,Correct

!*********************************
!	Input/Output Communications
!*********************************

	OPEN (10, FILE="INPUT.TXT")
	OPEN (20, FILE="TOTAL.DAT")
	OPEN (30, FILE="Boundary.TXT",status="unknown")	
	OPEN (40, FILE="WET.TXT",status="unknown")
	OPEN (50, FILE="Island.TXT",status="unknown")
	OPEN (60, FILE="DRY.TXT",status="unknown")
	OPEN (70, FILE="Bathymetry.DAT",status="unknown")
	OPEN (80, FILE="Bathymetry.TXT",status="unknown")

	READ (10,*) AngLat,Ke,Ce,Grav,Pi,Time,X,Y,Dx,Dy,Dt,Z0,Pa,Pw,Po,Thick,VolOIL,Cfilm,Karman
	READ (10,*) MINoil,CENTERx,CENTERy,Diameter,OILSource,OILSink,HWave,Period,Slope,Rip,HPl,HPc,CL,UB,FBd

	M=2*X+2
	N=2*Y+2

    ALLOCATE (XX(MAX(M-1,N-1)),BX(MAX(M-1,N-1)),AX(MAX(M-1,N-1),MAX(M-1,N-1)),F0(MAX(M,N)),G0(MAX(M,N)),H0(MAX(M,N)))
	ALLOCATE (CX1(M,N),MX1(M,N),MX2(M,N),MX0(MAX(M,N)),ADVv(M,N),DIFx(M,N),DIFy(M,N),CXX(MAX(M,N)))
	ALLOCATE (Vw(M,N),Twx(M,N),Twy(M,N),WINDteta(M,N),LANDA(M,N),EVAP(M,N),Coriolis(M,N))
	ALLOCATE (FVp(M,N),FUp(M,N),WINDx(M,N),WINDy(M,N),Friction(M,N),S(M,N),SX(M,N))
	ALLOCATE (Vegetation(M,N),Cvg(M,N),Mvg(M,N),Dvg(M,N),FUStar(M,N),FVStar(M,N),CorrectTest(M,N),TECz(M,N),TECs(M,N))
	ALLOCATE (Z(M,N),H(M,N),U(M,N),V(M,N),D(M,N),ZZ(M,N),Ahh(M,N),UU(M,N),VV(M,N)) 
	ALLOCATE (AAA(MAX(M,N)),BBB(MAX(M,N)),CCC(MAX(M,N)),DDD(MAX(M,N)))
	ALLOCATE (A1(M,N),B1(M,N),C1(M,N),D1(M,N),E1(M,N),F1(M,N),G1(M,N),H1(M,N))
	ALLOCATE (WET(M,N),Ah(M,N),Cchezy(M,N),Bath(M,N),FluxU(M,N),FluxV(M,N),Beta(M,N),UDirection(M,N),VDirection(M,N))

!*********************************
!	TecPlot Setup
!*********************************

	WRITE(20,*), 'Title="Hydrodynamic"'
	WRITE(20,*), 'variables="i","j","z","s","h","d","u","v","w"'
	
	WRITE(70,*), 'Title="Hydrodynamic"'
	WRITE(70,*), 'variables="i","j","b"'
    WRITE(70,*), 'zone t="zone number-' ,1,'",i=',M,',j=',N,''
    WRITE(70,*), 'f=point'

!*********************************
!	Wetting & Drying Setup
!*********************************



!*********************************
!	Bathymetry Setup
!*********************************

	DO J=1,N
	DO I=1,M

	IF (I<=(2*INT(CL/Dx)+1)) THEN 
	Bath(I,J)=((I-(2*INT(CL/Dx)+1))*UB/(2*INT(CL/Dx)+1))
    ELSE IF (I>=(2*(HPl+CL)/Dx)) THEN
	Bath(I,J)=FBd 	  
	ELSE 
	Bath(I,J)=Slope*((I-(2*INT(CL/Dx)))*(Dx/2)-HPc*SIN(Pi*(I-(2*INT(CL/Dx)))*(Dx/2)/(HPl))*SIN(2*Pi*(J)*(Dy/2)/(Rip)))
	END IF

	WRITE(70,*),I,J,Bath(I,J)

    END DO
	END DO

!*********************************
!	Initial & Boundary Condition
!*********************************

	D=Bath
	P=0
	Count=0
	U=0.0
	V=0.0
	FluxU=U*H
	UU=0.0
	VV=0.0
	Z=Z0
	H=D+Z
	S=0 
	SX=0
	WET=1
	EVAP=0
	Radius=INT(Diameter/2)
    SIGMA=SQRT(-(Radius*Dx)**2/(2*LOG(MINoil)))
	Omega=2*Pi/(24.*3600.) 
	Cori=2*Omega*SIN(Pi*AngLat)
		
	DO I=1,2*X+1
	DO J=1,2*Y+1
	    
	IF (H(I,J)<0.00) H(I,J)=0

	END DO 
	END DO

	DO I=1,2*X+1
	DO J=1,2*Y+1

	Cd=(0.63+0.066*Vw(I,J))*1.E-3
	Twx(I,J)=Cd*Pa*ABS(Vw(I,J))*Vw(I,J)*COS(Pi*WINDteta(I,J))
	Twy(I,J)=Cd*Pa*ABS(Vw(I,J))*Vw(I,J)*SIN(Pi*WINDteta(I,J))

	END DO 
	END DO

!*********************************
!	Initial OIL Condition
!*********************************

	DO I=62,64,2
	DO J=8,10,2

	S(I,J)=0.01
	  
	END DO
	END DO

!*********************************
!	Starting the Main Loop (Time/dT)+1
!*********************************

	DO P=1,(Time/dT)+1

!*********************************
!	Updating TecPlot
!*********************************

	Q=MOD(P,1)
	IF (Q.eq.0.or.P.eq.1) 	THEN
	
	Count=Count+1

	WRITE(20,*), 'zone t="zone number-',Count,'",i=',2*X,',j=',2*Y,''
	WRITE(20,*), 'f=point'

	TECz=Z
	TECs=S

	DO J=1,2*Y
	DO I=1,2*X

	IF (I==1 .or. J==1) THEN

	TECz(I,J)=(Z(I+1,J+1))
	IF (MOD(J,2)==0)	TECz(I,J)=(Z(I+1,J))
	IF (MOD(I,2)==0)	TECz(I,J)=(Z(I,J+1)+Z(I+2,J+1))/2

	ELSE

	IF (MOD(I,2)/=0 .and. MOD(J,2)/=0)	TECz(I,J)=(Z(I+1,J+1)+Z(I-1,J+1)+Z(I+1,J-1)+Z(I-1,J-1))/4
	IF (MOD(I,2)/=0 .and. MOD(J,2)==0)	TECz(I,J)=(Z(I+1,J)+Z(I-1,J))/2
	IF (MOD(I,2)==0 .and. MOD(J,2)/=0)	TECz(I,J)=(Z(I,J+1)+Z(I,J-1))/2

	END IF

	IF (I==1 .or. J==1) THEN

	TECs(I,J)=(S(I+1,J+1))
	IF (MOD(J,2)==0)	TECs(I,J)=(S(I+1,J))
	IF (MOD(I,2)==0)	TECs(I,J)=(S(I,J+1))

	ELSE

	IF (MOD(I,2)/=0 .and. MOD(J,2)/=0)	TECs(I,J)=(S(I+1,J+1)+S(I-1,J+1)+S(I+1,J-1)+S(I-1,J-1))/4
	IF (MOD(I,2)/=0 .and. MOD(J,2)==0)	TECs(I,J)=(S(I+1,J)+S(I-1,J))/2
	IF (MOD(I,2)==0 .and. MOD(J,2)/=0)	TECs(I,J)=(S(I,J+1)+S(I,J-1))/2

	END IF

	WRITE(20,*),I,J,TECz(I,J),TECs(I,J),H(I,J),D(I,J),U(I,J),V(I,J),0

	END DO
	END DO

	END IF

!*********************************
!	Alternating Directin Implicit
!*********************************

!*********************************
!	Subroutine for Average Values
!*********************************

	ZZ=0
	VV=0

	CALL ZxShapeFunction (UDirection,Dx,Dy,M-1,N-1,Z,ZZ)

	CONTINUE

	CALL VxShapeFunction (UDirection,Dx,Dy,M-1,N-1,V,VV)

	CONTINUE

!*********************************
!	Starting 1st Half Time Step
!*********************************

	FUStar=FluxU

!*********************************
!	Starting Correction Loop
!*********************************

	DO Correct=1,10

	UU=0

	CALL FxShapeFunction (VDirection,Dx,Dy,M-1,N-1,FUStar,UU)

	CONTINUE

!*********************************
!	Updating Boundary Condition
!*********************************

	DO J=2,2*Y+1,2

	Z(2*X,J)=HWave/2*COS(-2*Pi*P*dT/Period)
	ZZ(2*X+1,J)=HWave/2*COS(-2*Pi*P*dT/Period)

	DO I=1,2*X+1

	Z(I,2)=Z(I,4)
	Z(I,2*Y)=Z(I,2*Y-2)

	IF (MOD(I,2)==0)	H(I,J)=D(I,J)+Z(I,J)
	IF (MOD(I,2)/=0)	H(I,J)=D(I,J)+ZZ(I,J)

	IF (H(I,J)<0.00) H(I,J)=0
			
!	U(I,2)=0.1*U(I,4)
!	U(I,2*Y)=0.1*U(I,2*Y-2)
!	V(I,2)=0.1*V(I,4)
!	V(I,2*Y)=0.1*V(I,2*Y-2)
!	V(I,2*Y-2)=0 !.1*V(I,2*Y-3)


	IF (H(I,J)==0)   WET(I,J)=0 
	IF (H(I,J)>0)    WET(I,J)=1
!	IF (WET(I,J)==0) H(I,J)=0
!	IF (WET(I,J)==0) S(I,J)=0
!	IF (WET(I,J)==0) U(I,J)=0
	IF (WET(I-2,J)==0.and.WET(I,J)==1) U(I-2,J)=U(I,J)
	IF (WET(I-2,J)==0.and.WET(I,J)==1) H(I-2,J)=H(I,J)
!	IF (WET(I,J)==0) V(I,J)=0
!	IF (H(I,J)==0)   Ah(I,J)=0  

	END DO
	END DO

!*********************************
!	Updating Temporary Coefficient
!*********************************

	DO J=2,2*Y+1,2
	DO I=2,2*X+1,2

	IF (H(I,J)>0.025) 	Cchezy(I,J)=-18*LOG10(Ke/(H(I,J)*12))
	IF (H(I,J)>0.025)	Ah(I,J)=Ce*H(I,J)*SQRT(Grav*(U(I+1,J)**2+V(I,J+1)**2))/Cchezy(I,J)
	Beta(I,J)=1
	IF (H(I,J)>0.025)	Beta(I,J)=1.+Grav/(Cchezy(I,J)**2+Karman**2)

	END DO
	END DO

	ADVv=0
	FVp=0
	FUp=0
	WINDx=0
	WINDy=0
	Coriolis=0
	DIFx=0
	DIFy=0
	Friction=0
	A1=0
	B1=0
	C1=0
	D1=0
	E1=0
	F1=0
	G1=0
	H1=0
	F0=0
	G0=0
	H0=0

!*********************************
!	Subroutine for Average Values
!*********************************

	AHH=0

	CALL ZxShapeFunction (UDirection,Dx,Dy,M-1,N-1,AH,AHH)

	CONTINUE

	CALL FxShapeFunction (VDirection,Dx,Dy,M-1,N-1,AHH,AHH)

	CONTINUE

!*********************************
!	Updating Advection Coefficient
!*********************************

	DO J=2,2*Y,2
	DO I=0,2*X,2

!	IF (WET(I,J)==1) THEN

	ADVv(I+1,J)=Beta(I+1,J)*(VV(I+1,J+1)*UU(I+1,J+1)-VV(I+1,J-1)*UU(I+1,J-1))/Dy

!	END IF

	END DO
	END DO

!*********************************
!	Updating Coriolis Coefficient
!*********************************

	DO J=2,2*Y,2
	DO I=2,2*X,2

!	IF (WET(I,J)==1) THEN

	FVp(1,J)=(FluxV(2,J+1)+FluxV(2,J-1))/4.
	FVp(I+1,J)=(FluxV(I+2,J+1)+FluxV(I+2,J-1)+FluxV(I,J-1)+FluxV(I,J+1))/4.

	Coriolis(1,J)=Cori*FVp(1,J)
	Coriolis(I+1,J)=Cori*FVp(I+1,J)

!	END IF
	  
	END DO
	END DO

!*********************************
!	Updating Wind Coefficient
!*********************************

	DO J=2,2*Y,2
	DO I=0,2*X,2

!	IF (WET(I,J)==1) THEN

	WINDx(I+1,J)=Twx(I+1,J)/Pw

!	END IF

	END DO
	END DO

!*********************************
!	Updating Friction Coefficient
!*********************************

	DO J=2,2*Y,2
	DO I=1,2*X,2

	IF (H(I+1,J)>0.025) THEN

	Friction(I,J)=Grav*SQRT(FUStar(I,J)**2+FVp(I,J)**2)/((Cchezy(I+1,J)*H(I,J))**2)

	END IF

	END DO
	END DO

!*********************************
!	Updating Vegetive Coefficient
!*********************************

	DO J=2,2*Y,2
	DO I=0,2*X,2

	IF (H(I+1,J)/=0) THEN

	Vegetation(I+1,J)=(Cvg(I+1,J)*Mvg(I+1,J)*Dvg(I+1,J))*SQRT(FUStar(I+1,J)**2+FVp(I+1,J)**2)/(H(I+1,J))

	END IF

	END DO
	END DO

!*********************************
!	Updating Diffusion Coefficient
!*********************************

	DO J=2,2*Y,2
	DO I=0,2*X,2

	IF (I==0 .or. J==2) THEN

	DIFy(I+1,J)=-(Ahh(I+1,J+1)*(FluxV(I+2,J+1))-Ahh(I+1,J-1)*(FluxV(I+2,J-1)))/(Pw*Dx*Dy)+ &
	&		  -(Ahh(I+1,J+1)*FUStar(I+1,J+2))/(Pw*Dy**2)

	ELSE

	DIFy(I+1,J)=-(Ahh(I+1,J+1)*(FluxV(I+2,J+1)-FluxV(I,J+1))-Ahh(I+1,J-1)*(FluxV(I+2,J-1)-FluxV(I,J-1)))/(Pw*Dx*Dy)+ &
	&		  -(Ahh(I+1,J+1)*FUStar(I+1,J+2)+Ahh(I+1,J-1)*FUStar(I+1,J-2))/(Pw*Dy**2)

	END IF

	END DO
	END DO

!*********************************
!	Matrix Coefficients Set Up [A]
!*********************************

	F0=0
	G0=0
	H0=0
	A1=0
	B1=0
	C1=0
	D1=0
	E1=0
	F1=0
	G1=0
	H1=0

	DO J=2,2*Y,2

	IF (H(1,J)/=0.) THEN

	F0(J)=2/Dt+Beta(1,J)*FUStar(1,J)/(Dx*H(1,J))+2*Ah(2,J)/(Pw*Dx**2)+ &
	&		(Ah(1,J+1)+Ah(1,J-1))/(Pw*Dy**2)+Friction(1,J)+Vegetation(1,J)
	G0(J)=+Grav*H(1,J)/Dx
	H0(J)=-2*Ah(2,J)/(Pw*Dx**2)

	ELSE 

	F0(J)=2/Dt+2*Ah(2,J)/(Pw*Dx**2)+(Ah(1,J+1)+Ah(1,J-1))/(Pw*Dy**2)+Friction(1,J)+Vegetation(1,J)

	END IF

	DO I=2,2*X+1,2

	A1(I,J)=-1/(Dx)
	B1(I,J)=+2/(Dt)
	C1(I,J)=+1/(Dx)

	IF (H(I-1,J)/=0 .and. H(I+1,J)/=0) THEN
	
	D1(I,J)=-Beta(I+1,J)*FUStar(I-1,J)/(Dx*H(I-1,J))-2*Ah(I,J)/(Pw*Dx**2)
	F1(I,J)=2/Dt+Beta(I+1,J)*FUStar(I+1,J)/(Dx*H(I+1,J))+2*(Ah(I+2,J)+Ah(I,J))/(Pw*Dx**2)+ &
	&		(Ahh(I+1,J+1)+Ahh(I+1,J-1))/(Pw*Dy**2)+Friction(I+1,J)+Vegetation(I+1,J)

	ELSE

	D1(I,J)=-2*Ah(I,J)/(Pw*Dx**2)
	F1(I,J)=2/Dt+2*(Ah(I+2,J)+Ah(I,J))/(Pw*Dx**2)+(Ahh(I+1,J+1)+Ahh(I+1,J-1))/(Pw*Dy**2)+Friction(I+1,J)+Vegetation(I+1,J)

	END IF

	E1(I,J)=-Grav*H(I+1,J)/Dx
	G1(I,J)=+Grav*H(I+1,J)/Dx
	H1(I,J)=-2*Ah(I+2,J)/(Pw*Dx**2)

	END DO
	END DO

!*********************************
!	Matrix Answers Set Up [B]
!*********************************

	MX1=0
	CX1=0

	DO J=2,2*Y,2

	MX0(J)=2*FluxU(1,J)/Dt-ADVv(1,J)+Coriolis(1,J)+WINDx(1,J)-DIFy(1,J)

	DO I=2,2*X,2

	MX1(I,J)=2*FluxU(I+1,J)/Dt-ADVv(I+1,J)+Coriolis(I+1,J)+WINDx(I+1,J)-DIFy(I+1,J)

	END DO 
	END DO

	DO J=2,2*Y,2
	DO I=2,2*X,2

	CX1(I,J)=2*Z(I,J)/Dt-(FluxV(I,J+1)-FluxV(I,J-1))/Dy

	END DO 
	END DO

!*********************************
!	Starting X-Direction Sweep
!*********************************

	DEALLOCATE (XX,BX,AX)
	ALLOCATE (XX(MAX(M-1,N-1)),BX(MAX(M-1,N-1)),AX(MAX(M-1,N-1),MAX(M-1,N-1)))

	DO J=2,2*Y,2

	AX=0
	BX=0

!*********************************
!	Matrix Set Up [A][X]=[B]
!*********************************
		
	DO I=1,2*X+1

	IF (I==1) 					BX(I)=MX0(J)
	IF (MOD(I,2)==0)			BX(I)=CX1(I,J)
	IF (MOD(I,2)/=0 .and. I/=1)	BX(I)=MX1(I-1,J)
	IF (I==2*X+1)				BX(I)=MX1(I-1,J)-G1(I-1,J)*ZZ(2*X+1,J)

	END DO

	DO I=1,M-1

	IF (I==1) 	AX(I,I)=F0(J)
	IF (I==1) 	AX(I,I+1)=G0(J)
	IF (I==1) 	AX(I,I+2)=H0(J)

	IF (MOD(I,2)/=0 .and. I/=1)				AX(I,I)=F1(I-1,J)
	IF (MOD(I,2)/=0 .and. I/=1 .and. I<M-1)	AX(I,I+1)=G1(I-1,J)
	IF (MOD(I,2)/=0 .and. I/=1 .and. I<M-1)	AX(I,I+2)=H1(I-1,J)
	IF (MOD(I,2)/=0 .and. I/=1)				AX(I,I-1)=E1(I-1,J)
	IF (MOD(I,2)/=0 .and. I/=1)				AX(I,I-2)=D1(I-1,J)

	IF (MOD(I,2)==0)	AX(I,I)=B1(I,J)
	IF (MOD(I,2)==0)	AX(I,I+1)=C1(I,J)
	IF (MOD(I,2)==0)	AX(I,I-1)=A1(I,J)

	END DO

!*********************************
!	Matrix Solving Method 
!*********************************

!	CALL GuassSeidel	(M-1,AX,BX,XX)
!	CALL GuassJordan	(M-1,AX,BX,XX)
	CALL PentaDiagonal	(M-1,AX,BX,XX)

	CONTINUE

!*********************************
!	Updating Temporary Component
!*********************************

	DO R=1,2*X+1

	IF (MOD(R,2)/=0)	UU(R,J)=XX(R)
	IF (MOD(R,2)==0)	ZZ(R,J)=XX(R)

	END DO

!*********************************
!	Ending X-Direction Sweep 
!*********************************

	END DO

!*********************************
!	Ending Correction Loop
!*********************************

	FUStar=(FluxU+UU)/2
	CorrectTest=ABS(FluxU-FUStar)

	IF(MAXVAL(CorrectTest)<1.E-6) EXIT

	END DO

!*********************************
!	Updating Temporary Component
!*********************************

	UDirection=0
	VDirection=0
	
	DO I=1,M
	DO J=1,N

	FluxU(I,J)=UU(I,J)
	IF (H(I,J)/=0)	U(I,J)=FluxU(I,J)/H(I,J)
	IF (H(I,J)<.025)	U(I,J)=0
	IF (U(I,J)<0)	UDirection(I,J)=1
	IF (V(I,J)<0)	VDirection(I,J)=1
	Z(I,J)=ZZ(I,J)							!*************

	END DO
	END DO

!*********************************
!	Oil Dynamic
!*********************************

!*********************************
!	Source & Sink of Oil 
!*********************************

!*********************************
!	Evaporation of Oil 
!*********************************

!	DO I=2,2*X,2
!	DO J=2,2*Y,2

!	EVAP(I,J)=0

!	END DO
!	END DO

!*********************************
!	Updating Temporary OIL Coefficient
!*********************************

	ADVv=0
	DIFy=0
	CXX=0
	SX=0
	A1=0
	B1=0
	C1=0
	D1=0
	E1=0
	F1=0
	G1=0
	H1=0
	CX1=0

	DO I=1,M
	DO J=1,N

	LANDA(I,J)=(Grav*S(I,J)**2)*(Pw-Po)*Po/(Pw*Cfilm)

	END DO
	END DO

	DO I=2,2*X,2
	DO J=2,2*Y,2

!	IF (WET(I,J)==1) THEN

	IF (J==2) THEN
	
	ADVv(I,J)=((S(I,J)+S(I,J+2))*(V(I,J+1)+Twy(I,J+1)/Cfilm)-(S(I,J))*(V(I,J-1)+Twy(I,J-1)/Cfilm))/(2*Dy)
	DIFy(I,J)=(LANDA(I,J+1)*(S(I,J+2)+S(I,J))-LANDA(I,J-1)*(S(I,J)))/(Dy**2)

	ELSE

	ADVv(I,J)=((S(I,J)+S(I,J+2))*(V(I,J+1)+Twy(I,J+1)/Cfilm)-(S(I,J)+S(I,J-2))*(V(I,J-1)+Twy(I,J-1)/Cfilm))/(2*Dy)
	DIFy(I,J)=(LANDA(I,J+1)*(S(I,J+2)-S(I,J))-LANDA(I,J-1)*(S(I,J)-S(I,J-2)))/(Dy**2)

	END IF

!	END IF

	END DO
	END DO

!*********************************
!	Matrix Coefficients Set Up [A]
!*********************************

	DO I=2,2*X,2
	DO J=2,2*Y,2

	A1(I,J)=-(U(I-1,J)+Twx(I-1,J)/Cfilm)/(2*Dx)-(LANDA(I-1,J))/(Dx**2)

	B1(I,J)=2/Dt+(U(I+1,J)+Twx(I+1,J)/Cfilm)/(2*Dx)-(U(I-1,J)+Twx(I-1,J)/Cfilm)/(2*Dx)+(LANDA(I+1,J)+LANDA(I-1,J))/(Dx**2)

	C1(I,J)=+(U(I+1,J)+Twx(I+1,J)/Cfilm)/(2*Dx)-(LANDA(I+1,J))/(Dx**2)

	END DO
	END DO

!*********************************
!	Matrix Answers Set Up [B]
!*********************************

	DO I=2,2*X,2
	DO J=2,2*Y,2

	CX1(I,J)=2*S(I,J)/Dt+DIFy(I,J)-ADVv(I,J)+2.*EVAP(I,J)/Dt

	END DO
	END DO

!*********************************
!	Starting X-Direction OIL Sweep
!*********************************

	DEALLOCATE (CXX,AAA,BBB,CCC)
	ALLOCATE (CXX(MAX(X,Y)),AAA(MAX(X,Y)),BBB(MAX(X,Y)),CCC(MAX(X,Y)))

	!!!!!!!!!!MIN & MAX!!!!!!!!!!!!!!!!!

	DO J=2,2*Y,2

!*********************************
!	Matrix Set Up [A][X]=[B]
!*********************************
	  	    
	DO I=2,2*X,2

	IF (I==2) THEN
	DDD(I/2)=CX1(I,J)-A1(2,J)*OILSource
	ELSE IF (I==2*X) THEN
	DDD(I/2)=CX1(I,J)-C1(2*X,J)*OILSink
	ELSE
	DDD(I/2)=CX1(I,J)
	END IF

	AAA(I/2)=A1(I,J)
	BBB(I/2)=B1(I,J)
	CCC(I/2)=C1(I,J)

	END DO

!*********************************
!	Matrix Solving Method 
!*********************************

	CALL Thomas (X,AAA,BBB,CCC,DDD,CXX)

	CONTINUE  
	  
!*********************************
!	Updating Temporary OIL Component
!*********************************
	  
	DO M=2,2*X,2

	IF (CXX(M/2)<0) CXX(M/2)=0

	SX(M,J)=CXX(M/2) 

	END DO

!*********************************
!	Ending X-Direction OIL Sweep 
!*********************************

	END DO

!*********************************
!	Updating Temporary OIL Component
!*********************************
	  
	S=SX
	
!*********************************
!	Subroutine for Average Values
!*********************************

	ZZ=0
	UU=0

	CALL ZyShapeFunction (VDirection,Dx,Dy,M-1,N-1,Z,ZZ)

	CONTINUE

	CALL UyShapeFunction (VDirection,Dx,Dy,M-1,N-1,U,UU)

	CONTINUE		
	
!*********************************
!	Starting 2nd Half Time Step
!*********************************

	FVStar=FluxV

!*********************************
!	Starting Correction Loop
!*********************************

	DO Correct=1,10

	VV=0

	CALL FyShapeFunction (UDirection,Dx,Dy,M-1,N-1,FVStar,VV)

	CONTINUE

!*********************************
!	Updating Boundary Condition
!*********************************

	DO J=1,2*Y+1

	Z(2*X,J)=HWave/2*COS(-2*Pi*P*dT/Period)
	ZZ(2*X+1,J)=HWave/2*COS(-2*Pi*P*dT/Period)

	DO I=2,2*X+1,2

	Z(I,2)=Z(I,4)
	Z(I,2*Y)=Z(I,2*Y-2)

	IF (MOD(J,2)==0)	H(I,J)=D(I,J)+Z(I,J)
	IF (MOD(J,2)/=0)	H(I,J)=D(I,J)+ZZ(I,J)

	IF (H(I,J)<0.00) H(I,J)=0

!	U(I,2)=0.1*U(I,4)
!	U(I,2*Y)=0.1*U(I,2*Y-2)
!	V(I,2)=0.1*V(I,4)
!	V(I,2*Y)=0.1*V(I,2*Y-2)
!	V(I,2*Y-2)=0 !.1*V(I,2*Y-3)


	IF (H(I,J)==0)   WET(I,J)=0 
	IF (H(I,J)>0)    WET(I,J)=1
!	IF (WET(I,J)==0) H(I,J)=0
!	IF (WET(I,J)==0) S(I,J)=0
!	IF (WET(I,J)==0) U(I,J)=0
	IF (WET(I-2,J)==0.and.WET(I,J)==1) V(I-2,J)=V(I,J)
	IF (WET(I-2,J)==0.and.WET(I,J)==1) H(I-2,J)=H(I,J)
!	IF (WET(I,J)==0) V(I,J)=0
!	IF (H(I,J)==0)   Ah(I,J)=0  

	END DO
	END DO

!*********************************
!	Updating Temporary Coefficient
!*********************************

	DO J=2,2*Y+1,2
	DO I=2,2*X+1,2

	IF (H(I,J)>0.025) 	Cchezy(I,J)=-18*LOG10(Ke/(H(I,J)*12))
	IF (H(I,J)>0.025)	Ah(I,J)=Ce*H(I,J)*SQRT(Grav*(U(I+1,J)**2+V(I,J+1)**2))/Cchezy(I,J)
	Beta(I,J)=1
	IF (H(I,J)>0.025)	Beta(I,J)=1.+Grav/(Cchezy(I,J)**2+Karman**2)

	END DO
	END DO

	ADVv=0
	FVp=0
	FUp=0
	WINDx=0
	WINDy=0
	Coriolis=0
	DIFx=0
	DIFy=0
	Friction=0
	A1=0
	B1=0
	C1=0
	D1=0
	E1=0
	F1=0
	G1=0
	H1=0
	F0=0
	G0=0
	H0=0

!*********************************
!	Subroutine for Average Values
!*********************************

	AHH=0

	CALL ZyShapeFunction (VDirection,Dx,Dy,M-1,N-1,AH,AHH)

	CONTINUE

	CALL FyShapeFunction (UDirection,Dx,Dy,M-1,N-1,AHH,AHH)

	CONTINUE

!*********************************
!	Updating Advection Coefficient
!*********************************

	DO J=0,2*Y,2
	DO I=2,2*X,2

!	IF (WET(I,J)==1) THEN

	ADVv(I,J+1)=Beta(I,J+1)*(UU(I+1,J+1)*VV(I+1,J+1)-UU(I-1,J+1)*VV(I-1,J+1))/Dx

!	END IF

	END DO
	END DO

!*********************************
!	Updating Coriolis Coefficient
!*********************************

	DO J=2,2*Y,2
	DO I=2,2*X,2

!	IF (WET(I,J)==1) THEN

	FUp(I,1)=(FluxU(I+1,2)+FluxU(I-1,2))/4.
	FUp(I,J+1)=(FluxU(I+1,J+2)+FluxU(I-1,J+2)+FluxU(I+1,J)+FluxU(I-1,J))/4.

	Coriolis(I,1)=-Cori*FUp(I,1)
	Coriolis(I,J+1)=-Cori*FUp(I,J+1)

!	END IF
	  
	END DO
	END DO

!*********************************
!	Updating Wind Coefficient
!*********************************

	DO J=0,2*Y,2
	DO I=2,2*X,2

!	IF (WET(I,J)==1) THEN

	WINDx(I,J+1)=Twx(I,J+1)/Pw

!	END IF

	END DO
	END DO

!*********************************
!	Updating Friction Coefficient
!*********************************

	DO J=1,2*Y,2
	DO I=2,2*X,2

	IF (H(I,J+1)>0.025) THEN

	Friction(I,J)=Grav*SQRT(FVStar(I,J)**2+FUp(I,J)**2)/((Cchezy(I,J+1)*H(I,J))**2)

	END IF

	END DO
	END DO

!*********************************
!	Updating Vegetive Coefficient
!*********************************

	DO J=0,2*Y,2
	DO I=2,2*X,2

	IF (H(I,J+1)/=0) THEN

	Vegetation(I,J+1)=(Cvg(I,J+1)*Mvg(I,J+1)*Dvg(I,J+1))*SQRT(FVStar(I,J+1)**2+FUp(I,J+1)**2)/(H(I,J+1))

	END IF

	END DO
	END DO

!*********************************
!	Updating Diffusion Coefficient
!*********************************

	DO J=0,2*Y,2
	DO I=2,2*X,2

	IF (J==0 .or. I==2) THEN

	DIFx(I,J+1)=-(Ahh(I+1,J+1)*(FluxU(I+1,J+2))-Ahh(I-1,J+1)*(FluxU(I-1,J+2)))/(Pw*Dx*Dy)+ &
	&		  -(Ahh(I+1,J+1)*FVStar(I+2,J+1))/(Pw*Dx**2)

	ELSE

	DIFx(I,J+1)=-(Ahh(I+1,J+1)*(FluxU(I+1,J+2)-FluxU(I+1,J))-Ahh(I-1,J+1)*(FluxU(I-1,J+2)-FluxU(I-1,J)))/(Pw*Dx*Dy)+ &
	&		  -(Ahh(I+1,J+1)*FVStar(I+2,J+1)+Ahh(I-1,J+1)*FVStar(I-2,J+1))/(Pw*Dx**2)

	END IF

	END DO
	END DO

!*********************************
!	Matrix Coefficients Set Up [A]
!*********************************

	F0=0
	G0=0
	H0=0
	A1=0
	B1=0
	C1=0
	D1=0
	E1=0
	F1=0
	G1=0
	H1=0

	DO I=2,2*X,2

	IF (H(I,1)/=0.) THEN

	F0(I)=2/Dt+Beta(I,1)*FVStar(I,1)/(Dy*H(I,1))+2*Ah(I,2)/(Pw*Dy**2)+ &
	&		(Ahh(I+1,1)+Ahh(I-1,1))/(Pw*Dx**2)+Friction(I,1)+Vegetation(I,1)
	G0(I)=+Grav*H(I,1)/Dy
	H0(I)=-2*Ah(I,2)/(Pw*Dy**2)

	ELSE 

	F0(I)=2/Dt+2*Ah(I,2)/(Pw*Dy**2)+(Ahh(I+1,1)+Ahh(I-1,1))/(Pw*Dx**2)+Friction(I,1)+Vegetation(I,1)

	END IF

	DO J=2,2*Y+1,2

	A1(I,J)=-1/(Dy)
	B1(I,J)=+2/(Dt)
	C1(I,J)=+1/(Dy)

	IF (H(I,J-1)/=0 .and. H(I,J+1)/=0) THEN
	
	D1(I,J)=-Beta(I,J+1)*FVStar(I,J-1)/(Dy*H(I,J-1))-2*Ah(I,J)/(Pw*Dy**2)
	F1(I,J)=2/Dt+Beta(I,J+1)*FVStar(I,J+1)/(Dy*H(I,J+1))+2*(Ah(I,J+2)+Ah(I,J))/(Pw*Dy**2)+ &
	&		(Ahh(I+1,J+1)+Ahh(I-1,J+1))/(Pw*Dx**2)+Friction(I,J+1)+Vegetation(I,J+1)

	ELSE

	D1(I,J)=-2*Ah(I,J)/(Pw*Dy**2)
	F1(I,J)=2/Dt+2*(Ah(I,J+2)+Ah(I,J))/(Pw*Dy**2)+(Ahh(I+1,J+1)+Ahh(I-1,J+1))/(Pw*Dx**2)+Friction(I,J+1)+Vegetation(I,J+1)

	END IF

	E1(I,J)=-Grav*H(I,J+1)/Dy
	G1(I,J)=+Grav*H(I,J+1)/Dy
	H1(I,J)=-2*Ah(I,J+2)/(Pw*Dy**2)

	END DO
	END DO

!*********************************
!	Matrix Answers Set Up [B]
!*********************************

	MX1=0
	CX1=0

	DO I=2,2*X,2

	MX0(I)=2*FluxV(I,1)/Dt-ADVv(I,1)+Coriolis(I,1)+WINDx(I,1)-DIFy(I,1)

	DO J=2,2*Y,2

	MX1(I,J)=2*FluxV(I,J+1)/Dt-ADVv(I,J+1)+Coriolis(I,J+1)+WINDx(I,J+1)-DIFy(I,J+1)

	END DO 
	END DO

	DO J=2,2*Y,2
	DO I=2,2*X,2

	CX1(I,J)=2*Z(I,J)/Dt-(FluxU(I+1,J)-FluxU(I-1,J))/Dx

	END DO 
	END DO

!*********************************
!	Starting Y-Direction Sweep
!*********************************

	DEALLOCATE (XX,BX,AX)
	ALLOCATE (XX(MIN(M-1,N-1)),BX(MIN(M-1,N-1)),AX(MIN(M-1,N-1),MIN(M-1,N-1)))

	DO I=2,2*X,2

	AX=0
	BX=0

!*********************************
!	Matrix Set Up [A][X]=[B]
!*********************************
		
	DO J=1,2*Y+1

	IF (J==1) 					BX(J)=MX0(I)
	IF (MOD(J,2)==0)			BX(J)=CX1(I,J)
	IF (MOD(J,2)/=0 .and. J/=1)	BX(J)=MX1(I,J-1)


	END DO

	DO J=1,N-1

	IF (J==1) 	AX(J,J)=F0(I)
	IF (J==1) 	AX(J,J+1)=G0(I)
	IF (J==1) 	AX(J,J+2)=H0(I)

	IF (MOD(J,2)/=0 .and. J/=1)				AX(J,J)=F1(I,J-1)
	IF (MOD(J,2)/=0 .and. J/=1 .and. J<N-1)	AX(J,J+1)=G1(I,J-1)
	IF (MOD(J,2)/=0 .and. J/=1 .and. J<N-1)	AX(J,J+2)=H1(I,J-1)
	IF (MOD(J,2)/=0 .and. J/=1)				AX(J,J-1)=E1(I,J-1)
	IF (MOD(J,2)/=0 .and. J/=1)				AX(J,J-2)=D1(I,J-1)

	IF (MOD(J,2)==0)	AX(J,J)=B1(I,J)
	IF (MOD(J,2)==0)	AX(J,J+1)=C1(I,J)
	IF (MOD(J,2)==0)	AX(J,J-1)=A1(I,J)

	END DO

!*********************************
!	Matrix Solving Method 
!*********************************

!	CALL GuassSeidel	(N-1,AX,BX,XX)
!	CALL GuassJordan	(N-1,AX,BX,XX)
	CALL PentaDiagonal	(N-1,AX,BX,XX)

	CONTINUE

!*********************************
!	Updating Temporary Component
!*********************************

	DO R=1,2*Y+1

	IF (MOD(R,2)/=0)	VV(I,R)=XX(R)
	IF (MOD(R,2)==0)	ZZ(I,R)=XX(R)

	END DO

!*********************************
!	Ending Y-Direction Sweep 
!*********************************

	END DO

!*********************************
!	Ending Correction Loop
!*********************************

	FVStar=(FluxV+VV)/2
	CorrectTest=ABS(FluxV-FVStar)

	IF(MAXVAL(CorrectTest)<1.E-6) EXIT

	END DO

!*********************************
!	Updating Temporary Component
!*********************************

	UDirection=0
	VDirection=0
	
	DO I=1,M
	DO J=1,N

	FluxV(I,J)=VV(I,J)
	IF (H(I,J)/=0)	V(I,J)=FluxV(I,J)/H(I,J)
	IF (H(I,J)<.025)	V(I,J)=0
	IF (V(I,J)<0)	VDirection(I,J)=1
	IF (U(I,J)<0)	UDirection(I,J)=1
	Z(I,J)=ZZ(I,J)							!*************

	END DO
	END DO	
	
!*********************************
!	Oil Dynamic
!*********************************

!*********************************
!	Source & Sink of Oil 
!*********************************

!*********************************
!	Evaporation of Oil 
!*********************************

!	DO I=2,2*X,2
!	DO J=2,2*Y,2

!	EVAP(I,J)=0

!	END DO
!	END DO

!*********************************
!	Updating Temporary OIL Coefficient
!*********************************

	ADVv=0
	DIFy=0
	CXX=0
	SX=0
	A1=0
	B1=0
	C1=0
	D1=0
	E1=0
	F1=0
	G1=0
	H1=0
	CX1=0

	DO I=1,M
	DO J=1,N

	LANDA(I,J)=(Grav*S(I,J)**2)*(Pw-Po)*Po/(Pw*Cfilm)

	END DO
	END DO

	DO I=2,2*X,2
	DO J=2,2*Y,2

!	IF (WET(I,J)==1) THEN

	IF (I==2) THEN
	
	ADVv(I,J)=((S(I,J)+S(I+2,J))*(U(I+1,J)+Twx(I+1,J)/Cfilm)-(S(I,J))*(U(I-1,J)+Twx(I-1,J)/Cfilm))/(2*Dx)
	DIFy(I,J)=(LANDA(I+1,J)*(S(I+2,J)-S(I,J))-LANDA(I-1,J)*(S(I,J)))/(Dx**2)

	ELSE

	ADVv(I,J)=((S(I,J)+S(I+2,J))*(U(I+1,J)+Twx(I+1,J)/Cfilm)-(S(I,J)+S(I-2,J))*(U(I-1,J)+Twx(I-1,J)/Cfilm))/(2*Dx)
	DIFy(I,J)=(LANDA(I+1,J)*(S(I+2,J)-S(I,J))-LANDA(I-1,J)*(S(I,J)-S(I-2,J)))/(Dx**2)

	END IF

!	END IF

	END DO
	END DO

!*********************************
!	Matrix Coefficients Set Up [A]
!*********************************

	DO I=2,2*X,2
	DO J=2,2*Y,2

	A1(I,J)=-(V(I,J-1)+Twy(I,J-1)/Cfilm)/(2*Dy)-(LANDA(I,J-1))/(Dy**2)

	B1(I,J)=2/Dt+(V(I,J+1)+Twy(I,J+1)/Cfilm)/(2*Dy)-(V(I,J-1)+Twy(I,J-1)/Cfilm)/(2*Dy)+(LANDA(I,J+1)+LANDA(I,J-1))/(Dy**2)

	C1(I,J)=+(V(I,J+1)+Twy(I,J+1)/Cfilm)/(2*Dy)-(LANDA(I,J+1))/(Dy**2)

	END DO
	END DO

!*********************************
!	Matrix Answers Set Up [B]
!*********************************

	DO I=2,2*X,2
	DO J=2,2*Y,2

	CX1(I,J)=2*S(I,J)/Dt+DIFy(I,J)-ADVv(I,J)+2.*EVAP(I,J)/Dt

	END DO
	END DO

!*********************************
!	Starting X-Direction OIL Sweep
!*********************************

	DEALLOCATE (CXX,AAA,BBB,CCC)
	ALLOCATE (CXX(MAX(X,Y)),AAA(MAX(X,Y)),BBB(MAX(X,Y)),CCC(MAX(X,Y)))

	!!!!!!!!!!MIN & MAX!!!!!!!!!!!!!!!!!

	DO I=2,2*X,2

!*********************************
!	Matrix Set Up [A][X]=[B]
!*********************************
	  	    
	DO J=2,2*Y,2

	IF (J==2) THEN
	DDD(J/2)=CX1(I,J)-A1(I,2)*OILSource
	ELSE IF (J==2*Y) THEN
	DDD(J/2)=CX1(I,J)-C1(I,2*Y)*OILSink
	ELSE
	DDD(J/2)=CX1(I,J)
	END IF

	AAA(J/2)=A1(I,J)
	BBB(J/2)=B1(I,J)
	CCC(J/2)=C1(I,J)

	END DO

!*********************************
!	Matrix Solving Method 
!*********************************

	CALL Thomas (Y,AAA,BBB,CCC,DDD,CXX)

	CONTINUE  
	  
!*********************************
!	Updating Temporary OIL Component
!*********************************
	  
	DO N=2,2*Y,2

	IF (CXX(N/2)<0) CXX(N/2)=0

	SX(I,N)=CXX(N/2) 

	END DO

!*********************************
!	Ending X-Direction OIL Sweep 
!*********************************

	END DO

!*********************************
!	Updating Temporary OIL Component
!*********************************
	  
	S=SX								

!*********************************
!	Ending the Main Loop
!*********************************

	PRINT *, P

	END DO

!	PRINT *, P

!*********************************
!	Subroutine Contains
!*********************************

	CONTAINS

!*********************************
!	Subroutine (Zx)-Shape Function
!*********************************

!	A Computational Algorithm For Estimating of Average Water Surface Elevation (Z) in X-Direction 

	SUBROUTINE ZxShapeFunction (K,DeltaX,DeltaY,P,Q,HH,HHX)

	INTEGER I,J,P,Q,K(P,Q)
	REAL  HH(P,Q),DeltaX,DeltaY,X,Y
	REAL  AX(P,Q),BX(P,Q),CX(P,Q),DX(P,Q),EE(P,Q),AY(P,Q),BY(P,Q),CY(P,Q),DY(P,Q)
	REAL, INTENT(OUT) :: HHX(P,Q)

	X=(DeltaX/2)
	Y=(DeltaY/2+DeltaY/2)

	DO J=6,Q-4,2
	DO I=6,P-4,2

	AX(I,J)=(HH(I+4,J)+HH(I-4,J)-4*(HH(I+2,J)+HH(I-2,J))+6*HH(I,J))/(24*DeltaX**4)
	AY(I,J)=(HH(I,J+4)+HH(I,J-4)-4*(HH(I,J+2)+HH(I,J-2))+6*HH(I,J))/(24*DeltaY**4)
	BX(I,J)=(HH(I+4,J)-HH(I-4,J)-2*(HH(I+2,J)-HH(I-2,J)))/(12*DeltaX**3)
	BY(I,J)=(HH(I,J+4)-HH(I,J-4)-2*(HH(I,J+2)-HH(I,J-2)))/(12*DeltaY**3)
	CX(I,J)=(-3*(HH(I+4,J)+HH(I-4,J))+36*(HH(I+2,J)+HH(I-2,J))-66*HH(I,J))/(48*DeltaX**2)
	CY(I,J)=(-3*(HH(I,J+4)+HH(I,J-4))+36*(HH(I,J+2)+HH(I,J-2))-66*HH(I,J))/(48*DeltaY**2)
	DX(I,J)=(-5*(HH(I+4,J)-HH(I-4,J))+34*(HH(I+2,J)-HH(I-2,J)))/(48*DeltaX)
	DY(I,J)=(-5*(HH(I,J+4)-HH(I,J-4))+34*(HH(I,J+2)-HH(I,J-2)))/(48*DeltaY)
	EE(I,J)=(9*(HH(I+4,J)+HH(I-4,J)+HH(I,J+4)+HH(I,J-4))-116*(HH(I+2,J)+HH(I-2,J)+HH(I,J+2)+HH(I,J-2))+2348*HH(I,J))/(1920)
	
	END DO
	END DO

	DO J=2,Q,2
	DO I=1,P,2

	HHX(3,J)=(HH(2,J)+HH(4,J))/2
	HHX(5,J)=(HH(4,J)+HH(6,J))/2
	HHX(P-2,J)=(HH(P-1,J)+HH(P-3,J))/2
	HHX(P-4,J)=(HH(P-3,J)+HH(P-5,J))/2
	IF (I/=1 .and. I/=P)	HHX(I,2)=(HH(I-1,2)+HH(I+1,2))/2
	IF (I/=1 .and. I/=P)	HHX(I,4)=(HH(I-1,4)+HH(I+1,4))/2

	IF (I/=1 .and. I/=P)	HHX(I,Q-1)=(HH(I-1,Q-1)+HH(I+1,Q-1))/2
	IF (I/=1 .and. I/=P)	HHX(I,Q-3)=(HH(I-1,Q-3)+HH(I+1,Q-3))/2

	HHX(1,J)=HH(2,J)
	HHX(P,J)=HH(P-1,J)

	END DO
	END DO

	DO J=6,Q-4,2
	DO I=4,P-4,2

	IF (K(I+1,J)==0 .and. I/=4)		HHX(I+1,J)=(Y*(AX(I,J)*(X**4)+BX(I,J)*(X**3)+CX(I,J)*(X**2)+DX(I,J)*(X))+AY(I,J)*(Y**5)/5+BY(I,J)*(Y**4)/4+CY(I,J)*(Y**3)/3+DY(I,J)*(Y**2)/2+EE(I,J)*Y)/(DeltaY)
	IF (K(I+1,J)==1 .and. I/=P-5)	HHX(I+1,J)=(Y*(AX(I+2,J)*(X**4)+BX(I+2,J)*(X**3)+CX(I+2,J)*(X**2)+DX(I+2,J)*(X))+AY(I+2,J)*(Y**5)/5+BY(I+2,J)*(Y**4)/4+CY(I+2,J)*(Y**3)/3+DY(I+2,J)*(Y**2)/2+EE(I+2,J)*Y)/(DeltaY)

	END DO
	END DO

	RETURN

	END SUBROUTINE ZxShapeFunction

!*********************************
!	Subroutine (Vx)-Shape Function
!*********************************

!	A Computational Algorithm For Estimating of Average Velocity (V) in X-Direction With U>0 and U<0

	SUBROUTINE VxShapeFunction (K,DeltaX,DeltaY,P,Q,HH,HHX)

	INTEGER I,J,P,Q,K(P,Q)
	REAL  HH(P,Q),DeltaX,DeltaY,X,Y
	REAL  AX(P,Q),BX(P,Q),CX(P,Q),DX(P,Q),EE(P,Q),AY(P,Q),BY(P,Q),CY(P,Q),DY(P,Q)
	REAL, INTENT(OUT) :: HHX(P,Q)

	X=(3*DeltaX/2+DeltaX/2)
	Y=(0)

	DO J=5,Q-4,2
	DO I=6,P-4,2

	AX(I,J)=(HH(I+4,J)+HH(I-4,J)-4*(HH(I+2,J)+HH(I-2,J))+6*HH(I,J))/(24*DeltaX**4)
	AY(I,J)=(HH(I,J+4)+HH(I,J-4)-4*(HH(I,J+2)+HH(I,J-2))+6*HH(I,J))/(24*DeltaY**4)
	BX(I,J)=(HH(I+4,J)-HH(I-4,J)-2*(HH(I+2,J)-HH(I-2,J)))/(12*DeltaX**3)
	BY(I,J)=(HH(I,J+4)-HH(I,J-4)-2*(HH(I,J+2)-HH(I,J-2)))/(12*DeltaY**3)
	CX(I,J)=(-3*(HH(I+4,J)+HH(I-4,J))+36*(HH(I+2,J)+HH(I-2,J))-66*HH(I,J))/(48*DeltaX**2)
	CY(I,J)=(-3*(HH(I,J+4)+HH(I,J-4))+36*(HH(I,J+2)+HH(I,J-2))-66*HH(I,J))/(48*DeltaY**2)
	DX(I,J)=(-5*(HH(I+4,J)-HH(I-4,J))+34*(HH(I+2,J)-HH(I-2,J)))/(48*DeltaX)
	DY(I,J)=(-5*(HH(I,J+4)-HH(I,J-4))+34*(HH(I,J+2)-HH(I,J-2)))/(48*DeltaY)
	EE(I,J)=(9*(HH(I+4,J)+HH(I-4,J)+HH(I,J+4)+HH(I,J-4))-116*(HH(I+2,J)+HH(I-2,J)+HH(I,J+2)+HH(I,J-2))+2348*HH(I,J))/(1920)
	
	END DO
	END DO

	DO J=1,Q,2
	DO I=1,P,2

	HHX(3,J)=(HH(2,J)+HH(4,J))/2
	HHX(5,J)=(HH(4,J)+HH(6,J))/2
	HHX(P-2,J)=(HH(P-1,J)+HH(P-3,J))/2
	HHX(P-4,J)=(HH(P-3,J)+HH(P-5,J))/2
	IF (I/=1 .and. I/=P)HHX(I,3)=(HH(I-1,3)+HH(I+1,3))/2
	IF (I/=1 .and. I/=P)HHX(I,Q-2)=(HH(I-1,Q-2)+HH(I+1,Q-2))/2
	HHX(1,J)=0
	HHX(P,J)=0
	HHX(I,1)=0
	HHX(I,Q)=0

	END DO
	END DO

	DO J=5,Q-4,2
	DO I=4,P-4,2

	IF (K(I+1,J-1)==0 .and. I/=4)	HHX(I+1,J)=(AX(I,J)*(X**5)/5+BX(I,J)*(X**4)/4+CX(I,J)*(X**3)/3+DX(I,J)*(X**2)/2+EE(I,J)*X)/(2*DeltaX)
	IF (K(I+1,J-1)==1 .and. I/=P-5)	HHX(I+1,J)=(AX(I+2,J)*(X**5)/5+BX(I+2,J)*(X**4)/4+CX(I+2,J)*(X**3)/3+DX(I+2,J)*(X**2)/2+EE(I+2,J)*X)/(2*DeltaX)

	END DO
	END DO

	RETURN

	END SUBROUTINE VxShapeFunction

!*********************************
!	Subroutine FluxUx Shape Function
!*********************************

!	A Computational Algorithm For Estimating of Average Flux (FluxU) in X-Direction With V>0 and V<0

	SUBROUTINE FxShapeFunction (K,DeltaX,DeltaY,P,Q,HH,HHX)

	INTEGER I,J,P,Q,K(P,Q)
	REAL  HH(P,Q),DeltaX,DeltaY,X,Y
	REAL  AX(P,Q),BX(P,Q),CX(P,Q),DX(P,Q),EE(P,Q),AY(P,Q),BY(P,Q),CY(P,Q),DY(P,Q)
	REAL, INTENT(OUT) :: HHX(P,Q)

	X=(DeltaX+DeltaX)
	Y=(DeltaY/2)

	DO J=6,Q-4,2
	DO I=5,P-4,2

	AX(I,J)=(HH(I+4,J)+HH(I-4,J)-4*(HH(I+2,J)+HH(I-2,J))+6*HH(I,J))/(24*DeltaX**4)
	AY(I,J)=(HH(I,J+4)+HH(I,J-4)-4*(HH(I,J+2)+HH(I,J-2))+6*HH(I,J))/(24*DeltaY**4)
	BX(I,J)=(HH(I+4,J)-HH(I-4,J)-2*(HH(I+2,J)-HH(I-2,J)))/(12*DeltaX**3)
	BY(I,J)=(HH(I,J+4)-HH(I,J-4)-2*(HH(I,J+2)-HH(I,J-2)))/(12*DeltaY**3)
	CX(I,J)=(-3*(HH(I+4,J)+HH(I-4,J))+36*(HH(I+2,J)+HH(I-2,J))-66*HH(I,J))/(48*DeltaX**2)
	CY(I,J)=(-3*(HH(I,J+4)+HH(I,J-4))+36*(HH(I,J+2)+HH(I,J-2))-66*HH(I,J))/(48*DeltaY**2)
	DX(I,J)=(-5*(HH(I+4,J)-HH(I-4,J))+34*(HH(I+2,J)-HH(I-2,J)))/(48*DeltaX)
	DY(I,J)=(-5*(HH(I,J+4)-HH(I,J-4))+34*(HH(I,J+2)-HH(I,J-2)))/(48*DeltaY)
	EE(I,J)=(9*(HH(I+4,J)+HH(I-4,J)+HH(I,J+4)+HH(I,J-4))-116*(HH(I+2,J)+HH(I-2,J)+HH(I,J+2)+HH(I,J-2))+2348*HH(I,J))/(1920)

	END DO
	END DO

	DO J=1,Q,2
	DO I=1,P,2

	HHX(I,3)=(HH(I,2)+HH(I,4))/2
	HHX(I,5)=(HH(I,4)+HH(I,6))/2
	HHX(I,Q-2)=(HH(I,Q-1)+HH(I,Q-3))/2
	HHX(I,Q-4)=(HH(I,Q-3)+HH(I,Q-5))/2
	IF (J/=1 .and. J/=Q)	HHX(3,J)=(HH(3,J-1)+HH(3,J+1))/2
	IF (J/=1 .and. J/=Q)	HHX(P-1,J)=(HH(P-1,J-1)+HH(P-1,J+1))/2
	HHX(1,J)=0
	HHX(P,J)=0
	HHX(I,1)=0
	HHX(I,Q)=0

	END DO
	END DO

	DO J=4,Q-4,2
	DO I=5,P-4,2

	IF (K(I-1,J+1)==0 .and. J/=4)	HHX(I,J+1)=(AX(I,J)*(X**5)/5+BX(I,J)*(X**4)/4+CX(I,J)*(X**3)/3+DX(I,J)*(X**2)/2+EE(I,J)*X+X*(AY(I,J)*(Y**4)+BY(I,J)*(Y**3)+CY(I,J)*(Y**2)+DY(I,J)*(Y)))/(2*DeltaX)
	IF (K(I-1,J+1)==1 .and. J/=Q-5)	HHX(I,J+1)=(AX(I,J+2)*(X**5)/5+BX(I,J+2)*(X**4)/4+CX(I,J+2)*(X**3)/3+DX(I,J+2)*(X**2)/2+EE(I,J+2)*X+X*(AY(I,J+2)*(Y**4)+BY(I,J+2)*(Y**3)+CY(I,J+2)*(Y**2)+DY(I,J+2)*(Y)))/(2*DeltaX)

	END DO
	END DO

	END SUBROUTINE FxShapeFunction

!*********************************
!	Subroutine (Zy)-Shape Function
!*********************************

!	A Computational Algorithm For Estimating of Average Water Surface Elevation (Z) in Y-Direction 

	SUBROUTINE ZyShapeFunction (K,DeltaX,DeltaY,P,Q,HH,HHX)

	INTEGER I,J,P,Q,K(P,Q)
	REAL  HH(P,Q),DeltaX,DeltaY,X,Y
	REAL  AX(P,Q),BX(P,Q),CX(P,Q),DX(P,Q),EE(P,Q),AY(P,Q),BY(P,Q),CY(P,Q),DY(P,Q)
	REAL, INTENT(OUT) :: HHX(P,Q)

	X=(DeltaX/2+DeltaX/2)
	Y=(DeltaY/2)

	DO J=6,Q-4,2
	DO I=6,P-4,2

	AX(I,J)=(HH(I+4,J)+HH(I-4,J)-4*(HH(I+2,J)+HH(I-2,J))+6*HH(I,J))/(24*DeltaX**4)
	AY(I,J)=(HH(I,J+4)+HH(I,J-4)-4*(HH(I,J+2)+HH(I,J-2))+6*HH(I,J))/(24*DeltaY**4)
	BX(I,J)=(HH(I+4,J)-HH(I-4,J)-2*(HH(I+2,J)-HH(I-2,J)))/(12*DeltaX**3)
	BY(I,J)=(HH(I,J+4)-HH(I,J-4)-2*(HH(I,J+2)-HH(I,J-2)))/(12*DeltaY**3)
	CX(I,J)=(-3*(HH(I+4,J)+HH(I-4,J))+36*(HH(I+2,J)+HH(I-2,J))-66*HH(I,J))/(48*DeltaX**2)
	CY(I,J)=(-3*(HH(I,J+4)+HH(I,J-4))+36*(HH(I,J+2)+HH(I,J-2))-66*HH(I,J))/(48*DeltaY**2)
	DX(I,J)=(-5*(HH(I+4,J)-HH(I-4,J))+34*(HH(I+2,J)-HH(I-2,J)))/(48*DeltaX)
	DY(I,J)=(-5*(HH(I,J+4)-HH(I,J-4))+34*(HH(I,J+2)-HH(I,J-2)))/(48*DeltaY)
	EE(I,J)=(9*(HH(I+4,J)+HH(I-4,J)+HH(I,J+4)+HH(I,J-4))-116*(HH(I+2,J)+HH(I-2,J)+HH(I,J+2)+HH(I,J-2))+2348*HH(I,J))/(1920)
	
	END DO
	END DO

	DO J=1,Q,2
	DO I=2,P,2

	HHX(I,3)=(HH(I,2)+HH(I,4))/2
	HHX(I,5)=(HH(I,4)+HH(I,6))/2
	HHX(I,Q-2)=(HH(I,Q-1)+HH(I,Q-3))/2
	HHX(I,Q-4)=(HH(I,Q-3)+HH(I,Q-5))/2
	IF (J/=1 .and. J/=Q)	HHX(2,J)=(HH(2,J-1)+HH(2,J+1))/2
	IF (J/=1 .and. J/=Q)	HHX(4,J)=(HH(4,J-1)+HH(4,J+1))/2

	IF (J/=1 .and. J/=Q)	HHX(P-1,J)=(HH(P-1,J-1)+HH(P-1,J+1))/2
	IF (J/=1 .and. J/=Q)	HHX(P-3,J)=(HH(P-3,J-1)+HH(P-3,J+1))/2

	HHX(I,1)=HH(I,2)
	HHX(I,Q)=HH(I,Q-1)

	END DO
	END DO

	DO J=4,Q-4,2
	DO I=6,P-4,2

	IF (K(I,J+1)==0 .and. J/=4)		HHX(I,J+1)=(AX(I,J)*(X**5)/5+BX(I,J)*(X**4)/4+CX(I,J)*(X**3)/3+DX(I,J)*(X**2)/2+EE(I,J)*X+X*(AY(I,J)*(Y**4)+BY(I,J)*(Y**3)+CY(I,J)*(Y**2)+DY(I,J)*(Y)))/(DeltaX)
	IF (K(I,J+1)==1 .and. J/=Q-5)	HHX(I,J+1)=(AX(I,J+2)*(X**5)/5+BX(I,J+2)*(X**4)/4+CX(I,J+2)*(X**3)/3+DX(I,J+2)*(X**2)/2+EE(I,J+2)*X+X*(AY(I,J+2)*(Y**4)+BY(I,J+2)*(Y**3)+CY(I,J+2)*(Y**2)+DY(I,J+2)*(Y)))/(DeltaX)

	END DO
	END DO

	RETURN

	END SUBROUTINE ZyShapeFunction

!*********************************
!	Subroutine (Uy)-Shape Function
!*********************************

!	A Computational Algorithm For Estimating of Average Velocity (U) in Y-Direction With V>0 and V<0

	SUBROUTINE UyShapeFunction (K,DeltaX,DeltaY,P,Q,HH,HHX)

	INTEGER I,J,P,Q,K(P,Q)
	REAL  HH(P,Q),DeltaX,DeltaY,X,Y
	REAL  AX(P,Q),BX(P,Q),CX(P,Q),DX(P,Q),EE(P,Q),AY(P,Q),BY(P,Q),CY(P,Q),DY(P,Q)
	REAL, INTENT(OUT) :: HHX(P,Q)

	X=(0)
	Y=(3*DeltaY/2+DeltaY/2)

	DO J=6,Q-4,2
	DO I=5,P-4,2

	AX(I,J)=(HH(I+4,J)+HH(I-4,J)-4*(HH(I+2,J)+HH(I-2,J))+6*HH(I,J))/(24*DeltaX**4)
	AY(I,J)=(HH(I,J+4)+HH(I,J-4)-4*(HH(I,J+2)+HH(I,J-2))+6*HH(I,J))/(24*DeltaY**4)
	BX(I,J)=(HH(I+4,J)-HH(I-4,J)-2*(HH(I+2,J)-HH(I-2,J)))/(12*DeltaX**3)
	BY(I,J)=(HH(I,J+4)-HH(I,J-4)-2*(HH(I,J+2)-HH(I,J-2)))/(12*DeltaY**3)
	CX(I,J)=(-3*(HH(I+4,J)+HH(I-4,J))+36*(HH(I+2,J)+HH(I-2,J))-66*HH(I,J))/(48*DeltaX**2)
	CY(I,J)=(-3*(HH(I,J+4)+HH(I,J-4))+36*(HH(I,J+2)+HH(I,J-2))-66*HH(I,J))/(48*DeltaY**2)
	DX(I,J)=(-5*(HH(I+4,J)-HH(I-4,J))+34*(HH(I+2,J)-HH(I-2,J)))/(48*DeltaX)
	DY(I,J)=(-5*(HH(I,J+4)-HH(I,J-4))+34*(HH(I,J+2)-HH(I,J-2)))/(48*DeltaY)
	EE(I,J)=(9*(HH(I+4,J)+HH(I-4,J)+HH(I,J+4)+HH(I,J-4))-116*(HH(I+2,J)+HH(I-2,J)+HH(I,J+2)+HH(I,J-2))+2348*HH(I,J))/(1920)
	
	END DO
	END DO

	DO J=1,Q,2
	DO I=1,P,2

	HHX(I,3)=(HH(I,2)+HH(I,4))/2
	HHX(I,5)=(HH(I,4)+HH(I,6))/2
	HHX(I,Q-2)=(HH(I,Q-1)+HH(I,Q-3))/2
	HHX(I,Q-4)=(HH(I,Q-3)+HH(I,Q-5))/2
	IF (J/=1 .and. J/=Q)HHX(3,J)=(HH(3,J-1)+HH(3,J+1))/2
	IF (J/=1 .and. J/=Q)HHX(P-2,J)=(HH(P-2,J-1)+HH(P-2,J+1))/2
	HHX(1,J)=0
	HHX(P,J)=0
	HHX(I,1)=0
	HHX(I,Q)=0

	END DO
	END DO

	DO J=4,Q-4,2
	DO I=5,P-4,2

	IF (K(I-1,J+1)==0 .and. J/=4)	HHX(I,J+1)=(AY(I,J)*(Y**5)/5+BY(I,J)*(Y**4)/4+CY(I,J)*(Y**3)/3+DY(I,J)*(Y**2)/2+EE(I,J)*Y)/(2*DeltaY)
	IF (K(I-1,J+1)==1 .and. J/=Q-5)	HHX(I,J+1)=(AY(I,J+2)*(Y**5)/5+BY(I,J+2)*(Y**4)/4+CY(I,J+2)*(Y**3)/3+DY(I,J+2)*(Y**2)/2+EE(I,J+2)*Y)/(2*DeltaY)

	END DO
	END DO

	RETURN

	END SUBROUTINE UyShapeFunction

!*********************************
!	Subroutine FluxVy Shape Function
!*********************************

!	A Computational Algorithm For Estimating of Average Flux (FluxV) in Y-Direction With U>0 and U<0

	SUBROUTINE FyShapeFunction (K,DeltaX,DeltaY,P,Q,HH,HHX)

	INTEGER I,J,P,Q,K(P,Q)
	REAL  HH(P,Q),DeltaX,DeltaY,X,Y
	REAL  AX(P,Q),BX(P,Q),CX(P,Q),DX(P,Q),EE(P,Q),AY(P,Q),BY(P,Q),CY(P,Q),DY(P,Q)
	REAL, INTENT(OUT) :: HHX(P,Q)

	X=(DeltaX/2)
	Y=(DeltaY+DeltaY)

	DO J=5,Q-4,2
	DO I=6,P-4,2

	AX(I,J)=(HH(I+4,J)+HH(I-4,J)-4*(HH(I+2,J)+HH(I-2,J))+6*HH(I,J))/(24*DeltaX**4)
	AY(I,J)=(HH(I,J+4)+HH(I,J-4)-4*(HH(I,J+2)+HH(I,J-2))+6*HH(I,J))/(24*DeltaY**4)
	BX(I,J)=(HH(I+4,J)-HH(I-4,J)-2*(HH(I+2,J)-HH(I-2,J)))/(12*DeltaX**3)
	BY(I,J)=(HH(I,J+4)-HH(I,J-4)-2*(HH(I,J+2)-HH(I,J-2)))/(12*DeltaY**3)
	CX(I,J)=(-3*(HH(I+4,J)+HH(I-4,J))+36*(HH(I+2,J)+HH(I-2,J))-66*HH(I,J))/(48*DeltaX**2)
	CY(I,J)=(-3*(HH(I,J+4)+HH(I,J-4))+36*(HH(I,J+2)+HH(I,J-2))-66*HH(I,J))/(48*DeltaY**2)
	DX(I,J)=(-5*(HH(I+4,J)-HH(I-4,J))+34*(HH(I+2,J)-HH(I-2,J)))/(48*DeltaX)
	DY(I,J)=(-5*(HH(I,J+4)-HH(I,J-4))+34*(HH(I,J+2)-HH(I,J-2)))/(48*DeltaY)
	EE(I,J)=(9*(HH(I+4,J)+HH(I-4,J)+HH(I,J+4)+HH(I,J-4))-116*(HH(I+2,J)+HH(I-2,J)+HH(I,J+2)+HH(I,J-2))+2348*HH(I,J))/(1920)

	END DO
	END DO

	DO J=1,Q,2
	DO I=1,P,2

	HHX(3,J)=(HH(2,J)+HH(4,J))/2
	HHX(5,J)=(HH(4,J)+HH(6,J))/2
	HHX(P-2,J)=(HH(P-1,J)+HH(P-3,J))/2
	HHX(P-4,J)=(HH(P-3,J)+HH(P-5,J))/2
	IF (I/=1 .and. I/=P)	HHX(I,3)=(HH(I-1,3)+HH(I+1,3))/2
	IF (I/=1 .and. I/=P)	HHX(I,Q-1)=(HH(I-1,Q-1)+HH(I+1,Q-1))/2
	HHX(1,J)=0
	HHX(P,J)=0
	HHX(I,1)=0
	HHX(I,Q)=0

	END DO
	END DO

	DO J=5,Q-4,2
	DO I=4,P-4,2

	IF (K(I+1,J-1)==0 .and. I/=4)	HHX(I+1,J)=(Y*(AX(I,J)*(X**4)+BX(I,J)*(X**3)+CX(I,J)*(X**2)+DX(I,J)*(X))+AY(I,J)*(Y**5)/5+BY(I,J)*(Y**4)/4+CY(I,J)*(Y**3)/3+DY(I,J)*(Y**2)/2+EE(I,J)*Y)/(2*DeltaY)
	IF (K(I+1,J-1)==1 .and. I/=P-5)	HHX(I+1,J)=(Y*(AX(I+2,J)*(X**4)+BX(I+2,J)*(X**3)+CX(I+2,J)*(X**2)+DX(I+2,J)*(X))+AY(I+2,J)*(Y**5)/5+BY(I+2,J)*(Y**4)/4+CY(I+2,J)*(Y**3)/3+DY(I+2,J)*(Y**2)/2+EE(I+2,J)*Y)/(2*DeltaY)

	END DO
	END DO

	END SUBROUTINE FyShapeFunction

!*********************************
!	Subroutine Tri-Diagonal
!*********************************

!	A Computational Algorithm For Solving Tri-Diagonal Matrix

	SUBROUTINE Thomas (Arrays,AAA,BBB,CCC,DDD,CXX)

	INTEGER M,R,Arrays
	REAL ZX
	REAL  AAA(Arrays),BBB(Arrays),CCC(Arrays),DDD(Arrays),C(Arrays),CX(Arrays)
	REAL, INTENT(OUT) :: CXX(Arrays)

	IF (BBB(1)==0.0) STOP 666
	ZX=BBB(1)
	C(1)=DDD(1)/ZX

	DO  M=2,Arrays
	CX(M)=CCC(M-1)/ZX
	ZX=BBB(M)-AAA(M)*CX(M)
	IF (ZX==0.0) EXIT
	C(M)=(DDD(M)-AAA(M)*C(M-1))/ZX
	END DO

	DO  M=Arrays-1,1,-1
          CXX(M)=C(M)-CX(M+1)*C(M+1)
	END DO	  

	RETURN

	END SUBROUTINE Thomas

!*********************************
!	Subroutine Penta-Diagonal
!*********************************

!	A Computational Algorithm For Solving Nearly Penta-Diagonal Matrix 

	SUBROUTINE PentaDiagonal (Arrays,A,BX,XX)

	INTEGER M,R,Arrays
	REAL  S,T
	REAL  BX(Arrays),A(Arrays,Arrays)
	REAL  C(Arrays),E(Arrays),F(Arrays)
	REAL  H(Arrays),V(Arrays),Z(Arrays),SIGMA(Arrays)
	REAL, INTENT(OUT) :: XX(Arrays)

	S=0.
	T=0.
	C(1)=A(1,1)
	F(1)=0.
	E(1)=A(1,2)
	V(1)=T
	H(1)=S/C(1)
	F(2)=A(2,1)/C(1)
	C(2)=A(2,2)-F(2)*E(1)
	V(2)=-F(2)*V(1)
	H(2)=-E(1)*H(1)/C(2)
	E(2)=A(2,3)-F(2)*A(1,3)

	DO M=3,Arrays-1
	F(M)=A(M,M-1)/C(M-1)-A(M,M-2)*E(M-2)/(C(M-2)*C(M-1))
	IF (M<Arrays-1) THEN
	E(M)=A(M,M+1)-F(M)*A(M-1,M+1)
	END IF
	C(M)=A(M,M)-F(M)*E(M-1)-A(M,M-2)*A(M-2,M)/C(M-2)
	END DO

	DO M=3,Arrays-1
	IF (M==Arrays-2) THEN
	V(M)=A(M,M+2)-A(M,M-2)*V(M-2)/C(M-2)-F(M)*V(M-1)
	H(M)=A(M+2,M)/C(M)-(A(M-2,M)*H(M-2)/C(M)+E(M-1)*H(M-1)/C(M))
	ELSE IF (M==Arrays-1) THEN
	V(M)=A(M,M+1)-A(M,M-2)*V(M-2)/C(M-2)-F(M)*V(M-1)
	H(M)=A(M+1,M)/C(M)-(A(M-2,M)*H(M-2)/C(M)+E(M-1)*H(M-1)/C(M))
	ELSE 
	V(M)=A(M,M-2)*V(M-2)/C(M-2)-F(M)*V(M-1)
	H(M)=-(A(M-2,M)*H(M-2)/C(M)+E(M-1)*H(M-1)/C(M))
	END IF
	END DO

	SIGMA=H*V
	C(Arrays)=A(Arrays,Arrays)-(SUM(SIGMA)-H(Arrays)*V(Arrays))

	Z(1)=BX(1)
	Z(2)=BX(2)-F(2)*Z(1)

	DO M=3,Arrays-1
	Z(M)=BX(M)-F(M)*Z(M-1)-A(M,M-2)*Z(M-2)/C(M-2)
	END DO

	SIGMA=H*Z
	Z(Arrays)=BX(Arrays)-(SUM(SIGMA)-H(Arrays)*Z(Arrays))

	XX(Arrays)=Z(Arrays)/C(Arrays)
	XX(Arrays-1)=(Z(Arrays-1)-V(Arrays-1)*XX(Arrays))/C(Arrays-1)
	XX(Arrays-2)=(Z(Arrays-2)-E(Arrays-2)*XX(Arrays-1)-V(Arrays-2)*XX(Arrays))/C(Arrays-2)

	DO M=Arrays-3,1,-1
	XX(M)=(Z(M)-E(M)*XX(M+1)-A(M,M+2)*XX(M+2)-V(M)*XX(Arrays))/C(M)
	END DO

	RETURN

	END SUBROUTINE PentaDiagonal

!*********************************
!	Subroutine Guass-Jordan
!*********************************

!	A Computational Algorithm For Solving Matrix 

	SUBROUTINE GuassJordan (Arrays,A,BX,XX)

	INTEGER I,J,L,K,NN,Arrays
	REAL A(Arrays,Arrays),BX(Arrays),C(Arrays,Arrays),SUM
	REAL, INTENT(OUT) :: XX(Arrays)

!	X=0
	NN=Arrays-1
	DO I=1,Arrays
	DO J=1,Arrays
	C(I,J)=A(I,J)
	END DO
	END DO

	DO K=1,NN
	L=K+1
	DO I=L,Arrays
	DO J=1,Arrays
	C(I,J)=A(I,J)-A(K,J)*A(I,K)/A(K,K)
	END DO
	BX(I)=BX(I)-BX(K)*A(I,K)/A(K,K)
	END DO
	DO I=1,Arrays
	DO J=1,Arrays
	A(I,J)=C(I,J)
	END DO 
	END DO
	END DO
	DO K=Arrays,1,-1.
	SUM=0.
	DO I=1,Arrays
	IF (I==K) EXIT
	SUM=SUM+A(K,I)*XX(I)
	END DO
	XX(K)=(BX(K)-SUM)/A(K,K)
	END DO

	RETURN

	END SUBROUTINE GuassJordan

!*********************************
!	Subroutine Guass-Seidel
!*********************************

!	A Computational Algorithm For Solving Matrix

	SUBROUTINE GuassSeidel (Arrays,A,BX,XX)

	INTEGER M,R,Iteration,Arrays
	REAL  DiffMax,Diff,Temp,SUM,Test
	REAL  BX(Arrays),A(Arrays,Arrays)
	REAL, INTENT(OUT) :: XX(Arrays)

	TEST=0.0001
	XX=0
	Iteration=0
	DO
	DiffMax=0
	Iteration=Iteration+1
	IF (Iteration > 100) STOP 1005
	DO M=1,Arrays
	Temp=XX(M)
	SUM=0
	DO R=1,Arrays
	IF (R /= M) THEN
	SUM=SUM+XX(R)*AX(M,R)
	END IF
	END DO
	XX(M)=(BX(M)-SUM)/AX(M,M)
	Diff=ABS(Temp-XX(M))
	IF (Diff > DiffMax) DiffMax=Diff
	END DO
	IF (DiffMax < Test) EXIT
	END DO

	RETURN

	END SUBROUTINE GuassSeidel

!*********************************
!	Ending the Program
!*********************************

	END PROGRAM Hydrodynamic
