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
	DIFy(I,J)=(LANDA(I,J+1)*(S(I,J+2)+S(I,J))-LANDA(I,J-1)*(S(I,J)-S(I,J-2)))/(Dy**2)

	END IF

!	END IF

	END DO
	END DO

!*********************************
!	Matrix Coefficients Set Up [A]
!*********************************

	DO I=2,2*X,2
	DO J=2,2*Y,2

	A1(I,J)=-(U(I-1,J)+Twy(I-1,J)/Cfilm)/(2*Dx)-(LANDA(I-1,J))/(Dx**2)

	B1(I,J)=2/Dt+(U(I+1,J)+Twy(I+1,J)/Cfilm)/(2*Dx)-(U(I-1,J)+Twy(I-1,J)/Cfilm)/(2*Dx)+(LANDA(I+1,J)+LANDA(I-1,J))/(Dx**2)

	C1(I,J)=+(U(I+1,J)+Twy(I+1,J)/Cfilm)/(2*Dx)-(LANDA(I+1,J))/(Dx**2)

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