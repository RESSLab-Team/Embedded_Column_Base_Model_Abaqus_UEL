C --------------------------------------------------------------------------------------------------------------
C UEL subroutine to simulate embedded column base, ECB cyclic behavior with a rotational spring
C --------------------------------------------------------------------------------------------------------------
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4     PERIOD)
C     
      INCLUDE 'ABA_PARAM.INC'
C
      DOUBLE PRECISION Fy  
      DOUBLE PRECISION Fyn 
      DOUBLE PRECISION Fy2 
      DOUBLE PRECISION Fyn2
      DOUBLE PRECISION D  
      DOUBLE PRECISION Dn 
      DOUBLE PRECISION D2               
      DOUBLE PRECISION Dn2 
      DOUBLE PRECISION AK
      DOUBLE PRECISION AK2 
      DOUBLE PRECISION AK3 
      DOUBLE PRECISION Crp
      DOUBLE PRECISION Crn
      DOUBLE PRECISION px
      DOUBLE PRECISION py
      DOUBLE PRECISION Etot
      DOUBLE PRECISION DU_REPLACE
      DOUBLE PRECISION U_REPLACE 
C      
      PARAMETER ( ZERO = 0.D0, HALF = 0.5D0, ONE = 1.D0 )
C
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
     1     SVARS(NSVARS),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
     2     U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
     3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     5     JPROPS(*)
      DIMENSION SRESID(6)
C
C ##### INPUT PARAMETERS #####
C Get following information from .inp file
C
      Fy   =  PROPS(1) !Yield moment 
      Fyn  =  -Fy
      Fy2  =  PROPS(2) !Peak moment 
      Fyn2 =  -Fy2
      D    =  PROPS(3) !Yield rotation 
      Dn   =  -D
      D2   =  PROPS(4) !Rotation at max moment 
      Dn2  =  -D2
      AK   =  Fy/D ! Initial stiffness
      AK2  =  (Fy2-Fy)/(D2-D) ! Post-yield stiffness
      AK3  =  0.01*AK2 ! Stiffness after reaching peak moment (>0)
      px   = PROPS(5) !Pinching parameter px
      py   = PROPS(6) !Pinching parameter py
      Etot = PROPS(7) !Reference hysteretic energy dissipation capacity
C      
C -------------------------------------------------------
C        
      DO K1 = 1, NDOFEL                      
        SRESID(K1) = ZERO
        DO KRHS = 1, NRHS
          RHS(K1,KRHS) = ZERO
        END DO
        DO K2 = 1, NDOFEL
          AMATRX(K2,K1) = ZERO
        END DO
      END DO
C
      IF (LFLAGS(3).EQ.1) THEN
C
C ##### Obtain stiffness and force #####
C
C ----------------------------------- Notation ----------------------------------- 
C SVARS(#)  - solution dependent variable used to store the variable for the next step
C SVARS(1)  - store the calculated force
C SVARS(2)  - store the calculated stiffness
C SVARS(3)  - store the maximum force reached during the analysis (>Fy)
C SVARS(4)  - store the minimum force reached during the analysis (<Fyn)
C SVARS(5)  - store the maximum displacement reached during the analysis (>D)
C SVARS(6)  - store the minimum displacement reached during the analysis (<Dn)
C SVARS(7)  - store the initial displacement point to define the pinching in positive side 
C SVARS(8)  - store the initial displacement point to define the pinching in negative side
C SVARS(9)  - used to distinguish the step of unloading
C SVARS(10) - store the dissipated energy
C Force_pre - Force in previous step read from SVARS(1)
C Crp       - Maximum force reached until the previous step (>Fy) read from SVARS(3)
C Crn       - Minimum force reached until the previous step (<Fyn) read from SVARS(4)
C Fy        - Positive yield force (user input) 
C Fyn       - Negative yield force (=-Fy)
C Fy2       - Positive maximum force (user input)
C Fyn2      - Negative minimum force (=-Fy2)
C D         - Positive yield displacement (user input)
C Dn        - Negative yield displacement (=-D)
C D2        - Positive displacement at maximum force (user input) 
C Dn2       - Negative displacement at minimum force (=-Dn2)
C AK        - Initial stiffness
C AK2       - Post-yield stiffness
C AK3       - Stiffness after maximum force
C px        - Pinching parameter to define ratio of the displacements between part 1 and part 2 (user input)
C py        - Pinching parameter to define ratio of the forces between part 1 and part 2 (user input) 
C Etot      - Reference total energy used for the unloading stiffness degradation (user input) 
C ------------------------------------------------------------------------------------ 
C
C Read variables from previous step
        Force_pre = SVARS(1) 
        Crp = SVARS(3)
        Crn = SVARS(4)
C Subscript of DU and U depend on the aixs of the rotational spring
        DU_REPLACE = DU(8,1)-DU(4,1)
        U_REPLACE = U(8)-U(4) 
C         1 Displacement increment => 0
          IF ((DU_REPLACE) .GE. 0) THEN
C           Elastic
            IF ((Force_pre+SVARS(2)*((DU_REPLACE))) 
     +                .LE. Crp .and. Force_pre .GE. 0) THEN
C           1.1 Elastic before experiencing yielding
              IF (SVARS(3) .LE. Fy .and. SVARS(4) .GE. Fyn) THEN
                SVARS(2) =max(Fy,SVARS(3))/(max(D,SVARS(5))-SVARS(7))
                FORCE = Force_pre+SVARS(2)*((DU_REPLACE))
                SVARS(3)  = max(Fy,SVARS(3))
                SVARS(4)  = min(Fyn,SVARS(4))   
                SVARS(5)  = max(D,SVARS(5))
                SVARS(6) = min(Dn,SVARS(6))
                SVARS(7) = SVARS(7)   
                SVARS(8) = (U_REPLACE)-Force_pre/(AK*(1-(SVARS(10)
     +                    -HALF*(Force_pre**2)/AK)/Etot))   
                SVARS(9) = ZERO
                SVARS(10) = (DU_REPLACE)*(FORCE+Force_pre)*HALF
     +                       +SVARS(10)
C           Elastic after experiencing yielding                
              ELSE IF (SVARS(3) .GT. Fy .or. SVARS(4) .LT. Fyn) THEN
C               1.2 Pinching part 1              
                IF ((Force_pre+SVARS(2)*((DU_REPLACE))) .LE. 
     +             py*Crp) THEN
                  SVARS(2) = py*max(Fy,SVARS(3))/px/(max(D,SVARS(5))-
     +             SVARS(7))
                  FORCE = Force_pre+SVARS(2)*((DU_REPLACE))
                  SVARS(3)  = max(Fy,SVARS(3))
                  SVARS(4)  = min(Fyn,SVARS(4))   
                  SVARS(5)  = max(D,SVARS(5))
                  SVARS(6) = min(Dn,SVARS(6))
                  SVARS(7) = SVARS(7)   
                  SVARS(8) = (U_REPLACE)-Force_pre/(AK*(1-(SVARS(10)
     +                    -HALF*(Force_pre**2)/AK)/Etot))   
                  SVARS(9) = ZERO
                  SVARS(10) = (DU_REPLACE)*(FORCE+Force_pre)*HALF
     +                       +SVARS(10)
C               1.3 Pinching part 1 to part 2                  
                ELSE IF ((Force_pre+SVARS(2)*((DU_REPLACE))) .GE. 
     +             py*Crp .and. Force_pre .LE. py*Crp) THEN
                  FORCE = py*Crp+(((DU_REPLACE))-(py*Crp-
     +                Force_pre)/SVARS(2))*((ONE-py)*max(Fy,SVARS(3))/
     +                (ONE-px)*(max(D,SVARS(5))-SVARS(7)))                  
                  SVARS(2) = (ONE-py)*max(Fy,SVARS(3))/(ONE-px)/
     +             (max(D,SVARS(5))-SVARS(7))
                  SVARS(3)  = max(Fy,SVARS(3))
                  SVARS(4)  = min(Fyn,SVARS(4))   
                  SVARS(5)  = max(D,SVARS(5))
                  SVARS(6) = min(Dn,SVARS(6))
                  SVARS(7) = SVARS(7)   
                  SVARS(8) = (U_REPLACE)-Force_pre/(AK*(1-(SVARS(10)
     +                    -HALF*(Force_pre**2)/AK)/Etot))   
                  SVARS(9) = ZERO
                  SVARS(10) = (DU_REPLACE)*(FORCE+Force_pre)*HALF
     +                       +SVARS(10)
C               1.4 Pinching part 2                  
                ELSE IF ((Force_pre+SVARS(2)*((DU_REPLACE))) .GT. 
     +             py*Crp) THEN
                  SVARS(2) = (ONE-py)*max(Fy,SVARS(3))/(ONE-px)/
     +             (max(D,SVARS(5))-SVARS(7))
                  FORCE = Force_pre+SVARS(2)*((DU_REPLACE))
                  SVARS(3)  = max(Fy,SVARS(3))
                  SVARS(4)  = min(Fyn,SVARS(4))   
                  SVARS(5)  = max(D,SVARS(5))
                  SVARS(6) = min(Dn,SVARS(6))
                  SVARS(7) = SVARS(7)   
                  SVARS(8) = (U_REPLACE)-Force_pre/(AK*(1-(SVARS(10)
     +                    -HALF*(Force_pre**2)/AK)/Etot))   
                  SVARS(9) = ZERO
                  SVARS(10) = (DU_REPLACE)*(FORCE+Force_pre)*HALF
     +                       +SVARS(10)
                END IF
              END IF 
C           1.5 Yielding 1             
            ELSE IF ((Force_pre+SVARS(2)*((DU_REPLACE))) .GT. 
     +                  Crp .and. Force_pre .LE. 
     +                  Crp .and. Force_pre .GE. 0 .and.
     +                  Crp .LE. Fy2) THEN
              FORCE = Crp+(((DU_REPLACE))-(Crp-Force_pre)
     +                /SVARS(2))*AK2
              SVARS(2)  = AK2
              SVARS(3)  = max(FORCE,SVARS(3))
              SVARS(4)  = min(FORCE,SVARS(4))   
              SVARS(5)  = max((U_REPLACE),SVARS(5))
              SVARS(6) = min((U_REPLACE),SVARS(6))
              SVARS(7) = SVARS(7)
              SVARS(8) = SVARS(8)
              SVARS(9) = ZERO              
              SVARS(10) = (DU_REPLACE)*(FORCE+Force_pre)*HALF
     +                       +SVARS(10)
C           1.6 Hardening 1
            ELSE IF (Force_pre .GT. Crp .and. Crp+SVARS(2)*((DU_REPLACE
     +                )) .LE. Fy2) THEN 
              FORCE = Force_pre+AK2*(DU_REPLACE)
              SVARS(2)  = AK2
              SVARS(3)  = max(FORCE,SVARS(3))
              SVARS(4)  = min(FORCE,SVARS(4))   
              SVARS(5)  = max((U_REPLACE),SVARS(5))
              SVARS(6) = min((U_REPLACE),SVARS(6))
              SVARS(7) = SVARS(7)
              SVARS(8) = SVARS(8)
              SVARS(9) = ZERO
              SVARS(10) = (DU_REPLACE)*(FORCE+Force_pre)*HALF
     +                       +SVARS(10)
C           1.55 Yielding 2             
            ELSE IF ((Force_pre+SVARS(2)*((DU_REPLACE))) .GT. 
     +                  Crp .and. Force_pre .LE. 
     +                  Crp .and. Force_pre .GE. 0 .and.
     +                  Crp .GT. Fy2) THEN
              FORCE = Crp+(((DU_REPLACE))-(Crp-Force_pre)
     +                /SVARS(2))*AK3
              SVARS(2)  = AK3
              SVARS(3)  = max(FORCE,SVARS(3))
              SVARS(4)  = min(FORCE,SVARS(4))   
              SVARS(5)  = max((U_REPLACE),SVARS(5))
              SVARS(6) = min((U_REPLACE),SVARS(6))
              SVARS(7) = SVARS(7)
              SVARS(8) = SVARS(8)
              SVARS(9) = ZERO              
              SVARS(10) = (DU_REPLACE)*(FORCE+Force_pre)*HALF
     +                       +SVARS(10)
C           1.65 Hardening 2
            ELSE IF (Force_pre .GT. Crp .and. Crp .GT. Fy2) THEN
              FORCE = Force_pre+AK3*(DU_REPLACE)
              SVARS(2)  = AK3
              SVARS(3)  = max(FORCE,SVARS(3))
              SVARS(4)  = min(FORCE,SVARS(4))   
              SVARS(5)  = max((U_REPLACE),SVARS(5))
              SVARS(6) = min((U_REPLACE),SVARS(6))
              SVARS(7) = SVARS(7)
              SVARS(8) = SVARS(8)
              SVARS(9) = ZERO              
              SVARS(10) = (DU_REPLACE)*(FORCE+Force_pre)*HALF
     +                       +SVARS(10)
C           Unloading in Negative side  
            ELSE IF ((Force_pre)+(AK*(1-(SVARS(10)-HALF*(Force_pre**2)
     +                   /AK)/Etot))*((DU_REPLACE)) .LE. 0) THEN
C             1.7 first step of unloading   
              IF (SVARS(9) .EQ. ZERO) THEN
              	FORCE = Force_pre+(AK*(1-(SVARS(10)-HALF*(Force_pre**2) 
     +        	       /AK)/Etot))*((DU_REPLACE))
              	SVARS(2)  = HALF*((AK-AK*SVARS(10)/Etot)+((AK-AK*SVARS(10)/
     +        	            Etot)**2+4*HALF*(Force_pre**2)*AK/Etot)**(HALF))
              	SVARS(3)  = max(FORCE,SVARS(3))
              	SVARS(4)  = min(FORCE,SVARS(4))   
              	SVARS(5)  = max((U(7)-U(4)),SVARS(5))
              	SVARS(6) = min((U(7)-U(4)),SVARS(6))
              	SVARS(7) = (U_REPLACE)-Force_pre/(AK*(1-(SVARS(10)-HALF
     +        	            *(Force_pre**2)/AK)/Etot))   
              	SVARS(8) = (U_REPLACE)
              	SVARS(9) = ONE
              	SVARS(10) = (DU_REPLACE)*(FORCE+Force_pre)*HALF
     +        	               +SVARS(10)
C             1.8 following steps of unloading              
              ELSE IF (SVARS(9) .NE. ZERO) THEN
              	FORCE = Force_pre+SVARS(2)*((DU_REPLACE))
              	SVARS(2)  = SVARS(2) 
              	SVARS(3)  = max(FORCE,SVARS(3))
              	SVARS(4)  = min(FORCE,SVARS(4))   
              	SVARS(5)  = max((U_REPLACE),SVARS(5))
              	SVARS(6) = min((U_REPLACE),SVARS(6))
              	SVARS(7) = (U_REPLACE)-Force_pre/SVARS(2)   
              	SVARS(8) = (U_REPLACE)
              	SVARS(9) = HALF
              	SVARS(10) = (DU_REPLACE)*(FORCE+Force_pre)*HALF
     +        	               +SVARS(10)
              END IF              
C           1.9 Reversal loading from Negative side
            ELSE IF ((Force_pre)+AK*((DU_REPLACE)) 
     +        .GT. 0 .and. Force_pre .LT. 0) THEN
              FORCE = (((DU_REPLACE))-(0-Force_pre)/(AK*(1-
     +          (SVARS(10)-HALF*(Force_pre**2)/AK)/Etot)))*(py*
     +          max(Fy,SVARS(3))/px/(max(D,SVARS(5))-SVARS(7)))
              SVARS(2)  = py*max(Fy,SVARS(3))/px/(max(D,SVARS(5))-
     +                    SVARS(7))
              SVARS(3)  = max(FORCE,SVARS(3))
              SVARS(4)  = min(FORCE,SVARS(4))   
              SVARS(5)  = max((U_REPLACE),SVARS(5))
              SVARS(6) = min((U_REPLACE),SVARS(6))
              IF (SVARS(5) .LT. 2.0) THEN
                SVARS(7) = SVARS(7)
              ELSE IF (SVARS(5) .GE. 2.0) THEN
                SVARS(7) = (U_REPLACE)-Force_pre/(AK*(1-(SVARS(10)
     +                      -HALF*(Force_pre**2)/AK)/Etot))
              END IF   
              SVARS(8) = SVARS(8)
              SVARS(10) = (DU_REPLACE)*(FORCE+Force_pre)*HALF
     +                       +SVARS(10)
            END IF
C            
C         2 Displacement increment < 0                 
          ELSE IF ((DU_REPLACE) .LT. 0) THEN
C           Elastic
            IF ((Force_pre+SVARS(2)*((DU_REPLACE))) 
     +                .GE. Crn .and. Force_pre .LE. 0) THEN
C             2.1 Elastic before experiencing yielding     
              IF (SVARS(3) .LE. Fy .and. SVARS(4) .GE. Fyn) THEN
                SVARS(2) =min(Fyn,SVARS(4))/(min(Dn,SVARS(6))-
     +                    SVARS(8))
                FORCE = Force_pre+SVARS(2)*((DU_REPLACE))
                SVARS(3)  = max(Fy,SVARS(3))
                SVARS(4)  = min(Fyn,SVARS(4))   
                SVARS(5)  = max(D,SVARS(5))
                SVARS(6) = min(Dn,SVARS(6))
                SVARS(7) = (U_REPLACE)-Force_pre/(AK*(1-(SVARS(10)
     +                    -HALF*(Force_pre**2)/AK)/Etot))      
                SVARS(8) = SVARS(8)
                SVARS(9) = ZERO
                SVARS(10) = (DU_REPLACE)*(FORCE+Force_pre)*HALF
     +                       +SVARS(10)
C           Elastic after experiencing yielding      
              ELSE IF (SVARS(3) .GT. Fy .or. SVARS(4) .LT. Fyn) THEN
C               2.2 Pinching part 1              
                IF ((Force_pre+SVARS(2)*((DU_REPLACE))) .GE. 
     +             py*Crn) THEN
                  SVARS(2) = py*min(Fyn,SVARS(4))/px/
     +             (min(Dn,SVARS(6))-SVARS(8))
                  FORCE = Force_pre+SVARS(2)*((DU_REPLACE))
                  SVARS(3)  = max(Fy,SVARS(3))
                  SVARS(4)  = min(Fyn,SVARS(4))   
                  SVARS(5)  = max(D,SVARS(5))
                  SVARS(7) = (U_REPLACE)-Force_pre/(AK*(1-(SVARS(10)
     +                    -HALF*(Force_pre**2)/AK)/Etot))     
                  SVARS(8) = SVARS(8)
                  SVARS(9) = ZERO
                  SVARS(10) = (DU_REPLACE)*(FORCE+Force_pre)*HALF
     +                       +SVARS(10)
C               2.3 Pinching part 1 to part 2     
                ELSE IF ((Force_pre+SVARS(2)*((DU_REPLACE))) .LE. 
     +             py*Crn .and. Force_pre .GE. py*Crn) THEN
                  FORCE = py*Crn+(((DU_REPLACE))-(py*Crn-
     +                Force_pre)/SVARS(2))*((ONE-py)*min(Fyn,SVARS(4))
     +                /(ONE-px)*(min(Dn,SVARS(6))-SVARS(8)))                  
                  SVARS(2) = (ONE-py)*min(Fyn,SVARS(4))/(ONE-px)/
     +             (min(Dn,SVARS(6))-SVARS(8))
                  SVARS(3)  = max(Fy,SVARS(3))
                  SVARS(4)  = min(Fyn,SVARS(4))   
                  SVARS(5)  = max(D,SVARS(5))
                  SVARS(6) = min(Dn,SVARS(6))
                  SVARS(7) = (U_REPLACE)-Force_pre/(AK*(1-(SVARS(10)-HALF
     +                    *(Force_pre**2)/AK)/Etot))     
                  SVARS(8) = SVARS(8)
                  SVARS(9) = ZERO
                  SVARS(10) = (DU_REPLACE)*(FORCE+Force_pre)*HALF
     +                       +SVARS(10)
C               2.4 Pinching part 2     
                ELSE IF ((Force_pre+SVARS(2)*((DU_REPLACE))) .LT. 
     +             py*Crn) THEN
                  SVARS(2) = (ONE-py)*min(Fyn,SVARS(4))/(ONE-px)/
     +             (min(Dn,SVARS(6))-SVARS(8))
                  FORCE = Force_pre+SVARS(2)*((DU_REPLACE))
                  SVARS(3)  = max(Fy,SVARS(3))
                  SVARS(4)  = min(Fyn,SVARS(4))   
                  SVARS(5)  = max(D,SVARS(5))
                  SVARS(6) = min(Dn,SVARS(6))
                  SVARS(7) = (U_REPLACE)-Force_pre/(AK*(1-(SVARS(10)-HALF
     +                    *(Force_pre**2)/AK)/Etot))     
                  SVARS(8) = SVARS(8)
                  SVARS(9) = ZERO
                  SVARS(10) = (DU_REPLACE)*(FORCE+Force_pre)*HALF
     +                       +SVARS(10)
                END IF
              END IF 
C           2.5 Yielding 1             
            ELSE IF ((Force_pre+SVARS(2)*((DU_REPLACE))) .LT. 
     +                  Crn .and. Force_pre .GE. 
     +                  Crn .and. Force_pre .LE. 0 .and.
     +                  Crn .GE. Fyn2) THEN
              FORCE = Crn+(((DU_REPLACE))-(Crn-Force_pre)
     +                /SVARS(2))*AK2
              SVARS(2)  = AK2
              SVARS(3)  = max(FORCE,SVARS(3))
              SVARS(4)  = min(FORCE,SVARS(4))   
              SVARS(5)  = max((U_REPLACE),SVARS(5))
              SVARS(6) = min((U_REPLACE),SVARS(6))
              SVARS(7) = SVARS(7)
              SVARS(8) = SVARS(8)
              SVARS(9) = ZERO              
              SVARS(10) = (DU_REPLACE)*(FORCE+Force_pre)*HALF
     +                       +SVARS(10)
C           2.6 Hardening 1 
            ELSE IF (Force_pre .LT. Crn .and. Crn+SVARS(2)*((DU_REPLACE
     +                )) .GE. Fyn2) THEN 
              FORCE = Force_pre+AK2*(DU_REPLACE)
              SVARS(2)  = AK2
              SVARS(3)  = max(FORCE,SVARS(3))
              SVARS(4)  = min(FORCE,SVARS(4))   
              SVARS(5)  = max((U_REPLACE),SVARS(5))
              SVARS(6) = min((U_REPLACE),SVARS(6))
              SVARS(7) = SVARS(7)
              SVARS(8) = SVARS(8)
              SVARS(9) = ZERO
              SVARS(10) = (DU_REPLACE)*(FORCE+Force_pre)*HALF
     +                       +SVARS(10)
C           2.55 Yielding 2             
            ELSE IF ((Force_pre+SVARS(2)*((DU_REPLACE))) .LT. 
     +                  Crn .and. Force_pre .GE. 
     +                  Crn .and. Force_pre .LE. 0 .and.
     +                  Crn .LT. Fyn2) THEN
              FORCE = Crn+(((DU_REPLACE))-(Crn-Force_pre)
     +                /SVARS(2))*AK3
              SVARS(2)  = AK3
              SVARS(3)  = max(FORCE,SVARS(3))
              SVARS(4)  = min(FORCE,SVARS(4))   
              SVARS(5)  = max((U_REPLACE),SVARS(5))
              SVARS(6) = min((U_REPLACE),SVARS(6))
              SVARS(7) = SVARS(7)
              SVARS(8) = SVARS(8)
              SVARS(9) = ZERO              
              SVARS(10) = (DU_REPLACE)*(FORCE+Force_pre)*HALF
     +                       +SVARS(10)
C           2.65 Hardening 2               
            ELSE IF (Force_pre .LT. Crn .and. Crn .LT. Fyn2) THEN
              FORCE = Force_pre+AK3*(DU_REPLACE)
              SVARS(2)  = AK3
              SVARS(3)  = max(FORCE,SVARS(3))
              SVARS(4)  = min(FORCE,SVARS(4))   
              SVARS(5)  = max((U_REPLACE),SVARS(5))
              SVARS(6) = min((U_REPLACE),SVARS(6))
              SVARS(7) = SVARS(7)
              SVARS(8) = SVARS(8)
              SVARS(9) = ZERO              
              SVARS(10) = (DU_REPLACE)*(FORCE+Force_pre)*HALF
     +                       +SVARS(10)
C           Unloading in Positive side
            ELSE IF ((Force_pre)+(AK*(1-(SVARS(10)-HALF*(Force_pre**2)
     +                   /AK)/Etot))*((DU_REPLACE)) .GE. 0) THEN
C             2.7 first step of unloading        
              IF (SVARS(9) .EQ. ZERO) THEN
              	FORCE = Force_pre+(AK*(1-(SVARS(10)-HALF*(Force_pre**2)
     +        	       /AK)/Etot))*((DU_REPLACE))
              	SVARS(2)  = HALF*((AK-AK*SVARS(10)/Etot)+((AK-AK*SVARS(10)/
     +        	            Etot)**2+4*HALF*(Force_pre**2)*AK/Etot)**(HALF))
              	SVARS(3)  = max(FORCE,SVARS(3))
              	SVARS(4)  = min(FORCE,SVARS(4))   
              	SVARS(5)  = max((U_REPLACE),SVARS(5))
              	SVARS(6) = min((U_REPLACE),SVARS(6))
              	SVARS(7) = (U_REPLACE) 
              	SVARS(8) = (U_REPLACE)-Force_pre/(AK*(1-(SVARS(10)
     +        	            -HALF*(Force_pre**2)/AK)/Etot))   
              	SVARS(9) = ONE
              	SVARS(10) = (DU_REPLACE)*(FORCE+Force_pre)*HALF
     +        	               +SVARS(10)
C             2.8 following steps of unloading                 
              ELSE IF (SVARS(9) .NE. ZERO) THEN
              	FORCE = Force_pre+SVARS(2)*((DU_REPLACE))
              	SVARS(2)  = SVARS(2) 
              	SVARS(3)  = max(FORCE,SVARS(3))
              	SVARS(4)  = min(FORCE,SVARS(4))   
              	SVARS(5)  = max((U_REPLACE),SVARS(5))
              	SVARS(6) = min((U_REPLACE),SVARS(6))
              	SVARS(7) = (U_REPLACE)
              	SVARS(8) = (U_REPLACE)-Force_pre/SVARS(2)   
              	SVARS(9) = HALF
              	SVARS(10) = (DU_REPLACE)*(FORCE+Force_pre)*HALF
     +        	               +SVARS(10)
              END IF              
C           2.9 Reversal loading from Positive side
            ELSE IF ((Force_pre)+AK*((DU_REPLACE))
     +         .LT. 0 .and. Force_pre .GT. 0) THEN
              FORCE = (((DU_REPLACE))-(0-Force_pre)/(AK*(1-
     +                (SVARS(10)-HALF*(Force_pre**2)/AK)/Etot))
     +       )*(py*min(Fyn,SVARS(4))/px/(min(Dn,SVARS(6))-SVARS(8)))
              SVARS(2)  = py*min(Fyn,SVARS(4))/px/(min(Dn,SVARS(6))-
     +          SVARS(8))          
              SVARS(3)  = max(FORCE,SVARS(3))
              SVARS(4)  = min(FORCE,SVARS(4))   
              SVARS(5)  = max((U_REPLACE),SVARS(5))
              SVARS(6) = min((U_REPLACE),SVARS(6))
              SVARS(7) = SVARS(7)   
              IF (SVARS(5) .GT. 2.0) THEN
                SVARS(8) = SVARS(8)
              ELSE IF (SVARS(5) .LE. 2.0) THEN
                SVARS(8) = (U_REPLACE)-Force_pre/(AK*(1-(SVARS(10)
     +                      -HALF*(Force_pre**2)/AK)/Etot))
              END IF   
              SVARS(10) = (DU_REPLACE)*(FORCE+Force_pre)*HALF
     +                       +SVARS(10)
            END IF
          END IF
C Subscript of AMATRX, SRESID and RHS depend on the aixs of the rotational spring
C Provide stiffness and force 
        AMATRX(4,4) = SVARS(2)  
        AMATRX(8,8) = SVARS(2)  
        AMATRX(4,8) = -SVARS(2)  
        AMATRX(8,4) = -SVARS(2)
        SRESID(4) = -FORCE
        SRESID(8) =  FORCE
        RHS(4,1) = RHS(4,1)-SRESID(4)
        RHS(8,1) = RHS(8,1)-SRESID(8)
      ELSE IF (LFLAGS(3).EQ.2) THEN
C       Define the current stiffness matrix (AMATRX) 
C       Stiffness matrix
        AMATRX(4,4) =  SVARS(2)  
        AMATRX(8,8) =  SVARS(2)  
        AMATRX(4,8) = -SVARS(2)  
        AMATRX(8,4) = -SVARS(2)
      END IF
C Store variables to be used in next step      
      SVARS(1) = FORCE
      RETURN
      END         