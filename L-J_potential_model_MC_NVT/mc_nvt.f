!mc simulation

    PROGRAM mc_nvt







    IMPLICIT NONE
    
    INTEGER iseed, equil, prod, nsamp, ii, icycl, ndispl, attempt, nacc, ncycl, nmoves, imove
   !iseed : seed random number generator
    !equil : number of Monte Carlo cycles during equilibration
    !prod  : number of Monte Carlo cycles during production
    !nsamp : number of Monte Carlo cycles between two sampling periods
    !ii    : 
    !icycl :
    !ndispl: number of attemps to displace a particle per MC cycle
    !attempt : number of attemps that have been performed to displace a particle
    !nacc : number of successful attemps to displace a particle
    !ncycl : number of monte carlo cycles
    !imove :
    
    DOUBLE PRECISION en, ent, vir, virt, dr
    ! en : total energy
    ! ent
    ! Dr         : maximum displacement

    WRITE (6, *) '**************** MC_NVT ***************'
    !initialize sysem
    CALL READDAT(equil, prod, nsamp, ndispl, dr, iseed) ! getting data
    nmoves = ndispl
    ! total energy of the system
    CALL TOTERG(en, vir)
    WRITE (6, 99001) en, vir
    !start mc cycle
    DO ii = 1, 2 !ii=1 is equilibration and ii=2 is production
        IF (ii.EQ.1)THEN
            ncycl = equil
            IF (ncycl.NE.0) WRITE (6, *) ' Start equilibration '
        ELSE
            IF (ncycl.NE.0) WRITE (6, *) ' Start production '
            ncycl = prod
        END IF
        attempt = 0
        nacc = 0
        !intialize the subroutine that adjust the maximum displacement
        CALL ADJUST(attempt, nacc, dr)
        DO icycl = 1, ncycl
            DO imove = 1, nmoves
                !attempt to displace a particle
                CALL MCMOVE(en, vir, attempt, nacc, dr, iseed)
            END DO
            IF (ii.EQ.2) THEN
                !sample averages
                IF (MOD(icycl,nsamp).EQ.0) CALL SAMPLE(icycl, en, vir)
            END IF
            IF (MOD(icycl,ncycl/5).EQ.0) THEN
                WRITE (6, *) '======>> Done ', icycl, ' out of ', ncycl
                !write intermidiate configuration to file
                CALL STORE(8, dr)
                !adjust maximum displacements
                CALL ADJUST(attempt, nacc, dr)
            END IF
        END DO
        IF (ncycl.NE.0) THEN
            IF (attempt.NE.0) WRITE (6, 99003) attempt, nacc, 
     &                               100.*FLOAT(nacc)/FLOAT(attempt)
!          ---test total energy
            CALL TOTERG(ent, virt)
            IF (ABS(ent-en).GT.1.D-6) THEN
               WRITE (6, *) 
     &                    ' ######### PROBLEMS ENERGY ################ '
            END IF
            IF (ABS(virt-vir).GT.1.D-6) THEN
               WRITE (6, *) 
     &                    ' ######### PROBLEMS VIRIAL ################ '
            END IF
            WRITE (6, 99002) ent, en, ent - en, virt, vir, virt - vir
         END IF
    END DO
    CALL STORE(21, dr)
    STOP



99001 FORMAT (' Total energy initial configuration: ', f12.5, /, 
     &        ' Total virial initial configuration: ', f12.5)
99002 FORMAT (' Total energy end of simulation    : ', f12.5, /, 
     &        '       running energy              : ', f12.5, /, 
     &        '       difference                  :  ', e12.5, /, 
     &        ' Total virial end of simulation    : ', f12.5, /, 
     &        '       running virial              : ', f12.5, /, 
     &        '       difference                  :  ', e12.5)
99003 FORMAT (' Number of att. to displ. a part.  : ', i10, /, 
     &        ' success: ', i10, '(= ', f5.2, '%)')
     

    END

