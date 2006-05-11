!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module for analyzing the charge with a voronoi analysis
!
! By Andri Arnaldsson and Graeme Henkelman
! Last modified by 
!
!-----------------------------------------------------------------------------------!
MODULE voronoi_mod
  USE kind_mod
  USE matrix_mod
  USE ions_mod
  USE charge_mod
  IMPLICIT NONE

! Public, allocatable variables
  TYPE voronoi_obj
    REAL(q2),ALLOCATABLE,DIMENSION(:) :: vorchg
  END TYPE

  PRIVATE
  PUBLIC :: voronoi_obj
  PUBLIC :: voronoi

  CONTAINS

!-----------------------------------------------------------------------------------!
!  voronoi:  Calculate the charge density populations based upon Voronoi polyhedra.
!    In this scheme each element of charge density is associated with the atom that
!    it is closest to.
!-----------------------------------------------------------------------------------!

  SUBROUTINE voronoi(vor,ions,chg)

    TYPE(voronoi_obj) :: vor
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    REAL(q2),DIMENSION(3) :: r_lat,dr_lat,dr_car,ngf,ngf_2
    REAL(q2) :: dist,min_dist,shift
    INTEGER :: i,n1,n2,n3,closest,tenths_done,cr,count_max,t1,t2

    CALL system_clock(t1,cr,count_max)

    WRITE(*,'(/,2x,A)') 'CALCULATING VORONOI CHARGE DISTRIBUTION'
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'

    ALLOCATE(vor%vorchg(ions%nions))

    ngf=REAL(chg%npts,q2)
    ngf_2=ngf/2.0_q2

    vor%vorchg=0.0_q2
    tenths_done=0
    DO n1=1,chg%npts(1)
      r_lat(1)=REAL(n1,q2)
      IF ((n1*10/chg%npts(1)) > tenths_done) THEN
        tenths_done=(n1*10/chg%npts(1))
        WRITE(*,'(A,$)') '**'
      END IF
      DO n2=1,chg%npts(2)
        r_lat(2)=REAL(n2,q2)
        DO n3=1,chg%npts(3)
          r_lat(3)=REAL(n3,q2)
          closest=1
          dr_lat=r_lat-ions%r_lat(1,:)
          CALL dpbc(dr_lat,ngf,ngf_2)
          dr_car=lat2car(chg,dr_lat)
          min_dist=DOT_PRODUCT(dr_car,dr_car)
          DO i=2,ions%nions
            dr_lat=r_lat-ions%r_lat(i,:)
            CALL dpbc(dr_lat,ngf,ngf_2)
            dr_car=lat2car(chg,dr_lat)
            dist=DOT_PRODUCT(dr_car,dr_car)
            IF (dist < min_dist) THEN
              min_dist=dist
              closest=i
            END IF
          END DO
          vor%vorchg(closest)=vor%vorchg(closest)+rho_val(chg,n1,n2,n3)
        END DO
      END DO
    END DO

    ! Don't have this normalization for MONDO
    vor%vorchg(:)=vor%vorchg(:)/REAL(chg%nrho,q2)

    CALL system_clock(t2,cr,count_max)
    WRITE(*,'(/,1A12,1F6.2,1A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

  RETURN
  END SUBROUTINE voronoi

!-----------------------------------------------------------------------------------!

END MODULE voronoi_mod
