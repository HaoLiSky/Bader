!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module for analyzing the charge
!
! By Andri Arnaldsson and Graeme Henkelman
! Last modified by GH on May 11, 2006
!-----------------------------------------------------------------------------------!

MODULE charge_mod

  USE kind_mod
  USE matrix_mod
  USE ions_mod
  IMPLICIT NONE

  TYPE :: charge_obj
    REAL(q2),ALLOCATABLE,DIMENSION(:,:,:) :: rho
    REAL(q2),DIMENSION(3,3) :: lat2car,car2lat
    REAL(q2),DIMENSION(-1:1,-1:1,-1:1) :: lat_dist,lat_i_dist
    REAL(q2),DIMENSION(3) :: org_lat,org_dir,org_car
    REAL(q2),DIMENSION(3) :: i_npts
    INTEGER,DIMENSION(3) :: npts
    INTEGER :: nrho
  END TYPE

  PRIVATE
  PUBLIC :: charge_obj
  PUBLIC :: rho_val,rho_grad
  PUBLIC :: pbc,dpbc_dir,dpbc
  PUBLIC :: is_max
  PUBLIC :: lat2car,car2lat,lat2dir,dir2lat

  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE assign_charge
  END INTERFACE

  CONTAINS

!-----------------------------------------------------------------------------------!
!  assign: copy one charge object to another
!-----------------------------------------------------------------------------------!

  SUBROUTINE assign_charge(chg1,chg2)

    TYPE(charge_obj), INTENT(INOUT) :: chg1
    TYPE(charge_obj), INTENT(IN) :: chg2
    
    chg1%lat2car=chg2%lat2car
    chg1%car2lat=chg2%car2lat
    chg1%org_lat=chg2%org_lat
    chg1%org_dir=chg2%org_dir
    chg1%org_car=chg2%org_car
    chg1%lat_dist=chg2%lat_dist
    chg1%lat_i_dist=chg2%lat_i_dist
    chg1%i_npts=chg2%i_npts
    chg1%npts=chg2%npts
    chg1%nrho=chg2%nrho

    ALLOCATE(chg1%rho(chg1%npts(1),chg1%npts(2),chg1%npts(3)))
    chg1%rho=chg2%rho

    END SUBROUTINE

!-----------------------------------------------------------------------------------!
!  rho_val:  Return the density at the point (p1,p2,p3) taking into account the
!    boundary conditions.  This function is used to address points outside the
!    charge density array without a bunch of if statements at the place the value
!    is needed.
!-----------------------------------------------------------------------------------!

  FUNCTION rho_val(chg,p1,p2,p3)

    TYPE(charge_obj) :: chg
    INTEGER,INTENT(IN) :: p1,p2,p3
    REAL(q2) :: rho_val

    INTEGER,DIMENSION(3) :: p
    INTEGER :: i

    p=(/p1,p2,p3/)
    DO i=1,3
      DO
        IF(p(i) >= 1) EXIT
        p(i)=p(i)+chg%npts(i)
      END DO
      DO
        IF(p(i) <= chg%npts(i)) EXIT
        p(i)=p(i)-chg%npts(i)
      END DO
    END DO

    rho_val=chg%rho(p(1),p(2),p(3))

  RETURN
  END FUNCTION rho_val

!-----------------------------------------------------------------------------------!
!  rho_grad:  Return the density and gradient at the point p
!-----------------------------------------------------------------------------------!

  FUNCTION rho_grad(chg,p,rho)

    TYPE(charge_obj) :: chg
    REAL(q2),DIMENSION(3),INTENT(IN) :: p
    REAL(q2),INTENT(OUT) :: rho
    REAL(q2),DIMENSION(3) :: rho_grad

    INTEGER :: p1,p2,p3
    REAL(q2),DIMENSION(3) :: rho_grad_lat
    REAL(q2) :: f1,f2,f3,g1,g2,g3
    REAL(q2) :: rho000,rho001,rho010,rho100,rho011,rho101,rho110,rho111
    REAL(q2) :: rho00_,rho01_,rho10_,rho11_
    REAL(q2) :: rho0__,rho1__,rho_0_,rho_1_,rho__0,rho__1
    REAL(q2) :: rho_00,rho_01,rho_10,rho_11

    p1=FLOOR(p(1))
    p2=FLOOR(p(2))
    p3=FLOOR(p(3))

    f1=p1-REAL(p1,q2)
    f2=p2-REAL(p2,q2)
    f3=p3-REAL(p3,q2)

    g1=1.0_q2-f1
    g2=1.0_q2-f2
    g3=1.0_q2-f3

    rho000=rho_val(chg,p1,p2,p3)
    rho001=rho_val(chg,p1,p2,p3+1)
    rho010=rho_val(chg,p1,p2+1,p3)
    rho100=rho_val(chg,p1+1,p2,p3)
    rho011=rho_val(chg,p1,p2+1,p3+1)
    rho101=rho_val(chg,p1+1,p2,p3+1)
    rho110=rho_val(chg,p1+1,p2+1,p3)
    rho111=rho_val(chg,p1+1,p2+1,p3+1)

    rho00_=rho000*g3+rho001*f3
    rho01_=rho010*g3+rho011*f3
    rho10_=rho100*g3+rho101*f3
    rho11_=rho110*g3+rho111*f3

    rho0__=rho00_*g2+rho01_*f2
    rho1__=rho10_*g2+rho11_*f2

    rho=rho0__*g1+rho1__*f1

! More work for gradients

    rho_0_=rho00_*g1+rho10_*f1
    rho_1_=rho01_*g1+rho11_*f1

    rho_00=rho000*g1+rho100*f1
    rho_01=rho001*g1+rho101*f1
    rho_10=rho010*g1+rho110*f1
    rho_11=rho011*g1+rho111*f1

    rho__0=rho_00*g2+rho_10*f2
    rho__1=rho_01*g2+rho_11*f2

    rho_grad_lat(1)=(rho1__-rho0__)
    rho_grad_lat(2)=(rho_1_-rho_0_)
    rho_grad_lat(3)=(rho__1-rho__0)

    CALL vector_matrix(rho_grad_lat,chg%car2lat,rho_grad)

  RETURN
  END FUNCTION rho_grad

!-----------------------------------------------------------------------------------!
! pbc: Wrap the point (p(1),p(2),p(3)) to the boundary conditions [0,pmax].
!-----------------------------------------------------------------------------------!

  SUBROUTINE pbc(p,pmax)

    INTEGER,DIMENSION(3),INTENT(INOUT) :: p
    INTEGER,DIMENSION(3),INTENT(IN) :: pmax

    INTEGER :: i

    DO i=1,3
      DO
        IF(p(i) > 0) EXIT
        p(i)=p(i)+pmax(i)
      END DO
      DO
        IF(p(i) <= pmax(i)) EXIT
        p(i)=p(i)-pmax(i)
      END DO
    END DO

  RETURN
  END SUBROUTINE pbc

!-----------------------------------------------------------------------------------!
! dpbc_dir:  Wrap the vector dr to the boundary conditions [-1/2,1/2].
!-----------------------------------------------------------------------------------!

  SUBROUTINE dpbc_dir(dr)

    REAL(q2),INTENT(INOUT),DIMENSION(3) :: dr

    INTEGER :: i

    DO i=1,3
      DO
        IF(dr(i) > -0.5_q2) EXIT
        dr(i)=dr(i)+1.0_q2
      END DO
      DO
        IF(dr(i) < 0.5_q2) EXIT
        dr(i)=dr(i)-1.0_q2
      END DO
    END DO
  RETURN
  END SUBROUTINE dpbc_dir

!-----------------------------------------------------------------------------------!
! dpbc:  Wrap the vector dr to the boundary conditions [-ngf/2,ngf/2].
!-----------------------------------------------------------------------------------!

  SUBROUTINE dpbc(dr,nf,nf_2)

    REAL(q2),INTENT(IN),DIMENSION(3) :: nf,nf_2
    REAL(q2),INTENT(INOUT),DIMENSION(3) :: dr

    INTEGER :: i

    DO i=1,3
      DO
        IF(dr(i) > -nf_2(i)) EXIT
        dr(i)=dr(i)+nf(i)
      END DO
      DO
        IF(dr(i) < nf_2(i)) EXIT
        dr(i)=dr(i)-nf(i)
      END DO
    END DO

  RETURN
  END SUBROUTINE dpbc

!-----------------------------------------------------------------------------------!
! is_max: return .true. if the grid point is a maximum of charge density
!-----------------------------------------------------------------------------------!

  FUNCTION is_max(chg,p)

    TYPE(charge_obj) :: chg
    INTEGER,DIMENSION(3),INTENT(IN) :: p
    LOGICAL :: is_max

    REAL(q2) :: rho
    INTEGER :: d1,d2,d3,p1,p2,p3
  
    is_max=.TRUE. 
    p1=p(1)
    p2=p(2)
    p3=p(3)
    rho=rho_val(chg,p1,p2,p3)
    DO d1=-1,1
      p1=p(1)+d1
      DO d2=-1,1
        p2=p(2)+d2
        DO d3=-1,1
          p3=p(3)+d3
          IF(rho_val(chg,p1,p2,p3) > rho) THEN
            is_max=.FALSE.
          END IF
        END DO
      END DO
    END DO

  RETURN 
  END FUNCTION is_max

!-----------------------------------------------------------------------------------!
! lat2dir: convert from a lattice grid point to a direct coordinate vector
!-----------------------------------------------------------------------------------!

  FUNCTION lat2dir(chg,p)

    TYPE(charge_obj) :: chg
    REAL(q2),DIMENSION(3),INTENT(IN) :: p
    REAL(q2),DIMENSION(3) :: lat2dir

    lat2dir=p-chg%org_lat
    lat2dir=lat2dir*chg%i_npts
    lat2dir=lat2dir+chg%org_dir

  RETURN
  END FUNCTION lat2dir

!-----------------------------------------------------------------------------------!
! lat2car: convert from a lattice grid point to a Cartesian coordinate vector
!-----------------------------------------------------------------------------------!

  FUNCTION lat2car(chg,p)

    TYPE(charge_obj) :: chg
    REAL(q2),DIMENSION(3),INTENT(IN) :: p
    REAL(q2),DIMENSION(3) :: lat2car,v

    v=p-chg%org_lat
    CALL matrix_vector(chg%lat2car,v,lat2car)
    lat2car=lat2car+chg%org_car

  RETURN
  END FUNCTION lat2car

!-----------------------------------------------------------------------------------!
! dir2lat: convert from a direct coordinate vector to a lattice grid point
!-----------------------------------------------------------------------------------!

  FUNCTION dir2lat(chg,p)

    TYPE(charge_obj) :: chg
    REAL(q2),DIMENSION(3),INTENT(IN) :: p
    REAL(q2),DIMENSION(3) :: dir2lat

    dir2lat=p-chg%org_dir
    dir2lat=dir2lat*REAL(chg%npts,q2)
    dir2lat=dir2lat+chg%org_dir

  RETURN
  END FUNCTION dir2lat

!-----------------------------------------------------------------------------------!
! car2lat: convert from a Cartesian coordinate vector to a lattice grid point
!-----------------------------------------------------------------------------------!

  FUNCTION car2lat(chg,p)

    TYPE(charge_obj) :: chg
    REAL(q2),DIMENSION(3),INTENT(IN) :: p
    REAL(q2),DIMENSION(3) :: car2lat,v

    v=p-chg%org_car
    CALL matrix_vector(chg%car2lat,v,car2lat)
    car2lat=car2lat+chg%org_lat

  RETURN
  END FUNCTION car2lat

!-----------------------------------------------------------------------------------!

END MODULE charge_mod
