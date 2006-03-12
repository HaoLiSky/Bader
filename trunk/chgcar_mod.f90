!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!    Module for reading and writing VASP CHGCAR files
!
! By Andri Arnaldson and Graeme Henkelman
! Last modified by 
!-----------------------------------------------------------------------------------!

MODULE chgcar_mod
  USE kind_mod , ONLY : q1,q2
  USE ions_mod 
  USE charge_mod 
  USE matrix_mod , ONLY : transpose_matrix,matrix_vector
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_charge_chgcar,write_charge_chgcar
  CONTAINS

!-----------------------------------------------------------------------------------!
! read_charge_chgcar: Reads the charge density from a file in vasp format
!-----------------------------------------------------------------------------------!

  SUBROUTINE read_charge_chgcar(ions,chg,chargefile)

    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg
    CHARACTER(64) :: chargefile

    REAL(q2),DIMENSION(3,3) :: B
    REAL(q2),DIMENSION(3) :: box,v
    REAL(q2) :: side,vol,t
    INTEGER :: i
    INTEGER,DIMENSION(110) :: nionlist
    CHARACTER(LEN=7) :: text
    INTEGER :: nx,ny,nz

    OPEN(100,FILE=chargefile,STATUS='old',ACTION='read',BLANK='null',PAD='yes')
    WRITE(*,'(/,1A11,1A20)') 'OPEN ... ',chargefile
    elements=0
    WRITE(*,'(2x,A)') 'VASP-STYLE INPUT FILE'
    READ(100,'(/,1F20.16)') side
    READ(100,'(3F13.6)') (ions%lattice(i,1:3) , i=1,3)
    READ(100,'(110I4)') nionlist
    READ(100,*)
    DO i=1,110
     if(nionlist(i).eq.0) exit
    ENDDO
    ions%natypes=i-1
    ALLOCATE(ions%num_ions(ions%niontypes))
    DO i=1,ions%niontypes
      ions%num_ions(i)=nionlist(i)
    END DO
    ions%nions=SUM(nionlist)
    ions%lattice=side*ions%lattice
    CALL transpose_matrix(ions%lattice,B,3,3)
    wdim=ions%nions
    ALLOCATE(ions%r_car(ions%nions,3),ions%r_dir(ions%nions,3))
    DO i=1,ions%nions
!   Shouldn't r_dir be multiplied by side?
      READ(100,'(3(2X,1F8.6))') ions%r_dir(i,:)
      CALL matrix_vector(B,ions%r_dir(i,:),v,3,3)
      ions%r_car(i,:)=v
    END DO
    READ(100,*) 
    READ(100,*) chg%nxf,chg%nyf,chg%nzf
    WRITE(*,'(1A12,1I5,1A2,1I4,1A2,1I4)') 'FFT-grid: ',chg%nxf,'x',chg%nyf,'x',chg%nzf
    WRITE(*,'(2x,A,1A20)') 'CLOSE ... ', chargefile
    CLOSE(100)

  RETURN
  END SUBROUTINE read_charge_chgcar

!-----------------------------------------------------------------------------------!
! write_charge_chgcar: Write the charge density from a file in vasp format
!-----------------------------------------------------------------------------------!
    
  SUBROUTINE write_charge_chgcar(ions,chg,chargefile)
    
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg
    CHARACTER(64) :: chargefile

    OPEN(100,FILE=chargefile)
    WRITE(100,*)'Bader charge density file'
    WRITE(100,*)'1.00'
    WRITE(100,'(3F13.6)') (ions%lattice(i,1:3) , i=1,3)
    WRITE(100,'(110I4)') ions%num_ions
    WRITE(100,*)'DIRECT'
    WRITE(100,'(3(2X,1F8.6))') (ions%r_dir(i,:) , i=1,ions%nions)
    WRITE(100,*)
    WRITE(100,*) chg%nxf,chg%nyf,chg%nzf
    rho_tmp=0.0_q1
!    WHERE(bader_num == BaderCur) rho_tmp=rho
    WRITE(100,'(5E18.11)') (((chg%rho(nx,ny,nz),nx=1,chg%nxf),ny=1,chg%nyf),nz=1,chg%nzf)
    CLOSE(100)

  RETURN
  END SUBROUTINE write_charge_chgcar

!------------------------------------------------------------------------------------!

END MODULE chgcar_mod
