!   Program to calculate the individual MSDs and rotational correlation functions
!    and average them to get the total.
!

Program msd_rot_calc

    implicit none
    integer :: i, j, k, it
    integer :: ntimes
    integer :: natms, nmol, atms_per_mol
    real :: dt, eavg, volume
    real :: dp_O1, dp_O2, dOsq, d1sq, d2sq

    real, dimension(3) :: L
    real, dimension(500,3) :: rO_zero, r1_zero, r2_zero, rO, r1, r2
    real, dimension(500,3) :: rO_old, r1_old, r2_old, shiftO, shift1, shift2
    real, dimension(500,3) :: eO1, eO2, eO1_zero, eO2_zero
    real, dimension(0:30000,3) :: msd
    real, dimension(0:30000) :: c1, c2,c3

    real, dimension(5000) :: ener, dener

    character(len=10) nfile, mol_name
    character(len=2) ctmp

    atms_per_mol = 3  ! Set # of atoms per mol !!!!!!!!!!!!

    ! Read in the input parameters for the calculation
    open(10,file='corr_calc.in',status='old')
    
    read(10,*)
    read(10,*) nfile
    read(10,*) 
    read(10,*) ntimes, dt
    read(10,*) 
    read(10,*) volume
    read(10,*) 
    read(10,*) mol_name
    ! mol_type = 1 for water, acn
    ! mol_type = 2 for co2
    close(10)

    L(1)=volume ** (1.0/3.0)
    L(2)=L(1)
    L(3)=L(1)

    open(11,file='traj_'//trim(nfile)//'_'//trim(mol_name)//'.xyz',status='old')  !open traj file
     
    ! Zero the results for NVE trajectory
    msd = 0.0; c1 = 0.0; c2 = 0.0; c3 = 0.0

    ! Read in the first configuration and set it as the zero.
    read(11,*) natms
    read(11,*)
    nmol = natms/atms_per_mol
    do j = 1, nmol
        read(11,*) ctmp, (rO_zero(j,k),k=1,3)
        read(11,*) ctmp, (r1_zero(j,k),k=1,3)
        read(11,*) ctmp, (r2_zero(j,k),k=1,3)
        r1_zero(j,:) = r1_zero(j,:) - L(:)*anint(( r1_zero(j,:) - rO_zero(j,:) )/L(:))
        r2_zero(j,:) = r2_zero(j,:) - L(:)*anint(( r2_zero(j,:) - rO_zero(j,:) )/L(:))
    
        ! Define the bond unit vectors
        call bond_vec(nmol, rO_zero, r1_zero, eO1_zero)
        call bond_vec(nmol, rO_zero, r2_zero, eO2_zero)
        c1(0) = c1(0) + 2.0
        c2(0) = c2(0) + 2.0
        c3(0) = c3(0) + 2.0
    enddo

    ! Set the "old" coordinates for unwrapping to the initial ones
    rO_old = rO_zero; r1_old = r1_zero; r2_old = r2_zero
    shiftO = 0.0; shift1 = 0.0; shift2 = 0.0

    ! Loop over the timesteps in each trajectory
    do it = 2, ntimes

    ! Read in configuration at timestep it
        read(11,*)
        read(11,*)
        do j = 1, nmol
            read(11,*) ctmp, (rO(j,k),k=1,3)
            read(11,*) ctmp, (r1(j,k),k=1,3)
            read(11,*) ctmp, (r2(j,k),k=1,3)
                r1(j,:) = r1(j,:) - L(:)*anint(( r1(j,:) - rO(j,:) )/L(:))
                r2(j,:) = r2(j,:) - L(:)*anint(( r2(j,:) - rO(j,:) )/L(:))
         
            ! Define the bond unit vectors                                                                                                                   
            
            call bond_vec(nmol, rO, r1, eO1)
            call bond_vec(nmol, rO, r2, eO2)
         
         
            ! Calculate the contribution to the MSDs
            !   First m
            shiftO(j,:) =  shiftO(j,:) - L(:)*anint(( rO(j,:) - rO_old(j,:) )/L(:))
            shift1(j,:) =  shift1(j,:) - L(:)*anint(( r1(j,:) - r1_old(j,:) )/L(:))
            shift2(j,:) =  shift2(j,:) - L(:)*anint(( r2(j,:) - r2_old(j,:) )/L(:))

            !   Next the squared displacements
            dOsq = ( rO(j,1) + shiftO(j,1) - rO_zero(j,1) )**2 + ( rO(j,2) + shiftO(j,2) - rO_zero(j,2) )**2 &
                 + ( rO(j,3) + shiftO(j,3) - rO_zero(j,3) )**2
            d1sq = ( r1(j,1) + shift1(j,1) - r1_zero(j,1) )**2 + ( r1(j,2) + shift1(j,2) - r1_zero(j,2) )**2 &
                 + ( r1(j,3) + shift1(j,3) - r1_zero(j,3) )**2
            d2sq = ( r2(j,1) + shift2(j,1) - r2_zero(j,1) )**2 + ( r2(j,2) + shift2(j,2) - r2_zero(j,2) )**2 &
                 + ( r2(j,3) + shift2(j,3) - r2_zero(j,3) )**2

            !   Then add them to the MSDs
            msd(it-1,1) = msd(it-1,1) + dOsq
            msd(it-1,2) = msd(it-1,2) + d1sq
            msd(it-1,3) = msd(it-1,3) + d2sq
            
            ! Calculate the contribution to the correlation functions
            dp_O1 = eO1(j,1)*eO1_zero(j,1) + eO1(j,2)*eO1_zero(j,2) + eO1(j,3)*eO1_zero(j,3)
            dp_O2 = eO2(j,1)*eO2_zero(j,1) + eO2(j,2)*eO2_zero(j,2) + eO2(j,3)*eO2_zero(j,3)
            c1(it-1) = c1(it-1) + dp_O1 + dp_O2
            c2(it-1) = c2(it-1) + 0.5*(3.0*dp_O1**2 + 3.0*dp_O2**2 - 2.0)
            c3(it-1) = c3(it-1) + 0.5*(5.0*dp_O1**3 + 5.0*dp_O2**3 - 3.0*dp_O1 - 3.0*dp_O2)
            
        enddo ! End configuration read in

        ! Set the "old" coordinates for unwrapping to the last ones
        rO_old = rO; r1_old = r1; r2_old = r2

    enddo ! end loop over timesteps
    close(11)
    
    ! Write out the MSD, C1, C2
    open(21,file='c1_'//trim(nfile)//'_'//trim(mol_name)//'.dat')  !open C1(t) file for this trajectory
    open(22,file='c2_'//trim(nfile)//'_'//trim(mol_name)//'.dat')  !open C2(t) file for this trajectory
    open(24,file='c3_'//trim(nfile)//'_'//trim(mol_name)//'.dat')  !open C3(t) file for this trajectory
    open(23,file='msd_'//trim(nfile)//'_'//trim(mol_name)//'.dat') !open MSD(t) file for this trajectory

    do it = 0, ntimes - 1
        write(21,'(2F12.5)') real(it), c1(it)/real(2*nmol)
        write(22,'(2F12.5)') real(it), c2(it)/real(2*nmol)
        write(24,'(2F12.5)') real(it), c3(it)/real(2*nmol)
        write(23,'(5F12.5)') real(it), (msd(it,k)/real(nmol), k=1,3)
    enddo
    close(21)
    close(22)
    close(23)


End Program msd_rot_calc


! Subroutine to calculate the normalized bond vector from two atomic positions
!  e = rb - ra normed

Subroutine bond_vec(nmol, ra, rb, e)

implicit none
integer :: j, k, nmol
real, dimension(500,3) :: ra, rb, e
real, dimension(3) :: etmp
real :: norm

do j = 1, nmol
    etmp(:) = rb(j,:) - ra(j,:)
    norm = dot_product(etmp, etmp)

    e(j,:) = etmp(:)/sqrt(norm)

enddo


End Subroutine bond_vec
