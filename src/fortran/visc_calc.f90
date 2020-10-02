! This is a program to calculate the shear and bulk viscosity

Program visc_calc

      implicit none
      integer :: i, k, t, itmp
      integer :: ntimes
      real :: dt, volume
      
      real :: pxxo, pyyo, pzzo
      real :: pxyo, pxzo, pyzo
      real :: p_xyo, p_yzo

      real :: Pxy, Pxz, Pyz
      real :: P_xy, P_yz

      real, dimension(3) :: L
      real, dimension(6) :: press
      real, dimension(30000) :: P_av

      character(len=10) nfile, mol_name

      open(10,file='corr_calc.in',status='old')

      read(10,*)
      read(10,*) nfile
      read(10,*) 
      read(10,*) ntimes, dt
      read(10,*)
      read(10,*) volume
      read(10,*)
      read(10,*) mol_name
      close(10)

      L(1) = volume ** (1.0/3.0)
      L(2) = L(1)
      L(3) = L(1)
      
      open(11, file='pressures_out.log',status='old')
      
      ! Zero the results for NVE trajectory
      pxxo = 0.0; pyyo = 0.0; pzzo = 0.0
      pxyo = 0.0; pxzo = 0.0; pyzo = 0.0
      p_xyo = 0.0; p_yzo = 0.0
      
      read(11,*) pxxo, pyyo, pzzo, pxyo, pxzo, pyzo
      p_xyo = (pxxo-pyyo)/2.0
      p_yzo = (pyyo-pzzo)/2.0

      ! Loops over times
      do t = 1, ntimes
           read(11,*)  (press(k),k=1,6) 
           Pxy = pxyo * press(4) 
           Pxz = pxzo * press(5) 
           Pyz = pyzo * press(6)
           P_xy = p_xyo * (press(1) - press(2))/2.0
           P_yz = p_yzo * (press(2) - press(3))/2.0
           P_av(t) = (Pxy + Pxz + Pyz + P_xy + P_yz)/5.0
      end do
      
      open(15, file='shear_'//trim(nfile)//'_'//trim(mol_name)//'.dat')
      do t = 1, ntimes
        write (15,*) t, P_av(t)
      end do


end program visc_calc
