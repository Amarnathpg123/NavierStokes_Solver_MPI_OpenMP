PROGRAM main
    ! This program numerically solves the 3D incompressible Navier-Stokes
	! on a Cubic Domain [0,2pi]x[0,2pi]x[0,2pi] using pseudo-spectral methods and
	! Implicit Midpoint rule timestepping. The numerical solution is compared to
	! an exact solution reported by Shapiro 
	!
	! Analytical Solution:
	!	u(x,y,z,t)=-0.25*(cos(x)sin(y)sin(z)+sin(x)cos(y)cos(z))exp(-t/Re)  
	!	v(x,y,z,t)= 0.25*(sin(x)cos(y)sin(z)-cos(x)sin(y)cos(z))exp(-t/Re)
	!	w(x,y,z,t)= 0.5*cos(x)cos(y)sin(z)exp(-t/Re)
    
    ! mpifort nse3dMpi.f90 -I./2decomp_fft/include -L./2decomp_fft/lib -l2decomp_fft -lm
    ! mpirun -np <no.of cores/mpi tasks> ./a.out

    USE decomp_2d
	USE decomp_2d_fft
	!USE decomp_2d_io
	
	IMPLICIT NONE	
	INCLUDE 'mpif.h'
    
    INTEGER(kind=4), PARAMETER 		:: Nx=128, Ny=128, Nz=128
    INTEGER(kind=4), PARAMETER 		:: Lx=1, Ly=1, Lz=1, Nt=20
    REAL(kind=8), PARAMETER			:: dt=0.2d0/Nt, Re=1.0d0, tol=0.1d0**10, theta=0.0d0
    REAL(kind=8), PARAMETER	        ::  pi=3.14159265358979323846264338327950288419716939937510d0
    REAL(kind=8), PARAMETER		    ::	ReInv=1.0d0/REAL(Re,kind(0d0)), dtInv=1.0d0/REAL(dt,kind(0d0)) 
    REAL(kind=8)									:: chg, factor, s, e, scalemodes=1.0d0/REAL(Nx*Ny*Nz,kind(0d0))

    REAL(kind=8), DIMENSION(:), ALLOCATABLE			:: x, y, z, time, mychg, allchg
    COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE		:: kx, ky, kz

    COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE	:: u, v, w, ux, uy, uz, vx, vy, vz, wx, wy, wz, uold, uxold, uyold, &
                                                        uzold, vold, vxold, vyold, vzold, wold, wxold, wyold, wzold,&
                                                        utemp, vtemp, wtemp, temp_r
                                                                    						
    COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE	:: uhat, vhat, what, rhsuhatfix, rhsvhatfix, rhswhatfix, nonlinuhat,&
                                                        nonlinvhat, nonlinwhat, phat,temp_c
    REAL(kind=8), DIMENSION(:,:,:), ALLOCATABLE 	:: realtemp

    
    TYPE(DECOMP_INFO)                       ::  decomp,sp
    INTEGER(kind=4)				:: i, j, k, n, t, allocatestatus, ind, ierr, p_row=0, p_col=0, numprocs, myid

    CALL MPI_INIT(ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr) 
	CALL decomp_2d_init(Nx,Ny,Nz,p_row,p_col)
	CALL decomp_info_init(Nx,Ny,Nz,decomp)
	CALL decomp_2d_fft_init
    
    IF (myid == 0) THEN
        PRINT *,'Grid:',Nx,'X',Ny,'Y',Nz,'Z'
        PRINT *,'dt:',dt	
    END IF

    ALLOCATE(x(1:Nx),y(1:Ny),z(1:Nz),time(1:Nt+1),mychg(1:3),allchg(1:3),&
				u(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),& 
 				v(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&
 				w(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&
 				ux(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&
				uy(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&
				uz(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&
				vx(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&
				vy(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&
				vz(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&
				wx(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&
				wy(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&
				wz(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&

				uold(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&
				uxold(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&
				uyold(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&
				uzold(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&
				vold(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&
				vxold(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&
				vyold(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&
				vzold(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&
				wold(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&
				wxold(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&
				wyold(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&
				wzold(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&

				utemp(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&
 				vtemp(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&
 				wtemp(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&
 				temp_r(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&
				kx(1:Nx),ky(1:Ny),kz(1:Nz),&

				uhat(decomp%zst(1):decomp%zen(1), decomp%zst(2):decomp%zen(2), decomp%zst(3):decomp%zen(3)),&
				vhat(decomp%zst(1):decomp%zen(1), decomp%zst(2):decomp%zen(2), decomp%zst(3):decomp%zen(3)),&
	 			what(decomp%zst(1):decomp%zen(1), decomp%zst(2):decomp%zen(2), decomp%zst(3):decomp%zen(3)),&
	 			rhsuhatfix(decomp%zst(1):decomp%zen(1), decomp%zst(2):decomp%zen(2), decomp%zst(3):decomp%zen(3)),&
 				rhsvhatfix(decomp%zst(1):decomp%zen(1), decomp%zst(2):decomp%zen(2), decomp%zst(3):decomp%zen(3)),&
 				rhswhatfix(decomp%zst(1):decomp%zen(1), decomp%zst(2):decomp%zen(2), decomp%zst(3):decomp%zen(3)),&
 				nonlinuhat(decomp%zst(1):decomp%zen(1), decomp%zst(2):decomp%zen(2), decomp%zst(3):decomp%zen(3)),&
 				nonlinvhat(decomp%zst(1):decomp%zen(1), decomp%zst(2):decomp%zen(2), decomp%zst(3):decomp%zen(3)),&
 				nonlinwhat(decomp%zst(1):decomp%zen(1), decomp%zst(2):decomp%zen(2), decomp%zst(3):decomp%zen(3)),&
 				phat(decomp%zst(1):decomp%zen(1), decomp%zst(2):decomp%zen(2), decomp%zst(3):decomp%zen(3)),&
 				temp_c(decomp%zst(1):decomp%zen(1), decomp%zst(2):decomp%zen(2), decomp%zst(3):decomp%zen(3)),&
				realtemp(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)),&
   					stat=AllocateStatus)	
	IF (AllocateStatus .ne. 0) STOP
	IF (myid == 0) PRINT *,'allocated space'
    
    CALL cpu_time(s)

    DO i=1,Nx/2+1
        kx(i)= cmplx(0.0d0,1.0d0)*REAL(i-1,kind(0d0))/Lx
    END DO
    kx(1+Nx/2)=0.0d0
    
    DO i = 1,Nx/2 -1
        kx(i+1+Nx/2)=-kx(1-i+Nx/2)
    END DO	
    
    ind=1
    DO i=-Nx/2,Nx/2-1
        x(ind)=2.0d0*pi*REAL(i,kind(0d0))*Lx/REAL(Nx,kind(0d0))
        ind=ind+1
    END DO
    
    DO j=1,Ny/2+1
        ky(j)= cmplx(0.0d0,1.0d0)*REAL(j-1,kind(0d0))/Ly
    END DO
    ky(1+Ny/2)=0.0d0
    
    DO j = 1,Ny/2 -1
        ky(j+1+Ny/2)=-ky(1-j+Ny/2)
    END DO	
    
    ind=1
    DO j=-Ny/2,Ny/2-1
        y(ind)=2.0d0*pi*REAL(j,kind(0d0))*Ly/REAL(Ny,kind(0d0))
        ind=ind+1
    END DO
    
    DO k=1,Nz/2+1
        kz(k)= cmplx(0.0d0,1.0d0)*REAL(k-1,kind(0d0))/Lz
    END DO
    kz(1+Nz/2)=0.0d0
    
    DO k = 1,Nz/2 -1
        kz(k+1+Nz/2)=-kz(1-k+Nz/2)
    END DO	
    
    ind=1
    DO k=-Nz/2,Nz/2-1
        z(ind)=2.0d0*pi*REAL(k,kind(0d0))*Lz/REAL(Nz,kind(0d0))
        ind=ind+1
    END DO
    
    time(1)=0.0d0
    factor=sqrt(3.0d0)
    DO k=decomp%xst(3),decomp%xen(3)
        DO j=decomp%xst(2),decomp%xen(2)
            DO i=decomp%xst(1),decomp%xen(1)
                u(i,j,k)=-0.5*( factor*cos(x(i))*sin(y(j))*sin(z(k))&
                        +sin(x(i))*cos(y(j))*cos(z(k)) )*exp(-(factor**2)*time(1)/Re)
            END DO
        END DO
    END DO

    DO k=decomp%xst(3),decomp%xen(3)
        DO j=decomp%xst(2),decomp%xen(2)
            DO i=decomp%xst(1),decomp%xen(1)
                v(i,j,k)=0.5*( factor*sin(x(i))*cos(y(j))*sin(z(k))&
                        -cos(x(i))*sin(y(j))*cos(z(k)) )*exp(-(factor**2)*time(1)/Re)
            END DO 
        END DO 
    END DO

    DO k=decomp%xst(3),decomp%xen(3)
        DO j=decomp%xst(2),decomp%xen(2)
            DO i=decomp%xst(1),decomp%xen(1)
                w(i,j,k)=cos(x(i))*cos(y(j))*sin(z(k))*exp(-(factor**2)*time(1)/Re)
            END DO 
        END DO 
    END DO
    
    CALL decomp_2d_fft_3d(u,uhat,DECOMP_2D_FFT_FORWARD)
	CALL decomp_2d_fft_3d(v,vhat,DECOMP_2D_FFT_FORWARD)
	CALL decomp_2d_fft_3d(w,what,DECOMP_2D_FFT_FORWARD)
    
    DO k=decomp%zst(3),decomp%zen(3)
        DO j=decomp%zst(2),decomp%zen(2)
            DO i=decomp%zst(1),decomp%zen(1)
                temp_c(i,j,k)=uhat(i,j,k)*kx(i)*scalemodes
            END DO 
        END DO 
    END DO

    CALL decomp_2d_fft_3d(temp_c,ux,DECOMP_2D_FFT_BACKWARD)
    
    DO k=decomp%zst(3),decomp%zen(3)
        DO j=decomp%zst(2),decomp%zen(2)
            DO i=decomp%zst(1),decomp%zen(1)
                temp_c(i,j,k)=uhat(i,j,k)*ky(j)*scalemodes
            END DO 
        END DO 
    END DO

    CALL decomp_2d_fft_3d(temp_c,uy,DECOMP_2D_FFT_BACKWARD)	
    
    DO k=decomp%zst(3),decomp%zen(3)
        DO j=decomp%zst(2),decomp%zen(2)
            DO i=decomp%zst(1),decomp%zen(1)
                temp_c(i,j,k)=uhat(i,j,k)*kz(k)*scalemodes
            END DO 
        END DO 
    END DO

    CALL decomp_2d_fft_3d(temp_c,uz,DECOMP_2D_FFT_BACKWARD)
    
    DO k=decomp%zst(3),decomp%zen(3)
        DO j=decomp%zst(2),decomp%zen(2)
            DO i=decomp%zst(1),decomp%zen(1)
                temp_c(i,j,k)=vhat(i,j,k)*kx(i)*scalemodes
            END DO 
        END DO 
    END DO

    CALL decomp_2d_fft_3d(temp_c,vx,DECOMP_2D_FFT_BACKWARD)

    DO k=decomp%zst(3),decomp%zen(3)
        DO j=decomp%zst(2),decomp%zen(2)
            DO i=decomp%zst(1),decomp%zen(1)
                temp_c(i,j,k)=vhat(i,j,k)*ky(j)*scalemodes
            END DO 
        END DO 
    END DO

    CALL decomp_2d_fft_3d(temp_c,vy,DECOMP_2D_FFT_BACKWARD)

    DO k=decomp%zst(3),decomp%zen(3)
        DO j=decomp%zst(2),decomp%zen(2)
            DO i=decomp%zst(1),decomp%zen(1)
                temp_c(i,j,k)=vhat(i,j,k)*kz(k)*scalemodes
            END DO 
        END DO 
    END DO

    CALL decomp_2d_fft_3d(temp_c,vz,DECOMP_2D_FFT_BACKWARD)
    
    DO k=decomp%zst(3),decomp%zen(3)
        DO j=decomp%zst(2),decomp%zen(2)
            DO i=decomp%zst(1),decomp%zen(1)
                temp_c(i,j,k)=what(i,j,k)*kx(i)*scalemodes
            END DO 
        END DO 
    END DO

    CALL decomp_2d_fft_3d(temp_c,wx,DECOMP_2D_FFT_BACKWARD)

    DO k=decomp%zst(3),decomp%zen(3)
        DO j=decomp%zst(2),decomp%zen(2)
            DO i=decomp%zst(1),decomp%zen(1)
                temp_c(i,j,k)=what(i,j,k)*ky(j)*scalemodes
            END DO 
        END DO 
    END DO

    CALL decomp_2d_fft_3d(temp_c,wy,DECOMP_2D_FFT_BACKWARD)

    DO k=decomp%zst(3),decomp%zen(3)
        DO j=decomp%zst(2),decomp%zen(2)
            DO i=decomp%zst(1),decomp%zen(1)
                temp_c(i,j,k)=what(i,j,k)*kz(k)*scalemodes
            END DO 
        END DO 
    END DO

    CALL decomp_2d_fft_3d(temp_c,wz,DECOMP_2D_FFT_BACKWARD)

    time(1)=0.0
    n=0

    DO n=1,Nt

        DO k=decomp%xst(3),decomp%xen(3)
            DO j=decomp%xst(2),decomp%xen(2)
                DO i=decomp%xst(1),decomp%xen(1)
                    uold(i,j,k)=u(i,j,k)
                    uxold(i,j,k)=ux(i,j,k)
                    uyold(i,j,k)=uy(i,j,k)
                    uzold(i,j,k)=uz(i,j,k)
                END DO 
            END DO 
        END DO

        DO k=decomp%xst(3),decomp%xen(3)
            DO j=decomp%xst(2),decomp%xen(2)
                DO i=decomp%xst(1),decomp%xen(1)
                    vold(i,j,k)=v(i,j,k)
                    vxold(i,j,k)=vx(i,j,k)
                    vyold(i,j,k)=vy(i,j,k)
                    vzold(i,j,k)=vz(i,j,k)
                END DO 
            END DO 
        END DO

        DO k=decomp%xst(3),decomp%xen(3)
            DO j=decomp%xst(2),decomp%xen(2)
                DO i=decomp%xst(1),decomp%xen(1)
                    wold(i,j,k)=w(i,j,k)
                    wxold(i,j,k)=wx(i,j,k)
                    wyold(i,j,k)=wy(i,j,k)
                    wzold(i,j,k)=wz(i,j,k)
                END DO 
            END DO 
        END DO

        DO k=decomp%zst(3),decomp%zen(3)
            DO j=decomp%zst(2),decomp%zen(2)
                DO i=decomp%zst(1),decomp%zen(1)
                    rhsuhatfix(i,j,k) = (dtInv+(0.5d0*ReInv)*&
                    (kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))*uhat(i,j,k) 
                END DO 
            END DO 
        END DO

        DO k=decomp%zst(3),decomp%zen(3)
            DO j=decomp%zst(2),decomp%zen(2)
                DO i=decomp%zst(1),decomp%zen(1)
                    rhsvhatfix(i,j,k) = (dtInv+(0.5d0*ReInv)*&
                    (kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))*vhat(i,j,k) 
                END DO 
            END DO 
        END DO

        DO k=decomp%zst(3),decomp%zen(3)
            DO j=decomp%zst(2),decomp%zen(2)
                DO i=decomp%zst(1),decomp%zen(1)
                    rhswhatfix(i,j,k) = (dtInv+(0.5d0*ReInv)*&
                    (kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))*what(i,j,k) 
                END DO 
            END DO 
        END DO
        
        chg=1
        DO WHILE (chg > tol)

            DO k=decomp%xst(3),decomp%xen(3)
                DO j=decomp%xst(2),decomp%xen(2)
                    DO i=decomp%xst(1),decomp%xen(1)
                        temp_r(i,j,k)=0.25d0*((u(i,j,k)+uold(i,j,k))*(ux(i,j,k)+uxold(i,j,k))&
                                        +(v(i,j,k)+vold(i,j,k))*(uy(i,j,k)+uyold(i,j,k))&
                                        +(w(i,j,k)+wold(i,j,k))*(uz(i,j,k)+uzold(i,j,k)))
                    END DO 
                END DO 
            END DO

            CALL decomp_2d_fft_3d(temp_r,nonlinuhat,DECOMP_2D_FFT_FORWARD)

            DO k=decomp%xst(3),decomp%xen(3)
                DO j=decomp%xst(2),decomp%xen(2)
                    DO i=decomp%xst(1),decomp%xen(1)
                        temp_r(i,j,k)=0.25d0*((u(i,j,k)+uold(i,j,k))*(vx(i,j,k)+vxold(i,j,k))&
                                        +(v(i,j,k)+vold(i,j,k))*(vy(i,j,k)+vyold(i,j,k))&
                                        +(w(i,j,k)+wold(i,j,k))*(vz(i,j,k)+vzold(i,j,k)))
                    END DO 
                END DO 
            END DO

            CALL decomp_2d_fft_3d(temp_r,nonlinvhat,DECOMP_2D_FFT_FORWARD)

            DO k=decomp%xst(3),decomp%xen(3)
                DO j=decomp%xst(2),decomp%xen(2)
                    DO i=decomp%xst(1),decomp%xen(1)
                        temp_r(i,j,k)=0.25d0*((u(i,j,k)+uold(i,j,k))*(wx(i,j,k)+wxold(i,j,k))&
                                        +(v(i,j,k)+vold(i,j,k))*(wy(i,j,k)+wyold(i,j,k))&
                                        +(w(i,j,k)+wold(i,j,k))*(wz(i,j,k)+wzold(i,j,k)))
                    END DO 
                END DO 
            END DO

            CALL decomp_2d_fft_3d(temp_r,nonlinwhat,DECOMP_2D_FFT_FORWARD)

            DO k=decomp%zst(3),decomp%zen(3)
                DO j=decomp%zst(2),decomp%zen(2)
                    DO i=decomp%zst(1),decomp%zen(1)
                        phat(i,j,k)=-1.0d0*( kx(i)*nonlinuhat(i,j,k)+ky(j)*nonlinvhat(i,j,k)&
                                    +kz(k)*nonlinwhat(i,j,k))/(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)+0.1d0**13)
                    END DO 
                END DO 
            END DO
    
            DO k=decomp%zst(3),decomp%zen(3)
                DO j=decomp%zst(2),decomp%zen(2)
                    DO i=decomp%zst(1),decomp%zen(1)
				        uhat(i,j,k)=(rhsuhatfix(i,j,k)-nonlinuhat(i,j,k)-kx(i)*phat(i,j,k))/&
							(dtInv-(0.5d0*ReInv)*(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k))) 
			        END DO 
                END DO 
            END DO

			DO k=decomp%zst(3),decomp%zen(3)
                DO j=decomp%zst(2),decomp%zen(2)
                    DO i=decomp%zst(1),decomp%zen(1)
				        vhat(i,j,k)=(rhsvhatfix(i,j,k)-nonlinvhat(i,j,k)-ky(j)*phat(i,j,k))/&
							(dtInv-(0.5d0*ReInv)*(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k))) 
			        END DO 
                END DO 
            END DO

			DO k=decomp%zst(3),decomp%zen(3)
                DO j=decomp%zst(2),decomp%zen(2)
                    DO i=decomp%zst(1),decomp%zen(1)
				        what(i,j,k)=(rhswhatfix(i,j,k)-nonlinwhat(i,j,k)-kz(k)*phat(i,j,k))/&
							(dtInv-(0.5d0*ReInv)*(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k))) 
			        END DO 
                END DO 
            END DO

			DO k=decomp%zst(3),decomp%zen(3)
                DO j=decomp%zst(2),decomp%zen(2)
                    DO i=decomp%zst(1),decomp%zen(1)
				        temp_c(i,j,k)=uhat(i,j,k)*kx(i)*scalemodes
			        END DO 
                END DO 
            END DO

			CALL decomp_2d_fft_3d(temp_c,ux,DECOMP_2D_FFT_BACKWARD)	

			DO k=decomp%zst(3),decomp%zen(3)
                DO j=decomp%zst(2),decomp%zen(2)
                    DO i=decomp%zst(1),decomp%zen(1)
				        temp_c(i,j,k)=uhat(i,j,k)*ky(j)*scalemodes
			        END DO 
                END DO 
            END DO

			CALL decomp_2d_fft_3d(temp_c,uy,DECOMP_2D_FFT_BACKWARD)	

			DO k=decomp%zst(3),decomp%zen(3)
                DO j=decomp%zst(2),decomp%zen(2)
                    DO i=decomp%zst(1),decomp%zen(1)
				        temp_c(i,j,k)=uhat(i,j,k)*kz(k)*scalemodes
			        END DO 
                END DO 
            END DO

			CALL decomp_2d_fft_3d(temp_c,uz,DECOMP_2D_FFT_BACKWARD)	

			DO k=decomp%zst(3),decomp%zen(3)
                DO j=decomp%zst(2),decomp%zen(2)
                    DO i=decomp%zst(1),decomp%zen(1)
				        temp_c(i,j,k)=vhat(i,j,k)*kx(i)*scalemodes
			        END DO 
                END DO 
            END DO
            
			CALL decomp_2d_fft_3d(temp_c,vx,DECOMP_2D_FFT_BACKWARD)	

			DO k=decomp%zst(3),decomp%zen(3)
                DO j=decomp%zst(2),decomp%zen(2)
                    DO i=decomp%zst(1),decomp%zen(1)
				        temp_c(i,j,k)=vhat(i,j,k)*ky(j)*scalemodes
			        END DO 
                END DO 
            END DO

			CALL decomp_2d_fft_3d(temp_c,vy,DECOMP_2D_FFT_BACKWARD)	

			DO k=decomp%zst(3),decomp%zen(3)
                DO j=decomp%zst(2),decomp%zen(2)
                    DO i=decomp%zst(1),decomp%zen(1)
				        temp_c(i,j,k)=vhat(i,j,k)*kz(k)*scalemodes
			        END DO 
                END DO 
            END DO

			CALL decomp_2d_fft_3d(temp_c,vz,DECOMP_2D_FFT_BACKWARD)	

			DO k=decomp%zst(3),decomp%zen(3)
                DO j=decomp%zst(2),decomp%zen(2)
                    DO i=decomp%zst(1),decomp%zen(1)
				        temp_c(i,j,k)=what(i,j,k)*kx(i)*scalemodes
			        END DO 
                END DO 
            END DO

			CALL decomp_2d_fft_3d(temp_c,wx,DECOMP_2D_FFT_BACKWARD)	

			DO k=decomp%zst(3),decomp%zen(3)
                DO j=decomp%zst(2),decomp%zen(2)
                    DO i=decomp%zst(1),decomp%zen(1)
				        temp_c(i,j,k)=what(i,j,k)*ky(j)*scalemodes
			        END DO 
                END DO 
            END DO

			CALL decomp_2d_fft_3d(temp_c,wy,DECOMP_2D_FFT_BACKWARD)	

			DO k=decomp%zst(3),decomp%zen(3)
                DO j=decomp%zst(2),decomp%zen(2)
                    DO i=decomp%zst(1),decomp%zen(1)
				        temp_c(i,j,k)=what(i,j,k)*kz(k)*scalemodes
			        END DO 
                END DO 
            END DO

			CALL decomp_2d_fft_3d(temp_c,wz,DECOMP_2D_FFT_BACKWARD)	

			DO k=decomp%xst(3),decomp%xen(3)
                DO j=decomp%xst(2),decomp%xen(2)
                    DO i=decomp%xst(1),decomp%xen(1)
				        utemp(i,j,k)=u(i,j,k)
			        END DO 
                END DO 
            END DO

			DO k=decomp%xst(3),decomp%xen(3)
                DO j=decomp%xst(2),decomp%xen(2)
                    DO i=decomp%xst(1),decomp%xen(1)
				        vtemp(i,j,k)=v(i,j,k)
			        END DO 
                END DO 
            END DO

			DO k=decomp%xst(3),decomp%xen(3)
                DO j=decomp%xst(2),decomp%xen(2)
                    DO i=decomp%xst(1),decomp%xen(1)
				        wtemp(i,j,k)=w(i,j,k)
			        END DO 
                END DO 
            END DO

			CALL decomp_2d_fft_3d(uhat,u,DECOMP_2D_FFT_BACKWARD)	
			CALL decomp_2d_fft_3d(vhat,v,DECOMP_2D_FFT_BACKWARD)	
			CALL decomp_2d_fft_3d(what,w,DECOMP_2D_FFT_BACKWARD)

			DO k=decomp%xst(3),decomp%xen(3)
                DO j=decomp%xst(2),decomp%xen(2)
                    DO i=decomp%xst(1),decomp%xen(1)
				        u(i,j,k)=u(i,j,k)*scalemodes
			        END DO 
                END DO 
            END DO

			DO k=decomp%xst(3),decomp%xen(3)
                DO j=decomp%xst(2),decomp%xen(2)
                    DO i=decomp%xst(1),decomp%xen(1)
				        v(i,j,k)=v(i,j,k)*scalemodes
			        END DO 
                END DO 
            END DO

			DO k=decomp%xst(3),decomp%xen(3)
                DO j=decomp%xst(2),decomp%xen(2)
                    DO i=decomp%xst(1),decomp%xen(1)
				        w(i,j,k)=w(i,j,k)*scalemodes
			        END DO 
                END DO 
            END DO
						
			mychg(1) =maxval(abs(utemp-u))
			mychg(2) =maxval(abs(vtemp-v))
			mychg(3) =maxval(abs(wtemp-w))
			CALL MPI_ALLREDUCE(mychg,allchg,3,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
			chg=allchg(1)+allchg(2)+allchg(3)
			!PRINT *,'chg:',chg
		END DO
		time(n+1)=n*dt
		IF (myid == 0) PRINT *,'time',n*dt
    END DO
    
    ! DO k=decomp%xst(3),decomp%xen(3)
    !     DO j=decomp%xst(2),decomp%xen(2)
    !         DO i=decomp%xst(1),decomp%xen(1)
    !             utemp(i,j,k)=u(i,j,k)-(-0.5*( factor*cos(x(i))*sin(y(j))*sin(z(k))&
    !                             +sin(x(i))*cos(y(j))*cos(z(k)) )*exp(-(factor**2)*time(Nt+1)/Re))
    !         END DO
    !     END DO
    ! END DO

    ! DO k=decomp%xst(3),decomp%xen(3)
    !     DO j=decomp%xst(2),decomp%xen(2)
    !         DO i=decomp%xst(1),decomp%xen(1)
    !             vtemp(i,j,k)=v(i,j,k)-(0.5*(  factor*sin(x(i))*cos(y(j))*sin(z(k))&
    !                             -cos(x(i))*sin(y(j))*cos(z(k)) )*exp(-(factor**2)*time(Nt+1)/Re))
    !         END DO 
    !     END DO 
    ! END DO

    ! DO k=decomp%xst(3),decomp%xen(3)
    !     DO j=decomp%xst(2),decomp%xen(2)
    !         DO i=decomp%xst(1),decomp%xen(1)
    !             wtemp(i,j,k)=w(i,j,k)-(cos(x(i))*cos(y(j))*sin(z(k))*exp(-(factor**2)*time(Nt+1)/Re))
    !         END DO 
    !     END DO 
    ! END DO
    
    chg=maxval(abs(utemp))+maxval(abs(vtemp))+maxval(abs(wtemp))
    IF (myid == 0) PRINT*,'The error at the final timestep is',chg

    CALL cpu_time(e)
    IF (myid == 0) PRINT *, "Time took to execute: ", (e-s), "s"
    
    CALL decomp_2d_fft_finalize
  	CALL decomp_2d_finalize
    DEALLOCATE(x,y,z,time,mychg,allchg,u,v,w,ux,uy,uz,vx,vy,vz,wx,wy,wz,uold,uxold,uyold,uzold,&
                    vold,vxold,vyold,vzold,wold,wxold,wyold,wzold,utemp,vtemp,wtemp,&
                    temp_r,kx,ky,kz,uhat,vhat,what,rhsuhatfix,rhsvhatfix,&
                    rhswhatfix,phat,nonlinuhat,nonlinvhat,nonlinwhat,temp_c,&
                    realtemp,stat=AllocateStatus)		
    IF (AllocateStatus .ne. 0) STOP
    IF (myid == 0) PRINT *,'Program execution complete'
    CALL MPI_FINALIZE(ierr)	

END PROGRAM main