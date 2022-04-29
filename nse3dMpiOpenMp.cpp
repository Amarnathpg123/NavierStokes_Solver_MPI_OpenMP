#include <fftw3-mpi.h>
#include <complex>
#include <vector>
#include <math.h>
#include <omp.h>

using namespace std;


//mpic++ -std=c++11 -O3 nse3dMpiOpenMp.cpp -o nse3dMpiOpenMp -lfftw3_mpi -lfftw3 -fopenmp -lm
//mpirun -np 4 ./nse3dMpiOpenMp

#define vd vector<double>
#define vcd vector<complex<double>>


// pseudospectral Fourier-Galerkin method
// uing 1d slab decompositon

int main() {

    int rank, N, Nf;
    double L, dx, nu, dt, T;

    double t1, start_time;

    MPI::Init();
    fftw_mpi_init();

    rank = MPI::COMM_WORLD.Get_rank();
    double pi = atan(1.0)*4.0;;
    long alloc_local, local_n0, local_0_start, local_n1, local_1_start, i, j, k;

    vd s_in(1), s_out(1);  

    omp_set_num_threads(omp_get_max_threads());
    
    N = pow(2,6);    // --> change N here
    nu = 0.000625;
    T = 0.1;
    dt = 0.001;

    

    L = 2*pi;
    dx = L / N;

    if (rank==0)
        std::cout << "N = " << N << std::endl;

    vd a(4), b(3);

    //coefficients for RK4 method
    a[0] = 1./6.; a[1] = 1./3.; a[2] = 1./3.; a[3] = 1./6.;
    b[0] = 0.5; b[1] = 0.5; b[2] = 1.0;

    int tot = N*N*N;
    Nf = N/2+1;
    alloc_local = fftw_mpi_local_size_3d_transposed(N, N, Nf, MPI::COMM_WORLD,
                                            &local_n0, &local_0_start,
                                            &local_n1, &local_1_start);

    vd U(2*alloc_local), V(2*alloc_local), W(2*alloc_local), U_tmp(2*alloc_local), V_tmp(2*alloc_local), 
        W_tmp(2*alloc_local), CU(2*alloc_local), CV(2*alloc_local), CW(2*alloc_local), 
        Fx(2*alloc_local),Fy(2*alloc_local),Fz(2*alloc_local);

    vector<int> dealias(2*alloc_local);   // for dealiasing
    vd kk(2*alloc_local);

    vcd U_hat(alloc_local), V_hat(alloc_local), W_hat(alloc_local), P_hat(alloc_local), U_hat0(alloc_local), 
        V_hat0(alloc_local), W_hat0(alloc_local), U_hat1(alloc_local), V_hat1(alloc_local), W_hat1(alloc_local), 
        dU(alloc_local), dV(alloc_local), dW(alloc_local), curlX(alloc_local), curlY(alloc_local), curlZ(alloc_local);

    // Starting time
    MPI::COMM_WORLD.Barrier();
    start_time = MPI::Wtime();

    vd kx(N), kz(Nf);

    //fourier grid

    for (i=0; i<N/2; i++)
    {
        kx[i] = i;
        kz[i] = i;
    }

    kz[N/2] = N/2;
    for (i=-N/2; i<0; i++)
        kx[i+N] = i;

    //fftw_plan plan_backward;
    fftw_plan rfftn, irfftn;
    rfftn = fftw_mpi_plan_dft_r2c_3d(N, N, N, U.data(), reinterpret_cast<fftw_complex*>(U_hat.data()),
                                    MPI::COMM_WORLD, FFTW_MPI_TRANSPOSED_OUT);
    irfftn = fftw_mpi_plan_dft_c2r_3d(N, N, N, reinterpret_cast<fftw_complex*>(U_hat.data()),  U.data(),
                                    MPI::COMM_WORLD, FFTW_MPI_TRANSPOSED_IN);

    //Constant Force values
    #pragma omp for
    for (i=0; i<local_n0; i++)
        for (j=0; j<N; j++)
            for (k=0; k<N; k++) {
                const int z = (i*N+j)*2*Nf+k;
                Fx[z] = 0.01;   
                Fy[z] = 0.001;
                Fz[z] = 0.0001;
            }

    // boundary conditions
    #pragma omp for
    for (i=0; i<local_n0; i++)
        for (j=0; j<N; j++)
            for (k=0; k<N; k++) {
                const int z = (i*N+j)*2*Nf+k;
                U[z] = sin(dx*(i+local_0_start))*cos(dx*j)*cos(dx*k);   //---> initial boundary conditions
                V[z] = -cos(dx*(i+local_0_start))*sin(dx*j)*cos(dx*k);
                W[z] = 0.0;
            }

    //taking velocity field from physical to fourier
    fftw_mpi_execute_dft_r2c( rfftn, U.data(), reinterpret_cast<fftw_complex*>(U_hat.data()));
    fftw_mpi_execute_dft_r2c( rfftn, V.data(), reinterpret_cast<fftw_complex*>(V_hat.data()));
    fftw_mpi_execute_dft_r2c( rfftn, W.data(), reinterpret_cast<fftw_complex*>(W_hat.data()));


    double kmax = (2.0/3)*Nf;

    //dealiasing using 2/3 rule
    #pragma omp for
    for (i=0; i<local_n1; i++)
        for (j=0; j<N; j++)
            for (k=0; k<Nf; k++)
            {
                const int z = (i*N+j)*Nf+k;
                dealias[z] = (abs(kx[i+local_1_start])<kmax)*(abs(kx[j])<kmax)*(abs(kx[k])<kmax);
            }

    #pragma omp for
    for (i=0; i<local_n1; i++)
        for (j=0; j<N; j++)
            for (k=0; k<Nf; k++)
            {
                const int z = (i*N+j)*Nf+k;
                int m = kx[i+local_1_start]*kx[i+local_1_start] + kx[j]*kx[j] + kx[k]*kx[k];
                kk[z] = m > 0 ? m : 1;
            }

    complex<double> one(0,1);

    double t=0.0;
    int tstep = 0;

    //time stepping
    while (t < T-1e-8) {
        t += dt;
        tstep++;

        //copying fourier space velocity to mpi work arrays
        std::copy(U_hat.begin(), U_hat.end(), U_hat0.begin());
        std::copy(V_hat.begin(), V_hat.end(), V_hat0.begin());
        std::copy(W_hat.begin(), W_hat.end(), W_hat0.begin());
        std::copy(U_hat.begin(), U_hat.end(), U_hat1.begin());
        std::copy(U_hat.begin(), U_hat.end(), U_hat1.begin());
        std::copy(U_hat.begin(), U_hat.end(), U_hat1.begin());

        //RK-4 Method
        for (int rk=0; rk<4; rk++) {
            if (rk > 0) {

                fftw_mpi_execute_dft_c2r(irfftn, reinterpret_cast<fftw_complex*>(U_hat.data()), U.data());
                fftw_mpi_execute_dft_c2r(irfftn, reinterpret_cast<fftw_complex*>(V_hat.data()), V.data());
                fftw_mpi_execute_dft_c2r(irfftn, reinterpret_cast<fftw_complex*>(W_hat.data()), W.data());
                
                for (k=0; k<U.size(); k++) {
                    U[k] /= tot;
                    V[k] /= tot;
                    W[k] /= tot;
                }
            }

            // Compute curl
            #pragma omp for
            for (i=0; i<local_n1; i++)
                for (j=0; j<N; j++)
                    for (k=0; k<Nf; k++) {
                        const int z = (i*N+j)*Nf+k;
                        curlZ[z] = one*(kx[i+local_1_start]*V_hat[z]-kx[j]*U_hat[z]);
                        curlY[z] = one*(kz[k]*U_hat[z]-kx[i+local_1_start]*W_hat[z]);
                        curlX[z] = one*(kx[j]*W_hat[z]-kz[k]*V_hat[z]);
                    }

            // taking curl into fourier space
            fftw_mpi_execute_dft_c2r(irfftn, reinterpret_cast<fftw_complex*>(curlX.data()), CU.data());
            fftw_mpi_execute_dft_c2r(irfftn, reinterpret_cast<fftw_complex*>(curlY.data()), CV.data());
            fftw_mpi_execute_dft_c2r(irfftn, reinterpret_cast<fftw_complex*>(curlZ.data()), CW.data());

            for (k=0; k<CU.size(); k++)
            {
                CU[k] /= tot;
                CV[k] /= tot;
                CW[k] /= tot;
            }

            // Cross
            #pragma omp for
            for (i=0; i<local_n0; i++)
                for (j=0; j<N; j++)
                    for (k=0; k<N; k++) {
                        const int z = (i*N+j)*2*Nf+k;
                        U_tmp[z] = V[z]*CW[z]-W[z]*CV[z]+Fx[z];
                        V_tmp[z] = W[z]*CU[z]-U[z]*CW[z]+Fy[z];
                        W_tmp[z] = U[z]*CV[z]-V[z]*CU[z]+Fz[z];
                    }
            
            //taking cross into fourier space
            fftw_mpi_execute_dft_r2c( rfftn, U_tmp.data(), reinterpret_cast<fftw_complex*>(dU.data()));
            fftw_mpi_execute_dft_r2c( rfftn, V_tmp.data(), reinterpret_cast<fftw_complex*>(dV.data()));
            fftw_mpi_execute_dft_r2c( rfftn, W_tmp.data(), reinterpret_cast<fftw_complex*>(dW.data()));

            //dealiasing
            #pragma omp for
            for (i=0; i<local_n1; i++)
                for (j=0; j<N; j++)
                    for (k=0; k<Nf; k++) {
                        const int z = (i*N+j)*Nf+k;
                        dU[z] *= (dealias[z]*dt);
                        dV[z] *= (dealias[z]*dt);
                        dW[z] *= (dealias[z]*dt);
                    }
            
            //taking care of pressure and other linear terms
            #pragma omp for
            for (i=0; i<local_n1; i++)
                for (j=0; j<N; j++)
                    for (k=0; k<Nf; k++) {
                        const int z = (i*N+j)*Nf+k;
                        P_hat[z] = (dU[z]*kx[i+local_1_start] + dV[z]*kx[j] + dW[z]*kz[k])/kk[z];
                        dU[z] -= (P_hat[z]*kx[i+local_1_start] + nu*dt*kk[z]*U_hat[z]);
                        dV[z] -= (P_hat[z]*kx[j] + nu*dt*kk[z]*V_hat[z]);
                        dW[z] -= (P_hat[z]*kz[k] + nu*dt*kk[z]*W_hat[z]);
                    }
            
            if (rk < 3) {
                #pragma omp for
                for (i=0; i<local_n1; i++)
                    for (j=0; j<N; j++)
                        for (k=0; k<Nf; k++) {
                            const int z = (i*N+j)*Nf+k;
                            U_hat[z] = U_hat0[z] + b[rk]*dU[z];
                            V_hat[z] = V_hat0[z] + b[rk]*dV[z];
                            W_hat[z] = W_hat0[z] + b[rk]*dW[z];
                        }

            }
            
            #pragma omp for
            for (i=0; i<local_n1; i++)
                for (j=0; j<N; j++)
                    for (k=0; k<Nf; k++) {
                        const int z = (i*N+j)*Nf+k;
                        U_hat1[z] += a[rk]*dU[z];
                        V_hat1[z] += a[rk]*dV[z];
                        W_hat1[z] += a[rk]*dW[z];
                    }
        }

        //copying the temp fourier velocity into the hat arrays
        #pragma omp for
        for (i=0; i<local_n1; i++)
            for (j=0; j<N; j++)
                for (k=0; k<Nf; k++) {
                    const int z = (i*N+j)*Nf+k;
                    U_hat[z] = U_hat1[z];
                    V_hat[z] = V_hat1[z];
                    W_hat[z] = W_hat1[z];
                }

        //printing the total energy of the system at every even step
        if (tstep % 2 == 0) {
            s_in[0] = 0.0;
            #pragma omp for
            for (i=0; i<local_n0; i++)
                for (j=0; j<N; j++)
                    for (k=0; k<N; k++) {
                        int z = (i*N+j)*2*Nf+k;
                        s_in[0] += (U[z]*U[z] + V[z]*V[z] + W[z]*W[z]);
                    }

            s_in[0] *= (0.5*dx*dx*dx/L/L/L);

            MPI::COMM_WORLD.Reduce(s_in.data(), s_out.data(), 1, MPI::DOUBLE, MPI::SUM, 0);

            // if (rank==0)
            //     std::cout << " k = " << s_out[0] << std::endl;
        }
    }

    s_in[0] = 0.0;


    //get final velocities from U,V and W
    //calculating total final energy
    #pragma omp for
    for (i=0; i<local_n0; i++)
        for (j=0; j<N; j++)
            for (k=0; k<N; k++) {
                int z = (i*N+j)*2*Nf+k;
                s_in[0] += (U[z]*U[z] + V[z]*V[z] + W[z]*W[z]);
            }

    s_in[0] *= (0.5*dx*dx*dx/L/L/L);

    MPI::COMM_WORLD.Reduce(s_in.data(), s_out.data(), 1, MPI::DOUBLE, MPI::SUM, 0);

    if (rank==0) 
        std::cout << "\nFinal Energy: \n k = " << s_out[0] << std::endl;

    MPI::COMM_WORLD.Barrier();

    t1 = MPI::Wtime();
    
    if (rank == 0)
        std::cout << "Total time taken = " << t1 - start_time  << "s" << std::endl;

    fftw_destroy_plan(rfftn);
    fftw_destroy_plan(irfftn);

    MPI::Finalize ( );
}
