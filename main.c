#include "helper.h"
#include "visual.h"
#include "init.h"
#include "uvp.h"
#include "boundary_val.h"
#include "sor.h"
#include "parallel.h"
#include "string.h"
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>


/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the 
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */	
int main(int argn, char** args){
			double xlength;           /* length of the domain x-dir.*/
                    	double ylength;           /* length of the domain y-dir.*/
                    	int  imax;                /* number of cells x-direction*/
                    	int  jmax;                /* number of cells y-direction*/
                    	double dx;                /* length of a cell x-dir. */
                    	double dy;                /* length of a cell y-dir. */
                   	double t ;                /* gravitation y-direction */
                    	double t_end;             /* end time */
                    	double dt;                /* time step */
                    	double tau;               /* safety factor for time step*/
                    	double dt_value;          /* time interval for writing visualization data in afile*/
                    	int n;                    /* current time iteration step*/
                    	int itermax;              /* max. number of iterations  */
                    	int it;                   /* SOR iteration counter*/
                    	double res;               /*residual norm of the pressure equation*/
                    	double eps;               /* accuracy bound for pressure*/
                    	double alpha;             /* uppwind differencing factor*/
                    	double omg;               /* relaxation factor */
                    	double Re;                /* reynolds number   */
                    	double UI;                /* velocity x-direction */
                    	double VI;                /* velocity y-direction */
                   	double PI;                /* pressure */
                    	double GX;                /* gravitation x-direction */
                    	double GY;
		  	const char *szFileName="cavity100.dat";       /* name of the file */
			char *problemOutput = malloc(strlen(szFileName) + 10);
			char buffer[4] = {0};
			/* MPI variables */
			int iproc = 0;
			int jproc = 0;
			int myrank = 0;
			int il = 0;
			int ir = 0;
			int jb = 0;
			int jt = 0;
			int rank_l = 0;
			int rank_r = 0;
			int rank_b = 0;
			int rank_t = 0;
			int omg_i = 0;
			int omg_j = 0;
			int num_proc = 0;
			int comm_size = 0;
			int x_dim, y_dim;

			/* Initialize MPI and begin parallel processes */
			MPI_Init(&argn, &args);
			MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
			MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
			
			/* Append name for output vtk file */
			strcpy(problemOutput, szFileName);
			strcat(problemOutput, "_output.");
			sprintf(buffer, "%i", myrank);
			strcat(problemOutput, buffer);

		        //Reading the program configuration file using read_parameters()
			if(myrank==0){
		 		int param = read_parameters(szFileName, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value, &iproc, &jproc);
				//printf("dx = :%d and dy = %d", dx,dy);
				printf("Parameter Read: (main.c)  %d \n",myrank);
				param++; 		// Just using param so that C does not throw an error message
				
				param=comm_size;
				/*if(param != 1){
					Programm_Stop("Input parameters potentially corrupt!");
					return -1;
				}
				MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
				if(comm_size != (iproc*jproc)){
					printf("The number of processes entered via 'mpirun -np %u' does not equal the size specified by the input file (iproc * jproc = %u)\n", comm_size, iproc*jproc);
					return -1;
				}*/
			}
			
			printf("Broadcasting all variables :rank = %d : (main.c)\n ",myrank);

			/* Broadcast parameters to the remaining processes */
			MPI_Bcast(&Re, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&UI, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&VI, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&PI, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&GX, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&GY, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&t_end, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&xlength, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&ylength, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&dx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&dy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&imax, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&jmax, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&omg, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&itermax, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&eps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&dt_value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&iproc, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&jproc, 1, MPI_INT, 0, MPI_COMM_WORLD);

			printf("Init_parallel : rank = %d : (main.c)\n ",myrank);

			/* Initialize process dependent variables */
			init_parallel(iproc, jproc, imax, jmax, &myrank, &il, &ir, &jb, &jt, &rank_l,
			    &rank_r, &rank_b, &rank_t, &omg_i, &omg_j, num_proc);
			printf("Rank : %d, Omg : %dx%d, ilxir : %dx%d, jbxjt : %dx%d ",myrank, omg_i, omg_j, il,ir,jb,jt);
			x_dim = ir - il + 1;
			y_dim = jt - jb + 1;
			printf("Creating matrices : rank = %d : (main.c)\n ",myrank);
			double **P = matrix(0, x_dim+1, 0, y_dim+1);
			double **U = matrix(0, x_dim+1, 0, y_dim+1); //Dynamically allocating memmory for matrices U,V,P, RS, F and G
			double **V = matrix(0, x_dim+1, 0, y_dim+1); 
		 	double **F = matrix(0, x_dim+1, 0, y_dim+1); 
			double **G = matrix(0, x_dim+1, 0, y_dim+1);    	
		 	double **RS = matrix(0, x_dim+1, 0, y_dim+1);
			printf("Init_uvp by rank = %d in domain:%dx%d (main.c)\n ",myrank, omg_i, omg_j);
			/* Initialize matrices for velocity, pressure, rhs, etc. */
			init_uvp(UI, VI, PI, x_dim, y_dim, U, V, P);
			printf("Init_uvp completed by rank = %d in domain:%dx%d (main.c)\n \n",myrank, omg_i, omg_j);
		 	
		 	t=0;
			n=0;
			int n1 = 0;
			int chunk = 0;
			//int residue;
			double* bufSend = 0;
			double* bufRecv = 0;
			MPI_Status status;
			//init_uvp(UI,VI,PI,imax,jmax,U,V,P);    //Initializing U, V and P
			printf("Starting while loop for simulation : rank = %d : (main.c)\n ",myrank);
			while (t<t_end)
			{
				printf("boundaryvalues : rank = %d : (main.c)\n ",myrank);
				boundaryvalues(x_dim, y_dim, U, V, rank_l, rank_r, rank_b, rank_t);		 //Setting the boundary values for the next time step.
				printf("calculating Fg : rank = %d : (main.c)\n ",myrank);
				calculate_fg(Re,GX,GY,alpha,dt,dx,dy,x_dim,y_dim,U,V,F,G);			 //Determining the values of F and G (diffusion and confection).
				printf("calculating RS : rank = %d : (main.c)\n ",myrank);
				calculate_rs(dt,dx,dy,x_dim,y_dim,F,G,RS);					 //Calculating the right hand side of the pressure equation.
				
				it=0;
				res = 1;
				printf("Entering while loop of SOR : rank = %d in domain:%dx%d (main.c)\n \n ",myrank, omg_i, omg_j);
				while (it<itermax && res>eps)  //Iterating the pressure poisson equation until the residual becomes smaller than eps or the maximal number of iterations is performed. 
				{
					sor(omg, dx, dy, P, RS, &res, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, imax, jmax);
					it++;
				}
				if (myrank == 0)
				{
					printf("res=%f, it=%u, t=%f, dt=%f\n", res, it, t, dt);
				}
					/*printf("Entering SOR.c : rank = %d in domain:%dx%d (main.c)\n \n",myrank, omg_i, omg_j);
					sor(omg,dx,dy,imax,jmax,P,RS,&res,il,ir,jb,jt,rank_l,rank_r,rank_b,rank_t);  */    				//Within the iteration loop the operation sor() is used.

 					/* Communicate between processes regarding pressure boundaries */
					/*printf("Entering pressure_comm.c : rank = %d in domain:%dx%d (main.c)\n \n",myrank, omg_i, omg_j);
					pressure_comm(P, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, bufSend, bufRecv, &status, chunk);*/
			 		 /* Sum the squares of all local residuals then square root that sum for global residual */
				
					/*printf("Exiting pressure_comm.c : rank = %d in domain:%dx%d (main.c)\n \n",myrank, omg_i, omg_j);
					Programm_Sync("Synchronising all the threads for calculating residue");
					printf("Begin ALLREDUCE : rank = %d in domain:%dx%d (main.c)\n ",myrank, omg_i, omg_j);
					MPI_Allreduce(&res, &residue, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
					residue = sqrt((residue)/(imax*jmax));
					it=it+1;
					printf("Ending ALLREDUCE : rank = %d in domain:%dx%d (main.c)\n \n",myrank, omg_i, omg_j);*/
				
				printf("Exiting while loop of SOR : rank = %d in domain: %dx%d (main.c)\n \n",myrank, omg_i, omg_j);
				printf("Calculate UV : rank = %d in domain:%dx%d (main.c)\n ",myrank, omg_i, omg_j);
				calculate_uv(dt,dx,dy,x_dim,y_dim,U,V,F,G,P);   //Calculating the velocity at the next time step.
				
				printf("UV_comm : rank = %d in domain: %dx%d (main.c)\n \n",myrank, omg_i, omg_j);
				uv_comm(U, V, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, bufSend, bufRecv, &status, chunk);
				printf("Calculating dt : rank = %d in domain:%dx%d (main.c)\n \n",myrank, omg_i, omg_j);
				Programm_Sync("Synchronising all the threads for calculating dt");
				calculate_dt(Re,tau,&dt,dx,dy,x_dim,y_dim,U,V); 

				if(t>=n1*dt_value){
					output_uvp(U, V, P, il, ir, jb, jt, omg_i, omg_j, problemOutput, n1, dx, dy);
					n1 = n1 + 1;
				}
				
				t=t+dt;
				n=n+1;
				
			}
			printf("Exiting while loop for simulation : rank = %d in domain:%dx%d (main.c)\n \n",myrank, omg_i, omg_j);

			/* Deallocate heap memory */
			free_matrix(U, 0, x_dim+1, 0, y_dim+1);
			free_matrix(V, 0, x_dim+1, 0, y_dim+1);
			free_matrix(P, 0, x_dim+1, 0, y_dim+1);
			free_matrix(F, 0, x_dim+1, 0, y_dim+1);
			free_matrix(G, 0, x_dim+1, 0, y_dim+1);
			free_matrix(RS, 0, x_dim+1, 0, y_dim+1);
			
			/* Finalize MPI */
			Programm_Stop("End of simulation.");
			
  return -1;
}
