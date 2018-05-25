#include "parallel.h"


void Program_Message(char *txt)
/* produces a stderr text output  */
{
    int myrank;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
    fflush(stdout);
    fflush(stderr);
}


void Programm_Sync(char *txt)
/* produces a stderr textoutput and synchronize all processes */
{
    int myrank;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Barrier(MPI_COMM_WORLD);                             /* synchronize output */  
    fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
    fflush(stdout);
    fflush(stderr);
    MPI_Barrier(MPI_COMM_WORLD);
}



void Programm_Stop(char *txt)
/* all processes will produce a text output, be synchronized and finished */
{
    int myrank;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Barrier(MPI_COMM_WORLD);                           /* synchronize output */
    fprintf(stderr,"-STOP- P:%2d : %s\n",myrank,txt);
    fflush(stdout);
    fflush(stderr);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    exit(1);
}


/* thread-wise assignment of MPI variables */
void init_parallel (int iproc, int jproc, 
                    int imax, int jmax, int *myrank,
                    int *il, int *ir, int *jb, int *jt, 
                    int *rank_l, int *rank_r, 
                    int *rank_b, int *rank_t,
                    int *omg_i, int *omg_j, int num_proc)
{
	int i_per_iproc, j_per_jproc, i_rem, j_rem;

	/* Set sub-domain coordinates */
	*omg_i = ((*myrank) % iproc) + 1;
	*omg_j = (((*myrank) - (*omg_i) + 1) / iproc) + 1;

	/* Compute il, ir for each sub-domain*/
	i_per_iproc = (imax / iproc);
	j_per_jproc = (jmax / jproc);
	i_rem = (imax % iproc);
	j_rem = (jmax % jproc);

	/* for il ir*/
	if ((*omg_i) == 1)   /* to rank zero assign the remainder*/
	{
	  *il = (((*omg_i) -1) * i_per_iproc) + 1;
	}
	else
	{
	  *il = (((*omg_i) -1) * i_per_iproc) + i_rem + 1;
	}
	*ir = ((*omg_i) * i_per_iproc) + i_rem;

	/* for jb jt*/
	if ((*omg_j) == 1)
	{
	  *jb = (((*omg_j) -1) * j_per_jproc) + 1;
	}
	else
	{
	  *jb = (((*omg_j) -1) * j_per_jproc) + j_rem + 1;
	}
	*jt = ((*omg_j) * j_per_jproc) + j_rem;

    /* Assign rank of neighbour to each sub-domain: rank_l, rank_r, rank_b, rank_t*/
	/* Left boundary*/
	if ((*il) == 1)
	{
	  *rank_l = MPI_PROC_NULL;
	}
	else
	{
	  *rank_l = (*myrank) - 1;
	}

	/* Right boundary*/
	if ((*ir) == imax)
	{
	  *rank_r = MPI_PROC_NULL;
	}
	else
	{
	  *rank_r = (*myrank) + 1;
	}

	/* Bottom boundary*/
	if ((*jb) == 1)
	{
	  *rank_b = MPI_PROC_NULL;
	}
	else
	{
	  *rank_b = (*myrank) - iproc;
	}

	/* Top boundary*/
	if ((*jt) == jmax)
	{
	  *rank_t = MPI_PROC_NULL;
	}
	else
	{
	  *rank_t = (*myrank) + iproc;
	}

}



/* communication of pressure values along sub-domain boundaries */
void pressure_comm( double **P,
                    int il, int ir,
                    int jb,int jt,
                    int rank_l,int rank_r,
                    int rank_b,int rank_t,
                    double *bufSend,double *bufRecv, 
                    MPI_Status *status, int chunk){
    int xdim;
    int ydim;
    int count;

    //Sending and receiving pressure values from the left and the right sub-domains
     ydim = jt - jb + 1;
    
    /* Step 1: Send to the left - receive from the right */
    bufSend = (double*)malloc(ydim*sizeof(double));
    bufRecv = (double*)malloc(ydim*sizeof(double));

    if ( MPI_PROC_NULL != rank_l ){
        count = 0;
        for(int i = jb; i<=jt; i++, count++){
          bufSend[count] = P[i][il];
        }

        MPI_Send( bufSend, ydim, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD );
        MPI_Recv( bufRecv, ydim, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD, &status );

        count = 0;
        for(int i = jb; i<=jt; i++, count++){
          P[i][ir] = bufRecv[count];
        }
    }
    

    /* Step 2: Send to the right - receive from the left */
    if ( MPI_PROC_NULL != rank_r ){
        count = 0;
        for(int i = jb; i<=jt; i++, count++){
          bufSend[count] = P[i][ir];
        }

        MPI_Send( bufSend, ydim , MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD );
        MPI_Recv( bufRecv, ydim, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD, &status );

        count = 0;
        for(int i = jb; i<=jt; i++, count++){
          P[i][il] = bufRecv[count];
        }
    }


    //Sending and receiving pressure values from the top and the bottom sub-domains
    xdim = ir - il + 1;

    /* Step 3: Send to the top - receive from the bottom */
    bufSend = (double*)malloc(xdim*sizeof(double));
    bufRecv = (double*)malloc(xdim*sizeof(double));

    if ( MPI_PROC_NULL != rank_t ){
        count = 0;
        for(int j = il; j<=ir; j++, count++){
          bufSend[count] = P[jt][j];
        }

        MPI_Send( bufSend, xdim, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD );
        MPI_Recv( bufRecv, xdim, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD, &status );

        count = 0;
        for(int j = il; i<=ir; j++, count++){
          P[jb][j] = bufRecv[count];
        }
    }


    /* Step 4: Send to the bottom - receive from the top */
    if ( MPI_PROC_NULL != rank_b ){
        count = 0;
        for(int j = il; j<=ir; j++, count++){
          bufSend[count] = P[jb][j];
        }

        MPI_Send( bufSend, xdim, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD );
        MPI_Recv( bufRecv, xdim, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD, &status );

        count = 0;
        for(int j = il; i<=ir; j++, count++){
          P[jt][j] = bufRecv[count];
        }
    }
}



/* communication of u and v velocity values along sub-domain boundaries */
void uv_comm   (double **U,double **V,
                int il,int ir, 
                int jb,int jt,
                int rank_l,int rank_r,
                int rank_b,int rank_t,
                double *bufSend, double *bufRecv, 
                MPI_Status *status, int chunk){

    int xdim;
    int ydim;
    int count;

    //Sending and receiving u and v velocity values from the left and the right sub-domains
     ydim = 2*(jt - jb + 1);
    
    /* Step 1: Send to the left - receive from the right */
    bufSend = (double*)malloc(ydim*sizeof(double));
    bufRecv = (double*)malloc(ydim*sizeof(double));

    if ( MPI_PROC_NULL != rank_l ){
        count = 0;
        for(int i = jb; i<=jt; i++, count++){
          bufSend[count] = U[i][il];
        }
        for(int i = jb; i<=jt; i++, count++){
          bufSend[count] = V[i][il];
        }

        MPI_Send( bufSend, ydim, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD );
        MPI_Recv( bufRecv, ydim, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD, &status );

        count = 0;
        for(int i = jb; i<=jt; i++, count++){
          U[i][ir] = bufRecv[count];
        }
        for(int i = jb; i<=jt; i++, count++){
          V[i][ir] = bufRecv[count];
        }
    }
    

    /* Step 2: Send to the right - receive from the left */
    if ( MPI_PROC_NULL != rank_r ){
        count = 0;
        for(int i = jb; i<=jt; i++, count++){
          bufSend[count] = U[i][ir];
        }
        for(int i = jb; i<=jt; i++, count++){
          bufSend[count] = V[i][ir];
        }

        MPI_Send( bufSend, ydim , MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD );
        MPI_Recv( bufRecv, ydim, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD, &status );

        count = 0;
        for(int i = jb; i<=jt; i++, count++){
          U[i][il] = bufRecv[count];
        }
        for(int i = jb; i<=jt; i++, count++){
          V[i][il] = bufRecv[count];
        }
    }


    //Sending and receiving u and v velocity values from the top and the bottom sub-domains
    xdim = 2*(ir - il + 1);

    /* Step 3: Send to the top - receive from the bottom */
    bufSend = (double*)malloc(xdim*sizeof(double));
    bufRecv = (double*)malloc(xdim*sizeof(double));

    if ( MPI_PROC_NULL != rank_t ){
        count = 0;
        for(int j = il; j<=ir; j++, count++){
          bufSend[count] = U[jt][j];
        }
        for(int j = il; j<=ir; j++, count++){
          bufSend[count] = V[jt][j];
        }

        MPI_Send( bufSend, xdim, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD );
        MPI_Recv( bufRecv, xdim, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD, &status );

        count = 0;
        for(int j = il; i<=ir; j++, count++){
          U[jb][j] = bufRecv[count];
        }
        for(int j = il; i<=ir; j++, count++){
          V[jb][j] = bufRecv[count];
        }
    }


    /* Step 4: Send to the bottom - receive from the top */
    if ( MPI_PROC_NULL != rank_b ){
        count = 0;
        for(int j = il; j<=ir; j++, count++){
          bufSend[count] = U[jb][j];
        }
        for(int j = il; j<=ir; j++, count++){
          bufSend[count] = V[jb][j];
        }

        MPI_Send( bufSend, xdim, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD );
        MPI_Recv( bufRecv, xdim, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD, &status );

        count = 0;
        for(int j = il; i<=ir; j++, count++){
          U[jt][j] = bufRecv[count];
        }
        for(int j = il; i<=ir; j++, count++){
          V[jt][j] = bufRecv[count];
        }
    }

}


