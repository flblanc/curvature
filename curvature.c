/******************************************************************************
 *
 *  CURVATURE Module for WORDOM
 *  Calculates the mean and gaussian quadratic curvatures of a selection of atoms
 *
 *  Copyright (C) 2015  Florian Blanc
 *
 *  CURVATURE Module for WORDOM is free software: you can redistribute it and/or 
 *  modify it under the terms of the GNU General Public License as published 
 *  by the Free Software Foundation, either version 3 of the License, or 
 *  (at your option) any later version.
 *
 *  CURVATURE Module for WORDOM is distributed in the hope that it will be
 *  useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/


#include "wordom.h"
#include "tools.h"
#include "analysis.h"
#include "curvature.h"

/* Overview of the algorithm:
 * --------------------------
 * 
 * Input: a structure, a selection of atoms, a trajectory
 * Ouput: the time series of the curvature of the selection of atoms over the trajectory
 * 
 * Procedure:
 * ----------
 * Initialization: Read in the structure and the selection
 * For each frame:
 *      get the coordinates of the atoms of the selection
 *      fit a quadratic surface over the set of atoms using the Levenberg-Marquardt algorithm
 *      compute the curvature 
 *      print the frame number and the curvatures (mean and gaussian)
 */


int 
Read_iCurvature ( char **input, int inp_index, struct inp_curvature *inp_curvature , char *printout, Molecule *molecule )
{
    /* I have to define:
     * a inp_curvature structure
     * the data structure for the levmar optimizer
     */
      
    char          buffer[256];
    char          title[64]="";
    int           gotit;
    
    memset ( buffer, '\0', sizeof(buffer));
    
    // Read the input script
    
    while ( strncmp(buffer,"END",3) ) {
        
        gotit = 0;
        memset ( buffer, '\0', sizeof(buffer));
        sprintf( buffer, "%s", input[inp_index]);
        
        if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#') {
            gotit = 1;
        } 
        
        else if ( !strncmp(buffer, "--TITLE", 7)) {
            memset( title, '\0', 64);
            sscanf( buffer, "--TITLE %s", title);
            gotit = 1;
        }
        
        else if ( !strncmp(buffer, "--SELE", 6) ) {
            sscanf( buffer, "--SELE %[^:]", inp_curvature->sele.selestring);
            GetSele ( inp_curvature->sele.selestring, &inp_curvature->sele, molecule);
            if( inp_curvature->sele.nselatm > 1 ) {
                fprintf( stderr, "Warning: GC doesn't work with PBC yet; centers will be computed w/o PBC, distances with PBC\n");
            }
            if( inp_curvature->sele.nselatm < 7 ) {
                fprintf( stderr, "\nERROR! Selection %s must have more than 7 atoms (have %d)\n\n", inp_curvature->sele.selestring, inp_curvature->sele.nselatm);
                exit(0);
            }
            gotit = 1;
        }
        
        if( gotit==0 ) {
            fprintf( stderr, "Could not understand option: %s\n", buffer);
            exit(5);
        }
        inp_index++;
    }
    
    
    if (title[0]=='\0') { //if no title is provided
        sprintf( printout, "%10s %10s", "mean", "gaussian");
    }
    else {
        sprintf( printout, "%4.2s-mean %4.2-gaussian", title, title);
    }
    
    return 24;
}


/* Auxiliary functions:
 * quadsurf is the analytic expression of the quadratic surface
 * jacquadsurf is the analytic jacobian matrix of the quadratic surface;
 * using the analytic jacobian avoids the use of finite difference approximations
 * and should improve the convergence properties.
 */
 
// Analytic expression of the quadratic surface; it has 6 parameters that are
// to be fitted over the atomic coordinates of the selection.


void quadsurf(double *p,  double *z, int m, int n, void *adata)
{
    register int i; //I use a register integer because I'm smarter than the compiler. And the compiler will never know because it's so dumb,
                    // it cannot read comments. You're so lame, converter.
    double *x, *y;
    struct MyData *dptr;
    dptr=(struct mydata *)adata; // adata is passed as void *, dptr = data pointer

    x=dptr->XD;
    y=dptr->YD;


    for (i=0; i<n; i++){
        z[i] = p[0] + p[1]*x[i] + p[2]*(x[i]*x[i]) + p[3]*y[i] + p[4]*(y[i]*y[i]) + p[5]*x[i]*y[i];
        } // i loops over the number of data points
}

/* Analytic jacobian matrix of the model */
void jacquadsurf(double *p, double *jac,int m, int n, void *adata)
{
    register int i, j;
    double *x, *y;

    struct MyData *dptr;
    dptr=(struct mydata *)adata; // adata is passed as void *, dptr = data pointer

    x=dptr->XD;
    y=dptr->YD;


    /* The jacobian is filled row by row */
    for (i=j=0; i<n; i++){
        jac[j++] = 1.0 ;      // df/dp0
        jac[j++] = x[i];      // df/dp1
        jac[j++] = x[i]*x[i]; // df/dp2
        jac[j++] = y[i];      // df/dp3
        jac[j++] = y[i]*y[i]; // df/dp4
        jac[j++] = x[i]*y[i]; // df/dp5


        }
}


int Compute_Curvature ( struct inp_curvature *inp_curvature, CoorSet *trj_crd, char *outstring)
{
    int n; //n is the number of points and m is the number of parameters
    const int m=6; // there are m=6 parameters to fit
    
    register int i; //here comes the register int again ; )
    
    int ret;

    double *X = NULL, *Y = NULL, *Z = NULL; //data arrays
    struct MyData data;
    double p[m]; //parameter array
    
    double mean_c, gauss_c; //mean and gaussian curvatures
    
    double opts[LM_OPTS_SZ], info[LM_INFO_SZ]; //levmar related variables
    int niter=1000; //max number of iterations

    Selection sele = inp_curvature->sele;
    n = sele.nselatm; 
    
    /* Memory allocation for the coordinate arrays*/
    X = malloc(n*sizeof(double));
    Y = malloc(n*sizeof(double));
    Z = malloc(n*sizeof(double));
    
    /* Read atomic coordinates and store them in the X,Y, Z arrays */
    
    for (i=0;i<n;i++) {
        X[i] = trj_crd->xcoor[sele.selatm[i]];
        Y[i] = trj_crd->ycoor[sele.selatm[i]];
        Z[i] = trj_crd->zcoor[sele.selatm[i]];
    }
    
    
    /*X & Y are stored in the special purpose MyData structure;
     * this is the correct way to pass them to the optimizer.
     * Please let me know if you find a clever way to explain this.
     */
     
    data.XD = X;
    data.YD = Y; 
    
    /* We need a first guess of the parameters to initialize the
     * iterative procedure.
     * 
     * For now, they are actual guesses: I have no f**king idea of the
     * expected order of magnitude you expect for the curvature of 
     * a beta-sheet (yeah let's be honest, this function is all about
     * beta sheets).
     * 
     * Once I have a working code I will run it on selected beta-sheets 
     * from the PDB (read: I will randomly grab a couple of myosin structures)
     * to have better initial estimates.
     *
     * Another, probably better solution -but slightly more expensive- 
     * would be to find a clever way to do a quick estimate of the parameters
     * based on the atomic coordinates.
     * 
     * I don't have enough experience with the LM algorithm to evaluate the
     * influence of the guess. My feeling is that a bad guess will essentially
     * slow down convergence without affecting too much the final result. A boy can dream.
     */
     
     p[0] = 0.1;
     p[1] = 0.1;
     p[2] = 0.1;
     p[3] = 0.1;
     p[4] = 0.1;
     p[5] = 0.1;
     
    /* All right folks, coordinates are stored, 
     * fasten your seat belts for HERE COMES THE OPTIMIZER ! */
     
     
     /* optimization control parameters; passing to levmar NULL instead of opts reverts to defaults */
    opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
    opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used
    
    ret=dlevmar_der(quadsurf, jacquadsurf, p, Z, m, n, niter, opts, info, NULL, NULL, &data); // with analytic Jacobian
    
    // The relevant parameters for calculating the curvatures are p[2], p[4] and p[5]
    
    // The mean curvature is half the trace of the hessian matrix
    // The gaussian curvature is the determinant of the hessian matrix 
    
    mean_c = 0.5*(p[2]+p[4]);
    gauss_c = 4.*p[2]*p[4] - (p[5]*p[5]);
    
    //Please note that mean_c has dimension [1/L] (here, probably 1/Å) 
    // whereas gauss_c has dimension [1/L²] 
    
    
    // Print the results
    sprintf( outstring, " %10.5f  %10.5f", mean_c, gauss_c);


    return 24;
}
