/*
 *
 * Copyright 2015 Florian Blanc <flblanc@cecchini5>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 *
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "/usr/local/include/levmar.h"

#define MAX_LINE_LENGTH 1000

void quadsurf(double *p,  double *z, int m, int n, void *adata);
void jacquadsurf(double *p, double *jac,int m, int n, void *adata);


/* This structure allows us to pass the experimental data points
 * to the optimizer. It is weird but actually, it gives flexibility
 * because the user can design a struct that is adapted to the particular
 * problem (s)he is aiming to solve. 
 */
 

struct MyData {
    double *XD;
    double *YD;
    
    };


/* Model to be fitted to observed points */
void quadsurf(double *p,  double *z, int m, int n, void *adata)
{
    register int i;
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



int main(int argc, char **argv)
{

    int n; //n is the number of points and m is the number of parameters
    const int m=6;
    int ret;
    //register int i;
    double *X = NULL, *Y = NULL, *Z = NULL; //data arrays
    struct MyData data;
    double p[m]; //parameter array
    double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
    int niter=1000; //max number of iterations

    char *line;
    int k=0;



    // First I need to read the data from the files
    FILE* myfile=fopen("series.dat","r");
//    FILE* yfile=fopen("Y-series.dat","r");
//    FILE* zfile=fopen("Z-series.dat","r");


    line=malloc(MAX_LINE_LENGTH*sizeof(char));
    fgets(line,MAX_LINE_LENGTH,myfile);
    sscanf(line,"%d", &n); // the number of points is now stored in n

    printf("Found %d data points.\n", n);


    X = malloc(n*sizeof(double));
    Y = malloc(n*sizeof(double));
    Z = malloc(n*sizeof(double));

//    for (k=0;k<n;k++) {
//        printf("%10.5f %10.5f %10.5f \n",X[k],Y[k],Z[k]);}
//
//    printf("It's better if you allocate!\n");
//    printf("\n\n\n\n\n");

    //Reading data and storing them in the arrays
    while (fgets(line,MAX_LINE_LENGTH,myfile) != NULL)
			{
            sscanf(line,"%lf %lf %lf ",&X[k],&Y[k],&Z[k]);
            printf("%10.5f %10.5f %10.5f \n",X[k],Y[k],Z[k]);
            k++;
			}


    k=0 ; //reinit the counter
//    k=0; //reinitialize the counter
//    if (Y != NULL){
//        while (fgets(line,MAX_LINE_LENGTH,yfile) != NULL)
//			{
//
//				sscanf(line,"%lf",&Y[k]);
//				k++;
//			}
//    }
//    k=0; //reinitialize the counter
//    if (Z != NULL){
//        while (fgets(line,MAX_LINE_LENGTH,zfile) != NULL)
//			{
//
//				sscanf(line,"%lf",&Z[k]);
//				k++;
//			}
//    }

//    if (X != NULL) {
//        for (k=0;k<n;k++) {
//            if (fgets(line,MAX_LINE_LENGTH,myfile) != NULL) {
//
//                sscanf(line,"%lf %lf %lf ",&X[k],&Y[k],&Z[k]);}
//
//
//        }
//    }
    printf("Closing the file.\n");
    printf("\n\n\n\n\n");
    fclose(myfile);
//    fclose(yfile);
//    fclose(zfile);



    for (k=0;k<n;k++) {
        printf("%10.5f %10.5f %10.5f \n",X[k],Y[k],Z[k]);
        }


    data.XD = X;
    data.YD = Y;

//    exit(0);


    //Guess of the parameters
//     p[0]=0.78743658;
//     p[1]=0.25167087;
//     p[2]=0.32906466;
//     p[3]=0.26746894;
//     p[4]=0.14558755;
//     p[5]=0.45377611;
//
     p[0]=0.9;
     p[1]=0.9;
     p[2]=0.9;
     p[3]=0.9;
     p[4]=0.9;
     p[5]=0.9;

     /* optimization control parameters; passing to levmar NULL instead of opts reverts to defaults */
    opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
    opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used


    /* invoke the optimization function */
//    ret=dlevmar_der(quadsurf, jacquadsurf, p, Z, m, n, 1000, opts, info, NULL, NULL, &data); // with analytic Jacobian
    ret=dlevmar_der(quadsurf, jacquadsurf, p, Z, m, n, niter, opts, info, NULL, NULL, &data); // with analytic Jacobian
    //ret=dlevmar_dif(quadsurf, p, x, m, n, 1000, opts, info, NULL, NULL, NULL); // without Jacobian
    printf("Levenberg-Marquardt returned in %g iter, reason %g, sumsq %g [%g]\n", info[5], info[6], info[1], info[0]);
    printf("Best fit parameters: %.7g %.7g %.7g %.7g %.7g %.7g\n", p[0], p[1], p[2], p[3], p[4], p[5]);

    free(X);
    free(Y);
    free(Z);

	return 0;
}

