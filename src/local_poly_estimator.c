#include <R_ext/Lapack.h>



#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Lapack.h>

int n, d;//, length_band;

double factorial (double x)
{
    if (x == 0)
        return 1;
    else
        return (x * factorial(x-1));
} 


// calculates beta = (a'a)^-1a'b and returns the result in b
void reg(const int* m, const int* n,
         const int* nrhs, double* a, const int* lda,
         double* b, const int* ldb,
         double* work, const int* lwork, int* info)
{
    const char* trans = "N";
    F77_CALL(dgels)(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
}



double K(char * kernel_type, double ponto)
{
    if (strcmp(kernel_type,"box") == 0)               // --------------------------------- box
        return(dunif(ponto,-1,1,0));
    else if (strcmp(kernel_type,"trun.normal") == 0)  // --------------------------------- Truncated Normal
        return((dnorm(ponto,0,1,0)/(pnorm(-3,0,1,0,0) - pnorm(3,0,1,0,0)))*2*dunif(ponto,-3,3,0));
    else if (strcmp(kernel_type,"gaussian") == 0)     // --------------------------------- Gaussian
        return(dnorm(ponto,0,1,0));
    else if (strcmp(kernel_type,"epanech") == 0)      // --------------------------------- Epanech
        return((2*3*dunif(ponto,-1,1,0)*(1-pow(ponto,2))/4));
    else if (strcmp(kernel_type,"biweight") == 0)     // --------------------------------- biweight (Quartic)
        return((2*15*dunif(ponto,-1,1,0)*pow(1-pow(ponto,2),2)/16));
    else if (strcmp(kernel_type,"triweight") == 0)    // --------------------------------- triweight
        return((2*35*dunif(ponto,-1,1,0)*pow(1-pow(ponto,2),3)/32));
    else if (strcmp(kernel_type,"triangular") == 0)   // --------------------------------- triangular
        return((2*dunif(ponto,-1,1,0)*(1-fabs(ponto))));
    
}


// generalized cross validation min_{h} nsum[Y - Y_hat(i)_{-i}]^2/(n-tr(H))^2
// or simple leave-one-out cross-validation
double GCV( double *X, double *Y, int n, int d, char * kernel_type, double *bands, int grid_size, int degree_pol, int deriv, double *p_bandwidth)   
{
    int i, j, k;
    
    // variables used to solve (X'X)^-1X'Y
    const int m_c = n; 
    int n2_c;
    if (degree_pol == 1)
        n2_c = 1 + 1;
    else if (degree_pol == 2)
        n2_c = 1+1 + 1*(1+1)/2; 
    else
        n2_c = degree_pol + 1;  
    double a_c[n2_c*n]; // this will be X
    double b_c[n]; // this will be Y
    double a_copy[n2_c*n]; 
    double b_copy[n]; 
    const int nrhs_c = 1; 
    const int lda_c = n;
    const int ldb_c = n;
    int lwork_c;
    int mn_c = m_c;
    if (n2_c < m_c)
        mn_c = n2_c;
    if (mn_c == 1)
        lwork_c = mn_c + 1;
    else
        lwork_c = mn_c + mn_c;
    int info_c = 0;
    double work_c[lwork_c]; 
    for(i = 0; i < lwork_c; i++)
        work_c[i] = 0;
    
    int ind_cov;
    int method = p_bandwidth[0];
    int ind_band;    
    double res[grid_size];
    double aux;
    double H_ii;
    for (ind_cov = 0; ind_cov < d; ind_cov++)
    {
        for (ind_band = 0; ind_band < grid_size; ind_band++)
        {
            H_ii = 0;
            res[ind_band] = 0;
            for (i = 0; i < n; i++) // for each row/observation
            {
                    aux = X[i*d + ind_cov];
            
                // construct a = sqrt(W)XX                                 
                for (j = 0; j < n; j++)                                            
                {                                                                  
                    a_c[j] = sqrt(K(kernel_type, (X[j*d + ind_cov]-aux)/bands[ind_band*d + ind_cov]));
                    a_copy[j] = sqrt(K(kernel_type, (X[j*d + ind_cov]-aux)/bands[ind_band*d + ind_cov]));
                
                    if ((degree_pol >= 1)) // works only for d == 1
                        for (k = 1; k <= degree_pol; k++)                          
                        {
                            a_c[j+n*k] = pow((X[j*d + ind_cov]-aux),k)*a_c[j];                
                            a_copy[j+n*k] = pow((X[j*d+ind_cov]-aux),k)*a_copy[j];                
                        }
                
                    b_c[j] = Y[j]*a_c[j]; // b = sqrt(W)Y                           
                    if (i == j)
                        b_copy[j] = 1*a_copy[j]; // b_copy and a_copy are used to find H_ii: setting Y_copy = (0,0...,1,0,...0)
                    else
                        b_copy[j] = 0; 
                }
            
                for(j = 0; j < lwork_c; j++)
                    work_c[j] = 0;
                reg(&m_c, &n2_c, &nrhs_c, a_copy, &lda_c, b_copy, &ldb_c, work_c, &lwork_c, &info_c); // just to get H_ii 
                H_ii = H_ii + factorial(deriv)*b_copy[deriv];
            
            
                for(j = 0; j < lwork_c; j++)
                    work_c[j] = 0;
                reg(&m_c, &n2_c, &nrhs_c, a_c, &lda_c, b_c, &ldb_c, work_c, &lwork_c, &info_c);
            
                if (method == -1) // GCV
                    res[ind_band] = res[ind_band] + pow(Y[i] - factorial(deriv)*b_c[deriv],2);
                else // CV: method == 0
                    res[ind_band] = res[ind_band] + pow(Y[i] - factorial(deriv)*b_c[deriv],2)/pow(1-factorial(deriv)*b_copy[deriv],2);
            }
            if (method == -1) // GCV
                res[ind_band] = res[ind_band]/pow(1-H_ii/n,2);
        }        
    
        int ans;
        ans = 0;
        // select the bandwidth that has minimum res, and not NAN
        for (ind_band = 1; ind_band < grid_size; ind_band++)
            if (isnan(res[ind_band-1]))
                ans = ind_band;
            else
                if (res[ind_band] < res[ans])
                    ans = ind_band;
    
        p_bandwidth[ind_cov] = bands[ans*d + ind_cov];
    }        
}


// generalized cross validation min_{h} nsum[Y - Y_hat(i)_{-i}]^2/(n-tr(H))^2
// or simple leave-one-out cross-validation
void GCV_each_dimens( double *X, double *Y, int n, int d, char * kernel_type, double *bands, int grid_size, int degree_pol, int deriv, double* p_bandwidth)   
{
    int i, j, k;
    
    // variables used to solve (X'X)^-1X'Y
    const int m_c = n; 
    int n2_c;
    if (degree_pol == 1)
        n2_c = 1 + d;
    else if (degree_pol == 2)
        n2_c = 1+d + d*(d+1)/2; 
    else
        n2_c = degree_pol + 1;  
    double a_c[n2_c*n]; // this will be X
    double b_c[n]; // this will be Y
    double a_copy[n2_c*n]; 
    double b_copy[n]; 
    const int nrhs_c = 1; 
    const int lda_c = n;
    const int ldb_c = n;
    int lwork_c;
    int mn_c = m_c;
    if (n2_c < m_c)
        mn_c = n2_c;
    if (mn_c == 1)
        lwork_c = mn_c + 1;
    else
        lwork_c = mn_c + mn_c;
    int info_c = 0;
    double work_c[lwork_c]; 
    for(i = 0; i < lwork_c; i++)
        work_c[i] = 0;
    
    // this is for the 
    long limit = 1;
    for (i = 1; i <= d; i++)
        limit = limit*grid_size;    
    int combination_indices[d];
    long temp;
    
    long ind_band;    
    double res[limit];
    double aux[d];
    double H_ii;
    for (ind_band = 0; ind_band < limit; ind_band++)
    {
        temp = ind_band;                               // dynamic number of loops implemented by finding the position of the 
        for (i = 0; i < d; ++i)                        // matrix by remainder and division.
        {                                              // Each position of "combination_indices" has the line(entry) of the 
            combination_indices[i] = temp % grid_size; // corresponding column, making the combination of bandwidths for
            temp = temp/grid_size;                     // the kernel regression. Objective: find the combination with min MSE
        }
        
        H_ii = 0;
        res[ind_band] = 0;
        for (i = 0; i < n; i++) // for each row/observation
        {
            // aux is the point where we estimate m1
            if (d == 1)
                aux[0] = X[i];
            else
            {
                for (j = 0; j < d; j++)
                    aux[j] = X[i*d + j];
            }
            
            // construct a = sqrt(W)XX                                 
            for (j = 0; j < n; j++)                                            
            {                                                                  
                
                a_c[j] = 1;
                a_copy[j] = 1;
                for (k = 1; k <= d; k++)                             
                {
                    a_c[j] = a_c[j]*sqrt(K(kernel_type, (X[j*d + k-1]-aux[k-1])/bands[combination_indices[k-1]*d + k-1]));
                    a_copy[j] = a_copy[j]*sqrt(K(kernel_type, (X[j*d + k-1]-aux[k-1])/bands[combination_indices[k-1]*d + k-1]));
                }
                
                
                if ((degree_pol == 1) || (degree_pol == 2))          // add columns X1-x, X2-x,... Xd-x
                    for (k = 1; k <= d; k++)                             
                    {
                        a_c[j+n*k] = (X[j*d + k-1]-aux[k-1])*a_c[j];    // note that a is transpose manner
                        a_copy[j+n*k] = (X[j*d + k-1]-aux[k-1])*a_copy[j];    // note that a is transpose manner
                    }                
                
                if (degree_pol == 2) // include columns of half vectorization: VECH
                {
                    int l, ind_vech;
                    ind_vech = 1;
                    for (k = 1; k <= d; k++)       
                        for (l = k; l <= d; l++)
                        {
                            a_c[j+n*d+n*ind_vech] = (X[j*d + k-1]-aux[k-1])*(X[j*d + l-1]-aux[l-1])*a_c[j];    
                            a_copy[j+n*d+n*ind_vech] = (X[j*d + k-1]-aux[k-1])*(X[j*d + l-1]-aux[l-1])*a_copy[j];    
                            ind_vech = ind_vech + 1;
                        }
                }
                
                if ((degree_pol > 2) && (d == 1)) // works only for d == 1
                    for (k = 1; k <= degree_pol; k++)                          
                    {
                        a_c[j+n*k] = pow((X[j]-aux[0]),k)*a_c[j];                
                        a_copy[j+n*k] = pow((X[j]-aux[0]),k)*a_copy[j];                
                    }
                
                b_c[j] = Y[j]*a_c[j]; // b = sqrt(W)Y                           
                if (i == j)
                    b_copy[j] = 1*a_copy[j]; // b_copy and a_copy are used to find H_ii: setting Y_copy = (0,0...,1,0,...0)
                else
                    b_copy[j] = 0; 
            }
            
            for(j = 0; j < lwork_c; j++)
                work_c[j] = 0;
            reg(&m_c, &n2_c, &nrhs_c, a_copy, &lda_c, b_copy, &ldb_c, work_c, &lwork_c, &info_c); // just to get H_ii 
            H_ii = H_ii + factorial(deriv)*b_copy[deriv];
            
            
            for(j = 0; j < lwork_c; j++)
                work_c[j] = 0;
            reg(&m_c, &n2_c, &nrhs_c, a_c, &lda_c, b_c, &ldb_c, work_c, &lwork_c, &info_c);
            
            if (p_bandwidth[0] == -3) // GCV
                res[ind_band] = res[ind_band] + pow(Y[i] - factorial(deriv)*b_c[deriv],2);
            else // CV: method == 0
                res[ind_band] = res[ind_band] + pow(Y[i] - factorial(deriv)*b_c[deriv],2)/pow(1-factorial(deriv)*b_copy[deriv],2);
        }
        if (p_bandwidth[0] == -3) // GCV
            res[ind_band] = res[ind_band]/pow(1-H_ii/n,2);
    }        
    
    long ans = 0;
    // select the bandwidth that has minimum res, and not NAN
    while (isnan(res[ans]))
        ans = ans + 1;
    
    for (ind_band = ans+1; ind_band < limit; ind_band++)
        if (!isnan(res[ind_band]))
            if (res[ind_band] < res[ans])
                ans = ind_band;

    // decode the index whose MSE is min, to matrix entry and put the corresponding band on p_bandwidth to return
    temp = ans;
    for (k = 1; k <= d; k++) 
    {
        combination_indices[k-1] = temp % grid_size;
        p_bandwidth[k-1] = bands[combination_indices[k-1]*d + k-1];
        temp = temp/grid_size;
    }

    
}






/*----------------------------------------------------------------------------------------------------------------------------
 ------------------------------------------------------------------------------------------------------------------------------
 MAIN FUNCTION
 ------------------------------------------------------------------------------------------------------------------------------
 ------------------------------------------------------------------------------------------------------------------------------ */
SEXP local_poly_estimator(SEXP X, SEXP Y, SEXP points, SEXP band, SEXP grid1, SEXP degree_poly, SEXP kernel_type1, SEXP deriv1)   
{
    int i, j;
    i = 0;j = 0;
    
    
    /* Digest the datastructures (SEXPs) from R */ 
    double *xptr, *yptr, *grid;
    
    int degree_pol = INTEGER_VALUE(degree_poly);
    int deriv = INTEGER_VALUE(deriv1);
    PROTECT(grid1 = coerceVector (grid1, REALSXP) ) ; 
    grid = REAL(grid1);
    SEXP dimgrid = coerceVector(getAttrib(grid1, R_DimSymbol), INTSXP);
    int n_grid = INTEGER(dimgrid)[1];
    
    char *kernel_type;
    PROTECT(kernel_type1 = AS_CHARACTER(kernel_type1));
    kernel_type = R_alloc(strlen(CHAR(STRING_ELT(kernel_type1, 0))), sizeof(char)); 
    strcpy(kernel_type, CHAR(STRING_ELT(kernel_type1, 0)));
        
    
    // get dimensions of matrix X
    SEXP dimX = coerceVector(getAttrib(X, R_DimSymbol), INTSXP);
    d = INTEGER(dimX)[0];
    n = INTEGER(dimX)[1];
    
    // get dimensions of matrix points
    double *pontos;
    int n_pontos, d_pontos;
    SEXP dimpoints = coerceVector(getAttrib(points, R_DimSymbol), INTSXP);
    d_pontos = INTEGER(dimpoints)[0];
    n_pontos = INTEGER(dimpoints)[1];
    
    if ((d > 1) && (d_pontos == 1)) // X is a matrix n by d and points is a vector
    {                               // then, points is one point of d dimension
        n_pontos = 1;
        d_pontos = d;
    }
    
    PROTECT(X = coerceVector (X, REALSXP) ) ;
    xptr = REAL(X);
    PROTECT(Y = coerceVector (Y, REALSXP) ) ; 
    yptr = REAL(Y);
    PROTECT(points = coerceVector (points, REALSXP) ) ; 
    pontos = REAL(points);
    
    
    // aux is at each step the point x at which we predict y
    double aux[d];
    int k;
    
    
    // pred is the predicted values that will be returned
    SEXP pred;
    double *p_pred;
    PROTECT(pred = NEW_NUMERIC(n_pontos)); 
    p_pred = NUMERIC_POINTER(pred);

    
    PROTECT(band = coerceVector (band, REALSXP) ) ; 
    double * banda = REAL(band);
    // banda must have dimensions: n_points by d
    
    //SEXP dimbanda = coerceVector(getAttrib(band, R_DimSymbol), INTSXP);
    //int d_banda = INTEGER(dimbanda)[0];
    //int n_banda = INTEGER(dimbanda)[1];

    SEXP bandwidth;
    double *p_bandwidth;
    PROTECT(bandwidth = NEW_NUMERIC(d*n_pontos)); 
    p_bandwidth = NUMERIC_POINTER(bandwidth);

    
    // ------------------------------------------------------------- Cross Validation or GCV
    if ((banda[0] == 0) || (banda[0] == -1))
    {
        GCV(xptr, yptr, n , d , kernel_type, grid, n_grid, degree_pol, deriv, p_bandwidth);
        //for (i = 1; i < d; i++)
        //    p_bandwidth[i] = p_bandwidth[0];

        for (i = 1; i < n_pontos; i++)
            for (j = 0; j < d; j++)
                p_bandwidth[i*d + j] = p_bandwidth[j];
    } else
    // ------------------------------------------------------------- Cross Validation or GCV multidimensional
    if ((banda[0] == -2) || (banda[0] == -3))
    {
        GCV_each_dimens(xptr, yptr, n , d , kernel_type, grid, n_grid, degree_pol, deriv, p_bandwidth);
        for (i = 1; i < n_pontos; i++)
            for (j = 0; j < d; j++)
                p_bandwidth[i*d + j] = p_bandwidth[j];
    } else
    // ------------------------------------------------------------- 
    { // if no cross-validation, I still need to fill the matrix of bandwidths
      // where each row correspond to a point in 'points' sent by the user here
       for (i = 0; i < n_pontos; i++)
          for (j = 0; j < d; j++)
             p_bandwidth[i*d + j] = banda[i*d + j];
    }

    
    // variables used to solve (X'X)^-1X'Y
    const int m = n; 
    int n2;
    if (degree_pol == 1)
        n2 = 1 + d;
    else if (degree_pol == 2)
        n2 = 1+d + d*(d+1)/2; 
    else
        n2 = degree_pol + 1;  
    
    
    double a[n2*n]; // this will be X
    double b[n]; // this will be Y
    const int nrhs = 1; 
    const int lda = n;
    const int ldb = n;
    int lwork;
    int mn = m;
    if (n2 < m)
        mn = n2;
    if (mn == 1)
        lwork = mn + 1;
    else
        lwork = mn + mn;
    int info = 0;
    double work[lwork]; 
    for(i = 0; i < lwork; i++)
        work[i] = 0;
    
    // ------------------------------------------------------------------------------------------------- Prediction
    for (i = 0; i < n_pontos; i++)
    {
        
        // ------------------------------------ construct aux
        //aux is the point where m1 is to be estimated
        if (d == 1)
            aux[0] = pontos[i];
        else
            if (n_pontos == 1)           // here, X is a matrix n by d (d>1) 
                for (j = 0; j < d; j++)  //and points is a vector size d, thus there is 1 point
                    aux[j] = pontos[j];
            else
            {
                for (j = 0; j < d; j++)
                    aux[j] = pontos[i*d + j];
            }
        
        // for each observation in X construct a and obtain beta_hat_0 = m_hat(aux)
        for (j = 0; j < n; j++)                                            
        {      
            // construct a = sqrt(W)XX                                  
            a[j] = 1;
            for (k = 1; k <= d; k++)                             
                a[j] = a[j]*sqrt(K(kernel_type, (xptr[j*d + k-1]-aux[k-1])/p_bandwidth[i*d + k-1])); // for a vector of bandwidths
            
            
            if ((degree_pol == 1) || (degree_pol == 2))          // add columns X1-x, X2-x,... Xd-x
                for (k = 1; k <= d; k++)                             
                    a[j+n*k] = (xptr[j*d + k-1]-aux[k-1])*a[j];    // note that a is transpose manner
            
            
            if (degree_pol == 2) // include columns of half vectorization: VECH
            {
                int l, ind_vech;
                ind_vech = 1;
                for (k = 1; k <= d; k++)       
                    for (l = k; l <= d; l++)
                    {
                        a[j+n*d+n*ind_vech] = (xptr[j*d + k-1]-aux[k-1])*(xptr[j*d + l-1]-aux[l-1])*a[j];    
                        ind_vech = ind_vech + 1;
                    }
            }
            
            if ((degree_pol > 2) && (d == 1)) // works only for d == 1
                for (k = 1; k <= degree_pol; k++)                          
                    a[j+n*k] = pow((xptr[j]-aux[0]),k)*a[j];                
            
            
            b[j] = yptr[j]*a[j]; // b = sqrt(W)Y                           
        }
        
        // reg does (a'a)^-1a'b
        reg(&m, &n2, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
        p_pred[i] = factorial(deriv)*b[deriv];
        
    }
    // -------------------------------------------------------------------------------------------------- Prediction
    
    
    
    SEXP list, list_names;
    char *names[2] = {"predicted", "bandwidth"};
    PROTECT(list_names = allocVector(STRSXP,2));    
    PROTECT(list = allocVector(VECSXP, 2)); 
    for(i = 0; i < 2; i++)   
        SET_STRING_ELT(list_names,i,mkChar(names[i])); 
    SET_VECTOR_ELT(list, 0, pred); 
    SET_VECTOR_ELT(list, 1, bandwidth); 
    setAttrib(list, R_NamesSymbol, list_names); 
    
    UNPROTECT( 10 ) ;
    return(list);
}




