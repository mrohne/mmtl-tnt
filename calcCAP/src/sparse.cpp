/******* functions, needed for sparse matrices  *******/

//#include "stdafx.h"
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>

#include "complex.h" 
#include "vector.h" 
#include "matrix.h" 
#include "cmplxvec.h" 
#include "cmplxmat.h"
#include "ivectorl.h"

#include "sparse.h"


void  sprsin_d( Matrix& a, double thresh, Vector& sa, IVectorl& ija ) {

    const int  ns = a.dim_i();
    int        i, j, k, nh;
     
    nh = ns + 1;

    for( i = 0; i < ns; i++ ) {

        for( j = 0; j < ns; j++ ) {

            if( fabs( a(i,j) ) > thresh && i != j )
               nh += 1;
        }
    }

    sa.resize( nh );
    ija.resize( nh );
    
    for( j = 0; j < ns; j++ )
        sa[j] = a(j,j);

    ija[0] = ns + 1;
    k = ns;

    for( i = 0; i < ns; i++ ) {

        for( j = 0; j < ns; j++ ) {

	    if( fabs( a(i,j) ) > thresh && i != j ) {

               k += 1;
               sa[k] = a(i,j);
               ija[k] = j;
            }               
        }

        ija[i+1] = k + 1;
    }

    a.resize(0,0);             
}    
   

void  sprsax_d( Vector& sa, IVectorl& ija, Vector& x, Vector& b ) {

    const int  ns = x.dim();
    int        i, k; 

    b.resize(ns); 

    if( ija[0] != ns + 1 ) {

       std::cerr << " Something is wrong in sprsax_d !!! " << std::endl;
       exit(1);
    }

    for( i = 0; i < ns; i++ ) {
        
      b[i] = sa[i] * x[i]; //std::cerr << i << std::endl;

        for( k = ija[i]; k <= ija[i+1]-1; k++ ) {

	  //std::cerr << k << std::endl;
 	    b[i] += sa[k] * x[ija[k]]; 
        }
    }
}


void  sprsin_c( CmplxMatrix& a, double thresh, CmplxVector& sa, IVectorl& ija ) {

    const int  ns = a.dim_i();
    int        i, j, k, nh;

    nh = ns + 1;

    for( i = 0; i < ns; i++ ) {

        for( j = 0; j < ns; j++ ) {

            if( cabs( a(i,j) ) > thresh && i != j )
               nh += 1;
        }
    }

    std::cout << " Matrix sparsity = "
         << double( nh ) / double( ns * ns ) * 100.0
         << " % "
         << std::endl; 
      
    sa.resize( nh );
    ija.resize( nh );
    
    for( j = 0; j < ns; j++ )
        sa[j] = a(j,j);

    ija[0] = ns + 1;
    k = ns;

    for( i = 0; i < ns; i++ ) {

        for( j = 0; j < ns; j++ ) {

	    if( cabs( a(i,j) ) > thresh && i != j ) {

               k += 1;
               sa[k] = a(i,j);
               ija[k] = j;
            }               
        }

        ija[i+1] = k + 1;
    }

    a.resize(0,0);             
}    
   

void  sprsax_c( CmplxVector& sa, IVectorl& ija, CmplxVector& x, CmplxVector& b ) {

    const int  ns = x.dim();
    int        i, k; 

    b.resize(ns); 

    if( ija[0] != ns + 1 ) {

       std::cerr << " Something is wrong in sprsax_c !!! " << std::endl;
       exit(1);
    }

    for( i = 0; i < ns; i++ ) {
        
        b[i] = sa[i] * x[i];

        for( k = ija[i]; k <= ija[i+1]-1; k++ ) {

	    b[i] += sa[k] * x[ija[k]];
        }
    }
}


Vector  sprsax_d_v( Vector& sa, IVectorl& ija, Vector& x ) {

    Vector res;

    sprsax_d( sa, ija, x, res );
 
    return( res );
}


CmplxVector  sprsax_c_v( CmplxVector& sa, IVectorl& ija, CmplxVector& x ) {

    CmplxVector res;

    sprsax_c( sa, ija, x, res );
 
    return( res );
}


void sprsin_simple( CmplxMatrix& a, double lev ) {

    for( int m = 0; m < a.dim_i(); m++ ) {

	for( int n = 0; n < a.dim_j(); n++ ) {
               
            if( cabs( a(m,n) ) < lev )
               a(m,n) = 0.0;
        }
    }
}


CmplxVector  sprsax_c_v_simple( CmplxMatrix& a, CmplxVector& v ) {

    const int    mi = a.dim_i();  
    const int    nj = a.dim_j();
    int          m, n;    
    CmplxVector  res(mi,0.0);

    for( m = 0; m < mi; m++ ) {

        for( n = 0; n < nj; n++ ) {

            if( a(m,n) != 0.0 )
               res[m] += a(m,n) * v[n];
        }
    }
     
    return( res );  
}
