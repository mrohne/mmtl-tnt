/***********************************************************/
/***                                                     ***/ 
/***    Capacitance Calculation using Coifman wavelets   ***/
/***            for the method of moments                ***/
/***                                                     ***/
/***                  By Zhichao Zhang                   ***/
/***                      2001.6                         ***/
/***                                                     ***/
/***********************************************************/

#include "stdafx.h"
#ifdef HAVE_SYS_TIMES_H
#include <sys/times.h>
#endif
#include <time.h>
#include <stdio.h>

/********************** Global variables *******************/

   int     Nc, Nd;
   int     Die, Rec, Cir, Gnd, Tra;

   Matrix  Diedi;
   Matrix  Gnddi;
   Matrix  Recdi;
   Matrix  Cirdi;
   Matrix  Tradix;
   Matrix  Tradiy;

   Vector  a, b, c;
   Vector  at, bt, ct;

   Vector  Circum;

   Vector  Diec;
   Vector  Recc;
   Vector  Circ;
   Vector  Trac;
   Vector  Gndc;

   double  eps0;
   Vector  eps;
   Vector  d;

   int     matr, Np;
   int     Nh, Nc4, Nit, J;
   int     Nw, Nwx, Nwy, Nws;
   double  Pi, Pih, Pit;
   double  step_w, hw, hwh;
   double  EPS;
   Matrix  uns;
   Vector  Iq_sq;
   Complex Ic;

   int     od;
   double  maxga;
   double  maxerr = -0.0;


/********************** Main program ***********************/

int main ( int argc, char *argv[] )
{
  
  /************ to setup global parameters *********/
    
  char fileIn[1024];
  char fileName[1024], fileOut[1024], strg[20];
   
  //---------------------------------------------------------------
  // Get the input and output file names from the command-line.
  //---------------------------------------------------------------
  if ( argc > 0)
    strcpy (fileName, argv[1]);    
  else
    strcpy (fileName, "data");
  
  if ( strlen(fileName) < 2 )
    strcpy (fileName, "data");
  
  sprintf (fileIn, "%s.in", fileName);
  sprintf (fileOut, "%s.out", fileName);
  std::cout << "Input file : " << fileIn << std::endl;
  std::cout << "Output file: " << fileOut << std::endl;
  std::cout << " " << std::endl;
    
    /************ to setup global parameters *********/
    
    getparam(fileIn);
    
    /************ global variables *******************/

    Matrix       Cap(Nc,Nc);

    CmplxVector  sumc(Nc);    
    CmplxVector  Rhsw(Nc*Nw), Jnj(Nc*Nw), xc;
    CmplxMatrix  Pw(Nc*Nw,Nc*Nw);
    
    int          i, j, i1, j1, ix, jx, iy, m, n, k;
    double       tst, T;
    double       dlu;
    double       t, tp, xt, yt;
    Complex      tmp, tmp0, sum, sum0;
    Vector       xcnt(Nws), ycnt(Nws);
    Vector       twx(Nwx), wwx(Nwx), twy(Nwy), wwy(Nwy);
    Vector       ptwx(Nwx), whwx(Nwx), ptwy(Nwy), whwy(Nwy);
    Vector       sxc4(Nc4), sc4(Nc4), wxc4(Nc4), wc4(Nc4);
    Vector       tn(Nw);
    Vector       Indx;
    Matrix       uns(3,Nw);
    CmplxMatrix  lbd, si;

    lbd.resize(Nd*(Nd+1),4*od);
    si.resize(Nd*(Nd+1),4*od);
    
    std::ofstream  outdata( fileOut, std::ios::out );

    /*
    ofstream  out1( "charge1.m", ios::out );
    ofstream  out2( "charge2.m", ios::out );
    ofstream  out3( "charge3.m", ios::out );
    */

    std::ofstream  sysmat( "sysmat.dat",std::ios::out );

    getcoif4( sxc4, sc4, wxc4, wc4, Nit );
    getsbs( uns );
    getcoe( lbd, si );
    //    std::cout << "Approximation Coefficients Determined." << std::endl;

    /** points and weights for Gauss the method */
    gauleg( -1.0, 1.0, ptwx, whwx );
    gauleg( -1.0, 1.0, ptwy, whwy );
    
    for( i = 0; i < Nw; i++ ) {
      
      tn[i] = i * step_w;
      
    }
    
    /***** to create system matrix **************/
    
    for( i1 = 0; i1 < Nc; i1 ++ ) {
      
      //      std::cout << i1 << std::endl;
      
      for( i = 0; i < Nw; i ++ ) {
	
	for( j1 = 0; j1 < Nc; j1 ++ ) {
	  
	  for( j = 0; j < Nw; j ++ ) {
	    
	    if ( i1 == j1 && i == j ) {
	      
	      for( m = 0; m < Nws; m++ ) {
		
		xcnt[m] = uns(1,j) + hw * (double(m) + 0.5);
		ycnt[m] = uns(1,i) + hw * (double(m) + 0.5);
		
	      }
	      
	      sum = 0.0;
	      
	      for( iy = 0; iy < Nws; iy++ ) {
		
		for( jx = 0; jx < Nws; jx++ ) {
		  
		  intpnt( xcnt[jx], twx, wwx, ptwx, whwx,
			  ycnt[iy], twy, wwy, ptwy, whwy,
			  hwh );
		  
		  tmp = 0.0;
		  
		  for( m = 0; m < Nwy; m++ ) {
		    
		    tmp0 = 0.0;
		    
		    for( n = 0; n < Nwx; n++ ) {
		      
		      tmp0 += kernel( twy[m], twx[n], i1, j1, lbd, si ) *
			wwx[n] * scalcoif4( twx[n], J, j, sc4 );
		      
		    }
		    
		    tmp += tmp0 * wwy[m] * scalcoif4( twy[m], J, i, sc4 );
		    
		  }
		  
		  sum += tmp;
		  
		}
		
	      }
	      
	      sum = sum * Circum[j1] * Circum[i1];
	      
	    } else {
	      
	      sum = kernel( tn[i], tn[j], i1, j1, lbd, si ) * step_w 
		* Circum[j1] *Circum[i1];
	      
	    }
	    
	    Pw(i1*Nw+i,j1*Nw+j) = sum;

	    sysmat << i1*Nw+i << " " << j1*Nw+j << " " << Pw(i1*Nw+i,j1*Nw+j) << std::endl;
	  }
	}
      }
    }

    std::cout << Nc*Nw
	 << "X"
	 << Nc*Nw
	 << " System Matrix Created !!! "
         << std::endl;
    
    /* Solve equation for different voltage combination */
    
    for( i = 0; i < Nc; i++ ) {

      for ( j = 0; j < Nc; j++ ) {

	if( j == i ) {
	  Iq_sq[j] = 1.0;
	} else {
	  Iq_sq[j] = 0.0;
	}

	for ( k = 0; k < Nw; k++ ) {
	  Rhsw[j*Nw+k] = Circum[j] * Iq_sq[j] * sqrt(step_w);
	}

      }

      /*  By Using LU Decomposition  */
      if( matr == 1 ) {	
	ludcmp( Pw, Indx, dlu );
	lubksb( Pw, Rhsw, Indx );
      }
      
      /*  By Using Bi_Cgstab */
      if( matr == 2 ) {
	xc.resize(Nc*Nw);
	bi_cgstab_c( Pw, xc, Rhsw, EPS );
	Rhsw = xc;
      }
      
      std::cout << i+1
	   << " of "
	   << Nc
	   << " linear equations for N = "
	   << Nc*Nw
	   << " solved !!! "
	   << std::endl;
 
      /* To calculate the capacitance ****/
      
      for( ix = 0; ix < Nc; ix++ ) {
	
	for( jx = 0; jx < Nw; jx++ ) {
	  
	  tmp = 0.0;
	  
	  for( j = jx-3; j <= jx+3; j++ ) {
	    
	    if( j < 0 ) {
	      tst = tn[jx] + 1.0;
	      j1 = j + Nw;
	    } else if ( j > Nw-1) {
	      tst = tn[jx] - 1.0;
	      j1 = j - Nw;
	    } else {
	      tst = tn[jx];
	      j1 = j;
	    }
	    
	    tmp += scalcoif4( tst, J, j1, sc4 ) * Rhsw[ix*Nw+j1];
	    
	  }
	  
	  Jnj[ix*Nw+jx] = tmp;
	  
	}
	
      }
      
      for ( i1 = 0; i1 < Nc; i1 ++ ) {
	
	sumc[i1] = cmplx( 0.0, 0.0 );

	
	for ( j1 = 0; j1 < Nw; j1 ++ ) {

	  /*	    
	  if( i == 0 && i1 == 1 ) {

	    out1 << "c(" << j1+1 << ")=" << cabs(Jnj[i1*Nw+j1]) << ";" <<std::endl;
	    
	  }

	  if( i == 1 && i1 == 1 ) {

	    out2 << "c(" << j1+1 << ")=" << cabs(Jnj[i1*Nw+j1]) << ";" <<std::endl;
	    
	  }

	  if( i == 2 && i1 == 1 ) {

	    out3 << "c(" << j1+1 << ")=" << cabs(Jnj[i1*Nw+j1]) << ";" << std::endl;
	    
	  }
	  */

	  sumc[i1] += Jnj[i1*Nw+j1];
	  
	}
	
	Cap(i1,i) = real(sumc[i1]) * Circum[i1] / double(Nw);
	
      }
      
    }

    //----------------------------------------------------------------
    // Output the capacitance values.
    //----------------------------------------------------------------
    outdata << "---------------------------------------" << std::endl;
    outdata << "Capacitance Values" << std::endl;
    outdata << "Input file : " << fileIn << std::endl;
    outdata << "---------------------------------------" << std::endl;
    //    outdata << Cap*1.0e12 << std::endl;
    for( i = 0; i < Nc; i++ ) {
      
      for( j = 0; j < Nc; j++ ) {
	
	outdata << Cap(i,j)*1.0e12 << "     ";
	std::cout << Cap(i,j)*1.0e12 << "     ";

      }
      
      outdata <<std::endl;
      std::cout <<std::endl;
      
    }

    //    std::cout << Cap*1.0e12 << std::endl;
    
    std::cout << "Capacitance Calculation Completed." << std::endl;
    std::cout << "Estimated Maximum Computational Error is: " << maxerr << std::endl;
    outdata << "Estimated Maximum Computational Error is: " << maxerr << std::endl;

#ifdef HAVE_SYS_TIMES_H
    tms run_times;

    times(&run_times);
#endif

    return 0;
    
}







