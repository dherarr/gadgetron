/** \file fatWaterOperator.h
    \brief Implement fat-water signal model operator (no R2* decay).
    \author  Diego Hernando
*/

#pragma once

#include "curveFittingOperator.h"
#include <boost/math/special_functions/sign.hpp>

namespace Gadgetron {

    // y = bi[0] * exp(-x/bi[1])

    template <class ARRAY> class fatWaterOperator : public curveFittingOperator<ARRAY>
    {
    public:

        typedef curveFittingOperator<ARRAY> BaseClass;
        typedef typename BaseClass::ELEMENT_TYPE ELEMENT_TYPE;
        typedef typename BaseClass::REAL REAL;

        fatWaterOperator();
        virtual ~fatWaterOperator();

        virtual void gradient(const ELEMENT_TYPE& xi, const ARRAY& b, ARRAY& grad);
        virtual void magnitude(const ARRAY& x, const ARRAY& b, ARRAY& y); //Note that magnitude in our case means real and imaginary components

    protected:
    };

    template <class ARRAY>
    fatWaterOperator<ARRAY>::fatWaterOperator() : BaseClass()
    {
    }

    template <class ARRAY>
    fatWaterOperator<ARRAY>::~fatWaterOperator()
    {
    }

    template <class ARRAY>
    void fatWaterOperator<ARRAY>::gradient(const ELEMENT_TYPE& curte, const ARRAY& b, ARRAY& grad)
    {
        try
        {

            size_t num = b.size();
            if(grad.size()!=num) grad.resize(num, 0);

            int sign_b1 = boost::math::sign(b[1]);
            ELEMENT_TYPE rb = 1.0 / ( (std::abs(b[1])<FLT_EPSILON) ? sign_b1*FLT_EPSILON : b[1] );

            ELEMENT_TYPE val = std::exp(-1 * xi * rb);
            grad[0] = val;
            grad[1] = b[0] * val * xi * rb * rb;


	    //	    expr2 = exp(-te[kt]*r2);
	    expr2 = 1.0;
	    sinfm = sin(2*PI*fieldmap*curte);
	    cosfm = cos(2*PI*fieldmap*curte);
	    
	    shatr=cosfm*expr2*(Wr*swr + Fr*sfr - Wi*swi - Fi*sfi) - sinfm*expr2*(Wr*swi + Fr*sfi + Wi*swr + Fi*sfr);
	    shati=sinfm*expr2*(Wr*swr + Fr*sfr - Wi*swi - Fi*sfi) + cosfm*expr2*(Wr*swi + Fr*sfi + Wi*swr + Fi*sfr);
	    
	    grad[0] = cosfm*expr2*swr - sinfm*expr2*swi; 
	    grad[1] = sinfm*expr2*swr + cosfm*expr2*swi;
	    

	    grad[2] = -cosfm*expr2*swi - sinfm*expr2*swr;
	    grad[3] = -sinfm*expr2*swi + cosfm*expr2*swr;
	    
	    grad[4] = cosfm*expr2*sfr - sinfm*expr2*sfi;
	    grad[5] = sinfm*expr2*sfr + cosfm*expr2*sfi;
	    
	    grad[6] = -cosfm*expr2*sfi[kt] - sinfm*expr2*sfr[kt];
	    grad[7] = -sinfm*expr2*sfi[kt] + cosfm*expr2*sfr[kt];
	    
	    /*     curJ5 = -te[kt]*shatr; */
	    /*     gsl_matrix_set (J, kt, 4, curJ5);  */
	    /*     curJ5 = -te[kt]*shati; */
	    /*     gsl_matrix_set (J, kt+nte, 4, curJ5);  */
	    
	    grad[8] = -2*PI*te*(sinfm*expr2*(Wr*swr + Fr*sfr - Wi*swi - Fi*sfi) + cosfm*expr2*(Wr*swi + Fr*sfi + Wi*swr + Fi*sfr));
	    grad[9] =  2*PI*te*(cosfm*expr2*(Wr*swr + Fr*sfr - Wi*swi - Fi*sfi) - sinfm*expr2*(Wr*swi + Fr*sfi + Wi*swr + Fi*sfr));
	    
	    /*    mexPrintf("%s%f%s%f%s%f\n", "j1: ", curJ1, ", j2: ", curJ2 , ", j3: ", curJ3);*/
	    
	    }


        }
        catch(...)
        {
            GADGET_THROW("Errors happened in fatWaterOperator<ARRAY>::gradient(...) ... ");
        }
    }

    template <class ARRAY>
    void fatWaterOperator<ARRAY>::magnitude(const ARRAY& te, const ARRAY& b, ARRAY& y)
    {
      
      ELEMENT_TYPE Wr = b[0];
      ELEMENT_TYPE Wi = b[1];
      ELEMENT_TYPE Fr = b[2];
      ELEMENT_TYPE Fi = b[3];
      ELEMENT_TYPE fieldmap = b[4];

      
      size_t num = x.size();
      if(y.size()!=x.size()) y.resize(num, 0);
      
      size_t kt;
      for kt=0; kt<num/2; kt++)
        {
	  //            y[ii] = b[0] * exp( -1 * x[ii] * rb);
	  // Need to write this code still - lots of missing stuff including reserving memory
	  y[kt] = cos(2*PI*fieldmap*te[kt])*(Wr*swr[kt] + Fr*sfr[kt] - Wi*swi[kt] - Fi*sfi[kt]) - sin(2*PI*fieldmap*te[kt])*(Wr*swi[kt] + Fr*sfi[kt] + Wi*swr[kt] + Fi*sfr[kt]);
	  y[num/2+kt] = sin(2*PI*fieldmap*te[kt])*(Wr*swr[kt] + Fr*sfr[kt] - Wi*swi[kt] - Fi*sfi[kt]) + cos(2*PI*fieldmap*te[kt])*(Wr*swi[kt] + Fr*sfi[kt] + Wi*swr[kt] + Fi*sfr[kt]);

	}
    }
}
