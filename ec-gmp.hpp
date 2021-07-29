
/*The implementation is based in the Python code written in 
https://github.com/Chia-Network/bls-signatures/tree/main/python-impl*/


#ifndef EC_BNUM_HPP
#define EC_BNUM_HPP






#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "gmp.h"
#include "gmpxx.h"

class AffinePoint
{

    //private:

        //static unsigned x;
    //    mpz_t x;
        //static unsigned y;
    //   mpz_t y;

    //    bool is_infinity;

    
    public:

        mpz_t x;
        mpz_t y;
        bool is_infinity;
    
    
        AffinePoint()
        {
            mpz_init(x); //x = 0;
            mpz_init(y); //y = 0;

            is_infinity = true;
        }
        
        //AffinePoint(bint X, bint Y)
        AffinePoint(mpz_t  X, mpz_t Y)
        {

	    mpz_t q;
            mpz_init_set_str(q,"4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10);

	        //bint r = "52435875175126190479447740508185965837690552500527637822603658699938581184513";
	        //mpz_t r;
            //mpz_init_set_str(r,"52435875175126190479447740508185965837690552500527637822603658699938581184513",10);
            
            //mpz_t x;
            mpz_init_set(x,X);
            
            mpz_mod(x, X, q); //X % q; //it must be reduced modulus R
            
            //mpz_t y;
            mpz_init_set(y,Y);
            mpz_mod(y, Y, q); //Y % q; //it must be reduced modulus R

            if (mpz_cmp_si(x, 0) == 0 && mpz_cmp_si(y, 0) == 0)
            {
                is_infinity = true;
            }
            else is_infinity = false;
        }
        
        void Set_x(mpz_t X)
        {
            mpz_set(x,X);
        }
        
        void Set_y(mpz_t Y)
        {
            mpz_set(y, Y);
        }

        /*mpz_t Get_x(void)
        {
            //string x_str;; 
            //mpz_get_str(x_str, 10, x);
            return x;
        }*/
        
        //string Get_y(void){return y;}
        /*mpz_t Get_y(void)
        {
            //string y_str;; 
            //mpz_get_str(y_str, 10, x);
            return y;
        }*/
        
        static bool is_on_curve(AffinePoint &P)
        {
            //Check that y^2 = x^3 + ax + b
            if (P.is_infinity)
                return true;

            mpz_t x_coordinate;
            mpz_init_set(x_coordinate, P.x);
            
            mpz_t y_coordinate;
            mpz_init_set(y_coordinate, P.y);

            //void mpz_mul(mpz_t rop, mpz_t x, mpz_t y)  // x*y ⇒ rop
            mpz_t left;
            mpz_init(left);
            mpz_pow_ui(left, y_coordinate, 2);
            //bint left = y_coordinate * y_coordinate;
            //bint left = y_coordinate * y_coordinate;
            
            mpz_t right;
            mpz_init(right);
            mpz_pow_ui(right, x_coordinate, 3);
            
            mpz_add_ui(right, right, 4);
            //unsigned right = P.x * P.x * P.x + P.ec.a * P.x + P.ec.b;
            //bint right = x_coordinate * x_coordinate * x_coordinate + 4;
            //bint right = x_coordinate * x_coordinate * x_coordinate + 4;

            return mpz_cmp(left, right) == 0;;
        }

        static AffinePoint AddPoints(AffinePoint &P1 , AffinePoint &P2)
        {
            if (P1.is_infinity)
                return P2;

            if (P2.is_infinity)
                return P1;

            mpz_t x1; 
            mpz_init_set(x1, P1.x);
            mpz_t y1;
            mpz_init_set(y1, P1.y);
            mpz_t x2;
            mpz_init_set(x2, P2.x);
            mpz_t y2;
            mpz_init_set(y2, P2.y);
            
            //if (mpz_cmp(x1,x2) == 0  && mpz_cmp(y1,y2) == 0)
            //    return AffineDoublePoint(P1);

            if (mpz_cmp(x1,x2) == 0)
            {
                if(mpz_cmp(y1,y2) == 0)
                    return AffineDoublePoint(P1);
                else
                    return AffinePoint();
            }


            mpz_t s;
            mpz_init(s);
            mpz_t sn;
            mpz_init(sn);
            mpz_t sd;
            mpz_init(sd);
            
            mpz_sub(sn, y2, y1);
            mpz_sub(sd, x2, x1);
            mpz_fdiv_q(s, sn, sd);
            //bint s = (y2 - y1) / (x2 - x1);

            mpz_t new_x;
            mpz_init(new_x);
            mpz_mul(new_x, s, s);
            mpz_sub(new_x, new_x, x1);
            mpz_sub(new_x, new_x, x2);
            //s * s - x1 - x2;
            
            mpz_t new_y;
            mpz_init(new_y);
            mpz_sub(new_y, x1, new_x);
            mpz_mul(new_y, new_y, s);
            mpz_sub(new_y, new_y, y1);
            //bint new_y = s * (x1 - new_x) - y1;

            return AffinePoint(new_x, new_y);

        }

        static AffinePoint AffineDoublePoint(AffinePoint &P)
        {

            //bint x_coordinate = P.Get_x();
            mpz_t x_coordinate;
            mpz_init_set(x_coordinate, P.x);
            //bint y_coordinate = P.Get_y();
            mpz_t y_coordinate;
            mpz_init_set(y_coordinate, P.y);
            //unsigned left = Fq(ec.q, 3) * x * x;//Fq class
            //bint left = 3 * x_coordinate * x_coordinate;
            mpz_t left;
            mpz_init(left);
            mpz_mul(left, x_coordinate, x_coordinate);
            mpz_mul_ui(left, left, 3);
            //left = left + ec.a; in this case a = 0
            //s = left / (Fq(ec.q, 2) * y)
            //bint s = left / 2 * y_coordinate;
            mpz_t s;
            mpz_init(s);
            mpz_mul_ui(s, y_coordinate,2);
            mpz_fdiv_q(s, left, s);
            

            //bint new_x = s * s - x_coordinate - x_coordinate;
            mpz_t new_x;
            mpz_init(new_x);
            mpz_mul(new_x, s, s);
            mpz_sub(new_x, new_x, x_coordinate);
            mpz_sub(new_x, new_x, x_coordinate);
            
            //bint new_y = s * (x_coordinate - new_x) - y_coordinate;
            mpz_t new_y; //= s * (x_coordinate - new_x) - y_coordinate;
            mpz_init(new_y);
            mpz_sub(new_y, x_coordinate, new_x);
            mpz_mul(new_y, new_y, s);
            mpz_sub(new_y, new_y, y_coordinate);


	        mpz_clear(x_coordinate);
	        mpz_clear(y_coordinate);
	        mpz_clear(s);
	        mpz_clear(left);
	        mpz_clear(new_x);
	        mpz_clear(new_y);
            return AffinePoint(new_x, new_y);

        }


        
        
        /*Little-endian version
        Q := 0
        R := P
        for i = 0 .. k:
            if(d_i = 1)
                Q = add(Q, R)
            R = double(R)
        return Q
        */
        
        static AffinePoint ScalarMult_Little_endian(mpz_t k , AffinePoint &P)
        {
            
            mpz_t scalar;
            mpz_init_set(scalar, k);
            
            if (P.is_infinity || mpz_cmp_si(scalar, 0) == 0)//k % q == 0)
                return AffinePoint();
            
            AffinePoint Result = AffinePoint();
            AffinePoint Addend = AffinePoint(P.x,P.y);
            while(mpz_cmp_si(scalar, 0) > 0 )  //scalar > 0)
            {
                
                //mpz_t scalar_mod;
                //mpz_init(scalar_mod);
                //mpz_mod_ui(scalar_mod, scalar, 2);
                //if(mpz_cmp_ui(scalar_mod, 1) == 0)
                
                //if(scalar % 2 == 1)
                if(mpz_tstbit(scalar, 0) == 1)// % 2 == 1)
                {
                    Result = AddPoints(Result, Addend);
                }
                Addend = AffineDoublePoint(Addend);
                //scalar = scalar / 2;
                mpz_fdiv_q_ui(scalar, scalar, 2);//scalar = scalar / 2;
            }
            mpz_clear(scalar);
            return Result;
        }
        
        

};



#endif //EC_HPP
