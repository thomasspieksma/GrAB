// -------------------------------------------------------------------
// -------------------------------------------------------------------
// Compute the energy lost to other bound states on a parabolic orbit.
// -------------------------------------------------------------------
// -------------------------------------------------------------------


// ----------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------
// Main reference:
//
// Ref. [1]: G. M. Tomaselli, T. F. M. Spieksma, and G. Bertone, "Dynamical friction in gravitational atoms",
// JCAP 07 (2023) 070, arXiv:2305.15460 [gr-qc]
//
// Specifically, Section 3, Eqs. (3.6).
// ----------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------

// Compile with
// -----------------------------------------------------------------------------------------------
// 1. g++ -Wall -I/opt/homebrew/Cellar/gsl/2.7.1/include -c CaptureCrossSection_bound_PUBLIC.cpp
// 2. g++ CaptureCrossSection_bound_PUBLIC.o -L/opt/homebrew/Cellar/gsl/2.7.1/lib -lgsl -lgslcblas
// 3. ./a.out
// -----------------------------------------------------------------------------------------------


// ----------------
// Import packages.
// ----------------

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <cmath>
#include <iostream>
#include <float.h>
#include <fstream>
#include "gsl/gsl_sf_coulomb.h"
#include "gsl/gsl_sf_coupling.h"
#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_sf_legendre.h"

// ----------------------------------------------------------
// We normalize by setting the Schwarzschild radius, r_s = 1.
// ----------------------------------------------------------


inline double factorial (int n)
{
    if (n == 0) return 1;
    else return n * factorial(n-1);
}


// -----------------------------------------
// Wigner small d-matrix (Eq. (6.2) in [1]).
// -----------------------------------------

double Wigner_small_d (int j, int mprime, int m, double beta)
{
    double x = 0;
    int smin = std::max(0, m-mprime);
    int smax = std::min(j+m, j-mprime);
    
    if (smin > smax) return 0;
    
    for (int s = smin; s <= smax; s++)
        x += (pow(-1,mprime-m+s) * pow(cos(beta/2),2*j+m-mprime-2*s) * pow(sin(beta/2),mprime-m+2*s) ) / (factorial(j+m-s)*factorial(s)*factorial(mprime-m+s)*factorial(j-mprime-s));
    
    return x * sqrt(factorial(j+mprime)*factorial(j-mprime)*factorial(j+m)*factorial(j-m));
}


// -------------------------------------------------------------------------------------------------------------------
// Radial function of the bound states (Eq. (2.7) in [1]). Normalized such that integral_0^\infty Rbound^2 r^2 dr = 1.
// -------------------------------------------------------------------------------------------------------------------

double Rbound (int n, int l, double r, double alpha, double rs)
{
    double rbohr = rs / (2 * pow(alpha,2));
            
    return gsl_sf_hydrogenicR(n,l,1/rbohr,r);
}


// --------------------------------------------------------------------------------------------------------------------------------------------
// Energy of bound states of the cloud (Eq. (2.9) in [1]).
// m_set_spin is meant to be the m of the most populated state. It is used to set the BH spin on the threshold of the superradiant instability.
// --------------------------------------------------------------------------------------------------------------------------------------------

double Energy_bound (int n, int l, int m, double alpha, double rs, int m_set_spin)
{
    double mu = 2*alpha/rs;
    
    double fnl = (double) 2/n - (double) 6/(2*l+1);
    double hl = 0;
    if (l > 0) hl = (double) 16/(2*l*(2*l+1)*(2*l+2));
    double atilde = (double) (m_set_spin/alpha) / (1+ pow(m_set_spin/(2*alpha),2));
    
    if (2*alpha > m_set_spin) atilde = 1;
    
    // Currently turned off the 5th order correction so we can use the same time integral for different values of mprime.
        
    return mu * ( - pow(alpha/n,2)/2 - pow(alpha/n,4)/8 + ( fnl*alpha + hl*atilde*m*alpha*alpha) * pow(alpha/n,3) );
}


// -------------------------------------------------------
//  Integrand for the radial integral (Eq. (2.18) in [1]).
// -------------------------------------------------------

struct integrandRadial {int lprime; int l_star; int l; int nprime; int n; double alpha; double rs;};

// ---------------------------------------------------------------------------------------------
// "Inner part" of the potential (Eqs. (2.15)-(2.16) in [1]), excluding the dipole (l_star = 1).
// ---------------------------------------------------------------------------------------------

double integrandInner (double r, void * p)
{
    struct integrandRadial * params = (struct integrandRadial *)p;
        int lprime = (params->lprime);
        int l_star = (params->l_star);
        int l = (params->l);
        int nprime = (params->nprime);
        int n = (params->n);
        double alpha = (params->alpha);
        double rs = (params->rs);
    
    return pow(r, l_star+2) * Rbound(nprime, lprime, r, alpha, rs) * Rbound(n, l, r, alpha, rs);
}


// ---------------------------------------------------------------------------------------------
// "Outer part" of the potential (Eqs. (2.15)-(2.16) in [1]), excluding the dipole (l_star = 1).
// ---------------------------------------------------------------------------------------------

double integrandOuter (double r, void * p)
{
    struct integrandRadial * params = (struct integrandRadial *)p;
        int lprime = (params->lprime);
        int l_star = (params->l_star);
        int l = (params->l);
        int nprime = (params->nprime);
        int n = (params->n);
        double alpha = (params->alpha);
        double rs = (params->rs);
    
    return Rbound(nprime, lprime, r, alpha, rs) * Rbound(n, l, r, alpha, rs) / pow(r, l_star-1);
}


// ------------------------------------------------------------------------------------------
// Dipole contribution of the potential (only has an outer part) (Eqs. (2.15)-(2.16) in [1]).
// ------------------------------------------------------------------------------------------

struct integrandRadial_dipole {int lprime; int l_star; int l; int nprime; int n; double R_star; double alpha; double rs;};

double integrandOuter_dipole (double r, void * p)
{
    struct integrandRadial_dipole * params = (struct integrandRadial_dipole *)p;
        int lprime = (params->lprime);
        int l_star = (params->l_star);
        int l = (params->l);
        int nprime = (params->nprime);
        int n = (params->n);
        double R_star = (params->R_star);
        double alpha = (params->alpha);
        double rs = (params->rs);
    
    return Rbound(nprime, lprime, r, alpha, rs) * Rbound(n, l, r, alpha, rs) * (- pow(r,l_star)/pow(R_star,l_star+1) + pow(R_star,l_star)/pow(r,l_star+1)) * pow(r,2);
}


// -------------------------------------------------------------------------------------------------------------
// Compute radial integral (Eq. (2.18) in [1]).
// We use the QAG adaptive integration algorithm in the gsl library.
// For more information, see https://www.gnu.org/software/gsl/doc/html/integration.html#qag-adaptive-integration
// -------------------------------------------------------------------------------------------------------------

struct Fourier_integrandRadial {gsl_integration_workspace * w; int lprime; int l_star; int l; int nprime; int n; int msecond; double R_p; double Omega_p; double alpha; double rs;};

double I_r_sine (double t, void * p)
{
    struct Fourier_integrandRadial * params = (struct Fourier_integrandRadial *)p;
        gsl_integration_workspace * w = (params->w);
        int lprime = (params->lprime);
        int l_star = (params->l_star);
        int l = (params->l);
        int nprime = (params->nprime);
        int n = (params->n);
        int msecond = (params->msecond);
        double R_p = (params->R_p);
        double Omega_p = (params->Omega_p);
        double alpha = (params->alpha);
        double rs = (params->rs);
    
    double y = pow( ( 3*Omega_p*t + sqrt(pow(3*Omega_p*t,2) + 4) ) / 2, 1./3 ); // See Eq. (3.3) in [1].
    
    double phi_star = 2*atan(y-1/y); // See Eq. (3.2) in [1].
    
    double R_star = R_p * (1 + pow(y-1/y,2)); // See Eq. (3.2) in [1].
    
    int max_subdivisions = 1000; // Sets the size of the subintervals in which the integration region is divided.
    
    double relativeError = 1e-6;
        
    if (l_star == 1) // Dipole contribution.
    {
        double outerIntegral, outerError;
        integrandRadial_dipole parameters = {lprime, l_star, l, nprime, n, R_star, alpha, rs,};
        
        gsl_function Outer;
        Outer.function = &integrandOuter_dipole;
        Outer.params = &parameters;
        
        gsl_integration_qagiu(&Outer, R_star, 0, relativeError, max_subdivisions, w, &outerIntegral, &outerError);
            
        // For large R_star, there is a huge exponential suppression making the integration numerically unfeasible. We set it to zero by hand if necessary.
        if(outerIntegral!=outerIntegral) outerIntegral = 0;
        
        return ( sqrt(rs) * outerIntegral ) * sin(msecond * phi_star);
    }
    
    else
    {
        double innerIntegral, innerError, outerIntegral, outerError;
        integrandRadial parameters = {lprime, l_star, l, nprime, n, alpha, rs,};
        
        gsl_function Inner;
        Inner.function = &integrandInner;
        Inner.params = &parameters;
        
        // The integration rule is chosen to be the 61 point Gauss-Kronrod rule.
        gsl_integration_qag(&Inner, 0, R_star, 0, relativeError, max_subdivisions, 6, w, &innerIntegral, &innerError);
        
        gsl_function Outer;
        Outer.function = &integrandOuter;
        Outer.params = &parameters;
        
        gsl_integration_qagiu(&Outer, R_star, 0, relativeError, max_subdivisions, w, &outerIntegral, &outerError);
            
        if(innerIntegral!=innerIntegral) innerIntegral = 0;
        if(outerIntegral!=outerIntegral) outerIntegral = 0;

        return ( sqrt(rs) * innerIntegral / pow(R_star, l_star+1) + sqrt(rs) * outerIntegral * pow(R_star, l_star) ) * sin(msecond * phi_star);
    }
}


double I_r_cosine (double t, void * p)
{
    struct Fourier_integrandRadial * params = (struct Fourier_integrandRadial *)p;
        gsl_integration_workspace * w = (params->w);
        int lprime = (params->lprime);
        int l_star = (params->l_star);
        int l = (params->l);
        int nprime = (params->nprime);
        int n = (params->n);
        int msecond = (params->msecond);
        double R_p = (params->R_p);
        double Omega_p = (params->Omega_p);
        double alpha = (params->alpha);
        double rs = (params->rs);
    
    double y = pow( ( 3*Omega_p*t + sqrt(pow(3*Omega_p*t,2) + 4) ) / 2, 1./3 );
    
    double phi_star = 2*atan(y-1/y);
    
    double R_star = R_p * (1 + pow(y-1/y,2));
    
    int max_subdivisions = 1000;
    
    double relativeError = 1e-6;
    
    if (l_star == 1)
    {
        double outerIntegral, outerError;
        integrandRadial_dipole parameters = {lprime, l_star, l, nprime, n, R_star, alpha, rs,};
        
        gsl_function Outer;
        Outer.function = &integrandOuter_dipole;
        Outer.params = &parameters;
        
        gsl_integration_qagiu(&Outer, R_star, 0, relativeError, max_subdivisions, w, &outerIntegral, &outerError);
            
        if(outerIntegral!=outerIntegral) outerIntegral = 0;
                    
        return ( sqrt(rs) * outerIntegral ) * cos(msecond * phi_star);
    }
    
    else
    {
        double innerIntegral, innerError, outerIntegral, outerError;
        integrandRadial parameters = {lprime, l_star, l, nprime, n, alpha, rs,};
        
        gsl_function Inner;
        Inner.function = &integrandInner;
        Inner.params = &parameters;
        
        gsl_integration_qag(&Inner, 0, R_star, 0, relativeError, max_subdivisions, 6, w, &innerIntegral, &innerError);
        
        gsl_function Outer;
        Outer.function = &integrandOuter;
        Outer.params = &parameters;
        
        gsl_integration_qagiu(&Outer, R_star, 0, relativeError, max_subdivisions, w, &outerIntegral, &outerError);
                    
        if(innerIntegral!=innerIntegral) innerIntegral = 0;
        if(outerIntegral!=outerIntegral) outerIntegral = 0;
                    
        return ( sqrt(rs) * innerIntegral / pow(R_star, l_star+1) + sqrt(rs) * outerIntegral * pow(R_star, l_star) ) * cos(msecond * phi_star);
    }
}


// ---------------------------------------------------------------------------
// Compute angular integral (Eq. (2.19) in [1]) in terms of Wigner-3j symbols.
// ---------------------------------------------------------------------------

double I_Omega (int lprime, int l_star, int l, int mprime, int m_star, int m)
{
    return pow(-1,mprime+m_star) * sqrt(((2*lprime+1)*(2*l_star+1)*(2*l+1))/(4*M_PI)) * gsl_sf_coupling_3j(2*lprime, 2*l_star, 2*l, 0, 0, 0) * gsl_sf_coupling_3j(2*lprime, 2*l_star, 2*l, -2*mprime, -2*m_star, 2*m);
}


// -------------------------------------------------------------------------------------------------------
// Relation binary's frequency and orbital radius. (There is an additional 2 compared to Kepler's formula.
// -------------------------------------------------------------------------------------------------------

double Omega_p (double R_p, double q, double rs)
{
    return sqrt(rs*(1+q)/(4*R_p*R_p*R_p));
}


// ------------------------------------------------------------------------------------------------------------------------------------------------------------
// Main function to compute the energy loss to other bound states (first term in Eq. (3.8) in [1]).
// This function returns the sum of energies lost towards states with same nprime and lprime (i.e., over all mprime). Later we then sum over nprime and lprime.
// The result is given divided by (1/2)q/(1+q).
// We use the QAWF adaptive integration algorithm for Foureir integrals in the gsl library.
// For more information, see https://www.gnu.org/software/gsl/doc/html/integration.html#qawo-adaptive-integration-for-oscillatory-functions
// ------------------------------------------------------------------------------------------------------------------------------------------------------------

double E_lost_bound (gsl_integration_workspace * w, gsl_integration_workspace * cycle_w, int nprime, int lprime, int n, int l, int m, double R_p, double alpha, double McOverM, double q, double chi, double rs)
{
    double epsabs = 1e-6;
        
    double Omega = Omega_p(R_p, q, rs);
    
    double Im_c_nlm[2*lprime+1]; // Excluding a factor 4 \pi \alpha q.
    
    double tmp[2*lprime+1];
    
    for (int i = 0; i < 2*lprime+1; i++)
        Im_c_nlm[i] = 0;
    
    bool do_integral;
    
    double result_sin, abserr_sin, result_cos, abserr_cos;
    
    gsl_integration_qawo_table * wf_sin, * wf_cos;
        
    if (nprime == n && lprime == l)
    {
        for (int l_star = abs(lprime-l); l_star <= l+lprime; l_star++) // Selection rules, see Eqs. (2.20)-(2.22) in [1].
        {
            if ( (lprime+l_star+l)%2 == 0 )
            {
                for (int msecond = -l_star; msecond <= l_star; msecond++)
                {
                    for (int mprime = -lprime; mprime <= lprime; mprime++)
                    {
                        do_integral = 0;

                        tmp[mprime+lprime] = I_Omega(lprime,l_star,l,mprime,m-mprime,m) * Wigner_small_d(l_star, m-mprime, msecond, chi);
                        
                        if (abs(tmp[mprime+lprime]) > 1e-15) // Check if the angular part is zero, to speed up the code.
                            do_integral = 1;
                        
                        if (do_integral)
                        {
                            double DeltaE = Energy_bound(nprime, lprime, mprime, alpha, rs, m) - Energy_bound(n, l, m, alpha, rs, m);
                            
                            wf_sin = gsl_integration_qawo_table_alloc(DeltaE, 0, GSL_INTEG_SINE, 10);
                            wf_cos = gsl_integration_qawo_table_alloc(DeltaE, 0, GSL_INTEG_COSINE, 10);
                            
                            Fourier_integrandRadial parameters = {w, lprime, l_star, l, nprime, n, msecond, R_p, Omega, alpha, rs,};
                            gsl_function Fourier_sine;
                            Fourier_sine.function = &I_r_sine;
                            Fourier_sine.params = &parameters;
                            
                            gsl_function Fourier_cosine;
                            Fourier_cosine.function = &I_r_cosine;
                            Fourier_cosine.params = &parameters;
                                        
                            gsl_integration_qawf(&Fourier_sine, 0, epsabs, 1000, w, cycle_w, wf_sin, &result_sin, &abserr_sin);
                                        
                            gsl_integration_qawf(&Fourier_cosine, 0, epsabs, 1000, w, cycle_w, wf_cos, &result_cos, &abserr_cos);
                            
                            Im_c_nlm[mprime+lprime] += - tmp[mprime+lprime] * gsl_sf_legendre_sphPlm(l_star,abs(msecond),0) * 2 * (result_cos - result_sin) / (2*l_star+1); // See Eq. (3.7) in [1].
                        }
                    }
                }
            }
        }
        
        double final_sum = 0;
        for (int mprime = -lprime; mprime <= lprime; mprime++)
            final_sum += pow(Im_c_nlm[mprime+lprime],2) * (Energy_bound(nprime, lprime, mprime, alpha, rs, m) - Energy_bound(n, l, m, alpha, rs, m)); // See Eq. (3.8) in [1].
        
        return pow(4*M_PI*alpha*q,2) * ( McOverM/( (2*alpha/rs) * 0.5 * q / (1+q) ) ) * final_sum; // See Eq. (3.8) in [1].
    }
    else
    {
        double DeltaE = Energy_bound(nprime, lprime, 0, alpha, rs, m) - Energy_bound(n, l, m, alpha, rs, m);
        
        wf_sin = gsl_integration_qawo_table_alloc(DeltaE, 0, GSL_INTEG_SINE, 10);
        wf_cos = gsl_integration_qawo_table_alloc(DeltaE, 0, GSL_INTEG_COSINE, 10);
                    
        for (int l_star = abs(lprime-l); l_star <= l+lprime; l_star++)
        {
            if ( (lprime+l_star+l)%2 == 0 )
            {
                for (int msecond = -l_star; msecond <= l_star; msecond++)
                {
                    do_integral = 0;

                    for (int mprime = -lprime; mprime <= lprime; mprime++)
                    {
                        tmp[mprime+lprime] = I_Omega(lprime,l_star,l,mprime,m-mprime,m) * Wigner_small_d(l_star, m-mprime, msecond, chi);
                        if (abs(tmp[mprime+lprime]) > 1e-15)
                            do_integral = 1;
                    }
                    
                    if (do_integral)
                    {
                        Fourier_integrandRadial parameters = {w, lprime, l_star, l, nprime, n, msecond, R_p, Omega, alpha, rs,};
                        gsl_function Fourier_sine;
                        Fourier_sine.function = &I_r_sine;
                        Fourier_sine.params = &parameters;
                        
                        gsl_function Fourier_cosine;
                        Fourier_cosine.function = &I_r_cosine;
                        Fourier_cosine.params = &parameters;
                                    
                        gsl_integration_qawf(&Fourier_sine, 0, epsabs, 1000, w, cycle_w, wf_sin, &result_sin, &abserr_sin);
                                    
                        gsl_integration_qawf(&Fourier_cosine, 0, epsabs, 1000, w, cycle_w, wf_cos, &result_cos, &abserr_cos);
                        
                        for (int mprime = -lprime; mprime <= lprime; mprime++)
                        {
                            Im_c_nlm[mprime+lprime] += - tmp[mprime+lprime] * gsl_sf_legendre_sphPlm(l_star,abs(msecond),0) * 2 * (result_cos - result_sin) / (2*l_star+1);
                        }
                    }
                }
            }
        }
        
        double final_sum = 0;
        for (int mprime = -lprime; mprime <= lprime; mprime++)
            final_sum += pow(Im_c_nlm[mprime+lprime],2);
        
        return pow(4*M_PI*alpha*q,2) * ( McOverM/( (2*alpha/rs) * 0.5 * q / (1+q) ) ) * DeltaE * final_sum;
    }
}


// ----------------------------------------------------------------------------------
// Specifying parameters, summing over (n_prime, l_prime) and outputting to txt file.
// ----------------------------------------------------------------------------------

int main(int argc, const char * argv[]) {
        
    gsl_set_error_handler_off();
        
    int max_subdivisions_K = 1000;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (max_subdivisions_K);
    gsl_integration_workspace * cycle_w = gsl_integration_workspace_alloc (max_subdivisions_K);
    
    std::ofstream file;
    file.open ("Elost_211_count_bound.txt");
    
    
    // -----------
    // Parameters.
    // -----------
    
    
    // ---------------------------------------------------------------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------------------------

    double tmp1, tmp2, tmp3;
    
    double E_lost_tot_bound = 0;
    
    double Rp = 0; // Periapsis of the orbit.

    double q = 0.001; // Mass ratio of the binary (Default: 1e-3).
    
    double alpha = 0.2; // Gravitational fine-structure constant (Default: 0.2).
        
    double chi = 0; // Set to 0 (M_PI) for co (counter)-rotating orbits.
    
    double McOverM = 0.01; // Cloud's mass (Default: 1% of the BH mass).
    
    int n_b = 2; // Initial (bound) state of the cloud, defined by the principal, angular momentum and azimuthal quantum numbers (Default: nlm = 211).
    
    int l_b = 1;
    
    int m_b = 1;
    
    int min_n_max = 4; // We loop over n_prime until one of the conditions* is satisfied.
    
    int max_n_max = 40;
    
    double max_fractional_increment_n = 1e-3;
    
    double max_fractional_increment_l = 1e-4;
    
    // ---------------------------------------------------------------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------------------------

         
    
    // -----------------------------------
    // Loop over the desired radial range.
    // -----------------------------------
    
    for (int i = 10; i < 4000; i++)
    {
        Rp = (double) i / 10;
        
        E_lost_tot_bound = 0;
        
        tmp1 = 0;
        tmp2 = 0;
        tmp3 = 0;
        
        for (int n = 1; n < max_n_max; n++)
        {
            tmp3 = tmp1;
            
            for (int l = 0; l < n; l++)
            {
                tmp2 = E_lost_bound(w,cycle_w,n,l,n_b,l_b,m_b,Rp,alpha,McOverM,q,chi,1);
                
                tmp1 += tmp2;
                
                // ---------------------------------------------------------------------------------------------------------------------
                // We break the loop over l_prime (per n_prime),
                // when enough modes have been summed such that the increment reaches "max_fractional_increment_lprime" (Default: 1e-4).
                // ---------------------------------------------------------------------------------------------------------------------
                
                if (l >= 3 and tmp2 != 0 and tmp1 != 0 and tmp3 != 0)
                    if ( abs(tmp2/(tmp1-tmp3)) < max_fractional_increment_l)
                        l = n;
            }
            
            if (n >= min_n_max and tmp3 != 0 and tmp1 != 0)
            {
                // ------------------------------------------------------------------------------------------------------------------------------
                // Specify the conditions* to break the loop over n_prime,
                // which is when enough modes have been summed such that the increment reaches "max_fractional_increment_nprime" (Default: 1e-3).
                // ------------------------------------------------------------------------------------------------------------------------------
                
                if ( abs(tmp1/tmp3-1) < max_fractional_increment_n )
                {
                    std::cout << "n = " << n << std::endl;
                    break;
                }
            }
        }
        E_lost_tot_bound = tmp1;
        
        std::cout << 2*Rp << " " << E_lost_tot_bound << std::endl;
        
        file << 2*Rp << " " << E_lost_tot_bound << std::endl;
    }
    
    gsl_integration_workspace_free(w);
    gsl_integration_workspace_free(cycle_w);
    
    return 0;
}
