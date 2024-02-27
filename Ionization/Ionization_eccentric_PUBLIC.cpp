// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Compute ionization rate, power and torque on eccentric, co-rotating orbits.
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


// ----------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------
// Main reference:
//
// Ref. [1]: G. M. Tomaselli, T. F. M. Spieksma, and G. Bertone, "Dynamical friction in gravitational atoms",
// JCAP 07 (2023) 070, arXiv:2305.15460 [gr-qc]
//
// Specifically, Section 5, Eqs. (5.6) - (5.9).
// ----------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------

// --------------------------------------------------------------------------------------------------------------------------------
// Further details on ionization provided in:
//
// Ref. [2]: D. Baumann, G. Bertone, J. Stout, and G. M. Tomaselli, "Ionization of gravitational atoms",
// Phys. Rev. D 105 no. 11, (2022) 115036, arXiv:2112.14777 [gr-qc].
//
// Ref. [3]: D. Baumann, G. Bertone, J. Stout, and G. M. Tomaselli, "Sharp Signals of Boson Clouds in Black Hole Binary Inspirals",
// Phys. Rev. Lett. 128 no. 22, (2022) 221102, arXiv:2206.01212 [gr-qc].
// --------------------------------------------------------------------------------------------------------------------------------


// Compile with
// ------------------------------------------------------------------------------------------
// 1. g++ -Wall -I/opt/homebrew/Cellar/gsl/2.7.1/include -c Ionization_eccentric_PUBLIC.cpp
// 2. g++ Ionization_eccentric_PUBLIC.o -L/opt/homebrew/Cellar/gsl/2.7.1/lib -lgsl -lgslcblas
// 3. ./a.out
// ------------------------------------------------------------------------------------------


// ----------------
// Import packages.
// ----------------

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <cmath>
#include <iostream>
#include <float.h>
#include <fstream>
#include <vector>
#include "gsl/gsl_sf_coulomb.h"
#include "gsl/gsl_sf_coupling.h"
#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_roots.h"


// ----------------------------------------------------------
// We normalize by setting the Schwarzschild radius, r_s = 1.
// ----------------------------------------------------------


// -------------------------------------------------------------------------------------------------------------------
// Radial function of the bound states (Eq. (2.7) in [1]). Normalized such that integral_0^\infty Rbound^2 r^2 dr = 1.
// -------------------------------------------------------------------------------------------------------------------

double Rbound (int n, int l, double r, double alpha, double rs)
{
    double rbohr = rs / (2 * pow(alpha,2));
            
    return gsl_sf_hydrogenicR(n,l,1/rbohr,r);
}


// ---------------------------------------------------------------------------------------------------------------
// Radial function of the unbound states (Eq. (2.11) in [1]). Normalized such that Runbound ~ (2/r) * np.sin(rho).
// ---------------------------------------------------------------------------------------------------------------

double Runbound (double k, int l, double r, double alpha, double rs)
{
    double mu = 2*alpha/rs;
    double eta = - mu * alpha/ k;
    double rho = k*r;
    
    gsl_sf_result result, dummy1, dummy2, dummy3;
    double dummy4, dummy5;
    
    gsl_sf_coulomb_wave_FG_e(eta, rho, l, 0, &result, &dummy1, &dummy2, &dummy3, &dummy4, &dummy5);
        
    return result.val * 2 / r;
}


// -------------------------------------------------------
//  Integrand for the radial integral (Eq. (2.18) in [1]).
// -------------------------------------------------------

struct integrandRadial {int lprime; int l_star; int l; double k; int n; double alpha; double rs;};

// ---------------------------------------------------------------------------------------------
// "Inner part" of the potential (Eqs. (2.15)-(2.16) in [1]), excluding the dipole (l_star = 1).
// ---------------------------------------------------------------------------------------------

double integrandInner (double r, void * p)
{
    struct integrandRadial * params = (struct integrandRadial *)p;
        int lprime = (params->lprime);
        int l_star = (params->l_star);
        int l = (params->l);
        double k = (params->k);
        int n = (params->n);
        double alpha = (params->alpha);
        double rs = (params->rs);
    
    return pow(r, l_star+2) * Runbound(k, lprime, r, alpha, rs) * Rbound(n, l, r, alpha, rs);
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
        double k = (params->k);
        int n = (params->n);
        double alpha = (params->alpha);
        double rs = (params->rs);
    
    return Runbound(k, lprime, r, alpha, rs) * Rbound(n, l, r, alpha, rs) / pow(r, l_star-1);
}


// ------------------------------------------------------------------------------------------
// Dipole contribution of the potential (only has an outer part) (Eqs. (2.15)-(2.16) in [1]).
// ------------------------------------------------------------------------------------------

struct integrandRadial_dipole {int lprime; int l_star; int l; double k; int n; double R_star; double alpha; double rs;};

double integrandOuter_dipole (double r, void * p)
{
    struct integrandRadial_dipole * params = (struct integrandRadial_dipole *)p;
        int lprime = (params->lprime);
        int l_star = (params->l_star);
        int l = (params->l);
        double k = (params->k);
        int n = (params->n);
        double R_star = (params->R_star);
        double alpha = (params->alpha);
        double rs = (params->rs);
    
    return Runbound(k, lprime, r, alpha, rs) * Rbound(n, l, r, alpha, rs) * (- pow(r,l_star)/pow(R_star,l_star+1) + pow(R_star,l_star)/pow(r,l_star+1)) * pow(r,2);
}


// -------------------------------------------------------------------------------------------------------------
// Compute radial integral (Eq. (2.18) in [1]), excluding the dipole.
// We use the QAG adaptive integration algorithm in the gsl library.
// For more information, see https://www.gnu.org/software/gsl/doc/html/integration.html#qag-adaptive-integration
// -------------------------------------------------------------------------------------------------------------

double I_r (gsl_integration_workspace * w, int lprime, int l_star, int l, double k, int n, double R_star, double alpha, double rs)
{
    int max_subdivisions = 1000; // Sets the size of the subintervals in which the integration region is divided.
    
    double relativeError = 1e-6;
        
    double innerIntegral, innerError, outerIntegral, outerError;
    integrandRadial parameters = {lprime, l_star, l, k, n, alpha, rs,};
    
    gsl_function Inner;
    Inner.function = &integrandInner;
    Inner.params = &parameters;
    
    // The integration rule is chosen to be the 61 point Gauss-Kronrod rule.
    gsl_integration_qag(&Inner, 0, R_star, 0, relativeError, max_subdivisions, 6, w, &innerIntegral, &innerError);
    
    gsl_function Outer;
    Outer.function = &integrandOuter;
    Outer.params = &parameters;
    
    // Similar integration method as before, yet specialized to semi-infinite intervals (R_star, \infty).
    gsl_integration_qagiu(&Outer, R_star, 0, relativeError, max_subdivisions, w, &outerIntegral, &outerError);
    
    // For large R_star, there is a huge exponential suppression making the integration numerically unfeasible. We set it to zero by hand if necessary.
    if(innerIntegral!=innerIntegral) innerIntegral = 0;
    if(outerIntegral!=outerIntegral) outerIntegral = 0;
            
    return sqrt(rs) * innerIntegral / pow(R_star, l_star+1) + sqrt(rs) * outerIntegral * pow(R_star, l_star);
}


// -------------------------------------------------------------------
// Compute dipole contribution to radial integral (Eq. (2.18) in [1]).
// -------------------------------------------------------------------

double I_r_dipole (gsl_integration_workspace * w, int lprime, int l_star, int l, double k, int n, double R_star, double alpha, double rs)
{
    assert (l_star == 1);
    
    int max_subdivisions = 1000;
    
    double relativeError = 1e-6;
        
    double outerIntegral, outerError;
    integrandRadial_dipole parameters = {lprime, l_star, l, k, n, R_star, alpha, rs,};
    
    gsl_function Outer;
    Outer.function = &integrandOuter_dipole;
    Outer.params = &parameters;
    
    gsl_integration_qagiu(&Outer, R_star, 0, relativeError, max_subdivisions, w, &outerIntegral, &outerError);
    
    if(outerIntegral!=outerIntegral) outerIntegral = 0;
            
    return sqrt(rs) * outerIntegral;
}


// ---------------------------------------------------------------------------
// Compute angular integral (Eq. (2.19) in [1]) in terms of Wigner-3j symbols.
// ---------------------------------------------------------------------------

double I_Omega (int lprime, int l_star, int l, int mprime, int m_star, int m)
{
    return pow(-1,mprime+m_star) * sqrt(((2*lprime+1)*(2*l_star+1)*(2*l+1))/(4*M_PI)) * gsl_sf_coupling_3j(2*lprime, 2*l_star, 2*l, 0, 0, 0) * gsl_sf_coupling_3j(2*lprime, 2*l_star, 2*l, -2*mprime, -2*m_star, 2*m);
}


// -----------------------------------------------
// Relation binary's frequency and orbital radius.
// -----------------------------------------------

double OmegaKepler (double R_star, double q, double rs)
{
    return sqrt(rs/(2*R_star*R_star*R_star));
}


// -------------------------------------
// Kepler's equation (Eq. (5.3) in [1]).
// -------------------------------------

struct KeplerRHS_params {double e; double M;};

double KeplerRHS (double E, void * p)
{
    struct KeplerRHS_params * params = (struct KeplerRHS_params *)p;
        double e = (params->e);
        double M = (params->M);
    
    return E - e * sin(E) - M;
}


// --------------------------------------------------------------------------------------------------------------------------
// Solve for the eccetric anomaly in Kepler's equation: M = E − e sin E. Here, M is in [0,2π].
// Makes use of the gsl Root-Finding library, see https://www.gnu.org/software/gsl/doc/html/roots.html?highlight=gsl_function
// --------------------------------------------------------------------------------------------------------------------------

double EccentricAnomaly(gsl_root_fsolver * s, double M, double e)
{
    gsl_function F;
    KeplerRHS_params parameters = {e, M};
    
    F.function = &KeplerRHS;
    F.params = &parameters;
    
    gsl_root_fsolver_set (s, &F, M-M_PI, M+M_PI);
    
    double epsabs = 1e-8, epsrel = 1e-6;
    
    double a, b, r;
    int iter = 0, max_iter = 100, status;
    
    do {
        iter++;
        status = gsl_root_fsolver_iterate(s);
        r = gsl_root_fsolver_root(s);
        a = gsl_root_fsolver_x_lower(s);
        b = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(a, b, epsabs, epsrel);
    } while (status == GSL_CONTINUE && iter < max_iter);
    
    return r;
}


// --------------------------------------------------------------------
// Integrand for the Fourier coefficients \eta^(g) (Eq. (5.8) in [1]).
// --------------------------------------------------------------------

struct integrand_eta_g_Params {gsl_integration_workspace * w; gsl_root_fsolver * s; double k_g; int lprime; int mprime; int n; int l; int m; double a; double alpha; double q; double e; int g; double rs;};


// --------------------------------------------------------------
// Split integrand into cosine and sine part (Eq. (5.9) in [1]).
// --------------------------------------------------------------

double integrand_eta_g_cos (double M, void * p)
{
    struct integrand_eta_g_Params * params = (struct integrand_eta_g_Params *)p;
        gsl_integration_workspace * w = (params->w);
        gsl_root_fsolver * s = (params->s);
        double k_g = (params->k_g);
        int lprime = (params->lprime);
        int mprime = (params->mprime);
        int n = (params->n);
        int l = (params->l);
        int m = (params->m);
        double a = (params->a);
        double alpha = (params->alpha);
        double e = (params->e);
        double rs = (params->rs);
    
    double I_Omega_result, R_star, phi_star, I_r_result = 0;
    
    double result = 0;
    
    double E = EccentricAnomaly(s,M,e);
            
    for (int l_star = abs(lprime-l); l_star <= l+lprime; l_star++)
    {
        I_Omega_result = I_Omega(lprime,l_star,l,mprime,m-mprime,m);
        
        if (abs(I_Omega_result) > 1e-15) // Check if the angular integral is zero, to speed up the code.
        {
            R_star = a * (1 - e * cos(E));
            if (E < M_PI) phi_star = 2 * atan( sqrt( (1+e)*pow(tan(E/2),2)/(1-e) ) );
            else phi_star = 2 * M_PI - 2 * atan( sqrt( (1+e)*pow(tan(E/2),2)/(1-e) ) );
                        
            if (l_star != 1)
            {
                I_r_result = I_r(w,lprime,l_star,l,k_g,n,R_star,alpha,rs);
            }
            if (l_star == 1)
            {
                I_r_result = I_r_dipole(w,lprime,l_star,l,k_g,n,R_star,alpha,rs);
            }
            
            // See Eqs. (5.1) and (5.8) in [1].
            result += I_Omega_result * gsl_sf_legendre_sphPlm(l_star,abs(m-mprime),0) * cos((m-mprime)*(phi_star-M)) * I_r_result / (2*l_star+1);
        }
    }
        
    return result;
}

double integrand_eta_g_sin (double M, void * p)
{
    struct integrand_eta_g_Params * params = (struct integrand_eta_g_Params *)p;
        gsl_integration_workspace * w = (params->w);
        gsl_root_fsolver * s = (params->s);
        double k_g = (params->k_g);
        int lprime = (params->lprime);
        int mprime = (params->mprime);
        int n = (params->n);
        int l = (params->l);
        int m = (params->m);
        double a = (params->a);
        double alpha = (params->alpha);
        double e = (params->e);
        double rs = (params->rs);
    
    double I_Omega_result, R_star, phi_star, I_r_result = 0;
    
    double result = 0;
    
    double E = EccentricAnomaly(s,M,e);
            
    for (int l_star = abs(lprime-l); l_star <= l+lprime; l_star++)
    {
        I_Omega_result = I_Omega(lprime,l_star,l,mprime,m-mprime,m);
        
        if (abs(I_Omega_result) > 1e-15)
        {
            R_star = a * (1 - e * cos(E));
            if (E < M_PI) phi_star = 2 * atan( sqrt( (1+e)*pow(tan(E/2),2)/(1-e) ) );
            else phi_star = 2 * M_PI - 2 * atan( sqrt( (1+e)*pow(tan(E/2),2)/(1-e) ) );
                        
            if (l_star != 1)
            {
                I_r_result = I_r(w,lprime,l_star,l,k_g,n,R_star,alpha,rs);
            }
            if (l_star == 1)
            {
                I_r_result = I_r_dipole(w,lprime,l_star,l,k_g,n,R_star,alpha,rs);
            }
            
            result += I_Omega_result * gsl_sf_legendre_sphPlm(l_star,abs(m-mprime),0) * sin((m-mprime)*(phi_star-M)) * I_r_result / (2*l_star+1);

        }
    }
        
    return result;
}


// ----------------------------------------------------------------------------------------------------------------------------------------
// Compute the Fourier coefficients \eta^(g) (Eq. (5.8) in [1]).
// We use the QAWO adaptive integration algorithm for oscillatory functions in the gsl library, which makes use of Chebyschev moments.
// For more information, see https://www.gnu.org/software/gsl/doc/html/integration.html#qawo-adaptive-integration-for-oscillatory-functions
// ----------------------------------------------------------------------------------------------------------------------------------------

double eta_g (gsl_integration_workspace * w, gsl_root_fsolver * s, double k_g, int lprime, int mprime, int n, int l, int m, double a, double alpha, double q, double e, int g, double rs)
{
    double relativeError = 1e-5; // For high eccentricities, these error bounds need to be adjusted.
    double absoluteError = 1e-5;
    
    gsl_integration_qawo_table * wf_sin, * wf_cos;
    
    wf_cos = gsl_integration_qawo_table_alloc(m-mprime+g, M_PI, GSL_INTEG_COSINE, 10);
    wf_sin = gsl_integration_qawo_table_alloc(m-mprime+g, M_PI, GSL_INTEG_SINE, 10);
        
    double integral_cos, error_cos, integral_sin, error_sin;
    integrand_eta_g_Params parameters = {w, s, k_g, lprime, mprime, n, l, m, a, alpha, q, e, g, rs,};
    
    gsl_function integrand_cos;
    integrand_cos.function = &integrand_eta_g_cos;
    integrand_cos.params = &parameters;
    
    gsl_integration_qawo(&integrand_cos, 0, absoluteError, relativeError, 1000, w, wf_cos, &integral_cos, &error_cos);
    
    gsl_function integrand_sin;
    integrand_sin.function = &integrand_eta_g_sin;
    integrand_sin.params = &parameters;
    
    gsl_integration_qawo(&integrand_sin, 0, absoluteError, relativeError, 1000, w, wf_sin, &integral_sin, &error_sin);
     
     return 2*(integral_cos-integral_sin);
}

// -----------------------------------------------------------------------------------------
// Main function to compute the ionization rate, power and torque (Eqs. (5.5)-(5.7) in [1]).
// -----------------------------------------------------------------------------------------

struct IonizationRates {double IonizationRate; double IonizationPower; double IonizationAngularMomentumRate_z;};

IonizationRates Rates (gsl_integration_workspace * w, gsl_root_fsolver * s, int n, int l, int m, double a, double alpha, double McOverM, double q, double e, double rs)
{
    double OmegaMean = OmegaKepler(a,q,rs);
    
    double mu = 2*alpha/rs;
    
    double E_n = - mu * alpha*alpha / (2*n*n); // Energy of bound state to leading order (Eq. (2.9) in [1]).
    
    int min_lprime_max = l+6; // Minimum amount of modes that should be summed over.
    
    int max_lprime_max = 40; // Maximum amount of modes that should be summed over.
    
    double max_fractional_increment_lprime = 1e-4;
    
    IonizationRates result;
    
    result.IonizationRate = 0;
    result.IonizationPower = 0;
    result.IonizationAngularMomentumRate_z = 0;
    
    int lprime, mprime;
    
    double k_g;
    
    double tmp1, tmp2;
    
    int how_many_g = 24; // For high eccentricities, needs to be increased.
    
    double etag;
        
    for (lprime = 0; ; lprime++) // We loop over l_prime until one of the conditions* is satisfied.
    {
        std::cout << "lprime = " << lprime << std::endl;
        tmp1 = result.IonizationRate;
        tmp2 = result.IonizationPower;
        
        for (int g = 1; g <= how_many_g; g++)
        {
            if (E_n + g*OmegaMean > 0)
            {
                k_g = sqrt(2*mu*(E_n+g*OmegaMean));
                
                for (mprime = -lprime; mprime <= lprime; mprime++)
                {
                    etag = eta_g(w, s, k_g, lprime, mprime, n, l, m, a, alpha, q, e, g, rs);
                    
                    etag /= 2*M_PI; // Normalization constant
                    
                    std::cout << "lprime = " << lprime << ", g = " << g << ", mprime = " << mprime << ", eta_g = " << etag << std::endl;
                    
                    result.IonizationRate += mu * pow(etag,2) / k_g; // See Eq. (5.5) in [1].
                    result.IonizationPower += g * OmegaMean * mu * pow(etag,2) / k_g; // See Eq. (5.6) in [1].
                    result.IonizationAngularMomentumRate_z += (mprime-m) * mu * pow(etag,2) / k_g; // See Eq. (5.7) in [1].
                }
            }
        }
        
        // ----------------------------------------------------------------------------------------------------------------------------------------
        // Specify the conditions* to break the loop over l_prime, which is either:
        // 1. Enough modes have been summed such that the increment of the rate or power reaches "max_fractional_increment_lprime" (Default: 1e-4).
        // 2. l_prime has reached "max_lprime_max", which is set by hand (Default: 40) (if set too high, code slows down).
        // ----------------------------------------------------------------------------------------------------------------------------------------

        
        if (lprime >= min_lprime_max and tmp1 != 0 and tmp2 != 0)
            if ( fmax(abs(result.IonizationRate/tmp1), abs(result.IonizationPower/tmp2)) < 1+max_fractional_increment_lprime or lprime >= max_lprime_max)
            {
                result.IonizationRate *= pow(4*M_PI*alpha*q,2) / rs; // See Eq. (5.5) in [1].
                result.IonizationPower *= McOverM * pow(4*M_PI*alpha*q,2) / (2*alpha); // See Eq. (5.6) in [1].
                result.IonizationAngularMomentumRate_z *= McOverM * pow(4*M_PI*alpha*q,2) / (2*alpha); // See Eq. (5.7) in [1].
                
                std::cout << "R = " << 2*a << ", lprime = " << lprime << std::endl;
                
                return result;
            }
    }
         
    return result;
}


// -------------------------------------------------
// Specifying parameters and outputting to txt file.
// -------------------------------------------------

int main(int argc, const char * argv[]) {
    
    srand(static_cast<unsigned int>(time(NULL)));
    
    gsl_set_error_handler_off();
        
    int max_subdivisions_K = 1000;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (max_subdivisions_K);
    
    const gsl_root_fsolver_type * T = gsl_root_fsolver_brent;
    gsl_root_fsolver * s = gsl_root_fsolver_alloc(T);
    
    std::ofstream file;
    
    file.open ("211_e_0.txt"); // Path to where you want to save your file.
    
    IonizationRates Rate;
    
    
    // -----------
    // Parameters.
    // -----------
    
    
    // ---------------------------------------------------------------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------------------------

    double P_GW; // Gravitational wave power.
    
    double a; // Semi-major axis.
    
    double q = 0.001; // Mass ratio of the binary (Default: 1e-3).
    
    double alpha = 0.2; // Gravitational fine-structure constant (Default: 0.2).
    
    double e = 0.0; // Eccentricity.
        
    double McOverM = 0.01; // Cloud's mass (Default: 1% of the BH mass).
    
    int n = 2; // Initial (bound) state of the cloud, defined by the principal, angular momentum and azimuthal quantum numbers (Default: nlm = 211).
    
    int l = 1;
    
    int m = 1;
    
    // ---------------------------------------------------------------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------------------------

    
    // -----------------------------------
    // Loop over the desired radial range.
    // -----------------------------------

    for (int i = 10; i < 400; i++)
    {
        a = i*0.5 / (pow(alpha,2)/pow(0.2,2));
        
        P_GW = 16 * pow(q,2) * pow(a, 4) * pow(OmegaKepler(a,q,1), 6) / 5;
        
        P_GW *= (1 + (73*e*e/24) + (37*e*e*e*e/96)) / pow(1-e*e,3.5); // Peter's formula (Eq. 2.13) in [1].
        
        std::cout << "a = " << 2*a << ", e = " << e << std::endl;;
        
        Rate = Rates(w, s, n, l, m, a, alpha, McOverM, q, e, 1);
        
        file << 2*a << " " << Rate.IonizationRate << " " << Rate.IonizationPower << " " << -Rate.IonizationAngularMomentumRate_z << " " << Rate.IonizationPower/P_GW;
        
        file << std::endl;
    }
        
    file.close();
    
    gsl_integration_workspace_free(w);
    
    gsl_root_fsolver_free(s);
        
    return 0;
}
