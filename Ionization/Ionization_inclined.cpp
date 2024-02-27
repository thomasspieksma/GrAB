// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// Compute ionization rate, power and torque on inclined, quasi-circular orbits.
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------


// ----------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------
// Main reference:
//
// Ref. [1]: G. M. Tomaselli, T. F. M. Spieksma, and G. Bertone, "Dynamical friction in gravitational atoms",
// JCAP 07 (2023) 070, arXiv:2305.15460 [gr-qc]
//
// Specifically, Section 6, Eqs. (6.5) - (6.7) and Figure 8 for the conventions.
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
// 1. g++ -Wall -I/opt/homebrew/Cellar/gsl/2.7.1/include -c Ionization_inclined_PUBLIC.cpp
// 2. g++ Ionization_inclined_PUBLIC.o -L/opt/homebrew/Cellar/gsl/2.7.1/lib -lgsl -lgslcblas
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
#include "gsl/gsl_sf_coulomb.h"
#include "gsl/gsl_sf_coupling.h"
#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_sf_legendre.h"


// ----------------------------------------------------------
// We normalize by setting the Schwarzschild radius, r_s = 1.
// ----------------------------------------------------------


inline long double factorial (int n)
{
    if (n == 0) return 1;
    else return n * factorial(n-1);
}


// -----------------------------------------
// Wigner small d-matrix (Eq. (6.2) in [1]).
// -----------------------------------------

double Wigner_small_d (int j, int mprime, int m, double beta)
{
    long double x = 0;
    int smin = std::max(0, m-mprime);
    int smax = std::min(j+m, j-mprime);
    
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


// -----------------------------------------------------------------------------------------
// Main function to compute the ionization rate, power and torque (Eqs. (6.5)-(6.7) in [1]).
// -----------------------------------------------------------------------------------------

struct IonizationRates {double IonizationRate; double IonizationPower; double IonizationAngularMomentumRate_z; double IonizationAngularMomentumRate_x;};

IonizationRates Rates (gsl_integration_workspace * w, int n, int l, int m, double R_star, double alpha, double McOverM, double q, double beta, double rs)
{
    double Omega = OmegaKepler(R_star,q,rs);
    
    double mu = 2*alpha/rs;
    
    double E_n = - mu * alpha*alpha / (2*n*n); // Energy of bound state to leading order (Eq. (2.9) in [1]).
    
    int min_lprime_max = l+6; // Minimum amount of modes that should be summed over.
    
    int max_lprime_max = 40; // Maximum amount of modes that should be summed over.
    
    double max_fractional_increment_lprime = 1e-4;
    
    IonizationRates result;

    result.IonizationRate = 0;
    result.IonizationPower = 0;
    result.IonizationAngularMomentumRate_z = 0;
    result.IonizationAngularMomentumRate_x = 0;
    
    int lprime, mprime;
    
    int mtilde;
    
    double d_matrix[2*l+1];
    
    for (mtilde = -l; mtilde <= l; mtilde++)
        d_matrix[mtilde+l] = Wigner_small_d(l, m, mtilde, beta);
    
    double I_Omega_result, I_r_result, k_g, L_x, L_z;
    
    double eta[2*l+1];
    
    double tmp1, tmp2;
    
    int how_many_g = 1 + max_lprime_max + l;
    
    int how_many_l_star = 1 + max_lprime_max + l;
    
    double stored_values_I_r[how_many_g][how_many_l_star];
    
    double dumb_I_r_value = 0.123456789; // Dummy variable
    
    for (lprime = 0; ; lprime++) // We loop over l_prime until one of the conditions* is satisfied.
    {
        tmp1 = result.IonizationRate;
        tmp2 = result.IonizationPower;
        
        for (int g = -floor(E_n/Omega); g <= l+lprime; g++)
        {
            
            k_g = sqrt(2*mu*(E_n+g*Omega));
            
            for (int gg = 0; gg < 1 + l + lprime; gg++)
                for (int ll_star = 0; ll_star < 1 + lprime + l; ll_star++)
                    stored_values_I_r[gg][ll_star] = dumb_I_r_value; // Initially, "stored_values_I_r" has the dummy variable "dumb_I_r_value" in all its entries.
            
            for (mtilde = -l; mtilde <= l; mtilde++)
                eta[mtilde+l] = 0;
                                    
            for (mtilde = std::max(-l, -lprime - g); mtilde <= std::min(l, lprime - g); mtilde++)
            {
                mprime = g + mtilde;
                
                for (int l_star = std::max(g,abs(lprime-l)); l_star <= l+lprime; l_star++)
                {
                    if (l_star != 1)
                    {
                        I_Omega_result = I_Omega(lprime,l_star,l,mprime,-g,mtilde);
                        
                        if (abs(I_Omega_result) > 1e-15) // Check if the angular integral is zero, to speed up the code.
                        {
                            
                            if (abs(stored_values_I_r[g][l_star] - dumb_I_r_value) < 1e-15) // Only compute values for I_r that haven't been computed before, to speed up the code.
                            {
                                I_r_result = I_r(w,lprime,l_star,l,k_g,n,R_star,alpha,rs);
                                stored_values_I_r[g][l_star] = I_r_result;
                            }
                            else
                                I_r_result = stored_values_I_r[g][l_star]; // Otherwise, used the value that was already stored (no need to recompute it).
                            
                            eta[mtilde+l] += I_r_result * I_Omega_result * gsl_sf_legendre_sphPlm(l_star,g,0) / (2*l_star+1); // The same form as Eq. (5.1) in [1].
                        }
                    }
                    if (l_star == 1)
                    {
                        I_Omega_result = I_Omega(lprime,l_star,l,mprime,-g,mtilde);
                        
                        if (abs(I_Omega_result) > 1e-15)
                        {
                            if (abs(stored_values_I_r[g][l_star] - dumb_I_r_value) < 1e-15)
                            {
                                I_r_result = I_r_dipole(w,lprime,l_star,l,k_g,n,R_star,alpha,rs);
                                stored_values_I_r[g][l_star] = I_r_result;
                            }
                            else
                                I_r_result = stored_values_I_r[g][l_star];
                            
                            eta[mtilde+l] += I_r_result * I_Omega_result * gsl_sf_legendre_sphPlm(l_star,g,0) / (2*l_star+1);
                        }
                    }
                }
                
                                                                
                result.IonizationRate += mu * pow(eta[mtilde+l]*d_matrix[mtilde+l],2) / k_g; // See Eq. (6.5) in [1].
                result.IonizationPower += g * Omega * mu * pow(eta[mtilde+l]*d_matrix[mtilde+l],2) / k_g; // See Eq. (6.6) in [1].
                result.IonizationAngularMomentumRate_z += g * mu * pow(eta[mtilde+l]*d_matrix[mtilde+l],2) / k_g; // See Eq. (6.12) in [1].
                
                if (mtilde > -l)
                {
                    result.IonizationAngularMomentumRate_x += sqrt(lprime*(lprime+1)-mprime*(mprime-1)) * mu * eta[mtilde+l-1]*eta[mtilde+l] * d_matrix[mtilde+l-1]*d_matrix[mtilde+l]/ k_g; // See Eq. (6.17) in [1].
                    result.IonizationAngularMomentumRate_x -= sqrt(l*(l+1)-mtilde*(mtilde-1)) * d_matrix[mtilde+l-1]*d_matrix[mtilde+l] * mu * (pow(eta[mtilde+l],2) + pow(eta[mtilde+l-1],2)) / (2*k_g) ; // See Eq. (6.12) in [1].
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
                result.IonizationRate *= pow(4*M_PI * alpha * q, 2) / rs; // See Eq. (6.5) in [1].
                result.IonizationPower *= McOverM * pow(4*M_PI * alpha * q, 2) / (2*alpha); // See Eq. (6.6) in [1].
                result.IonizationAngularMomentumRate_z *= McOverM * pow(4*M_PI * alpha * q, 2) / (2*alpha); // See Eq. (6.7) in [1].
                result.IonizationAngularMomentumRate_x *= McOverM * pow(4*M_PI * alpha * q, 2) / (2*alpha); // See Eq. (6.12) in [1].
                
                L_z = result.IonizationAngularMomentumRate_z;
                L_x = result.IonizationAngularMomentumRate_x;
                
                result.IonizationAngularMomentumRate_z = cos(beta) * L_z - sin(beta) * L_x; // See Eq. (6.21) in [1].
                result.IonizationAngularMomentumRate_x = cos(beta) * L_x + sin(beta) * L_z; // See Eq. (6.21) in [1].
                
                std::cout << "R = " << 2*R_star << ", lprime = " << lprime << ", beta = " << beta << std::endl;
                
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
    
    std::ofstream file;
    
    file.open ("211_inclined.txt"); // Path to where you want to save your file.
    
    IonizationRates Rate;
    
    // -----------
    // Parameters.
    // -----------
    
    
    // ---------------------------------------------------------------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------------------------
    
    double P_GW; // Gravitational wave power.
    
    double R; // Binary separation.
    
    double q = 0.001; // Mass ratio of the binary (Default: 1e-3).
    
    double alpha = 0.2; // Gravitational fine-structure constant (Default: 0.2).
    
    double McOverM = 0.01; // Cloud's mass (Default: 1% of the BH mass).
    
    int n = 2; // Initial (bound) state of the cloud, defined by the principal, angular momentum and azimuthal quantum numbers (Default: nlm = 211).
    
    int l = 1;
    
    int m = 1;
    
    // ---------------------------------------------------------------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------------------------
            
    
    // -----------------------------------
    // Loop over the desired radial range.
    // -----------------------------------
    
    for (int i = 10; i < 4000; i++)
    {
        R = i*0.05 / (pow(alpha,2)/pow(0.2,2));
        
        P_GW = 16 * pow(q,2) * pow(R, 4) * pow(OmegaKepler(R,q,1), 6) / 5;
        
        // --------------------------------------------------
        // Loop over desired angular range, beta is in [0,π].
        // --------------------------------------------------
        
        for (int j = 0; j < 181; j++)
        {
            Rate = Rates(w, n, l, m, R, alpha, McOverM, q, j*M_PI/180, 1);
            file << 2*R << " " << j*M_PI/180 <<" " << Rate.IonizationRate << " " << Rate.IonizationPower << " " << Rate.IonizationAngularMomentumRate_z << " " << Rate.IonizationAngularMomentumRate_x << " " << Rate.IonizationPower/P_GW;
            
            // --------------------------------------------------------------------------------------------------------------------------
            // Output format: for each R, it prints beta (looped from [0,π]), IonRate, IonPower, Torque_z, Torque_x and IonPower/GWPower.
            // So each R will have j times 6 entries.
            // --------------------------------------------------------------------------------------------------------------------------

        }
        
        std::cout << std::endl;
        
        file << std::endl;
    }
    
    file.close();
    
    gsl_integration_workspace_free(w);
        
    return 0;
}
