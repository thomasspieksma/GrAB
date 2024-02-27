// --------------------------------------------------------------
// --------------------------------------------------------------
// Compute resonance strength on inclined, quasi-circular orbits.
// --------------------------------------------------------------
// --------------------------------------------------------------


// -----------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------
// Main reference:
//
// Ref. [1]: G. M. Tomaselli, T. F. M. Spieksma, and G. Bertone, "The resonant history of gravitational atoms in black hole binaries".
//
// Specifically, Section 3.
//
// Ref. [2]: G. M. Tomaselli, T. F. M. Spieksma, and G. Bertone, "Dynamical friction in gravitational atoms",
// JCAP 07 (2023) 070, arXiv:2305.15460 [gr-qc]
// -----------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------------
// Further details on resonances provided in:
//
// Ref. [3]: D. Baumann, H. S. Chia, and R. A. Porto, “Probing Ultralight Bosons with Binary Black Holes”,
// Phys. Rev. D 99 no. 4, (2019) 044001, arXiv:1804.03208 [gr-qc].
//
// Ref. [4]: D. Baumann, H. S. Chia, R. A. Porto, and J. Stout, “Gravitational Collider Physics”,
// Phys. Rev. D 101 no. 8, (2020) 083019, arXiv:1912.04932 [gr-qc].
// -------------------------------------------------------------------------------------------------------

// --------------------------------------------------------------------------------------------------------------------------------
// Additional references on ionization:
//
// Ref. [5]: D. Baumann, G. Bertone, J. Stout, and G. M. Tomaselli, "Ionization of gravitational atoms",
// Phys. Rev. D 105 no. 11, (2022) 115036, arXiv:2112.14777 [gr-qc].
//
// Ref. [6]: D. Baumann, G. Bertone, J. Stout, and G. M. Tomaselli, "Sharp Signals of Boson Clouds in Black Hole Binary Inspirals",
// Phys. Rev. Lett. 128 no. 22, (2022) 221102, arXiv:2206.01212 [gr-qc].
// --------------------------------------------------------------------------------------------------------------------------------


// Compile with
// -----------------------------------------------------------------------------------------
// 1. g++ -Wall -I/opt/homebrew/Cellar/gsl/2.7.1/include -c Resonances_inclined_PUBLIC.cpp
// 2. g++ Resonances_inclined_PUBLIC.o -L/opt/homebrew/Cellar/gsl/2.7.1/lib -lgsl -lgslcblas
// 3. ./a.out
// -----------------------------------------------------------------------------------------


// ----------------
// Import packages.
// ----------------


#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <cmath>
#include <iostream>
#include <sstream>
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
// Wigner small d-matrix (Eq. (6.2) in [2]).
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
// Radial function of the bound states (Eq. (2.7) in [2]). Normalized such that integral_0^\infty Rbound^2 r^2 dr = 1.
// -------------------------------------------------------------------------------------------------------------------

double Rbound (int n, int l, double r, double alpha, double rs)
{
    double rbohr = rs / (2 * pow(alpha,2));
            
    return gsl_sf_hydrogenicR(n,l,1/rbohr,r);
}


// --------------------------------------------------------------------------------------------------------------------------------------------
// Energy of bound states of the cloud (Eq. (2.3) in [1]).
// m_set_spin is meant to be the m of the most populated state. It is used to set the BH spin on the threshold of the superradiant instability.
// --------------------------------------------------------------------------------------------------------------------------------------------

double Energy_bound (int n, int l, int m, double alpha, double rs, int m_set_spin)
{
    double mu = 2*alpha/rs;
    
    double fnl = (double) 2/n - (double) 6/(2*l+1);
    double hl = 0;
    if (l > 0) hl = (double) 16/(2*l*(2*l+1)*(2*l+2));
    double atilde = 0.5;
    
    if (2*alpha > m_set_spin) atilde = 1;
        
    return mu * ( - pow(alpha/n,2)/2 - pow(alpha/n,4)/8 + ( fnl*alpha + hl*atilde*m*alpha*alpha) *pow(alpha/n,3) );
}


// -------------------------------------------------------
//  Integrand for the radial integral (Eq. (2.18) in [2]).
// -------------------------------------------------------

struct integrandRadial {int lprime; int l_star; int l; int nprime; int n; double alpha; double rs;};

// -------------------------------------------------------------------------------------------
// "Inner part" of the potential (Eqs. (2.4)-(2.5) in [1]), excluding the dipole (l_star = 1).
// -------------------------------------------------------------------------------------------

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


// -------------------------------------------------------------------------------------------
// "Outer part" of the potential (Eqs. (2.4)-(2.5) in [1]), excluding the dipole (l_star = 1).
// -------------------------------------------------------------------------------------------

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


// ----------------------------------------------------------------------------------------
// Dipole contribution of the potential (only has an outer part) (Eqs. (2.4)-(2.5) in [1]).
// ----------------------------------------------------------------------------------------

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
// Compute radial integral (Eq. (2.18) in [2]), excluding the dipole.
// We use the QAG adaptive integration algorithm in the gsl library.
// For more information, see https://www.gnu.org/software/gsl/doc/html/integration.html#qag-adaptive-integration
// -------------------------------------------------------------------------------------------------------------

double I_r (gsl_integration_workspace * w, int lprime, int l_star, int l, int k, int n, double R_star, double alpha, double rs)
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
// Compute dipole contribution to radial integral (Eq. (2.18) in [2]).
// -------------------------------------------------------------------

double I_r_dipole (gsl_integration_workspace * w, int lprime, int l_star, int l, int k, int n, double R_star, double alpha, double rs)
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
// Compute angular integral (Eq. (2.19) in [2]) in terms of Wigner-3j symbols.
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


// --------------------------------------------------------------
// Main function to calculate coupling strength (Eq. 2.6 in [1]).
// --------------------------------------------------------------

void Resonance (gsl_integration_workspace * w, int nprime, int lprime, int mprime, int n, int l, int m, double alpha, double McOverM, double q, double chi, double rs, bool include_ionization)
{
    
    double DeltaE = Energy_bound(nprime,lprime,mprime,alpha,rs,m) - Energy_bound(n,l,m,alpha,rs,m);
    
    double Omega_resonance, R_resonance;
    
    double I_r_result, I_Omega_result, D;
        
    double eta, gamma, z;
        
    double B, Delta_t_float, exponent;
    
    std::cout << "chi = " << chi << " ";
    
    for (int g = -lprime-l; g <= l+lprime; g++)
    {
        eta = 0;
        I_r_result = 0;
        
        Omega_resonance = DeltaE/g; // See Eq. (3.5) in [1].
                        
        if (Omega_resonance > 0) // Makes sure the g = 0 entries are filtered out.
        {
            std::cout << "g = " << g << " ";

            R_resonance = pow( rs*(1+q)/(2*Omega_resonance*Omega_resonance), 1./3);
                        
            gamma = (96/5) * q * pow(0.5*Omega_resonance, 5./3) * pow(Omega_resonance,2) / pow(1+q, 1./3); // See Eq. (3.6) in [1].
            
            // ---------------------------------------------------------------------------------------------------
            // To include ionization in \gamma, needs two .txt files with saved energy losses.
            
            if (include_ionization == true && R_resonance > 0)
            {
                std::ifstream file;
                file.open("/Data/211_ion.txt.txt");
                
                std::string line1, line2;
                double R1, chi1, chi2, Pion_Over_Pgw11, Pion_Over_Pgw12, Pion_Over_Pgw21, Pion_Over_Pgw;
                
                double R2 = 0.0;

                do {
                    line1 = line2;
                    R1 = R2;
                    
                    std::getline(file,line2);
                    R2 = std::stod(line2) * pow(alpha/0.2,-2) / 2;
                } while (!(R1 <= R_resonance && R2 > R_resonance));
                                
                chi1 = floor(chi*180/M_PI)*M_PI/180;
                chi2 = floor(1+chi*180/M_PI)*M_PI/180;
                
                std::istringstream ss1(line1);
                std::istringstream ss2(line2);
                
                for (int i = 0; i < (3 + 4*floor(chi*180/M_PI)); ++i)
                {
                    ss1 >> Pion_Over_Pgw11;
                    ss2 >> Pion_Over_Pgw12;
                }
                ss1 >> Pion_Over_Pgw21; ss1 >> Pion_Over_Pgw21; ss1 >> Pion_Over_Pgw21; ss1 >> Pion_Over_Pgw21;
                
                Pion_Over_Pgw11 *= pow(alpha/0.2,-5) * McOverM/0.01; // Scaling of ionization, see Eq. (3.31) in [5].
                Pion_Over_Pgw12 *= pow(alpha/0.2,-5) * McOverM/0.01;
                Pion_Over_Pgw21 *= pow(alpha/0.2,-5) * McOverM/0.01;
                                
                Pion_Over_Pgw = Pion_Over_Pgw11 + (Pion_Over_Pgw12-Pion_Over_Pgw11)*(R_resonance-R1)/(R2-R1) + (Pion_Over_Pgw21-Pion_Over_Pgw11)*(chi-chi1)/(chi2-chi1);
                                
                gamma *= 1 + Pion_Over_Pgw;
                
                file.close();
            }
            
            // ---------------------------------------------------------------------------------------------------

            
            B = - 3 * McOverM * pow(Omega_resonance,4./3) * pow(0.5*(1+q),1./3) * g / (q * alpha * sqrt(gamma/abs(g))); // See (3.23) in [1].
            
            Delta_t_float = B / sqrt(gamma*abs(g)); // See (3.24) in [1].
                        
            exponent = gamma * Delta_t_float / Omega_resonance;

            std::cout << "Delta_t_float (yrs) = " << Delta_t_float * 3.12436 * pow(10,-9) << " "; // In years with M = 10^4 M_{\odot}.
                                    
            for (int l_star = std::max(abs(g),abs(lprime-l)); l_star <= l+lprime; l_star++)
            {
                if (l_star != 1)
                {
                    I_Omega_result = I_Omega(lprime,l_star,l,mprime,m-mprime,m);
                    D = Wigner_small_d(l_star, m-mprime, -g, chi);
                    
                    if (l_star == 4 && g == -2) I_Omega_result=0;
                    
                    
                    if (abs(I_Omega_result*D) > 1e-15) // Check if the angular integral is zero, to speed up the code.
                        I_r_result = I_r(w, lprime,l_star,l,nprime,n,R_resonance,alpha,rs);
                                        
                    eta += (4*M_PI * alpha * q) * I_r_result * I_Omega_result * D * gsl_sf_legendre_sphPlm(l_star,abs(g),0) / (2*l_star+1); // See e.g. Eq. (3.8) in [2].
                    
                }
                
                if (l_star == 1)
                {
                    I_Omega_result = I_Omega(lprime,l_star,l,mprime,m-mprime,m);
                    D = Wigner_small_d(l_star, m-mprime, -g, chi);
                    
                    I_Omega_result=0;
                                        
                    if (abs(I_Omega_result*D) > 1e-15)
                    {
                        I_r_result = I_r_dipole(w, lprime,l_star,l,nprime,n,R_resonance,alpha,rs);
                    }
                                        
                    eta += (4*M_PI * alpha * q) * I_r_result * I_Omega_result * D * gsl_sf_legendre_sphPlm(l_star,abs(g),0) / (2*l_star+1);
                }
            }
            
            z = eta*eta/(gamma*abs(g));
            
            std::cout << "B=" << B << " ";
                        
            std::cout << "2Pi*z*B=" << 2*M_PI*z*B << " ";
            
            std::cout << "D=" << exponent << " ";
            
            std::cout << " " << " " << " " << " ";
        }
        
    }
    
    std::cout << std::endl;
}


// -------------------------------------------------
// Specifying parameters and outputting to txt file.
// -------------------------------------------------

int main(int argc, const char * argv[]) {
    
    gsl_set_error_handler_off();
        
    int max_subdivisions_K = 1000;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (max_subdivisions_K);
        
    std::ofstream file;
    file.open ("Testdata.txt");
    
    // -----------
    // Parameters.
    // -----------
    
    
    // ---------------------------------------------------------------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------------------------
                    
        double q = 0.001; // Mass ratio of the binary (Default: 1e-3).
        
        double alpha = 0.2; // Gravitational fine-structure constant (Default: 0.2).
                        
        double McOverM = 0.01; // Cloud's mass (Default: 1% of the BH mass).
        
        int n = 2; // Initial (bound) state of the cloud, defined by the principal, angular momentum and azimuthal quantum numbers (Default: nlm = 211).

        int l = 1;

        int m = 1;
    
        int nprime = 2; // Second bound state that participates in the resonance.
    
        int lprime = 1;
    
        int mprime = -1;
    
        bool include_ionization = false; // Including ionization energy losses or not. NOTE: hyperfine and fine resonances happen at large radii at which ionization energy losses essentially become zero. Therefore, we have not tabulated ionization energy losses at very large radii. One should thus use "false" when computing hyperfine or fine resonances and only use "true" for Bohr resonances.
    
    // ---------------------------------------------------------------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------------------------
    
    // ------------------------------------
    // Loop over the desired angular range.
    // -------------------------------------
    
    for (int chi = 0; chi <= 10; chi++)
    {
        Resonance(w,nprime,lprime,mprime,n,l,m,alpha,McOverM,q,chi*M_PI/180,1,include_ionization);
        file << chi*M_PI/180 << " " << mprime;
        
        file << std::endl;

    }
    
    file.close();

    return 0;
}
