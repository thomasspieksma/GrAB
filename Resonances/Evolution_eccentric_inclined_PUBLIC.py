# ---------------------------------------------
# ---------------------------------------------
# Evolve cloud-binary system with backreaction.
# ---------------------------------------------
# ---------------------------------------------


# -----------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------
# Main reference:
#
# Ref. [1]: G. M. Tomaselli, T. F. M. Spieksma, and G. Bertone, "The resonant history of gravitational atoms in black hole binaries".
#
# Specifically, TBA
# -----------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------


# ----------------
# Import packages.
# ----------------

import math
import random

import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
import scipy.special as sc
from scipy import integrate
from scipy.linalg import cosm, expm, sinm, inv
from sympy.physics.quantum.cg import Wigner3j
from sympy.physics.quantum.spin import Rotation


# -----------------------------------------------
# Dimensionless parameters, see Eq. (3.8) in [1].

Z = 0.01 # See Eq. (3.8) in [1].

B = 1000 # See Eq. (3.23) in [1].

C = 3*B * 0.25 # See Eq. (3.23) in [1].

Deltam_g = 1 # See Fig. 6 in [1].

Gamma = 0 # Equals 1/\tau_decay, see Fig. 9 in [1].

init_epsilon_squared = 0.3 # Initial eccentricity.

init_beta = 1.3 # Initial inclination angle.

init_omega = 600

omega = np.array([-init_omega*np.sqrt(Z)])

beta = np.array([init_beta])

epsilon_squared = np.array([init_epsilon_squared])

# -----------------------------------------------


# ----------------------------------------------------------------------------------
# Occupation densities

d1 = omega[0]/2 - np.sqrt((omega[0]/2)**2+np.sqrt(Z)**2) # Initially occupied state.
d2 = np.sqrt(Z) # Initially deoccupied state.
norm = np.sqrt(d1**2 + d2**2)
d1 /= norm
d2 /= norm
d = np.array([[d1, d2]])

# ----------------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# Function that will enable adaptive timestepping.

def check (x_old, x_new):
    threshold = 0.2
    threshold_2 = 1e-1
    if abs(x_new-x_old) < threshold * (abs(x_old) + abs(x_new)) + threshold_2:
        return True
    else:
        return False
        
# ----------------------------------------------------------------------------


t = np.array([0])

dt = 0.01 # Initial timestep.

# Peter's formula, see Eq. (3.15) in [1].
def f(x):
    return (1+73*x/24+37*x**2/96) / (1-x)**(7/2)

# See Eq. (3.19) in [1].
def h(x):
    return (1+7*x/8) / (1-x)**(2)

# -----------------------------------------------
# Useful fuctions that will be needed down below.

def H(x):
    return (1-x)**(1/2) * f(x) - h(x)

def k(x):
    return (1-x)**0.5-(1-x)

# -----------------------------------------------



# ---------------------------------------------------------------------------------------------------
# We use a combination of Crank-Nicoloson (for the occupations) and Heun methods for the integration.
# Therefore, it is accurate to second order and unitary.
# ---------------------------------------------------------------------------------------------------

# We continue the evolution until some predefined value for omega.
while (omega[-1] < 1*init_omega*np.sqrt(Z)):

    
    print("Completion: ", "{:.2f}".format(0.5 + omega[-1] / (2*init_omega*np.sqrt(Z))), ", |c_2|^2 =", "{:.2f}".format(abs(d[-1][1])**2))

    
    B_repl = B * 2 * np.sqrt(Z) # An additional factor of 2 \sqrt(Z) comes from rewriting B d|c_b|^{2}/dt using the Schrödinger eqn, as done in the lines below.
    
    domega = (f(epsilon_squared[-1]) - B_repl * (1j*np.conj(d[-1][0])*d[-1][1]).real) * dt/2 # See Eq. (3.20) in [1].
    
    depsilon_squared = - 2 * (1-epsilon_squared[-1])**(1/2) * ( H(epsilon_squared[-1]) + ( Deltam_g * np.cos(beta[-1]) - (1-epsilon_squared[-1])**(1/2) ) * B_repl * (1j*np.conj(d[-1][0])*d[-1][1]).real ) * (dt/2) / C # See Eq. (3.21) in [1].
    
    dbeta = - B_repl * Deltam_g * (1j*np.conj(d[-1][0])*d[-1][1]).real * np.sin(beta[-1]) * (dt/2) / (C * (1-epsilon_squared[-1])**(1/2)) # See Eq. (3.22) in [1].
    
    
    # Heun method. First calulcate intermediate value.
    omega_mid = omega[-1] + domega
    
    epsilon_squared_mid = epsilon_squared[-1] + depsilon_squared
    
    beta_mid = beta[-1] + dbeta


    # Crank-Nicolson method for the occupation densities.
    H_dressed = np.array([[omega_mid/2, np.sqrt(Z)],[np.sqrt(Z), -omega_mid/2 - Gamma*1j]])
    
    matrix1 = np.identity(2) - 1j * H_dressed * dt/2
    
    matrix2 = np.identity(2) + 1j * H_dressed * dt/2
    
    
    d_new = matrix1.dot(d[-1])
    
    d_new = inv(matrix2).dot(d_new)


    # Heun method. Use the intermediate value to calculate the new values.
    omega_new = omega[-1] + domega + (f(epsilon_squared[-1]+2*depsilon_squared) - B_repl * (1j*np.conj(d_new[0])*d_new[1]).real) * dt/2
    
    epsilon_squared_new = epsilon_squared[-1] + depsilon_squared - 2 * (1-epsilon_squared[-1]-2*depsilon_squared)**(1/2) * ( H(epsilon_squared[-1]+2*depsilon_squared) + ( Deltam_g * np.cos(beta[-1]+2*dbeta) - (1-epsilon_squared[-1]-2*depsilon_squared)**(1/2) ) * B_repl * (1j*np.conj(d_new[0])*d_new[1]).real ) * (dt/2) / C
    
    beta_new = beta[-1] + dbeta - B_repl * Deltam_g * (1j*np.conj(d_new[0])*d_new[1]).real * np.sin(beta[-1]+2*dbeta) * (dt/2) / (C * (1-epsilon_squared[-1]-2*depsilon_squared)**(1/2))
    
    
    # ------------------------------------------------------------------------------------------------
    # Adaptive time steps. If the change in occupation densities, frequency and eccentricity is small,
    # we increase the step size (by 0.1%). If it's big, we decrease the step size by 2.
    # ------------------------------------------------------------------------------------------------
    
    if (check(d[-1][0]+d[-1][1], d_new[0]+d_new[1]) and check(omega[-1],omega_new) and check(epsilon_squared_new,epsilon_squared[-1])):
        t = np.append(t, t[-1]+dt)
        omega = np.append(omega, omega_new)
        epsilon_squared = np.append(epsilon_squared, epsilon_squared_new)
        beta = np.append(beta, beta_new)
        d = np.append(d,[d_new], axis=0)
        dt *= 1.001
    else:
        dt /= 2
        

# ------------------
# Write to txt file.
# ------------------

file = open("data.txt", "w")
for i in range(len(omega)):
    if (i%5==0):
        file.write(str(t[i])+" "+str(omega[i])+" "+str(epsilon_squared[i])+" "+str(beta[i])+" "+str(abs(d[i][0])**2)+" "+str(abs(d[i][1])**2)+"\n")
file.close()


# -------------------------
# Plot relevant quantities.
# -------------------------

fig,ax=plt.subplots(1,1, figsize = (12,8))

ax.plot(t,abs(d[:,0])**2, label = r'$|c_{a}|^{2}$')
ax.plot(t,abs(d[:,1])**2, label = r'$|c_{b}|^{2}$')
ax.plot(t,-omega/(2*omega[0])+0.5,label=r'$\omega$')
ax.plot(t, omega * np.exp(-2*np.pi*np.sqrt(Z)**2) / omega, label = 'LZ prediction')
ax.plot(t, epsilon_squared**0.5,label=r'$\varepsilon$')
ax.plot(t, beta,label=r'$\beta$')
ax.set_xlabel(r'$\tau$', fontsize = 20)
plt.legend(fontsize = 20)
plt.show()
