# GrAB
# DESCRIPTION
This is a repository containing various codes used in arXiv:2305.15460 [1] and arXiv:XXXX.XXXXX [2] to study gravitational atoms in binary systems. 

We include three main folders:\
\
(i) Capture cross section, based on Section 3 in [1]. Two separate codes to calculate energy lost to other bound states and unbound states.\
\
(ii) Ionization, based on Section 5 and 6 in [1], yet see earlier work (arXiv:2112.14777, arXiv:2206.01212). Two separate codes to calculate the energy losses by ionization on eccentric or inclined orbits.\
\
(iii) Resonances, based on [2]. Two separate codes to calculate the overlap integrals on eccentric or inclined orbits. In addition, a python script to calculate the evolution of the system using dimensionless variables.\
# COMPILING
\
All c++ codes can be compiled as:\
\
g++ -Wall -I/opt/homebrew/Cellar/gsl/2.7.1/include -c XXX.cpp\
\
g++ XXX.o -L/opt/homebrew/Cellar/gsl/2.7.1/lib -lgsl -lgslcblas\
\
./a.out
# QUESTIONS
For any questions, feedback or suggestions, please send an email to <a href="mailto:thomas.spieksma@nbi.ku.dk">thomas.spieksma@nbi.ku.dk</a> and/or <a href="mailto:gimmytomas@gmail.com">gimmytomas@gmail.com</a>
# CITING
