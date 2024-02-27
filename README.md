# GrAB (Gravitational Atoms in Binaries)
# Description
This is a repository containing various codes used in <a href="https://arxiv.org/abs/2305.15460">arXiv:2305.15460</a> [1] and <a href="https://arxiv.org/abs/XXXX.XXXXX">arXiv:XXXX.XXXXX</a> [2] to study gravitational atoms in binary systems. 

We include three main folders:\
\
(i) *Capture cross section*, see Section 3 in [1].\
Contains two codes to calculate energy losses to other bound and unbound states when the secondary black hole is on a parabolic orbit.\
\
\
(ii) *Ionization*, see Section 5 and 6 in [1] (also earlier work: <a href="https://arxiv.org/abs/2112.14777">arXiv:2112.14777</a> and <a href="https://arxiv.org/abs/2206.01212">arXiv:2206.01212</a>).\
Contains two codes to calculate energy losses by ionization when the secondary black hole is on eccentric and equatorial or inclined and quasi-circular orbits.\
\
\
(iii) *Resonances*, see [2] (also earlier work: <a href="https://arxiv.org/abs/1804.03208">arXiv:1804.03208</a> and <a href="https://arxiv.org/abs/1912.04932">arXiv:1912.04932</a>)).\
Contains two codes to calculate the overlap integrals when the secondary black hole is on eccentric and equatorial or inclined and quasi-circular orbits. To include energy losses by ionization, we make use of the files in the "Data" folder, which have precomputed values for the 211 and 322 state.\
Additionally, we include a python script that evolves the nonlinear system using dimensionsless variables (see Section 3 of [2]).
# Compiling
All c++ codes can be compiled and run as follows, replacing `XXX.cpp` with the name of the file. The directory of the gsl library should be changed to its location on the user's machine.
<pre><code>g++ -Wall -I/opt/homebrew/Cellar/gsl/2.7.1/include -c XXX.cpp
  
g++ XXX.o -L/opt/homebrew/Cellar/gsl/2.7.1/lib -lgsl -lgslcblas

./a.out
</code></pre>
# Questions
For any questions, feedback or suggestions, please send an email to <a href="mailto:thomas.spieksma@nbi.ku.dk">thomas.spieksma@nbi.ku.dk</a> and <a href="mailto:gimmytomas@gmail.com">gimmytomas@gmail.com</a>.
# Citing
If you make use of this code, please consider citing the corresponding papers,
<pre><code> @article{Tomaselli:2023ysb,
    author = "Tomaselli, Giovanni Maria and Spieksma, Thomas F. M. and Bertone, Gianfranco",
    title = "{Dynamical friction in gravitational atoms}",
    eprint = "2305.15460",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    doi = "10.1088/1475-7516/2023/07/070",
    journal = "JCAP",
    volume = "07",
    pages = "070",
    year = "2023"
}
  
@article{Tomaselli:2024XXX,
    author = "Tomaselli, Giovanni Maria and Spieksma, Thomas F. M. and Bertone, Gianfranco",
    title = "{The resonant history of gravitational atoms in black hole binaries}",
    eprint = "24XX.XXXXX",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    doi = "XXXXXXXXX",
    journal = "XXXX",
    volume = "XX",
    pages = "XXX",
    year = "2024"
}</code></pre>
