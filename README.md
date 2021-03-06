# SNO
**Spin Natural Orbital**
This is a script to compute SNO based on Gaussian calculation results. Gaussian itself has a keyword to do the same thing which however is not very convenient.

## Requirement
- Python 3
- Numpy 1.18.5
- Gaussian 09 or 16

## Tutorial
1. Run Gaussian calculation and generate 47 file

    **Input**

    - *CH3.gjf*
    ```
    %chk=CH3.chk
    #p opt freq b3lyp/6-31G* pop=nboread

    CH3

    0 2
     C                 -0.00000000    0.00000000    0.00000000
     H                 -0.00000000   -1.00880567   -0.35666667
     H                 -0.87365134    0.50440284   -0.35666667
     H                  0.87365134    0.50440284   -0.35666667

     $nbo archive FILENAME=CH3 $end
     ```

    **Output**
    - *CH3.chk*
    - *CH3.47*

    Using methyl radical as an example, shown above are a typical Gaussian input file and the corresponding output files needed for SNO analysis. In addition to ordinary Gaussian output files, there is an extra CH3.47 file that contains necessary info for NBO analysis. However, NBO analysis is not needed for SNO analysis. The 47 file serves only as an interface for SNO analysis.

2. Run SNO analysis
    First, convert chk file to fchk file for later visualization.
    ```
    formchk CH3.chk CH3.fchk
    ```
    Then, run the following command to perform SNO analysis
    ```
    python SNO.py CH3.47 [CH3.fchk] [SNO|s|COT|c]
    ```
    If fchk file is supplied, the SNOs will be saved as a fchk file named `CH3_SNO.fchk`. Otherwise, only SNO eigenvalues will be printed.

    When saving SNO to fchk file, a warning message `Warning: Spin symmetry not match!` might be printed which can be safely ignored. This warning arises only to remind that MOs are spin-polarized (alpha MOs are not idential to beta MOs) but SNOs are spinless.

    Corresponding orbital transformation (COT) is also implemented and can be called with a parameter 'COT' or simply 'c'. Without any extra parameter, the program computes SNO by default.

3. Analyze SNO results
    If you open the SNO fchk file in GaussView, please remember that the occupations (shown as up/down arrows in the orbital visualization window) are NOT meaningful any more. SNOs are not ordered in orbital energy (as ordinary MOs do). Instead, they are ordered by absolute SNO eigenvalues (i.e. how much spin each orbital contributes). Hence the 1st orbital (with largest absolute eigenvalue) is usually what you should first take a look.

## Related publication
For an extensive discussion on the chemical implication of SNO and comparison against COT, see doi.org/10.1002/jcc.25762


