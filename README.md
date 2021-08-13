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

    Using methyl radical as an example, shown above is a typical Gaussian input file for SNO analysis. There would be an extra output file named CH3.47 in addition to ordinary Gaussian output files.

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

    Corresponding orbital transformation (COT) is also implemented and can be called with a parameter 'COT' or simply 'c'. Without any extra parameter, the program computes SNO by default.


