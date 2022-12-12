===================
READ ME
===================

REPLICATION CODE FOR
"The Risk-Premium Channel of Uncertainty: Implications for Unemployment and Inflation"
Freund, Lee, and Rendahl (2022)


# REPLICATION CODE IS COMPOSED OF FOUR PARTS:

1) MAIN: "mainFLR.m" 
 - MAIN CODE FOR THE REPLICATION
 - SIMPLY RUNNING THIS FILE COMPLETES THE REPLICATION
 - CALLS THE OTHER THREE PARTS BELOW

2) BASELINE MODEL: "mainFLR_baseline.m"
 - CALIBRATED BASELINE MODEL
 - SAVES THE CALIBRATED PARAMETER LEVEL "Parameters.mat"
 - CALLS "dynareFLR.mod"
 - DYNARE OUTPUTS ARE SAVED IN "./Output/IRFs"
 - PRINTS TABLE 2 ON MATLAB COMMAND WINDOW

3) OTHER SPECIFICATIONS: "./mainFLR_others/"
 - FOLDER OF DYNARE FILES WITH OTHER SPECIFICATIONS FOR THE ANALYSIS
 - DYNARE OUTPUTS ARE SAVED IN "./Output/IRFs"

4) FIGURE PLOTTER: "mainFLR_plot.m"
 - LOAD SAVED DYNARE OUTPUTS FROM "./Output/IRFs"
 - CALLS "./mainFLR_others/mainFLR_plot_decomp" FOR DECOMPOSITION ANALYSIS (FIGURE 3)
 - CALLS "./mainFLR_others/mainFLR_print_IRF" FOR FORMATTING FIGURES
 - FIGURES ARE SAVED IN "./Output/Figures"


# SOFTWARE
The code was written using "Matlab R2020b" AND "Dynare 4.4.3" and tested both on Windows and Mac devices.

# RUNTIME
The total runtime is approx. 400 seconds (Windows, Intel 11th gen i7 3Ghz, 32GB RAM).

# CONTACT
For any questions or comments, please contact
Lukas Freund (lukas.beat.freund@gmail.com) & Hanbaek Lee (hanbaeklee1@gmail.com)