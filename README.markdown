BMAGWA software v1.0 README
===========================

Introduction
------------

This software implements a method for computing posterior association 
probabilities of SNPs (and other quantities) in genome-wide association studies 
using Bayesian variable selection and model averaging.


Availability
------------

Binary versions for 64bit Linux and Windows and source code can be downloaded at
the [homepage][] of the software or from their [GitHub repositories][]. The
program may potentially run faster when compiled from source, but installing the
required libraries may require some work (although often they may already be
available if your computer is set up for development for numerical computing).

The software has been developed (and tested) on 64bit Linux.

  [homepage]:            http://www.lce.hut.fi/research/mm/bmagwa/
  [GitHub repositories]: https://github.com/to-mi


Usage
-----

The software is used from command line. On Linux use terminal client and on 
Windows run `cmd` to access the command line. Then, change to the directory of 
the binary and execute it with the name of your configuration file as an 
attribute:

    cd /directory/of/the/binary/
    ./bmagwa config.ini


Configuration file
------------------

~~~~~ ini

; Example configuration file
; (lines starting with ; are comments)

[datafiles]
; data files are in PLINK format (use only standard formats;
; not long or transposed etc.)
; (see http://pngu.mgh.harvard.edu/~purcell/plink/)
;
; fam contains information on individuals (only phenotype is used)
file_fam = testdata/testdata.fam
; genotypes in binary format (SNP-major mode; this is default in PLINK)
file_g = testdata/testdata.bed
; covariates (excluding constant, which is added automatically)
file_e = testdata/testdata.e
; optional phenotype file (overruns phenotypes of fam if provided)
file_y = testdata/testdata.y
; set following to 1 to force coding of genotypes to 0,1,2 according to the
; number of minor alleles (i.e. if computing allele freq with current coding
; gives > 0.5, 0s are swapped to 2s and 2s to 0s).
; Usually you probably use 0.
recode_g_to_minor_allele_count = 1

[sizes]
; number of individuals
n = 1000
; number of SNPs
m_g = 10000
; number of covariates (excluding constant)
m_e = 2

[sampler]
; number of iterations
do_n_iter = 1000000
; skip (in iterations) between computing rao-blackwellization (RB)
n_rao = 500
; number of RB-steps to discard as burnin (affects _rao.dat outputfiles)
; relative to n_rao (so these settings start computing RB posterior
; estimates after 500*1000 iterations)
n_rao_burnin = 1000
; keep at 0 (experimental)
adaptation = 0
; proposal distribution flattening parameter
t_proposal = 0.5
; skip between printing status to _verbose.dat file
verbosity = 10000
; skip between saving MCMC samples
thin = 100
; set 1 to save samples of regression coefficients
save_beta = 0
; skip between updating tau2 and missing genotypes
n_sample_tau2_and_missing = 10

[thread]
; number of chains to run simultanously
n_threads = 2
; prefix for result files (directory should exist; in this case
; 'testdata/results/')
basename = testdata/results/chain
; seed for random number generator for each chain (comma separated)
seeds = 1234,2345

[model]
; allowed effect types (possible values are A,H,D,R,AH; comma separated)
types = A,AH

[prior]
; prior weights on different effect types (t parameter prior)
type_A = 1
type_H = 1
type_D = 1
type_R = 1
type_AH = 1
; model size prior (mean and variance) used to set prior for gamma parameter
e_qg = 20
var_qg = 400
; residual variance (sigma2 parameter) prior
; give either R2mode_sigma2 (expected mode for proportion of variance
; explained) or s2_sigma2 (directly the parameter of invchi2-prior)
R2mode_sigma2 = 0.40
; s2_sigma2 = 1.0
nu_sigma2 = 1
; effect size (tau2_t parameters) priors
; give either eh_tau2_X (expected heritability contribution) or s2_tau2_X
nu_tau2_A = 3
eh_tau2_A = 0.02
nu_tau2_H = 3
eh_tau2_H = 0.02
nu_tau2_D = 3
eh_tau2_D = 0.02
nu_tau2_R = 3
eh_tau2_R = 0.02
; contant term effect size prior parameter (1/sigma^2_alpha_0)
inv_tau2_e_const_val = 0
; covariate term effect size prior parameter (1/sigma^2_alpha_l)
inv_tau2_e_val = 1

~~~~~


Output files and formats
------------------------

Most of the output files are in binary format. A Python script/module to convert
the results to text file format is included (run `python postprocess.py` for 
help).

The output file names are prefixed with the `basename`-setting and the number of
the chain.

Note: thinning affects the output. If `do_n_iter = 100` and `thin = 10` only 10 
numbers will be in the (MCMC trace) output files.

### Text files

 * `basename_prior.dat`

   Contains the values of prior parameters.

 * `basenameX_verbose.dat`

   The state of the sampler is printed to this file every verbosity iterations.
   You can check this file to monitor the progress.


### Binary files

Format of a single value is given after filename.

 * `basenameX_accepted.dat [char ('b')]`

   Contains 1 if move was accepted and 0 if not for each MCMC iteration.

 * `basenameX_loci.dat [unsigned int32 ('I')]`

   Contains the SNPs included in the model for each MCMC iteration. You need
   `basenameX_modelsize.dat` to parse this file.

 * `basenameX_log_likelihood.dat [double ('d')]`

   Contains the log likelihood of the linear model for each MCMC iteration.

 * `basenameX_log_prior.dat [double ('d')]`

   Contains the log of model prior (model size and effect type prior) for each
   MCMC iteration.

 * `basenameX_modelsize.dat [unsigned int32 ('I')]`

   Contains the number of SNPs in the model for each MCMC iteration.

 * `basenameX_move_type.dat [char ('b')]`

   Contains the type of move proposed for each MCMC iteration.

 * `basenameX_pve.dat [double ('d')]`

   Contains the total proportion of variance explained and the proportions of
   variance explained by genetic effects and covariates for each MCMC iteration.

 * `basenameX_rao.dat [double ('d')]`

   Contains the Rao-Blackwellised posterior association probability for each 
   SNP.

 * `basenameX_rao_types.dat [double ('d')]`

   Contains the Rao-Blackwellised posterior probabilities of effect types for
   each SNP (conditional on this SNP being included in the model). Note the
   ordering of types. Not created if only one possible type is specified in
   types.

 * `basenameX_sigma2.dat [double ('d')]`

   Contains the residual variance for each MCMC iteration.

 * `basenameX_singlesnp_post.dat [double ('d')]`

   Contains the posterior association probability for each SNP computed with
   only that SNP in the model (and using only one value of tau2).

 * `basenameX_singlesnp_post_types.dat [double ('d')]`

   Contains the posterior probabilities of effect types for each SNP conditional
   on only that SNP being in the model (again using only one value of tau2).
   Note the ordering of types. Not created if only one possible type is
   specified in types.

 * `basenameX_tau2.dat [double ('d')]`

   Contains the value of tau2 for each MCMC iteration.

 * `basenameX_types.dat [char ('b')]`

   Contains the effect type of each SNP included in the model for each MCMC 
   iteration. Use `basenameX_loci.dat` to find the corresponding SNPs. Note the
   ordering of types. Not created if only one possible type is specified in 
   types.

Ordering of types: Effect types are in order (with corresponding numeric value 
in `basenameX_types.dat`) A=0, H=1, D=2, R=3, AH=4 (excluding those not defined 
in types; however, the order in types is insignificant).


Example
-------

A simple dataset generated with PLINK is included in 'testdata' directory. It 
has 10000 SNPs and 1000 individuals. The last 20 SNPs are causal with 
heritability contribution of 0.02 each. To run the BMAGWA software on the 
dataset, on the command line run:

    cd /directory/of/the/binary/
    ./bmagwa testdata/testdata.ini

Running the sampler for the dataset took under 10 minutes on 2.8GHz Intel Core2 
Quad CPU. You can monitor the progress by looking into 
`testdata/results/chain0_verbose.dat` (you can tweak the sampling options in the
configuration file if it seems to take too long...).

Print out from the program should look something like this:

    -------------------------------------------------------------
    BMAGWA software version 1.0
    http://www.lce.hut.fi/research/mm/bmagwa/
    -------------------------------------------------------------
    
    Warning: using default value of 0 for option prior.s2_sigma2
    Warning: using default value of 0 for option prior.s2_tau2_A
    Warning: using default value of 0 for option prior.s2_tau2_H
    var y = 0.990225
    var x = 0.353579
    mean x = 0.530252
    Precomputing SNP covariances
    Initializing sampler 0
    Initializing sampler 1
    Creating thread 0
    Creating thread 1
    Completed join with thread 0 having a status of 0
    Completed join with thread 1 having a status of 0

You may then use the provided Python script to convert the binary files to text 
files and to compute posterior association probabilities from the MCMC chain. On
command line:

    python postprocess.py convert testdata/results
    python postprocess.py mcmcpos testdata/results/chain 10000 5000 1

The .txt files are probably easy to open in your preferred numerical software 
for further analysis. For example, in Matlab:

~~~~~ matlab

p_rao1 = load('testdata/results/chain0_rao.dat.txt');
p_rao2 = load('testdata/results/chain1_rao.dat.txt');
p_mcmc = load('testdata/results/chain_mcmcpos.txt');

% manhattan plot
figure(1); clf;
plot(p_mcmc, '.');

% compare rao-blackwellized and mcmc frequency estimates
figure(2); clf;
plot(p_mcmc, 0.5 * (p_rao1 + p_rao2), '.');

% modelsize histogram for the first chain
figure(3); clf;
hist(load('testdata/results/chain0_modelsize.dat.txt'));

% see which SNPs have posterior assocation probability > 0.5
find(p_mcmc > 0.5)

~~~~~


Contact
-------

 * E-mail: tomi.peltola at aalto.fi
 * WWW: [http://www.lce.hut.fi/research/mm/bmagwa/][homepage]
 * GitHub: [https://github.com/to-mi][GitHub repositories]

Acknowledgements: The software uses [inih][], [Boost][], [BLAS][], [LAPACK][]
and [LINPACK][] libraries.

  [inih]:    http://code.google.com/p/inih/
  [Boost]:   http://www.boost.org/
  [BLAS]:    http://www.netlib.org/blas/
  [LAPACK]:  http://www.netlib.org/lapack/
  [LINPACK]: http://www.netlib.org/linpack/


Compilation from source
-----------------------

Requirements: C, C++ and Fortran compilers (use [GNU compiler collection][]) and
BLAS, LAPACK (e.g., [ATLAS][] or [GotoBLAS][]; at least with GotoBLAS separate
LAPACK is not needed) and [Boost][] libraries.

Compilation: use the makefile provided (i.e., run `make` on command line in the 
source directory). You may have to adjust the paths of the libraries in the 
makefile.

The software have been developed only on 64 bit Linux. I have successfully 
compiled it from source under [Cygwin][] for 64-bit Windows. However, I
recommend using the provided binary on Windows.

  [Gnu compiler collection]: http://gcc.gnu.org/
  [ATLAS]:                   http://math-atlas.sourceforge.net/
  [GotoBLAS]:                http://www.tacc.utexas.edu/tacc-projects/gotoblas2
  [Cygwin]:                  http://www.cygwin.com/


License
-------

GPL v3. See `licenses/gpl-3.0.txt`.

Licenses for the utilized libraries are also included in the `licenses`
directory.

