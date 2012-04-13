BMAGWA software v2.0 README
===========================

Introduction
------------

This software implements a method for computing posterior association 
probabilities of SNPs (and other quantities) in genome-wide association studies 
using Bayesian variable selection and model averaging (Peltola et al., 2012).

References:

 * Peltola T, Marttinen P, Jula A, Salomaa V, Perola M, Vehtari A (2012)
   Bayesian Variable Selection in Searching for Additive and Dominant
   Effects in Genome-Wide Data. PLoS ONE 7(1): e29115.
   http://www.plosone.org/article/info:doi/10.1371/journal.pone.0029115


Availability
------------

The source code is available at the [homepage][] of the software or from their
[GitHub repositories][]. See below for compilation. There are no ready-made
binary executables at the moment, but if you need one, send me an e-mail.

The software has been developed (and tested) on 64bit Linux.

  [homepage]:            http://becs.aalto.fi/en/research/bayes/bmagwa/
  [GitHub repositories]: https://github.com/to-mi


Compilation from source
-----------------------

Requirements: C, C++ and Fortran compilers (use [GNU compiler collection][]) and
BLAS, LAPACK (e.g., [ATLAS][] or [OpenBLAS][]) and [Boost][] libraries.

Compilation: use the makefile provided (i.e., run `make` on command line in the 
source directory). You may have to adjust the paths of the libraries in the 
makefile.

The software have been developed only on 64 bit Linux. Complitation under
Windows (using [Cygwin]) should be possible, albeit minor modification may be
neccessary.

  [Gnu compiler collection]: http://gcc.gnu.org/
  [ATLAS]:                   http://math-atlas.sourceforge.net/
  [OpenBLAS]:                https://github.com/xianyi/OpenBLAS/
  [Cygwin]:                  http://www.cygwin.com/


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
; Data files are in PLINK format (use only standard formats; not long or
; transposed etc.; see http://pngu.mgh.harvard.edu/~purcell/plink/).
;
; fam contains information on individuals (only phenotype is used).
file_fam = testdata/testdata.fam
; Genotypes in binary format (SNP-major mode; this is default in PLINK).
file_g = testdata/testdata.bed
; Covariates (excluding constant, which is added automatically).
file_e = testdata/testdata.e
; Optional phenotype file (overruns phenotypes of fam if provided).
file_y = testdata/testdata.y
; Set following to 1 to force coding of genotypes to 0,1,2 according to the
; number of minor alleles (i.e. if computing allele freq with current coding
; gives > 0.5, 0s are swapped to 2s and 2s to 0s). Usually you probably use 0.
recode_g_to_minor_allele_count = 1

[sizes]
; Number of individuals.
n = 1000
; Number of SNPs.
m_g = 10000
; Number of covariates (excluding constant).
m_e = 2

[sampler]
; Sampler types: PMV / NK / KSC / G.
; The latter three have limited implementation (no effect types other than A
; can be used).
type = PMV
; Number of iterations.
do_n_iter = 1000000
; Skip (in iterations) between computing rao-blackwellization (RB).
n_rao = 500
; Number of RB-steps to discard as burnin (affects _rao.dat outputfiles)
; relative to n_rao (so these settings start computing RB posterior estimates
; after 500*1000 iterations).
n_rao_burnin = 1000
; Keep at 0 (experimental; set to 1 to continue adaptations after burnin).
adaptation = 0
; Skip between printing status to _log.txt file.
verbosity = 100000
; Skip between saving MCMC samples.
thin = 10
; Skip between updating tau2 and missing genotypes.
n_sample_tau2_and_missing = 10
; Delayed rejection: 0 off, > 0 gives the threshold for restricting
; DR moves to move sizes less than or equal to the value.
delay_rejection = 10
; Initial value of the move size proposal distribution parameter (geom. dist).
p_move_size = 0.2
; The move size proposal distribution parameter for state changes of near-by
; snps.
p_move_size_nbc = 0.25
; The move size proposal distribution parameter for swaps of near-by snps.
p_move_size_nbs = 0.7
; SNP neighborhood size / 2 (i.e., to one direction as difference in the
; indices of snps) for state changes and swaps of near-by SNPs.
max_SNP_neighborhood_size = 20
; Adapt move size distribtion: 0 or 1.
adapt_p_move_size = 1
; Set to 0 to use expected jump distance optimization.
; Set to (0,1] to use acceptance rate coercion with the specified goal value
; for the rate.
; Has not effect is adapt_p_move_size = 0.
p_move_size_acpt_goal = 0
; Maximum move size: [1,255] (not tested with values > 20).
max_move_size = 20
; Don't adapt the proposal distribution for additions/removals (i.e., use
; uniform dist.): 0 or 1.
flat_proposal_dist = 0

[thread]
; Number of chains to run simultanously.
n_threads = 2
; Prefix for result files (directory should exist; in this case
; 'testdata/results/').
basename = testdata/results/chain
; Seed for random number generator for each chain (comma separated list).
seeds = 1234,2345

[model]
; Allowed effect types (possible values are A,H,D,R,AH; comma separated list).
types = A

[prior]
; Use individual tau2 parameters: 0 or 1.
use_individual_tau2 = 1
; Prior weights on different effect types (t parameter prior).
type_A = 1
type_H = 1
type_D = 1
type_R = 1
type_AH = 1
; Model size prior (mean and variance) used to set prior for gamma parameter.
e_qg = 20
var_qg = 300
; Residual variance (sigma2 parameter) prior.
; Give either R2mode_sigma2 (expected mode for proportion of variance
; explained) or s2_sigma2 (directly the parameter of invchi2-prior).
R2mode_sigma2 = 0.2
;s2_sigma2 = 1.0
nu_sigma2 = 1
; Effect size (tau2_t parameters) priors.
; _A is for additive effect. _H, _R, _D for others.
nu_tau2_A = 5
s2_tau2_A = 0.05
nu_tau2_H = 5
s2_tau2_H = 0.05
; Alpha prior mean value.
mu_alpha = 1
; Constant term effect size prior parameter (1/sigma^2_alpha_0).
inv_tau2_e_const_val = 0
; Covariate term effect size prior parameter (1/sigma^2_alpha_l).
inv_tau2_e_val = 0.001

~~~~~


Output files and formats
------------------------

Most of the output files are in binary format. A Python script/module to convert
the results to text file format is included (run `python postprocess.py` for 
help). Note that the binary format may not be portable across computers with
different architectures.

The output file names are prefixed with the `basename`-setting and the number of
the chain.

Note: thinning affects the output. If `do_n_iter = 100` and `thin = 10` only 10 
numbers will be in the (MCMC trace) output files.

### Text files

 * `basename_prior.txt`

   Contains the values of prior parameters.

 * `basenameX_log.txt`

   The state of the sampler is printed to this file every verbosity iterations.
   You can check this file to monitor the progress.

 * `basenameX_samplerstats.txt`

   Some statistics from the run (sampling times etc.).


### Binary files

Format of a single value is given after filename.

 * `basenameX_alpha.dat [double ('d')]`

   Contains the value of alpha for each MCMC iteration.

 * `basenameX_jumpdistance.dat [unsigned char ('B')]`

   Contains the realized number of changes made to the variable inclusion
   vector for each MCMC iteration (not thinned!).

 * `basenameX_loci.dat [unsigned int32 ('I')]`

   Contains the SNPs included in the model for each MCMC iteration. You need
   `basenameX_modelsize.dat` to parse this file.

 * `basenameX_log_likelihood.dat [double ('d')]`

   Contains the marginal log likelihood (conditional on beta variance
   parameters) of the linear model for each MCMC iteration.

 * `basenameX_log_prior.dat [double ('d')]`

   Contains the log of model prior (model size and effect type prior) for each
   MCMC iteration. Does not include prior contribution of the beta variance
   parameters.

 * `basenameX_modelsize.dat [unsigned int32 ('I')]`

   Contains the number of SNPs in the model for each MCMC iteration.

 * `basenameX_move_size.dat [unsigned char ('B')]`

   Contains the size of move proposed for each MCMC iteration. Not thinned.

 * `basenameX_move_type.dat [unsigned char ('B')]`

   Contains the type of move proposed for each MCMC iteration. Not thinned.

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

 * `basenameX_types.dat [unsigned char ('B')]`

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
    BMAGWA software version 2.0
    http://becs.aalto.fi/en/research/bayes/bmagwa/
    -------------------------------------------------------------
    
    Warning: using default value of 0 for option sampler.save_beta
    Warning: using default value of 0 for option prior.s2_sigma2
    Warning: using default value of 0 for option prior.eh_tau2_A
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
    python postprocess.py mcmcpos testdata/results/chain 10000 50000 1

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


Limitations in the implementation
---------------------------------

 * Other sampler types than PMV do not allow effect types other than A.
 * Near-by SNP state change proposal is not available when there are more than
   one effect type in use.
 * Acceptance rate coercion is experimental (not tested / will not work in all
   cases).


Contact
-------

 * E-mail: tomi.peltola at aalto.fi
 * WWW: [http://becs.aalto.fi/en/research/bayes/bmagwa/][homepage]
 * GitHub: [https://github.com/to-mi][GitHub repositories]

Acknowledgements: The software uses [inih][], [Boost][], [BLAS][], [LAPACK][]
and [LINPACK][] libraries.

  [inih]:    http://code.google.com/p/inih/
  [Boost]:   http://www.boost.org/
  [BLAS]:    http://www.netlib.org/blas/
  [LAPACK]:  http://www.netlib.org/lapack/
  [LINPACK]: http://www.netlib.org/linpack/


License
-------

GPL v3. See `licenses/gpl-3.0.txt`.

Licenses for the utilized libraries are also included in the `licenses`
directory.

