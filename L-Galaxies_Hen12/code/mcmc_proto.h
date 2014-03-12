

// Routines used in the MCMC sampling
void Senna (void);
int mcmc(int Nsteps);
double SAM(int filenr);
void initialize_mcmc_par (void);
void read_mcmc_par (int snapnum);
void read_sample_info(void);
void read_observations(void);
double propose_new_parameters(void);
void save_galaxy_for_mcmc(int gal_index);
void write_comparison_to_obs(void);
void reset_cosmology(void);

double get_likelihood (void);
void bin_function(int ObsNr, double *binsamdata, double *samdata, int snap);
void bin_colors(int ObsNr, double *binredfraction, int snap);
void bin_color_hist(int ObsNr, double *bincolorhist, int snap);
void bin_bhbm(double *binblackholeup, double *binblackholedown, int snap);
double chi_square_probability(int ObsNr, double *samdata, int snap);
double maximum_likelihood_probability(int ObsNr, double *samfract, int snap);
double binomial_probability(int ObsNr, double *samup, double *samdown, int snap);


double gammp (double a, double x);
double gammq (double a, double x);
double gser (double a, double x);
double gcf(double a, double x);
double gammpapprox(double a, double x, int psig);
double gammln(double xx);
double betai(double a, double b, double x);
double betacf(double a,double  b, double x);
double ran3(long *idum);
double gassdev(long *idum);

