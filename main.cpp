//////////////////////////
// Copyright (c) 2015-2019 Julian Adamek
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESSED OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//////////////////////////

//////////////////////////
// main.cpp
//////////////////////////
//
// main control sequence of Geneva N-body code with evolution of metric perturbations (gevolution)
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen Mary University of London)
//
// Last modified: November 2019
//
//////////////////////////

#include <stdlib.h>
#include <set>
#include <vector>
#ifdef HAVE_CLASS
#include "class.h"
#undef MAX			// due to macro collision this has to be done BEFORE including LATfield2 headers!
#undef MIN
#endif
#include "LATfield2.hpp"
#include "metadata.hpp"
#include "class_tools.hpp"
#include "tools.hpp"
#include "background.hpp"
#include "Particles_gevolution.hpp"
#include "gevolution.hpp"
#include "ic_basic.hpp"
#include "ic_read.hpp"
#ifdef ICGEN_PREVOLUTION
#include "ic_prevolution.hpp"
#endif
#ifdef ICGEN_FALCONIC
#include "fcn/togevolution.hpp"
#endif
#include "radiation.hpp"
#include "parser.hpp"
#include "output.hpp"
#include "hibernation.hpp"
#ifdef VELOCITY
#include "velocity.hpp"
#endif

// Quintessence
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "background_import.hpp"

using namespace std;
using namespace LATfield2;

int main(int argc, char **argv)
{
#ifdef BENCHMARK
	//benchmarking variables
	double ref_time, ref2_time, cycle_start_time;
  double a_quintessence;
	double initialization_time;
	double run_time;
  double quintessence_update_time = 0;
	double cycle_time=0;
	double projection_time = 0;
	double snapshot_output_time = 0;
	double spectra_output_time = 0;
	double lightcone_output_time = 0;
	double gravity_solver_time = 0;
	double fft_time = 0;
	int fft_count = 0;
	double update_q_time = 0;
	int update_q_count = 0;
	double moveParts_time = 0;
	int  moveParts_count =0;
#endif  //BENCHMARK

	int n = 0, m = 0;
	int io_size = 0;
	int io_group_size = 0;

	int i, j, cycle = 0, snapcount = 0, pkcount = 0, restartcount = 0, usedparams, numparam = 0, numspecies, done_hij;
	int numsteps_ncdm[MAX_PCL_SPECIES-2];
	long numpts3d;
	int box[3];
	double dtau, dtau_old, dx, tau, a, fourpiG, tmp, start_time;
	double maxvel[MAX_PCL_SPECIES];
	FILE * outfile;
	char filename[2*PARAM_MAX_LENGTH+24];
	string h5filename;
	char * settingsfile = NULL;
	char * precisionfile = NULL;
	parameter * params = NULL;
	metadata sim;
	cosmology cosmo;
	// quintessence
	mg_cosmology quintessence;
  double phi_bg, phi_p_bg, phi_pp_bg, V_prime_varphi;
	icsettings ic;
	double T00hom;

#ifndef H5_DEBUG
	H5Eset_auto2 (H5E_DEFAULT, NULL, NULL);
#endif

	for (i=1 ; i < argc ; i++ ){
		if ( argv[i][0] != '-' )
			continue;
		switch(argv[i][1]) {
			case 's':
				settingsfile = argv[++i]; //settings file name
				break;
			case 'n':
				n = atoi(argv[++i]); //size of the dim 1 of the processor grid
				break;
			case 'm':
				m =  atoi(argv[++i]); //size of the dim 2 of the processor grid
				break;
			case 'p':
#ifndef HAVE_CLASS
				cout << "HAVE_CLASS needs to be set at compilation to use CLASS precision files" << endl;
				exit(-100);
#endif
				precisionfile = argv[++i];
				break;
			case 'i':
#ifndef EXTERNAL_IO
				cout << "EXTERNAL_IO needs to be set at compilation to use the I/O server"<<endl;
				exit(-1000);
#endif
				io_size =  atoi(argv[++i]);
				break;
			case 'g':
#ifndef EXTERNAL_IO
				cout << "EXTERNAL_IO needs to be set at compilation to use the I/O server"<<endl;
				exit(-1000);
#endif
				io_group_size = atoi(argv[++i]);
		}
	}

#ifndef EXTERNAL_IO
	parallel.initialize(n,m);
#else
	if (!io_size || !io_group_size)
	{
		cout << "invalid number of I/O tasks and group sizes for I/O server (-DEXTERNAL_IO)" << endl;
		exit(-1000);
	}
	parallel.initialize(n,m,io_size,io_group_size);
	if(parallel.isIO()) ioserver.start();
	else
	{
#endif

	COUT << COLORTEXT_WHITE << endl;
	COUT << "  _  _                         " << endl;
	COUT << " |_ | |   _      _         __ ,  _" << endl;
	COUT << " |_ |_|\\ (-' \\/ (_) (_ (_| (  ( (_) /\\/	version 0.0         running on " << n*m << " cores." << endl;
	COUT << endl << " Based on gevolution version 1.2 " << endl << COLORTEXT_RESET << endl;
	//COUT << "  _   _      _         __ ,  _" << endl;
	//COUT << " (_| (-' \\/ (_) (_ (_| (  ( (_) /\\/	version 1.2         running on " << n*m << " cores." << endl;
	//COUT << "  -'" << endl << COLORTEXT_RESET << endl;

	if (settingsfile == NULL)
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": no settings file specified!" << endl;
		parallel.abortForce();
	}

	COUT << " initializing..." << endl;

	start_time = MPI_Wtime();

	numparam = loadParameterFile(settingsfile, params);

	usedparams = parseMetadata(params, numparam, sim, cosmo, quintessence, ic);

	COUT << " parsing of settings file completed. " << numparam << " parameters found, " << usedparams << " were used." << endl;
  if (quintessence.NL_quintessence == 1)
  {
    COUT<<"The quintessence is solved with non-linear corrections!"<<endl;
  }
  else COUT<<"The linear quintessence equations are being solved!"<<endl;


	sprintf(filename, "%s%s_settings_used.ini", sim.output_path, sim.basename_generic);
	saveParameterFile(filename, params, numparam);

	free(params);

#ifdef HAVE_CLASS
	background class_background;
	thermo class_thermo;
  	perturbs class_perturbs;

  	if (precisionfile != NULL)
	  	numparam = loadParameterFile(precisionfile, params);
	else
#endif
		numparam = 0;

	h5filename.reserve(2*PARAM_MAX_LENGTH);
	h5filename.assign(sim.output_path);

	box[0] = sim.numpts;
	box[1] = sim.numpts;
	box[2] = sim.numpts;

	Lattice lat(3,box,1);
	Lattice latFT;
	latFT.initializeRealFFT(lat,0);

	Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> pcls_cdm;
	Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> pcls_b;
	Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> pcls_ncdm[MAX_PCL_SPECIES-2];
	Field<Real> * update_cdm_fields[3];
	Field<Real> * update_b_fields[3];
	Field<Real> * update_ncdm_fields[3];
	double f_params[5];
	set<long> IDbacklog[MAX_PCL_SPECIES];
	// Quintessence
  Field<Real> pi;   // pi =: delta varphi
  Field<Real> V_pi; // V_pi =: conformal time derivative of delta varphi
	Field<Real> phi_old;
	Field<Real> chi_old;
	//Field<Real> phi_prime;
	//Field<Real> chi_prime;
  Field<Real> TiimT00;
	Field<Cplx> scalarFT_phi_old;
	Field<Cplx> scalarFT_chi_old;
	//Field<Cplx> scalarFT_phi_prime;
	//Field<Cplx> scalarFT_chi_prime;
  Field<Cplx> scalarFT_pi;
  Field<Cplx> scalarFT_V_pi;
  Field<Cplx> scalarFT_TiimT00;
  pi.initialize(lat,1);
	V_pi.initialize(lat,1);
	phi_old.initialize(lat,1);
	chi_old.initialize(lat,1);
	//phi_prime.initialize(lat,1);
	//chi_prime.initialize(lat,1);
  TiimT00.initialize(lat,1);
	scalarFT_pi.initialize(latFT,1);
	PlanFFT<Cplx> plan_pi(&pi, &scalarFT_pi);
  scalarFT_V_pi.initialize(latFT,1);
  PlanFFT<Cplx> plan_V_pi(&V_pi, &scalarFT_V_pi);
	scalarFT_phi_old.initialize(latFT,1);
	PlanFFT<Cplx> plan_phi_old(&phi_old, &scalarFT_phi_old);
	scalarFT_chi_old.initialize(latFT,1);
	PlanFFT<Cplx> plan_chi_old(&chi_old, &scalarFT_chi_old);
	//scalarFT_phi_prime.initialize(latFT,1);
	//PlanFFT<Cplx> plan_phi_prime(&phi_prime, &scalarFT_phi_prime);
	//scalarFT_chi_prime.initialize(latFT,1);
	//PlanFFT<Cplx> plan_chi_prime(&chi_prime, &scalarFT_chi_prime);
  scalarFT_TiimT00.initialize(latFT,1);
  PlanFFT<Cplx> plan_TiimT00(&TiimT00, &scalarFT_TiimT00);
	Field<Real> phi;
	Field<Real> source;
	Field<Real> chi;
	Field<Real> Sij;
	Field<Real> Bi;
	Field<Cplx> scalarFT;
	Field<Cplx> SijFT;
	Field<Cplx> BiFT;
	source.initialize(lat,1);
	phi.initialize(lat,1);
	chi.initialize(lat,1);
	scalarFT.initialize(latFT,1);
	PlanFFT<Cplx> plan_source(&source, &scalarFT);
	PlanFFT<Cplx> plan_phi(&phi, &scalarFT);
	PlanFFT<Cplx> plan_chi(&chi, &scalarFT);
	Sij.initialize(lat,3,3,symmetric);
	SijFT.initialize(latFT,3,3,symmetric);
	PlanFFT<Cplx> plan_Sij(&Sij, &SijFT);
	Bi.initialize(lat,3);
	BiFT.initialize(latFT,3);
	PlanFFT<Cplx> plan_Bi(&Bi, &BiFT);
#ifdef CHECK_B
	Field<Real> Bi_check;
	Field<Cplx> BiFT_check;
	Bi_check.initialize(lat,3);
	BiFT_check.initialize(latFT,3);
	PlanFFT<Cplx> plan_Bi_check(&Bi_check, &BiFT_check);
#endif
#ifdef VELOCITY
	Field<Real> vi;
	Field<Cplx> viFT;
	vi.initialize(lat,3);
	viFT.initialize(latFT,3);
	PlanFFT<Cplx> plan_vi(&vi, &viFT);
	double a_old;
#endif

	update_cdm_fields[0] = &phi;
	update_cdm_fields[1] = &chi;
	update_cdm_fields[2] = &Bi;

	update_b_fields[0] = &phi;
	update_b_fields[1] = &chi;
	update_b_fields[2] = &Bi;

	update_ncdm_fields[0] = &phi;
	update_ncdm_fields[1] = &chi;
	update_ncdm_fields[2] = &Bi;

	Site x(lat);
	rKSite kFT(latFT);

	dx = 1.0 / (double) sim.numpts;
	numpts3d = (long) sim.numpts * (long) sim.numpts * (long) sim.numpts;

	for (i = 0; i < 3; i++) // particles may never move farther than to the adjacent domain
	{
		if (lat.sizeLocal(i)-1 < sim.movelimit)
			sim.movelimit = lat.sizeLocal(i)-1;
	}
	parallel.min(sim.movelimit);

	fourpiG = 1.5 * sim.boxsize * sim.boxsize / C_SPEED_OF_LIGHT / C_SPEED_OF_LIGHT;
	a = 1. / (1. + sim.z_in);
	tau = particleHorizon(a, fourpiG, cosmo);
  #ifdef HAVE_CLASS_BG
  #ifdef HAVE_CLASS
  initializeCLASSstructures(sim, ic, cosmo, quintessence, class_background, class_thermo, class_perturbs, params, numparam);
  loadBGFunctions(class_background, quintessence, quintessence.spline_H, quintessence.acc_H, "H [1/Mpc]", sim.z_in, fourpiG);
  loadBGFunctions(class_background, quintessence, quintessence.spline_H_prime, quintessence.acc_H_prime, "H_prime", sim.z_in, fourpiG);
  loadBGFunctions(class_background, quintessence, quintessence.spline_w_mg,quintessence.acc_w_mg,"w_mg",sim.z_in, fourpiG);
  loadBGFunctions(class_background, quintessence, quintessence.spline_Omega_m,quintessence.acc_Omega_m,"Omega_m",sim.z_in, fourpiG);
  loadBGFunctions(class_background, quintessence, quintessence.spline_Omega_rad,quintessence.acc_Omega_rad,"Omega_rad",sim.z_in, fourpiG);
  loadBGFunctions(class_background, quintessence, quintessence.spline_Omega_mg,quintessence.acc_Omega_mg,"Omega_mg",sim.z_in, fourpiG);
  loadBGFunctions(class_background, quintessence, quintessence.spline_particleHorizon,quintessence.acc_particleHorizon,"conf. time [Mpc]",sim.z_in, fourpiG);
  loadBGFunctions(class_background, quintessence, quintessence.spline_mg_field,quintessence.acc_mg_field,"phi_smg",sim.z_in, fourpiG);
  loadBGFunctions(class_background, quintessence, quintessence.spline_mg_field_p,quintessence.acc_mg_field_p,"phi_prime_smg",sim.z_in, fourpiG);
  loadBGFunctions(class_background, quintessence, quintessence.spline_mg_field_pp,quintessence.acc_mg_field_pp,"phi_prime_prime_smg",sim.z_in, fourpiG);
  loadBGFunctions(class_background, quintessence, quintessence.spline_a,quintessence.acc_a,"scale factor",sim.z_in, fourpiG);
  #else
  COUT << " Quintessence background from hiclass requested while CLASS is not defined: an error occurred while importing background data." << endl;
  parallel.abortForce();
  #endif

  #else
	if( mg_import(a, fourpiG, &quintessence, quintessence.mg_bkg_file) ){
    // The import function should be in this position because it relies on fourpiG value to correctly normalize input values
    COUT << endl << " Reading quintessence background file." << endl;
    // Calling function to import background data and testing the correct import
		COUT << " File correctly imported." << endl;
		COUT << " Size of imported background vector: " << quintessence.a_vec.size() << endl;
		// Definitions for interpolation (spline method), with accelerators
    quintessence.acc_H = gsl_interp_accel_alloc();
    quintessence.spline_H = gsl_spline_alloc(gsl_interp_cspline,quintessence.last_int);
    quintessence.acc_H_prime = gsl_interp_accel_alloc();
    quintessence.spline_H_prime = gsl_spline_alloc(gsl_interp_cspline,quintessence.last_int);
    quintessence.acc_w_mg = gsl_interp_accel_alloc();
    quintessence.spline_w_mg = gsl_spline_alloc(gsl_interp_cspline,quintessence.last_int);
    quintessence.acc_Omega_m = gsl_interp_accel_alloc();
    quintessence.spline_Omega_m = gsl_spline_alloc(gsl_interp_cspline,quintessence.last_int);
    quintessence.acc_Omega_rad = gsl_interp_accel_alloc();
    quintessence.spline_Omega_rad = gsl_spline_alloc(gsl_interp_cspline,quintessence.last_int);
    quintessence.acc_Omega_mg = gsl_interp_accel_alloc();
    quintessence.spline_Omega_mg = gsl_spline_alloc(gsl_interp_cspline,quintessence.last_int);
    quintessence.acc_particleHorizon = gsl_interp_accel_alloc();
    quintessence.spline_particleHorizon = gsl_spline_alloc(gsl_interp_cspline,quintessence.last_int);
    quintessence.acc_mg_field = gsl_interp_accel_alloc();
    quintessence.spline_mg_field = gsl_spline_alloc(gsl_interp_cspline,quintessence.last_int);
    quintessence.acc_mg_field_p = gsl_interp_accel_alloc();
    quintessence.spline_mg_field_p = gsl_spline_alloc(gsl_interp_cspline,quintessence.last_int);
    quintessence.acc_a = gsl_interp_accel_alloc();
    quintessence.spline_a = gsl_spline_alloc(gsl_interp_cspline,quintessence.last_int);
		// Initialization of interpolation structures
		gsl_spline_init(quintessence.spline_H,quintessence.a,quintessence.H,quintessence.last_int);
		gsl_spline_init(quintessence.spline_H_prime,quintessence.a,quintessence.H_prime,quintessence.last_int);
		gsl_spline_init(quintessence.spline_w_mg,quintessence.a,quintessence.w_mg,quintessence.last_int);
		gsl_spline_init(quintessence.spline_Omega_m,quintessence.a,quintessence.Omega_m,quintessence.last_int);
		gsl_spline_init(quintessence.spline_Omega_rad,quintessence.a,quintessence.Omega_rad,quintessence.last_int);
		gsl_spline_init(quintessence.spline_Omega_mg,quintessence.a,quintessence.Omega_mg,quintessence.last_int);
		gsl_spline_init(quintessence.spline_particleHorizon,quintessence.a,quintessence.particleHorizon,quintessence.last_int);
	  gsl_spline_init(quintessence.spline_mg_field,quintessence.a,quintessence.mg_field,quintessence.last_int);
	  gsl_spline_init(quintessence.spline_mg_field_p,quintessence.a,quintessence.mg_field_p,quintessence.last_int);
		gsl_spline_init(quintessence.spline_a,quintessence.particleHorizon,quintessence.a,quintessence.last_int);
	}
	// If the import function fails, interrupt the software
	else{
		COUT << " Quintessence background ERROR: an error occurred while importing background data." << endl;
		parallel.abortForce();
	}
  #endif
  // Interpolation test
  if(quintessence.mg_verbose>0)
  {
    double a_eval = 1./(1+50);
    COUT << " Quintessence TEST: Hubble parameter at z_in: " << gsl_spline_eval(quintessence.spline_H,a_eval,quintessence.acc_H) << " , while GR value " << Hconf(a_eval, fourpiG, cosmo) << endl;
    COUT << " Quintessence TEST: Hubble parameter derivative at z_in: " << gsl_spline_eval(quintessence.spline_H_prime,a_eval,quintessence.acc_H_prime) << endl;
    COUT << " Quintessence TEST: Omega_m at z_in: " << gsl_spline_eval(quintessence.spline_Omega_m,a_eval,quintessence.acc_Omega_m) << " , while GR value " << Omega_m(a_eval, cosmo) << endl;
    COUT << " Quintessence TEST: Omega_rad at z_in: " << gsl_spline_eval(quintessence.spline_Omega_rad,a_eval,quintessence.acc_Omega_rad) << " , while GR value " << Omega_rad(a_eval, cosmo) << endl;
    COUT << " Quintessence TEST: Omega_Lambda at z_in: " << gsl_spline_eval(quintessence.spline_Omega_mg,a_eval,quintessence.acc_Omega_mg) << " , while GR value " << Omega_Lambda(a_eval, cosmo) << endl;
    COUT << " Quintessence TEST: particle horizon at z_in: " << gsl_spline_eval(quintessence.spline_particleHorizon, a_eval,quintessence.acc_particleHorizon) << " , while GR value " << particleHorizon(a_eval, fourpiG, cosmo) << endl;
    COUT << " Quintessence TEST: scalar field at z_in: " << gsl_spline_eval(quintessence.spline_mg_field,a_eval,quintessence.acc_mg_field) << endl;
    COUT << " Quintessence TEST: scalar field velocity at z_in: " << gsl_spline_eval(quintessence.spline_mg_field_p,a_eval,quintessence.acc_mg_field_p) << endl;
  }

	//Quintessence
	double alpha = quintessence.mg_alpha;
  // LambdaVal
	double Lambda =pow(2.*fourpiG/3. * 2.197266,0.25); //quintessence.mg_Lambda * pow(2.*fourpiG/3.,0.25);//Lambda is dimensionful so we convert to gevolution's units.
	double sigma = quintessence.mg_sigma;
  double mg_field = gsl_spline_eval(quintessence.spline_mg_field, a, quintessence.acc_mg_field);
  double mg_field_prime = gsl_spline_eval(quintessence.spline_mg_field_p, a, quintessence.acc_mg_field_p);
  double mg_field_prime_prime = gsl_spline_eval(quintessence.spline_mg_field_pp, a, quintessence.acc_mg_field_pp);
	double Hconf_quintessence = gsl_spline_eval(quintessence.spline_H,a,quintessence.acc_H);
	double Hconf_prime_quintessence = gsl_spline_eval(quintessence.spline_H_prime,a,quintessence.acc_H_prime);

	COUT << " Extended Quintessence called with parameters:" << endl;
	COUT << "\t - alpha: " << alpha << endl;
	COUT << "\t - Lambda: " << Lambda << endl;
	COUT << "\t - sigma: " << sigma << endl;

	if( (Hconf(a, fourpiG, cosmo) - Hconf_quintessence)/Hconf(a, fourpiG, cosmo) > 0.1 )
	{
		COUT << " Quintessence WARNING: It seems that the tau/boxsize is very different from LambdaCDM at z_in." << endl;
	}

	COUT << endl;

	tau = gsl_spline_eval(quintessence.spline_particleHorizon,a,quintessence.acc_particleHorizon);

	if (sim.Cf * dx < sim.steplimit / gsl_spline_eval(quintessence.spline_H,a,quintessence.acc_H))
		dtau = sim.Cf * dx;
	else
		dtau = sim.steplimit / gsl_spline_eval(quintessence.spline_H,a,quintessence.acc_H);

	dtau_old = 0.;



	if (ic.generator == ICGEN_BASIC)
		generateIC_basic(sim, ic, cosmo, quintessence, fourpiG, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &pi, &V_pi, &phi, &chi, &Bi, &source, &Sij, &scalarFT_pi, &scalarFT_V_pi, &scalarFT, &BiFT, &SijFT, &plan_pi, &plan_V_pi, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, params, numparam); // generates ICs on the fly

  else if (ic.generator == ICGEN_READ_FROM_DISK)
		readIC(sim, ic, cosmo, quintessence, fourpiG, a, tau, dtau, dtau_old, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, cycle, snapcount, pkcount, restartcount, IDbacklog);
#ifdef ICGEN_PREVOLUTION
	else if (ic.generator == ICGEN_PREVOLUTION)
		generateIC_prevolution(sim, ic, cosmo, fourpiG, a, tau, dtau, dtau_old, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, params, numparam);
#endif
#ifdef ICGEN_FALCONIC
	else if (ic.generator == ICGEN_FALCONIC)
		maxvel[0] = generateIC_FalconIC(sim, ic, cosmo, fourpiG, dtau, &pcls_cdm, pcls_ncdm, maxvel+1, &phi, &source, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_source, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
#endif
	else
	{
		COUT << " error: IC generator not implemented!" << endl;
		parallel.abortForce();
	}
	if (sim.baryon_flag > 1)
	{
		COUT << " error: baryon_flag > 1 after IC generation, something went wrong in IC generator!" << endl;
		parallel.abortForce();
	}
	numspecies = 1 + sim.baryon_flag + cosmo.num_ncdm;
	parallel.max<double>(maxvel, numspecies);

	if (sim.gr_flag > 0)
	{
		for (i = 0; i < numspecies; i++)
			maxvel[i] /= sqrt(maxvel[i] * maxvel[i] + 1.0);
	}

#ifdef CHECK_B
	if (sim.vector_flag == VECTOR_ELLIPTIC)
	{
		for (kFT.first(); kFT.test(); kFT.next())
		{
			BiFT_check(kFT, 0) = BiFT(kFT, 0);
			BiFT_check(kFT, 1) = BiFT(kFT, 1);
			BiFT_check(kFT, 2) = BiFT(kFT, 2);
		}
	}
#endif
#ifdef VELOCITY
	a_old = a;
	projection_init(&vi);
#endif

#ifdef BENCHMARK
	initialization_time = MPI_Wtime() - start_time;
	parallel.sum(initialization_time);
	COUT << COLORTEXT_GREEN << " initialization complete." << COLORTEXT_RESET << " BENCHMARK: " << hourMinSec(initialization_time) << endl << endl;
#else
	COUT << COLORTEXT_GREEN << " initialization complete." << COLORTEXT_RESET << endl << endl;
#endif

#ifdef HAVE_CLASS
	if (sim.radiation_flag > 0 || sim.fluid_flag > 0)
	{
		initializeCLASSstructures(sim, ic, cosmo, quintessence, class_background, class_thermo, class_perturbs, params, numparam);
		if (sim.gr_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.) && (ic.generator == ICGEN_BASIC || (ic.generator == ICGEN_READ_FROM_DISK && cycle == 0)))
		{
			prepareFTchiLinear(class_background, class_perturbs, scalarFT, sim, ic, cosmo, quintessence, fourpiG, a);
			plan_source.execute(FFT_BACKWARD);
			for (x.first(); x.test(); x.next())
				chi(x) += source(x);
			chi.updateHalo();
		}
	}
	if (numparam > 0) free(params);
#endif

  sprintf(filename, "%s%s_background.dat", sim.output_path, sim.basename_generic);
  outfile = fopen(filename, "w+");
  fclose(outfile);

  // EQ
  double a2, rhoa2, Pa2, f_varphi, f_varphi_prime, V_varphi;
  phi_bg = gsl_spline_eval(quintessence.spline_mg_field, a, quintessence.acc_mg_field);
  phi_p_bg = gsl_spline_eval(quintessence.spline_mg_field_p, a, quintessence.acc_mg_field_p);

	while (true)    // main loop
	{
    for (x.first(); x.test(); x.next())
    {
			phi_old(x)=phi(x);
			chi_old(x)=chi(x);
     }

#ifdef BENCHMARK
		cycle_start_time = MPI_Wtime();
#endif
		// construct stress-energy tensor
		projection_init(&source);
#ifdef HAVE_CLASS
		if (sim.radiation_flag > 0 || sim.fluid_flag > 0)
			projection_T00_project(class_background, class_perturbs, source, scalarFT, &plan_source, sim, ic, cosmo, quintessence, fourpiG, a);
#endif
		if (sim.gr_flag > 0)
		{
			projection_T00_project(&pcls_cdm, &source, a, &phi);
			if (sim.baryon_flag)
				projection_T00_project(&pcls_b, &source, a, &phi);
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (a >= 1. / (sim.z_switch_deltancdm[i] + 1.) && sim.numpcl[1+sim.baryon_flag+i] > 0)
					projection_T00_project(pcls_ncdm+i, &source, a, &phi);
				else if (sim.radiation_flag == 0 || (a >= 1. / (sim.z_switch_deltancdm[i] + 1.) && sim.numpcl[1+sim.baryon_flag+i] == 0))
				{
					tmp = bg_ncdm(a, cosmo, i);
					for(x.first(); x.test(); x.next())
						source(x) += tmp;
				}
			}
		}
		else
		{
			scalarProjectionCIC_project(&pcls_cdm, &source);
			if (sim.baryon_flag)
				scalarProjectionCIC_project(&pcls_b, &source);
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (a >= 1. / (sim.z_switch_deltancdm[i] + 1.) && sim.numpcl[1+sim.baryon_flag+i] > 0)
					scalarProjectionCIC_project(pcls_ncdm+i, &source);
			}
		}
		projection_T00_comm(&source);

#ifdef VELOCITY
		if ((sim.out_pk & MASK_VEL) || (sim.out_snapshot & MASK_VEL))
		{
			projection_init(&Bi);
            projection_Ti0_project(&pcls_cdm, &Bi, &phi, &chi);
            vertexProjectionCIC_comm(&Bi);
            compute_vi_rescaled(cosmo, &vi, &source, &Bi, a, a_old);
            a_old = a;
		}
#endif

		if (sim.vector_flag == VECTOR_ELLIPTIC)
		{
			projection_init(&Bi);
			projection_T0i_project(&pcls_cdm, &Bi, &phi);
			if (sim.baryon_flag)
				projection_T0i_project(&pcls_b, &Bi, &phi);
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (a >= 1. / (sim.z_switch_Bncdm[i] + 1.) && sim.numpcl[1+sim.baryon_flag+i] > 0)
					projection_T0i_project(pcls_ncdm+i, &Bi, &phi);
			}
			projection_T0i_comm(&Bi);
		}

		projection_init(&Sij);
		projection_Tij_project(&pcls_cdm, &Sij, a, &phi);
		if (sim.baryon_flag)
			projection_Tij_project(&pcls_b, &Sij, a, &phi);
		if (a >= 1. / (sim.z_switch_linearchi + 1.))
		{
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (sim.numpcl[1+sim.baryon_flag+i] > 0)
					projection_Tij_project(pcls_ncdm+i, &Sij, a, &phi);
			}
		}
		projection_Tij_comm(&Sij);

#ifdef BENCHMARK
		projection_time += MPI_Wtime() - cycle_start_time;
		ref_time = MPI_Wtime();
#endif

		if (sim.gr_flag > 0)
		{
			T00hom = 0.;
			for (x.first(); x.test(); x.next())
				T00hom += source(x);
			parallel.sum<double>(T00hom);
			T00hom /= (double) numpts3d;

			if (cycle % CYCLE_INFO_INTERVAL == 0)
			{
				COUT << " cycle " << cycle << ", background information: z = " << (1./a) - 1. << ", average T00 = " << T00hom << ", background model = " << cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo) << endl;
			}


			if (dtau_old > 0.)
			{
				prepareFTsource<Real>(phi, chi, source, cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo), source, 3. * Hconf(a, fourpiG, cosmo) * dx * dx / dtau_old, fourpiG * dx * dx / a, 3. * Hconf(a, fourpiG, cosmo) * Hconf(a, fourpiG, cosmo) * dx * dx);  // prepare nonlinear source for phi update
        // prepareFTsource_Quintessence(phi, chi, source, cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo), source, pi, V_pi, mg_field, mg_field_prime, alpha, Lambda, sigma, Hconf_quintessence, fourpiG, a, dx, dtau_old, quintessence.NL_quintessence); // prepare nonlinear source for phi update using extended quintessence expressions

#ifdef BENCHMARK
				ref2_time= MPI_Wtime();
#endif
				plan_source.execute(FFT_FORWARD);  // go to k-space
#ifdef BENCHMARK
				fft_time += MPI_Wtime() - ref2_time;
				fft_count++;
#endif

				solveModifiedPoissonFT(scalarFT, scalarFT, 1. / (dx * dx), 3. * Hconf(a, fourpiG, cosmo) / dtau_old);  // phi update (k-space)
        // solveModifiedPoissonFT_quintessence (scalarFT, scalarFT, 1. / (dx * dx), mg_field, mg_field_prime, alpha, Hconf_quintessence, dtau_old);


#ifdef BENCHMARK
				ref2_time= MPI_Wtime();
#endif
				plan_phi.execute(FFT_BACKWARD);	 // go back to position space
#ifdef BENCHMARK
				fft_time += MPI_Wtime() - ref2_time;
				fft_count++;
#endif
			}
		}
		else
		{
#ifdef BENCHMARK
			ref2_time= MPI_Wtime();
#endif
			plan_source.execute(FFT_FORWARD);  // Newton: directly go to k-space
#ifdef BENCHMARK
			fft_time += MPI_Wtime() - ref2_time;
			fft_count++;
#endif

			solveModifiedPoissonFT(scalarFT, scalarFT, fourpiG / a);  // Newton: phi update (k-space)

#ifdef BENCHMARK
			ref2_time= MPI_Wtime();
#endif
			plan_phi.execute(FFT_BACKWARD);	 // go back to position space
#ifdef BENCHMARK
			fft_time += MPI_Wtime() - ref2_time;
			fft_count++;
#endif
		}

		phi.updateHalo();  // communicate halo values

		// record some background data
		if (kFT.setCoord(0, 0, 0))
		{
			outfile = fopen(filename, "a");
			if (outfile == NULL)
			{
				cout << " error opening file for background output!" << endl;
			}
			else
			{
				if (cycle == 0)
					fprintf(outfile, "# background statistics\n# 0: cycle   1: tau/boxsize    2: a             3: conformal H/H0   4: scalar(phi)   5: scalar_p/H0   6: scalar_pp/H0^2    7: phi_bg(gevolution)* H0   8: phi_prime_bg(gevolution)      9: phi(k=0)       10: T00(k=0)\n");
				fprintf(outfile, " %6d   %e   %e   %e   %e    %e   %e    %e   %e    %e    %e\n", cycle, tau, a, Hconf_quintessence / gsl_spline_eval(quintessence.spline_H, 1., quintessence.acc_H), mg_field, mg_field_prime/gsl_spline_eval(quintessence.spline_H, 1., quintessence.acc_H), mg_field_prime_prime/gsl_spline_eval(quintessence.spline_H, 1., quintessence.acc_H)/gsl_spline_eval(quintessence.spline_H, 1., quintessence.acc_H), phi_bg , phi_p_bg/ gsl_spline_eval(quintessence.spline_H,1.,quintessence.acc_H), scalarFT(kFT).real(), T00hom);
				fclose(outfile);
			}
      double a2    = a * a;
      double rhoa2 = (cosmo.Omega_cdm + cosmo.Omega_b) / a + ( cosmo.Omega_Lambda * a2) + (cosmo.Omega_rad / a2);
      double Pa2   = - (cosmo.Omega_Lambda * a2) + 1./3. * (cosmo.Omega_rad / a2 );

      phi_p_bg = phi_p_bg + dtau * phi_pp_bg;
      V_prime_varphi = - sigma * Lambda * Lambda * Lambda * Lambda * pow(phi_bg, -sigma-1.);
      f_varphi =  alpha *  pow(phi_bg, 2.0);
      f_varphi_prime = 2.0 * alpha * phi_bg;
      V_varphi =  Lambda * Lambda * Lambda * Lambda * pow(phi_bg, -sigma);

      // phi_pp_bg = + ( rhoa2 - 3. * Pa2) * sqrt(2.*fourpiG/3.) * sqrt(2.*fourpiG/3.) /2. * f_varphi_prime / (1. + f_varphi + 3./4./fourpiG * f_varphi_prime *f_varphi_prime );
      // phi_pp_bg = -2 * gsl_spline_eval(quintessence.spline_H, a, quintessence.acc_H) * phi_p_bg +
      //  ( - V_prime_varphi * a2 * ( 1. +  f_varphi ) + ( rhoa2 - 3. * Pa2 + 4. * a2 * V_varphi - phi_p_bg * phi_p_bg * ( 1. + 3./2./fourpiG * 2. * alpha ) ) /2. * f_varphi_prime ) / (1. + f_varphi + 3./4./fourpiG * f_varphi_prime *f_varphi_prime );

      phi_pp_bg = -2 * gsl_spline_eval(quintessence.spline_H, a, quintessence.acc_H) * phi_p_bg- a*a * V_prime_varphi +  (gsl_spline_eval(quintessence.spline_H, a, quintessence.acc_H) * gsl_spline_eval(quintessence.spline_H, a, quintessence.acc_H) + gsl_spline_eval(quintessence.spline_H_prime, a, quintessence.acc_H_prime)) * 2. * alpha * phi_bg;

      phi_bg = phi_bg + dtau * phi_p_bg;
		}
		// done recording background data
		// prepareFTsource<Real>(phi, Sij, Sij, 2. * fourpiG * dx * dx / a);  // prepare nonlinear source for additional equations
    prepareFTsource_Tii<Real>(phi, source, Sij, Sij, TiimT00, 2. * fourpiG * dx * dx / a);  // prepare nonlinear source for additional equations
    // prepareFTsource_Quintessence(phi, Sij, Sij, chi, pi, V_pi, mg_field, mg_field_prime, alpha, Lambda, sigma, Hconf_quintessence, fourpiG, a, dx, dtau_old);


#ifdef BENCHMARK
		ref2_time= MPI_Wtime();
#endif
		plan_Sij.execute(FFT_FORWARD);  // go to k-space
#ifdef BENCHMARK
		fft_time += MPI_Wtime() - ref2_time;
		fft_count += 6;
#endif

#ifdef HAVE_CLASS
		if (sim.radiation_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.))
		{
			prepareFTchiLinear(class_background, class_perturbs, scalarFT, sim, ic, cosmo, quintessence, fourpiG, a);
			projectFTscalar(SijFT, scalarFT, 1);
		}
		else
#endif
		projectFTscalar(SijFT, scalarFT);  // construct chi by scalar projection (k-space)

#ifdef BENCHMARK
		ref2_time= MPI_Wtime();
#endif
		plan_chi.execute(FFT_BACKWARD);	 // go back to position space
#ifdef BENCHMARK
		fft_time += MPI_Wtime() - ref2_time;
		fft_count++;
#endif
		chi.updateHalo();  // communicate halo values

		if (sim.vector_flag == VECTOR_ELLIPTIC)
		{
#ifdef BENCHMARK
			ref2_time= MPI_Wtime();
#endif
			plan_Bi.execute(FFT_FORWARD);
#ifdef BENCHMARK
			fft_time += MPI_Wtime() - ref2_time;
			fft_count++;
#endif
			projectFTvector(BiFT, BiFT, fourpiG * dx * dx); // solve B using elliptic constraint (k-space)
#ifdef CHECK_B
			evolveFTvector(SijFT, BiFT_check, a * a * dtau_old);
#endif
		}
		else
			evolveFTvector(SijFT, BiFT, a * a * dtau_old);  // evolve B using vector projection (k-space)

		if (sim.gr_flag > 0)
		{
#ifdef BENCHMARK
			ref2_time= MPI_Wtime();
#endif
			plan_Bi.execute(FFT_BACKWARD);  // go back to position space
#ifdef BENCHMARK
			fft_time += MPI_Wtime() - ref2_time;
			fft_count += 3;
#endif
			Bi.updateHalo();  // communicate halo values
		}

#ifdef BENCHMARK
		gravity_solver_time += MPI_Wtime() - ref_time;
		ref_time = MPI_Wtime();
#endif

		// lightcone output
		if (sim.num_lightcone > 0)
			writeLightcones(sim, cosmo, quintessence, fourpiG, a, tau, dtau, dtau_old, maxvel[0], cycle, h5filename + sim.basename_lightcone, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &chi, &Bi, &Sij, &BiFT, &SijFT, &plan_Bi, &plan_Sij, done_hij, IDbacklog);
		else done_hij = 0;

#ifdef BENCHMARK
		lightcone_output_time += MPI_Wtime() - ref_time;
		ref_time = MPI_Wtime();
#endif

		// snapshot output
		if (snapcount < sim.num_snapshot && 1. / a < sim.z_snapshot[snapcount] + 1.)
		{
			COUT << COLORTEXT_CYAN << " writing snapshot" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;

			writeSnapshots(sim, cosmo, quintessence, fourpiG, a, dtau_old, done_hij, snapcount, h5filename + sim.basename_snapshot, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &pi, &V_pi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij
#ifdef CHECK_B
				, &Bi_check, &BiFT_check, &plan_Bi_check
#endif
#ifdef VELOCITY
				, &vi
#endif
			);

			snapcount++;
		}

#ifdef BENCHMARK
		snapshot_output_time += MPI_Wtime() - ref_time;
		ref_time = MPI_Wtime();
#endif

		// power spectra
		if (pkcount < sim.num_pk && 1. / a < sim.z_pk[pkcount] + 1.)
		{
			COUT << COLORTEXT_CYAN << " writing power spectra" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;

			writeSpectra(sim, cosmo, quintessence, fourpiG, a, pkcount,
#ifdef HAVE_CLASS
				class_background, class_perturbs, ic,
#endif
				&pcls_cdm, &pcls_b, pcls_ncdm, &phi, &pi, &V_pi, &chi, &Bi, &source, &Sij, &scalarFT, &scalarFT_pi, &scalarFT_V_pi,  &BiFT, &SijFT, &plan_phi, &plan_pi, &plan_V_pi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij
#ifdef CHECK_B
				, &Bi_check, &BiFT_check, &plan_Bi_check
#endif
#ifdef VELOCITY
				, &vi, &viFT, &plan_vi
#endif
			);

			pkcount++;
		}

#ifdef EXACT_OUTPUT_REDSHIFTS
		tmp = a;

		// Quintessence
		tmp = gsl_spline_eval(quintessence.spline_a , gsl_spline_eval(quintessence.spline_particleHorizon, tmp, quintessence.acc_particleHorizon) + 0.5 * dtau + 0.5 * dtau, quintessence.acc_a);

		//rungekutta4bg(tmp, fourpiG, cosmo, 0.5 * dtau);
		//rungekutta4bg(tmp, fourpiG, cosmo, 0.5 * dtau);

		if (pkcount < sim.num_pk && 1. / tmp < sim.z_pk[pkcount] + 1.)
		{
			writeSpectra(sim, cosmo, quintessence, fourpiG, a, pkcount,
#ifdef HAVE_CLASS
				class_background, class_perturbs, ic,
#endif
        &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &pi, &V_pi, &chi, &Bi, &source, &Sij, &scalarFT, &scalarFT_pi, &scalarFT_V_pi,  &BiFT, &SijFT, &plan_phi, &plan_pi, &plan_V_pi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij
#ifdef CHECK_B
				, &Bi_check, &BiFT_check, &plan_Bi_check
#endif
#ifdef VELOCITY
				, &vi, &viFT, &plan_vi
#endif
			);
		}
#endif // EXACT_OUTPUT_REDSHIFTS

#ifdef BENCHMARK
		spectra_output_time += MPI_Wtime() - ref_time;
#endif

		if (pkcount >= sim.num_pk && snapcount >= sim.num_snapshot)
		{
			for (i = 0; i < sim.num_lightcone; i++)
			{
				if (sim.lightcone[i].z + 1. < 1. / a)
					i = sim.num_lightcone + 1;
			}
			if (i == sim.num_lightcone) break; // simulation complete
		}

		// compute number of step subdivisions for ncdm particle updates
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			if (dtau * maxvel[i+1+sim.baryon_flag] > dx * sim.movelimit)
				numsteps_ncdm[i] = (int) ceil(dtau * maxvel[i+1+sim.baryon_flag] / dx / sim.movelimit);
			else numsteps_ncdm[i] = 1;
		}

		if (cycle % CYCLE_INFO_INTERVAL == 0)
		{
			COUT << " cycle " << cycle << ", time integration information: max |v| = " << maxvel[0] << " (cdm Courant factor = " << maxvel[0] * dtau / dx;
			if (sim.baryon_flag)
			{
				COUT << "), baryon max |v| = " << maxvel[1] << " (Courant factor = " << maxvel[1] * dtau / dx;
			}

			COUT << "), time step / Hubble time = " << Hconf_quintessence * dtau;

			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (i == 0)
				{
					COUT << endl << " time step subdivision for ncdm species: ";
				}
				COUT << numsteps_ncdm[i] << " (max |v| = " << maxvel[i+1+sim.baryon_flag] << ")";
				if (i < cosmo.num_ncdm-1)
				{
					COUT << ", ";
				}
			}

			COUT << endl;
		}


    //Kessence
#ifdef BENCHMARK
    ref_time = MPI_Wtime();
#endif

// Leap-frog quintessence - the Klein-Gordon equation.
    a_quintessence = gsl_spline_eval(quintessence.spline_a ,
                                     gsl_spline_eval(quintessence.spline_particleHorizon,a, quintessence.acc_particleHorizon) - dtau / 2.,
                                     quintessence.acc_a);
  //First we update zeta_half to have it at -1/2 just in the first loop

  if(cycle==0)
  {
      update_V_pi(phi, phi_old, chi, chi_old, pi, V_pi, TiimT00, source, mg_field, mg_field_prime, alpha, Lambda, sigma, Hconf_quintessence,  Hconf_prime_quintessence, fourpiG, a, dx, -dtau/ 2., quintessence.NL_quintessence);
      V_pi.updateHalo();
  }
 //Then we start the main loop V_pi is updated to get V_pi(n+1/2) from pi(n) and V_pi(n-1/2)
  for (i=0;i<sim.nq_numsteps;i++)
  {
      update_V_pi(phi, phi_old, chi, chi_old, pi, V_pi, TiimT00, source,
			 gsl_spline_eval(quintessence.spline_mg_field,a_quintessence,quintessence.acc_mg_field),
			 gsl_spline_eval(quintessence.spline_mg_field_p,a_quintessence,quintessence.acc_mg_field_p),
			 alpha, Lambda, sigma,
			 gsl_spline_eval(quintessence.spline_H,a_quintessence,quintessence.acc_H),
			 gsl_spline_eval(quintessence.spline_H_prime,a_quintessence,quintessence.acc_H_prime),
			 fourpiG, a_quintessence, dx, dtau/ sim.nq_numsteps, quintessence.NL_quintessence);
    V_pi.updateHalo();

		a_quintessence = gsl_spline_eval(quintessence.spline_a ,
			gsl_spline_eval(quintessence.spline_particleHorizon,a_quintessence, quintessence.acc_particleHorizon) + dtau / sim.nq_numsteps / 2.0,
			quintessence.acc_a);
    // rungekutta4bg(a_quintessence, fourpiG, cosmo,  dtau  / sim.nq_numsteps / 2.0);

		update_pi(dtau/ sim.nq_numsteps, pi, V_pi);
    pi.updateHalo();
    //********************************************************************************
    // Now we have pi(n+1) and a_kess(n+1/2) so we update background by halfstep to have a_kess(n+1)
    //********************************************************************************

		a_quintessence = gsl_spline_eval(quintessence.spline_a ,
			gsl_spline_eval(quintessence.spline_particleHorizon, a_quintessence, quintessence.acc_particleHorizon) + dtau / sim.nq_numsteps / 2.0,
			quintessence.acc_a);
    // rungekutta4bg(a_quintessence, fourpiG, cosmo,  dtau  / sim.nq_numsteps / 2.0 );
  }

#ifdef BENCHMARK
    quintessence_update_time += MPI_Wtime() - ref_time;
    ref_time = MPI_Wtime();
#endif

#ifdef BENCHMARK
		ref2_time = MPI_Wtime();
#endif
		for (i = 0; i < cosmo.num_ncdm; i++) // non-cold DM particle update
		{
			if (sim.numpcl[1+sim.baryon_flag+i] == 0) continue;

			tmp = a;

			for (j = 0; j < numsteps_ncdm[i]; j++)
			{
				f_params[0] = tmp;
				f_params[1] = tmp * tmp * sim.numpts;
				if (sim.gr_flag > 0)
					maxvel[i+1+sim.baryon_flag] = pcls_ncdm[i].updateVel(update_q, (dtau + dtau_old) / 2. / numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
				else
					maxvel[i+1+sim.baryon_flag] = pcls_ncdm[i].updateVel(update_q_Newton, (dtau + dtau_old) / 2. / numsteps_ncdm[i], update_ncdm_fields, ((sim.radiation_flag + sim.fluid_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);

#ifdef BENCHMARK
				update_q_count++;
				update_q_time += MPI_Wtime() - ref2_time;
				ref2_time = MPI_Wtime();
#endif

				// Quintessence
				tmp = gsl_spline_eval(quintessence.spline_a , gsl_spline_eval(quintessence.spline_particleHorizon, tmp, quintessence.acc_particleHorizon) + 0.5 * dtau / numsteps_ncdm[i], quintessence.acc_a);
				//rungekutta4bg(tmp, fourpiG, cosmo, 0.5 * dtau / numsteps_ncdm[i]);
				f_params[0] = tmp;
				f_params[1] = tmp * tmp * sim.numpts;

				if (sim.gr_flag > 0)
					pcls_ncdm[i].moveParticles(update_pos, dtau / numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
				else
					pcls_ncdm[i].moveParticles(update_pos_Newton, dtau / numsteps_ncdm[i], NULL, 0, f_params);
#ifdef BENCHMARK
				moveParts_count++;
				moveParts_time += MPI_Wtime() - ref2_time;
				ref2_time = MPI_Wtime();
#endif
				//quintessence
				tmp = gsl_spline_eval(quintessence.spline_a , gsl_spline_eval(quintessence.spline_particleHorizon, tmp, quintessence.acc_particleHorizon) + 0.5 * dtau / numsteps_ncdm[i], quintessence.acc_a);
				//rungekutta4bg(tmp, fourpiG, cosmo, 0.5 * dtau / numsteps_ncdm[i]);
			}
		}

		// cdm and baryon particle update
		f_params[0] = a;
		f_params[1] = a * a * sim.numpts;
		if (sim.gr_flag > 0)
		{
			maxvel[0] = pcls_cdm.updateVel(update_q, (dtau + dtau_old) / 2., update_cdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
			if (sim.baryon_flag)
				maxvel[1] = pcls_b.updateVel(update_q, (dtau + dtau_old) / 2., update_b_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
		}
		else
		{
			maxvel[0] = pcls_cdm.updateVel(update_q_Newton, (dtau + dtau_old) / 2., update_cdm_fields, ((sim.radiation_flag + sim.fluid_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);
			if (sim.baryon_flag)
				maxvel[1] = pcls_b.updateVel(update_q_Newton, (dtau + dtau_old) / 2., update_b_fields, ((sim.radiation_flag + sim.fluid_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);
		}

#ifdef BENCHMARK
		update_q_count++;
		update_q_time += MPI_Wtime() - ref2_time;
		ref2_time = MPI_Wtime();
#endif

		// Quintessence
		a = gsl_spline_eval(quintessence.spline_a , gsl_spline_eval(quintessence.spline_particleHorizon, a, quintessence.acc_particleHorizon) + 0.5 * dtau, quintessence.acc_a);

		mg_field = gsl_spline_eval(quintessence.spline_mg_field,a,quintessence.acc_mg_field);
		mg_field_prime = gsl_spline_eval(quintessence.spline_mg_field_p,a,quintessence.acc_mg_field_p);
    mg_field_prime_prime = gsl_spline_eval(quintessence.spline_mg_field_pp,a,quintessence.acc_mg_field_pp);
		Hconf_quintessence = gsl_spline_eval(quintessence.spline_H,a,quintessence.acc_H);
		Hconf_prime_quintessence = gsl_spline_eval(quintessence.spline_H_prime,a,quintessence.acc_H_prime);
		// rungekutta4bg(a, fourpiG, cosmo, 0.5 * dtau);  // evolve background by half a time step

		f_params[0] = a;
		f_params[1] = a * a * sim.numpts;
		if (sim.gr_flag > 0)
		{
			pcls_cdm.moveParticles(update_pos, dtau, update_cdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 0), f_params);
			if (sim.baryon_flag)
				pcls_b.moveParticles(update_pos, dtau, update_b_fields, (1. / a < ic.z_relax + 1. ? 3 : 0), f_params);
		}
		else
		{
			pcls_cdm.moveParticles(update_pos_Newton, dtau, NULL, 0, f_params);
			if (sim.baryon_flag)
				pcls_b.moveParticles(update_pos_Newton, dtau, NULL, 0, f_params);
		}

#ifdef BENCHMARK
		moveParts_count++;
		moveParts_time += MPI_Wtime() - ref2_time;
#endif

		// Quintessence
		a = gsl_spline_eval(quintessence.spline_a , gsl_spline_eval(quintessence.spline_particleHorizon, a, quintessence.acc_particleHorizon) + 0.5 * dtau, quintessence.acc_a);

		mg_field = gsl_spline_eval(quintessence.spline_mg_field,a,quintessence.acc_mg_field);
		mg_field_prime = gsl_spline_eval(quintessence.spline_mg_field_p,a,quintessence.acc_mg_field_p);
    mg_field_prime_prime = gsl_spline_eval(quintessence.spline_mg_field_pp,a,quintessence.acc_mg_field_pp);
		Hconf_quintessence = gsl_spline_eval(quintessence.spline_H,a,quintessence.acc_H);
		Hconf_prime_quintessence = gsl_spline_eval(quintessence.spline_H_prime,a,quintessence.acc_H_prime);
		// rungekutta4bg(a, fourpiG, cosmo, 0.5 * dtau);  // evolve background by half a time step

		parallel.max<double>(maxvel, numspecies);

		if (sim.gr_flag > 0)
		{
			for (i = 0; i < numspecies; i++)
				maxvel[i] /= sqrt(maxvel[i] * maxvel[i] + 1.0);
		}
		// done particle update

		tau += dtau;

		if (sim.wallclocklimit > 0.)   // check for wallclock time limit
		{
			tmp = MPI_Wtime() - start_time;
			parallel.max(tmp);
			if (tmp > sim.wallclocklimit)   // hibernate
			{
				COUT << COLORTEXT_YELLOW << " reaching hibernation wallclock limit, hibernating..." << COLORTEXT_RESET << endl;
				COUT << COLORTEXT_CYAN << " writing hibernation point" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;
				if (sim.vector_flag == VECTOR_PARABOLIC && sim.gr_flag == 0)
					plan_Bi.execute(FFT_BACKWARD);
#ifdef CHECK_B
				if (sim.vector_flag == VECTOR_ELLIPTIC)
				{
					plan_Bi_check.execute(FFT_BACKWARD);
					hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi_check, a, tau, dtau, cycle);
				}
				else
#endif
				hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi, a, tau, dtau, cycle);
				break;
			}
		}

		if (restartcount < sim.num_restart && 1. / a < sim.z_restart[restartcount] + 1.)
		{
			COUT << COLORTEXT_CYAN << " writing hibernation point" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;
			if (sim.vector_flag == VECTOR_PARABOLIC && sim.gr_flag == 0)
				plan_Bi.execute(FFT_BACKWARD);
#ifdef CHECK_B
			if (sim.vector_flag == VECTOR_ELLIPTIC)
			{
				plan_Bi_check.execute(FFT_BACKWARD);
				hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi_check, a, tau, dtau, cycle, restartcount);
			}
			else
#endif
			hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi, a, tau, dtau, cycle, restartcount);
			restartcount++;
		}

		dtau_old = dtau;

		// Quintessence
		if(quintessence.mg_verbose>0)
		{
			COUT << " Quintessence VERBOSE: Hubble parameter check at a = " << a << " : " << Hconf_quintessence<< " , while GR value " << Hconf(a, fourpiG, cosmo) << endl;
		}

		if (sim.Cf * dx < sim.steplimit / Hconf_quintessence)
			dtau = sim.Cf * dx;
		else
			dtau = sim.steplimit / Hconf_quintessence;

		cycle++;

#ifdef BENCHMARK
		cycle_time += MPI_Wtime()-cycle_start_time;
#endif
	}

	COUT << COLORTEXT_GREEN << " simulation complete." << COLORTEXT_RESET << endl;

#ifdef BENCHMARK
		ref_time = MPI_Wtime();
#endif

#ifdef HAVE_CLASS
	if (sim.radiation_flag > 0 || sim.fluid_flag > 0)
		freeCLASSstructures(class_background, class_thermo, class_perturbs);
#endif

// Quintessence
// Free interpolation structures
gsl_spline_free(quintessence.spline_H);gsl_interp_accel_free(quintessence.acc_H);
gsl_spline_free(quintessence.spline_H_prime);gsl_interp_accel_free(quintessence.acc_H_prime);
gsl_spline_free(quintessence.spline_w_mg);gsl_interp_accel_free(quintessence.acc_w_mg);
gsl_spline_free(quintessence.spline_Omega_m);gsl_interp_accel_free(quintessence.acc_Omega_m);
gsl_spline_free(quintessence.spline_Omega_rad);gsl_interp_accel_free(quintessence.acc_Omega_rad);
gsl_spline_free(quintessence.spline_Omega_mg);gsl_interp_accel_free(quintessence.acc_Omega_mg);
gsl_spline_free(quintessence.spline_mg_field);gsl_interp_accel_free(quintessence.acc_mg_field);
gsl_spline_free(quintessence.spline_mg_field_p);gsl_interp_accel_free(quintessence.acc_mg_field_p);
gsl_spline_free(quintessence.spline_particleHorizon);gsl_interp_accel_free(quintessence.acc_particleHorizon);

gsl_spline_free(quintessence.spline_a);gsl_interp_accel_free(quintessence.acc_a);

COUT << endl << " Quintessence interpolation structures correctly deallocated." << endl << endl;

#ifdef BENCHMARK
	lightcone_output_time += MPI_Wtime() - ref_time;
	run_time = MPI_Wtime() - start_time;

	parallel.sum(run_time);
	parallel.sum(cycle_time);
	parallel.sum(projection_time);
	parallel.sum(snapshot_output_time);
	parallel.sum(spectra_output_time);
	parallel.sum(lightcone_output_time);
	parallel.sum(gravity_solver_time);
	parallel.sum(fft_time);
	parallel.sum(update_q_time);
	parallel.sum(moveParts_time);

	COUT << endl << "BENCHMARK" << endl;
	COUT << "total execution time  : "<<hourMinSec(run_time) << endl;
	COUT << "total number of cycles: "<< cycle << endl;
	COUT << "time consumption breakdown:" << endl;
	COUT << "initialization   : "  << hourMinSec(initialization_time) << " ; " << 100. * initialization_time/run_time <<"%."<<endl;
	COUT << "main loop        : "  << hourMinSec(cycle_time) << " ; " << 100. * cycle_time/run_time <<"%."<<endl;

	COUT << "----------- main loop: components -----------"<<endl;
	COUT << "projections                : "<< hourMinSec(projection_time) << " ; " << 100. * projection_time/cycle_time <<"%."<<endl;
	COUT << "snapshot outputs           : "<< hourMinSec(snapshot_output_time) << " ; " << 100. * snapshot_output_time/cycle_time <<"%."<<endl;
	COUT << "lightcone outputs          : "<< hourMinSec(lightcone_output_time) << " ; " << 100. * lightcone_output_time/cycle_time <<"%."<<endl;
	COUT << "power spectra outputs      : "<< hourMinSec(spectra_output_time) << " ; " << 100. * spectra_output_time/cycle_time <<"%."<<endl;
	COUT << "update momenta (count: "<<update_q_count <<"): "<< hourMinSec(update_q_time) << " ; " << 100. * update_q_time/cycle_time <<"%."<<endl;
	COUT << "move particles (count: "<< moveParts_count <<"): "<< hourMinSec(moveParts_time) << " ; " << 100. * moveParts_time/cycle_time <<"%."<<endl;
	COUT << "gravity solver             : "<< hourMinSec(gravity_solver_time) << " ; " << 100. * gravity_solver_time/cycle_time <<"%."<<endl;
	COUT << "-- thereof Fast Fourier Transforms (count: " << fft_count <<"): "<< hourMinSec(fft_time) << " ; " << 100. * fft_time/gravity_solver_time <<"%."<<endl;
#endif

#ifdef EXTERNAL_IO
		ioserver.stop();
	}
#endif

	return 0;
}
