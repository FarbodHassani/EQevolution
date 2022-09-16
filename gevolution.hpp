//////////////////////////
// gevolution.hpp
//////////////////////////
//
// Geneva algorithms for evolution of metric perturbations
// and relativistic free-streaming particles (gevolution)
//
// 1. Suite of Fourier-based methods for the computation of the
//    relativistic scalar (Phi, Phi-Psi) and vector modes [see J. Adamek,
//    R. Durrer, and M. Kunz, Class. Quant. Grav. 31, 234006 (2014)]
//
// 2. Collection of "update position" and "update velocity/momentum" methods
//    [see J. Adamek, D. Daverio, R. Durrer, and M. Kunz, JCAP 1607, 053 (2016)]
//
// 3. Collection of projection methods for the construction of the
//    stress-energy-tensor
//
// 4. Fourier-space projection methods for the computation of the
//    curl and divergence of the velocity field
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen Mary University of London)
//
// Last modified: April 2019
//
//////////////////////////

#ifndef GEVOLUTION_HEADER
#define GEVOLUTION_HEADER

#ifndef Cplx
#define Cplx Imag
#endif

using namespace std;
using namespace LATfield2;



//////////////////////////
// prepareFTsource (1)
//////////////////////////
// Description:
//   construction of real-space source tensor for Fourier-based solvers
//
// Arguments:
//   phi        reference to field configuration
//   Tij        reference to symmetric tensor field containing the space-space
//              components of the stress-energy tensor (rescaled by a^3)
//   Sij        reference to allocated symmetric tensor field which will contain
//              the source tensor (may be identical to Tji)
//   coeff      scaling coefficient for Tij ("8 pi G dx^2 / a")
//
// Returns:
//
//////////////////////////

template <class FieldType>
void prepareFTsource_Tii(Field<FieldType> & phi, Field<FieldType> & source, Field<FieldType> & Tij, Field<FieldType> & Sij, Field<FieldType> & source_klein_gordon, const double coeff)
{
	Site x(phi.lattice());

	for (x.first(); x.test(); x.next())
	{
    source_klein_gordon (x) = Tij(x, 0, 0) + Tij(x, 1, 1) +  Tij(x, 2, 2)  -  source (x);
		// 0-0-component:
		Sij(x, 0, 0) = coeff * Tij(x, 0, 0);
#ifdef PHINONLINEAR
#ifdef ORIGINALMETRIC
		Sij(x, 0, 0) -= 4. * phi(x) * (phi(x-0) + phi(x+0) - 2. * phi(x));
		Sij(x, 0, 0) -= 0.5 * (phi(x+0) - phi(x-0)) * (phi(x+0) - phi(x-0));
#else
		Sij(x, 0, 0) += 0.5 * (phi(x+0) - phi(x-0)) * (phi(x+0) - phi(x-0));
#endif
#endif

		// 1-1-component:
		Sij(x, 1, 1) = coeff * Tij(x, 1, 1);
#ifdef PHINONLINEAR
#ifdef ORIGINALMETRIC
		Sij(x, 1, 1) -= 4. * phi(x) * (phi(x-1) + phi(x+1) - 2. * phi(x));
		Sij(x, 1, 1) -= 0.5 * (phi(x+1) - phi(x-1)) * (phi(x+1) - phi(x-1));
#else
		Sij(x, 1, 1) += 0.5 * (phi(x+1) - phi(x-1)) * (phi(x+1) - phi(x-1));
#endif
#endif

		// 2-2-component:
		Sij(x, 2, 2) = coeff * Tij(x, 2, 2);
#ifdef PHINONLINEAR
#ifdef ORIGINALMETRIC
		Sij(x, 2, 2) -= 4. * phi(x) * (phi(x-2) + phi(x+2) - 2. * phi(x));
		Sij(x, 2, 2) -= 0.5 * (phi(x+2) - phi(x-2)) * (phi(x+2) - phi(x-2));
#else
		Sij(x, 2, 2) += 0.5 * (phi(x+2) - phi(x-2)) * (phi(x+2) - phi(x-2));
#endif
#endif

		// 0-1-component:
		Sij(x, 0, 1) = coeff * Tij(x, 0, 1);
#ifdef PHINONLINEAR
		Sij(x, 0, 1) += phi(x+0) * phi(x+1) - phi(x) * phi(x+0+1);
#ifdef ORIGINALMETRIC
		Sij(x, 0, 1) -= 1.5 * phi(x) * phi(x);
		Sij(x, 0, 1) += 1.5 * phi(x+0) * phi(x+0);
		Sij(x, 0, 1) += 1.5 * phi(x+1) * phi(x+1);
		Sij(x, 0, 1) -= 1.5 * phi(x+0+1) * phi(x+0+1);
#else
		Sij(x, 0, 1) += 0.5 * phi(x) * phi(x);
		Sij(x, 0, 1) -= 0.5 * phi(x+0) * phi(x+0);
		Sij(x, 0, 1) -= 0.5 * phi(x+1) * phi(x+1);
		Sij(x, 0, 1) += 0.5 * phi(x+0+1) * phi(x+0+1);
#endif
#endif

		// 0-2-component:
		Sij(x, 0, 2) = coeff * Tij(x, 0, 2);
#ifdef PHINONLINEAR
		Sij(x, 0, 2) += phi(x+0) * phi(x+2) - phi(x) * phi(x+0+2);
#ifdef ORIGINALMETRIC
		Sij(x, 0, 2) -= 1.5 * phi(x) * phi(x);
		Sij(x, 0, 2) += 1.5 * phi(x+0) * phi(x+0);
		Sij(x, 0, 2) += 1.5 * phi(x+2) * phi(x+2);
		Sij(x, 0, 2) -= 1.5 * phi(x+0+2) * phi(x+0+2);
#else
		Sij(x, 0, 2) += 0.5 * phi(x) * phi(x);
		Sij(x, 0, 2) -= 0.5 * phi(x+0) * phi(x+0);
		Sij(x, 0, 2) -= 0.5 * phi(x+2) * phi(x+2);
		Sij(x, 0, 2) += 0.5 * phi(x+0+2) * phi(x+0+2);
#endif
#endif

		// 1-2-component:
		Sij(x, 1, 2) = coeff * Tij(x, 1, 2);
#ifdef PHINONLINEAR
		Sij(x, 1, 2) += phi(x+1) * phi(x+2) - phi(x) * phi(x+1+2);
#ifdef ORIGINALMETRIC
		Sij(x, 1, 2) -= 1.5 * phi(x) * phi(x);
		Sij(x, 1, 2) += 1.5 * phi(x+1) * phi(x+1);
		Sij(x, 1, 2) += 1.5 * phi(x+2) * phi(x+2);
		Sij(x, 1, 2) -= 1.5 * phi(x+1+2) * phi(x+1+2);
#else
		Sij(x, 1, 2) += 0.5 * phi(x) * phi(x);
		Sij(x, 1, 2) -= 0.5 * phi(x+1) * phi(x+1);
		Sij(x, 1, 2) -= 0.5 * phi(x+2) * phi(x+2);
		Sij(x, 1, 2) += 0.5 * phi(x+1+2) * phi(x+1+2);
#endif
#endif
	}
}



//////////////////////////
// prepareFTsource (1)
//////////////////////////
// Description:
//   construction of real-space source tensor for Fourier-based solvers
//
// Arguments:
//   phi        reference to field configuration
//   Tij        reference to symmetric tensor field containing the space-space
//              components of the stress-energy tensor (rescaled by a^3)
//   Sij        reference to allocated symmetric tensor field which will contain
//              the source tensor (may be identical to Tji)
//   coeff      scaling coefficient for Tij ("8 pi G dx^2 / a")
//
// Returns:
//
//////////////////////////

template <class FieldType>
void prepareFTsource(Field<FieldType> & phi, Field<FieldType> & Tij, Field<FieldType> & Sij, const double coeff)
{
	Site x(phi.lattice());

	for (x.first(); x.test(); x.next())
	{
		// 0-0-component:
		Sij(x, 0, 0) = coeff * Tij(x, 0, 0);
#ifdef PHINONLINEAR
#ifdef ORIGINALMETRIC
		Sij(x, 0, 0) -= 4. * phi(x) * (phi(x-0) + phi(x+0) - 2. * phi(x));
		Sij(x, 0, 0) -= 0.5 * (phi(x+0) - phi(x-0)) * (phi(x+0) - phi(x-0));
#else
		Sij(x, 0, 0) += 0.5 * (phi(x+0) - phi(x-0)) * (phi(x+0) - phi(x-0));
#endif
#endif

		// 1-1-component:
		Sij(x, 1, 1) = coeff * Tij(x, 1, 1);
#ifdef PHINONLINEAR
#ifdef ORIGINALMETRIC
		Sij(x, 1, 1) -= 4. * phi(x) * (phi(x-1) + phi(x+1) - 2. * phi(x));
		Sij(x, 1, 1) -= 0.5 * (phi(x+1) - phi(x-1)) * (phi(x+1) - phi(x-1));
#else
		Sij(x, 1, 1) += 0.5 * (phi(x+1) - phi(x-1)) * (phi(x+1) - phi(x-1));
#endif
#endif

		// 2-2-component:
		Sij(x, 2, 2) = coeff * Tij(x, 2, 2);
#ifdef PHINONLINEAR
#ifdef ORIGINALMETRIC
		Sij(x, 2, 2) -= 4. * phi(x) * (phi(x-2) + phi(x+2) - 2. * phi(x));
		Sij(x, 2, 2) -= 0.5 * (phi(x+2) - phi(x-2)) * (phi(x+2) - phi(x-2));
#else
		Sij(x, 2, 2) += 0.5 * (phi(x+2) - phi(x-2)) * (phi(x+2) - phi(x-2));
#endif
#endif

		// 0-1-component:
		Sij(x, 0, 1) = coeff * Tij(x, 0, 1);
#ifdef PHINONLINEAR
		Sij(x, 0, 1) += phi(x+0) * phi(x+1) - phi(x) * phi(x+0+1);
#ifdef ORIGINALMETRIC
		Sij(x, 0, 1) -= 1.5 * phi(x) * phi(x);
		Sij(x, 0, 1) += 1.5 * phi(x+0) * phi(x+0);
		Sij(x, 0, 1) += 1.5 * phi(x+1) * phi(x+1);
		Sij(x, 0, 1) -= 1.5 * phi(x+0+1) * phi(x+0+1);
#else
		Sij(x, 0, 1) += 0.5 * phi(x) * phi(x);
		Sij(x, 0, 1) -= 0.5 * phi(x+0) * phi(x+0);
		Sij(x, 0, 1) -= 0.5 * phi(x+1) * phi(x+1);
		Sij(x, 0, 1) += 0.5 * phi(x+0+1) * phi(x+0+1);
#endif
#endif

		// 0-2-component:
		Sij(x, 0, 2) = coeff * Tij(x, 0, 2);
#ifdef PHINONLINEAR
		Sij(x, 0, 2) += phi(x+0) * phi(x+2) - phi(x) * phi(x+0+2);
#ifdef ORIGINALMETRIC
		Sij(x, 0, 2) -= 1.5 * phi(x) * phi(x);
		Sij(x, 0, 2) += 1.5 * phi(x+0) * phi(x+0);
		Sij(x, 0, 2) += 1.5 * phi(x+2) * phi(x+2);
		Sij(x, 0, 2) -= 1.5 * phi(x+0+2) * phi(x+0+2);
#else
		Sij(x, 0, 2) += 0.5 * phi(x) * phi(x);
		Sij(x, 0, 2) -= 0.5 * phi(x+0) * phi(x+0);
		Sij(x, 0, 2) -= 0.5 * phi(x+2) * phi(x+2);
		Sij(x, 0, 2) += 0.5 * phi(x+0+2) * phi(x+0+2);
#endif
#endif

		// 1-2-component:
		Sij(x, 1, 2) = coeff * Tij(x, 1, 2);
#ifdef PHINONLINEAR
		Sij(x, 1, 2) += phi(x+1) * phi(x+2) - phi(x) * phi(x+1+2);
#ifdef ORIGINALMETRIC
		Sij(x, 1, 2) -= 1.5 * phi(x) * phi(x);
		Sij(x, 1, 2) += 1.5 * phi(x+1) * phi(x+1);
		Sij(x, 1, 2) += 1.5 * phi(x+2) * phi(x+2);
		Sij(x, 1, 2) -= 1.5 * phi(x+1+2) * phi(x+1+2);
#else
		Sij(x, 1, 2) += 0.5 * phi(x) * phi(x);
		Sij(x, 1, 2) -= 0.5 * phi(x+1) * phi(x+1);
		Sij(x, 1, 2) -= 0.5 * phi(x+2) * phi(x+2);
		Sij(x, 1, 2) += 0.5 * phi(x+1+2) * phi(x+1+2);
#endif
#endif
	}
}


//////////////////////////
// prepareFTsource (2)
//////////////////////////
// Description:
//   construction of real-space source field for Fourier-based solvers
//
// Arguments:
//   phi        reference to field configuration (first Bardeen potential)
//   chi        reference to field configuration (difference between Bardeen potentials, phi-psi)
//   source     reference to fully dressed source field (rescaled by a^3)
//   bgmodel    background model of the source (rescaled by a^3) to be subtracted
//   result     reference to allocated field which will contain the result (may be identical to source)
//   coeff      diffusion coefficient ("3 H_conformal dx^2 / dtau")
//   coeff2     scaling coefficient for the source ("4 pi G dx^2 / a")
//   coeff3     scaling coefficient for the psi-term ("3 H_conformal^2 dx^2")
// Returns:
//
//////////////////////////

template <class FieldType>
void prepareFTsource(Field<FieldType> & phi, Field<FieldType> & chi, Field<FieldType> & source, const FieldType bgmodel, Field<FieldType> & result, const double coeff, const double coeff2, const double coeff3)
{

	Site x(phi.lattice());

	for (x.first(); x.test(); x.next())
	{
		result(x) = coeff2 * (source(x) - bgmodel);
#ifdef PHINONLINEAR
#ifdef ORIGINALMETRIC
    result(x) *= 1. - 4. * phi(x);
    result(x) -= 0.375 * (phi(x-0) - phi(x+0)) * (phi(x-0) - phi(x+0));
    result(x) -= 0.375 * (phi(x-1) - phi(x+1)) * (phi(x-1) - phi(x+1));
    result(x) -= 0.375 * (phi(x-2) - phi(x+2)) * (phi(x-2) - phi(x+2));
#else
		result(x) *= 1. - 2. * phi(x);
		result(x) += 0.125 * (phi(x-0) - phi(x+0)) * (phi(x-0) - phi(x+0));
		result(x) += 0.125 * (phi(x-1) - phi(x+1)) * (phi(x-1) - phi(x+1));
		result(x) += 0.125 * (phi(x-2) - phi(x+2)) * (phi(x-2) - phi(x+2));
#endif
#endif
		result(x) += (coeff3 - coeff) * phi(x) - coeff3 * chi(x);
	}
}


//////////////////////////
// prepareFTsource (1)
//////////////////////////
// Description:
//   construction of real-space source tensor for Fourier-based solvers
//
// Arguments:
//   phi        reference to field configuration (first Bardeen potential)
//   chi        reference to field configuration (difference between Bardeen potentials, phi-psi)
//   Tij        reference to symmetric tensor field containing the space-space
//              components of the stress-energy tensor (rescaled by a^3)
//   Sij        reference to allocated symmetric tensor field which will contain
//              the source tensor (may be identical to Tji)
//   QUINTESSENCE
//   pi              reference to scalar field in extended quintessence
//   V_pi            reference to time derivative (conformal time) of the scalar field in e-quintessence
//   varphi_bg       reference to backgroudn value of the scalar field in extended quintessence
//   varphi_prime_bg reference to conformal time derivative of varphi extended quintessence
//   alpha           alpha in f(varphi) = alpha varphi^2
//   Lambda          Lambda in f(varphi) = alpha varphi^2
//   sigma           sigma in V(varphi) = Lambda^4 varphi^(-sigma)
//   Hcon          Hubble function in the code unit
//   fourpiG         1/2 Mpl^2 in the code unit
//   a               scale factor (a)
//   dx              dx lattice resolution
//   dtau            dtau time step

// Returns:
//
//////////////////////////

template <class FieldType>
void prepareFTsource_Quintessence(Field<FieldType> & phi, Field<FieldType> & Tij, Field<FieldType> & Sij, Field<FieldType> & chi, Field<FieldType> & pi, Field<FieldType> & V_pi, const double varphi_bg, const double varphi_prime_bg, const double alpha, const double Lambda, const double sigma, const double Hcon, const double fourpiG, const double a, const double dx, const double dtau)
{
	Site x(phi.lattice());
  double f_varphi = alpha * varphi_bg * varphi_bg;
  double f_prime_varphi = 2. * alpha * varphi_bg;
  double f_ddprime_varphi = 2. * alpha;
  double V_varphi = Lambda * Lambda * Lambda * Lambda * pow(varphi_bg, -sigma) ;
  double V_prime_varphi = - sigma * Lambda * Lambda * Lambda * Lambda * pow(varphi_bg, -sigma-1.);

  // Coefficients:
  double beta_coeff = -3. * ( 1. + f_varphi + .5 * varphi_prime_bg * f_prime_varphi) ;
  #ifndef PHINONLINEAR
  COUT << COLORTEXT_RED << " ERROR" << COLORTEXT_RESET << ": You asked for LINEAR PHI while the Equintessence equations are written up using the PHINONLINEAR!" << endl;
  parallel.abortForce();
  #endif

  #ifndef ORIGINALMETRIC
  COUT << COLORTEXT_RED << " ERROR" << COLORTEXT_RESET << ": You asked for expenential metric while the Equintessence equations are written up using the original metric!" << endl;
  parallel.abortForce();
  #endif

  #ifdef PHINONLINEAR
  #ifdef ORIGINALMETRIC
	for (x.first(); x.test(); x.next())
	{
    // 0-0-component:
    Sij(x, 0, 0) =  2. * fourpiG * Tij(x, 0, 0) / a / beta_coeff;
		Sij(x, 0, 0) += (2. * chi(x) - 4. * phi(x)) * (phi(x-0) + phi(x+0) - 2. * phi(x)) / dx / dx;
		Sij(x, 0, 0) += -.5 * (phi(x+0) - phi(x-0)) * (phi(x+0) - phi(x-0)) / dx / dx;
    Sij(x, 0, 0) += .25 * (f_prime_varphi / beta_coeff) * (phi(x+0) - phi(x-0)) * (pi(x+0) - pi(x-0)) / dx / dx;
    Sij(x, 0, 0) += .25 * ((2.*fourpiG + f_ddprime_varphi)/beta_coeff) * (pi(x+0) - pi(x-0)) * (pi(x+0) - pi(x-0)) / dx / dx;
    Sij(x, 0, 0) += ((f_prime_varphi + f_ddprime_varphi)/beta_coeff) * (pi(x-0) + pi(x+0) - 2. * pi(x)) / dx / dx;
    Sij(x, 0, 0) *= (dx * dx);

		// 1-1-component:
    Sij(x, 1, 1) =  2. * fourpiG * Tij(x, 1, 1) / a / beta_coeff;
    Sij(x, 1, 1) += (2. * chi(x) - 4. * phi(x)) * (phi(x-1) + phi(x+1) - 2. * phi(x)) / dx / dx;
    Sij(x, 1, 1) += -.5 * (phi(x+1) - phi(x-1)) * (phi(x+1) - phi(x-1)) / dx / dx;
    Sij(x, 1, 1) += .25 * (f_prime_varphi / beta_coeff) * (phi(x+1) - phi(x-1)) * (pi(x+1) - pi(x-1)) / dx / dx;
    Sij(x, 1, 1) += .25 * ((2.*fourpiG + f_ddprime_varphi)/beta_coeff) * (pi(x+1) - pi(x-1)) * (pi(x+1) - pi(x-1)) / dx / dx;
    Sij(x, 1, 1) += ((f_prime_varphi + f_ddprime_varphi)/beta_coeff) * (pi(x-1) + pi(x+1) - 2. * pi(x)) / dx / dx;
    Sij(x, 1, 1) *= dx * dx;

  	// 2-2-component:
    Sij(x, 2, 2) =  2. * fourpiG * Tij(x, 2, 2) / a / beta_coeff;
    Sij(x, 2, 2) += (2. * chi(x) - 4. * phi(x)) * (phi(x-2) + phi(x+2) - 2. * phi(x)) / dx / dx;
    Sij(x, 2, 2) += -.5 * (phi(x+2) - phi(x-2)) * (phi(x+2) - phi(x-2)) / dx / dx;
    Sij(x, 2, 2) += .25 * (f_prime_varphi / beta_coeff) * (phi(x+2) - phi(x-2)) * (pi(x+2) - pi(x-2)) / dx / dx;
    Sij(x, 2, 2) += .25 * ((2.*fourpiG + f_ddprime_varphi)/beta_coeff) * (pi(x+2) - pi(x-2)) * (pi(x+2) - pi(x-2)) / dx / dx;
    Sij(x, 2, 2) += ((f_prime_varphi + f_ddprime_varphi)/beta_coeff) * (pi(x-2) + pi(x+2) - 2. * pi(x)) / dx / dx;
    Sij(x, 2, 2) *= dx * dx;

		// 0-1-component:
		Sij(x, 0, 1) =   2. * fourpiG * Tij(x, 0, 1) / a / beta_coeff;
		Sij(x, 0, 1) += -2. * (phi(x+0+1) - phi(x+1) + phi(x+0) - phi(x)) * (phi(x+0+1) - phi(x+0) + phi(x+1) - phi(x)) / (2. * dx) / (2. * dx);
		Sij(x, 0, 1) +=  2. * (chi(x+0+1) + chi(x+0) + chi(x+1) + chi(x)) * (phi(x+0+1) - phi(x+1) - phi(x+0) + phi(x)) / dx / dx / 4.;
    Sij(x, 0, 1) +=  2. * (chi(x+0+1) + chi(x+0) + chi(x+1) + chi(x)) * (phi(x+0+1) - phi(x+1) - phi(x+0) + phi(x)) / dx / dx / 4.;
    Sij(x, 0, 1) += -4. * (phi(x+0+1) + phi(x+0) + phi(x+1) + phi(x)) * (phi(x+0+1) - phi(x+1) - phi(x+0) + phi(x)) / dx / dx / 4.;
    Sij(x, 0, 1) +=  (f_prime_varphi/beta_coeff) * (pi(x+0+1) - pi(x+1) + pi(x+0) - pi(x)) * (phi(x+0+1) - phi(x+0) + phi(x+1) - phi(x)) / (2. * dx) / (2. * dx);
    Sij(x, 0, 1) +=  (f_prime_varphi/beta_coeff) * (phi(x+0+1) - phi(x+1) + phi(x+0) - phi(x)) * (pi(x+0+1) - pi(x+0) + pi(x+1) - pi(x)) / (2. * dx) / (2. * dx);
    Sij(x, 0, 1) +=  ((2.*fourpiG + f_ddprime_varphi)/beta_coeff) * (pi(x+0+1) - pi(x+1) + pi(x+0) - pi(x)) * (pi(x+0+1) - pi(x+0) + pi(x+1) - pi(x)) / (2. * dx) / (2. * dx);
    Sij(x, 0, 1) +=  ((f_prime_varphi + f_ddprime_varphi)/beta_coeff) * (pi(x+0+1) - pi(x+1) - pi(x+0) + pi(x)) / dx / dx;
    Sij(x, 0, 1) *= dx * dx;

		// 0-2-component:
		Sij(x, 0, 2) =   2. * fourpiG * Tij(x, 0, 2) / a / beta_coeff;
		Sij(x, 0, 2) += -2. * (phi(x+0+2) - phi(x+2) + phi(x+0) - phi(x)) * (phi(x+0+2) - phi(x+0) + phi(x+2) - phi(x)) / (2. * dx) / (2. * dx);
		Sij(x, 0, 2) +=  2. * (chi(x+0+2) + chi(x+0) + chi(x+2) + chi(x)) * (phi(x+0+2) - phi(x+2) - phi(x+0) + phi(x)) / dx / dx / 4.;
    Sij(x, 0, 2) +=  2. * (chi(x+0+2) + chi(x+0) + chi(x+2) + chi(x)) * (phi(x+0+2) - phi(x+2) - phi(x+0) + phi(x)) / dx / dx / 4.;
    Sij(x, 0, 2) += -4. * (phi(x+0+2) + phi(x+0) + phi(x+2) + phi(x)) * (phi(x+0+2) - phi(x+2) - phi(x+0) + phi(x)) / dx / dx / 4.;
    Sij(x, 0, 2) +=  (f_prime_varphi/beta_coeff) * (pi(x+0+2) - pi(x+2) + pi(x+0) - pi(x)) * (phi(x+0+2) - phi(x+0) + phi(x+2) - phi(x)) / (2. * dx) / (2. * dx);
    Sij(x, 0, 2) +=  (f_prime_varphi/beta_coeff) * (phi(x+0+2) - phi(x+2) + phi(x+0) - phi(x)) * (pi(x+0+2) - pi(x+0) + pi(x+2) - pi(x)) / (2. * dx) / (2. * dx);
    Sij(x, 0, 2) +=  ((2.*fourpiG + f_ddprime_varphi)/beta_coeff) * (pi(x+0+2) - pi(x+2) + pi(x+0) - pi(x)) * (pi(x+0+2) - pi(x+0) + pi(x+2) - pi(x)) / (2. * dx) / (2. * dx);
    Sij(x, 0, 2) +=  ((f_prime_varphi + f_ddprime_varphi)/beta_coeff) * (pi(x+0+2) - pi(x+2) - pi(x+0) + pi(x)) / dx / dx;
    Sij(x, 0, 2) *= dx * dx;

  // 1-2-component:
	Sij(x, 1, 2) =   2. * fourpiG * Tij(x, 1, 2) / a / beta_coeff;
	Sij(x, 1, 2) += -2. * (phi(x+1+2) - phi(x+2) + phi(x+1) - phi(x)) * (phi(x+1+2) - phi(x+1) + phi(x+2) - phi(x)) / (2. * dx) / (2. * dx);
	Sij(x, 1, 2) +=  2. * (chi(x+1+2) + chi(x+1) + chi(x+2) + chi(x)) * (phi(x+1+2) - phi(x+2) - phi(x+1) + phi(x)) / dx / dx /4.;
  Sij(x, 1, 2) +=  2. * (chi(x+1+2) + chi(x+1) + chi(x+2) + chi(x)) * (phi(x+1+2) - phi(x+2) - phi(x+1) + phi(x)) / dx / dx /4.;
  Sij(x, 1, 2) += -4. * (phi(x+1+2) + phi(x+1) + phi(x+2) + phi(x)) * (phi(x+1+2) - phi(x+2) - phi(x+1) + phi(x)) / dx / dx /4.;
  Sij(x, 1, 2) +=  (f_prime_varphi/beta_coeff) * (pi(x+1+2) - pi(x+2) + pi(x+1) - pi(x)) * (phi(x+1+2) - phi(x+1) + phi(x+2) - phi(x)) / (2. * dx) / (2. * dx);
  Sij(x, 1, 2) +=  (f_prime_varphi/beta_coeff) * (phi(x+1+2) - phi(x+2) + phi(x+1) - phi(x)) * (pi(x+1+2) - pi(x+1) + pi(x+2) - pi(x)) / (2. * dx) / (2. * dx);
  Sij(x, 1, 2) +=  ((2.*fourpiG + f_ddprime_varphi)/beta_coeff) * (pi(x+1+2) - pi(x+2) + pi(x+1) - pi(x)) * (pi(x+1+2) - pi(x+1) + pi(x+2) - pi(x)) / (2. * dx) / (2. * dx);
  Sij(x, 1, 2) +=  ((f_prime_varphi + f_ddprime_varphi)/beta_coeff) * (pi(x+1+2) - pi(x+2) - pi(x+1) + pi(x)) / dx / dx;

  Sij(x, 1, 2) *= dx * dx;
	}
  #endif
  #endif
}

//////////////////////////
// prepareFTsource (2)
//////////////////////////
// Description:
//   construction of real-space source field for Fourier-based solvers
//
// Arguments:
//   phi        reference to field configuration (first Bardeen potential)
//   chi        reference to field configuration (difference between Bardeen potentials, phi-psi)
//   source     reference to fully dressed source field (rescaled by a^3)
//   bgmodel    background model of the source (rescaled by a^3) to be subtracted
//   result     reference to allocated field which will contain the result (may be identical to source)
//   QUINTESSENCE
//   pi              reference to scalar field in extended quintessence
//   V_pi            reference to time derivative (conformal time) of the scalar field in e-quintessence
//   varphi_bg       reference to backgroudn value of the scalar field in extended quintessence
//   varphi_prime_bg reference to conformal time derivative of varphi extended quintessence
//   alpha           alpha in f(varphi) = alpha varphi^2
//   Lambda          Lambda in f(varphi) = alpha varphi^2
//   sigma           sigma in V(varphi) = Lambda^4 varphi^(-sigma)
//   Hcon          Hubble function in the code unit
//   fourpiG         1/2 Mpl^2 in the code unit
//   a               scale factor (a)
//   dx              dx lattice resolution
//   dtau            dtau time step
// Returns:
//
//////////////////////////

template <class FieldType>
void prepareFTsource_Quintessence(Field<FieldType> & phi, Field<FieldType> & chi, Field<FieldType> & source, const FieldType bgmodel, Field<FieldType> & result, Field<FieldType> & pi, Field<FieldType> & V_pi, const double varphi_bg, const double varphi_prime_bg, const double alpha, const double Lambda, const double sigma, const double Hcon, const double fourpiG, const double a, const double dx, const double dtau, int non_linearity)
{
	Site x(phi.lattice());

  double f_varphi = alpha * varphi_bg * varphi_bg;
  double f_prime_varphi = 2. * alpha * varphi_bg;
  double f_ddprime_varphi = 2. * alpha;
  double V_varphi = Lambda * Lambda * Lambda * Lambda * pow(varphi_bg, -sigma) ;
  double V_prime_varphi = - sigma * Lambda * Lambda * Lambda * Lambda * pow(varphi_bg, -sigma-1.);

  // Coefficients:
  double alpha_coeff;
  double lambda_coeff = 1.5 * (1. + f_varphi);
  double beta_coeff = - 3. * ( Hcon * (1. + f_varphi) + .5 * varphi_prime_bg * f_prime_varphi) ;
  double eta_coeff = 1.5 * Hcon * f_prime_varphi - varphi_prime_bg * fourpiG ;
  double gamma_coeff = 3. * ( Hcon * Hcon * (1. + f_varphi) - varphi_prime_bg * varphi_prime_bg * fourpiG/3. + Hcon * varphi_prime_bg * f_prime_varphi) ;// 4 \pi G = 1/(2 M_pl^2)
  double mu_coeff = -.25 * (1. + 2. * f_ddprime_varphi);
  double nu_coeff;
  double sigma_coeff = .5 * f_prime_varphi;
  double zeta_coeff = .5 * (3. * Hcon * (Hcon * f_prime_varphi + varphi_prime_bg * f_ddprime_varphi) - a * a *  V_prime_varphi);

  // Terms:
  double term1, term2 = 0., term3, term4, term5, term6 = 0., term7;

  #ifndef PHINONLINEAR
  COUT << COLORTEXT_RED << " ERROR" << COLORTEXT_RESET << ": You asked for LINEAR PHI while the Equintessence equations are written up using the PHINONLINEAR!" << endl;
  parallel.abortForce();
  #endif

  #ifndef ORIGINALMETRIC
  COUT << COLORTEXT_RED << " ERROR" << COLORTEXT_RESET << ": You asked for expenential metric while the Equintessence equations are written up using the original metric!" << endl;
  parallel.abortForce();
  #endif

  #ifdef PHINONLINEAR
  #ifdef ORIGINALMETRIC

	for (x.first(); x.test(); x.next())
	{
    alpha_coeff = 6. * phi(x) * (1. + f_varphi) + f_varphi + f_prime_varphi * pi(x); // Alpha should be kept at full order because  (fourpiG/a) * (source(x) - bgmodel)/(1.+alpha_coeff)
    nu_coeff = -.5 * f_prime_varphi;

    if (non_linearity == 1)
    {
      nu_coeff += -2. * f_prime_varphi * phi(x) - .5 * f_ddprime_varphi * pi(x) ;

      term2 = lambda_coeff * (phi(x-0) - phi(x+0)) * (phi(x-0) - phi(x+0));
      term2 += lambda_coeff * (phi(x-1) - phi(x+1)) * (phi(x-1) - phi(x+1));
      term2 += lambda_coeff * (phi(x-2) - phi(x+2)) * (phi(x-2) - phi(x+2));
      term2 /= 4. * dx * dx;

      term6 = mu_coeff * (pi(x-0) - pi(x+0)) * (pi(x-0) - pi(x+0));
      term6 += mu_coeff * (pi(x-1) - pi(x+1)) * (pi(x-1) - pi(x+1));
      term6 += mu_coeff * (pi(x-2) - pi(x+2)) * (pi(x-2) - pi(x+2));
      term6 /= 4. * dx * dx;

      term7 = sigma_coeff * (phi(x-0) - phi(x+0)) * (pi(x-0) - pi(x+0));
      term7 += sigma_coeff * (phi(x-1) - phi(x+1)) * (pi(x-1) - pi(x+1));
      term7 += sigma_coeff * (phi(x-2) - phi(x+2)) * (pi(x-2) - pi(x+2));
      term7 /= 4. * dx * dx;
    }

    // Terms:
    term1 = gamma_coeff * (chi(x) - phi(x));
    term3 = eta_coeff * V_pi(x);
    term4 = zeta_coeff * pi(x);
    term5 = nu_coeff * (pi(x-0) + pi(x+0) - 2. * pi(x));
    term5 += nu_coeff * (pi(x-1) + pi(x+1) - 2. * pi(x));
    term5 += nu_coeff * (pi(x-2) + pi(x+2) - 2. * pi(x));
    term5 /= dx * dx;

    result(x) = (1. - 4. * phi(x)) * (fourpiG/a) * (source(x) - bgmodel);
    result(x) += beta_coeff * phi(x) /dtau;
    result(x) +=  - (term1 + term2 + term3 + term4 + term5 + term6 + term7);
    result(x) *= dx * dx/(1.+alpha_coeff); // It's multiplied to dx^2 as it is done in geovlution GR case!
	}
  #endif
  #endif
}



#ifdef BACKREACTION_TEST
  //Checking field
    template <class FieldType>
    void check_field(Field<FieldType> & field, double constant , string field_name, long n3, string message = "")
    {
    Site x(field.lattice());
    ios::fmtflags f(cout.flags());
    double max = 0., hom = 0., sum = 0., temp;
    for(x.first(); x.test(); x.next())
    {
    temp = field(x);
    hom += temp;
    sum += fabs(temp);
    if(fabs(temp) >= max)
    {
    max = fabs(temp);
    }
    }
    parallel.max(max);
    parallel.sum(sum);
    parallel.sum(hom);
    sum /= n3;
    hom /= n3;
    COUT << scientific << setprecision(6);
    COUT << message
    << setw(17) << field_name
    << " Max = " << setw(9) << max * constant
    << " hom = " << setw(9) << hom * constant
    << endl;
    cout.flags(f);
    }


    // Printing out the average
    template <class FieldType>
    double average(Field<FieldType> & field, double constant , long n3)
    {
    Site x(field.lattice());
    ios::fmtflags f(cout.flags());
    double max = 0., hom = 0., sum = 0., temp;
    for(x.first(); x.test(); x.next())
    {
    temp = field(x);
    hom += temp;
    sum += fabs(temp);
    if(fabs(temp) >= max)
    {
    max = fabs(temp);
    }
    }
    parallel.max(max);
    parallel.sum(sum);
    parallel.sum(hom);
    sum /= n3;
    hom /= n3;
    return hom * constant;
    }

#endif


//////////////////////////
// Update quintessence field (pi)
//////////////////////////
// Description:
//   Updating quintessence field pi = deltaphi field based on equation obtained by energy momentum conservation
//
// Arguments:
//   pi              reference to scalar field in extended quintessence
//   V_pi            reference to time derivative (conformal time) of the scalar field //   dtau            time step
// Returns:
//
//////////////////////////
template <class FieldType>
void update_pi( double dtau, Field<FieldType> & pi, Field<FieldType> & V_pi)
{
  Site x(pi.lattice());
  for (x.first(); x.test(); x.next())
    {
      pi(x) = pi(x)  + dtau * V_pi(x); //  pi(n+1)
    }
}

//////////////////////////
// Update quintessence velocity field (V_pi)
//////////////////////////
// Description:
//   Updating quintessence field V_pi = deltaphi' field based on equation obtained by Klein-Gordon equation
//
// Arguments:
//   phi        reference to field configuration (first Bardeen potential)
//   chi        reference to field configuration (difference between Bardeen potentials, phi-psi)
//   QUINTESSENCE
//   pi              reference to scalar field in extended quintessence
//   V_pi            reference to time derivative (conformal time) of the scalar field in e-quintessence
//   varphi_bg       reference to backgroudn value of the scalar field in extended quintessence
//   varphi_prime_bg reference to conformal time derivative of varphi extended quintessence
//   alpha           alpha in f(varphi) = alpha varphi^2
//   Lambda          Lambda in f(varphi) = alpha varphi^2
//   sigma           sigma in V(varphi) = Lambda^4 varphi^(-sigma)
//   Hcon            conformal Hubble function in the code unit
//   H_prime         time derivative of conformal Hubble function in the code unit
//   fourpiG          1/(2*Mpl^2) in the code unit
//   a               scale factor (a)
//   dx              dx lattice resolution
//   dtau            dtau time step
//   non_linearity   refered to whether the system qoes non-linear: 0-> linear solver ; 1-> non-linear solver
// Returns:
//
//////////////////////////

template <class FieldType>
void update_V_pi(Field<FieldType> & phi, Field<FieldType> & phi_old, Field<FieldType> & chi, Field<FieldType> & chi_old, Field<FieldType> & pi, Field<FieldType> & V_pi, Field<FieldType> & TiimT00, Field<FieldType> & source, const double varphi_bg, const double varphi_prime_bg, const double alpha, const double Lambda, const double sigma, const double Hcon, const double H_prime, const double fourpiG, const double a, const double dx, const double dtau, int non_linearity)
  {
    double phi_prime, chi_prime, pi_prime_prime;
    double f_varphi = alpha * varphi_bg * varphi_bg;
    double f_prime_varphi = 2. * alpha * varphi_bg;
    double f_ddprime_varphi = 2. * alpha;
    double f_dddprime_varphi = 0.;
    double V_varphi = Lambda * Lambda * Lambda * Lambda * pow(varphi_bg, -sigma);
    double V_prime_varphi = - sigma * Lambda * Lambda * Lambda * Lambda * pow(varphi_bg, -sigma-1.);
    double V_ddprime_varphi = sigma * (sigma + 1.0) * Lambda * Lambda * Lambda * Lambda * pow(varphi_bg, -sigma-2.);
    double M_pl2 = .5/fourpiG;
    double gamma = (2./3.) * (1. + f_varphi + 3. *f_prime_varphi * f_prime_varphi/2.);
    double coeff_C = -(2. * Hcon + 2. * (varphi_prime_bg * f_prime_varphi/ (3. * gamma)) * (1.+3. * f_ddprime_varphi) );
    double coeff_D = - (1./gamma) * ( (Hcon*Hcon+H_prime) * (f_prime_varphi*f_prime_varphi * (2.+3. *f_ddprime_varphi) - 2.* f_ddprime_varphi * (1.+f_varphi))+ (a*a/3.) * (2.*V_ddprime_varphi * (1.+f_varphi)  -  f_prime_varphi*V_prime_varphi * (4.+3. *f_ddprime_varphi) )  +  varphi_prime_bg * varphi_prime_bg * f_prime_varphi * f_dddprime_varphi);
    double coeff_E = 1.;
    double coeff_F =0.;
    if (non_linearity == 1) coeff_F = (f_prime_varphi/(3.*gamma)) * (1. + 3.*f_ddprime_varphi);
    double coeff_G = (6.*f_prime_varphi*(Hcon*Hcon+H_prime) - 2.*a*a*V_prime_varphi + (1.+3.*f_ddprime_varphi) * (2.*varphi_prime_bg*varphi_prime_bg*f_prime_varphi)/(3.0*gamma));
    double coeff_H = 0.;
    double Laplacian_pi, Laplacian_phi = 0., Gradpi_Gradpi=0.;
    Site x(phi.lattice());

    for (x.first(); x.test(); x.next())
      {
        phi_prime = (phi(x)-phi_old(x))/dtau;
        chi_prime = (chi(x)-chi_old(x))/dtau;
        Laplacian_pi= pi(x-0) + pi(x+0) - 2. * pi(x);
        Laplacian_pi+=pi(x+1) + pi(x-1) - 2. * pi(x);
        Laplacian_pi+=pi(x+2) + pi(x-2) - 2. * pi(x);
        Laplacian_pi/= dx*dx;

        if (non_linearity == 1)
          {
            coeff_E = (1. + 2. * phi(x) + (2.0 * f_prime_varphi * f_prime_varphi/3./gamma) * (phi(x) - chi(x)) + (f_prime_varphi * f_ddprime_varphi/gamma) * pi(x));

            coeff_H = - (2./(3.*gamma)) * (2.*f_prime_varphi*(1.+f_varphi) * (phi(x) - chi(x)) + (f_prime_varphi*f_prime_varphi - f_ddprime_varphi * (1.+f_varphi)) * pi(x));

            Laplacian_phi= phi(x-0) + phi(x+0) - 2. * phi(x);
            Laplacian_phi+=phi(x+1) + phi(x-1) - 2. * phi(x);
            Laplacian_phi+=phi(x+2) + phi(x-2) - 2. * phi(x);
            Laplacian_phi/= dx * dx;

            Gradpi_Gradpi= .25 * (pi(x + 0)  - pi(x - 0)) * (pi(x + 0) - pi(x - 0));
            Gradpi_Gradpi+=.25 * (pi(x + 1)  - pi(x - 1)) * (pi(x + 1) - pi(x - 1));
            Gradpi_Gradpi+=.25 * (pi(x + 2)  - pi(x - 2)) * (pi(x + 2) - pi(x - 2));
            Gradpi_Gradpi/= dx * dx;
          }
        // TEST:
        pi_prime_prime =  + (Hcon * V_pi(x) + Hcon * Hcon * pi(x)) + 1 * Laplacian_pi + TiimT00(x); //+ 1 * Laplacian_pi ); //+ coeff_D * pi(x));
        V_pi(x) = V_pi(x) + pi_prime_prime * dtau;
      // V_pi(x) = (1. + coeff_C * dtau/2.) * V_pi(x)
      // + (dtau/2.) * (coeff_D * pi(x)
      //         + coeff_E * Laplacian_pi
      //         + coeff_F * Gradpi_Gradpi
      //         + varphi_prime_bg *  (4. * phi_prime - chi_prime)
      //         + coeff_G * (phi(x) - chi(x))
      //         + coeff_H * Laplacian_phi
      // + a * a * f_prime_varphi/gamma/3. * TiimT00(x));
      // V_pi(x) /= (1. - coeff_C * dtau/2.);
  }
}

#ifdef FFT3D
//////////////////////////
// projectFTscalar
//////////////////////////
// This equation do not change for the extended-quintessence case.
// Description:
//   projection of the Fourier image of a tensor field on the trace-free
//   longitudinal (scalar) component
//
// Arguments:
//   SijFT      reference to the Fourier image of the input tensor field
//   chiFT      reference to allocated field which will contain the Fourier
//              image of the trace-free longitudinal (scalar) component
//
// Returns:
//
//////////////////////////

void projectFTscalar(Field<Cplx> & SijFT, Field<Cplx> & chiFT, const int add = 0)
{
	const int linesize = chiFT.lattice().size(1);
	int i;
	Real * gridk2;
	Cplx * kshift;
	rKSite k(chiFT.lattice());

	gridk2 = (Real *) malloc(linesize * sizeof(Real));
	kshift = (Cplx *) malloc(linesize * sizeof(Cplx));

	for (i = 0; i < linesize; i++)
	{
		gridk2[i] = 2. * (Real) linesize * sin(M_PI * (Real) i / (Real) linesize);
		kshift[i] = gridk2[i] * Cplx(cos(M_PI * (Real) i / (Real) linesize), -sin(M_PI * (Real) i / (Real) linesize));
		gridk2[i] *= gridk2[i];
	}

	k.first();
	if (k.coord(0) == 0 && k.coord(1) == 0 && k.coord(2) == 0)
	{
		chiFT(k) = Cplx(0.,0.);
		k.next();
	}

	if (add)
	{
		for (; k.test(); k.next())
		{
			chiFT(k) += ((gridk2[k.coord(1)] + gridk2[k.coord(2)] - 2. * gridk2[k.coord(0)]) * SijFT(k, 0, 0) +
						(gridk2[k.coord(0)] + gridk2[k.coord(2)] - 2. * gridk2[k.coord(1)]) * SijFT(k, 1, 1) +
						(gridk2[k.coord(0)] + gridk2[k.coord(1)] - 2. * gridk2[k.coord(2)]) * SijFT(k, 2, 2) -
						6. * kshift[k.coord(0)] * kshift[k.coord(1)] * SijFT(k, 0, 1) -
						6. * kshift[k.coord(0)] * kshift[k.coord(2)] * SijFT(k, 0, 2) -
						6. * kshift[k.coord(1)] * kshift[k.coord(2)] * SijFT(k, 1, 2)) /
						(2. * (gridk2[k.coord(0)] + gridk2[k.coord(1)] + gridk2[k.coord(2)]) * (gridk2[k.coord(0)] + gridk2[k.coord(1)] + gridk2[k.coord(2)]) * linesize);
		}
	}
	else
	{
		for (; k.test(); k.next())
		{
			chiFT(k) = ((gridk2[k.coord(1)] + gridk2[k.coord(2)] - 2. * gridk2[k.coord(0)]) * SijFT(k, 0, 0) +
						(gridk2[k.coord(0)] + gridk2[k.coord(2)] - 2. * gridk2[k.coord(1)]) * SijFT(k, 1, 1) +
						(gridk2[k.coord(0)] + gridk2[k.coord(1)] - 2. * gridk2[k.coord(2)]) * SijFT(k, 2, 2) -
						6. * kshift[k.coord(0)] * kshift[k.coord(1)] * SijFT(k, 0, 1) -
						6. * kshift[k.coord(0)] * kshift[k.coord(2)] * SijFT(k, 0, 2) -
						6. * kshift[k.coord(1)] * kshift[k.coord(2)] * SijFT(k, 1, 2)) /
						(2. * (gridk2[k.coord(0)] + gridk2[k.coord(1)] + gridk2[k.coord(2)]) * (gridk2[k.coord(0)] + gridk2[k.coord(1)] + gridk2[k.coord(2)]) * linesize);
		}
	}

	free(gridk2);
	free(kshift);
}


//////////////////////////
// evolveFTvector
//////////////////////////
// Description:
//   projects the Fourier image of a tensor field on the spin-1 component
//   used as a source for the evolution of the vector perturbation
//
// Arguments:
//   SijFT      reference to the Fourier image of the input tensor field
//   BiFT       reference to the Fourier image of the vector perturbation
//   a2dtau     conformal time step times scale factor squared (a^2 * dtau)
//
// Returns:
//
//////////////////////////

void evolveFTvector(Field<Cplx> & SijFT, Field<Cplx> & BiFT, const Real a2dtau)
{
	const int linesize = BiFT.lattice().size(1);
	int i;
	Real * gridk2;
	Cplx * kshift;
	rKSite k(BiFT.lattice());
	Real k4;

	gridk2 = (Real *) malloc(linesize * sizeof(Real));
	kshift = (Cplx *) malloc(linesize * sizeof(Cplx));

	for (i = 0; i < linesize; i++)
	{
		gridk2[i] = 2. * (Real) linesize * sin(M_PI * (Real) i / (Real) linesize);
		kshift[i] = gridk2[i] * Cplx(cos(M_PI * (Real) i / (Real) linesize), -sin(M_PI * (Real) i / (Real) linesize));
		gridk2[i] *= gridk2[i];
	}

	k.first();
	if (k.coord(0) == 0 && k.coord(1) == 0 && k.coord(2) == 0)
	{
		BiFT(k, 0) = Cplx(0.,0.);
		BiFT(k, 1) = Cplx(0.,0.);
		BiFT(k, 2) = Cplx(0.,0.);
		k.next();
	}

	for (; k.test(); k.next())
	{
		k4 = gridk2[k.coord(0)] + gridk2[k.coord(1)] + gridk2[k.coord(2)];
		k4 *= k4;

		BiFT(k, 0) += Cplx(0.,-2.*a2dtau/k4) * (kshift[k.coord(0)].conj() * ((gridk2[k.coord(1)] + gridk2[k.coord(2)]) * SijFT(k, 0, 0)
				- gridk2[k.coord(1)] * SijFT(k, 1, 1) - gridk2[k.coord(2)] * SijFT(k, 2, 2) - 2. * kshift[k.coord(1)] * kshift[k.coord(2)] * SijFT(k, 1, 2))
				+ (gridk2[k.coord(1)] + gridk2[k.coord(2)] - gridk2[k.coord(0)]) * (kshift[k.coord(1)] * SijFT(k, 0, 1) + kshift[k.coord(2)] * SijFT(k, 0, 2)));
		BiFT(k, 1) += Cplx(0.,-2.*a2dtau/k4) * (kshift[k.coord(1)].conj() * ((gridk2[k.coord(0)] + gridk2[k.coord(2)]) * SijFT(k, 1, 1)
				- gridk2[k.coord(0)] * SijFT(k, 0, 0) - gridk2[k.coord(2)] * SijFT(k, 2, 2) - 2. * kshift[k.coord(0)] * kshift[k.coord(2)] * SijFT(k, 0, 2))
				+ (gridk2[k.coord(0)] + gridk2[k.coord(2)] - gridk2[k.coord(1)]) * (kshift[k.coord(0)] * SijFT(k, 0, 1) + kshift[k.coord(2)] * SijFT(k, 1, 2)));
		BiFT(k, 2) += Cplx(0.,-2.*a2dtau/k4) * (kshift[k.coord(2)].conj() * ((gridk2[k.coord(0)] + gridk2[k.coord(1)]) * SijFT(k, 2, 2)
				- gridk2[k.coord(0)] * SijFT(k, 0, 0) - gridk2[k.coord(1)] * SijFT(k, 1, 1) - 2. * kshift[k.coord(0)] * kshift[k.coord(1)] * SijFT(k, 0, 1))
				+ (gridk2[k.coord(0)] + gridk2[k.coord(1)] - gridk2[k.coord(2)]) * (kshift[k.coord(0)] * SijFT(k, 0, 2) + kshift[k.coord(1)] * SijFT(k, 1, 2)));
	}

	free(gridk2);
	free(kshift);
}


//////////////////////////
// projectFTvector
//////////////////////////
// Description:
//   projects the Fourier image of a vector field on the transverse component
//   and solves the constraint equation for the vector perturbation
//
// Arguments:
//   SiFT       reference to the Fourier image of the input vector field
//   BiFT       reference to the Fourier image of the vector perturbation (can be identical to input)
//   coeff      rescaling coefficient (default 1)
//   modif      modification k^2 -> k^2 + modif (default 0)
//
// Returns:
//
//////////////////////////

void projectFTvector(Field<Cplx> & SiFT, Field<Cplx> & BiFT, const Real coeff = 1., const Real modif = 0.)
{
	const int linesize = BiFT.lattice().size(1);
	int i;
	Real * gridk2;
	Cplx * kshift;
	rKSite k(BiFT.lattice());
	Real k2;
	Cplx tmp(0., 0.);

	gridk2 = (Real *) malloc(linesize * sizeof(Real));
	kshift = (Cplx *) malloc(linesize * sizeof(Cplx));

	for (i = 0; i < linesize; i++)
	{
		gridk2[i] = 2. * (Real) linesize * sin(M_PI * (Real) i / (Real) linesize);
		kshift[i] = gridk2[i] * Cplx(cos(M_PI * (Real) i / (Real) linesize), -sin(M_PI * (Real) i / (Real) linesize));
		gridk2[i] *= gridk2[i];
	}

	k.first();
	if (k.coord(0) == 0 && k.coord(1) == 0 && k.coord(2) == 0)
	{
		BiFT(k, 0) = Cplx(0.,0.);
		BiFT(k, 1) = Cplx(0.,0.);
		BiFT(k, 2) = Cplx(0.,0.);
		k.next();
	}

	for (; k.test(); k.next())
	{
		k2 = gridk2[k.coord(0)] + gridk2[k.coord(1)] + gridk2[k.coord(2)];

		tmp = (kshift[k.coord(0)] * SiFT(k, 0) + kshift[k.coord(1)] * SiFT(k, 1) + kshift[k.coord(2)] * SiFT(k, 2)) / k2;

		BiFT(k, 0) = (SiFT(k, 0) - kshift[k.coord(0)].conj() * tmp) * 4. * coeff / (k2 + modif);
		BiFT(k, 1) = (SiFT(k, 1) - kshift[k.coord(1)].conj() * tmp) * 4. * coeff / (k2 + modif);
		BiFT(k, 2) = (SiFT(k, 2) - kshift[k.coord(2)].conj() * tmp) * 4. * coeff / (k2 + modif);
	}

	free(gridk2);
	free(kshift);
}


//////////////////////////
// projectFTtensor
//////////////////////////
// Note that assuming the quasi static approximation, the results would not change.
// Description:
//   projection of the Fourier image of a tensor field on the transverse
//   trace-free tensor component
//
// Arguments:
//   SijFT      reference to the Fourier image of the input tensor field
//   hijFT      reference to allocated field which will contain the Fourier
//              image of the transverse trace-free tensor component
//
// Returns:
//
//////////////////////////

void projectFTtensor(Field<Cplx> & SijFT, Field<Cplx> & hijFT)
{
	const int linesize = hijFT.lattice().size(1);
	int i;
	Real * gridk2;
	Cplx * kshift;
	rKSite k(hijFT.lattice());
	Cplx SxxFT, SxyFT, SxzFT, SyyFT, SyzFT, SzzFT;
	Real k2, k6;

	gridk2 = (Real *) malloc(linesize * sizeof(Real));
	kshift = (Cplx *) malloc(linesize * sizeof(Cplx));

	for (i = 0; i < linesize; i++)
	{
		gridk2[i] = 2. * (Real) linesize * sin(M_PI * (Real) i / (Real) linesize);
		kshift[i] = gridk2[i] * Cplx(cos(M_PI * (Real) i / (Real) linesize), -sin(M_PI * (Real) i / (Real) linesize));
		gridk2[i] *= gridk2[i];
	}

	k.first();
	if (k.coord(0) == 0 && k.coord(1) == 0 && k.coord(2) == 0)
	{
		for (i = 0; i < hijFT.components(); i++)
			hijFT(k, i) = Cplx(0.,0.);

		k.next();
	}

	for (; k.test(); k.next())
	{
		SxxFT = SijFT(k, 0, 0);
		SxyFT = SijFT(k, 0, 1);
		SxzFT = SijFT(k, 0, 2);
		SyyFT = SijFT(k, 1, 1);
		SyzFT = SijFT(k, 1, 2);
		SzzFT = SijFT(k, 2, 2);

		k2 = gridk2[k.coord(0)] + gridk2[k.coord(1)] + gridk2[k.coord(2)];
		k6 = k2 * k2 * k2 * linesize;

		hijFT(k, 0, 0) = ((gridk2[k.coord(0)] - k2) * ((gridk2[k.coord(0)] - k2) * SxxFT + 2. * kshift[k.coord(0)] * (kshift[k.coord(1)] * SxyFT + kshift[k.coord(2)] * SxzFT))
				+ ((gridk2[k.coord(0)] + k2) * (gridk2[k.coord(1)] + k2) - 2. * k2 * k2) * SyyFT
				+ ((gridk2[k.coord(0)] + k2) * (gridk2[k.coord(2)] + k2) - 2. * k2 * k2) * SzzFT
				+ 2. * (gridk2[k.coord(0)] + k2) * kshift[k.coord(1)] * kshift[k.coord(2)] * SyzFT) / k6;

		hijFT(k, 0, 1) = (2. * (gridk2[k.coord(0)] - k2) * (gridk2[k.coord(1)] - k2) * SxyFT + (gridk2[k.coord(2)] + k2) * kshift[k.coord(0)].conj() * kshift[k.coord(1)].conj() * SzzFT
				+ (gridk2[k.coord(0)] - k2) * kshift[k.coord(1)].conj() * (kshift[k.coord(0)].conj() * SxxFT + 2. * kshift[k.coord(2)] * SxzFT)
				+ (gridk2[k.coord(1)] - k2) * kshift[k.coord(0)].conj() * (kshift[k.coord(1)].conj() * SyyFT + 2. * kshift[k.coord(2)] * SyzFT)) / k6;

		hijFT(k, 0, 2) = (2. * (gridk2[k.coord(0)] - k2) * (gridk2[k.coord(2)] - k2) * SxzFT + (gridk2[k.coord(1)] + k2) * kshift[k.coord(0)].conj() * kshift[k.coord(2)].conj() * SyyFT
				+ (gridk2[k.coord(0)] - k2) * kshift[k.coord(2)].conj() * (kshift[k.coord(0)].conj() * SxxFT + 2. * kshift[k.coord(1)] * SxyFT)
				+ (gridk2[k.coord(2)] - k2) * kshift[k.coord(0)].conj() * (kshift[k.coord(2)].conj() * SzzFT + 2. * kshift[k.coord(1)] * SyzFT)) / k6;

		hijFT(k, 1, 1) = ((gridk2[k.coord(1)] - k2) * ((gridk2[k.coord(1)] - k2) * SyyFT + 2. * kshift[k.coord(1)] * (kshift[k.coord(0)] * SxyFT + kshift[k.coord(2)] * SyzFT))
				+ ((gridk2[k.coord(1)] + k2) * (gridk2[k.coord(0)] + k2) - 2. * k2 * k2) * SxxFT
				+ ((gridk2[k.coord(1)] + k2) * (gridk2[k.coord(2)] + k2) - 2. * k2 * k2) * SzzFT
				+ 2. * (gridk2[k.coord(1)] + k2) * kshift[k.coord(0)] * kshift[k.coord(2)] * SxzFT) / k6;

		hijFT(k, 1, 2) = (2. * (gridk2[k.coord(1)] - k2) * (gridk2[k.coord(2)] - k2) * SyzFT + (gridk2[k.coord(0)] + k2) * kshift[k.coord(1)].conj() * kshift[k.coord(2)].conj() * SxxFT
				+ (gridk2[k.coord(1)] - k2) * kshift[k.coord(2)].conj() * (kshift[k.coord(1)].conj() * SyyFT + 2. * kshift[k.coord(0)] * SxyFT)
				+ (gridk2[k.coord(2)] - k2) * kshift[k.coord(1)].conj() * (kshift[k.coord(2)].conj() * SzzFT + 2. * kshift[k.coord(0)] * SxzFT)) / k6;

		hijFT(k, 2, 2) = ((gridk2[k.coord(2)] - k2) * ((gridk2[k.coord(2)] - k2) * SzzFT + 2. * kshift[k.coord(2)] * (kshift[k.coord(0)] * SxzFT + kshift[k.coord(1)] * SyzFT))
				+ ((gridk2[k.coord(2)] + k2) * (gridk2[k.coord(0)] + k2) - 2. * k2 * k2) * SxxFT
				+ ((gridk2[k.coord(2)] + k2) * (gridk2[k.coord(1)] + k2) - 2. * k2 * k2) * SyyFT
				+ 2. * (gridk2[k.coord(2)] + k2) * kshift[k.coord(0)] * kshift[k.coord(1)] * SxyFT) / k6;
	}

	free(gridk2);
	free(kshift);
}



//////////////////////////
// solveModifiedPoissonFT_quintessence
//////////////////////////
// Description:
//   Modified Poisson solver using the standard Fourier method for quintessence
//
// Arguments:
//   sourceFT   reference to the Fourier image of the source field
//   potFT      reference to the Fourier image of the potential
//   dx2        lattice resolution 'dx' squared.
// in e-quintessence
//   varphi_bg       reference to backgroudn value of the scalar field in extended quintessence
//   varphi_prime_bg reference to conformal time derivative of varphi extended quintessence
//   alpha           alpha in f(varphi) = alpha varphi^2
//   Hcon            conformal Hubble function in the code unit
//   dtau            dtau time step
//   modif      modification k^2 -> k^2 + modif which is -beta(varphi)/(1+f(varphi)) dtau where we can neglect the first order terms )
// Returns:
//
//////////////////////////


void solveModifiedPoissonFT_quintessence(Field<Cplx> & sourceFT, Field<Cplx> & potFT, Real coeff, const double varphi_bg, const double varphi_prime_bg, const double alpha, const double Hcon, const double dtau)
{

  double alpha_phi = alpha * varphi_bg * varphi_bg;
  double f_prime_varphi = 2. * alpha * varphi_bg;
  double beta_coeff = - 3. * ( Hcon * (1. + alpha_phi) + .5 * varphi_prime_bg * f_prime_varphi) ;
  double modif = - beta_coeff/((1. + alpha_phi) * dtau);

	const int linesize = potFT.lattice().size(1);
	int i;
	Real * gridk2;
	Real * sinc;
	rKSite k(potFT.lattice());

	gridk2 = (Real *) malloc(linesize * sizeof(Real));
  coeff /= -((long) linesize * (long) linesize * (long) linesize);

	for (i = 0; i < linesize; i++)
	{
		gridk2[i] = 2. * (Real) linesize * sin(M_PI * (Real) i / (Real) linesize);
		gridk2[i] *= gridk2[i];
	}

	k.first();
  if (k.coord(0) == 0 && k.coord(1) == 0 && k.coord(2) == 0)
  {
    if (modif == 0.)
      potFT(k) = Cplx(0.,0.);
    else
      potFT(k) = Cplx(0.,0.);
    k.next();
  }

	for (; k.test(); k.next())
	{
		potFT(k) = sourceFT(k) * coeff / (gridk2[k.coord(0)] + gridk2[k.coord(1)] +  gridk2[k.coord(2)] + modif);
	}

	free(gridk2);
}


//////////////////////////
// solveModifiedPoissonFT
//////////////////////////
// Description:
//   Modified Poisson solver using the standard Fourier method
//
// Arguments:
//   sourceFT   reference to the Fourier image of the source field
//   potFT      reference to the Fourier image of the potential
//   coeff      coefficient applied to the source ("4 pi G / a")
//   modif      modification k^2 -> k^2 + modif (default 0 gives standard Poisson equation)
//
// Returns:
//
//////////////////////////

void solveModifiedPoissonFT(Field<Cplx> & sourceFT, Field<Cplx> & potFT, Real coeff, const Real modif = 0.)
{
	const int linesize = potFT.lattice().size(1);
	int i;
	Real * gridk2;
	Real * sinc;
	rKSite k(potFT.lattice());

	gridk2 = (Real *) malloc(linesize * sizeof(Real));

	coeff /= -((long) linesize * (long) linesize * (long) linesize);

	for (i = 0; i < linesize; i++)
	{
		gridk2[i] = 2. * (Real) linesize * sin(M_PI * (Real) i / (Real) linesize);
		gridk2[i] *= gridk2[i];
	}

	k.first();
	if (k.coord(0) == 0 && k.coord(1) == 0 && k.coord(2) == 0)
	{
		if (modif == 0.)
			potFT(k) = Cplx(0.,0.);
		else
			potFT(k) = sourceFT(k) * coeff / modif;
		k.next();
	}

	for (; k.test(); k.next())
	{
		potFT(k) = sourceFT(k) * coeff / (gridk2[k.coord(0)] + gridk2[k.coord(1)] + gridk2[k.coord(2)] + modif);
	}

	free(gridk2);
}
#endif


//////////////////////////
// update_q
//////////////////////////
// Description:
//   Update momentum method (arbitrary momentum)
//   Note that vel[3] in the particle structure is used to store q[3] in units
//   of the particle mass, such that as q^2 << m^2 a^2 the meaning of vel[3]
//   is ~ v*a.
//
// Arguments:
//   dtau       time step
//   dx         lattice unit
//   part       pointer to particle structure
//   ref_dist   distance vector to reference point
//   partInfo   global particle properties (unused)
//   fields     array of pointers to fields appearing in geodesic equation
//              fields[0] = phi
//              fields[1] = chi
//              fields[2] = Bi
//   sites      array of sites on the respective lattices
//   nfield     number of fields
//   params     array of additional parameters
//              params[0] = a
//              params[1] = scaling coefficient for Bi
//   outputs    array of reduction variables
//   noutputs   number of reduction variables
//
// Returns: squared velocity of particle after update
//
//////////////////////////

Real update_q(double dtau, double dx, part_simple * part, double * ref_dist, part_simple_info partInfo, Field<Real> ** fields, Site * sites, int nfield, double * params, double * outputs, int noutputs)
{
#define phi (*fields[0])
#define chi (*fields[1])
#define Bi (*fields[2])
#define xphi (sites[0])
#define xchi (sites[1])
#define xB (sites[2])

	Real gradphi[3]={0,0,0};
	Real pgradB[3]={0,0,0};
	Real v2 = (*part).vel[0] * (*part).vel[0] + (*part).vel[1] * (*part).vel[1] + (*part).vel[2] * (*part).vel[2];
	Real e2 = v2 + params[0] * params[0];

	gradphi[0] = (1.-ref_dist[1]) * (1.-ref_dist[2]) * (phi(xphi+0) - phi(xphi));
	gradphi[1] = (1.-ref_dist[0]) * (1.-ref_dist[2]) * (phi(xphi+1) - phi(xphi));
	gradphi[2] = (1.-ref_dist[0]) * (1.-ref_dist[1]) * (phi(xphi+2) - phi(xphi));
	gradphi[0] += ref_dist[1] * (1.-ref_dist[2]) * (phi(xphi+1+0) - phi(xphi+1));
	gradphi[1] += ref_dist[0] * (1.-ref_dist[2]) * (phi(xphi+1+0) - phi(xphi+0));
	gradphi[2] += ref_dist[0] * (1.-ref_dist[1]) * (phi(xphi+2+0) - phi(xphi+0));
	gradphi[0] += (1.-ref_dist[1]) * ref_dist[2] * (phi(xphi+2+0) - phi(xphi+2));
	gradphi[1] += (1.-ref_dist[0]) * ref_dist[2] * (phi(xphi+2+1) - phi(xphi+2));
	gradphi[2] += (1.-ref_dist[0]) * ref_dist[1] * (phi(xphi+2+1) - phi(xphi+1));
	gradphi[0] += ref_dist[1] * ref_dist[2] * (phi(xphi+2+1+0) - phi(xphi+2+1));
	gradphi[1] += ref_dist[0] * ref_dist[2] * (phi(xphi+2+1+0) - phi(xphi+2+0));
	gradphi[2] += ref_dist[0] * ref_dist[1] * (phi(xphi+2+1+0) - phi(xphi+1+0));

	gradphi[0] *= (v2 + e2) / e2;
	gradphi[1] *= (v2 + e2) / e2;
	gradphi[2] *= (v2 + e2) / e2;

	if (nfield>=2 && fields[1] != NULL)
	{
		gradphi[0] -= (1.-ref_dist[1]) * (1.-ref_dist[2]) * (chi(xchi+0) - chi(xchi));
		gradphi[1] -= (1.-ref_dist[0]) * (1.-ref_dist[2]) * (chi(xchi+1) - chi(xchi));
		gradphi[2] -= (1.-ref_dist[0]) * (1.-ref_dist[1]) * (chi(xchi+2) - chi(xchi));
		gradphi[0] -= ref_dist[1] * (1.-ref_dist[2]) * (chi(xchi+1+0) - chi(xchi+1));
		gradphi[1] -= ref_dist[0] * (1.-ref_dist[2]) * (chi(xchi+1+0) - chi(xchi+0));
		gradphi[2] -= ref_dist[0] * (1.-ref_dist[1]) * (chi(xchi+2+0) - chi(xchi+0));
		gradphi[0] -= (1.-ref_dist[1]) * ref_dist[2] * (chi(xchi+2+0) - chi(xchi+2));
		gradphi[1] -= (1.-ref_dist[0]) * ref_dist[2] * (chi(xchi+2+1) - chi(xchi+2));
		gradphi[2] -= (1.-ref_dist[0]) * ref_dist[1] * (chi(xchi+2+1) - chi(xchi+1));
		gradphi[0] -= ref_dist[1] * ref_dist[2] * (chi(xchi+2+1+0) - chi(xchi+2+1));
		gradphi[1] -= ref_dist[0] * ref_dist[2] * (chi(xchi+2+1+0) - chi(xchi+2+0));
		gradphi[2] -= ref_dist[0] * ref_dist[1] * (chi(xchi+2+1+0) - chi(xchi+1+0));
	}

	e2 = sqrt(e2);

	if (nfield>=3 && fields[2] != NULL)
	{
		pgradB[0] = ((1.-ref_dist[2]) * (Bi(xB+0,1) - Bi(xB,1)) + ref_dist[2] * (Bi(xB+2+0,1) - Bi(xB+2,1))) * (*part).vel[1];
		pgradB[0] += ((1.-ref_dist[1]) * (Bi(xB+0,2) - Bi(xB,2)) + ref_dist[1] * (Bi(xB+1+0,2) - Bi(xB+1,2))) * (*part).vel[2];
		pgradB[0] += (1.-ref_dist[1]) * (1.-ref_dist[2]) * ((ref_dist[0]-1.) * Bi(xB-0,0) + (1.-2.*ref_dist[0]) * Bi(xB,0) + ref_dist[0] * Bi(xB+0,0)) * (*part).vel[0];
		pgradB[0] += ref_dist[1] * (1.-ref_dist[2]) * ((ref_dist[0]-1.) * Bi(xB+1-0,0) + (1.-2.*ref_dist[0]) * Bi(xB+1,0) + ref_dist[0] * Bi(xB+1+0,0)) * (*part).vel[0];
		pgradB[0] += (1.-ref_dist[1]) * ref_dist[2] * ((ref_dist[0]-1.) * Bi(xB+2-0,0) + (1.-2.*ref_dist[0]) * Bi(xB+2,0) + ref_dist[0] * Bi(xB+2+0,0)) * (*part).vel[0];
		pgradB[0] += ref_dist[1] * ref_dist[2] * ((ref_dist[0]-1.) * Bi(xB+2+1-0,0) + (1.-2.*ref_dist[0]) * Bi(xB+2+1,0) + ref_dist[0] * Bi(xB+2+1+0,0)) * (*part).vel[0];

		pgradB[1] = ((1.-ref_dist[0]) * (Bi(xB+1,2) - Bi(xB,2)) + ref_dist[0] * (Bi(xB+1+0,2) - Bi(xB+0,2))) * (*part).vel[2];
		pgradB[1] += ((1.-ref_dist[2]) * (Bi(xB+1,0) - Bi(xB,0)) + ref_dist[2] * (Bi(xB+1+2,0) - Bi(xB+2,0))) * (*part).vel[0];
		pgradB[1] += (1.-ref_dist[0]) * (1.-ref_dist[2]) * ((ref_dist[1]-1.) * Bi(xB-1,1) + (1.-2.*ref_dist[1]) * Bi(xB,1) + ref_dist[1] * Bi(xB+1,1)) * (*part).vel[1];
		pgradB[1] += ref_dist[0] * (1.-ref_dist[2]) * ((ref_dist[1]-1.) * Bi(xB+0-1,1) + (1.-2.*ref_dist[1]) * Bi(xB+0,1) + ref_dist[1] * Bi(xB+0+1,1)) * (*part).vel[1];
		pgradB[1] += (1.-ref_dist[0]) * ref_dist[2] * ((ref_dist[1]-1.) * Bi(xB+2-1,1) + (1.-2.*ref_dist[1]) * Bi(xB+2,1) + ref_dist[1] * Bi(xB+2+1,1)) * (*part).vel[1];
		pgradB[1] += ref_dist[0] * ref_dist[2] * ((ref_dist[1]-1.) * Bi(xB+2+0-1,1) + (1.-2.*ref_dist[1]) * Bi(xB+2+0,1) + ref_dist[1] * Bi(xB+2+0+1,1)) * (*part).vel[1];

		pgradB[2] = ((1.-ref_dist[1]) * (Bi(xB+2,0) - Bi(xB,0)) + ref_dist[1] * (Bi(xB+2+1,0) - Bi(xB+1,0))) * (*part).vel[0];
		pgradB[2] += ((1.-ref_dist[0]) * (Bi(xB+2,1) - Bi(xB,1)) + ref_dist[0] * (Bi(xB+2+0,1) - Bi(xB+0,1))) * (*part).vel[1];
		pgradB[2] += (1.-ref_dist[0]) * (1.-ref_dist[1]) * ((ref_dist[2]-1.) * Bi(xB-2,2) + (1.-2.*ref_dist[2]) * Bi(xB,2) + ref_dist[2] * Bi(xB+2,2)) * (*part).vel[2];
		pgradB[2] += ref_dist[0] * (1.-ref_dist[1]) * ((ref_dist[2]-1.) * Bi(xB+0-2,2) + (1.-2.*ref_dist[2]) * Bi(xB+0,2) + ref_dist[2] * Bi(xB+0+2,2)) * (*part).vel[2];
		pgradB[2] += (1.-ref_dist[0]) * ref_dist[1] * ((ref_dist[2]-1.) * Bi(xB+1-2,2) + (1.-2.*ref_dist[2]) * Bi(xB+1,2) + ref_dist[2] * Bi(xB+2+1,2)) * (*part).vel[2];
		pgradB[2] += ref_dist[0] * ref_dist[1] * ((ref_dist[2]-1.) * Bi(xB+1+0-2,2) + (1.-2.*ref_dist[2]) * Bi(xB+1+0,2) + ref_dist[2] * Bi(xB+1+0+2,2)) * (*part).vel[2];

		gradphi[0] += pgradB[0] / params[1] / e2;
		gradphi[1] += pgradB[1] / params[1] / e2;
		gradphi[2] += pgradB[2] / params[1] / e2;
	}

	v2 = 0.;
	for (int i=0;i<3;i++)
	{
		(*part).vel[i] -= dtau * e2 * gradphi[i] / dx;
		v2 += (*part).vel[i] * (*part).vel[i];
	}

	return v2 / params[0] / params[0];

#undef phi
#undef chi
#undef Bi
#undef xphi
#undef xchi
#undef xB
}


//////////////////////////
// update_q_Newton
//////////////////////////
// Description:
//   Update momentum method (Newtonian version)
//   Note that vel[3] in the particle structure is used to store q[3] in units
//   of the particle mass, such that the meaning of vel[3] is v*a.
//
// Arguments:
//   dtau       time step
//   dx         lattice unit
//   part       pointer to particle structure
//   ref_dist   distance vector to reference point
//   partInfo   global particle properties (unused)
//   fields     array of pointers to fields appearing in geodesic equation
//              fields[0] = psi
//              fields[1] = chi
//   sites      array of sites on the respective lattices
//   nfield     number of fields (should be 1)
//   params     array of additional parameters
//              params[0] = a
//   outputs    array of reduction variables
//   noutputs   number of reduction variables
//
// Returns: squared velocity of particle after update
//
//////////////////////////

Real update_q_Newton(double dtau, double dx, part_simple * part, double * ref_dist, part_simple_info partInfo, Field<Real> ** fields, Site * sites, int nfield, double * params, double * outputs, int noutputs)
{
#define psi (*fields[0])
#define xpsi (sites[0])
#define chi (*fields[1])
#define xchi (sites[1])

	Real gradpsi[3]={0,0,0};

	gradpsi[0] = (1.-ref_dist[1]) * (1.-ref_dist[2]) * (psi(xpsi+0) - psi(xpsi));
	gradpsi[1] = (1.-ref_dist[0]) * (1.-ref_dist[2]) * (psi(xpsi+1) - psi(xpsi));
	gradpsi[2] = (1.-ref_dist[0]) * (1.-ref_dist[1]) * (psi(xpsi+2) - psi(xpsi));
	gradpsi[0] += ref_dist[1] * (1.-ref_dist[2]) * (psi(xpsi+1+0) - psi(xpsi+1));
	gradpsi[1] += ref_dist[0] * (1.-ref_dist[2]) * (psi(xpsi+1+0) - psi(xpsi+0));
	gradpsi[2] += ref_dist[0] * (1.-ref_dist[1]) * (psi(xpsi+2+0) - psi(xpsi+0));
	gradpsi[0] += (1.-ref_dist[1]) * ref_dist[2] * (psi(xpsi+2+0) - psi(xpsi+2));
	gradpsi[1] += (1.-ref_dist[0]) * ref_dist[2] * (psi(xpsi+2+1) - psi(xpsi+2));
	gradpsi[2] += (1.-ref_dist[0]) * ref_dist[1] * (psi(xpsi+2+1) - psi(xpsi+1));
	gradpsi[0] += ref_dist[1] * ref_dist[2] * (psi(xpsi+2+1+0) - psi(xpsi+2+1));
	gradpsi[1] += ref_dist[0] * ref_dist[2] * (psi(xpsi+2+1+0) - psi(xpsi+2+0));
	gradpsi[2] += ref_dist[0] * ref_dist[1] * (psi(xpsi+2+1+0) - psi(xpsi+1+0));

	if (nfield>=2 && fields[1] != NULL)
	{
		gradpsi[0] -= (1.-ref_dist[1]) * (1.-ref_dist[2]) * (chi(xchi+0) - chi(xchi));
		gradpsi[1] -= (1.-ref_dist[0]) * (1.-ref_dist[2]) * (chi(xchi+1) - chi(xchi));
		gradpsi[2] -= (1.-ref_dist[0]) * (1.-ref_dist[1]) * (chi(xchi+2) - chi(xchi));
		gradpsi[0] -= ref_dist[1] * (1.-ref_dist[2]) * (chi(xchi+1+0) - chi(xchi+1));
		gradpsi[1] -= ref_dist[0] * (1.-ref_dist[2]) * (chi(xchi+1+0) - chi(xchi+0));
		gradpsi[2] -= ref_dist[0] * (1.-ref_dist[1]) * (chi(xchi+2+0) - chi(xchi+0));
		gradpsi[0] -= (1.-ref_dist[1]) * ref_dist[2] * (chi(xchi+2+0) - chi(xchi+2));
		gradpsi[1] -= (1.-ref_dist[0]) * ref_dist[2] * (chi(xchi+2+1) - chi(xchi+2));
		gradpsi[2] -= (1.-ref_dist[0]) * ref_dist[1] * (chi(xchi+2+1) - chi(xchi+1));
		gradpsi[0] -= ref_dist[1] * ref_dist[2] * (chi(xchi+2+1+0) - chi(xchi+2+1));
		gradpsi[1] -= ref_dist[0] * ref_dist[2] * (chi(xchi+2+1+0) - chi(xchi+2+0));
		gradpsi[2] -= ref_dist[0] * ref_dist[1] * (chi(xchi+2+1+0) - chi(xchi+1+0));
	}

	Real v2 = 0.;
	for (int i=0;i<3;i++)
	{
		(*part).vel[i] -= dtau * params[0] * gradpsi[i] / dx;
		v2 += (*part).vel[i] * (*part).vel[i];
	}

	return v2 / params[0] / params[0];

#undef psi
#undef xpsi
#undef chi
#undef xchi
}


//////////////////////////
// update_pos
//////////////////////////
// Description:
//   Update position method (arbitrary momentum)
//   Note that vel[3] in the particle structure is used to store q[3] in units
//   of the particle mass, such that as q^2 << m^2 a^2 the meaning of vel[3]
//   is ~ v*a.
//
// Arguments:
//   dtau       time step
//   dx         lattice unit
//   part       pointer to particle structure
//   ref_dist   distance vector to reference point
//   partInfo   global particle properties (unused)
//   fields     array of pointers to fields appearing in geodesic equation
//              fields[0] = phi
//              fields[1] = chi
//              fields[2] = Bi
//   sites      array of sites on the respective lattices
//   nfield     number of fields
//   params     array of additional parameters
//              params[0] = a
//              params[1] = scaling coefficient for Bi
//   outputs    array of reduction variables
//   noutputs   number of reduction variables
//
// Returns:
//
//////////////////////////

void update_pos(double dtau, double dx, part_simple * part, double * ref_dist, part_simple_info partInfo, Field<Real> ** fields, Site * sites, int nfield, double * params, double * outputs, int noutputs)
{
	Real v[3];
	Real v2 = (*part).vel[0] * (*part).vel[0] + (*part).vel[1] * (*part).vel[1] + (*part).vel[2] * (*part).vel[2];
	Real e2 = v2 + params[0] * params[0];
	Real phi = 0;
	Real chi = 0;

	if (nfield >= 1)
	{
		phi = (*fields[0])(sites[0]) * (1.-ref_dist[0]) * (1.-ref_dist[1]) * (1.-ref_dist[2]);
		phi += (*fields[0])(sites[0]+0) * ref_dist[0] * (1.-ref_dist[1]) * (1.-ref_dist[2]);
		phi += (*fields[0])(sites[0]+1) * (1.-ref_dist[0]) * ref_dist[1] * (1.-ref_dist[2]);
		phi += (*fields[0])(sites[0]+0+1) * ref_dist[0] * ref_dist[1] * (1.-ref_dist[2]);
		phi += (*fields[0])(sites[0]+2) * (1.-ref_dist[0]) * (1.-ref_dist[1]) * ref_dist[2];
		phi += (*fields[0])(sites[0]+0+2) * ref_dist[0] * (1.-ref_dist[1]) * ref_dist[2];
		phi += (*fields[0])(sites[0]+1+2) * (1.-ref_dist[0]) * ref_dist[1] * ref_dist[2];
		phi += (*fields[0])(sites[0]+0+1+2) * ref_dist[0] * ref_dist[1] * ref_dist[2];
	}

	if (nfield >= 2)
	{
		chi = (*fields[1])(sites[1]) * (1.-ref_dist[0]) * (1.-ref_dist[1]) * (1.-ref_dist[2]);
		chi += (*fields[1])(sites[1]+0) * ref_dist[0] * (1.-ref_dist[1]) * (1.-ref_dist[2]);
		chi += (*fields[1])(sites[1]+1) * (1.-ref_dist[0]) * ref_dist[1] * (1.-ref_dist[2]);
		chi += (*fields[1])(sites[1]+0+1) * ref_dist[0] * ref_dist[1] * (1.-ref_dist[2]);
		chi += (*fields[1])(sites[1]+2) * (1.-ref_dist[0]) * (1.-ref_dist[1]) * ref_dist[2];
		chi += (*fields[1])(sites[1]+0+2) * ref_dist[0] * (1.-ref_dist[1]) * ref_dist[2];
		chi += (*fields[1])(sites[1]+1+2) * (1.-ref_dist[0]) * ref_dist[1] * ref_dist[2];
		chi += (*fields[1])(sites[1]+0+1+2) * ref_dist[0] * ref_dist[1] * ref_dist[2];
	}

	v2 = (1. + (3. - v2 / e2) * phi - chi) / sqrt(e2);

	v[0] = (*part).vel[0] * v2;
	v[1] = (*part).vel[1] * v2;
	v[2] = (*part).vel[2] * v2;

	if (nfield >= 3)
	{
		Real b[3];

		b[0] = (*fields[2])(sites[2], 0) * (1.-ref_dist[1]) * (1.-ref_dist[2]);
		b[1] = (*fields[2])(sites[2], 1) * (1.-ref_dist[0]) * (1.-ref_dist[2]);
		b[2] = (*fields[2])(sites[2], 2) * (1.-ref_dist[0]) * (1.-ref_dist[1]);
		b[1] += (*fields[2])(sites[2]+0, 1) * ref_dist[0] * (1.-ref_dist[2]);
		b[2] += (*fields[2])(sites[2]+0, 2) * ref_dist[0] * (1.-ref_dist[1]);
		b[0] += (*fields[2])(sites[2]+1, 0) * ref_dist[1] * (1.-ref_dist[2]);
		b[2] += (*fields[2])(sites[2]+1, 2) * (1.-ref_dist[0]) * ref_dist[1];
		b[0] += (*fields[2])(sites[2]+2, 0) * (1.-ref_dist[1]) * ref_dist[2];
		b[1] += (*fields[2])(sites[2]+2, 1) * (1.-ref_dist[0]) * ref_dist[2];
		b[1] += (*fields[2])(sites[2]+2+0, 1) * ref_dist[0] * ref_dist[2];
		b[0] += (*fields[2])(sites[2]+2+1, 0) * ref_dist[1] * ref_dist[2];
		b[2] += (*fields[2])(sites[2]+1+0, 2) * ref_dist[0] * ref_dist[1];

		for (int l=0;l<3;l++) (*part).pos[l] += dtau*(v[l] + b[l] / params[1]);
	}
	else
	{
		for (int l=0;l<3;l++) (*part).pos[l] += dtau*v[l];
	}
}


//////////////////////////
// update_pos_Newton
//////////////////////////
// Description:
//   Update position method (Newtonian version)
//   Note that vel[3] in the particle structure is used to store q[3] in units
//   of the particle mass, such that the meaning of vel[3] is v*a.
//
// Arguments:
//   dtau       time step
//   dx         lattice unit (unused)
//   part       pointer to particle structure
//   ref_dist   distance vector to reference point (unused)
//   partInfo   global particle properties (unused)
//   fields     array of pointers to fields appearing in geodesic equation (unused)
//   sites      array of sites on the respective lattices (unused)
//   nfield     number of fields (unused)
//   params     array of additional parameters
//              params[0] = a
//   outputs    array of reduction variables (unused)
//   noutputs   number of reduction variables (unused)
//
// Returns:
//
//////////////////////////

void update_pos_Newton(double dtau, double dx, part_simple * part, double * ref_dist, part_simple_info partInfo, Field<Real> ** fields, Site * sites, int nfield, double * params, double * outputs, int noutputs)
{
	for (int l=0;l<3;l++) (*part).pos[l] += dtau * (*part).vel[l] / params[0];
}


//////////////////////////
// projection_T00_project
//////////////////////////
// Description:
//   Particle-mesh projection for T00, including geometric corrections
//
// Arguments:
//   pcls       pointer to particle handler
//   T00        pointer to target field
//   a          scale factor at projection (needed in order to convert
//              canonical momenta to energies)
//   phi        pointer to Bardeen potential which characterizes the
//              geometric corrections (volume distortion); can be set to
//              NULL which will result in no corrections applied
//   coeff      coefficient applied to the projection operation (default 1)
//
// Returns:
//
//////////////////////////

template<typename part, typename part_info, typename part_dataType>
void projection_T00_project(Particles<part, part_info, part_dataType> * pcls, Field<Real> * T00, double a = 1., Field<Real> * phi = NULL, double coeff = 1.)
{
	if (T00->lattice().halo() == 0)
	{
		cout<< "projection_T00_project: target field needs halo > 0" << endl;
		exit(-1);
	}

	Site xPart(pcls->lattice());
	Site xField(T00->lattice());

	typename std::list<part>::iterator it;

	Real referPos[3];
	Real weightScalarGridUp[3];
	Real weightScalarGridDown[3];
	Real dx = pcls->res();

	double mass = coeff / (dx*dx*dx);
	mass *= *(double*)((char*)pcls->parts_info() + pcls->mass_offset());
	mass /= a;

	Real e = a, f = 0.;
	Real * q;
	size_t offset_q = offsetof(part,vel);

	Real localCube[8]; // XYZ = 000 | 001 | 010 | 011 | 100 | 101 | 110 | 111
	Real localCubePhi[8];

	for (int i=0; i<8; i++) localCubePhi[i] = 0.0;

	for (xPart.first(),xField.first(); xPart.test(); xPart.next(),xField.next())
	{
		if (pcls->field()(xPart).size != 0)
		{
			for(int i=0; i<3; i++) referPos[i] = xPart.coord(i)*dx;
			for(int i=0; i<8; i++) localCube[i] = 0.0;

			if (phi != NULL)
			{
				localCubePhi[0] = (*phi)(xField);
				localCubePhi[1] = (*phi)(xField+2);
				localCubePhi[2] = (*phi)(xField+1);
				localCubePhi[3] = (*phi)(xField+1+2);
				localCubePhi[4] = (*phi)(xField+0);
				localCubePhi[5] = (*phi)(xField+0+2);
				localCubePhi[6] = (*phi)(xField+0+1);
				localCubePhi[7] = (*phi)(xField+0+1+2);
			}

			for (it=(pcls->field())(xPart).parts.begin(); it != (pcls->field())(xPart).parts.end(); ++it)
			{
				for (int i=0; i<3; i++)
				{
					weightScalarGridUp[i] = ((*it).pos[i] - referPos[i]) / dx;
					weightScalarGridDown[i] = 1.0l - weightScalarGridUp[i];
				}

				if (phi != NULL)
				{
					q = (Real*)((char*)&(*it)+offset_q);

					f = q[0] * q[0] + q[1] * q[1] + q[2] * q[2];
					e = sqrt(f + a * a);
					f = 3. * e + f / e;
				}

				//000
				localCube[0] += weightScalarGridDown[0]*weightScalarGridDown[1]*weightScalarGridDown[2]*(e+f*localCubePhi[0]);
				//001
				localCube[1] += weightScalarGridDown[0]*weightScalarGridDown[1]*weightScalarGridUp[2]*(e+f*localCubePhi[1]);
				//010
				localCube[2] += weightScalarGridDown[0]*weightScalarGridUp[1]*weightScalarGridDown[2]*(e+f*localCubePhi[2]);
				//011
				localCube[3] += weightScalarGridDown[0]*weightScalarGridUp[1]*weightScalarGridUp[2]*(e+f*localCubePhi[3]);
				//100
				localCube[4] += weightScalarGridUp[0]*weightScalarGridDown[1]*weightScalarGridDown[2]*(e+f*localCubePhi[4]);
				//101
				localCube[5] += weightScalarGridUp[0]*weightScalarGridDown[1]*weightScalarGridUp[2]*(e+f*localCubePhi[5]);
				//110
				localCube[6] += weightScalarGridUp[0]*weightScalarGridUp[1]*weightScalarGridDown[2]*(e+f*localCubePhi[6]);
				//111
				localCube[7] += weightScalarGridUp[0]*weightScalarGridUp[1]*weightScalarGridUp[2]*(e+f*localCubePhi[7]);
			}

			(*T00)(xField)	   += localCube[0] * mass;
			(*T00)(xField+2)	 += localCube[1] * mass;
			(*T00)(xField+1)	 += localCube[2] * mass;
			(*T00)(xField+1+2)   += localCube[3] * mass;
			(*T00)(xField+0)	 += localCube[4] * mass;
			(*T00)(xField+0+2)   += localCube[5] * mass;
			(*T00)(xField+0+1)   += localCube[6] * mass;
			(*T00)(xField+0+1+2) += localCube[7] * mass;
		}
	}
}

#define projection_T00_comm scalarProjectionCIC_comm


//////////////////////////
// projection_T0i_project
//////////////////////////
// Description:
//   Particle-mesh projection for T0i, including geometric corrections
//
// Arguments:
//   pcls       pointer to particle handler
//   T0i        pointer to target field
//   phi        pointer to Bardeen potential which characterizes the
//              geometric corrections (volume distortion); can be set to
//              NULL which will result in no corrections applied
//   coeff      coefficient applied to the projection operation (default 1)
//
// Returns:
//
//////////////////////////

template<typename part, typename part_info, typename part_dataType>
void projection_T0i_project(Particles<part,part_info,part_dataType> * pcls, Field<Real> * T0i, Field<Real> * phi = NULL, double coeff = 1.)
{
	if (T0i->lattice().halo() == 0)
	{
		cout<< "projection_T0i_project: target field needs halo > 0" << endl;
		exit(-1);
	}

	Site xPart(pcls->lattice());
	Site xT0i(T0i->lattice());

	typename std::list<part>::iterator it;

	Real referPos[3];
	Real weightScalarGridDown[3];
	Real weightScalarGridUp[3];
	Real dx = pcls->res();

	double mass = coeff / (dx*dx*dx);
	mass *= *(double*)((char*)pcls->parts_info() + pcls->mass_offset());

    Real w;
	Real * q;
	size_t offset_q = offsetof(part,vel);

	Real  qi[12];
	Real  localCubePhi[8];

	for (int i=0; i<8; i++) localCubePhi[i] = 0;

	for(xPart.first(),xT0i.first();xPart.test();xPart.next(),xT0i.next())
	{
		if(pcls->field()(xPart).size!=0)
        {
        	for(int i=0; i<3; i++)
        		referPos[i] = xPart.coord(i)*dx;

            for(int i=0; i<12; i++) qi[i]=0.0;

			if (phi != NULL)
			{
				localCubePhi[0] = (*phi)(xT0i);
				localCubePhi[1] = (*phi)(xT0i+2);
				localCubePhi[2] = (*phi)(xT0i+1);
				localCubePhi[3] = (*phi)(xT0i+1+2);
				localCubePhi[4] = (*phi)(xT0i+0);
				localCubePhi[5] = (*phi)(xT0i+0+2);
				localCubePhi[6] = (*phi)(xT0i+0+1);
				localCubePhi[7] = (*phi)(xT0i+0+1+2);
			}

			for (it=(pcls->field())(xPart).parts.begin(); it != (pcls->field())(xPart).parts.end(); ++it)
			{
				for (int i =0; i<3; i++)
				{
					weightScalarGridUp[i] = ((*it).pos[i] - referPos[i]) / dx;
					weightScalarGridDown[i] = 1.0l - weightScalarGridUp[i];
				}

				q = (Real*)((char*)&(*it)+offset_q);

				w = mass * q[0];

				qi[0] +=  w * weightScalarGridDown[1] * weightScalarGridDown[2];
				qi[1] +=  w * weightScalarGridUp[1]   * weightScalarGridDown[2];
				qi[2] +=  w * weightScalarGridDown[1] * weightScalarGridUp[2];
				qi[3] +=  w * weightScalarGridUp[1]   * weightScalarGridUp[2];

				w = mass * q[1];

				qi[4] +=  w * weightScalarGridDown[0] * weightScalarGridDown[2];
				qi[5] +=  w * weightScalarGridUp[0]   * weightScalarGridDown[2];
				qi[6] +=  w * weightScalarGridDown[0] * weightScalarGridUp[2];
				qi[7] +=  w * weightScalarGridUp[0]   * weightScalarGridUp[2];

                w = mass * q[2];

				qi[8] +=  w * weightScalarGridDown[0] * weightScalarGridDown[1];
				qi[9] +=  w * weightScalarGridUp[0]   * weightScalarGridDown[1];
				qi[10]+=  w * weightScalarGridDown[0] * weightScalarGridUp[1];
				qi[11]+=  w * weightScalarGridUp[0]   * weightScalarGridUp[1];
			}

			(*T0i)(xT0i,0) += qi[0] * (1. + localCubePhi[0] + localCubePhi[4]);
			(*T0i)(xT0i,1) += qi[4] * (1. + localCubePhi[0] + localCubePhi[2]);
			(*T0i)(xT0i,2) += qi[8] * (1. + localCubePhi[0] + localCubePhi[1]);

            (*T0i)(xT0i+0,1) += qi[5] * (1. + localCubePhi[4] + localCubePhi[6]);
            (*T0i)(xT0i+0,2) += qi[9] * (1. + localCubePhi[4] + localCubePhi[5]);

            (*T0i)(xT0i+1,0) += qi[1] * (1. + localCubePhi[2] + localCubePhi[6]);
            (*T0i)(xT0i+1,2) += qi[10] * (1. + localCubePhi[2] + localCubePhi[3]);

            (*T0i)(xT0i+2,0) += qi[2] * (1. + localCubePhi[1] + localCubePhi[5]);
            (*T0i)(xT0i+2,1) += qi[6] * (1. + localCubePhi[1] + localCubePhi[3]);

            (*T0i)(xT0i+1+2,0) += qi[3] * (1. + localCubePhi[3] + localCubePhi[7]);
            (*T0i)(xT0i+0+2,1) += qi[7] * (1. + localCubePhi[5] + localCubePhi[7]);
            (*T0i)(xT0i+0+1,2) += qi[11] * (1. + localCubePhi[6] + localCubePhi[7]);
		}
	}
}

#define projection_T0i_comm vectorProjectionCICNGP_comm


//////////////////////////
// projection_Tij_project
//////////////////////////
// Description:
//   Particle-mesh projection for Tij, including geometric corrections
//
// Arguments:
//   pcls       pointer to particle handler
//   Tij        pointer to target field
//   a          scale factor at projection (needed in order to convert
//              canonical momenta to energies)
//   phi        pointer to Bardeen potential which characterizes the
//              geometric corrections (volume distortion); can be set to
//              NULL which will result in no corrections applied
//   coeff      coefficient applied to the projection operation (default 1)
//
// Returns:
//
//////////////////////////

template<typename part, typename part_info, typename part_dataType>
void projection_Tij_project(Particles<part, part_info, part_dataType> * pcls, Field<Real> * Tij, double a = 1., Field<Real> * phi = NULL, double coeff = 1.)
{
	if (Tij->lattice().halo() == 0)
	{
		cout<< "projection_Tij_project: target field needs halo > 0" << endl;
		exit(-1);
	}

	Site xPart(pcls->lattice());
	Site xTij(Tij->lattice());

	typename std::list<part>::iterator it;

	Real referPos[3];
	Real weightScalarGridDown[3];
	Real weightScalarGridUp[3];
	Real dx = pcls->res();

	double mass = coeff / (dx*dx*dx);
	mass *= *(double*)((char*)pcls->parts_info() + pcls->mass_offset());
	mass /= a;

	Real e, f, w;
	Real * q;
	size_t offset_q = offsetof(part,vel);

	Real  tij[6];           // local cube
	Real  tii[24];          // local cube
	Real  localCubePhi[8];

	for (int i=0; i<8; i++) localCubePhi[i] = 0;

	for (xPart.first(),xTij.first(); xPart.test(); xPart.next(),xTij.next())
	{
		if (pcls->field()(xPart).size != 0)
		{
			for (int i=0;i<3;i++)
				referPos[i] = (double)xPart.coord(i)*dx;

			for (int i=0; i<6; i++)  tij[i]=0.0;
			for (int i=0; i<24; i++) tii[i]=0.0;

			if (phi != NULL)
			{
				localCubePhi[0] = (*phi)(xTij);
				localCubePhi[1] = (*phi)(xTij+2);
				localCubePhi[2] = (*phi)(xTij+1);
				localCubePhi[3] = (*phi)(xTij+1+2);
				localCubePhi[4] = (*phi)(xTij+0);
				localCubePhi[5] = (*phi)(xTij+0+2);
				localCubePhi[6] = (*phi)(xTij+0+1);
				localCubePhi[7] = (*phi)(xTij+0+1+2);
			}

			for (it=(pcls->field())(xPart).parts.begin(); it != (pcls->field())(xPart).parts.end(); ++it)
			{
				for (int i =0; i<3; i++)
				{
					weightScalarGridUp[i] = ((*it).pos[i] - referPos[i]) / dx;
					weightScalarGridDown[i] = 1.0l - weightScalarGridUp[i];
				}

				q = (Real*)((char*)&(*it)+offset_q);
				f = q[0] * q[0] + q[1] * q[1] + q[2] * q[2];
				e = sqrt(f + a * a);
				f = 4. + a * a / (f + a * a);

				// diagonal components
				for (int i=0; i<3; i++)
				{
					w = mass * q[i] * q[i] / e;
					//000
					tii[0+i*8] += w * weightScalarGridDown[0] * weightScalarGridDown[1] * weightScalarGridDown[2] * (1. + f * localCubePhi[0]);
					//001
					tii[1+i*8] += w * weightScalarGridDown[0] * weightScalarGridDown[1] * weightScalarGridUp[2]   * (1. + f * localCubePhi[1]);
					//010
					tii[2+i*8] += w * weightScalarGridDown[0] * weightScalarGridUp[1]   * weightScalarGridDown[2] * (1. + f * localCubePhi[2]);
					//011
					tii[3+i*8] += w * weightScalarGridDown[0] * weightScalarGridUp[1]   * weightScalarGridUp[2]   * (1. + f * localCubePhi[3]);
					//100
					tii[4+i*8] += w * weightScalarGridUp[0]   * weightScalarGridDown[1] * weightScalarGridDown[2] * (1. + f * localCubePhi[4]);
					//101
					tii[5+i*8] += w * weightScalarGridUp[0]   * weightScalarGridDown[1] * weightScalarGridUp[2]   * (1. + f * localCubePhi[5]);
					//110
					tii[6+i*8] += w * weightScalarGridUp[0]   * weightScalarGridUp[1]   * weightScalarGridDown[2] * (1. + f * localCubePhi[6]);
					//111
					tii[7+i*8] += w * weightScalarGridUp[0]   * weightScalarGridUp[1]   * weightScalarGridUp[2]   * (1. + f * localCubePhi[7]);
				}

				w = mass * q[0] * q[1] / e;
				tij[0] +=  w * weightScalarGridDown[2] * (1. + f * 0.25 * (localCubePhi[0] + localCubePhi[2] + localCubePhi[4] + localCubePhi[6]));
				tij[1] +=  w * weightScalarGridUp[2] * (1. + f * 0.25 * (localCubePhi[1] + localCubePhi[3] + localCubePhi[5] + localCubePhi[7]));

				w = mass * q[0] * q[2] / e;
				tij[2] +=  w * weightScalarGridDown[1] * (1. + f * 0.25 * (localCubePhi[0] + localCubePhi[1] + localCubePhi[4] + localCubePhi[5]));
				tij[3] +=  w * weightScalarGridUp[1] * (1. + f * 0.25 * (localCubePhi[2] + localCubePhi[3] + localCubePhi[6] + localCubePhi[7]));

				w = mass * q[1] * q[2] / e;
				tij[4] +=  w * weightScalarGridDown[0] * (1. + f * 0.25 * (localCubePhi[0] + localCubePhi[1] + localCubePhi[2] + localCubePhi[3]));
				tij[5] +=  w * weightScalarGridUp[0] * (1. + f * 0.25 * (localCubePhi[4] + localCubePhi[5] + localCubePhi[6] + localCubePhi[7]));

			}


			for (int i=0; i<3; i++) (*Tij)(xTij,i,i) += tii[8*i];
			(*Tij)(xTij,0,1) += tij[0];
			(*Tij)(xTij,0,2) += tij[2];
			(*Tij)(xTij,1,2) += tij[4];

			for (int i=0; i<3; i++) (*Tij)(xTij+0,i,i) += tii[4+8*i];
			(*Tij)(xTij+0,1,2) += tij[5];

			for (int i=0; i<3; i++) (*Tij)(xTij+1,i,i) += tii[2+8*i];
			(*Tij)(xTij+1,0,2) += tij[3];

			for (int i=0; i<3; i++) (*Tij)(xTij+2,i,i) += tii[1+8*i];
			(*Tij)(xTij+2,0,1) += tij[1];

			for (int i=0; i<3; i++) (*Tij)(xTij+0+1,i,i) += tii[6+8*i];
			for (int i=0; i<3; i++) (*Tij)(xTij+0+2,i,i) += tii[5+8*i];
			for (int i=0; i<3; i++) (*Tij)(xTij+1+2,i,i) += tii[3+8*i];
			for (int i=0; i<3; i++) (*Tij)(xTij+0+1+2,i,i) += tii[7+8*i];
		}
	}
}

#ifndef projection_Tij_comm
#define projection_Tij_comm symtensorProjectionCICNGP_comm
#endif


//////////////////////////
// projection_Ti0_project
//////////////////////////
// Description:
//   Particle-mesh projection for Ti0, including geometric corrections
//
// Arguments:
//   pcls       pointer to particle handler
//   Ti0        pointer to target field
//   phi        pointer to Bardeen potential which characterizes the
//              geometric corrections (volume distortion); can be set to
//              NULL which will result in no corrections applied
//   chi        pointer to difference between the Bardeen potentials which
//              characterizes additional corrections; can be set to
//              NULL which will result in no corrections applied
//   coeff      coefficient applied to the projection operation (default 1)
//
// Returns:
//
//////////////////////////

template<typename part, typename part_info, typename part_dataType>
void projection_Ti0_project(Particles<part, part_info, part_dataType> * pcls, Field<Real> * Ti0, Field<Real> * phi = NULL, Field<Real> * chi = NULL, double coeff = 1.)
{
	if (Ti0->lattice().halo() == 0)
	{
		cout<< "projection_Ti0_project: target field needs halo > 0" << endl;
		exit(-1);
	}

	Site xPart(pcls->lattice());
	Site xField(Ti0->lattice());

	typename std::list<part>::iterator it;

	Real referPos[3];
	Real weightScalarGridUp[3];
	Real weightScalarGridDown[3];
	Real dx = pcls->res();

	Real * q;
	size_t offset_q = offsetof(part,vel);


	double mass = coeff / (dx*dx*dx);
	mass *= *(double*)((char*)pcls->parts_info() + pcls->mass_offset());


	Real localCube[24]; // XYZ = 000 | 001 | 010 | 011 | 100 | 101 | 110 | 111
	Real localCubePhi[8];
	Real localCubeChi[8];

	for (int i = 0; i < 8; i++) localCubePhi[i] = 0.0;
	for (int i = 0; i < 8; i++) localCubeChi[i] = 0.0;

	for (xPart.first(), xField.first(); xPart.test(); xPart.next(), xField.next())
	{
		if (pcls->field()(xPart).size != 0)
		{
			for(int i = 0; i < 3; i++) referPos[i] = xPart.coord(i)*dx;
			for(int i = 0; i < 24; i++) localCube[i] = 0.0;

			if (phi != NULL)
			{
				localCubePhi[0] = (*phi)(xField);
				localCubePhi[1] = (*phi)(xField+2);
				localCubePhi[2] = (*phi)(xField+1);
				localCubePhi[3] = (*phi)(xField+1+2);
				localCubePhi[4] = (*phi)(xField+0);
				localCubePhi[5] = (*phi)(xField+0+2);
				localCubePhi[6] = (*phi)(xField+0+1);
				localCubePhi[7] = (*phi)(xField+0+1+2);
			}
			if (chi != NULL)
			{
				localCubeChi[0] = (*chi)(xField);
				localCubeChi[1] = (*chi)(xField+2);
				localCubeChi[2] = (*chi)(xField+1);
				localCubeChi[3] = (*chi)(xField+1+2);
				localCubeChi[4] = (*chi)(xField+0);
				localCubeChi[5] = (*chi)(xField+0+2);
				localCubeChi[6] = (*chi)(xField+0+1);
				localCubeChi[7] = (*chi)(xField+0+1+2);
			}

			for (it = (pcls->field())(xPart).parts.begin(); it != (pcls->field())(xPart).parts.end(); ++it)
			{
				for (int i = 0; i < 3; i++)
				{
					weightScalarGridUp[i] = ((*it).pos[i] - referPos[i]) / dx;
					weightScalarGridDown[i] = 1.0l - weightScalarGridUp[i];
				}

				q = (Real*)((char*)&(*it)+offset_q);

                for (int i = 0; i < 3; i++){
                    //000
                    localCube[8*i] += weightScalarGridDown[0]*weightScalarGridDown[1]*weightScalarGridDown[2]*q[i]*(1.+6*localCubePhi[0]-localCubeChi[0]);
                    //001
                    localCube[8*i+1] += weightScalarGridDown[0]*weightScalarGridDown[1]*weightScalarGridUp[2]*q[i]*(1.+6*localCubePhi[1]-localCubeChi[1]);
                    //010
                    localCube[8*i+2] += weightScalarGridDown[0]*weightScalarGridUp[1]*weightScalarGridDown[2]*q[i]*(1.+6*localCubePhi[2]-localCubeChi[2]);
                    //011
                    localCube[8*i+3] += weightScalarGridDown[0]*weightScalarGridUp[1]*weightScalarGridUp[2]*q[i]*(1.+6*localCubePhi[3]-localCubeChi[3]);
                    //100
                    localCube[8*i+4] += weightScalarGridUp[0]*weightScalarGridDown[1]*weightScalarGridDown[2]*q[i]*(1.+6*localCubePhi[4]-localCubeChi[4]);
                    //101
                    localCube[8*i+5] += weightScalarGridUp[0]*weightScalarGridDown[1]*weightScalarGridUp[2]*q[i]*(1.+6*localCubePhi[5]-localCubeChi[5]);
                    //110
                    localCube[8*i+6] += weightScalarGridUp[0]*weightScalarGridUp[1]*weightScalarGridDown[2]*q[i]*(1.+6*localCubePhi[6]-localCubeChi[6]);
                    //111
                    localCube[8*i+7] += weightScalarGridUp[0]*weightScalarGridUp[1]*weightScalarGridUp[2]*q[i]*(1.+6*localCubePhi[7]-localCubeChi[7]);
                }
			}
			for (int i = 0; i < 3; i++)
            {
                (*Ti0)(xField,i)       += localCube[8*i] * mass;
                (*Ti0)(xField+2,i)     += localCube[8*i+1] * mass;
                (*Ti0)(xField+1,i)     += localCube[8*i+2] * mass;
                (*Ti0)(xField+1+2,i)   += localCube[8*i+3] * mass;
                (*Ti0)(xField+0,i)     += localCube[8*i+4] * mass;
                (*Ti0)(xField+0+2,i)   += localCube[8*i+5] * mass;
                (*Ti0)(xField+0+1,i)   += localCube[8*i+6] * mass;
                (*Ti0)(xField+0+1+2,i) += localCube[8*i+7] * mass;
            }
		}
	}
}


//////////////////////////
// projectFTtheta
//////////////////////////
// Description:
//   Compute the diverge of the velocity in Fourier space
//
// Arguments:
//   thFT       reference to the Fourier image of the divergence of the velocity field
//   viFT       reference to the Fourier image of the velocity field
//
// Returns:
//
//////////////////////////

void projectFTtheta(Field<Cplx> & thFT, Field<Cplx> & viFT)
{
	const int linesize = thFT.lattice().size(1);
	int i;
	Real * gridk;
	rKSite k(thFT.lattice());
	Cplx tmp(0., 0.);

	gridk = (Real *) malloc(linesize * sizeof(Real));

	for (i = 0; i < linesize; i++)
		gridk[i] = (Real) linesize * sin(M_PI * 2.0 * (Real) i / (Real) linesize);

	for (k.first(); k.test(); k.next())
		thFT(k) = Cplx(0.,1.)*(gridk[k.coord(0)] * viFT(k,0) + gridk[k.coord(1)] * viFT(k,1) + gridk[k.coord(2)] * viFT(k,2));

	free(gridk);
}


//////////////////////////
// projectFTomega
//////////////////////////
// Description:
//   Compute the curl part of the velocity field in Fourier space
//
// Arguments:
//   viFT      reference to the input Fourier image of the velocity field
//             the divergence part will be projected out
//
// Returns:
//
//////////////////////////

void projectFTomega(Field<Cplx> & viFT)
{
	const int linesize = viFT.lattice().size(1);
	int i;
	Real * gridk2;
	Cplx * kshift;
	Real * gridk;
	rKSite k(viFT.lattice());
	Cplx tmp(0., 0.);
	Cplx vr[3];

	gridk2 = (Real *) malloc(linesize * sizeof(Real));
	gridk = (Real *) malloc(linesize * sizeof(Real));

	for (i = 0; i < linesize; i++)
	{
		gridk[i] = (Real) linesize * sin(M_PI * 2.0 * (Real) i / (Real) linesize);
		gridk2[i] = gridk[i]*gridk[i];
    }

	k.first();
	if (k.coord(0) == 0 && k.coord(1) == 0 && k.coord(2) == 0)
	{
		viFT(k,0) = Cplx(0.,0.);
		viFT(k,1) = Cplx(0.,0.);
		viFT(k,2) = Cplx(0.,0.);
		k.next();
	}

	for (; k.test(); k.next())
	{
		if ((k.coord(0) == 0 || k.coord(0) == linesize/2) && (k.coord(1) == 0 || k.coord(1) == linesize/2) && (k.coord(2) == 0 || k.coord(2) == linesize/2))
		{
			viFT(k, 0) = Cplx(0.,0.);
			viFT(k, 1) = Cplx(0.,0.);
			viFT(k, 2) = Cplx(0.,0.);
		}
		else
		{
			tmp = (gridk[k.coord(0)] * viFT(k,0) + gridk[k.coord(1)] * viFT(k,1) + gridk[k.coord(2)] * viFT(k,2)) / (gridk2[k.coord(0)] + gridk2[k.coord(1)] + gridk2[k.coord(2)]);

			vr[0] = (viFT(k,0) - gridk[k.coord(0)] * tmp);
			vr[1] = (viFT(k,1) - gridk[k.coord(1)] * tmp);
			vr[2] = (viFT(k,2) - gridk[k.coord(2)] * tmp);

			viFT(k,0) = Cplx(0.,1.)*(gridk[k.coord(1)]*vr[2] - gridk[k.coord(2)]*vr[1]);
			viFT(k,1) = Cplx(0.,1.)*(gridk[k.coord(2)]*vr[0] - gridk[k.coord(0)]*vr[2]);
			viFT(k,2) = Cplx(0.,1.)*(gridk[k.coord(0)]*vr[1] - gridk[k.coord(1)]*vr[0]);
		}
	}

	free(gridk2);
}

#endif
