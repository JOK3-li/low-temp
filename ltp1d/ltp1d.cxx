/*******************************************************************
 * Low temperature plasma model
 *
 * J. Leddy, March 2017
 *******************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>
#include <invert_laplace.hxx>
#include <bout/invert/laplacexz.hxx>
#include <derivs.hxx>
#include <initialprofiles.hxx>

// Evolving variables
Field3D Ne, NeE;		 // electron density and energy and particle flux
Field3D Ng;				// Gas density
Field3D Ni; 			// ion density
Field3D Te, E; 			// Temeprature
Field3D phi, phiext, total_phi;		// potential
Field3D Sn, Se; 		// density and temperature sources due to interactions
Field3D linear;			// linear function in X

// Laplacian inversion
Laplacian *phiSolver; // Old Laplacian in X-Z

// parameters
BoutReal mu_e, mu_i, Dn, De, Ti, nu;
BoutReal w0, L0, n0, v0, T0;	// normalisations
BoutReal phimag, Efreq;  // imposed field
const BoutReal PI = acos(-1.);  // Pi
const BoutReal ee = 1.60217662e-19; // C
const BoutReal eps0 = 8.854187817e-12; // vacuum permitivity
const BoutReal Mn = 6.6335209e-26; // kg
const BoutReal Me = 9.10938356e-31; // kg
const BoutReal k_B = 1.0;//1.38064852e-23; // boltzmann constant
bool evolve_Ne, evolve_Ni, evolve_NeE, evolve_Ng;

Options *opt = Options::getRoot();

const Field3D Div_Perp_Lap_FV(const Field3D &a, const Field3D &f, const Field3D &bndry_flux);

int physics_init(bool restarting)
{
	// read options
	Coordinates *coords = mesh->coordinates();
	Options *options = Options::getRoot();
	options = options->getSection("ltp");
	OPTION(options, L0,  1.0);
	OPTION(options, w0,  1e4);
	OPTION(options, n0,  1e15);
	OPTION(options, T0,  1.0);
	OPTION(options, phimag,  0.0);
	OPTION(options, Efreq,  0.0);
	OPTION(options, evolve_Ne, true);
	OPTION(options, evolve_Ni, true);
	OPTION(options, evolve_NeE, true);
	OPTION(options, evolve_Ng, true);
	OPTION(options, mu_i, 0.5);// m^2 / Vs
	OPTION(options, mu_e, 156.0);// m^2 / Vs

	// constants
	Dn = 0.1;  // m2/s
	De = 0.3;  // m2/s
	Ti = 100;  // eV
	nu = 1.0; // charge exchange frequency

	// Create laplacian solver
	phiSolver = Laplacian::create(opt->getSection("phiSolver"));

	// Set evolving variables
	initial_profile("linear", linear);
	if(evolve_Ng) {
		bout_solve(Ng, "Ng");
	} else {
		initial_profile("Ng", Ng);
		SAVE_ONCE(Ng);
	}
	if(evolve_NeE) {
		bout_solve(NeE, "NeE");
	}
	if(evolve_Ne) {
		bout_solve(Ne, "Ne");
	}
	if(evolve_Ni) {
		bout_solve(Ni, "Ni");
	}
	SAVE_REPEAT4(phiext,total_phi,phi,E);
	SAVE_ONCE2(phimag,Efreq);
	SAVE_ONCE4(w0,T0,n0,L0);

	// Normalisations
	v0 = ee / eps0; // voltage normalisation
	Ti /= T0;
	mu_i *= v0 / (w0 * L0 * L0);
	mu_e *= v0 / (w0 * L0 * L0);
	Dn /= w0*L0*L0;
	De /= w0*L0*L0;
	coords->dx /= L0;
	coords->dz /= L0;
	coords->J = 1.0;
	phimag /= v0;
	Efreq *= 2. * PI / w0;
	Ne /= n0;
	Ni /= n0;
	Ng /= n0;
	NeE /= n0;

	// initial potential
	phi = 0.0;
	phiext = 0.0;
	total_phi = 0.0;
	E = NeE / Ne;

	return 0;
}

int physics_run(BoutReal t)
{
	mesh->communicate(Ne,Ni,NeE,Ng);//,Fex,Fix,Fez,Fiz);

	Ne = floor(Ne,1e-5);
	Ni = floor(Ni,1e-5);
	NeE = floor(Ne,1e-5);
	Ng = floor(Ng,1e-5);

	// invert Poisson equation
	phi = phiSolver->solve(Ni-Ne);
	mesh->communicate(phi);
	phiext = phimag * linear * sin(Efreq * t);
	total_phi = phi + phiext;
	mesh->communicate(total_phi);

	// calculate energy and temperature
	E = NeE / Ne;
	Te = 2.*E/3.;
	Field3D zero = 0.0;

	// Gas continuity equation
	ddt(Ng) = 0.0;
	if(evolve_Ng) {
		// Electron density
		ddt(Ng) += -Div_Perp_Lap_FV(Dn,Ng,zero);
		ddt(Ng) += Ne * Ni * 0.01 * NeE;
		ddt(Ng) += -Ne * Ng * 0.001 * NeE;
	}

	// Electron continuity equation
	ddt(Ne) = 0.0;
	if(evolve_Ne) {
		// Electron density
		ddt(Ne) += -Div_Perp_Lap_FV(-Ne*mu_e,total_phi,0.25*Ne*sqrt(Te));
		ddt(Ne) += -Div_Perp_Lap_FV(Dn,Ne,zero);
		// Density source/sink from interactions
		ddt(Ne) += -Ne * Ni * 0.01 * NeE; //-Ne * Ni * recombination(Ne*n0, Te*T0) * n0 / w0;
		ddt(Ne) += Ne * Ng * 0.001 * NeE; // Ne * Ng * ionisation(Te*Tnorm) * n0 / w0;
	}

	// Ion continuity equation
	ddt(Ni) = 0.0;
	if(evolve_Ni) {
		// Ion density
		ddt(Ni) += -Div_Perp_Lap_FV(Ni*mu_i,total_phi,0.25*Ni*sqrt(Ti));
		ddt(Ni) += -Div_Perp_Lap_FV(Dn,Ni,zero);
		// Density source/sink from interactions
		ddt(Ni) += -Ne * Ni * 0.01 * NeE; //-Ne * Ni * recombination(Ne*n0, Te*T0) * n0 / w0;
		ddt(Ni) += Ne * Ng * 0.001 * NeE; // Ne * Ng * ionisation(Te*Tnorm) * n0 / w0;
	}

	// Electron energy density equation
	ddt(NeE) = 0.0;
	if(evolve_NeE) {
		// Energy density
		ddt(NeE) += Div_Perp_Lap_FV(5./3.*NeE*mu_e,total_phi,5./12.*NeE*sqrt(Ti));
		ddt(NeE) += -Div_Perp_Lap_FV(5./3.*E*Dn,Ne,zero);
		ddt(NeE) += -Div_Perp_Lap_FV(5./3.*Ne*De,E,zero);

		// Energy source from interactions
		Se = 0.0; // Ne * Ng * excitation(Te*Tnorm) * n0 / w0;

		ddt(NeE) += 3. * Me/Mn * k_B * nu * Ne * (Te - Ti);
		ddt(NeE) += Se;
	}

	return 0;
}

// Div ( a * Grad_perp(f) )
const Field3D Div_Perp_Lap_FV(const Field3D &a, const Field3D &f, const Field3D &bndry_flux) {

  Coordinates *coords = mesh->coordinates();
  Field3D result = 0.0;

  //////////////////////////////////////////
  // X-Z diffusion
  //
  //			Z
  //			|
  //
  //	 o --- gU --- o
  //	 |	 nU	 |
  //	 |			|
  //	gL nL	  nR gR	-> X
  //	 |			|
  //	 |	 nD	 |
  //	 o --- gD --- o
  //


  Field3D fs = f;
  Field3D as = a;

  for(int i=mesh->xstart;i<=mesh->xend;i++)
	for(int j=mesh->ystart;j<=mesh->yend;j++)
	  for(int k=0;k<mesh->LocalNz;k++) {

		// Calculate gradients on cell faces
		BoutReal gRp = (coords->g11(i+1,j) + coords->g11(i+2,j)) * (fs(i+2,j,k) - fs(i+1,j,k))/(coords->dx(i+2,j) + coords->dx(i+1,j));
		BoutReal gR = (coords->g11(i,j) + coords->g11(i+1,j)) * (fs(i+1,j,k) - fs(i,j,k))/(coords->dx(i+1,j) + coords->dx(i,j));
		BoutReal gL = (coords->g11(i-1,j) + coords->g11(i,j))*(fs(i,j,k) - fs(i-1,j,k))/(coords->dx(i-1,j) + coords->dx(i,j));
		BoutReal gLm = (coords->g11(i-2,j) + coords->g11(i-1,j))*(fs(i-1,j,k) - fs(i-2,j,k))/(coords->dx(i-2,j) + coords->dx(i-1,j));

		// Flow right
		BoutReal flux = gR * 0.25*(coords->J(i+1,j) + coords->J(i,j)) *(as(i+1,j,k) + as(i,j,k));
		BoutReal fluxp = gRp * 0.25*(coords->J(i+2,j) + coords->J(i+1,j)) *(as(i+2,j,k) + as(i+1,j,k));
		if(i==mesh->xend) {
			if(flux<0.0){
				flux = 0.0;
			}
			if(fluxp<0.0){
				fluxp = 0.0;
			}
			flux += bndry_flux(i,j,k); //+
			fluxp += bndry_flux(i+1,j,k); //+
		}
		result(i,j,k) += (fluxp - flux) / (2.*coords->dx(i,j)*coords->J(i,j));
		// result(i,j,k) += flux / (coords->dx(i,j)*coords->J(i,j));
		//result(i+1,j,k) -= flux / (mesh->dx(i+1,j)*mesh->J(i+1,j));

		// Flow left
		flux = gL * 0.25*(coords->J(i-1,j) + coords->J(i,j)) *(as(i-1,j,k) + as(i,j,k));
		BoutReal fluxm = gLm * 0.25*(coords->J(i-2,j) + coords->J(i-1,j)) *(as(i-2,j,k) + as(i-1,j,k));
		if(i==mesh->xstart) {
			if(flux>0.0) {
				flux = 0.0;
			}
			if(fluxm>0.0) {
				fluxm = 0.0;
			}
			flux += -bndry_flux(i,j,k); //-
			fluxm += -bndry_flux(i-1,j,k); //-
		}
		result(i,j,k) -= (flux - fluxm) / (2.*coords->dx(i,j)*coords->J(i,j));
		// result(i,j,k) -= flux / (coords->dx(i,j)*coords->J(i,j));
		//result(i-1,j,k) += flux / (mesh->dx(i+1,j)*mesh->J(i+1,j));
	  }

  return result;
}
