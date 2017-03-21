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
BoutReal Eext, Efreq;  // imposed field
const BoutReal PI = acos(-1.);  // Pi
const BoutReal ee = 1.60217662e-19; // C
const BoutReal eps0 = 8.854187817e-12; // vacuum permitivity
const BoutReal Mn = 1.6726219e-27; // kg
const BoutReal Me = 9.10938356e-31; // kg
const BoutReal k_B = 1.0;//1.38064852e-23; // boltzmann constant
bool evolve_Ne, evolve_Ni, evolve_NeE, evolve_Ng;

Options *opt = Options::getRoot();

const Field3D Div_Perp_Lap_FV(const Field3D &a, const Field3D &f, bool xflux);

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
	OPTION(options, Eext,  0.0);
	OPTION(options, Efreq,  0.0);
	OPTION(options, evolve_Ne, true);
	OPTION(options, evolve_Ni, true);
	OPTION(options, evolve_NeE, true);
	OPTION(options, evolve_Ng, true);
	OPTION(options, mu_i, 1.0);// m^2 / Vs
	OPTION(options, mu_e, 10.0);// m^2 / Vs

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
	SAVE_ONCE2(Eext,Efreq);
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
	Eext *= L0 / v0;
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
	phiext = Eext * linear * cos(Efreq * t);
	total_phi = phi + phiext;
	mesh->communicate(total_phi);

	// calculate energy and temperature
	E = NeE / Ne;
	Te = 2.*E/3.;

	// Gas continuity equation
	ddt(Ng) = 0.0;
	if(evolve_Ng) {
		// Electron density
		ddt(Ng) += -Div_Perp_Lap_FV(Dn,Ng,false);
		ddt(Ng) += Ne * Ni * 0.01 * NeE;
		ddt(Ng) += -Ni * Ng * 0.01 * NeE;
		ddt(Ng) += -Ne * Ng * 0.001 * NeE;
	}

	// Electron continuity equation
	ddt(Ne) = 0.0;
	if(evolve_Ne) {
		// Electron density
		ddt(Ne) += Div_Perp_Lap_FV(Ne*mu_e,total_phi,false);
		ddt(Ne) += -Div_Perp_Lap_FV(Dn,Ne,false);
		ddt(Ne) += -Ne * Ni * 0.01 * NeE;
		ddt(Ne) += Ni * Ng * 0.01 * NeE;
		ddt(Ne) += Ne * Ng * 0.001 * NeE;
	    // Density source from interactions
		// BoutReal R_rc  = hydrogen.recombination(Ne*n0, Te*T0)*SQ(Ne) * n0 / w0;
		// BoutReal R_iz = Ne*Nn_C*hydrogen.ionisation(Te_C*Tnorm) * Nnorm / Omega_ci;
	}

	// Ion continuity equation
	ddt(Ni) = 0.0;
	if(evolve_Ni) {
		// Ion density
		ddt(Ni) += -Div_Perp_Lap_FV(Ni*mu_i,total_phi,false);
		ddt(Ni) += -Div_Perp_Lap_FV(Dn,Ni,false);
		ddt(Ni) += -Ne * Ni * 0.01 * NeE;
		ddt(Ni) += Ni * Ng * 0.01 * NeE;
		ddt(Ni) += Ne * Ng * 0.001 * NeE;
	    // Density source from interactions
		// BoutReal R_rc  = hydrogen.recombination(Ne*n0, Te*T0)*SQ(Ne) * n0 / w0;
		// BoutReal R_iz = Ne*Nn_C*hydrogen.ionisation(Te_C*Tnorm) * Nnorm / Omega_ci;
	}

	// Electron energy density equation
	ddt(NeE) = 0.0;
	if(evolve_NeE) {
		// Energy density
		ddt(NeE) += Div_Perp_Lap_FV(5./3.*NeE*mu_e,total_phi,false);
		ddt(NeE) += -Div_Perp_Lap_FV(5./3.*E*Dn,Ne,false);
		ddt(NeE) += -Div_Perp_Lap_FV(5./3.*Ne*De,E,false);

		// Energy source from interactions
		Se = 0.0; // nj * ni * ke_ij

		ddt(NeE) += 3. * Me/Mn * k_B * nu * Ne * (Te - Ti);
		ddt(NeE) += Se;
	}

 	///////////////////////////////////////////////
	// Set boundary conditions
	for(RangeIterator r=mesh->iterateBndryLowerY(); !r.isDone(); r++) {
  		for(int jz=0; jz<mesh->LocalNz; jz++) {
			// Zero-gradient density
			BoutReal nesheath = 0.5*( 3.*Ne(r.ind, mesh->ystart, jz) - Ne(r.ind, mesh->ystart+1, jz) );
			if(nesheath < 0.0)
				nesheath = 0.0;
			BoutReal nisheath = 0.5*( 3.*Ni(r.ind, mesh->ystart, jz) - Ne(r.ind, mesh->ystart+1, jz) );
			if(nisheath < 0.0)
				nisheath = 0.0;

			// Temperature at the sheath entrance
			BoutReal tesheath = Te(r.ind, mesh->ystart, jz);
			BoutReal tisheath = Ti;//(r.ind, mesh->ystart, jz);

			// Zero-gradient potential
			BoutReal phisheath = phi(r.ind, mesh->ystart, jz);

			// // Ion velocity goes to the sound speed
			// BoutReal visheath = -sqrt(tesheath + tisheath); // Sound speed outwards
			//
			// if(Vi(r.ind, mesh->ystart, jz) < visheath) {
			//   // If plasma is faster, go to plasma velocity
			// 	visheath = Vi(r.ind, mesh->ystart, jz);
			// }

			// Sheath current
			BoutReal phi_te = phisheath / Te(r.ind, mesh->ystart, jz);
			// BoutReal vesheath = -sqrt(tesheath) * (sqrt(mi_me)/(2.*sqrt(PI))) * exp(-phi_te);
			// J = n*(Vi - Ve)
			// BoutReal jsheath = nesheath * (visheath - vesheath);
			// if(nesheath < 1e-10) {
			// 	vesheath = visheath;
			// 	jsheath = 0.0;
			// }

			// Apply boundary condition half-way between cells
			for(int jy = mesh->ystart-1;jy >= 0; jy--) {
				// Neumann conditions
				phi(r.ind, jy, jz) = phisheath;
				// Vort(r.ind, jy, jz) = Vort(r.ind, mesh->ystart, jz);

				// Here zero-gradient Te, heat flux applied later
				Te(r.ind, jy, jz) = Te(r.ind, mesh->ystart, jz);

				// Dirichlet conditions
				Ne(r.ind, jy, jz) = 2.*nesheath - Ne(r.ind, mesh->ystart, jz);
				Ni(r.ind, jy, jz) = 2.*nisheath - Ni(r.ind, mesh->ystart, jz);
				//E(r.ind, jy, jz) = 2.*nesheath*tesheath - E(r.ind, mesh->ystart, jz);

				// Vi(r.ind, jy, jz)  = 2.*visheath - Vi(r.ind, mesh->ystart, jz);
				// Ve(r.ind, jy, jz)  = 2.*vesheath - Ve(r.ind, mesh->ystart, jz);
				// Jpar(r.ind, jy, jz) = 2.*jsheath - Jpar(r.ind, mesh->ystart, jz);
				// NVi(r.ind, jy, jz) = 2.*nesheath*visheath - NVi(r.ind, mesh->ystart, jz);
			}
		}
	}

	return 0;
}

// Div ( a * Grad_perp(f) )
const Field3D Div_Perp_Lap_FV(const Field3D &a, const Field3D &f, bool xflux) {

  Coordinates *coords = mesh->coordinates();
  Field3D result = 0.0;

  //////////////////////////////////////////
  // X-Z diffusion
  //
  //            Z
  //            |
  //
  //     o --- gU --- o
  //     |     nU     |
  //     |            |
  //    gL nL      nR gR    -> X
  //     |            |
  //     |     nD     |
  //     o --- gD --- o
  //


  Field3D fs = f;
  Field3D as = a;

  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++)
      for(int k=0;k<mesh->LocalNz;k++) {
        int kp = (k+1) % (mesh->LocalNz);
		int km = (k-1+mesh->LocalNz) % (mesh->LocalNz);
		//output<<kp<<" "<<km<<"\n";

        // Calculate gradients on cell faces
        BoutReal gR = (coords->g11(i,j) + coords->g11(i+1,j)) * (fs(i+1,j,k) - fs(i,j,k))/(coords->dx(i+1,j) + coords->dx(i,j));
        BoutReal gL = (coords->g11(i-1,j) + coords->g11(i,j))*(fs(i,j,k) - fs(i-1,j,k))/(coords->dx(i-1,j) + coords->dx(i,j));
        BoutReal gD = coords->g33(i,j)*(fs(i,j,k) - fs(i,j,km))/coords->dz;
        BoutReal gU = coords->g33(i,j)*(fs(i,j,kp) - fs(i,j,k))/coords->dz;

        // Flow right
        BoutReal flux = gR * 0.25*(coords->J(i+1,j) + coords->J(i,j)) *(as(i+1,j,k) + as(i,j,k));
        result(i,j,k)   += flux / (coords->dx(i,j)*coords->J(i,j));
        //result(i+1,j,k) -= flux / (mesh->dx(i+1,j)*mesh->J(i+1,j));

        // Flow left
        flux = gL * 0.25*(coords->J(i-1,j) + coords->J(i,j)) *(as(i-1,j,k) + as(i,j,k));
        result(i,j,k)   -= flux / (coords->dx(i,j)*coords->J(i,j));
        //result(i-1,j,k) += flux / (mesh->dx(i+1,j)*mesh->J(i+1,j));

        // Flow up
        flux = gU * 0.5*(as(i,j,k) + as(i,j,kp)) / coords->dz;
        result(i,j,k) += flux;
        //result(i,j,kp) -= flux;

		// Flow down
        flux = gD * 0.5*(as(i,j,k) + as(i,j,km)) / coords->dz;
        result(i,j,k) -= flux;
        //result(i,j,km) += flux;
      }

  return result;
}
