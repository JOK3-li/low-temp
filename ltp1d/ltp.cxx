/*******************************************************************
 * Low temperature plasma model
 *
 * J. Leddy, April 2017
 *******************************************************************/

#include <bout.hxx>
#include <bout/physicsmodel.hxx>
#include <bout/constants.hxx>
#include <boutmain.hxx>
#include <invert_laplace.hxx>
#include <bout/invert/laplacexz.hxx>
#include <derivs.hxx>
#include <initialprofiles.hxx>
#include "radiation.hxx"

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
bool evolve_Ne, evolve_Ni, evolve_NeE, evolve_Ng;

// neutral interaction
UpdatedRadiatedPower hydrogen; // Atomic rates (H.Willett)
Field3D S,F,Q,R;
Field3D Riz,Rrc,Rcx;
BoutReal Eionize;

const Field3D Div_Perp_Lap_FV(const Field3D &a, const Field3D &f, const Field3D &bndry_flux);
void neutral_rates(const Field3D &Ne, const Field3D &Te, const Field3D &Vi,const Field3D &Nn, const Field3D &Tn, const Field3D &Vnpar,
                      Field3D &S, Field3D &F, Field3D &Q, Field3D &R,Field3D &Riz, Field3D &Rrc, Field3D &Rcx);

const BoutReal ee = 1.60217662e-19; // C
const BoutReal eps0 = 8.854187817e-12; // vacuum permitivity
const BoutReal Mn = 6.6335209e-26; // kg
const BoutReal Me = 9.10938356e-31; // kg
const BoutReal k_B = 1.38064852e-23; // boltzmann constant

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
	phiSolver = Laplacian::create(options->getSection("phiSolver"));

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
		BoutReal gR = (coords->g11(i,j) + coords->g11(i+1,j)) * (fs(i+1,j,k) - fs(i,j,k))/(coords->dx(i+1,j) + coords->dx(i,j));
		BoutReal gL = (coords->g11(i-1,j) + coords->g11(i,j))*(fs(i,j,k) - fs(i-1,j,k))/(coords->dx(i-1,j) + coords->dx(i,j));

		// Flow right
		BoutReal flux = gR * 0.25*(coords->J(i+1,j) + coords->J(i,j)) *(as(i+1,j,k) + as(i,j,k));
		if(i==mesh->xend) {
			if(flux<0.0){
				flux = 0.0;
			}
			flux += bndry_flux(i,j,k); //+
		}
		result(i,j,k) += flux / (coords->dx(i,j)*coords->J(i,j));
		//result(i+1,j,k) -= flux / (mesh->dx(i+1,j)*mesh->J(i+1,j));

		// Flow left
		flux = gL * 0.25*(coords->J(i-1,j) + coords->J(i,j)) *(as(i-1,j,k) + as(i,j,k));
		if(i==mesh->xstart) {
			if(flux>0.0) {
				flux = 0.0;
			}
			flux += -bndry_flux(i,j,k); //-
		}
		result(i,j,k) -= flux / (coords->dx(i,j)*coords->J(i,j));
		//result(i-1,j,k) += flux / (mesh->dx(i+1,j)*mesh->J(i+1,j));
	  }

  return result;
}

void neutral_rates(const Field3D &Ne, const Field3D &Te, const Field3D &Vi,    // Plasma quantities
                      const Field3D &Nn, const Field3D &Tn, const Field3D &Vnpar, // Neutral gas
                      Field3D &S, Field3D &F, Field3D &Q, Field3D &R,  // Transfer rates
                      Field3D &Riz, Field3D &Rrc, Field3D &Rcx) {       // Rates

  // Allocate output fields
  /*
  S.allocate();
  F.allocate();
  Q.allocate();
  R.allocate();

  Riz.allocate();
  Rrc.allocate();
  Rcx.allocate();
  */
  S = 0.0;
  F = 0.0;
  Q = 0.0;
  R = 0.0;

  Riz = 0.0;
  Rrc = 0.0;
  Rcx = 0.0;

  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++)
      for(int k=0;k<mesh->LocalNz-1;k++) {
        // Integrate rates over each cell in Y
        // NOTE: This should integrate over (x,y,z)

        // Calculate cell centre (C), left (L) and right (R) values
        BoutReal Te_C = Te(i,j,k), Te_L = 0.5*(Te(i,j-1,k) + Te(i,j,k)), Te_R = 0.5*(Te(i,j,k) + Te(i,j+1,k));
        BoutReal Ne_C = Ne(i,j,k), Ne_L = 0.5*(Ne(i,j-1,k) + Ne(i,j,k)), Ne_R = 0.5*(Ne(i,j,k) + Ne(i,j+1,k));
        BoutReal Vi_C = Vi(i,j,k), Vi_L = 0.5*(Vi(i,j-1,k) + Vi(i,j,k)), Vi_R = 0.5*(Vi(i,j,k) + Vi(i,j+1,k));
        BoutReal Tn_C = Tn(i,j,k), Tn_L = 0.5*(Tn(i,j-1,k) + Tn(i,j,k)), Tn_R = 0.5*(Tn(i,j,k) + Tn(i,j+1,k));
        BoutReal Nn_C = Nn(i,j,k), Nn_L = 0.5*(Nn(i,j-1,k) + Nn(i,j,k)), Nn_R = 0.5*(Nn(i,j,k) + Nn(i,j+1,k));
        BoutReal Vn_C = Vnpar(i,j,k), Vn_L = 0.5*(Vnpar(i,j-1,k) + Vnpar(i,j,k)), Vn_R = 0.5*(Vnpar(i,j,k) + Vnpar(i,j+1,k));

        if(Ne_C < 0.) Ne_C = 0.0;
        if(Ne_L < 0.) Ne_L = 0.0;
        if(Ne_R < 0.) Ne_R = 0.0;
        if(Nn_C < 0.) Nn_C = 0.0;
        if(Nn_L < 0.) Nn_L = 0.0;
        if(Nn_R < 0.) Nn_R = 0.0;

        // Jacobian (Cross-sectional area)
		Coordinates *coords = mesh->coordinates();
        BoutReal J_C = coords->J(i,j), J_L = 0.5*(coords->J(i,j-1) + coords->J(i,j)), J_R = 0.5*(coords->J(i,j) + coords->J(i,j+1));

        ///////////////////////////////////////////
        // Charge exchange

        BoutReal R_cx_L = Ne_L*Nn_L*hydrogen.chargeExchange(Te_L*T0) * (n0 / w0);
        BoutReal R_cx_C = Ne_C*Nn_C*hydrogen.chargeExchange(Te_C*T0) * (n0 / w0);
        BoutReal R_cx_R = Ne_R*Nn_R*hydrogen.chargeExchange(Te_R*T0) * (n0 / w0);

        // Power transfer from plasma to neutrals
        // Factor of 3/2 to convert temperature to energy

        Q(i,j,k) =(3./2)* (
                                J_L * (Te_L - Tn_L)*R_cx_L
                           + 4.*J_C * (Te_C - Tn_C)*R_cx_C
                           +    J_R * (Te_R - Tn_R)*R_cx_R
                           ) / (6. * J_C);

        // Plasma-neutral friction
        F(i,j,k) =(
                        J_L * (Vi_L - Vn_L)*R_cx_L
                   + 4.*J_C * (Vi_C - Vn_C)*R_cx_C
                   +    J_R * (Vi_R - Vn_R)*R_cx_R
                   ) / (6. * J_C);

        // Cell-averaged rate
        Rcx(i,j,k) = (
                             J_L * R_cx_L
                      + 4. * J_C * R_cx_C
                      +      J_R * R_cx_R
                      ) / (6. * J_C);

        ///////////////////////////////////////
        // Recombination

        BoutReal R_rc_L  = hydrogen.recombination(Ne_L*n0, Te_L*T0)*SQ(Ne_L) * n0 / w0;
        BoutReal R_rc_C  = hydrogen.recombination(Ne_C*n0, Te_C*T0)*SQ(Ne_C) * n0 / w0;
        BoutReal R_rc_R  = hydrogen.recombination(Ne_R*n0, Te_R*T0)*SQ(Ne_R) * n0 / w0;

        // Radiated power from plasma
        // Factor of 1.09 so that recombination becomes an energy source at 5.25eV
        R(i,j,k) = (
                         J_L * (1.09*Te_L - 13.6/T0)*R_rc_L
                    + 4.*J_C * (1.09*Te_C - 13.6/T0)*R_rc_C
                    +    J_R * (1.09*Te_R - 13.6/T0)*R_rc_R
                    ) / (6. * J_C);

        // Plasma sink / neutral source
        S(i,j,k) = (
                          J_L * R_rc_L
                    + 4.* J_C * R_rc_C
                    +     J_R * R_rc_R
                    ) / (6. * J_C);

        // Transfer of ion momentum to neutrals
        F(i,j,k) += (
                           J_L * Vi_L * R_rc_L
                     + 4.* J_C * Vi_C * R_rc_C
                     +     J_R * Vi_R * R_rc_R
                     ) / (6. * J_C);

        // Transfer of ion energy to neutrals
        Q(i,j,k) +=(3./2) * (
                                  J_L * Te_L * R_rc_L
                             + 4.*J_C * Te_C * R_rc_C
                             +    J_R * Te_R * R_rc_R
                             ) / (6. * J_C);

        // Cell-averaged rate
        Rrc(i,j,k) = (
                             J_L * R_rc_L
                      + 4. * J_C * R_rc_C
                      +      J_R * R_rc_R
                      ) / (6. * J_C);

        ///////////////////////////////////////
        // Ionisation

        BoutReal R_iz_L = Ne_L*Nn_L*hydrogen.ionisation(Te_L*T0) * n0 / w0;
        BoutReal R_iz_C = Ne_C*Nn_C*hydrogen.ionisation(Te_C*T0) * n0 / w0;
        BoutReal R_iz_R = Ne_R*Nn_R*hydrogen.ionisation(Te_R*T0) * n0 / w0;

        // Neutral sink, plasma source
        S(i,j,k) -=  (
                            J_L * R_iz_L
                      + 4.* J_C * R_iz_C
                      +     J_R * R_iz_R
                      ) / (6. * J_C);

        // Transfer of neutral momentum to ions
        F(i,j,k) -= (
                            J_L * Vn_L * R_iz_L
                     + 4. * J_C * Vn_C * R_iz_C
                     +      J_R * Vn_R * R_iz_R
                     ) / (6. * J_C);

        // Transfer of neutral energy to ions
        Q(i,j,k) -= (3./2)* (
                                    J_L * Tn_L * R_iz_L
                             + 4. * J_C * Tn_C * R_iz_C
                             +      J_R * Tn_R * R_iz_R
                             ) / (6. * J_C);

        // Ionisation and electron excitation energy
        R(i,j,k) += (Eionize/T0) * (
                                            J_L * R_iz_L
                                       + 4.*J_C * R_iz_C
                                       +    J_R * R_iz_R
                                       ) / (6. * J_C);

        // Cell-averaged rate
        Riz(i,j,k) = (
                             J_L * R_iz_L
                      + 4. * J_C * R_iz_C
                      +      J_R * R_iz_R
                      ) / (6. * J_C);


      }
}
