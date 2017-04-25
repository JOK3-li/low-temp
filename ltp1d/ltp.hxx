
class ltp;

#ifndef __LTP_H__
#define __LTP_H__

#include <bout/physicsmodel.hxx>

#include <invert_laplace.hxx>
#include <bout/constants.hxx>

#include "radiation.hxx"

class ltp : public PhysicsModel {
public:
	virtual ~ltp() {}
protected:
	int init(bool restarting);
	int rhs(BoutReal t);

	int precon(BoutReal t, BoutReal gamma, BoutReal delta);
private:
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
	BoutReal mu_e, mu_i, Dn, De, Ti;
	BoutReal w0, L0, n0, v0, T0;	// normalisations
	int ZZ;
	BoutReal phimag, Efreq;  // imposed field
	bool evolve_Ne, evolve_Ni, evolve_NeE, evolve_Ng;

	// neutral interaction
	UpdatedRadiatedPower hydrogen; // Atomic rates (H.Willett)
	Field3D S,F,Q,Rp,nu;
	Field3D Riz,Rrc,Rcx;
	BoutReal Eionize;

	const Field3D Div_Perp_Lap_FV(const Field3D &a, const Field3D &f, const Field3D &bndry_flux);
	// const Field3D Div_par_FV(const Field3D &f, const Field3D &v);
	const Field3D Div_Perp_Lap_FV4(const Field3D &a, const Field3D &f, const Field3D &bndry_flux);
	void neutral_rates(const Field3D &Ne, const Field3D &Te, const Field3D &Nn, const Field3D &Tn,
		Field3D &S, Field3D &F, Field3D &Q, Field3D &R, Field3D &Riz, Field3D &Rrc, Field3D &Rcx);

};

/// Fundamental constants
const BoutReal Mn = 6.6335209e-26; // kg // 3.3435e-27;//
const BoutReal Me = 9.10938356e-31; // kg

#endif // __LTP_H__
