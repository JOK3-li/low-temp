/*******************************************************************
 * Low temperature plasma model
 *
 * J. Leddy, April 2017
 *******************************************************************/

 #include "ltp.hxx"

 #include <initialprofiles.hxx>
 #include <derivs.hxx>
 #include <field_factory.hxx>
 #include <invert_parderiv.hxx>
 #include <bout/constants.hxx>
 #include <bout/assert.hxx>

int ltp::init(bool restarting)
{
	// read options
	Coordinates *coords = mesh->coordinates();
	Options *options = Options::getRoot();
	options = options->getSection("ltp");
	OPTION(options, L0,  1.0);
	OPTION(options, w0,  1e4);
	OPTION(options, n0,  1e15);
	OPTION(options, T0,  1.0);
	OPTION(options, v0,  1.0);
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
	SAVE_REPEAT(nu);
	SAVE_ONCE3(phimag,Efreq,v0);
	SAVE_ONCE4(w0,T0,n0,L0);

	/////////////////////////
	// Normalisations
	// Evolving variables
	Ne /= n0;
	Ni /= n0;
	NeE /= n0*T0;

	// Derived variables
	Ti /= T0;
	Ng /= n0;

	// physics quantities
	mu_i *= v0 / (w0 * L0 * L0);
	mu_e *= v0 / (w0 * L0 * L0);
	Dn /= w0*L0*L0;
	De /= w0*L0*L0;

	// grid
	coords->dx /= L0;
	coords->dz /= L0;
	// coords->J = coords->dx;

	// source
	phimag /= v0;
	Efreq *= 2. * PI / w0;

	///////////////////////////////////////////
	// Initialise potential
	phi = 0.0;
	phiext = 0.0;
	total_phi = 0.0;
	E = NeE / Ne;

	return 0;
}

int ltp::rhs(BoutReal t) {
	mesh->communicate(Ne,Ni,NeE,Ng);//,Fex,Fix,Fez,Fiz);

	Ne = floor(Ne,1e-5);
	Ni = floor(Ni,1e-5);
	NeE = floor(Ne,1e-5);
	Ng = floor(Ng,1e-5);

	// invert Poisson equation
	TRACE("Potential inversion");
	phi = phiSolver->solve(Ne-Ni);
	mesh->communicate(phi);
	phiext = phimag * linear * sin(Efreq * t);
	total_phi = phi + phiext;
	mesh->communicate(total_phi);

	// calculate energy and temperature
	E = NeE / Ne;
	Te = 2.*E/3.;
	Field3D zero = 0.0;

	// Sheath boundary conditions
	// BoundaryRegion *bndry;
	// for(bndry->first(); !bndry->isDone(); bndry->nextX()) {
	// 	BoutReal dNedx = DDX(Ne)[bndry->x][bndry->y]
	// 	Ne[bndry->x - bndry->bx][bndry->y] = 0.25 * sqrt(Te) * dNedx
	// }

	// Gas continuity equation
	TRACE("Neutral density equation");
	ddt(Ng) = 0.0;
	if(evolve_Ng) {
		// Electron density
		ddt(Ng) += -Div_Perp_Lap_FV(Dn,Ng,zero);
	}

	// Electron continuity equation
	TRACE("Electron density equation");
	ddt(Ne) = 0.0;
	if(evolve_Ne) {
		// Electron density
		ddt(Ne) += -Div_Perp_Lap_FV(-Ne*mu_e,total_phi,0.25*Ne*sqrt(Te));
		ddt(Ne) += -Div_Perp_Lap_FV(Dn,Ne,zero);
	}

	// Ion continuity equation
	TRACE("Ion density equation");
	ddt(Ni) = 0.0;
	if(evolve_Ni) {
		// Ion density
		ddt(Ni) += -Div_Perp_Lap_FV(Ni*mu_i,total_phi,0.25*Ni*sqrt(Ti));
		ddt(Ni) += -Div_Perp_Lap_FV(Dn,Ni,zero);
	}

	// Electron energy density equation
	TRACE("Energy density equation");
	ddt(NeE) = 0.0;
	if(evolve_NeE) {
		// Energy density
		ddt(NeE) += Div_Perp_Lap_FV(5./3.*NeE*mu_e,total_phi,5./12.*NeE*sqrt(Ti));
		ddt(NeE) += -Div_Perp_Lap_FV(5./3.*E*Dn,Ne,zero);
		ddt(NeE) += -Div_Perp_Lap_FV(5./3.*Ne*De,E,zero);
		// ddt(NeE) += -Ne*mu_e*DDX(total_phi)*DDX(total_phi);

		// Energy source from interactions
		nu = 0.0;
		for(int i=mesh->xstart;i<=mesh->xend;i++) {
			for(int j=mesh->ystart;j<=mesh->yend;j++) {
				for(int k=0;k<mesh->LocalNz;k++) {
					nu(i,j,k) = Ng(i,j,k) * (2.336e-14*pow(Te(i,j,k)*T0,1.609) *
						exp(0.0618*pow(log(Te(i,j,k)*T0),2) - 0.1171*pow(log(Te(i,j,k)*T0),3))) * n0 / w0;
				}
			}
		}
		ddt(NeE) += 3. * Me/Mn * nu * Ne * (Te - Ti);
	}

	TRACE("Calculate neutral rates");
	Field3D Riz, Rrc, Rcx;
	neutral_rates(Ne, Te, Ng, Ti, S, F, Q, Rp, Riz, Rrc, Rcx);

	/////////////////////////////////////////////////////
	// Neutral density source
	TRACE("Neutral density");
	if(evolve_Ng) {
		ddt(Ng) += S; // Sink of plasma density
	}
	if(evolve_Ni) {
		ddt(Ni) -= S; // Sink of plasma density
	}
	if(evolve_Ne) {
		ddt(Ne) -= S; // Sink of plasma density
	}

	/////////////////////////////////////////////////////
	// Neutral pressure
	TRACE("Neutral energy density");
	if(evolve_NeE) {
		ddt(NeE) -= Q + Rp; // originally Pe (NeTe), need to figure out how to make this for NeE -> (2./3) * (Q + Rp) originally
	}

	return 0;
}

// Div ( a * Grad_perp(f) )
const Field3D ltp::Div_Perp_Lap_FV(const Field3D &a, const Field3D &f, const Field3D &bndry_flux) {

	Coordinates *coords = mesh->coordinates();
	Field3D result = 0.0;

	//////////////////////////////////////////
	// d/dx(a * d/dx(f))
	//
	//
	//			f_i				f_j				f_k
	//  		a_i				a_j				a_k
	//				  df_ij/dx		  df_jk/dx
	//					F_ij			F_jk
	//		  dF_ij/dx				  dF_ij/dx
	//   |-------o-------|-------o-------|-------o-------|
	//			 i				 j				 k
	//
	//

	Field3D fs = f;
	Field3D as = a;

	BoutReal Globaldx = mesh->GlobalX(1) - mesh->GlobalX(0);

	for(int i=mesh->xstart;i<=mesh->xend;i++) {
		for(int j=mesh->ystart;j<=mesh->yend;j++) {
			for(int k=0;k<mesh->LocalNz;k++) {

				// Calculate gradients on cell faces
				BoutReal gR = (fs(i+1,j,k) - fs(i,j,k)) / (coords->dx(i+1,j) + coords->dx(i,j));
				BoutReal gL = (fs(i,j,k) - fs(i-1,j,k)) / (coords->dx(i-1,j) + coords->dx(i,j));

				// Flow out of right face
				BoutReal fluxR = gR * 0.5 * (as(i+1,j,k) + as(i,j,k));
				if(mesh->GlobalX(i)>1-Globaldx) {
					// output<<mesh->GlobalX(i)<<"  "<<i<<"\n";
					if(fluxR<0.0){
						fluxR = 0.0;
					}
					fluxR += bndry_flux(i,j,k);
				}

				// Flow into left face
				BoutReal fluxL = gL * 0.5 * (as(i-1,j,k) + as(i,j,k));
				if(mesh->GlobalX(i)<Globaldx) {
					// output<<mesh->GlobalX(i)<<"  "<<i<<"\n";
					if(fluxL>0.0) {
						fluxL = 0.0;
					}
					fluxL += -bndry_flux(i,j,k);
				}

				// Divergence of flux on cell centre
				result(i,j,k) = (fluxR - fluxL) / coords->dx(i,j);
			}
		}
	}

	return result;
}

// Div ( a * Grad_perp(f) )
const Field3D ltp::Div_Perp_Lap_FV4(const Field3D &a, const Field3D &f, const Field3D &bndry_flux) {

	Coordinates *coords = mesh->coordinates();
	Field3D result = 0.0;

	//////////////////////////////////////////
	// d/dx(a * d/dx(f))
	//
	//
	//			f_i				f_j				f_k
	//  		a_i				a_j				a_k
	//				  df_ij/dx		  df_jk/dx
	//					F_ij			F_jk
	//		  dF_ij/dx				  dF_ij/dx
	//   |-------o-------|-------o-------|-------o-------|
	//			 i				 j				 k
	//
	//

	Field3D fs = f;
	Field3D as = a;

	BoutReal Globaldx = mesh->GlobalX(1) - mesh->GlobalX(0);

	for(int i=mesh->xstart;i<=mesh->xend;i++) {
		for(int j=mesh->ystart;j<=mesh->yend;j++) {
			for(int k=0;k<mesh->LocalNz;k++) {

				// Calculate gradients on cell faces
				BoutReal gR = 2. * ((fs(i-1,j,k)+fs(i,j,k))/24. - 2./3.*fs(i,j,k) + 2./3.*fs(i+1,j,k) - (fs(i+1,j,k)+fs(i+2,j,k))/24.) / coords->dx(i,j);
				BoutReal gL = 2. * ((fs(i-2,j,k)+fs(i-1,j,k))/24. - 2./3.*fs(i-1,j,k) + 2./3.*fs(i,j,k) - (fs(i,j,k)+fs(i+1,j,k))/24.) / coords->dx(i,j);

				// Flow out of right face
				BoutReal fluxR = gR * 0.5 * (as(i+1,j,k) + as(i,j,k));
				if(mesh->GlobalX(i)>1-Globaldx) {
					// output<<mesh->GlobalX(i)<<"  "<<i<<"\n";
					if(fluxR<0.0){
						fluxR = 0.0;
					}
					fluxR += bndry_flux(i,j,k);
				}

				// Flow into left face
				BoutReal fluxL = gL * 0.5 * (as(i-1,j,k) + as(i,j,k));
				if(mesh->GlobalX(i)<Globaldx) {
					// output<<mesh->GlobalX(i)<<"  "<<i<<"\n";
					if(fluxL>0.0) {
						fluxL = 0.0;
					}
					fluxL += -bndry_flux(i,j,k);
				}

				// Divergence of flux on cell centre
				result(i,j,k) = (fluxR - fluxL) / coords->dx(i,j);
			}
		}
	}

	return result;
}

// const Field3D Div_n_bxGrad_f_B_XPPM(const Field3D &n_in, const Field3D &f_in) {
// 	Field3D result = 0;
//
// 	//////////////////////////////////////////
// 	// X-Z advection.
// 	//
// 	//             Z
// 	//             |
// 	//
// 	//    fmp --- vU --- fpp
// 	//     |      nU      |
// 	//     |               |
// 	//    vL nL        nR vR    -> X
// 	//     |               |
// 	//     |      nD       |
// 	//    fmm --- vD --- fpm
// 	//
//
// 	Field3D n = n_in;  // Done in orthogonal X-Z coordinates
// 	Field3D f = f_in;
//
// 	for(int i=mesh->xstart;i<=mesh->xend;i++)
//       for(int j=mesh->ystart;j<=mesh->yend;j++)
//         for(int k=0;k<mesh->ngz-1;k++) {
// 		int kp = (k+1) % (mesh->ngz-1);
// 		int kpp = (kp+1) % (mesh->ngz-1);
// 		int km = (k-1+mesh->ngz-1) % (mesh->ngz-1);
// 		int kmm = (km-1+mesh->ngz-1) % (mesh->ngz-1);
//
// 		// 1) Interpolate stream function f onto corners fmp, fpp, fpm
//
// 		BoutReal fmm = 0.25*(f(i,j,k) + f(i-1,j,k) + f(i,j,km) + f(i-1,j,km));
// 		BoutReal fmp = 0.25*(f(i,j,k) + f(i,j,kp) + f(i-1,j,k) + f(i-1,j,kp)); // 2nd order accurate
// 		BoutReal fpp = 0.25*(f(i,j,k) + f(i,j,kp) + f(i+1,j,k) + f(i+1,j,kp));
// 		BoutReal fpm = 0.25*(f(i,j,k) + f(i+1,j,k) + f(i,j,km) + f(i+1,j,km));
//
// 		// 2) Calculate velocities on cell faces
//
// 		BoutReal vU = mesh->J(i,j)*(fmp - fpp)/mesh->dx(i,j); // -J*df/dx
// 		BoutReal vD = mesh->J(i,j)*(fmm - fpm)/mesh->dx(i,j); // -J*df/dx
//
// 		BoutReal vR = 0.5*(mesh->J(i,j)+mesh->J(i+1,j))*(fpp - fpm)/mesh->dz; // J*df/dz
// 		BoutReal vL = 0.5*(mesh->J(i,j)+mesh->J(i-1,j))*(fmp - fmm)/mesh->dz; // J*df/dz
//
// 	        //output.write("NEW: (%d,%d,%d) : (%e/%e, %e/%e)\n", i,j,k,vL,vR, vU,vD);
//
// 		// 3) Calculate n on the cell faces. The sign of the
// 		//    velocity determines which side is used.
//
// 		// X direction
// 		Stencil1D s;
// 		s.c  = n(i,  j,k);
// 		s.m  = n(i-1,j,k);
// 		s.mm = n(i-2,j,k);
// 		s.p  = n(i+1,j,k);
// 		s.pp = n(i+2,j,k);
//
// 		//Upwind(s, mesh->dx(i,j));
// 		//XPPM(s, mesh->dx(i,j));
// 		Fromm(s, mesh->dx(i,j));
//
// 		// Right side
// 		if((i==mesh->xend) && (mesh->lastX())) {
// 			// At right boundary in X
//
// 			if(bndry_flux) {
// 				BoutReal flux;
// 				if(vR > 0.0) {
// 					// Flux to boundary
// 					flux = vR * s.R;
// 				}else {
// 					// Flux in from boundary
// 					flux = vR * 0.5*(n(i+1,j,k) + n(i,j,k));
// 				}
// 				result(i,j,k)   += flux / (mesh->dx(i,j) * mesh->J(i,j));
// 				result(i+1,j,k) -= flux / (mesh->dx(i+1,j) * mesh->J(i+1,j));
// 			}
// 		}else {
// 			// Not at a boundary
// 			if(vR > 0.0) {
// 				// Flux out into next cell
// 				BoutReal flux = vR * s.R;
// 				result(i,j,k)   += flux / (mesh->dx(i,j) * mesh->J(i,j));
// 				result(i+1,j,k) -= flux / (mesh->dx(i+1,j) * mesh->J(i+1,j));
//
// 				//if(i==mesh->xend)
// 				//  output.write("Setting flux (%d,%d) : %e\n", j,k,result(i+1,j,k));
// 			}
// 		}
//
//         // Left side
//
// 		if((i==mesh->xstart) && (mesh->firstX())) {
// 			// At left boundary in X
// 			if(bndry_flux) {
// 				BoutReal flux;
// 				if(vL < 0.0) {
// 					// Flux to boundary
// 					flux = vL * s.L;
// 				}else {
// 					// Flux in from boundary
// 					flux = vL * 0.5*(n(i-1,j,k) + n(i,j,k));
// 				}
// 				result(i,j,k)   -= flux / (mesh->dx(i,j) * mesh->J(i,j));
// 				result(i-1,j,k) += flux / (mesh->dx(i-1,j) * mesh->J(i-1,j));
// 			}
// 		}else {
// 			// Not at a boundary
//
// 			if(vL < 0.0) {
// 				BoutReal flux = vL * s.L;
// 				result(i,j,k)   -= flux / (mesh->dx(i,j) * mesh->J(i,j));
// 				result(i-1,j,k) += flux / (mesh->dx(i-1,j) * mesh->J(i-1,j));
// 			}
// 		}
//
// 		/// NOTE: Need to communicate fluxes
//
// 		// Z direction
// 		s.m  = n(i,j,km);
// 		s.mm = n(i,j,kmm);
// 		s.p  = n(i,j,kp);
// 		s.pp = n(i,j,kpp);
//
// 		//Upwind(s, mesh->dz);
// 		//XPPM(s, mesh->dz);
// 		Fromm(s, mesh->dz);
//
// 		if(vU > 0.0) {
// 			BoutReal flux = vU * s.R / (mesh->J(i,j)*mesh->dz);
// 			result(i,j,k)   += flux;
// 			result(i,j,kp)  -= flux;
// 		}
// 		if(vD < 0.0) {
// 			BoutReal flux = vD * s.L / (mesh->J(i,j)*mesh->dz);
// 			result(i,j,k)   -= flux;
// 			result(i,j,km)  += flux;
// 		}
// 	}
//
// 	return result;
// }

BOUTMAIN(ltp);
