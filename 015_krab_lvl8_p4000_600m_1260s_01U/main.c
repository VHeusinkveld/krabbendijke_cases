/** Include required libraries */ 
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "grid/octree.h" 		// For 3D
#include "view.h"			// For bview
#include "navier-stokes/centered.h"     // Navier stokes 
#include "tracer.h"			// Tracers
#include "diffusion.h"			// Diffusion 

/** Global variables */
int minlevel, maxlevel;         	// Grid depths
double meps, eps;			// Maximum error and error in u fields
double TEND = 1260.;

char sim_ID[] = "krab";		        // Simulation identifier
char sim_var[] = "None";  		// Notes if a variable is varied over runs

#include "physics.h"			// Physics of the simulation 
#include "fan.h"			// Include a fan
#include "diagnostics.h"		// Perform diagnostics

/** Initialisation */
int main() {	
    minlevel = 4;
    maxlevel = 8;

    L0 = 600.;
    X0 = Y0 = Z0 = 0.;

    // Possibility to run for variable changes
    for(rot.theta=97.*M_PI/180.; rot.theta<=110.*M_PI/180.; rot.theta+=100.*M_PI/180.) {
        init_grid(1<<6);
	a = av; 

        foreach_dimension() {
	    u.x.refine = refine_linear;  		// Momentum conserved 
	}

	fan.prolongation = fraction_refine;		// Fan is a volume fraction
	p.refine = p.prolongation = refine_linear;
	b.gradient = minmod2; 				// Flux limiter 

  	meps = 10.;					// Maximum adaptivity criterion
	DT = 10E-5;					// For poisson solver 
        TOLERANCE=10E-6;				// For poisson solver 
	CFL = 0.5;					// CFL condition

	sim_dir_create();				// Create relevant dir's
	out.sim_i++;					// Simulation iteration
 
    	run();						// Start simulation 

    }
}

/** Initialisation */
event init(t=0) {
    rot.fan = true;		// Yes we want a fan
    rot.rotate = true;		// If we want it to rotate 
    rot.start = 60;
    rot.stop = 1140;

    rot.phi = 0;		// Reset for different runs
    eps = .5;
    
    init_physics();

    if(rot.fan) {
        init_rotor();
	rotor_coord();
    }
        
    while(adapt_wavelet((scalar *){u,b},(double []){eps,eps,eps,0.38*9.81/273},maxlevel,minlevel).nf) {
	foreach() {
	    b[] = STRAT(y);
            u.x[] = WIND(y);
	}
	rotor_coord();
    }
 

}

/** Return to standard tolerances and DTs for poisson solver */ 
event init_change(i=10) {
    TOLERANCE=10E-3;
    DT = 0.5;
}

/** Adaptivity */
event adapt(i++) {
    adapt_wavelet((scalar *){fan,u,b},(double []){0.01,eps,eps,eps,0.38*9.81/273},maxlevel,minlevel);
}

/** Progress event */
event progress(t+=5) {
    fprintf(stderr, "i=%d t=%g p=%d u=%d b=%d \n", i, t, mgp.i, mgu.i, mgb.i);
}

event dumpfields(t=120; t+=120) {
    char nameDump[90];
    snprintf(nameDump, 90, "./%s/fielddump", out.dir);
    dump(file = nameDump, list = all);
}

/** End the simulation */
event end(t=TEND) {
}
