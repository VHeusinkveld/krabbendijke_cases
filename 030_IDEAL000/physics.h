#if dimension == 3
	#include "SGS.h"
#endif 

#define CP 1005.	// C_p for air 
#define gCONST 9.81	// Gravitational constant
#define TREF 273.	// Kelvin
#define INVERSION .5 	// Kelvin per meter
#define karman 0.4      // von Karman constant 

#define roughY0u 0.1    // roughness wind length 
#define roughY0h 0.1     // roughness heat length 
#define WIND(s)  (-1.) 
#define Lambda 1.
#define STRAT(s) gCONST*(INVERSION)*s/TREF //+ gCONST/CP


scalar b[];
scalar * tracers = {b};

double crho = 1.;	

face vector av[]; 
struct sCase def;

struct sCase {
	double wind;
	double wphi;
};	

void init_physics(){
 	def.wind = WIND(1.);
        def.wphi = 0.;

	b.nodump = false; // TODO

        u.n[bottom] = dirichlet(0.);
        u.t[bottom] = dirichlet(WIND(s)); 
        u.n[top] = dirichlet(0);
        u.t[top] = dirichlet(WIND(s));

        periodic (left);      
	
	b[bottom] = dirichlet(STRAT(y));
	b[top] = dirichlet(STRAT(y));

	#if dimension == 3
		u.r[bottom] = dirichlet(WIND(s));
                u.r[top] = dirichlet(WIND(s));
	        u.t[bottom] = dirichlet(0.); 
		u.t[top] = dirichlet(0.); 

		periodic(front);
	#endif  
	foreach() {
		b[] = STRAT(y);
 	        u.x[] = WIND(y);
	}
}

/* Gravity forcing */
event acceleration(i++){
	foreach_face(y){
		av.y[] = (b[] + b[0,-1])/2.;
	}
}

event inflow(i++){
    double sides = 25;
    double relaxtime = dt/25.;
    foreach(){
	if((x < sides || x > L0-sides || 
	    z < sides || z > L0-sides ||
 	    y < sides || y > L0-sides )) {
	    u.x[] = u.x[] + (WIND(s)-u.x[])*relaxtime;
 	    b[] = b[] + (STRAT(y) - b[])*relaxtime;
	    u.y[] = u.y[] - u.y[]*relaxtime;
            u.z[] = u.z[] - u.z[]*relaxtime;
	}
    }
}

mgstats mgb;
/* Diffusion */
event tracer_diffusion(i++){
    mgb = diffusion(b, dt, mu);
}
