#if dimension == 3
	#include "SGS.h"
#endif 

#define CP 1005.	// C_p for air 
#define gCONST 9.81	// Gravitational constant
#define TREF 273.	// Kelvin
#define INVERSION .2 	// Kelvin per meter
#define karman 0.4      // von Karman constant 

#define roughY0u 0.1    // roughness wind length 
#define roughY0h 0.01     // roughness heat length 
#define WIND(s) -max((0.30*log(2.*(s-roughY0u)+1.)),0.)   // log Wind profile 

#define QFLX 0. 	// 0 (0.001 = 20wm2)
#define BSURF ((b[0,1]-b[]*lut2[level])/(1.-lut2[level]))  // log estimate of surface b
//linear (1.5*b[] - 0.5*b[0,1])   // Estimation of surface b
#define GFLX (-Lambda*(BSURF - bd))
double Lambda = 0.005, bd = 0.;   // Grass coupling
#define STRAT(s) gCONST/TREF*(log(30*s + a1 + 1.) - log(a1 + 1.)) + (QFLX/Lambda + bd)
double a1 = 0.;


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
 	def.wind = 1.;
        def.wphi = 0.;

	b.nodump = false; // TODO

	if(def.wind == 0.){    	
	    u.n[bottom] = dirichlet(0.);
	    u.t[bottom] = dirichlet(0.);
	    u.n[top] = dirichlet(0.);
 	    u.t[top] = neumann(0.);

            periodic (left);
	
	} else if(fabs(def.wind) > 0.) {
	    u.n[bottom] = dirichlet(0.);
	    u.t[bottom] = dirichlet(0.);
	    u.n[top] = dirichlet(0);
 	    u.t[top] = dirichlet(WIND(y));

            periodic (left);
		
	}        
	
	b[bottom] = BSURF;
	b[top] = dirichlet(STRAT(y));

	#if dimension == 3
		u.r[bottom] = dirichlet(0.); 
		u.r[top] = neumann(0.); 
		
		Evis[bottom] = dirichlet(0.); // Flux is explicitly calculated
        	Evis[top] = dirichlet(0.);

		periodic(front);
	#endif  
	foreach() {
		b[] = STRAT(y);
		if(fabs(def.wind) > 0.) {
       	            u.x[] = WIND(y);
		}
	}
}

double lut[20];
double lut2[20];
event init (t = 0) {
    for (int m = 0; m <= 19; m++) {
	double d  = (L0/((double)(1 << m)))/roughY0u;
	double d2 = (L0/((double)(1 << m)))/roughY0h;

        if (m==0) {
	    lut[0] = 0;
	    lut2[0] = 0;
	} else {
            lut[m] = sq(karman/log(L0/d - 1.));
	}
        lut2[m] = (log(4.*d2) - 1.)/(log(d2) - 1.);
    }
}

#define dvdt(u) (sign(u)*lut[level]*sq(u)/Delta)

event law_of_the_wall(i++){
    double ui; //scratch for the mid-point-value estimate 
    foreach_boundary(bottom){
        ui = u.x[] - (dt*dvdt(u.x[])/2.); 
        u.x[] -= dt*dvdt(ui);
        ui = u.z[] - (dt*dvdt(u.z[])/2.);
        u.z[] -= dt*dvdt(ui);
    }
}


/* Gravity forcing */
event acceleration(i++){
	foreach_face(y){
		av.y[] = (b[] + b[0,-1])/2.;
	}
}

event inflow(i++){
    double sides = 0.08;
    double relaxtime = dt/40.;
    foreach(){
	if((x < sides*L0 || x > (1-sides)*L0 || 
	    z < sides*L0 || z > (1-sides)*L0 ||
 	    y > (1-2*sides)*L0 )) {
	    u.x[] = u.x[] + (WIND(y)-u.x[])*relaxtime;
	    u.y[] = u.y[] - u.y[]*relaxtime;
            u.z[] = u.z[] - u.z[]*relaxtime;
	}
	if((x < sides*L0 || x > (1-1.5*sides)*L0 || 
	    z < sides*L0 || z > (1-1.5*sides)*L0 ||
	    y > (1-3*sides)*L0 )) {
 	    b[] = b[] + (STRAT(y) - b[])*relaxtime/1.5;
	}

    }
}

mgstats mgb;
//face vector muz; This is taken care of by Evis
/* Diffusion */
event tracer_diffusion(i++){
    scalar r[];
    /*
    foreach_face(){
	if(y < Y0 + 1E-10){
            muz.x[] = 1.*mu.x[];
   	} else {
	    muz.x[] = 1.*mu.x[];
        }
    }*/
    
    foreach() {
        r[] = 0;
        if (y < Delta)
            r[] = (QFLX + GFLX)/sq(Delta); // div needed as normalization 
    }
    /*
    double flx = 0, bt = 0;
    double fctr = CP*TREF/gCONST;
    foreach_boundary(bottom reduction(+:flx) reduction(+:bt)) {
        flx = flx + (QFLX + GFLX) * sq(Delta);
         bt = bt + BSURF * sq(Delta);
    }
    bt = bt/sq(L0);
    flx = flx/sq(L0);
    fprintf(stderr, "%g %g %g %d\n", t, fctr*flx, fctr*bt/CP, i);  
    */
    mgb = diffusion(b, dt, mu, r = r);
}
