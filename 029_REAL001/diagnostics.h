#include "utils.h"
#if dimension == 3
	#include "lambda2.h"
#endif
#include "output_slices.h"
#include "profile5c.h"
#include "equi_data.h"

/** Define variables and structures to do: diagnostics, ouput data, output movies. */

struct sDiag dia; 			// Diagnostics

struct sDiag {
	double Ekin;			// Total kinetic energy
	double EkinOld;			// Track changes in kin energy 
	double WdoneOld;		// Track changes in work done 
	double rotVol;			// Diagnosed rotor volum
	double bE;			// Buoyancy energy
	double bEold;			// Track changes in buoyancy energy
	double diss;			// Diagnose dissipation
};

struct sEquiDiag {
    int level;			// Level for which diagnostics should be done
    int ii; 			// Keep track of how many additions are done
    double dtDiag;
    double startDiag;
    double endDiag;
    double dtOutput;
};

struct sOutput {
    double dtDiag;
    double dtVisual;
    double dtSlices;
    double dtProfile;
    double startAve;
    double dtAve;
    char main_dir[12];
    char dir[30];
    char dir_profiles[60];
    char dir_slices[60];
    char dir_equifields[60];
    char dir_strvel[60];
    char dir_refvel[60];
    char dir_diffbins[60];
    char dir_dts[60];
    int sim_i;
};
	
struct sbViewSettings {
    double phi; 			// Phi for 3D bview movie
    double theta;			// Theta for 3D bview movie
    double sphi; 			// Polar angle for sliced image RED
    double stheta;			// Azimuthal angle for sliced image RED
};

/** Initialize structures */
struct sOutput out = {.dtDiag = 1., .dtVisual=1., .dtSlices=5., .dtProfile=60., .main_dir="results", .sim_i=0};

struct sEquiDiag ediag = {.level = 5, .ii = 0, .startDiag = 0., .dtDiag = 0., .dtOutput = 0.};

struct sbViewSettings bvsets = {.phi=0., .theta=0., .sphi=0., .stheta=0.};

event init(i = 0){
	bvsets.phi = 0.;
	bvsets.theta = -M_PI/6.;
	bvsets.sphi = 0.;
	bvsets.stheta = 0.;
}

/** Diagnosing: kinetic energy, diagnosed rotor volume, buoyancy energy, ammount of cells used.*/
event diagnostics (t+=out.dtDiag){
	int n = 0.;
	scalar ekin[]; 		// Kinetic energy field 
	double tempVol = 0.;    // Temp volume 
	double tempEkin = 0.;   // Temp kinetic energy
	double tempDiss = 0.;   // Temp dissipation 
	double maxVel = 0.;     // Maximum velocity in fan
	double bEnergy = 0.;    // Buoyant energy
		
	/** Loop over cells to get diagnostics */ 
	foreach(reduction(+:n) reduction(+:tempVol) reduction(+:tempEkin) 
		reduction(max:maxVel) reduction(+:bEnergy) reduction(+:tempDiss)) {

		tempVol += dv()*fan[];
		if(y + Delta/2. <= rot.y0){
		     bEnergy += dv()*y*(b[] - STRAT(y));
		}
		foreach_dimension() {
			ekin[] += sq(u.x[]);
		}
		maxVel = max(maxVel, sq(ekin[]));
		ekin[] *= 0.5*rho[]*dv();	
		tempEkin += ekin[];
		n++;
	}
	/** Assign values to respective global sturcture vars */ 
	dia.diss = 1.*tempDiss;
	dia.bE = 1.*bEnergy;
	dia.rotVol = 1.*tempVol;
	dia.Ekin = 1.*tempEkin;
	
	/** Check if fan volume is within twenty percent of definition */
	// This piece of code got semi redundant, fan also performs if Vdia == 0
	//if(fabs(dia.rotVol/rot.V - 1) > 0.20){
	//	fprintf(stderr, "ERROR Check fan volume, V=%g, Vr=%g\n",rot.V, dia.rotVol);
	//}
	
	if (pid() == 0){
	/** Write away simulation data and case setup for main thread */ 
	char nameOut[90];
	char nameCase[90];
     	snprintf(nameOut, 90, "./%s/output", out.dir);
     	snprintf(nameCase, 90, "./%s/case", out.dir);
	static FILE * fpout = fopen(nameOut, "w");
	static FILE * fpca = fopen(nameCase, "w");

	if(t==0.){
		fprintf(fpca,"L0\tinversion\thubU\tTref\tLambda\txr\tyr\tzr\ttheta\tphit\tr\tW\tP\tcu\trampT\tmaxlvl\tminlvl\teps\n");
		fprintf(fpca, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%g\n", 
				L0,TREF/gCONST*STRAT(rot.y0)-STRAT(1.5), WIND(rot.y0),TREF, Lambda, rot.x0, rot.y0, rot.z0, rot.theta, rot.phit, rot.R, rot.W, rot.P, rot.cu, rot.rampT, maxlevel, minlevel, eps);
		
	        fprintf(stderr,"n\tred\tEkin\tWork\tbE\n");
		fprintf(fpout,"i\tt\tn\tred\tEkin\tWork\tbE\n");
	}
	fprintf(fpout, "%d\t%g\t%d\t%g\t%g\t%g\t%g\n",
		i,t,n,(double)((1<<(maxlevel*3))/n),dia.Ekin,rot.Work, dia.bE);
	
	fprintf(stderr, "%d\t%g\t%g\t%g\t%g\n",n,(double)((1<<(maxlevel*dimension))/n),dia.Ekin,rot.Work,dia.bE);
	
	fflush(fpout);
	fflush(fpca);	

	}

	dia.EkinOld = 1.*dia.Ekin;
	dia.WdoneOld = 1.*rot.Work;
	dia.bEold = 1.*dia.bE;
}

#if dimension == 3
/** Profiles in height */
event profiles(t += out.dtProfile) {
	char nameProf[90];
	snprintf(nameProf, 90, "./%s/t=%05g", out.dir_profiles, t);
	field_profile((scalar *){b, u}, nameProf, n=100);
}

/** Average in time */
event equidiags(t = ediag.startDiag; t += ediag.dtDiag) {
	ediag.ii = equi_diag(b, ediag.level, ediag.ii);
}

/** Eauidistant field output */
event equioutputs(t = (ediag.dtOutput + ediag.startDiag); t += ediag.dtOutput) {
    char nameEquif[90];
    snprintf(nameEquif, 90, "%s/%st=%05g", out.dir_equifields, "equifield", t);
    equi_output_binary(b, nameEquif, ediag.level, ediag.ii);

    ediag.ii = 0;
	
    free(equifield);				// We dont want memory leaks 
    equifield = NULL;				// Reset equifield pointer
}
/*
event tempbinning(t+=1) {
    double Tmax =  2;
    double Tmin = -2;

    int ntot = 51;

    double Tstep = (fabs(Tmax) + fabs(Tmin))/ntot;

    double binsp[ntot];
    memset(&binsp, 0., ntot*sizeof(double)); 
    double binnorm = 0;

#if _OPENMP
    foreach_leaf() {
#else 
    foreach(reduction(+:binnorm)) {
#endif
	double Tdif = TREF/gCONST*(b[] - STRAT(y)); // TODO think about removing linear term
	int place = round((Tdif + fabs(Tmin))/Tstep);

	place = place < 0 ? 0 : place;
	place = place > ntot ? ntot : place;	

	binsp[place] += dv();
	binnorm += dv();
    }
#if _MPI
    if(pid() == 0) {
	MPI_Reduce(MPI_IN_PLACE, &binsp[0], ntot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(&binsp[0], NULL, ntot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif     
    if(pid() == 0) {
	char nameDiffbins[90];
    	snprintf(nameDiffbins, 90, "%st=%05g", out.dir_diffbins, t);
        FILE * fpstr = fopen(nameDiffbins, "w");
       
	fprintf(fpstr, "dT\t Prob\n");
       
        for(int ip=0; ip<ntot; ip+=1) {
            binsp[ip] /= binnorm;
	    fprintf(fpstr, "%g\t %g\n", Tmin + ip*Tstep, binsp[ip]);            
	}
	
	fclose(fpstr);
    }

}
*/
event fanvelocity(t+=1) {
    if(pid() == 0) {
    char nameStrvel[90];
    snprintf(nameStrvel, 90, "%st=%05g", out.dir_strvel, t);
    FILE * fpstr = fopen(nameStrvel, "w");
    fprintf(fpstr, "R,v,vx,vy,vz\n");
        
    double length = 200;
    int ntot = 200;
    double vels1[ntot];	

    double xf0 = rot.x0;
    double yf0 = rot.y0;
    double zf0 = rot.z0;

    for(int n; n < ntot; n++) {
	double dist = length*n/ntot - 100;
	double phi = dist < 0 ? M_PI/360. : 0.;
	double theta = 97*M_PI/180;

	double nfx1 = sin(theta)*cos(phi);
	double nfz1 = sin(theta)*sin(phi);
	double nfy1 = cos(theta);

	double xx1 = xf0 + dist*nfx1; 
	double yy1 = yf0 + dist*nfz1; 
	double zz1 = zf0 + dist*nfy1; 
        // double rr = sqrt(sq(xx-xf0)+sq(yy-yf0)+sq(zz-zf0));	

	double valx1 = interpolate(u.x, xx1, yy1, zz1);
	double valy1 = interpolate(u.y, xx1, yy1, zz1);
	double valz1 = interpolate(u.z, xx1, yy1, zz1);
         
        vels1[n] = sqrt(sq(valx1) + sq(valy1) + sq(valz1));
	fprintf(fpstr, "%g,%g,%g,%g,%g\n", dist, vels1[n], valx1, valy1, valz1);     
    }  
    fclose(fpstr); 
    }
}

event refvelocity(t+=1) {
    if(pid() == 0) {
    char nameRefvel[90];
    snprintf(nameRefvel, 90, "%st=%05g", out.dir_refvel, t);
    FILE * fpstr = fopen(nameRefvel, "w");
    fprintf(fpstr, "x,v,vx,vy,vz\n");
        
    double length = L0;
    int ntot = L0/2;
    double vels1[ntot];	

    double xf0 = 0;
    double yf0 = 3;
    double zf0 = L0/2;

    for(int n; n < ntot; n++) {
	double dist = length*n/ntot;
	double xx = xf0 + dist; 
	double yy = yf0; 
	double zz = zf0; 

	double valx1 = interpolate(u.x, xx, yy, zz);
	double valy1 = interpolate(u.y, xx, yy, zz);
	double valz1 = interpolate(u.z, xx, yy, zz);
         
        vels1[n] = sqrt(sq(valx1) + sq(valy1) + sq(valz1));
	fprintf(fpstr, "%g,%g,%g,%g,%g\n", xx, vels1[n], valx1, valy1, valz1);     
    }  
    fclose(fpstr); 
    }
}

event dts_meas(t += 1) {
    if(pid() == 0) {
	char nameDtshor[90];
    	snprintf(nameDtshor, 90, "%shor_t=%05g", out.dir_dts, t);
        FILE * fpstrhor = fopen(nameDtshor, "w");
        fprintf(fpstrhor, "x,y,z,b\n");
        
	double lengthhor = L0;
        int ntothor = L0;
	
	double xf0 = rot.x0;
	double yf0 = rot.y0;
	double zf0 = rot.z0;

	for(int n = 0; n <= ntothor; n++) {
	    double dist = lengthhor*n/ntothor;

	    double xx = xf0; 
	    double yy = 1.; 
	    double zz = dist; 
       
	    double valb = interpolate(b, xx, yy, zz);
          
	    fprintf(fpstrhor, "%g,%g,%g,%g\n", xx, yy, zz, valb);     
	}  
 	fclose(fpstrhor); 

	char nameDtsver[90];
    	snprintf(nameDtsver, 90, "%sver_t=%05g", out.dir_dts, t);
        FILE * fpstrver = fopen(nameDtsver, "w");
        fprintf(fpstrver, "x,y,z,b\n");
        
	double lengthver = 20;
        int ntotver = 200;

	for(int n = 0; n <= ntotver; n++) {
	    double dist = lengthver*n/ntotver;

	    double xx = xf0 + 30; 
	    double yy = dist; 
	    double zz = zf0; 
       
	    double valb = interpolate(b, xx, yy, zz);
          
	    fprintf(fpstrver, "%g,%g,%g,%g\n", xx, yy, zz, valb);     
	}  
 	fclose(fpstrver); 
    }
    
}


#endif

/** Ouputting movies in 2- or 3-D*/

#if dimension == 2
event movies(t += 0.5) {
    vertex scalar omega[]; 	// Vorticity
    scalar lev[];	 	// Grid depth
    scalar ekinRho[]; 		// Kinetic energy

    foreach() {
        omega[] = ((u.y[1,0] - u.y[-1,0]) - (u.x[0,1] - u.x[0,-1]))/(2*Delta); // Curl(u) 
        ekinRho[] = 0.5*rho[]*(sq(u.x[]) + sq(u.y[]));
        lev[] = level;
    }

    boundary ({b, lev, omega, ekinRho});
    output_ppm (b, file = "ppm2mp4 ./results/buoyancy.mp4", n = 1<<maxlevel, linear = true, max=STRAT(L0), min=STRAT(0.));
    output_ppm (ekinRho, file = "ppm2mp4 ./results/ekin.mp4", n = 1<<maxlevel, min = 0, max = 1.*sq(rot.cu));
    output_ppm (omega, file = "ppm2mp4 ./results/vort.mp4", n = 1<<maxlevel, linear = true); 
    output_ppm (lev, file = "ppm2mp4 ./results/grid_depth.mp4", n = 1<<maxlevel, min = minlevel, max = maxlevel);
}
#elif dimension == 3
event movies(t += out.dtVisual) {
    scalar l2[];

    lambda2(u,l2);
    boundary({l2});
    
    view(fov=25, tx = 0., ty = 0., phi=bvsets.phi, theta=bvsets.theta, width = 800, height = 800);
    
    translate(-rot.x0,-rot.y0,-rot.z0) {
        box(notics=false);
        isosurface("l2", v=-0.02, color="b", min=STRAT(roughY0h), max=STRAT(2.*rot.y0));
	draw_vof("fan", fc = {1,0,0});
    }
    translate(-rot.z0,-rot.y0, -L0){
      	squares("u.x", n = {0,0,1}, alpha=rot.z0, min=-fabs(WIND(1.5*rot.y0)), max=fabs(WIND(1.5*rot.y0)));
        cells(n = {0,0,1}, alpha = rot.z0);
    }
    
    translate(0.,-rot.y0,-rot.z0){
        squares("u.x", n = {1,0,0}, alpha=rot.x0, min=-fabs(WIND(1.5*rot.y0)), max=fabs(WIND(1.5*rot.y0)));
    }

    /** Save file with certain fps*/
    char nameVid1[90];
    snprintf(nameVid1, 90, "ppm2mp4 -r %g ./%s/visual_3d_vel.mp4", 10., out.dir);
    save(nameVid1);
    clear();

    view(fov=25, tx = 0., ty = 0., phi=bvsets.phi, theta=bvsets.theta, width = 800, height = 800);
    
    translate(-rot.x0,-rot.y0,-rot.z0) {
        box(notics=false);
        isosurface("l2", v=-0.02, color="b", min=STRAT(0.), max=STRAT(1.5*rot.y0));
	draw_vof("fan", fc = {1,0,0});
    }
    translate(-rot.z0,-rot.y0, -L0){
      	squares("b", n = {0,0,1}, alpha=rot.z0, min=STRAT(0.), max=STRAT(1.5*rot.y0));
        cells(n = {0,0,1}, alpha = rot.z0);
    }
    
    translate(0.,-rot.y0,-rot.z0){
        squares("b", n = {1,0,0}, alpha=rot.x0, min=STRAT(0.), max=STRAT(1.5*rot.y0));
    }

    /** Save file with certain fps*/
    char nameVid2[90];
    snprintf(nameVid2, 90, "ppm2mp4 -r %g ./%s/visual_3d_b.mp4", 10., out.dir);
    save(nameVid2);
    clear();
}

/** Take relevant field slices and write away */
event slices(t=out.dtSlices; t+=out.dtSlices) {
    char nameSlice[90];
    coord slice = {1., 0., 1.};
    int res = L0/2;

    for(double yTemp = 0.5; yTemp<=1; yTemp+=0.5) {
	slice.y = yTemp/L0;
   
    	snprintf(nameSlice, 90, "%st=%05gy=%03g", out.dir_slices, t, yTemp);
    	FILE * fpsli = fopen(nameSlice, "w");
    	output_slice(list = (scalar *){b}, fp = fpsli, n = res, linear = true, plane=slice);
    	fclose(fpsli);
    }

    for(double yTemp = 2; yTemp<=4.; yTemp+=2.) {
	slice.y = yTemp/L0;
   
    	snprintf(nameSlice, 90, "%st=%05gy=%03g", out.dir_slices, t, yTemp);
    	FILE * fpsli = fopen(nameSlice, "w");
    	output_slice(list = (scalar *){b}, fp = fpsli, n = res, linear = true, plane=slice);
    	fclose(fpsli);
    }

/*
    slice.y = 1.;
    slice.z = 1.;

    for(double xTemp = rot.x0-40; xTemp<=rot.x0+40.; xTemp+=20.) {
	slice.x = xTemp/L0;
	
	snprintf(nameSlice, 90, "%st=%05gx=%03g", out.dir_slices, t, xTemp);
	FILE * fpsli = fopen(nameSlice, "w");
	output_slice(list = (scalar *){b}, fp = fpsli, n = 256, linear = true, plane=slice);
	fclose(fpsli);
	}

    slice.x = 1.;
    slice.y = 1.;
    slice.z = rot.z0/L0;
    	
    snprintf(nameSlice, 90, "%st=%05gz=%03g", out.dir_slices, t, rot.z0);
    FILE * fpsli = fopen(nameSlice, "w");
    output_slice(list = (scalar *){b}, fp = fpsli, n = 256, linear = true, plane=slice);
    fclose(fpsli);
*/
}
#endif

/** Usefull functions */ 

/** Checks if required folders exists, if not they get created. */
void sim_dir_create(){
   
    sprintf(out.dir, "./%s/%s%02d", out.main_dir, sim_ID, out.sim_i);
    sprintf(out.dir_profiles, "%s/profiles/", out.dir);
    sprintf(out.dir_slices, "%s/slices/", out.dir);
    sprintf(out.dir_equifields, "%s/equifields/", out.dir);
    sprintf(out.dir_strvel, "%s/strvel/", out.dir);
    sprintf(out.dir_refvel, "%s/refvel/", out.dir);
    sprintf(out.dir_diffbins, "%s/diffbins/", out.dir);
    sprintf(out.dir_dts, "%s/dts/", out.dir);

    
    if (pid() == 0){
    struct stat st = {0};
    if (stat(out.main_dir, &st) == -1) {
        mkdir(out.main_dir, 0777);
    }
   
    if (stat(out.dir, &st) == -1) {
        mkdir(out.dir, 0777);
    }
    if (stat(out.dir_slices, &st) == -1) {
        mkdir(out.dir_slices, 0777);
    }
    if (stat(out.dir_profiles, &st) == -1) {
        mkdir(out.dir_profiles, 0777);
    }  
    if (stat(out.dir_equifields, &st) == -1) {
        mkdir(out.dir_equifields, 0777);
    }  
    if (stat(out.dir_strvel, &st) == -1) {
	mkdir(out.dir_strvel, 0777);
    }
    if (stat(out.dir_refvel, &st) == -1) {
	mkdir(out.dir_refvel, 0777);
    }
    if (stat(out.dir_diffbins, &st) == -1) {
	mkdir(out.dir_diffbins, 0777);
    }
    if (stat(out.dir_dts, &st) == -1) {
	mkdir(out.dir_dts, 0777);
    }

    }
}


