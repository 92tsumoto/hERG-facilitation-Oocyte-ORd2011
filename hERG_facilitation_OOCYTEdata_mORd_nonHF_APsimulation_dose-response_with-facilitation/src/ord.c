/* produced by Tsumoto. K 2008.10.27 */

#include <string.h>
#include <stdlib.h>
#include "syspara.h"

FILE *fopen(), *fpin, *fp0, *fp1, *fp2, *fp3, *fp4;
int mode = 1;
int P = 2;
int beats = 300;

typedef double Number;

main(argc,argv)
int argc;
char **argv;
{
	int i,w;
	int ii=0;
	double x[NN];
	double t = 0.0;
	double time=0.0;
	double h;
	double v_old,dvdt,dvdt_new;
	double t_stok;
	char *tmpname;
	char cmd[BUFSIZ];
	double tend;
	double ave_apd90,ave_apd50,ave_apd20;

/* Action Potential Duration and Max. Info */
	double *vmax ; // Max. Voltage (mV)
	double *dvdtmax ; // Max. dv/dt (mV/ms)
	double *apd90,*apd50,*apd20; // Action Potential Duration
	double *toneapd; // Time of dv/dt Max.
	double *ttwoapd90,*ttwoapd50,*ttwoapd20; // Time of 90% Repolarization
	double *rmbp; // Resting Membrane Potential
	double *nair; // Intracellular Na At Rest
	double *cair; // Intracellular Ca At Rest
	double *kir ; // Intracellular K At Rest
	double caimax [beats] ; // Peak Intracellular Ca

	vmax=(Number *)calloc(beats,sizeof(Number));
	dvdtmax=(Number *)calloc(beats,sizeof(Number));
	apd90=(Number *)calloc(beats,sizeof(Number));
	apd50=(Number *)calloc(beats,sizeof(Number));
	apd20=(Number *)calloc(beats,sizeof(Number));
	toneapd=(Number *)calloc(beats,sizeof(Number));
	ttwoapd90=(Number *)calloc(beats,sizeof(Number));
	ttwoapd50=(Number *)calloc(beats,sizeof(Number));
	ttwoapd20=(Number *)calloc(beats,sizeof(Number));
	rmbp=(Number *)calloc(beats,sizeof(Number));
	nair=(Number *)calloc(beats,sizeof(Number));
	cair=(Number *)calloc(beats,sizeof(Number));
	kir=(Number *)calloc(beats,sizeof(Number));
	if(vmax==NULL || dvdtmax==NULL || apd90==NULL || apd50==NULL || apd20==NULL
		|| toneapd==NULL || ttwoapd90==NULL || ttwoapd50==NULL || ttwoapd20==NULL
		|| rmbp==NULL || nair==NULL || cair==NULL || kir==NULL
		) exit(1);

//int i; // Stimulation Counter

	tmpname = "temp";

	sprintf(cmd, "/usr/bin/cpp -P %s > %s", argv[1],tmpname);
	if(system(cmd) == -1){
		fprintf(stderr,"cannot open %s\n",argv[1]);
		exit(1);
	}
	if((fpin=fopen(tmpname,"r"))==NULL){
		fprintf(stderr,"cannot open %s\n",argv[1]);
		exit(1);
	}
	if ((fp1 = fopen("para.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp2 = fopen("data.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp3 = fopen("ndata.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp4 = fopen("initdata.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}

// parameter inputs
	input_para(fpin);

	if (var.write){
		if ((fp0 = fopen(argv[2],"w"))==NULL){
			fprintf(stderr, "%s cannot open.\n",argv[2]);
			exit(-1);
		}
	}

	for (ii = 0; ii < var.datas; ii++){
		long j;
		time = 0.0;
		tend = var.tend[ii];
		
		h = 1.0 / var.m;
		h *= var.tsign[ii];
		
		// initial values input.
		val_consts(fp1);
		printf("exit consts\n");
	
	// initial values input.
		initial_mem();
		printf("exit memory initialization\n");

		printf("Istim=%lf\n",var.Istim_base);

	// Tablize exp functions.	
		printf("start tablization\n");
		make_ExpTable();
		printf("finished tablization\n");

//for(var.drug_concentration=0;var.drug_concentration<=300E-6;var.drug_concentration+=1E-6){
for(var.block_rate=1;var.block_rate>0;var.block_rate-=0.001){

	// Initialization time
		time = 0.0;
		time -= h;
		var.dt = h;
		var.beat = 0;
		ave_apd90 = 0.0; ave_apd50 = 0.0; ave_apd20 = 0.0;

		for (i = 0; i < NN; i++){ 
			x[i] = var.x0[ii][i];
		}

	// unaffected fraction (without facilitation)
	//var.block_rate = 1/(1+pow((var.drug_concentration/var.ic50),var.hillc));
	var.drug_concentration = var.ic50*pow(1.0/var.block_rate-1.0,1/var.hillc);

	// affected fraction (with facilitation)
	var.fraction_facilitation = 1/(1+pow((var.drug_concentration/var.fic50),var.fhillc));
	printf("noaffected=%lf,affected=%lf,total=%lf\n",var.block_rate*var.fraction_facilitation,var.block_rate*(1-var.fraction_facilitation),var.block_rate);

		for (var.beat=0; var.beat < beats; var.beat++){
			if(var.beat == beats-1){
				var.l = 3.0*var.BCL;
			} else {
				var.l = var.BCL;
			}

			for (j = 0; j< (var.m * var.l ); j++){
				t = h*j;
				time += h;

				if ( time-(var.BCL*var.beat+10.0) >= 0.0 && time-(var.BCL*var.beat+10.0) < h ){
					apd90[var.beat] =0; apd50[var.beat] =0; apd20[var.beat] =0;
					toneapd[var.beat] =0; ttwoapd90[var.beat] =0; ttwoapd50[var.beat] =0; ttwoapd20[var.beat] =0;
					rmbp[var.beat] =x[0];
					nair[var.beat] = x[33];	kir[var.beat] = x[35];	cair[var.beat] = x[37];	caimax[var.beat] = x[37];
					vmax[var.beat] = -90.0;
					dvdtmax[var.beat] = 0.0;
				}

				if (time-(var.BCL*var.beat+10.0) >= 0.0 && time-(var.BCL*var.beat+10.0) < 0.5){
					var.Istim = var.Istim_base;
				} else {
					var.Istim = 0;
				}

				if (fabs(time) > tend &&  tend != 0.0) break;

				v_old = x[0];

				eular(NN,h,x,t);
				
				dvdt_new = (x[0]-v_old)/h;

				if(var.beat>=0){
					if (x[0] > vmax[var.beat] )
						vmax[var.beat] = x[0];
					if (x[37] > caimax[var.beat] )
						caimax[var.beat] = x[37];
					if (dvdt_new > dvdtmax[var.beat] ){
						dvdtmax[var.beat] = dvdt_new;
						toneapd[var.beat] = time;
					}
					if (dvdt_new < 0 && x[0] >= (vmax[var.beat] -0.9*(vmax[var.beat]-rmbp[var.beat]) ) )
						ttwoapd90[var.beat] = time;
					if (dvdt_new < 0 && x[0] >= (vmax[var.beat] -0.5*(vmax[var.beat]-rmbp[var.beat]) ) )
						ttwoapd50[var.beat] = time;
					if (dvdt_new < 0 && x[0] >= (vmax[var.beat] -0.2*(vmax[var.beat]-rmbp[var.beat]) ) )
						ttwoapd20[var.beat] = time;
				}


				//if (time>= (beats-5)*var.BCL && time < beats*var.BCL){
				//	fprintf(fp2,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
				//		time,x[0],ikr.ik,ikr2.ik,ikr.xrss,ikr2.xrss,ina.mss,ina.hss,ina.jss,x[1],ina.h,x[4],ikr.xr,ikr2.xr);
				//}
				
				dvdt = dvdt_new;

			} // end j-loop

			//fprintf(fp4,"#beats=%d\n",var.beat-1);
			//for(w=0;w<NN;w++){
			//	fprintf(fp4,"%e\n",x[w]);
			//}

			if (fabs(time) > tend &&  tend != 0.0) break;

		} // end for var.beat-loop

		// Data output
		for(i=0;i<5;i++){
			ave_apd90 += ttwoapd90 [beats-(i+1)] -toneapd [beats-(i+1)] ;
			ave_apd50 += ttwoapd50 [beats-(i+1)] -toneapd [beats-(i+1)] ;
			ave_apd20 += ttwoapd20 [beats-(i+1)] -toneapd [beats-(i+1)] ;
		}
		apd90[beats-1] = ave_apd90/5.0 ;
		apd50[beats-1] = ave_apd50/5,0;
		apd20[beats-1] = ave_apd20/5.0;
			fprintf(fp3,"%e\t%e\t%g\t%g\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",var.drug_concentration,var.block_rate,
				vmax[beats-1],dvdtmax[beats-1],apd90[beats-1],apd50[beats-1],apd20[beats-1],toneapd[beats-1],
				ttwoapd90[beats-1],ttwoapd50[beats-1],ttwoapd20[beats-1],cair[beats-1],caimax[beats-1],rmbp[beats-1]);
			printf("%e\t%e\t%g\t%g\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",var.drug_concentration,var.block_rate,
				vmax[beats-1],dvdtmax[beats-1],apd90[beats-1],apd50[beats-1],apd20[beats-1],toneapd[beats-1],
				ttwoapd90[beats-1],ttwoapd50[beats-1],ttwoapd20[beats-1],cair[beats-1],caimax[beats-1],rmbp[beats-1]);

	} // end for drug_concentration
	
	} // end for ii-loop

	/* -------------------------------------------------------------------- */
	/* .. Termination and release of memory. */
	/* -------------------------------------------------------------------- */
		fclose(fp1);
		fclose(fp2);
		fclose(fp3);
		fclose(fp4);
		free(vmax);free(dvdtmax);free(apd90);free(apd50);free(apd20);
		free(toneapd);free(ttwoapd90);free(ttwoapd50);free(ttwoapd20);
		free(rmbp);free(nair);free(cair);free(kir);
}

