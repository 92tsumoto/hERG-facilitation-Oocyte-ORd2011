/* produced by Tsumoto. K 2008.10.27 */

#include <string.h>
#include <stdlib.h>
#include "syspara.h"

FILE *fopen(), *fpin, *fp0, *fp1, *fp2, *fp3;
FILE *fp4,*fp5,*fp6,*fp7,*fp8,*fp9,*fp10,*fp11,*fp12,*fp13;
int mode = 1;
int P = 2;
int beats = 30;

typedef double Number;

main(argc,argv)
int argc;
char **argv;
{
	int i,w;
	int ii=0;
	double x[NN];
	double t = 0.0;
	double tt,ttt;
	double time=0.0;
	double ttime=0.0;
	double h;
	double v_old,dvdt,dvdt_new;
	double t_stok;
	char *tmpname;
	char cmd[BUFSIZ];
	double tend;
	double apd,rmbp,vmax,toneapd,ttwoapd,dvdtmax;

/* Action Potential Duration and Max. Info */
/*	double *vmax ; // Max. Voltage (mV)
	double *dvdtmax ; // Max. dv/dt (mV/ms)
	double *apd; // Action Potential Duration
	double *toneapd; // Time of dv/dt Max.
	double *ttwoapd; // Time of 90% Repolarization
	double *rmbp; // Resting Membrane Potential
	double *nair; // Intracellular Na At Rest
	double *cair; // Intracellular Ca At Rest
	double *kir ; // Intracellular K At Rest
	double caimax [beats] ; // Peak Intracellular Ca

	vmax=(Number *)calloc(beats,sizeof(Number));
	dvdtmax=(Number *)calloc(beats,sizeof(Number));
	apd=(Number *)calloc(beats,sizeof(Number));
	toneapd=(Number *)calloc(beats,sizeof(Number));
	ttwoapd=(Number *)calloc(beats,sizeof(Number));
	rmbp=(Number *)calloc(beats,sizeof(Number));
	nair=(Number *)calloc(beats,sizeof(Number));
	cair=(Number *)calloc(beats,sizeof(Number));
	kir=(Number *)calloc(beats,sizeof(Number));
	if(vmax==NULL || dvdtmax==NULL || apd==NULL || toneapd==NULL || ttwoapd==NULL 
		|| rmbp==NULL || nair==NULL || cair==NULL || kir==NULL
		) exit(1);
*/
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
	if ((fp3 = fopen("nstate.out","w")) == NULL){
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

	if(var.out_data){
		if ((fp4 = fopen("ikr_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp5 = fopen("iks_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp6 = fopen("ik1_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp7 = fopen("ito_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp8 = fopen("inak_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp9 = fopen("incx_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp10 = fopen("ina_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp11 = fopen("ical_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp12 = fopen("itot_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp13 = fopen("irel_data.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
	}
	xhplot(WINDOW, 700.0, 700.0, WHITE);
	xhplot(DIRECT, 0.0, 0.0, WHITE);

	for (ii = 0; ii < var.datas; ii++){
		long j;
		time = 0.0;
		tend = var.tend[ii];
		for (i = 0; i < NN; i++){ 
			x[i] = var.x0[ii][i];
		}

		//tt = var.ndis*(double)var.m;
		h = 1.0 / (double)var.m;
		//h = 2.0*M_PI / (double)var.m;
		h *= var.tsign[ii];

		xddp.line_wid = var.line_wid[ii];
		xhplot(LINEATT,0,0,WHITE);

		// initial values input.
		val_consts(fp1);
		//var.pswitch=1;
		printf("exit consts\n");

		//var.block_rate = 1.0/(1.0+pow((var.drug_concentration/var.ic50),var.hillc));
		var.drug_concentration = var.ic50*pow(1.0/var.block_rate-1.0,1/var.hillc);
		// affected fraction (with facilitation)
		if(var.model_type==0){ // facilitation
			var.fraction_facilitation = 1.0/(1.0+pow((var.drug_concentration/var.fic50),var.fhillc));
			printf("with facilitation model\n");
		} else {
			var.fraction_facilitation = 1.0;
			printf("without facilitation model\n");
		}
		printf("fraction=%lf,noaffected=%lf,affected=%lf,total=%lf\n",
			var.fraction_facilitation,var.block_rate*var.fraction_facilitation,var.block_rate*(1.0-var.fraction_facilitation),var.block_rate);
		fprintf(fp1,"fraction=%lf,noaffected=%lf,affected=%lf,total=%lf\n",
			var.fraction_facilitation,var.block_rate*var.fraction_facilitation,var.block_rate*(1.0-var.fraction_facilitation),var.block_rate);
	
		// initial values input.
		initial_mem();
		printf("exit memory initialization\n");

		printf("Istim=%lf\n",var.Istim_base);

		// Tablize exp functions.	
		printf("start tablization\n");
		make_ExpTable();
		printf("finished tablization\n");

		// Initialization time
		var.beat = 0;

		tt = var.ndis*(double)var.m*var.BCL;
		ttt = (1.0-var.ndis)*(double)var.m*var.BCL;
		printf("tt=%lf,ttt=%lf\n",tt,ttt);

		while (1){
			eventloop(fp1,&mode,&P,x);
			time=var.beat*var.BCL;
			
			apd = 0.0; toneapd = 0.0; ttwoapd = 0.0; rmbp = x[0]; vmax = -90.0; dvdtmax = 0.0, v_old = x[0];
			for (j = 0; j< (int)tt; j++){
				t = h*(double)j;
				v_old = x[0];
				var.Istim = var.Istim_base;
				runge(NN,h,x,t);
				dvdt_new = (x[0]-v_old)/(h);
				if(x[0] > vmax){vmax = x[0];}
				if(dvdt_new > dvdtmax){
					dvdtmax = dvdt_new;
					toneapd = time;
				}
				if(dvdt_new < 0 && x[0] >= (vmax - 0.9*(vmax - rmbp))) ttwoapd = time;

				if (var.pflag) orbit(&mode,x,dvdt_new);
				if (var.pswitch==1){
					data_out(fp2,ttime,x);
					if(var.out_data){
						current(fp4,fp5,fp6,fp7,fp8,fp9,fp10,fp11,fp12,fp13,ttime,x);
					}
				}
				time += h;
				ttime=time;
			} // 1st-loop

			for (j = 0; j< (int)ttt; j++){
				t = h*(double)j;
				v_old = x[0];
				var.Istim = 0.0;
				runge(NN,h,x,t);
				dvdt_new = (x[0]-v_old)/(h); // LRd -> dvdtnew
				if(x[0] > vmax){vmax = x[0];}
				if(dvdt_new > dvdtmax){
					dvdtmax = dvdt_new;
					toneapd = time;
				}
				if(dvdt_new < 0 && x[0] >= (vmax - 0.9*(vmax - rmbp))) ttwoapd = time;

				if (var.pflag) orbit(&mode,x,dvdt_new);
				if (var.pswitch==1){
					data_out(fp2,ttime,x);
					if(var.out_data){
						current(fp4,fp5,fp6,fp7,fp8,fp9,fp10,fp11,fp12,fp13,ttime,x);
					}
				}
				time += h;
				ttime=time;
			} // 2nd-loop

			fprintf(fp3,"#beats=%d\n",var.beat);
			for(w=0;w<NN;w++){
				fprintf(fp3,"%16.15e\n",x[w]);
			}
		
			if(var.sswitch==1){	
				printf("%d %lf ",var.beat,time);
				for(w=0;w<NN+1;w++){
					if(w!=NN){
						printf("%10.9lf ",x[w]);
					} else {
						printf("%10.9lf %10.9e %10.9lf\n",ttwoapd - toneapd,var.drug_concentration,var.block_rate);
					}
				}
			}
			draw_p(&mode,P,x,dvdt);
			mouse(&mode,x,dvdt);
			if (fabs(time) > tend &&  tend != 0.0) break;
			var.beat++;

		} // end for while loop

	} // end for ii-loop


	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	if(var.write0){
		fclose(fp4);fclose(fp5);fclose(fp6);fclose(fp7);fclose(fp8);
		fclose(fp9);fclose(fp10);fclose(fp11);fclose(fp12);fclose(fp13);
	}
	//free(vmax);free(dvdtmax);free(apd);free(toneapd);free(ttwoapd);
	//free(rmbp);free(nair);free(cair);free(kir);

}

