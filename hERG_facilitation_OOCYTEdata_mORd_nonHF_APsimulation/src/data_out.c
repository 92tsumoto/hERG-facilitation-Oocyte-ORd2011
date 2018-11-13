#include "syspara.h"

void data_out(FILE *fp2, double t, double u[])
{
	int i;

	fprintf(fp2,"%lf %lf ",t,t/var.omega);
	for(i=0;i<NN;i++){
		fprintf(fp2,"%10.9lf ",u[i]);
	}
	fprintf(fp2,"\n");
}

void current(FILE *fp4, FILE *fp5, FILE *fp6, FILE *fp7, FILE *fp8, FILE *fp9, FILE *fp10, FILE *fp11, FILE *fp12, FILE *fp13, double t, double u[])
{

	current_ikr(fp4,t,u);
	current_iks(fp5,t,u);
	current_ik1(fp6,t,u);
	current_ito(fp7,t,u);
	current_inak(fp8,t,u);
	current_incx(fp9,t,u);
	current_ina(fp10,t,u);
	current_ical(fp11,t,u);
	current_it(fp12,t,u);
	current_irel(fp13,t,u);

	//printf("t=%lf\n",t);

}

// Ikr, IKr2 
void current_ikr (FILE *fp4, double time, double p[])
{
	fprintf(fp4,"%lf %lf %lf\n",time,ikr.ik,ikr2.ik);

}

// IKs 
void current_iks (FILE *fp5, double time, double p[])
{
	
	fprintf(fp5,"%lf %lf\n",time,iks.ik);

}

// IK1 
void current_ik1 (FILE *fp6, double time, double p[])
{
	
	fprintf(fp6,"%lf %lf\n",time,ik1.ik);

}

// Ito 
void current_ito (FILE *fp7, double time, double p[])
{
	
	fprintf(fp7,"%lf %lf\n",time,ito.ik);

}

// INak

void current_inak (FILE *fp8, double time, double p[])
{
	fprintf(fp8,"%lf %lf\n",time,inak.inak);

}
// Incx

void current_incx (FILE *fp9, double time, double p[])
{
	fprintf(fp9,"%lf %lf %lf %lf\n",time,var.inaca_i,var.inaca_ss,var.inaca);

}

// INa

void current_ina (FILE *fp10, double time, double p[])
{
	fprintf(fp10,"%lf %lf %lf %lf\n",time,ina.fast,ina.late,ina.total);

}

// L-type calcium current
void current_ical(FILE *fp11, double time, double p[])
{

	fprintf(fp11,"%lf %lf %lf %lf\n",time,ical.ica,ical.icana,ical.icak);

}

// Itotal 
void current_it (FILE *fp12, double time, double p[])
{
	
	fprintf(fp12,"%10.9lf %10.9lf\n",time,var.Itotal);

}

// Jrel
void current_irel (FILE *fp13, double time, double p[])
{
	
	fprintf(fp13,"%10.9lf %10.9lf\n",time,jrel.ca);

}

