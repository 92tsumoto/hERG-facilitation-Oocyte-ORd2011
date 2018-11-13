#include <string.h>
#include "syspara.h"

void input_para(FILE *fpin)
{
	int i,ii;

	fscanf(fpin,"%d",&var.fail);
	fscanf(fpin,"%d",&var.celltype);
	fscanf(fpin,"%lf",&var.BCL);
	fscanf(fpin,"%lf",&var.Istim_base);
	fscanf(fpin,"%lf",&var.drug_concentration);
	fscanf(fpin,"%lf",&var.ic50);
	fscanf(fpin,"%lf",&var.hillc);
	fscanf(fpin,"%lf",&var.ko);
	fscanf(fpin,"%d",&var.datas);
	for (ii = 0; ii < var.datas; ii++){
		fscanf(fpin, "%d", &var.line_wid[ii]);
		for (i=0;i<NN;i++){
			fscanf(fpin,"%lf",&var.x0[ii][i]);
			printf("x[%d]=%e\n",i,var.x0[ii][i]);
		}
		fscanf(fpin, "%lf", &var.tsign[ii]);
		fscanf(fpin, "%lf", &var.tend[ii]);
	}
	fscanf(fpin,"%lf",&var.dIstim);
	fscanf(fpin,"%d",&var.l);
	fscanf(fpin,"%d",&var.m);
	fscanf(fpin,"%d",&var.pflag);
	fscanf(fpin,"%d",&var.write);
	fscanf(fpin,"%d",&var.write0);
	fscanf(fpin,"%d",&var.half);

}

