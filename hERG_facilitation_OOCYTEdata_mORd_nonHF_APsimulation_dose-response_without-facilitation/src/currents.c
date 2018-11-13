#include "syspara.h"

void comp_ina(double x[])
{
	//MKL_INT iV=0;
	int iV=0;
	double V1,V2,d1,d2;
	
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;
	//printf("iV=%d,V1=%f,V2=%f,d1=%f,d2=%f\n",iV,V1,V2,d1,d2);

	ina.mss = ina.Tmss[iV]*d2 + ina.Tmss[iV+1]*d1;
	ina.taum = ina.Ttaum[iV]*d2 + ina.Ttaum[iV+1]*d1;
	ina.hss = ina.Thss[iV]*d2 + ina.Thss[iV+1]*d1;
	ina.tauh_fast = ina.Ttauh_fast[iV]*d2 + ina.Ttauh_fast[iV+1]*d1;
	ina.tauh_slow = ina.Ttauh_slow[iV]*d2 + ina.Ttauh_slow[iV+1]*d1;
	ina.h = ina.Ah_fast*x[2] + ina.Ah_slow*x[3];
	ina.jss = ina.hss;
	ina.tauj = ina.Ttauj[iV]*d2 + ina.Ttauj[iV+1]*d1;

	//ina.jss = ina.hss;
	//ina.tauj = 2.038+1.0/(0.02136*exp(-(x[0]+100.6)/8.281)+0.3052*exp((x[0]+0.9941)/38.45));
	
	ina.hCaMKss = ina.ThCaMKss[iV]*d2 + ina.ThCaMKss[iV+1]*d1;
	ina.tauh_CaMK_slow = 3.0*ina.tauh_slow;
	//ina.hCaMK = ina.Ah_CaMK_fast*x[2] + ina.Ah_CaMK_slow*x[5];	// h_CaMK_fast = h_fast (x[2])
	ina.hCaMK = ina.Ah_fast*x[2] + ina.Ah_slow*x[5];	// h_CaMK_fast = h_fast (x[2])

	//ina.jCaMKss = ina.jss; //jss=hss;
	ina.tauj_CaMK = 1.46*ina.tauj;
	
	ina.fast_pCaMK = 1.0/(1.0 + CaMK.Km/CaMK.active);
	ina.fast = ina.Gna_fast*(x[0]-var.Ena)*x[1]*x[1]*x[1]*((1.0-ina.fast_pCaMK)*ina.h*x[4]+ina.fast_pCaMK*ina.hCaMK*x[6]);

	ina.mlss = ina.Tmlss[iV]*d2 + ina.Tmlss[iV+1]*d1;
	ina.tauml = ina.taum;
	ina.hlss = ina.Thlss[iV]*d2 + ina.Thlss[iV+1]*d1;
	ina.hlCaMKss = ina.ThlCaMKss[iV]*d2 + ina.ThlCaMKss[iV+1]*d1;

	ina.late_pCaMK = 1.0/(1.0 + CaMK.Km/CaMK.active);
	ina.late = ina.Gna_late*(x[0]-var.Ena)*x[7]*((1.0-ina.late_pCaMK)*x[8] + ina.late_pCaMK*x[9]);
	
	ina.total = ina.fast + ina.late;
}

// Ito Transient Outward Current
void comp_ito (double x[])
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ito.ass = ito.Tass[iV]*d2 + ito.Tass[iV+1]*d1;
	ito.taua = ito.Ttaua[iV]*d2 + ito.Ttaua[iV+1]*d1;

	ito.iss = ito.Tiss[iV]*d2 + ito.Tiss[iV+1]*d1;
	if(var.celltype==1){
		ito.depi = ito.Tdepi[iV]*d2 + ito.Tdepi[iV+1]*d1;
	} else {
		ito.depi = 1.0;
	}
	ito.taui_fast = ito.depi*(ito.Ttaui_fast[iV]*d2 + ito.Ttaui_fast[iV+1]*d1);
	ito.taui_slow = ito.depi*(ito.Ttaui_slow[iV]*d2 + ito.Ttaui_slow[iV+1]*d1);
	ito.Ai_fast = ito.TAi_fast[iV]*d2 + ito.TAi_fast[iV+1]*d1;

	if(var.fail == 1){
		ito.Ai_fast *= 0.64;
	}

	ito.Ai_slow = 1.0 - ito.Ai_fast;
	ito.i = ito.Ai_fast*x[11] + ito.Ai_slow*x[12];

	ito.aCaMKss = ito.TaCaMKss[iV]*d2 + ito.TaCaMKss[iV+1]*d1;
	//ito.taua_CaMK = ito.taua;

	//ito.iCaMKss = ito.iss;
	ito.deltaCaMK_dev = ito.TdeltaCaMK_dev[iV]*d2 + ito.TdeltaCaMK_dev[iV+1]*d1;
	ito.deltaCaMK_rec = ito.TdeltaCaMK_rec[iV]*d2 + ito.TdeltaCaMK_rec[iV+1]*d1;
	ito.taui_CaMK_fast = ito.taui_fast*ito.deltaCaMK_dev*ito.deltaCaMK_rec;
	ito.taui_CaMK_slow = ito.taui_slow*ito.deltaCaMK_dev*ito.deltaCaMK_rec;
	ito.Ai_CaMK_fast = ito.Ai_fast;
	ito.Ai_CaMK_slow = ito.Ai_slow;
	if(var.fail == 1){
		ito.Ai_CaMK_slow *= 1.5;
		ito.taui_CaMK_slow *= 0.082;
	}
	ito.iCaMK = ito.Ai_CaMK_fast*x[14] + ito.Ai_CaMK_slow*x[15];
	//ito.iCaMK = ito.Ai_fast*x[14] + ito.Ai_slow*x[15];

	ito.pCaMK = 1.0/(1.0 + CaMK.Km/CaMK.active);
	ito.ik = ito.Gto*(x[0]-var.Ek)*((1.0-ito.pCaMK)*x[10]*ito.i + ito.pCaMK*x[13]*ito.iCaMK);

}


// L-type calcium current
void comp_ical(double x[])
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

// VDA
	ical.dss = ical.Tdss[iV]*d2 + ical.Tdss[iV+1]*d1;
	ical.taud = ical.Ttaud[iV]*d2 + ical.Ttaud[iV+1]*d1;
// VDI 
	ical.fss = ical.Tfss[iV]*d2 + ical.Tfss[iV+1]*d1;
	ical.tauf_fast = ical.Ttauf_fast[iV]*d2 + ical.Ttauf_fast[iV+1]*d1;
	ical.tauf_slow = ical.Ttauf_slow[iV]*d2 + ical.Ttauf_slow[iV+1]*d1;
	ical.f = ical.Af_fast*x[17] + ical.Af_slow*x[18];

// CDI 
	//ical.fcass = ical.fss;
	ical.taufca_fast = ical.Ttaufca_fast[iV]*d2 + ical.Ttaufca_fast[iV+1]*d1;
	ical.taufca_slow = ical.Ttaufca_slow[iV]*d2 + ical.Ttaufca_slow[iV+1]*d1;
	ical.Afca_fast = ical.TAfca_fast[iV]*d2 + ical.TAfca_fast[iV+1]*d1;
	ical.Afca_slow = 1.0 - ical.Afca_fast;
	ical.fca = ical.Afca_fast*x[19] + ical.Afca_slow*x[20];

	//ical.jcass = ical.fcass;
	//ical.taujca = 75.0; // (ms).

// CaMK(VDI)
    //ical.f_CaMKss = ical.fss;
	ical.tauf_CaMK_fast = 2.5*ical.tauf_fast;
	//ical.Af_CaMK_fast = ical.Af_fast;
	//ical.Af_CaMK_slow = ical.Af_slow;
	//ical.f_CaMK_slow = x[18]; //(x[18] = f_slow )
	ical.f_CaMK = ical.Af_fast*x[22] + ical.Af_slow*x[18];

// CaMK(CDI)
	
	//ical.fca_CaMKss = ical.fss;
	ical.taufca_CaMK_fast = 2.5*ical.taufca_fast;
	//ical.Afca_CaMK_fast = ical.Afca_fast;
	//ical.Afca_CaMK_slow = ical.Afca_slow;
	//var.fca_CaMK_slow = x[20]; //(x[20] = fca_slow )
	ical.fca_CaMK = ical.Afca_fast*x[23] + ical.Afca_slow*x[20];

	//ical.kmn = 0.002;
	//ical.kp2n = 1000.0;
	ical.km2n = 1.0*x[21]; //(jca: recovery from CDI for ICaL )
	ical.tmp = (1.0+ical.kmn/x[38])*(1.0+ical.kmn/x[38])*(1.0+ical.kmn/x[38])*(1.0+ical.kmn/x[38]);
	ical.alpha_n = 1.0/((ical.kp2n/ical.km2n)+ical.tmp);
	
	ical.exp_Ca = ical.Texp_Ca[iV]*d2 + ical.Texp_Ca[iV+1]*d1;
	ical.exp_Na = ical.Texp_Na[iV]*d2 + ical.Texp_Na[iV+1]*d1;
	ical.exp_K = ical.Texp_K[iV]*d2 + ical.Texp_K[iV+1]*d1;

	//if(fabs(x[0])>1E-8){
		ical.phi_ca = (4.0*F*F*x[0]/R/T)*(ical.gacai*x[38]*ical.exp_Ca-ical.gacao*var.cao)/(ical.exp_Ca-1.0);
		ical.phi_na = (1.0*F*F*x[0]/R/T)*(ical.ganai*x[34]*ical.exp_Na-ical.ganao*var.nao)/(ical.exp_Na-1.0);
		ical.phi_k = (1.0*F*F*x[0]/R/T)*(ical.gaki*x[36]*ical.exp_K-ical.gako*var.ko)/(ical.exp_K-1.0);
	//} else {
	//	ical.phi_ca = 0.0;//(4.0*F*F*x[0]/R/T)*(ical.gacai*x[38]*ical.exp_Ca-ical.gacao*var.cao)/(ical.exp_Ca-1.0);
	//	ical.phi_na = 0.0;//(1.0*F*F*x[0]/R/T)*(ical.ganai*x[34]*ical.exp_Na-ical.ganao*var.nao)/(ical.exp_Na-1.0);
	//	ical.phi_k = 0.0;//(1.0*F*F*x[0]/R/T)*(ical.gaki*x[36]*ical.exp_K-ical.gako*var.ko)/(ical.exp_K-1.0);
	//}
	ical.ibarcal = ical.pca*ical.phi_ca;
	ical.ibarcana = ical.pcana*ical.phi_na;
	ical.ibarcak = ical.pcak*ical.phi_k;

	ical.ibarcal_CaMK = ical.pca_CaMK*ical.phi_ca;
	ical.ibarcana_CaMK = ical.pcana_CaMK*ical.phi_na;
	ical.ibarcak_CaMK = ical.pcak_CaMK*ical.phi_k;

	ical.pCaMK = 1.0/(1.0 + CaMK.Km/CaMK.active);

	ical.ica =ical.ibarcal*x[16]*(1.0-ical.pCaMK)*(ical.f*(1.0-x[24])+ical.fca*x[24]*x[21])
				+ical.ibarcal_CaMK*x[16]*ical.pCaMK*(ical.f_CaMK*(1.0-x[24])+ical.fca_CaMK*x[24]*x[21]);
	ical.icana = ical.ibarcana*x[16]*(1.0- ical.pCaMK)*(ical.f*(1.0-x[24])+ical.fca*x[24]*x[21])
				+ ical.ibarcana_CaMK*x[16]*ical.pCaMK*(ical.f_CaMK*(1.0-x[24])+ical.fca_CaMK*x[24]*x[21]);
	ical.icak = ical.ibarcak*x[16]*(1.0- ical.pCaMK)*(ical.f*(1.0-x[24])+ical.fca*x[24]*x[21])
				+ ical.ibarcak_CaMK*x[16]*ical.pCaMK*(ical.f_CaMK*(1.0-x[24])+ical.fca_CaMK*x[24]*x[21]);

}

// Rapidly Activating Potassium Current (non-facilitation component)
void comp_ikr (double x[])
{
	MKL_INT iV=0;	
	double V1,V2,d1,d2;
	
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ikr.xrss = ikr.Txrss[iV]*d2 + ikr.Txrss[iV+1]*d1;
	ikr.tauxr_fast = ikr.Ttauxr_fast[iV]*d2 + ikr.Ttauxr_fast[iV+1]*d1;
	ikr.tauxr_slow = ikr.Ttauxr_slow[iV]*d2 + ikr.Ttauxr_slow[iV+1]*d1;
	ikr.Axr_fast = ikr.TAxr_fast[iV]*d2 + ikr.TAxr_fast[iV+1]*d1;
	ikr.Axr_slow = 1.0 - ikr.Axr_fast;
	ikr.rkr = ikr.Trkr[iV]*d2 + ikr.Trkr[iV+1]*d1;

	ikr.xr = ikr.Axr_fast*x[25] + ikr.Axr_slow*x[26];

	ikr.ik = var.block_rate*ikr.Gkr*ikr.rategkr*ikr.xr*ikr.rkr*(x[0]-var.Ek);

}

// Slowly Activating Potassium Current 
void comp_iks (double x[])
{
	
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	iks.xs1ss = iks.Txs1ss[iV]*d2 + iks.Txs1ss[iV+1]*d1;
	//iks.xs2ss = iks.xs1ss;
	iks.tauxs1 = iks.Ttauxs1[iV]*d2 + iks.Ttauxs1[iV+1]*d1;
	iks.tauxs2 = iks.Ttauxs2[iV]*d2 + iks.Ttauxs2[iV+1]*d1;
	iks.KsCa = 1.0+0.6/(1.0+pow(3.8e-5/x[37],1.4));
	iks.ik = iks.Gks*iks.KsCa*x[27]*x[28]*(x[0]-var.Eks);

}

// Inward rectifier potassium current (Ik1)
void comp_ik1 (double x[])
{
        
	MKL_INT iV=0;   
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ik1.k1ss = ik1.Tk1ss[iV]*d2 + ik1.Tk1ss[iV+1]*d1;
	ik1.tauk1 = ik1.Ttauk1[iV]*d2 + ik1.Ttauk1[iV+1]*d1;
	ik1.rk1 = ik1.Trk1[iV]*d2 + ik1.Trk1[iV+1]*d1;

	ik1.ik = ik1.Gk1*ik1.rategk1*x[29]*ik1.rk1*(x[0]-var.Ek);

}

// Sodium-Calcium Exchanger V-S

void comp_inaca (double x[])
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	var.hca=var.Thca[iV]*d2 + var.Thca[iV+1]*d1;
	var.hna=var.Thna[iV]*d2 + var.Thna[iV+1]*d1;

	// intracallular space
	ncxi.h1 = 1.0 + (x[33]/var.kna3)*(1.0 + var.hna);
	ncxi.h2 = x[33]*var.hna/(var.kna3*ncxi.h1);
	ncxi.h3 = 1.0/ncxi.h1;
	ncxi.h4 = 1.0 + (x[33]/var.kna1)*(1.0+x[33]/var.kna2);
	ncxi.h5 = x[33]*x[33]/(ncxi.h4*var.kna1*var.kna2);
	ncxi.h6 = 1.0/ncxi.h4;
	ncxi.h7 = 1.0+(var.nao/var.kna3)*(1.0+1.0/var.hna);
	ncxi.h8 = var.nao/(var.kna3*var.hna*ncxi.h7);
	ncxi.h9 = 1.0/ncxi.h7;
	ncxi.h10 = var.kasym+1.0+(var.nao/var.kna1)*(1.0+var.nao/var.kna2);
	ncxi.h11 = var.nao*var.nao/(ncxi.h10*var.kna1*var.kna2);
	ncxi.h12 = 1.0/ncxi.h10;

	ncxi.k1 = ncxi.h12*var.cao*var.kca_on;
	ncxi.k2 = var.kca_off;
	ncxi.k31 = ncxi.h9*var.omega_ca;
	ncxi.k32 = ncxi.h8*var.omega_naca;
	ncxi.k3 = ncxi.k31 + ncxi.k32;
	ncxi.k41 = ncxi.h3*var.omega_ca/var.hca;
	ncxi.k42 = ncxi.h2*var.omega_naca;
	ncxi.k4 = ncxi.k41 + ncxi.k42;
	ncxi.k5 = var.kca_off;
	ncxi.k6 = ncxi.h6*x[37]*var.kca_on;
	ncxi.k7 = ncxi.h5*ncxi.h2*var.omega_na;
	ncxi.k8 = ncxi.h8*ncxi.h11*var.omega_na;
	
	ncxi.x1 = ncxi.k2*ncxi.k4*(ncxi.k7+ncxi.k6)+ncxi.k5*ncxi.k7*(ncxi.k2+ncxi.k3);
	ncxi.x2 = ncxi.k1*ncxi.k7*(ncxi.k4+ncxi.k5)+ncxi.k4*ncxi.k6*(ncxi.k1+ncxi.k8);
	ncxi.x3 = ncxi.k1*ncxi.k3*(ncxi.k7+ncxi.k6)+ncxi.k8*ncxi.k6*(ncxi.k2+ncxi.k3);
	ncxi.x4 = ncxi.k2*ncxi.k8*(ncxi.k4+ncxi.k5)+ncxi.k3*ncxi.k5*(ncxi.k1+ncxi.k8);

	ncxi.E1 = ncxi.x1/(ncxi.x1+ncxi.x2+ncxi.x3+ncxi.x4);	
	ncxi.E2 = ncxi.x2/(ncxi.x1+ncxi.x2+ncxi.x3+ncxi.x4);	
	ncxi.E3 = ncxi.x3/(ncxi.x1+ncxi.x2+ncxi.x3+ncxi.x4);	
	ncxi.E4 = ncxi.x4/(ncxi.x1+ncxi.x2+ncxi.x3+ncxi.x4);	

	ncxi.allo = 1.0/(1.0+(var.km_ca_act/x[37])*(var.km_ca_act/x[37]));

	ncxi.jnaca_na = 3.0*(ncxi.E4*ncxi.k7 - ncxi.E1*ncxi.k8) + ncxi.E3*ncxi.k42 - ncxi.E2*ncxi.k32;
	ncxi.jnaca_ca = ncxi.E2*ncxi.k2 - ncxi.E1*ncxi.k1;
	var.inaca_i = var.Gnaca*0.8*ncxi.allo*(var.zna*ncxi.jnaca_na + var.zca*ncxi.jnaca_ca);
	
	// subspace
	ncxss.h1 = 1.0 + (x[34]/var.kna3)*(1.0 + var.hna);
	ncxss.h2 = x[34]*var.hna/(var.kna3*ncxss.h1);
	ncxss.h3 = 1.0/ncxss.h1;
	ncxss.h4 = 1.0+x[34]/var.kna1*(1.0+x[34]/var.kna2);
	ncxss.h5 = x[34]*x[34]/(ncxss.h4*var.kna1*var.kna2);
	ncxss.h6 = 1.0/ncxss.h4;
	ncxss.h7 = 1.0+var.nao/var.kna3*(1.0+1.0/var.hna);
	ncxss.h8 = var.nao/(var.kna3*var.hna*ncxss.h7);
	ncxss.h9 = 1.0/ncxss.h7;
	ncxss.h10 = var.kasym+1.0+var.nao/var.kna1*(1.0+var.nao/var.kna2);
	ncxss.h11 = var.nao*var.nao/(ncxss.h10*var.kna1*var.kna2);
	ncxss.h12 = 1.0/ncxss.h10;

	ncxss.k1 = ncxss.h12*var.cao*var.kca_on;
	ncxss.k2 = var.kca_off;
	ncxss.k31 = ncxss.h9*var.omega_ca;
	ncxss.k32 = ncxss.h8*var.omega_naca;
	ncxss.k3 = ncxss.k31 + ncxss.k32;
	ncxss.k41 = ncxss.h3*var.omega_ca/var.hca;
	ncxss.k42 = ncxss.h2*var.omega_naca;
	ncxss.k4 = ncxss.k41 + ncxss.k42;
	ncxss.k5 = var.kca_off;
	ncxss.k6 = ncxss.h6*x[38]*var.kca_on;
	ncxss.k7 = ncxss.h5*ncxss.h2*var.omega_na;
	ncxss.k8 = ncxss.h8*ncxss.h11*var.omega_na;

	ncxss.x1 = ncxss.k2*ncxss.k4*(ncxss.k7+ncxss.k6)+ncxss.k5*ncxss.k7*(ncxss.k2+ncxss.k3);
	ncxss.x2 = ncxss.k1*ncxss.k7*(ncxss.k4+ncxss.k5)+ncxss.k4*ncxss.k6*(ncxss.k1+ncxss.k8);
	ncxss.x3 = ncxss.k1*ncxss.k3*(ncxss.k7+ncxss.k6)+ncxss.k8*ncxss.k6*(ncxss.k2+ncxss.k3);
	ncxss.x4 = ncxss.k2*ncxss.k8*(ncxss.k4+ncxss.k5)+ncxss.k3*ncxss.k5*(ncxss.k1+ncxss.k8);

	ncxss.E1 = ncxss.x1/(ncxss.x1+ncxss.x2+ncxss.x3+ncxss.x4);	
	ncxss.E2 = ncxss.x2/(ncxss.x1+ncxss.x2+ncxss.x3+ncxss.x4);	
	ncxss.E3 = ncxss.x3/(ncxss.x1+ncxss.x2+ncxss.x3+ncxss.x4);	
	ncxss.E4 = ncxss.x4/(ncxss.x1+ncxss.x2+ncxss.x3+ncxss.x4);	

	ncxss.allo = 1.0/(1.0+(var.km_ca_act/x[38])*(var.km_ca_act/x[38]));

	ncxss.jnaca_na = 3.0*(ncxss.E4*ncxss.k7 - ncxss.E1*ncxss.k8) + ncxss.E3*ncxss.k42 - ncxss.E2*ncxss.k32;
	ncxss.jnaca_ca = ncxss.E2*ncxss.k2 - ncxss.E1*ncxss.k1;
	var.inaca_ss = var.Gnaca*0.2*ncxss.allo*(var.zna*ncxss.jnaca_na + var.zca*ncxss.jnaca_ca);

    var.inaca = var.inaca_i+var.inaca_ss;
    
}

// Sodium-Potassium Pump

void comp_inak (double x[])
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	inak.knai = inak.Tknai[iV]*d2 + inak.Tknai[iV+1]*d1;
	inak.knao = inak.Tknao[iV]*d2 + inak.Tknao[iV+1]*d1;

	inak.P = inak.SigP/(1.0 + (inak.H/inak.HP) + (x[33]/inak.nap) + (x[35]/inak.kp));

	inak.a1 = inak.kp1*(x[33]/inak.knai)*(x[33]/inak.knai)*(x[33]/inak.knai)/
				((1.0+x[33]/inak.knai)*(1.0+x[33]/inak.knai)*(1.0+x[33]/inak.knai)+(1+x[35]/inak.kki)*(1+x[35]/inak.kki)-1.0);
	inak.b2 = inak.km2*(var.nao/inak.knao)*(var.nao/inak.knao)*(var.nao/inak.knao)/
	 			((1.0+var.nao/inak.knao)*(1.0+var.nao/inak.knao)*(1.0+var.nao/inak.knao)+(1.0+var.ko/inak.kko)*(1.0+var.ko/inak.kko)-1.0);
	inak.a3 = inak.kp3*(var.ko/inak.kko)*(var.ko/inak.kko)/
				((1.0+var.nao/inak.knao)*(1.0+var.nao/inak.knao)*(1.0+var.nao/inak.knao)+(1.0+var.ko/inak.kko)*(1.0+var.ko/inak.kko)-1.0);
	inak.b3 = inak.km3*inak.P*inak.H/(1.0+inak.MgATP/inak.k_MgATP);
	inak.b4 = inak.km4*(x[35]/inak.kki)*(x[35]/inak.kki)/
				((1.0+x[33]/inak.knai)*(1.0+x[33]/inak.knai)*(1.0+x[33]/inak.knai)+(1.0+x[35]/inak.kki)*(1.0+x[35]/inak.kki)-1.0);

	inak.x1 = inak.a4*inak.a1*inak.a2+inak.b2*inak.b4*inak.b3+inak.a2*inak.b4*inak.b3+inak.b3*inak.a1*inak.a2;
	inak.x2 = inak.b2*inak.b1*inak.b4+inak.a1*inak.a2*inak.a3+inak.a3*inak.b1*inak.b4+inak.a2*inak.a3*inak.b4;
	inak.x3 = inak.a2*inak.a3*inak.a4+inak.b3*inak.b2*inak.b1+inak.b2*inak.b1*inak.a4+inak.a3*inak.a4*inak.b1;
	inak.x4 = inak.b4*inak.b3*inak.b2+inak.a3*inak.a4*inak.a1+inak.b2*inak.a4*inak.a1+inak.b3*inak.b2*inak.a1;

	inak.E1 = inak.x1/(inak.x1+inak.x2+inak.x3+inak.x4);
	inak.E2 = inak.x2/(inak.x1+inak.x2+inak.x3+inak.x4);
	inak.E3 = inak.x3/(inak.x1+inak.x2+inak.x3+inak.x4);
	inak.E4 = inak.x4/(inak.x1+inak.x2+inak.x3+inak.x4);

	inak.jna = 3.0*(inak.E1*inak.a3 - inak.E2*inak.b3);
	inak.jk = 2.0*(inak.E4*inak.b1 - inak.E3*inak.a1);

	inak.inak = inak.G*(var.zna*inak.jna + var.zk*inak.jk);

}

// Sarcolemmal Ca Pump 

void comp_ipca (double x[])
{

	ipca.ca = ipca.G*x[37]/(ipca.km + x[37]);

}

// K Background Current
void comp_ikb (double x[])
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ikb.xkb = ikb.Txkb[iV]*d2 + ikb.Txkb[iV+1]*d1;

	ikb.k = ikb.G*ikb.xkb*(x[0] - var.Ek);
}

// Ca Background Current 

void comp_icab (double x[])
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	icab.exp = icab.Texp[iV]*d2 + icab.Texp[iV+1]*d1;
	
	icab.ca = icab.pcab*4.0*x[0]*F*F/R/T*(icab.gacai*x[37]*icab.exp-icab.gacao*var.cao)/(icab.exp-1.0);

}

// Na Background Current 

void comp_inab (double x[])
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	inab.exp = inab.Texp[iV]*d2 + inab.Texp[iV+1]*d1;

	inab.na = inab.pnab*x[0]*F*F/R/T*(x[33]*inab.exp-var.nao)/(inab.exp-1.0);

}

void comp_CaMK (double x[])
{
	CaMK.bound = CaMK.z*(1.0-x[30])/(1.0+CaM.Km/x[38]);

	CaMK.active = CaMK.bound + x[30];

}

void comp_diff (double x[])
{

	jdiff.na = (x[34]-x[33])/jdiff.tau_na;
	jdiff.k = (x[36]-x[35])/jdiff.tau_k;
	jdiff.ca = (x[38]-x[37])/jdiff.tau_ca;

}

void comp_jrel (double x[])
{

	int flag1=0,flag2=0;

	jrel.NPss = jrel.p*(-jrel.a*ical.ica/(1.0+pow((1.5/x[40]),8.0)));
	jrel.tau_NP = jrel.b_tau/(1.0+(0.0123/x[40]));
	if (jrel.tau_NP < 0.001){
		jrel.tau_NP = 0.001;
		flag1=1;
	}

	jrel.CaMKss = jrel.p*(-jrel.a_CaMK*ical.ica/(1.0+pow((1.5/x[40]),8.0)));
	jrel.tau_CaMK = jrel.b_tau_CaMK/(1.0+(0.0123/x[40]));
	if (jrel.tau_CaMK < 0.001){
		jrel.tau_CaMK = 0.001;
		flag1=2;
	}

	jrel.pCaMK = 1.0/(1.0+CaMK.Km/CaMK.active);

	jrel.ca = (1.0-jrel.pCaMK)*x[31]+jrel.pCaMK*x[32];
	
	if(flag1!=0 || flag2!=0){
		if(flag1==1){
			printf("jrel tauNP is small\n");
		} else if(flag2==1){
			printf("jrel tauNP_CaMK is small\n");
		}
	}
		
}

void comp_jup (double x[])
{

	jup.np = jup.p*(0.004375*x[37]/(0.00092+x[37]));
	jup.CaMK = jup.p*((1.0+jup.dCaMK)*0.004375*x[37]/(0.00092 - jup.dKm_PLB + x[37]));
	
	jup.pCaMK = 1.0/(1.0 + CaMK.Km/CaMK.active);

	jup.leak = 0.0039375*x[39]/15.0;
	if(var.fail==1){
		jup.leak *= 3.5;
	}

	jup.ca = (1.0-jup.pCaMK)*jup.np + jup.pCaMK*jup.CaMK - jup.leak;

}

void comp_jtr (double x[])
{

	jtr.ca = (x[39]-x[40])/jtr.tau;

}

void comp_concentration (double x[])
{
	var.b_Ca_i = 1.0/(1.0+var.cmdnbar*var.kmcmdn/((var.kmcmdn+x[37])*(var.kmcmdn+x[37]))+var.trpnbar*var.kmtrpn/((var.kmtrpn+x[37])*(var.kmtrpn+x[37])));
	var.b_Ca_ss = 1.0/(1.0+var.bsrbar*var.kmbsr/((var.kmbsr+x[38])*(var.kmbsr+x[38]))+var.bslbar*var.kmbsl/((var.kmbsl+x[38])*(var.kmbsl+x[38])));
	var.b_Ca_jsr = 1.0/(1.0+var.csqnbar*var.kmcsqn/((var.kmcsqn+x[40])*(var.kmcsqn+x[40])));

}


// Reversal potentials */

void comp_reversal_potential(double x[])
{
	var.Ena = var.RTonF*log(var.nao/x[33]);
	var.Ek = var.RTonF*log(var.ko/x[35]);
	var.Eks = var.RTonF*log((var.ko+var.prnak*var.nao)/(x[35]+var.prnak*x[33]));
	
	//printf("Ena=%lf, Ek=%lf, Eks=%lf\n",var.Ena,var.Ek,var.Eks);
}

