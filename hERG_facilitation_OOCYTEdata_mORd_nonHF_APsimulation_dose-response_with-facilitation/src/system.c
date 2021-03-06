#include "syspara.h"

void function(double x[],double f[],double t)
{

	int i;
	comp_reversal_potential(x);
	comp_CaMK(x);
	comp_ina(x);
	comp_ito(x);
	comp_ical(x);
	comp_ikr(x);
	comp_ikr2(x);
	comp_iks(x);
	comp_ik1(x);
	comp_inaca(x);
	comp_inak(x);
	comp_ipca(x);
	comp_ikb(x);
	comp_icab(x);
	comp_inab(x);
	//comp_CaMK(x);
	comp_diff(x);
	comp_jrel(x);
	comp_jup(x);
	comp_jtr(x);
	comp_concentration(x);
	
	var.Ina_i_total = ina.total + inab.na + 3.0*inak.inak + 3.0*var.inaca_i;
	var.Ina_ss_total = ical.icana + 3.0*var.inaca_ss;
	var.Ik_i_total = ito.ik + ikr.ik + ikr2.ik + iks.ik + ik1.ik + ikb.k - 2.0*inak.inak + var.Istim;
	var.Ik_ss_total = ical.icak;
	var.Ica_i_total = ipca.ca + icab.ca - 2.0*var.inaca_i;
	var.Ica_ss_total = ical.ica - 2.0*var.inaca_ss;
	var.Itotal = ina.total + ito.ik + ical.ica + ical.icana + ical.icak + ikr.ik + ikr2.ik + iks.ik + ik1.ik 
					+ var.inaca + inak.inak + inab.na + icab.ca + ikb.k + ipca.ca + var.Istim;

	f[0] = -var.Itotal;
	//Fast sodium current
	f[1] = (ina.mss - x[1])/ina.taum; // m
	f[2] = (ina.hss - x[2])/ina.tauh_fast; // h_fast
	f[3] = (ina.hss - x[3])/ina.tauh_slow; // h_slow
	f[4] = (ina.hss - x[4])/ina.tauj; // j
	f[5] = (ina.hCaMKss - x[5])/ina.tauh_CaMK_slow; // h_CaMK_slow
	f[6] = (ina.hss - x[6])/ina.tauj_CaMK; // j_CaMK
	//late sodium current
	f[7] = (ina.mlss - x[7])/ina.tauml; // ml
	f[8] = (ina.hlss - x[8])/ina.tauhl; // hl
	f[9] = (ina.hlCaMKss - x[9])/ina.tauhl_CaMK; // hl
	//Transient outward current
	f[10] = (ito.ass - x[10])/ito.taua;
	f[11] = (ito.iss - x[11])/ito.taui_fast;
	f[12] = (ito.iss - x[12])/ito.taui_slow;
	f[13] = (ito.aCaMKss - x[13])/ito.taua;
	f[14] = (ito.iss - x[14])/ito.taui_CaMK_fast;
	f[15] = (ito.iss - x[15])/ito.taui_CaMK_slow;
	// LTCC
	f[16] = (ical.dss - x[16])/ical.taud;
	f[17] = (ical.fss - x[17])/ical.tauf_fast;
	f[18] = (ical.fss - x[18])/ical.tauf_slow;
	f[19] = (ical.fss - x[19])/ical.taufca_fast;
	f[20] = (ical.fss - x[20])/ical.taufca_slow;
	f[21] = (ical.fss - x[21])/ical.taujca;
	f[22] = (ical.fss - x[22])/ical.tauf_CaMK_fast;
	f[23] = (ical.fss - x[23])/ical.taufca_CaMK_fast;
	f[24] = ical.alpha_n*ical.kp2n - x[24]*ical.km2n;
	// Ikr
	f[25] = (ikr.xrss - x[25])/ikr.tauxr_fast;
	f[26] = (ikr.xrss - x[26])/ikr.tauxr_slow;
	// Iks
	f[27] = (iks.xs1ss - x[27])/iks.tauxs1;
	f[28] = (iks.xs1ss - x[28])/iks.tauxs2;
	// Ik1
	f[29] = (ik1.k1ss - x[29])/ik1.tauk1;
	// CaMK
	f[30] = CaMK.a*CaMK.bound*(CaMK.bound+x[30]) - CaMK.b*x[30];
	// Jrel
	f[31] = (jrel.NPss - x[31])/jrel.tau_NP; 
	f[32] = (jrel.CaMKss - x[32])/jrel.tau_CaMK; 
	// [Na]i
	f[33] = -var.Ina_i_total*var.vr1 + jdiff.na*var.vr2;
	// [Na]ss
	f[34] = -var.Ina_ss_total*var.vr3 - jdiff.na;
	// [K]i
	f[35] = -var.Ik_i_total*var.vr1 + jdiff.k*var.vr2;
	// [K]ss
	f[36] = -var.Ik_ss_total*var.vr3 - jdiff.k;
	// [Ca]i
	f[37] = var.b_Ca_i*(-var.Ica_i_total*var.vr4 - jup.ca*var.vr5 + jdiff.ca*var.vr2);
	// [Ca]ss
	f[38] = var.b_Ca_ss*(-var.Ica_ss_total*var.vr6 +jrel.ca*var.vr7 - jdiff.ca);
	// [Ca]nsr
	f[39] = jup.ca - jtr.ca*var.vr8;
	// [Ca]jsr
	f[40] = var.b_Ca_jsr*(jtr.ca-jrel.ca);
	// Ikr2 (facilitation component)
	f[41] = (ikr2.xrss - x[41])/ikr2.tauxr_fast;
	f[42] = (ikr2.xrss - x[42])/ikr2.tauxr_slow;

	//printf("NPss=%lf,NPp=%lf\n",jrel.NPss,jrel.CaMKss);
	//for(i=0;i<NN;i++){
	//	printf("x[%d]=%e\n",i,f[i]);
	//}
}
