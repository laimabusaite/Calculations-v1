// Saakam programmeet, iesaakumaa, meegjinaasim paarrakstiit mathematica programmu uz c++
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <complex>
#define Pi 3.1415926535897932384626433
#define k_Boltz 1.3806504e-23
#define e 2.718281828459
#define cc 2.99792458e8
#define mi_Bora_MHz 1.3996 //4 //E[MHz]=B[G]*g*mf*koef_E_MHz
#include <limits>
#include <algorithm>
#include "viela.cc"
#include "viela_probe.cc"
#include <sys/stat.h>
#define prec pow(10.0,-10)
//#include <dcomplex>
using namespace std;
typedef std::complex<double> dcomplex;
//imaginaaraa vieniiba
gsl_complex imagunit = gsl_complex_rect(0.0,1.0);
dcomplex dc_imagunit(0.0,1.0);
// mazaaki par so buus nulle
//
// lenjkjiskie momenti F_gr un F_ex
//
//dimensijas
int dimNum = FDeltaGr * (FStartGr + FEndGr + 1) * FDeltaGr * (FStartGr + FEndGr + 1) + FDeltaEx * (FStartEx + FEndEx + 1) * FDeltaEx * (FStartEx + FEndEx + 1);
int dimNum_probe = FDeltaGr_probe * (FStartGr_probe + FEndGr_probe + 1) * FDeltaGr_probe * (FStartGr_probe + FEndGr_probe + 1) + FDeltaEx_probe * (FStartEx_probe + FEndEx_probe + 1) * FDeltaEx_probe * (FStartEx_probe + FEndEx_probe + 1);
// 3j * 6j simbolu matricas
//vietas
int vietaGr[10][20][10][20];
int vietaEx[10][20][10][20];
int vietaMxGr[10][20];
int vietaMxEx[10][20];
// lvs matricas
gsl_matrix_complex * ppMatrix = gsl_matrix_complex_calloc(dimNum, dimNum);
gsl_vector_complex * resVector = gsl_vector_complex_calloc(dimNum);
gsl_vector_complex * ppVector = gsl_vector_complex_calloc(dimNum);
gsl_matrix_complex * rotated;
gsl_matrix * dMatrixE = gsl_matrix_calloc(4, 4);
gsl_matrix * dMatrixG = gsl_matrix_calloc(2, 2);
//dcomplex ppMatrix[dimNum][dimNum];
//dcomplex resVector[dimNum];
//
// matricas, kuraas glabaasim bliivuma matricas
#define Edim 16
dcomplex rhoE[Edim][Edim];
#define Ecnt 16
dcomplex rhEP[Ecnt][Edim][Edim];
#define Gdim 8
dcomplex rhoG[Gdim][Gdim];
#define Gcnt 2
dcomplex rhGP[Gcnt][Gdim][Gdim];
dcomplex JJmatrix[2][4][4];
double FtoJmx[2][4][4][4][9][4][9];
// ieviesham papildus relaksaaciju augshai
double Gamma_exc = Gamma_n * 10;
//
// parametri
//
double Rabi = 2; //MHz
double Temp = 298;
//Probe
double energija_probe = pow(10.0, -6.0) * (cc / lambda_probe), vilnu_vektors_probe = 2 * Pi / lambda_probe, lazera_frekvence_probe = energija_probe;
double DWidth_probe = sqrt((8 * k_Boltz * Temp * log(2))/(masa_probe * cc * cc)) * lazera_frekvence_probe;
double Sigma_probe = (1.0/1.0) * sqrt(k_Boltz * Temp / masa_probe) * lazera_frekvence_probe / cc;
//Pump
double energija = pow(10.0, -6.0) * (cc / lambda), vilnu_vektors = 2 * Pi / lambda, lazera_frekvence = energija;
double DWidth = sqrt((8 * k_Boltz * Temp * log(2))/(masa * cc * cc)) * lazera_frekvence;
double atrums = sqrt((8 * k_Boltz * Temp)/(masa * Pi));
//double Sigma = Gamma_n * 10;
double Sigma = (1.0/1.0) * sqrt(k_Boltz * Temp / masa) * lazera_frekvence / cc;
double constDstep = 1.00;
double DScan = 2*Sigma*constDstep;
//double DScan = 56;
// int DSteps = 150*constDstep; //Doplera solji
int DSteps = 150*constDstep; //Doplera solji
double DStep = DScan / DSteps;
double detuning = 0, detuning_probe = 0; //MHz
double E_p, liin_platums = 2, liin_platums_probe = 2; //MHz
double Gamma_mazs, Gamma_mazs_probe; //MHz
double Gamma_p;
double BRange = 3000; //300; //3000; //Maksimala lauka vertiba
double DRange = 6300;
#define steps 1000 //solu skaits magnetiskaja lauka
double step = BRange / steps;
//double step = DRange / steps;
double energySplitsMatrix[2][10][20][2 * (steps + 1)];
double energySplitsMatrixProbe[2][10][20][2 * (steps + 1)];
double levelMixingMatrix[2][10][10][20][2 * (steps + 1)];
double levelMixingMatrixProbe[2][10][10][20][2 * (steps + 1)];
double apdz[2][10][2 * (steps + 1)];
double sabruksana[2][10][20][10][20][2 * (steps + 1)];
double three_j_matrix[10][20][10][20];
double three_j_matrix_probe[10][20][10][20];
double reduced_matrix[10][10];
double reduced_matrix_probe[10][10];
double Btrans;
int tmpCnt;
//
// sadursmju parametri
//
double pressure = 5e-6 * 133.322; // Pa
double sigma_sad = 3e-18; // m^2 sadursmes skeersgriezums
double Delta_col;
//
//ierosmes gjeometrija	//SVARIGI!!!!!!!
double delta_teta_i = 0.0 * Pi/180;
//
int pol = 0; //0 - lineaara; 1 = pa labi; -1 = pa kreisi
double teta_i = Pi/2; //
double fii_i = 0;
//
// ierosmes noskanjoshana	//SVARIGI!!!!!!!
int Fg_d = 3;// 4; //3;//4; //1; //2; //1; //
int Fe_d = 3; //3; //3; //3;//1; //1; //3; //1; //
//
//
// zondes noskanjoshana	//SVARIGI!!!!!!!
int Fg_p = Fg_d; //4; //3;//2; //1; //
int Fe_p = Fg_d; //4;//4; //4; //1; //3; //1; //
//
//noveerosanas geometrija	//SVARIGI!!!!!!!
double delta_teta = 0.0 * Pi/180;
double delta_fii = 0.0 * Pi/180;
//
int pol_o[2] = {0, 0};  //{0, 0}; //0 - lineaara; 1 - pa labi; -1 - pa kreisi //pol_o[2] nozime ka output bus 2 komponetnes
double teta_o[2] = {Pi/2 + delta_teta, Pi/4 + delta_teta}; //ar komatu atdala pirmo un otro komponenti
double fii_o[2] = {Pi/4 + delta_fii, 0 + delta_fii}; //{0, Pi/2}; //
//
//zondējošā stara ģeometrija	//SVARIGI!!!!!!!
//
int pol_p[2] = {0, 0};  //{0, 0}; //0 - lineaara; 1 - pa labi; -1 - pa kreisi //pol_o[2] nozime ka output bus 2 komponetnes
double teta_p[2] = {Pi/2 + delta_teta, Pi/4 + delta_teta}; //{0, Pi/2}; // ar komatu atdala pirmo un otro komponenti
double fii_p[2] = {Pi/4 + delta_fii, 0 + delta_fii}; //{0, Pi/2}; //{0, Pi/2}; //{Pi, Pi}; //{Pi, Pi}; //

// E kompleksie vektori
//
// Matricaa U transformaacijai no Dekarta uz ciklisko koordinaatu sisteemu (Auzinsh, Budker, Rochester (3.9), kuru njem kompleksi saistiitu saskanjaa ar (3.15), p. 32), matrica nav atkariiga no E vektora izveeles, taapeec to defineejam globaali
gsl_matrix_complex * U = gsl_matrix_complex_calloc(3, 3);

void init_u(){ // Saliek pareizaas veertiibas U matricaa
    gsl_matrix_complex_set(U, 0, 0, gsl_complex_rect(-1/sqrt(2), 0));
    gsl_matrix_complex_set(U, 0, 1, gsl_complex_rect(0, 1/sqrt(2)));
    gsl_matrix_complex_set(U, 1, 2, gsl_complex_rect(1, 0));
    gsl_matrix_complex_set(U, 2, 0, gsl_complex_rect(1/sqrt(2), 0));
    gsl_matrix_complex_set(U, 2, 1, gsl_complex_rect(0, 1/sqrt(2)));
}

gsl_complex E(int q, int pol_u, double teta_u, double fii_u){
/*
    Shiis funkcijas interfeiss nav mainiits, ir taads pats kaa ieprieksh, tomeer izsaukt to katru reizi kad ievajagaas kaadu E komponenti sanaak neefektiivi, jo katru reizi tiek iziets pils transformaaciju celjsh
*/
    if(abs(q) < 1.5) // ja q par lielu, tad nav jeegas neko dariit
    {
        // inicializeejam saakotneejo vektoru uz nulleem, peec tam saliek veertiibas
        gsl_vector_complex * initial = gsl_vector_complex_calloc(3);
        if(pol_u == 1) // Pa kreisi polarizeets 1/sqrt(2)(x + iy)
        {
            gsl_vector_complex_set(initial, 0, gsl_complex_rect(1/sqrt(2), 0));
            gsl_vector_complex_set(initial, 1, gsl_complex_rect(0, 1/sqrt(2)));
        }
        else if(pol_u == 0) // Lineaars z
        {
            gsl_vector_complex_set(initial, 2, gsl_complex_rect(1, 0));
        }
        else if(pol_u == -1) // Pa labi polarizeets 1/sqrt(2)(x - iy)
        {
            gsl_vector_complex_set(initial, 0, gsl_complex_rect(1/sqrt(2), 0));
            gsl_vector_complex_set(initial, 1, gsl_complex_rect(0, -1/sqrt(2)));
        }

        // Pagrieziena matrica ap y asi par teta
        gsl_matrix_complex * R1 = gsl_matrix_complex_calloc(3, 3);
        gsl_matrix_complex_set(R1, 0, 0, gsl_complex_rect(cos(teta_u), 0));
        gsl_matrix_complex_set(R1, 0, 2, gsl_complex_rect(sin(teta_u), 0));
        gsl_matrix_complex_set(R1, 1, 1, gsl_complex_rect(1, 0));
        gsl_matrix_complex_set(R1, 2, 0, gsl_complex_rect(-sin(teta_u), 0));
        gsl_matrix_complex_set(R1, 2, 2, gsl_complex_rect(cos(teta_u), 0));

        // Pagrieziena matrica ap z asi par fii
        gsl_matrix_complex * R2 = gsl_matrix_complex_calloc(3, 3);
        gsl_matrix_complex_set(R2, 0, 0, gsl_complex_rect(cos(fii_u), 0));
        gsl_matrix_complex_set(R2, 0, 1, gsl_complex_rect(-sin(fii_u), 0));
        gsl_matrix_complex_set(R2, 1, 0, gsl_complex_rect(sin(fii_u), 0));
        gsl_matrix_complex_set(R2, 1, 1, gsl_complex_rect(cos(fii_u), 0));
        gsl_matrix_complex_set(R2, 2, 2, gsl_complex_rect(1, 0));

        // Kopeejaa transformaacijas matrica R = R2 * R1 (transformaacijas matrica iedarbojas uz vektoru no kreisaas puses, taapeec shaada reizinaataaju seciiba)
        gsl_matrix_complex * R = gsl_matrix_complex_calloc(3, 3);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1, 0), R2, R1, gsl_complex_rect(0, 0), R);

        // Dekarta vektora transformaacija
        gsl_vector_complex * rotated_Cart = gsl_vector_complex_calloc(3);
        gsl_blas_zgemv(CblasNoTrans, gsl_complex_rect(1, 0), R, initial, gsl_complex_rect(0, 0), rotated_Cart);

        // Paareja no Dekarta uz cikliskajaam koordinaataam
        gsl_vector_complex * cyclic = gsl_vector_complex_calloc(3);
        gsl_blas_zgemv(CblasNoTrans, gsl_complex_rect(1, 0), U, rotated_Cart, gsl_complex_rect(0, 0), cyclic);

        // Izveelamies veertiibu, ko atgriezt
        gsl_complex retVal = gsl_vector_complex_get(cyclic, -q+1);

        // Atbriivojam atminju
        gsl_vector_complex_free(initial);
        gsl_vector_complex_free(rotated_Cart);
        gsl_vector_complex_free(cyclic);

        gsl_matrix_complex_free(R1);
        gsl_matrix_complex_free(R2);
        gsl_matrix_complex_free(R);
        return retVal;
    }
    else
    {
        return gsl_complex_rect(0,0);
    }
}
//
// paarbauda triistuura likumus un atgriezh 3j simbolu vai nulli
//
double three_j_wrap(int ja, int jb, int jc, int ma, int mb, int mc){
	if(abs(ja - jb) <= jc && ja + jb >= jc &&  ma + mb + mc == 0){
		return gsl_sf_coupling_3j(ja, jb, jc, ma, mb, mc);
		}
	else{
		return 0;
		}
	}
//
// paarbauda triistuura likumus un atgriezh 6j simbolu vai nulli
//
double six_j_wrap(int ja, int jb, int jc, int jd, int je, int jf){
	if(abs(ja - jb) <= jc && ja + jb >= jc && abs(ja - je) <= jf && ja + je >= jf && abs(jd - jb) <= jf && jd + jb >= jf && abs(jd - je) <= jc && jd + je >= jc){
		return gsl_sf_coupling_6j(ja, jb, jc, jd, je, jf);
		}
	else{
		return 0;
		}
	}
//
// paareja uz Je - Jg reduceeto elementu
//
double reducedJ(int fe, int fg, double Je, double Jg){
	double retVal;
	retVal = pow(-1.0, Je + ik + fg + 1) * sqrt((2 * fe + 1) * (2 * fg + 1)) * six_j_wrap(int(Je * 2), fe * 2, int(ik * 2), fg * 2, int(Jg * 2), 1 * 2);
	return retVal;
	}
double reducedJ_probe(int fe, int fg, double Je, double Jg){
	double retVal;
	retVal = pow(-1.0, Je + ik_probe + fg + 1) * sqrt((2 * fe + 1) * (2 * fg + 1)) * six_j_wrap(int(Je * 2), fe * 2, int(ik_probe * 2), fg * 2, int(Jg * 2), 1 * 2);
	return retVal;
	}
//
// optimizeejam: f-ja, kas saliek visus iespeejamos 3j simbolus matricaa
//
void build_three_j_matrix(){
	int FStart = min(FStartEx, FStartGr);
	int FEnd = max(FEndEx, FEndGr);
	int F, Fp, M, Mp, q, i, j, ip, jp;
	for(i = 0; i < 10; i++){
		for(j = -9; j <= 9; j++){
			for(ip = 0; ip < 10; ip++){
				for(jp = -9; jp <= 9; jp++){
					three_j_matrix[i][j][ip][jp] = 0;
				}
			}
		}
	}
	for(F = FStart; F <= FEnd; F++){
		for(M = -F; M <= F; M++){
			for(Fp = FStart; Fp <= FEnd; Fp++){
				for(Mp = -Fp; Mp <= Fp; Mp++){
					q = - M - Mp;
					three_j_matrix[F][M][Fp][Mp] = three_j_wrap(2 * F, 2, 2 * Fp, 2 * M, 2 * q, 2 * Mp);
					//cout << "(F,M||d||Fp,Mp) (" << F << "," << M << "||" << Fp << "," << Mp << ") " << three_j_wrap(2 * F, 2, 2 * Fp, 2 * M, 2 * q, 2 * Mp) << " " << three_j_matrix[F][M][Fp][Mp] << endl;
				}
			}
		}
	}
}
void build_three_j_matrix_probe(){
	int FStart = min(FStartEx_probe, FStartGr_probe);
	int FEnd = max(FEndEx_probe, FEndGr_probe);
	int F, Fp, M, Mp, q, i, j, ip, jp;
	for(i = 0; i < 10; i++){
		for(j = -9; j <= 9; j++){
			for(ip = 0; ip < 10; ip++){
				for(jp = -9; jp <= 9; jp++){
					three_j_matrix_probe[i][j][ip][jp] = 0;
				}
			}
		}
	}
	for(F = FStart; F <= FEnd; F++){
		for(M = -F; M <= F; M++){
			for(Fp = FStart; Fp <= FEnd; Fp++){
				for(Mp = -Fp; Mp <= Fp; Mp++){
					q = - M - Mp;
					three_j_matrix_probe[F][M][Fp][Mp] = three_j_wrap(2 * F, 2, 2 * Fp, 2 * M, 2 * q, 2 * Mp);
					//cout << "(F,M||d||Fp,Mp) (" << F << "," << M << "||" << Fp << "," << Mp << ") " << three_j_wrap(2 * F, 2, 2 * Fp, 2 * M, 2 * q, 2 * Mp) << " " << three_j_matrix[F][M][Fp][Mp] << endl;
				}
			}
		}
	}
}
//
// optimizeejam: f-ja, kas saliek visas iespeejamaas paarejas uz reduceetajiem elementiem matricaa
//
void build_reduced_matrix(){
	int FStart = min(FStartEx, FStartGr);
	int FEnd = max(FEndEx, FEndGr);
	int F, Fp, i, ip;
	// cout << FStart << " " << FEnd <<endl;
	for(i = 0; i < 10; i++){
		for(ip = 0; ip < 10; ip++){
			reduced_matrix[i][ip] = 0;
		}
	}
	for(F = FStart; F <= FEnd; F++){
		for(Fp = FStart; Fp <= FEnd; Fp++){
			//if((F == 2 && Fp == 3) || (F == 3 && Fp == 2)) reduced_matrix[F][Fp] = reducedJ(F, Fp, JJ[1], JJ[0]);
			reduced_matrix[F][Fp] = reducedJ(F, Fp, JJ[1], JJ[0]);
			//cout << "(F||d||Fp) (" << F << "||" << Fp << ") " << reducedJ(F, Fp, JJ[1], JJ[0]) << " " << reduced_matrix[F][Fp] << endl;
		}
	}

	FStart = min(FStartEx_probe, FStartGr_probe);
	FEnd = max(FEndEx_probe, FEndGr_probe);
	// cout << FStart << " " << FEnd <<endl;
	// int F, Fp, i, ip;
	for(i = 0; i < 10; i++){
		for(ip = 0; ip < 10; ip++){
			reduced_matrix_probe[i][ip] = 0;
		}
	}
	for(F = FStart; F <= FEnd; F++){
		for(Fp = FStart; Fp <= FEnd; Fp++){
			//if((F == 2 && Fp == 3) || (F == 3 && Fp == 2)) reduced_matrix[F][Fp] = reducedJ(F, Fp, JJ[1], JJ[0]);
			reduced_matrix_probe[F][Fp] = reducedJ_probe(F, Fp, JJ_probe[1], JJ_probe[0]);
			//cout << "(F||d||Fp) (" << F << "||" << Fp << ") " << reducedJ(F, Fp, JJ[1], JJ[0]) << " " << reduced_matrix[F][Fp] << endl;
		}
	}
}
//
//spontaanaas relaksaacijas reekinu funkcija
//
double relax(int fe, int fep, int fg, int fgp, int M, int Mp, int mi, int mip, int step){
    //step = 0;
	double retVal = 0.0, temp1 = 0.0, temp2 = 0.0;
	int q = M - mi;
	if(mi - mip != M - Mp) return 0.0;
	else{
		for(int F1mix = FStartGr; F1mix <= FEndGr; F1mix++){
			for(int F2mix = FStartEx; F2mix <= FEndEx; F2mix++){
				for(int F3mix = FStartGr; F3mix <= FEndGr; F3mix++){
					for(int F4mix = FStartEx; F4mix <= FEndEx; F4mix++){
						retVal +=  FDeltaEx * pow(-1.0, F2mix + F4mix - M - Mp) * levelMixingMatrix[0][fg - FStartGr][F1mix - FStartGr][mi + FEndGr][step] * levelMixingMatrix[1][fe - FStartEx][F2mix - FStartEx][M + FEndEx][step] * levelMixingMatrix[0][fgp - FStartGr][F3mix - FStartGr][mip + FEndGr][step] * levelMixingMatrix[1][fep - FStartEx][F4mix - FStartEx][Mp + FEndEx][step] * three_j_matrix[F2mix][-M][F1mix][mi] * reduced_matrix[F2mix][F1mix] * three_j_matrix[F4mix][-Mp][F3mix][mip] * reduced_matrix[F4mix][F3mix];
						//retVal +=  FDeltaEx * pow(-1.0, F2mix + F4mix - M - Mp) * levelMixingMatrix[0][fg - FStartGr][F1mix - FStartGr][mi + FEndGr][1] * levelMixingMatrix[1][fe - FStartEx][F2mix - FStartEx][M + FEndEx][1] * levelMixingMatrix[0][fgp - FStartGr][F3mix - FStartGr][mip + FEndGr][1] * levelMixingMatrix[1][fep - FStartEx][F4mix - FStartEx][Mp + FEndEx][1] * three_j_matrix[F2mix][-M][F1mix][mi] * reduced_matrix[F2mix][F1mix] * three_j_matrix[F4mix][-Mp][F3mix][mip] * reduced_matrix[F4mix][F3mix];
						//if(three_j_wrap(F2mix * 2, 1 * 2, F1mix * 2, -M * 2, q * 2, mi * 2) != three_j_matrix[F2mix][-M][F1mix][mi]) cout << "(F,M||d||Fp,Mp) (" << F2mix << "," << -M << "||" << F1mix << "," << mi << ") " << three_j_wrap(2 * F2mix, 2, 2 * F1mix, 2 * (-M), 2 * q, 2 * mi) << " " << three_j_matrix[F2mix][-M][F1mix][mi] << endl;
						}
					}
				}
			}
		//if(fe == fep && M == Mp && fg == fgp && mi == mip){
		//sabruksana[1][fe - FStartEx][M + FStartEx][fep - FStartEx][Mp + FStartEx][step] += retVal;
		//}
		retVal *= Gamma_n;
		return retVal;
		}
	}
//
// vienkaarshoti baazes elementiem
//
gsl_complex baseElementSimple(int Fe, int Fg, int M, int m, int step, int pol_u = pol, double teta_u = teta_i, double fii_u = fii_i){
    //step = 0;
	gsl_complex retVal, temp;
	double tempR = 0.0, tempI = 0.0;
	int q = M - m;
	gsl_complex E_cikl = E(q, pol_u, teta_u, fii_u);
	for(int F1mix = FStartEx; F1mix <= FEndEx; F1mix++){
		for(int F2mix = FStartGr; F2mix <= FEndGr; F2mix++){
			temp = gsl_complex_mul_real(E_cikl,  pow(-1.0, F1mix - M) * levelMixingMatrix[1][Fe - FStartEx][F1mix - FStartEx][M + FEndEx][step] * levelMixingMatrix[0][Fg - FStartGr][F2mix - FStartGr][m + FEndGr][step] * three_j_matrix[F1mix][-M][F2mix][m] * reduced_matrix[F1mix][F2mix]);
			//temp = gsl_complex_mul_real(E_cikl,  pow(-1.0, F1mix - M) * levelMixingMatrix[1][Fe - FStartEx][F1mix - FStartEx][M + FEndEx][1] * levelMixingMatrix[0][Fg - FStartGr][F2mix - FStartGr][m + FEndGr][1] * three_j_matrix[F1mix][-M][F2mix][m] * reduced_matrix[F1mix][F2mix]);
			tempR += GSL_REAL(temp);
			tempI += GSL_IMAG(temp);
		}
	}
	//temp = gsl_complex_mul_real(E_cikl, pow(-1.0, fe - M) * three_j_wrap(2 * fe, 2 * 1, 2 * fg, -(2 * M), 2 * q, 2 * m) * reducedJ(fe, fg, JJ[1], JJ[0]));
	retVal = gsl_complex_rect(tempR, tempI);
	return retVal;
}
// induceetais
gsl_complex baseElementSimpleInd(int Fg, int Fe, int m, int M, int step, int pol_u = pol, double teta_u = teta_i, double fii_u = fii_i){
    //step = 0;
	gsl_complex retVal, temp;
	double tempR = 0.0, tempI = 0.0;
	int q = M - m;
	gsl_complex E_cikl = gsl_complex_mul_real(gsl_complex_conjugate(E(q, pol_u, teta_u, fii_u)), pow(-1.0, q));
	for(int F1mix = FStartGr; F1mix <= FEndGr; F1mix++){
		for(int F2mix = FStartEx; F2mix <= FEndEx; F2mix++){
			temp = gsl_complex_mul_real(E_cikl, pow(-1.0, F1mix - m + F1mix - F2mix) * levelMixingMatrix[0][Fg - FStartGr][F1mix - FStartGr][m + FEndGr][step] * levelMixingMatrix[1][Fe - FStartEx][F2mix - FStartEx][M + FEndEx][step] * three_j_matrix[F1mix][-m][F2mix][M] * reduced_matrix[F2mix][F1mix]);
			//temp = gsl_complex_mul_real(E_cikl, pow(-1.0, F1mix - m + F1mix - F2mix) * levelMixingMatrix[0][Fg - FStartGr][F1mix - FStartGr][m + FEndGr][1] * levelMixingMatrix[1][Fe - FStartEx][F2mix - FStartEx][M + FEndEx][1] * three_j_matrix[F1mix][-m][F2mix][M] * reduced_matrix[F2mix][F1mix]);
			//if(reducedJ(F2mix, F1mix, JJ[1], JJ[0]) == 100) cout << "(F,M||d||Fp,Mp) (" << F2mix << "," << -m << "||" << F1mix << "," << M <<") " << reducedJ(F2mix, F1mix, JJ[1], JJ[0]) * three_j_wrap(2 * F1mix, 2 * 1, 2 * F2mix, -(2 * m), -(2 * q), 2 * M) << " " << reduced_matrix[F2mix][F1mix] * three_j_matrix[F1mix][-m][F2mix][M] << endl;
			tempR += GSL_REAL(temp);
			tempI += GSL_IMAG(temp);
			}
		}
	retVal = gsl_complex_rect(tempR, tempI);
	return retVal;
}

// vienkaarshoti baazes elementiem zondejosam staram
//
gsl_complex baseElementSimpleProbe(int Fe, int Fg, int M, int m, int step, int pol_u = pol, double teta_u = teta_i, double fii_u = fii_i){
    //step = 0;
	gsl_complex retVal, temp;
	double tempR = 0.0, tempI = 0.0;
	int q = M - m;
	gsl_complex E_cikl = E(q, pol_u, teta_u, fii_u);
	for(int F1mix = FStartEx_probe; F1mix <= FEndEx_probe; F1mix++){
		for(int F2mix = FStartGr_probe; F2mix <= FEndGr_probe; F2mix++){
			temp = gsl_complex_mul_real(E_cikl,  pow(-1.0, F1mix - M) * levelMixingMatrixProbe[1][Fe - FStartEx_probe][F1mix - FStartEx_probe][M + FEndEx_probe][step] * levelMixingMatrixProbe[0][Fg - FStartGr_probe][F2mix - FStartGr_probe][m + FEndGr_probe][step] * three_j_matrix_probe[F1mix][-M][F2mix][m] * reduced_matrix_probe[F1mix][F2mix]);
			//temp = gsl_complex_mul_real(E_cikl,  pow(-1.0, F1mix - M) * levelMixingMatrix[1][Fe - FStartEx][F1mix - FStartEx][M + FEndEx][1] * levelMixingMatrix[0][Fg - FStartGr][F2mix - FStartGr][m + FEndGr][1] * three_j_matrix[F1mix][-M][F2mix][m] * reduced_matrix[F1mix][F2mix]);
			tempR += GSL_REAL(temp);
			tempI += GSL_IMAG(temp);
		}
	}
	//temp = gsl_complex_mul_real(E_cikl, pow(-1.0, fe - M) * three_j_wrap(2 * fe, 2 * 1, 2 * fg, -(2 * M), 2 * q, 2 * m) * reducedJ(fe, fg, JJ[1], JJ[0]));
	retVal = gsl_complex_rect(tempR, tempI);
	return retVal;
}
// induceetais zondejosam staram
gsl_complex baseElementSimpleIndProbe(int Fg, int Fe, int m, int M, int step, int pol_u = pol, double teta_u = teta_i, double fii_u = fii_i){
    //step = 0;
	gsl_complex retVal, temp;
	double tempR = 0.0, tempI = 0.0;
	int q = M - m;
	gsl_complex E_cikl = gsl_complex_mul_real(gsl_complex_conjugate(E(q, pol_u, teta_u, fii_u)), pow(-1.0, q));
	for(int F1mix = FStartGr_probe; F1mix <= FEndGr_probe; F1mix++){
		for(int F2mix = FStartEx_probe; F2mix <= FEndEx_probe; F2mix++){
			temp = gsl_complex_mul_real(E_cikl, pow(-1.0, F1mix - m + F1mix - F2mix) * levelMixingMatrixProbe[0][Fg - FStartGr_probe][F1mix - FStartGr_probe][m + FEndGr_probe][step] * levelMixingMatrixProbe[1][Fe - FStartEx_probe][F2mix - FStartEx_probe][M + FEndEx_probe][step] * three_j_matrix_probe[F1mix][-m][F2mix][M] * reduced_matrix_probe[F2mix][F1mix]);
			//temp = gsl_complex_mul_real(E_cikl, pow(-1.0, F1mix - m + F1mix - F2mix) * levelMixingMatrix[0][Fg - FStartGr][F1mix - FStartGr][m + FEndGr][1] * levelMixingMatrix[1][Fe - FStartEx][F2mix - FStartEx][M + FEndEx][1] * three_j_matrix[F1mix][-m][F2mix][M] * reduced_matrix[F2mix][F1mix]);
			//if(reducedJ(F2mix, F1mix, JJ[1], JJ[0]) == 100) cout << "(F,M||d||Fp,Mp) (" << F2mix << "," << -m << "||" << F1mix << "," << M <<") " << reducedJ(F2mix, F1mix, JJ[1], JJ[0]) * three_j_wrap(2 * F1mix, 2 * 1, 2 * F2mix, -(2 * m), -(2 * q), 2 * M) << " " << reduced_matrix[F2mix][F1mix] * three_j_matrix[F1mix][-m][F2mix][M] << endl;
			tempR += GSL_REAL(temp);
			tempI += GSL_IMAG(temp);
			}
		}
	retVal = gsl_complex_rect(tempR, tempI);
	return retVal;
}

//
// Funkcijas, kas reekjina Gamma_p koeficientus
//
double split(int Fe, int M, int Fg, int m, int step){
    return energySplitsMatrix[1][Fe - FStartEx][M + FEndEx][step] - energySplitsMatrix[0][Fg - FStartGr][m + FEndGr][step]; // normaala izskanjoshanaas
    // return energySplitsMatrix[1][Fe - FStartEx][M + FEndEx][0] - energySplitsMatrix[0][Fg - FStartGr][m + FEndGr][0]; // izskanjoshanaas efekti "izsleegti", visur izmantojam nulles lauka izskanjosahnos
}
gsl_complex GammaPCalcPlus(int Fe, int m, int Fg, int mi, int step){
	gsl_complex retVal;
	dcomplex temp;
	temp = pow(Rabi, 2.0) / (Gamma_mazs + (Gamma_n/2) + (Gamma_exc) + (liin_platums/2) + dc_imagunit * (detuning - split(Fe, m, Fg, mi, step)));
	retVal = gsl_complex_rect(real(temp), imag(temp));
	//if(Fg == 0) cout << "Fg: " << Fg << "; Fe: " << Fe << "; m: " << m << "; Gp(sauc): " << GSL_REAL(sauc) << " + " << GSL_IMAG(sauc) << endl;
	if(abs(GSL_REAL(retVal)) > prec or abs(GSL_IMAG(retVal)) > prec) return retVal;
	else return gsl_complex_rect(0.0, 0.0);
	}
gsl_complex GammaPCalcMinus(int Fe, int m, int Fg, int mi, int step){
	gsl_complex retVal;
	dcomplex temp;
	//cout << Gamma_mazs << endl;
	temp = pow(Rabi, 2.0) / (Gamma_mazs + (Gamma_n/2) + (Gamma_exc) + (liin_platums/2) - dc_imagunit * (detuning - split(Fe, m, Fg, mi, step)));
	retVal = gsl_complex_rect(real(temp), imag(temp));
	//if(m == 0 && mi == 0) cout << "Fg: " << Fg << "; Fe: " << Fe << "; m: " << m << "; Gp: " << temp << endl;
	//if(Fg == 2 && Fe == 2 && m == 0 && mi == 0) cout << "Fg: " << Fg << "; Fe: " << Fe << "; m: " << m << "; Gp: " << temp << endl;
	if(abs(GSL_REAL(retVal)) > prec or abs(GSL_IMAG(retVal)) > prec) return retVal;
	else return gsl_complex_rect(0.0, 0.0);
	}
//
	//
// Funkcijas, kas reekjina Gamma_p koeficientus Probe
//
double splitProbe(int Fe, int M, int Fg, int m, int step){
    return energySplitsMatrixProbe[1][Fe - FStartEx_probe][M + FEndEx_probe][step] - energySplitsMatrixProbe[0][Fg - FStartGr_probe][m + FEndGr_probe][step]; // normaala izskanjoshanaas
    // return energySplitsMatrix[1][Fe - FStartEx][M + FEndEx][0] - energySplitsMatrix[0][Fg - FStartGr][m + FEndGr][0]; // izskanjoshanaas efekti "izsleegti", visur izmantojam nulles lauka izskanjosahnos
}
gsl_complex GammaPCalcPlusProbe(int Fe, int m, int Fg, int mi, int step){
	gsl_complex retVal;
	dcomplex temp;

	temp = 1.0 / (Gamma_mazs_probe + (Gamma_n_probe/2) + (Gamma_exc) + (liin_platums_probe/2) + dc_imagunit * (detuning_probe - splitProbe(Fe, m, Fg, mi, step)));
	retVal = gsl_complex_rect(real(temp), imag(temp));
	//if(Fg == 0) cout << "Fg: " << Fg << "; Fe: " << Fe << "; m: " << m << "; Gp(sauc): " << GSL_REAL(sauc) << " + " << GSL_IMAG(sauc) << endl;
	if(abs(GSL_REAL(retVal)) > prec or abs(GSL_IMAG(retVal)) > prec) return retVal;
	else return gsl_complex_rect(0.0, 0.0);
	}
gsl_complex GammaPCalcMinusProbe(int Fe, int m, int Fg, int mi, int step){
	gsl_complex retVal;
	dcomplex temp;

	//cout << Gamma_mazs << endl;
	temp = 1.0 / (Gamma_mazs_probe + (Gamma_n_probe/2) + (Gamma_exc) + (liin_platums_probe/2) - dc_imagunit * (detuning_probe - splitProbe(Fe, m, Fg, mi, step)));
	retVal = gsl_complex_rect(real(temp), imag(temp));
	//if(m == 0 && mi == 0) cout << "Fg: " << Fg << "; Fe: " << Fe << "; m: " << m << "; Gp: " << temp << endl;
	//if(Fg == 2 && Fe == 2 && m == 0 && mi == 0) cout << "Fg: " << Fg << "; Fe: " << Fe << "; m: " << m << "; Gp: " << temp << endl;
	if(abs(GSL_REAL(retVal)) > prec or abs(GSL_IMAG(retVal)) > prec) return retVal;
	else return gsl_complex_rect(0.0, 0.0);
	}
gsl_complex GammaPCalcProbe(int Fe, int m, int Fg, int mi, int step){
	gsl_complex retVal;
	dcomplex temp;
	//cout << Gamma_mazs << endl;
	// temp = (Gamma_mazs + (Gamma_n/2) + (Gamma_exc) + (liin_platums_probe/2)) / (pow(Gamma_mazs + (Gamma_n/2) + (Gamma_exc) + (liin_platums/2),2) + pow(detuning_probe - splitProbe(Fe, m, Fg, mi, step),2));
	// double gamma_temp = 0.001;
	// temp = 1.0 / (pow(gamma_temp,2) + pow(detuning_probe - splitProbe(Fe, m, Fg, mi, step),2));
	double gamma_temp = Gamma_mazs_probe + (Gamma_n_probe/2) + (Gamma_exc) + (liin_platums_probe/2);
	temp = gamma_temp / (pow(gamma_temp,2) + pow(detuning_probe - splitProbe(Fe, m, Fg, mi, step),2));
	// temp = 1.0 / ((0.1,2) + pow(detuning_probe - splitProbe(Fe, m, Fg, mi, step),2));
	// cout << "Fe = " << Fe << "; M = " << m << "; Fg = " << Fg << "; m = " << mi << "; detuning_probe - splitProbe(Fe, m, Fg, mi, step) " << detuning_probe - splitProbe(Fe, m, Fg, mi, step) << endl;
	retVal = gsl_complex_rect(real(temp), imag(temp));
	//if(m == 0 && mi == 0) cout << "Fg: " << Fg << "; Fe: " << Fe << "; m: " << m << "; Gp: " << temp << endl;
	//if(Fg == 2 && Fe == 2 && m == 0 && mi == 0) cout << "Fg: " << Fg << "; Fe: " << Fe << "; m: " << m << "; Gp: " << temp << endl;
	if(abs(GSL_REAL(retVal)) > prec or abs(GSL_IMAG(retVal)) > prec) return retVal;
	else return gsl_complex_rect(0.0, 0.0);
	}
//
// sadursmju relatiivie aatrumi
double col_rates(int F, int m, int Fp, int mp, int step, double Delta_col = Delta_col){
	return pow(Delta_col, 2) / (pow(Delta_col, 2) + pow((energySplitsMatrix[1][F - FStartEx][m + FEndEx][step] - energySplitsMatrix[1][Fp - FStartEx][mp + FEndEx][step]), 2));
}
//
// sadursmju relaksaacija
/*double col_relax(int F, int m, step){
	for(int Fk = FStartEx; Fk
}*/
//
// sareekjinam vietas
void vietas(){
	int vieta = 0;
	for(int Fg = FStartGr; Fg <= FEndGr; Fg++){
		for(int m = -Fg; m <= Fg; m++){
			for(int Fgp = FStartGr; Fgp <= FEndGr; Fgp++){
				for(int mp = -Fgp; mp <= Fgp; mp++){
					vietaGr[Fg][m][Fgp][mp] = vieta;
					vieta++;
				}
			}
		}
	}
	for(int Fe = FStartEx; Fe <= FEndEx; Fe++){
		for(int M = -Fe; M <= Fe; M++){
			for(int Fep = FStartEx; Fep <= FEndEx; Fep++){
				for(int Mp = -Fep; Mp <= Fep; Mp++){
					vietaEx[Fe][M][Fep][Mp] = vieta;
					vieta++;
				}
			}
		}
	}

	vieta = 0;
	for(int Fg = FStartGr; Fg <= FEndGr; Fg++){
		for(int m = -Fg; m <= Fg; m++){
			vietaMxGr[Fg][m] = vieta;
			vieta++;
		}
	}
	vieta = 0;
	for(int Fe = FStartEx; Fe <= FEndEx; Fe++){
		for(int M = -Fe; M <= Fe; M++){
			vietaMxEx[Fe][M] = vieta;
			vieta++;
		}
	}
}
//
// Landē faktors supersīkstruktūras līmenim
//
double p1(double w){
	return w * (w + 1);
	}
double gF(int F, int level){
	double gf;
	gf = gj[level] * ((p1(F) + p1(JJ[level]) - p1(ik)) / (2 * p1(F))) + gi * ((p1(F) + p1(ik) - p1(JJ[level])) / (2 * p1(F)));
	//if(level == 0) return gf;
	//else return 0.0;
	return gf;
	}
//
// saliekam matricaa rate equation sisteemu
//
void buildLVSAlternative(int step){
	// mainiigo deklaraacijas
	gsl_complex temp, tempElement;
	gsl_complex Gp;
	gsl_complex indGr = gsl_complex_rect(0, 0), indEx = gsl_complex_rect(0, 0);
	double tempR, tempI;
	int i = 0, j = 0;
	//
	// visiem elementiem piešķiram nulles veertiibas
	for(int ii = 0; ii < dimNum; ii++){
		for(int jj = 0; jj < dimNum; jj++){
			gsl_matrix_complex_set(ppMatrix, ii, jj, gsl_complex_rect(0.0,0.0));
			}
		gsl_vector_complex_set(ppVector, ii, gsl_complex_rect(0.0,0.0));
		}
	//cout << "dimNum: " << dimNum << endl;
	// vienkaarsi ejam cauri visai matricai
	// saakam ar pamatliimeniem
	for(int Fg = FStartGr; Fg <= FEndGr; Fg++){
		for(int m = -Fg; m <= Fg; m++){
			for(int Fgp = FStartGr; Fgp <= FEndGr; Fgp++){
				for(int mp = -Fgp; mp <= Fgp; mp++){
					// pirmais elements - induceetaas paarejas uz leju
					temp = gsl_complex_rect(0, 0);
					for(int Fe = FStartEx; Fe <= FEndEx; Fe++){
						for(int M = -Fe; M <= Fe; M++){
							for(int Fep = FStartEx; Fep <= FEndEx; Fep++){
								for(int Mp = -Fep; Mp <= Fep; Mp++){
									Gp = gsl_complex_add(GammaPCalcPlus(Fep, Mp, Fg, m, step), GammaPCalcMinus(Fe, M, Fgp, mp, step));
									temp = gsl_complex_mul(Gp, gsl_complex_mul(baseElementSimpleInd(Fg, Fe, m, M, step), baseElementSimple(Fep, Fgp, Mp, mp, step)));
									j = vietaEx[Fe][M][Fep][Mp];
									tempElement = gsl_complex_add(temp, gsl_matrix_complex_get(ppMatrix, i, j));
									if(abs(GSL_REAL(temp)) > prec or abs(GSL_IMAG(temp)) > prec){
										gsl_matrix_complex_set(ppMatrix, i, j, tempElement);
										//if(abs(GSL_REAL(temp)) > prec) cout << "i: " << i << "; j: " << j << "; val: (" << GSL_REAL(temp) << "," << GSL_IMAG(temp) << ")" << endl;
									}
									}
								}
							}
						}
					// otrais elements - ierosme uz augshu
					for(int Fgpp = FStartGr; Fgpp <= FEndGr; Fgpp++){
						for(int mpp = -Fgpp; mpp <= Fgpp; mpp++){
							temp = gsl_complex_rect(0,0);
							tempI = 0.0;
							tempR = 0.0;
							for(int Fe = FStartEx; Fe <= FEndEx; Fe++){
								for(int M = -Fe; M <= Fe; M++){
									Gp = GammaPCalcMinus(Fe, M, Fgp, mp, step);
									temp = gsl_complex_mul(Gp, gsl_complex_mul(baseElementSimpleInd(Fg, Fe, m, M, step), baseElementSimple(Fe, Fgpp, M, mpp, step)));
									tempR += GSL_REAL(temp);
									tempI += GSL_IMAG(temp);
									}
								}
							temp = gsl_complex_rect(tempR, tempI);
							j = vietaGr[Fgpp][mpp][Fgp][mp];
							if(abs(GSL_REAL(temp)) > prec or abs(GSL_IMAG(temp)) > prec){
								temp = gsl_complex_sub(gsl_matrix_complex_get(ppMatrix, i, j), temp);
								gsl_matrix_complex_set(ppMatrix, i, j, temp);
								}
							}
						}
					// treshais elements - ierosme uz augshu
					for(int Fgpp = FStartGr; Fgpp <= FEndGr; Fgpp++){
						for(int mpp = -Fgpp; mpp <= Fgpp; mpp++){
							temp = gsl_complex_rect(0.0,0.0);
							tempI = 0.0;
							tempR = 0.0;
							for(int Fe = FStartEx; Fe <= FEndEx; Fe++){
								for(int M = -Fe; M <= Fe; M++){
									Gp = GammaPCalcPlus(Fe, M, Fg, m, step);
									temp = gsl_complex_mul(Gp,gsl_complex_mul(baseElementSimpleInd(Fgpp, Fe, mpp, M, step), baseElementSimple(Fe, Fgp, M, mp, step)));
									tempR += GSL_REAL(temp);
									tempI += GSL_IMAG(temp);
									}
								}
							temp = gsl_complex_rect(tempR, tempI);
							j = vietaGr[Fg][m][Fgpp][mpp];
							//cout << "mi = " << m << "; mip = " << mp << "; mipp = " << mpp << "; i -> " << i <<"; j -> " << j << "; pp -> " << GSL_REAL(temp) << " + " << GSL_IMAG(temp) << " * i " << endl;
							if(abs(GSL_REAL(temp)) > prec or abs(GSL_IMAG(temp)) > prec){
								temp = gsl_complex_sub(gsl_matrix_complex_get(ppMatrix, i, j), temp);
								gsl_matrix_complex_set(ppMatrix, i, j, temp);
								}
							}
						}
					// ceturtais loceklis - saskjelshanaas !!!
					temp = gsl_complex_mul_real(imagunit, energySplitsMatrix[0][Fg - FStartGr][m + FEndGr][step] - energySplitsMatrix[0][Fgp - FStartGr][mp + FEndGr][step]);
					temp = gsl_complex_sub(gsl_matrix_complex_get(ppMatrix, i, i), temp);
					// if(Fg != Fgp) gsl_matrix_complex_set(ppMatrix, i, i, temp); //izslegti pamatstavokla koherentie efekti (2013)
					gsl_matrix_complex_set(ppMatrix, i, i, temp); // pamatstavokla koherences ieslegtas
					// piektais loceklis - spontaanaas paarejas no augshas !!!
					for(int Fe = FStartEx; Fe <= FEndEx; Fe++){
						for(int M = -Fe; M <= Fe; M++){
							for(int Fep = FStartEx; Fep <= FEndEx; Fep++){
								for(int Mp = -Fep; Mp <= Fep; Mp++){
									j = vietaEx[Fe][M][Fep][Mp];
									temp = gsl_complex_rect(relax(Fe, Fep, Fg, Fgp, M, Mp, m, mp, step),0);
									//if(Fe == Fep && Fe == 2 && M == Mp && M == 1 && Fg == Fgp && Fg == 2 && m == mp && m == 2) cout << "Fg = " << Fg << "; Fgp = " << Fgp << "; m = " << m << "; mp = " << mp << "; Fe = " << Fe << "; Fep = " << Fep << "; M = " << M << "; Mp = " << Mp << "; i -> " << i <<"; j -> " << j << "; pp -> " << GSL_REAL(temp) << " + " << GSL_IMAG(temp) << " * i " << endl;
									if(abs(GSL_REAL(temp)) > prec or abs(GSL_IMAG(temp)) > prec){
										temp = gsl_complex_add(temp, gsl_matrix_complex_get(ppMatrix, i, j));
										gsl_matrix_complex_set(ppMatrix, i, j, temp);
										}
									}
								}
							}
						}
					// sestais loceklis - relaksaacija izlidojot no laazera stara
					temp = gsl_complex_rect(-Gamma_mazs,0); // Normaala caurlidoshanas relaksaacija
					// temp = gsl_complex_rect(m == mp && Fg == Fgp ? -Gamma_mazs : -1e9,0); // "izsleegti" pamatstaavokla efekti, visas pamatstaavokla koherences relaksee loti aatri
					temp = gsl_complex_add(temp, gsl_matrix_complex_get(ppMatrix, i, i));
					gsl_matrix_complex_set(ppMatrix, i, i, temp);
					// septiitais loceklis - ielidoshanas repopulaacija - paarnests uz labo pusi
					if(m == mp && Fg == Fgp){
						temp = gsl_complex_rect(-Gamma_mazs,0);
						gsl_vector_complex_set(ppVector, i, temp);
						//cout << "Fg -> " << Fg << "; m -> " << m << "; Fgp -> " << Fgp << "; mp -> " << m << " lambda: " << GSL_REAL(temp) << " + " << GSL_IMAG(temp) << " * i" << endl;
						}
					i++;
					}
				}
			}
		}
	// tagad kjeramies pie ierosinaataa liimena
	for(int Fe = FStartEx; Fe <= FEndEx; Fe++){
		for(int M = -Fe; M <= Fe; M++){
			for(int Fep = FStartEx; Fep <= FEndEx; Fep++){
				for(int Mp = -Fep; Mp <= Fep; Mp++){
					// pirmais loceklis - ierosme no pamatstaavokla
					for(int Fg = FStartGr; Fg <= FEndGr; Fg++){
						for(int m = -Fg; m <= Fg; m++){
							for(int Fgp = FStartGr; Fgp <= FEndGr; Fgp++){
								for(int mp = -Fgp; mp <= Fgp; mp++){
									j = vietaGr[Fg][m][Fgp][mp];
									temp = gsl_complex_mul(baseElementSimple(Fe, Fg, M, m, step), baseElementSimpleInd(Fgp, Fep, mp, Mp, step));
									Gp = gsl_complex_add(GammaPCalcMinus(Fe, M, Fgp, mp, step), GammaPCalcPlus(Fep, Mp, Fg, m, step));
									temp = gsl_complex_mul(temp, Gp);
									//cout << "mi = " << mi << "; mip = " << mip << "; M = " << M << "; Mp = " << Mp << "; i -> " << i <<"; j -> " << j << "; pp -> " << GSL_REAL(temp) << " + " << GSL_IMAG(temp) << " * i " << endl;
									temp = gsl_complex_add(temp, gsl_matrix_complex_get(ppMatrix, i, j));
									if(abs(GSL_REAL(temp)) > prec or abs(GSL_IMAG(temp)) > prec) gsl_matrix_complex_set(ppMatrix, i, j, temp);
									}
								}
							}
						}
					// otrais loceklis - induceetaas paarejas uz leju
					for(int Fepp = FStartEx; Fepp <= FEndEx; Fepp++){
						for(int Mpp = -Fepp; Mpp <= Fepp; Mpp++){
							j = vietaEx[Fepp][Mpp][Fep][Mp];
							//cout << "Fepp nr: " << Fepp - FStartEx << "; Mpp nr: " << Mpp + FEndEx << "; Fep nr: " << Fep - FStartEx << "; Mp nr: " << Mp + FEndEx << "; j: " << j << endl;
							temp = gsl_complex_rect(0.0, 0.0);
							tempI = 0.0;
							tempR = 0.0;
							for(int Fg = FStartGr; Fg <= FEndGr; Fg++){
								for(int m = -Fg; m <= Fg; m++){
									Gp = GammaPCalcPlus(Fep, Mp, Fg, m, step);
									temp = gsl_complex_mul(Gp, gsl_complex_mul(baseElementSimple(Fe, Fg, M, m, step), baseElementSimpleInd(Fg, Fepp, m, Mpp, step)));
									tempR += GSL_REAL(temp);
									tempI += GSL_IMAG(temp);
									}
								}
							temp = gsl_complex_rect(tempR, tempI);
							//cout << "M = " << M << "; Mp = " << Mp << "; Mpp = " << Mpp << "; i -> " << i <<"; j -> " << j << "; Gp -> " << GSL_REAL(Gp) << " + " << GSL_IMAG(Gp) << " * i ; pp -> " << GSL_REAL(temp) << " + " << GSL_IMAG(temp) << " * i " << endl;
							if(abs(GSL_REAL(temp)) > prec or abs(GSL_IMAG(temp)) > prec){
								temp = gsl_complex_sub(gsl_matrix_complex_get(ppMatrix, i, j), temp);
								gsl_matrix_complex_set(ppMatrix, i, j, temp);
								}
							}
						}
					// treshais loceklis - induceetaas paarejas uz leju
					for(int Fepp = FStartEx; Fepp <= FEndEx; Fepp++){
						for(int Mpp = -Fepp; Mpp <= Fepp; Mpp++){
							j = vietaEx[Fe][M][Fepp][Mpp];
							temp = gsl_complex_rect(0.0, 0.0);
							tempI = 0.0;
							tempR = 0.0;
							for(int Fg = FStartGr; Fg <= FEndGr; Fg++){
								for(int m = -Fg; m <= Fg; m++){
									Gp = GammaPCalcMinus(Fe, M, Fg, m, step);
									temp = gsl_complex_mul(Gp, gsl_complex_mul(baseElementSimple(Fepp, Fg, Mpp, m, step), baseElementSimpleInd(Fg, Fep, m, Mp, step)));
									tempR += GSL_REAL(temp);
									tempI += GSL_IMAG(temp);
									//Gp = gsl_complex_mul_real(Gp, Rabi * Rabi);
									}
								}
							temp = gsl_complex_rect(tempR, tempI);
							if(abs(GSL_REAL(temp)) > prec or abs(GSL_IMAG(temp)) > prec){
								temp = gsl_complex_sub(gsl_matrix_complex_get(ppMatrix, i, j), temp);
								gsl_matrix_complex_set(ppMatrix, i, j, temp);
								}
							}
						}
					// ceturtais loceklis - aareejais lauks !!!
					temp = gsl_complex_mul_real(imagunit, energySplitsMatrix[1][Fe - FStartEx][M + FEndEx][step] - energySplitsMatrix[1][Fep - FStartEx][Mp + FEndEx][step]);
					temp = gsl_complex_sub(gsl_matrix_complex_get(ppMatrix, i, i), temp);
					// if(Fe != Fep) gsl_matrix_complex_set(ppMatrix, i, i, temp); // OFF
					gsl_matrix_complex_set(ppMatrix, i, i, temp); // ON
					//gsl_matrix_complex_set(ppMatrix, i, i, temp);
					// piektais loceklis - spontaanaa sabruksana
					temp = gsl_complex_add_real(gsl_matrix_complex_get(ppMatrix, i, i), -Gamma_n - Gamma_exc);
					gsl_matrix_complex_set(ppMatrix, i, i, temp);
					// piektais loceklis - izlidoshana
					temp = gsl_complex_add_real(gsl_matrix_complex_get(ppMatrix, i, i), -Gamma_mazs); // normaala caurlidoshanas relaksaacija
					// temp = gsl_complex_add_real(gsl_matrix_complex_get(ppMatrix, i, i), Fe == Fep && M == Mp ? -Gamma_mazs : -1e9); // "izsleegti" ierosinaataa staavokla efekti - ierosinaataa staavokla kohereneces relaksee loti aatri
					gsl_matrix_complex_set(ppMatrix, i, i, temp);
					// nekompensētais lauks
					/*if(Fe == Fep){
    					for(int Mpp = max(M-1, -Fe); Mpp <= min(M+1, Fe); Mpp++){
    					    j = vietaEx[Fe][M][Fe][Mpp];
    					    temp = gsl_complex_add(gsl_matrix_complex_get(ppMatrix, i, j), gsl_complex_rect(gF(1, Fe) * Btrans * mi_Bora_MHz,0));
        					gsl_matrix_complex_set(ppMatrix, i, j, temp);
    					}
    					for(int Mpp = max(Mp-1, -Fe); Mpp <= min(Mp+1, Fe); Mpp++){
    					    j = vietaEx[Fe][Mp][Fe][Mpp];
    					    temp = gsl_complex_sub(gsl_matrix_complex_get(ppMatrix, i, j), gsl_complex_rect(gF(1, Fe) * Btrans * mi_Bora_MHz,0));
        					gsl_matrix_complex_set(ppMatrix, i, j, temp);
    					}
					}*/
					//cout << i << endl;
//					if(Fe == Fep && M == Mp){
//					    tempI = 0;
//    					for(int Fepp = FStartEx; Fepp <= FEndEx; Fepp++){
//    						for(int Mpp = -Fepp; Mpp <= Fepp; Mpp++){
//    						    tempR = 0;
//            					for(int Feppp = FStartEx; Feppp <= FEndEx; Feppp++){
//            						for(int Mppp = -Feppp; Mppp <= Feppp; Mppp++){
//            						    tempR += col_rates(Fepp, Mpp, Feppp, Mppp, step);
//    						        }
//    						    }
//    						    j = vietaEx[Fepp][Mpp][Fepp][Mpp];
//    						    temp = gsl_complex_add_real(gsl_matrix_complex_get(ppMatrix, i, j), Gamma_exc * col_rates(Fepp, Mpp, Fe, M, step)/tempR);
//					        }
//				        }
//					}
					i++;
					}
				}
			}
		}
	//cout << "pabeidzu vienaadojumu matricu ... " << endl;
	}
//
// pati galvenaa reekjinu funkcija
//
int BSign = 1;
ofstream dataFileEx, dataFileGr, dataFilePr, matrixFile, momFile, matrixFileGr, matrixFileEx, dataFileExProbe, dataFileGrProbe;


#include "energija_jacobi.cc"
#include "utilities.cc"
int main(){
	init_u();
    //momFile.open("pol/moments.dat", ios::out | ios::trunc);
    //momFile << "B Fe0 Fe1 Fe2 Fe3 FeS Esum Fg1 Fg2 FgS Gsum" << endl;
	ofstream dataFile, dataFileProbe;
	cout << "Fg = " << FStartGr << " -- " << FEndGr << "; " << "Fe = " << FStartEx << " -- " << FEndEx << "; dim = " << dimNum << endl;
	cout << "Probe Fg = " << FStartGr_probe << " -- " << FEndGr_probe << "; " << "Fe = " << FStartEx_probe << " -- " << FEndEx_probe << "; dim = " << dimNum_probe << endl;
	//const char base[] = "data/D2_87/g00145_D125_L10/2_3/";
	//dataFileEx.open("magMixRb87Ex.txt", ios::out | ios::trunc);
	//dataFileGr.open("magMixRb87Gr.txt", ios::out | ios::trunc);
	//dataFileEx.open("magSplitsRb87Ex.txt", ios::out | ios::trunc);
	//dataFileGr.open("magSplitsRb87Gr.txt", ios::out | ios::trunc);
	//dataFilePr.open("data/D2_85/levels/probsRb85_pi_G3.txt", ios::out | ios::trunc);
	//dataFile.open("dataInt2.txt", ios::out | ios::trunc);

    char filename [ FILENAME_MAX ];
    char fldname [ FILENAME_MAX ];
    char fldname2 [ FILENAME_MAX ];
    char dmfldname [ FILENAME_MAX ];

	double diff, der1, der1a, der1b, intensity_par = 0, intensity_per = 0, intensity_z = 0, intensity_plus = 0, intensity_minus = 0, polarity, polarity_abs, tRabi, sinus, cosinus, d_eff, bright = 0, dark = 0, bark = 0, right = 0, sight = 0, absorption_par = 0, absorption_per = 0;
	int dim, FMU, shift;
	int useI, sign_status;
	int usedIs[10];
	int ev_signs[2][10][10];
	int Bi = 0;
	double B, BB, Bpar, DShift, detBase, lShift, detBaseProbe, lShift_probe;
	int indexOfMin;
	double apdzGrCnt = 0, apdzExCnt = 0;
	// 3_3 double Rabi_matrix[11] = {0.866, 1.194, 1.936, 2.739, 3.873, 5.477, 7.5, 10.607, 15, 21.213, 28.723};
		/*double Rabi_matrix[25] = {5, 20
        };*/

	double Rabi_matrix[30] =  { 0.0,
		0.5,	// 1
		1.0,	// 2 (shie ir indeksi ar ko izsauc vajadzigos Rabi)
		5.0,	// 3
		10.0,	// 4
		20.0,	// 5
		50.0,	// 6
		100.0,	// 7
		60.0,	// 8
		80.0,	// 9
		150.0,	// 10
		200.0,	// 11
		40.0,	// 12
		30.0,	// 13
		17*0.67,	// 14
		18*0.67,	// 15
		19*0.67,	// 16
		32,	// 17
		45.25,	// 18
		53.82,	// 19
		64,	// 20
		76.11,	// 21
		90.51,	// 22
		100.,	// 23
		Gamma_n / 100, // 24
		Gamma_n / 10, // 25
		Gamma_n, // 26
		Gamma_n *10 // 27
	};

	//double Rabi_matrix[25] = { 0.447, 0.379, 0.7086011328, 1.0021133322, 1.2273331642, 1.4172022655, 1.5844803017, 1.7357112063, 1.874782376 };
	//double Rabi_matrix[20] = { 6.45, 11.2, 25, 50, 100, 2, 2.83, 4, 5.66, 8, 11.31, 16, 22.63, 32, 45.25, 64, 90.51, 128 };
	//double Rabi_matrix[15] = { 30.17, 25, 30, 40, 50, 60 };
	//double Rabi_matrix[15] = { 6.45, 11.18, 15.81, 25, 35.36, 50 };
	//double Rabi_matrix[15] = { 1.5 };
	double lWidth_matrix[10] = { 0.1, 0.2, 0.5, 1, 2, 5, 10, 20 };
	double gamma_matrix[12] = { 0.0155, 0.019, 0.0095, 0.014, 0.124, 0.186, 0.228, Gamma_n * 0.25, Gamma_n * 0.5, Gamma_n, Gamma_n * 2, Gamma_n / 100 };
	double gamma_exc_matrix[12] = { 0, 1 * Gamma_n, 2 * Gamma_n, 5 * Gamma_n, 6 * Gamma_n, 10 * Gamma_n, 30 * Gamma_n };
	double lambda_matrix[12] = {1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 4};
	//double det_matrix[16] = { 0, 1.15, 2.3, 3.45, 4.6, 5.75, 6.9, 8.05, 9.2, 10.35, 11.5, 12.65, 13.8, 14.95, 16.1, 18 };
	double mult_matrix[13] = { 0.1, 0.2, 0.4, 0.6, 0.8, 1.25, 1.5, 2, 3, 4, 6, 8, 10 };
	double Bpar_matrix[10] = {-0.1, 0, -0.05, -0.025, 0.025, 0.05, 0.1};
	double Temp_matrix[10] = {293, 298, 303, 308, 313, 323};
	double pres_matrix[10] = {2.2e-7, 4e-7, 6.5e-7, 1.2e-6, 2e-6, 5e-6};
	double dFact = 0.0;
	double det_matrix[16] = {0};
	//gsl_permutation * p = gsl_permutation_alloc (dimNum);

	gsl_complex tmp;
	vietas();
	build_three_j_matrix();
	build_three_j_matrix_probe();
	build_reduced_matrix();
	fill_rot_matrices();
	cout << "ierosme: pol -> " << pol << "; teta -> " << teta_i << "; fii -> " << fii_i << endl;
	cout << "E+1: (" << GSL_REAL(E(1, pol, teta_i, fii_i)) << "," << GSL_IMAG(E(1, pol, teta_i, fii_i)) << ")"
		<< "; E 0: (" << GSL_REAL(E(0, pol, teta_i, fii_i)) << "," << GSL_IMAG(E(0, pol, teta_i, fii_i)) << ")"
		<< "; E-1: (" << GSL_REAL(E(-1, pol, teta_i, fii_i)) << "," << GSL_IMAG(E(-1, pol, teta_i, fii_i)) << ")"
		<< endl;
	cout << "noveeroshana: pol -> " << pol_o[0] << "; teta -> " << teta_o[0] << "; fii -> " << fii_o[0] << endl;
	cout << "E+1: (" << GSL_REAL(E(1, pol_o[0], teta_o[0], fii_o[0])) << "," << GSL_IMAG(E(1, pol_o[0], teta_o[0], fii_o[0])) << ")"
		<< "; E 0: (" << GSL_REAL(E(0, pol_o[0], teta_o[0], fii_o[0])) << "," << GSL_IMAG(E(0, pol_o[0], teta_o[0], fii_o[0])) << ")"
		<< "; E-1: (" << GSL_REAL(E(-1, pol_o[0], teta_o[0], fii_o[0])) << "," << GSL_IMAG(E(-1, pol_o[0], teta_o[0], fii_o[0])) << ")"
		<< endl;
	cout << "E+1: (" << GSL_REAL(E(1, pol_o[1], teta_o[1], fii_o[1])) << "," << GSL_IMAG(E(1, pol_o[1], teta_o[1], fii_o[1])) << ")"
		<< "; E 0: (" << GSL_REAL(E(0, pol_o[1], teta_o[1], fii_o[1])) << "," << GSL_IMAG(E(0, pol_o[1], teta_o[1], fii_o[1])) << ")"
		<< "; E-1: (" << GSL_REAL(E(-1, pol_o[1], teta_o[1], fii_o[1])) << "," << GSL_IMAG(E(-1, pol_o[1], teta_o[1], fii_o[1])) << ")"
		<< endl;
	//for(double atrums = 0.0; atrums <= 0.0; atrums += 10){
	B = 0;
	/*for(detuning = -DRange; detuning <= DRange; detuning += 50){
    sprintf(filename, "%s%d", base, (int)detuning);
	dataFile.open(filename, ios::out | ios::trunc);*/
	/*dataFile.open("data/87/pop/1_1.dat", ios::out | ios::trunc);
	dataFile << "Rabi ";
	for(int Fg = FStartGr; Fg <= FEndGr; Fg++){
		for(int m = -Fg; m <= Fg; m++){
			for(int Fgp = FStartGr; Fgp <= FEndGr; Fgp++){
				for(int mp = -Fgp; mp <= Fgp; mp++){
					if(Fg == Fgp && m == mp) dataFile << "F" << Fg << "" << m << " ";
				}
			}
		}
	}
	dataFile << "sum" << endl;*/
	double Sigma_base = Sigma;
	// for(int gamma_i = 1; gamma_i < 2; gamma_i++){
	double probe_det_step = 50;
	// for(double probe_det = -1000+9200; probe_det <= 1000+9200; probe_det += probe_det_step){
	for(double probe_det = 0; probe_det <= 0; probe_det += probe_det_step){
	// for(double dshift1 = -100; dshift1 <= 100; dshift1 += probe_det_step){

	//for(lShift = -100; lShift <= -50; lShift += 2.5){
	double lambda_k = lambda_matrix[0];
	// Gamma_mazs = gamma_matrix[gamma_i];
	Gamma_mazs = gamma_matrix[1]; //gamma_matrix[1];
	Gamma_mazs_probe = gamma_matrix[1];
	//liin_platums = lwidth_matrix[gamma_i];
	double Gamma_mazs_base = Gamma_mazs;
	//Sigma = Sigma_base * (Gamma_n / (2 * Gamma_mazs));
	//if(Sigma > 210.14) Sigma = 210.14;
	//DScan = 2 * Sigma;
	//lShift = det_matrix[gamma_i];
	lShift = 0; //-5316; //-1919; // 0; // -2280; //0; // //Izskanosana no parejas, vienibas - MHz
	lShift_probe = 0; // 24.9; // crossover(3;5) nprobe_det;//0; //-5316; //-1919; // 0; // -2280; //0; // //Izskanosana no parejas, vienibas - MHz
	for(int ge_i = 0; ge_i < 1; ge_i++){
	//Gamma_exc = gamma_exc_matrix[ge_i] * exp ((-1) * (3 - Gamma_mazs) / Gamma_mazs);
	Gamma_exc = gamma_exc_matrix[ge_i];
	Delta_col = Gamma_exc / 4;
	liin_platums = 2; //lWidth_matrix[4]; //
	liin_platums_probe = 2; //2; //
	// liin_platums = 10;
	Bpar = 0;
	Btrans = 0;

	// maticas datu saglabaashanas fails
	/*sprintf(fldname, "%s%s%.4f%s%i%s%.1f", baseFolder, "gamma=", Gamma_mazs, "-DSteps=", DSteps, "-lWidth=", liin_platums);
	if(!FileExists(fldname)){ cout << "creating folder " << fldname << endl; mkdir(fldname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); }
	sprintf(fldname, "%s%s%i%s%i%s%i%s%.0f", fldname, "/", Fg_d, "-", Fe_d, "-pol=", pol, "-det=", lShift);
	if(!FileExists(fldname)){ cout << "creating folder " << fldname << endl; mkdir(fldname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); }
	sprintf(filename, "%s%s%s%s", fldname, "/", baseFilename, ".dat");
	cout << "filename: " << filename << "; sigma_base: " << Sigma_base << endl;
	matrixFile.open(filename, ios::out | ios::trunc);*/

	for(int Ri = 1; Ri <= 7; Ri += 1){ //sheit izvelas un ievieto vajdzigos rabi indeksus, kas atbilsti rabi frekvencem kas definetas augstak/ Taisa rabi ciklus sakas ar Ri= un beigu vertiba Ri< un Ri+= saka par cik palielinat katru nakamo
	//if(ge_i == 0 && Ri == 8) Ri = 17;
	Rabi = Rabi_matrix[Ri];
	//if(Rabi > 16) BRange = 15;
	//if(Rabi > 32) BRange = 20;
	tRabi = Rabi * 1000;
    //sprintf(filename, "%s%f%s%f", base, Gamma_mazs_base, "/Rb85_D1_22_Rabi=", Rabi);
	double konc = pressure / (k_Boltz * Temp); // m^-3 koncentraacija atmosfeeras spiedienaa
	double Gamma_sad = 1e-6 * konc * atrums * sigma_sad / (2 * Pi);

	// if(!FileExists(baseFolder)){ cout << "creating folder " << baseFolder << endl; mkdir(baseFolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); }
	sprintf(fldname, "%s%s", baseFolder, baseFolderProbe);
	if(!FileExists(fldname)){ cout << "creating folder " << fldname << endl; mkdir(fldname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); }
	// sprintf(fldname, "%s%s%.4f%s%i%s%.2f%s%s%.1f", baseFolder, "gamma=", Gamma_mazs_base, "-DSteps=", DSteps, "-DScan=", DScan/Sigma, "sigma", "-lWidth=", liin_platums);
	sprintf(fldname, "%s%s%.4f%s%i%s%.2f%s%s%.1f", fldname, "gamma=", Gamma_mazs_base, "-DSteps=", DSteps, "-DScan=", DScan/Sigma, "sigma", "-lWidth=", liin_platums);
	if(!FileExists(fldname)){ cout << "creating folder " << fldname << endl; mkdir(fldname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); }
	sprintf(fldname, "%s%s%i%s%i%s%.2f%s%.2f%s%i%s%.0f%s%i%s%i%s%.0f", fldname, "/Pump_", Fg_d, "-", Fe_d, "-theta=", teta_i, "-phi=", fii_i, "-pol=", pol, "-det=", lShift, "_Probe_", Fg_p, "-", Fe_p, "-det=", lShift_probe);
	if(!FileExists(fldname)){ cout << "creating folder " << fldname << endl; mkdir(fldname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); }

	sprintf(fldname2, "%s%s%.2f%s%.2f%s%.2f%s%.2f%s%i%s%i", fldname, "/LIF_theta=", teta_o[0],"-",teta_o[1], "_phi=", fii_o[0],"-",fii_o[1],"_pol=",pol_o[0],"-",pol_o[1]);
	if(!FileExists(fldname2)){ cout << "creating folder " << fldname2 << endl; mkdir(fldname2, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); }
	//cout << lShift << endl;
	sprintf(filename, "%s%s%s%.2f%s%.1f%s%.0f", fldname2, "/", baseFilename, Rabi, "-shift=", lShift, "-Ge=", Gamma_exc);
	cout << "filename: " << filename << "; sigma_base: " << Sigma_base << endl;
	dataFile.open(filename, ios::out | ios::trunc);
	//Absorption file
	sprintf(fldname2, "%s%s%.2f%s%.2f%s%.2f%s%.2f%s%i%s%i%s%.4f%s%.1f", fldname, "/Absorption_theta=", teta_p[0],"-",teta_p[1], "_phi=", fii_p[0],"-",fii_p[1],"_pol=",pol_p[0],"-",pol_p[1],"_gammaprobe=",Gamma_mazs_probe,"_lWidthpr=", liin_platums_probe);
	if(!FileExists(fldname2)){ cout << "creating folder " << fldname2 << endl; mkdir(fldname2, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); }
	sprintf(filename, "%s%s%s%.2f%s%.1f%s%.0f", fldname2, "/", baseFilenameProbe, Rabi, "-shift=", lShift_probe, "-Ge=", Gamma_exc);
	cout << "filename absorption: " << filename << "; sigma_base: " << Sigma_base << endl;
	dataFileProbe.open(filename, ios::out | ios::trunc);
	//cout << "!!! data file not open" << endl;
	//dataFile.open("data/D2_87/Rb87_D2_det_scan_G1", ios::out | ios::trunc);
	//energySplitsMatrix building
	for(int m = -FEndEx; m <= FEndEx; m++){
		if(abs(m) > FStartEx) FMU = abs(m);
		else FMU = FStartEx;
		shift = FMU - FStartEx;
		dim = FEndEx - FMU + 1;
		//cout << dim << "; " << shift << endl;
		//dim = FEndEx - FStartEx + 1;
		energySplits(m, FStartEx, FEndEx, 1, 0);
		for (int i = 0; i < dim; i++){
			energySplitsMatrix[1][i + shift][m + FEndEx][0] = eigenValues[i];
			// cout << "m " << m << " i " << i << " "  << energySplitsMatrix[1][i + shift][m + FEndEx_probe][0] << " " << eigenValues[i] << endl;
			for(int j = 0; j < dim; j++){
				levelMixingMatrix[1][i + shift][j + shift][m + FEndEx][0] = ev_signs[1][i][j] * eigenVectors[j][i];
			}
		}
	}
	for(int m = -FEndGr; m <= FEndGr; m++){
		if(abs(m) > FStartGr) FMU = abs(m);
		else FMU = FStartGr;
		shift = FMU - FStartGr;
		dim = FEndGr - FMU + 1;
		//dim = FEndGr - FStartGr + 1;
		energySplits(m, FStartGr, FEndGr, 0, 0);
		for (int i = 0; i < dim; i++){
			energySplitsMatrix[0][i + shift][m + FEndGr][0] = eigenValues[i];
			for(int j = 0; j < dim; j++){
				levelMixingMatrix[0][i + shift][j + shift][m + FEndGr][0] = ev_signs[0][i][j] * eigenVectors[j][i];
			}
		}
	}

	//energySplitsMatrixProbe building
	for(int m = -FEndEx_probe; m <= FEndEx_probe; m++){
		if(abs(m) > FStartEx_probe) FMU = abs(m);
		else FMU = FStartEx_probe;
		shift = FMU - FStartEx_probe;
		dim = FEndEx_probe - FMU + 1;
		//cout << dim << "; " << shift << endl;
		//dim = FEndEx - FStartEx + 1;
		energySplitsProbe(m, FStartEx_probe, FEndEx_probe, 1, 0);
		// energySplits(m, FStartEx, FEndEx, 1, 0);
		for (int i = 0; i < dim; i++){
			energySplitsMatrixProbe[1][i + shift][m + FEndEx_probe][0] = eigenValuesProbe[i];
			// energySplitsMatrixProbe[1][i + shift][m + FEndEx_probe][0] = eigenValues[i];
			// cout << "m " << m << " i " << i << " " << energySplitsMatrixProbe[1][i + shift][m + FEndEx_probe][0] << " " << eigenValuesProbe[i] << " " << eigenValues[i] << endl;
			for(int j = 0; j < dim; j++){
				levelMixingMatrixProbe[1][i + shift][j + shift][m + FEndEx_probe][0] = ev_signs[1][i][j] * eigenVectorsProbe[j][i];
			}
		}
	}
	for(int m = -FEndGr_probe; m <= FEndGr_probe; m++){
		if(abs(m) > FStartGr_probe) FMU = abs(m);
		else FMU = FStartGr_probe;
		shift = FMU - FStartGr_probe;
		dim = FEndGr_probe - FMU + 1;
		//dim = FEndGr - FStartGr + 1;
		energySplitsProbe(m, FStartGr_probe, FEndGr_probe, 0, 0);
		for (int i = 0; i < dim; i++){
			energySplitsMatrixProbe[0][i + shift][m + FEndGr_probe][0] = eigenValuesProbe[i];
			for(int j = 0; j < dim; j++){
				levelMixingMatrixProbe[0][i + shift][j + shift][m + FEndGr_probe][0] = ev_signs[0][i][j] * eigenVectorsProbe[j][i];
			}
		}
	}

	Bi = 1;
	for(BSign = 1; BSign <= 1; BSign += 2){
	B = 0;
	for(int i = 0; i < 10; i++){
		for(int j = 0; j < 10; j++){
			ev_signs[0][i][j] = 1;
			ev_signs[1][i][j] = 1;
		}
	}
		//for(lShift = -450.0; lShift <= 600.0; lShift += 2.5){
	detBase = energySplitsMatrix[1][Fe_d - FStartEx][0 + FEndEx][0] - energySplitsMatrix[0][Fg_d - FStartGr][0 + FEndGr][0] + lShift;
	detBaseProbe = energySplitsMatrixProbe[1][Fe_p - FStartEx_probe][0 + FEndEx_probe][0] - energySplitsMatrixProbe[0][Fg_p - FStartGr_probe][0 + FEndGr_probe][0] + lShift_probe;
	// detBaseProbe = energySplitsMatrix[1][Fe_d - FStartEx_probe][0 + FEndEx][0] - energySplitsMatrix[0][Fg_d - FStartGr][0 + FEndGr][0] + lShift;

	BRange = 200; //3000; //
	for(BB = -200; BB <= BRange; BB += step){
	// for(BB = 1420; BB <= 1560; BB += step){
	// for(BB = 2060; BB <= 2160; BB += step){
		//lShift = BSign * BB;
		// step = 10;
		step = 5; //10;//1; //
		B = BSign * BB;
        if(abs(BB) < 50) step = 1;
		if(abs(BB) < 5) step = 0.1;
		
		// if(BB < 0.0099) step = 0.001;
		// // else if (BB < 0.00099) step = 0.00001;
		// // else if(BB < 0.099) step = 0.01; // shis ir solis kads tiek izmantos, Lauks shis vienibas ir gausos
		// else if (BB < 0.199) step = 0.1;
		// // if (BB < 4.99) step = 0.5;
		// else if (BB < 4.99) step = 0.1;
		//if (BB < 1.999) step = 0.25;
		// if (BB < 9.999) step = 0.1;
		// else step = 1;
		//step = 1;
		//step = 5;
		/*if(BB < 0.01 && abs(Bpar) > 0.01){
			teta_o[0] = Pi/2;
			teta_o[1] = 0;
			B = Bpar;
		}
		else if(abs(Bpar) > 0.01){
			teta_o[0] = Pi/2;
			teta_o[1] = Pi/2 - atan(Bpar/B);
			B = sqrt(BB*BB + Bpar*Bpar) * BSign * Bpar / abs(Bpar);
		}
		else{
			teta_o[0] = Pi/2;
			teta_o[1] = Pi/2;
		}*/
		//for(DShift = -2000; DShift <= 2000; DShift += 2){
		//if(B >= 1) step = 0.25;
		//else step = 0.1;
		intensity_par = 0; intensity_z = 0; intensity_plus = 0; intensity_minus = 0;
		intensity_per = 0; dark = 0; bark = 0; bright = 0; right = 0; sight = 0;
		absorption_par = 0; absorption_per = 0;
		time_t rawtime;
          struct tm * timeinfo;
          time (&rawtime);
          timeinfo = localtime (&rawtime);
        cout << endl << "Time: " << asctime(timeinfo) << "Rabi: " << Rabi << ";" << endl << "B -> " << B << "; step: " << step << ";" << endl << "DScan: " << DScan << "; DStep: " << DStep << "; sigma: " << Sigma <<"; Sigma_probe: " << Sigma_probe << ";" << endl << "gamma: " << Gamma_mazs << "; gamma_exc: " << Gamma_exc << "; atrums: " << atrums << "; energija: " << energija << ";" << endl;		
        for(int m = -FEndEx; m <= FEndEx; m++){
			if(abs(m) > FStartEx) FMU = abs(m);
			else FMU = FStartEx;
			shift = FMU - FStartEx;
			dim = FEndEx - FMU + 1;
			//cout << dim << "; " << shift << endl;
			//dim = FEndEx - FStartEx + 1;
			energySplits(m, FStartEx, FEndEx, 1, B);
			// energySplitsProbe(m, FStartEx_probe, FEndEx_probe, 1, B);
			for (int i = 0; i < dim; i++){
				energySplitsMatrix[1][i + shift][m + FEndEx][Bi] = eigenValues[i];
				for(int j = 0; j < dim; j++){
					levelMixingMatrix[1][i + shift][j + shift][m + FEndEx][Bi] = ev_signs[1][i][j] * eigenVectors[j][i];
				}
			}
			if(Bi > 1){
				for(int i = 0; i < dim; i++){
					for(int j = 0; j < dim; j++){
						der1 = levelMixingMatrix[1][i + shift][j + shift][m + FEndEx][Bi-2] - levelMixingMatrix[1][i + shift][j + shift][m + FEndEx][Bi-1];
						der1a = levelMixingMatrix[1][i + shift][j + shift][m + FEndEx][Bi-1] - ev_signs[1][i][j] * eigenVectors[j][i];
						der1b = levelMixingMatrix[1][i + shift][j + shift][m + FEndEx][Bi-1] - (-1) * ev_signs[1][i][j] * eigenVectors[j][i];
						if(abs(der1 - der1a) > abs(der1 - der1b)){
							ev_signs[1][i][j] *= -1;
							levelMixingMatrix[1][i + shift][j + shift][m + FEndEx][Bi] = ev_signs[1][i][j] * eigenVectors[j][i];
						}
					}
				}
			}
		}
		for(int m = -FEndGr; m <= FEndGr; m++){
			if(abs(m) > FStartGr) FMU = abs(m);
			else FMU = FStartGr;
			shift = FMU - FStartGr;
			dim = FEndGr - FMU + 1;
			//dim = FEndGr - FStartGr + 1;
			energySplits(m, FStartGr, FEndGr, 0, B);
			// energySplitsProbe(m, FStartGr_probe, FEndGr_probe, 0, B);
			for (int i = 0; i < dim; i++){
				energySplitsMatrix[0][i + shift][m + FEndGr][Bi] = eigenValues[i];
				for(int j = 0; j < dim; j++){
					levelMixingMatrix[0][i + shift][j + shift][m + FEndGr][Bi] = ev_signs[0][i][j] * eigenVectors[j][i];
				}
			}
			if(Bi > 1){
				for(int i = 0; i < dim; i++){
					for(int j = 0; j < dim; j++){
						der1 = levelMixingMatrix[0][i + shift][j + shift][m + FEndGr][Bi-2] - levelMixingMatrix[0][i + shift][j + shift][m + FEndGr][Bi-1];
						der1a = levelMixingMatrix[0][i + shift][j + shift][m + FEndGr][Bi-1] - ev_signs[0][i][j] * eigenVectors[j][i];
						der1b = levelMixingMatrix[0][i + shift][j + shift][m + FEndGr][Bi-1] - (-1) * ev_signs[0][i][j] * eigenVectors[j][i];
						if(abs(der1 - der1a) > abs(der1 - der1b)){
							ev_signs[0][i][j] *= -1;
							levelMixingMatrix[0][i + shift][j + shift][m + FEndGr][Bi] = ev_signs[0][i][j] * eigenVectors[j][i];
						}
					}
				}
			}
		}
			// energySplitsMatrixProbe
		for(int m = -FEndEx_probe; m <= FEndEx_probe; m++){
			if(abs(m) > FStartEx_probe) FMU = abs(m);
			else FMU = FStartEx_probe;
			shift = FMU - FStartEx_probe;
			dim = FEndEx_probe - FMU + 1;
			energySplitsProbe(m, FStartEx_probe, FEndEx_probe, 1, B);
			// energySplits(m, FStartEx, FEndEx, 1, B);
			for (int i = 0; i < dim; i++){
				energySplitsMatrixProbe[1][i + shift][m + FEndEx_probe][Bi] = eigenValuesProbe[i];
				// energySplitsMatrixProbe[1][i + shift][m + FEndEx_probe][Bi] = eigenValues[i];
				for(int j = 0; j < dim; j++){
					levelMixingMatrixProbe[1][i + shift][j + shift][m + FEndEx_probe][Bi] = ev_signs[1][i][j] * eigenVectorsProbe[j][i];
				}
			}
			if(Bi > 1){
				for(int i = 0; i < dim; i++){
					for(int j = 0; j < dim; j++){
						der1 = levelMixingMatrixProbe[1][i + shift][j + shift][m + FEndEx_probe][Bi-2] - levelMixingMatrixProbe[1][i + shift][j + shift][m + FEndEx_probe][Bi-1];
						der1a = levelMixingMatrixProbe[1][i + shift][j + shift][m + FEndEx_probe][Bi-1] - ev_signs[1][i][j] * eigenVectorsProbe[j][i];
						der1b = levelMixingMatrixProbe[1][i + shift][j + shift][m + FEndEx_probe][Bi-1] - (-1) * ev_signs[1][i][j] * eigenVectorsProbe[j][i];
						if(abs(der1 - der1a) > abs(der1 - der1b)){
							ev_signs[1][i][j] *= -1;
							levelMixingMatrixProbe[1][i + shift][j + shift][m + FEndEx_probe][Bi] = ev_signs[1][i][j] * eigenVectorsProbe[j][i];
						}
					}
				}
			}
		}
		for(int m = -FEndGr_probe; m <= FEndGr_probe; m++){
			if(abs(m) > FStartGr_probe) FMU = abs(m);
			else FMU = FStartGr_probe;
			shift = FMU - FStartGr_probe;
			dim = FEndGr_probe - FMU + 1;
			//dim = FEndGr - FStartGr + 1;
			energySplitsProbe(m, FStartGr_probe, FEndGr_probe, 0, B); //void energySplitsProbe(int m, int FMin, int FMax, int level, double B)
			// energySplitsProbe(m, FStartGr_probe, FEndGr_probe, 0, B);
			for (int i = 0; i < dim; i++){
				energySplitsMatrixProbe[0][i + shift][m + FEndGr_probe][Bi] = eigenValuesProbe[i];
				for(int j = 0; j < dim; j++){
					levelMixingMatrixProbe[0][i + shift][j + shift][m + FEndGr_probe][Bi] = ev_signs[0][i][j] * eigenVectorsProbe[j][i];
				}
			}
			if(Bi > 1){
				for(int i = 0; i < dim; i++){
					for(int j = 0; j < dim; j++){
						der1 = levelMixingMatrixProbe[0][i + shift][j + shift][m + FEndGr_probe][Bi-2] - levelMixingMatrixProbe[0][i + shift][j + shift][m + FEndGr_probe][Bi-1];
						der1a = levelMixingMatrixProbe[0][i + shift][j + shift][m + FEndGr_probe][Bi-1] - ev_signs[0][i][j] * eigenVectorsProbe[j][i];
						der1b = levelMixingMatrixProbe[0][i + shift][j + shift][m + FEndGr_probe][Bi-1] - (-1) * ev_signs[0][i][j] * eigenVectorsProbe[j][i];
						if(abs(der1 - der1a) > abs(der1 - der1b)){
							ev_signs[0][i][j] *= -1;
							levelMixingMatrixProbe[0][i + shift][j + shift][m + FEndGr_probe][Bi] = ev_signs[0][i][j] * eigenVectorsProbe[j][i];
						}
					}
				}
			}
		}
		// cout << "energySplitsMatrixProbe[1][Fe_p - FStartEx_probe][0 + FEndEx_probe][0] " << energySplitsMatrixProbe[1][Fe_p - FStartEx_probe][0 + FEndEx_probe][0] << " energySplitsMatrixProbe[1][Fe_p - FStartEx_probe][0 + FEndEx_probe][1] " <<energySplitsMatrixProbe[1][Fe_p - FStartEx_probe][0 + FEndEx_probe][1] <<endl;

		//detuning = 0;
		cout << "detuning -> " << detBase << "; lShift: " << lShift << "; detuning probe -> " << detBaseProbe << "; lShift probe: " << lShift_probe << endl;
		DShift = 0;
		//print_energy_splits(Bi);
		for(DShift = -DScan; DShift <= DScan; DShift += DStep){ // Doppler start
		detuning = detBase - DShift;
		detuning_probe = detBaseProbe - DShift;
		//Gamma_mazs = Gamma_mazs_base + (cc * ((DStep / 2) + abs(DShift)) / energija) / (lambda * 1e6);
		//Gamma_mazs = Gamma_mazs_base + 0.3*abs(DShift);
		//cout << "Dshift: " << DShift << "; Gamma_mazs: " << Gamma_mazs << endl;

		// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		buildLVSAlternative(Bi);
		int s;
		gsl_permutation * p = gsl_permutation_alloc (dimNum);
		gsl_linalg_complex_LU_decomp (ppMatrix, p, &s);
		gsl_linalg_complex_LU_solve (ppMatrix, p, ppVector, resVector);
		// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//rotate_dens_matrix(0,2);
		//count_bright_states(Rabi);
		//print_apdz(2);
		//saglabaa blivuma matricas
		// sprintf(dmfldname, "%s%s%.2f", fldname, "/dms-Rabi=", Rabi);
		// if(!FileExists(dmfldname)){ cout << "creating folder " << dmfldname << endl; mkdir(dmfldname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); }
		// sprintf(filename, "%s%s%.2f", dmfldname, "/densMatrix-det=m5-B=", B);
		// cout << "dm filename -> " << filename  << endl;
		// matrixFile.open(filename, ios::out | ios::trunc);
		// sprintf(filename, "%s%s%.2f", dmfldname, "/GrFgd-densMatrix-det=m5-B=", B);
		// cout << "dm gr filename -> " << filename  << endl;
		// matrixFileGr.open(filename, ios::out | ios::trunc);
		// sprintf(filename, "%s%s%.2f", dmfldname, "/ExFed-densMatrix-det=m5-B=", B);
		// cout << "dm exp filename -> " << filename  << endl;
		// matrixFileEx.open(filename, ios::out | ios::trunc);
		// save_density_matrix(B);
		// matrixFile.close();
		// matrixFileGr.close();
		// matrixFileEx.close();

		//collect_dens_matrix(0,(DStep * pow(e,-((DShift * DShift)/(2 * Sigma * Sigma)))/(sqrt(2 * Pi) * Sigma)));
		//gsl_vector_complex_fprintf (stdout, resVector, "%g");
		//print_apdz(1);
		/*double ksii = GSL_REAL(gsl_complex_add(GammaPCalcMinus(Fe_d, 0, Fg_d, 0, 0), GammaPCalcPlus(Fe_d, 0, Fg_d, 0, 0))) / Gamma_mazs;
		double aug = 0, apa = 0;
		int J = Fg_d;*/
		/*// RRL
		for(int m = -J; m <= J; m++){
			aug += (((J + 1) * (J + 1) -m * m) * ((J + 1) * (J + 1) + (J + 1) - 3 * m * m)) / (1 + (ksii * ((J + 1) * (J + 1) - m * m)) / ((J +1) * (2 * J + 1)));
			apa += (((J + 1) * (J + 1) - m * m) * (3 * J * J + 5 * J - m * m + 2)) / (1 + (ksii * ((J + 1) * (J + 1) - m * m)) / ((J +1) * (2 * J + 1)));
			}*/
		/*// QQL
		for(int m = -J; m <= J; m++){
			aug += (m * m * (3 * m * m - J * J - J)) / (1 + (ksii * m * m)/(J * (J + 1)));
			apa += (m * m * (m * m + J * J + J)) / (1 + (ksii * m * m)/(J * (J + 1)));
			}*/
		/*// QRL
		for(int m = -J; m <= J; m++){
			aug += (m * m * (J * J + J - 3 * m * m)) / (1 + (ksii * m * m)/(J * (J + 1)));
			apa += (m * m * (3 * J * J - J - m * m)) / (1 + (ksii * m * m)/(J * (J + 1)));
			}*/
		// PPL
		/*for(int m = -J; m <= J; m++){
			aug += ((J * J - m * m) * (J * J - J - 3 * m * m))/(1 + (ksii * (J * J - m * m))/(2 * J * J + J));
			apa += ((J * J - m * m) * (3 * J * J - m * m + J))/(1 + (ksii * (J * J - m * m))/(2 * J * J + J));
			}*/
		/*//QQC
		for(int m = -J; m <= J; m++){
			aug += ((J * J + J - m * m - m) * (m +1))/(1 + (ksii * (J * J - m * m - m))/(2 * J * (J + 1)));
			apa += ((J * J + J - m * m - m) * (J * J + J - m * m - 2 * m - 1))/(1 + (ksii * (J * J - m * m - m))/(2 * J * (J + 1)));
			}*/
		/*//RRC
		for(int m = -J; m <= J; m++){
			aug += ((2 * J + 1) * (m + 1) * (J + m + 1) * (J + m + 2))/(1 + (ksii * (J + m + 1) * (J + m + 2))/(2 * (J + 1) * (2 * J + 1)));
			apa += ((J * (J + 1) + m * (m + 2) + 1) * (J + m + 1) * (J + m + 2))/(1 + (ksii * (J + m + 1) * (J + m + 2))/(2 * (J + 1) * (2 * J + 1)));
			}*/
		/*//PPC
		for(int m = -J; m <= J; m++){
			aug += (-(J - m - 1) * (J - m) * (2 * J + 1) * (m + 1))/(1 + (ksii * (J - m - 1) * (J - m))/(2 * J * (2 * J + 1)));
			apa += ((J - m - 1) * (J - m) * (J * J + J + m * m + 2 * m +1))/(1 + (ksii * (J - m - 1) * (J - m))/(2 * J * (2 * J + 1)));
			}*/

		intensity_par += intensity(Bi, pol_o[0], teta_o[0], fii_o[0]) * (DStep * pow(e,-((DShift * DShift)/(2 * Sigma * Sigma)))/(sqrt(2 * Pi) * Sigma));
		intensity_per += intensity(Bi, pol_o[1], teta_o[1], fii_o[1]) * (DStep * pow(e,-((DShift * DShift)/(2 * Sigma * Sigma)))/(sqrt(2 * Pi) * Sigma));
		intensity_plus += intensity(Bi, 1, Pi/2, 0) * (DStep * pow(e,-((DShift * DShift)/(2 * Sigma * Sigma)))/(sqrt(2 * Pi) * Sigma));
		intensity_minus += intensity(Bi, -1, Pi/2, 0) * (DStep * pow(e,-((DShift * DShift)/(2 * Sigma * Sigma)))/(sqrt(2 * Pi) * Sigma));
		absorption_par += absorption_intensity(Bi, pol_p[0], teta_p[0], fii_p[0]) * (DStep * pow(e,-((DShift * DShift)/(2 * Sigma_probe * Sigma_probe)))/(sqrt(2 * Pi) * Sigma_probe));
		absorption_per += absorption_intensity(Bi, pol_p[1], teta_p[1], fii_p[1]) * (DStep * pow(e,-((DShift * DShift)/(2 * Sigma_probe * Sigma_probe)))/(sqrt(2 * Pi) * Sigma_probe));
	
		//dark += count_dark_states() * (DStep * pow(e,-((DShift * DShift)/(2 * Sigma * Sigma)))/(sqrt(2 * Pi) * Sigma));
		//dark += count_dark_states(0,2);
		//bark += count_dark_states(1,-1);
		//bright += count_bright_states(2) * 2 * (DStep * pow(e,-((DShift * DShift)/(2 * Sigma * Sigma)))/(sqrt(2 * Pi) * Sigma));
		//bright += count_bright_states(0);
		//right += count_bright_states(1) * 2;
		//sight += count_bright_states(2) * 2;
		//intensity_par = intensity_sel(Bi, Fe_d, Fg_d, 0, teta_i, fii_i);
		//intensity_per = intensity_sel(Bi, Fe_d, Fg_d, 0, teta_i + Pi/2, fii_i);
		/*intensity_par = intensity(Bi, 0, 0);
		intensity_per = intensity(Bi, 0, Pi/2, 0);
		polarity = (intensity_par - intensity_per) / (intensity_par + intensity_per);
		cout << "intensity sum: " << intensity_par + intensity_per << "; polarizaacija: " << polarity << endl;
		dataFile << DShift << " " << intensity_par + intensity_per << endl;*/
		} //Doppler end

	//cout << "teoreetiski pol: " << aug/apa << endl;

	cout.precision(10);
	polarity = (intensity_par - intensity_per);
	cout << "circ_1: " << intensity_par << "; circ_2: " << intensity_per << " ;intensity sum: " << intensity_par + intensity_per << "; total: " << intensity_per + intensity_z + intensity_par << "; starpiba: " << polarity << endl;
	polarity_abs = (absorption_par - absorption_per);
	cout << "abs circ_1: " << absorption_par << "; abs circ_2: " << absorption_per << " ; abs intensity sum: " << absorption_par + absorption_per << "; abs starpiba: " << polarity_abs << endl;
	//cout << "dark: " << dark << "; bark: " << bark << "; bright: " << bright << "; right: " << right << "; sight: " << sight << endl;
	dataFile << BB * BSign << " " << setprecision (9) << intensity_par + intensity_per << " " << intensity_par << " " << intensity_per << " " << intensity_par - intensity_per << endl;
	dataFileProbe << BB * BSign << " " << setprecision (9) << absorption_par + absorption_per << " " << absorption_par << " " << absorption_per << " " << absorption_par - absorption_per << endl;
	//dataFile << lShift << " " << setprecision (9) << intensity_par + intensity_per << endl;
    //treat_collected(BB);
	Bi++;
	//matrixFile << Rabi << " " << dark << " " << bark << endl;
	//matrixFile << Rabi << " " << bright << " " << right << " " << sight << endl;
	}
	sprintf(filename, "%s%s%s%s", fldname, "/", baseFilename, "Breit-Rabi-ex.txt");
	dataFileEx.open(filename, ios::out | ios::trunc);
	sprintf(filename, "%s%s%s%s", fldname, "/", baseFilename, "Breit-Rabi-gr.txt");
	dataFileGr.open(filename, ios::out | ios::trunc);
	save_level_splits();
	dataFileEx.close();
	dataFileGr.close();

	sprintf(filename, "%s%s%s%s", fldname, "/", baseFilename, "Probe-Breit-Rabi-ex.txt");
	dataFileExProbe.open(filename, ios::out | ios::trunc);
	sprintf(filename, "%s%s%s%s", fldname, "/", baseFilename, "Probe-Breit-Rabi-gr.txt");
	dataFileGrProbe.open(filename, ios::out | ios::trunc);
	save_level_splits_probe();
	dataFileExProbe.close();
	dataFileGrProbe.close();

	//Saglabā līmeņu sajaukšanos
	sprintf(filename, "%s%s%s%s", fldname, "/", baseFilename, "level-mixing-ex.txt");
	dataFileEx.open(filename, ios::out | ios::trunc);
	sprintf(filename, "%s%s%s%s", fldname, "/", baseFilename, "level-mixing-gr.txt");
	dataFileGr.open(filename, ios::out | ios::trunc);
	save_level_mixing();
	dataFileEx.close();
	dataFileGr.close();

	int polPr;
	for(polPr = -1; polPr <= 1; polPr += 1){
	sprintf(filename, "%s%s%s%s%i%s", fldname, "/", baseFilename, "trans-prob-pol", polPr,".txt");
	dataFilePr.open(filename, ios::out | ios::trunc);
	save_trans_prob(polPr, Fg_d);
	dataFilePr.close();}

	}
	//save_level_splits();
	dataFile.close();
	dataFileProbe.close();
	}
	}
	}

	//dataFileEx.close();
	//dataFileGr.close();
	//matrixFile.close();
	cout << "tas arī viss ..." << endl;
	return 0;
}
