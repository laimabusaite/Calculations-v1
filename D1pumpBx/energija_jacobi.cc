#include <stdlib.h> 
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

using std::cout; using std::endl;
double eigenValues[10], eigenValuesProbe[10]; //iipashveertiibas uz aaru
double eigenValuesSorted[10]; //iipashveertiibas sarkaartoti uz aaru
double eigenVectors[10][10], eigenVectorsProbe[10][10]; //iipashvektori uz aaru
double eigenVectorsSorted[10][10]; //iipashvektori sakaartoti uz aaru
double EM[10][10], EMp[10][10], EMs[10][10]; //energiju matrica
double vl[10]; //iipashveertiibas ieksheejai lietoshanai
double vc[10][10], vcp[10][10]; //iipashvektori iekshhejai lietoshanai
double RM[10][10]; //rotaacijas matrica
int i, j, k, l, m, mx, n, state; //indexi
double s, c, t, p, y, theta, sum, sumE, sumV; //roteeshanas parametri
int ind[10], changed[10]; //indexu masiivi


// Diraka - delta f-ja
//
int Ddelta(double i, double j){
	if (i == j) return 1; return 0;
}

/*double evec_decomp(int dim){
	// saakam ar apaksheejaa triistuura nulleeshanu
	for(i = 0; i < dim; i++){ //rindinja, no kuras njem diagonaalo elementu
		for(j = i + 1; j < dim; j++){ //rindinja, no kuras njem nedaigonaalo
			if(fabs(vcp[j][i]) > prec){ //vai shis vispaar atshkjiras no nulles
				c = (-1) * vcp[j][i] / vcp[i][i];
				for(k = 0; k <= dim; k++){ //nediagonaalaas rindinjas elementi
					vcp[j][k] += c;
				}
			}
		}
	}
	// turpinaam ar augsheejo triistuuri
	for(i = dim - 1; i >= 0; i--){ //rindinja, no kuras njem diagonaalo elementu
		for(j = i - 1; j >= 0; j--){ //rindinja, no kuras njem nedaigonaalo
			if(fabs(vcp[j][i]) > prec){ //vai shis vispaar atshkjiras no nulles
				c = (-1) * vcp[j][i] / vcp[i][i];
				for(k = 0; k <= dim; k++){ //nediagonaalaas rindinjas elementi
					vcp[j][k] += c;
				}
			}
		vcp[j][dim] = vcp[j][dim] / vcp[j][j];
		}
	}
}*/


int maxind(int k){
	int m;
	m = k + 1;
	for(i = k + 2; i < n; i++){
		if(fabs(EM[k][i]) > fabs(EM[k][m])){ m = i; }
	}
	return m;
}

double fill_RM(int h, int v, int dim){
	theta = atan(2 * EM[h][v] / (EM[h][h] - EM[v][v])) / 2;
	//cout << "theta: " << theta << "; EM: " << EM[h][v] << "; vl_sub: " << vl[h] - vl[v] << endl;
	for(i = 0; i < dim; i++){
		for(j = 0; j < dim; j++){
			RM[i][j] = Ddelta(i,j);
		}
	}
	RM[h][h] =  cos(theta);	RM[h][v] = sin(theta);
	RM[v][h] = -sin(theta);	RM[v][v] = cos(theta);
}

void rotate(int dim){
	for(i = 0; i < dim; i ++){
		for(j = 0; j < dim; j++){
			sumE = 0; sumV = 0;
			for(k = 0; k < dim; k++){
				sumE += EM[i][k] * RM[j][k];
				sumV += vc[i][k] * RM[j][k];
			}
			EMp[i][j] = sumE; vcp[i][j] = sumV;
		}
	}
	for(i = 0; i < dim; i ++){
		for(j = 0; j < dim; j++){
			sumE = 0; sumV = 0;
			for(k = 0; k < dim; k++){
				sumE += RM[i][k] * EMp[k][j];
				//sumV += RM[i][k] * vcp[k][j];
			}
			EM[i][j] = sumE; vc[i][j] = vcp[i][j];
		}
	}
}

double el_en_J(int m, int f, int ff, int level, double B){
	double retVal;
	retVal = mi_Bora_MHz * gj[level] * B * pow(-1., JJ[level] + ik + f + ff - m + 1) * sqrt((2*f + 1)*(2*ff + 1) * JJ[level] * (JJ[level] + 1) * (2*JJ[level] + 1)) * three_j_wrap(2*f, 2, 2*ff, -2*m, 0, 2*m) * six_j_wrap(2*JJ[level], 2*f, 2*ik, 2*ff, 2*JJ[level], 2);
	return retVal;
}

double el_en_I(int m, int f, int ff, int level, double B){
	double retVal;
	retVal = mi_Bora_MHz * gi * B * pow(-1., JJ[level] + ik + f + ff - m + 1) * sqrt((2*f + 1)*(2*ff + 1) * ik * (ik + 1) * (2*ik + 1)) * three_j_wrap(2*f, 2, 2*ff, -2*m, 0, 2*m) * six_j_wrap(2*ik, 2*f, LL[level], 2*ff, 2*ik, 2);
	return retVal;
}

void energySplits(int m, int FMin, int FMax, int level, double B){
    //if(level == 0) B = 0;
	int FMU;
	if(abs(m) > FMin) FMU = abs(m);
	else FMU = FMin;
	int diff = FMU - FMin;
	int dim = FMax - FMU + 1;
	n = dim;
	//cout << " m: " << m << "; dim: " << dim << endl;
	state = dim;
	int Fp, ii, jj, kk;
	double C;
	double energySelf;
	double energySide;
	double norm;
	double swap;
	double sides[10];
	/*gsl_matrix * EMtrans = gsl_matrix_calloc((int)(2*JJ[level]+1), (int)(2*JJ[level]+1));
	gsl_matrix * tmp = gsl_matrix_calloc((int)(2*JJ[level]+1), (int)(2*JJ[level]+1));
	for(int F = FMU; F <= FMax; F++){
		for(int FF = FMU; FF <= FMax; FF++){
			gsl_matrix_set(EMtrans, F - FMU, FF - FMU, el_en_J(m, F, FF, level, Btrans) + el_en_I(m, F, FF, level, Btrans));
		}
	}
	if(JJ[level] > 1){
    	gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1, dMatrixE, EMtrans, 0, tmp);
    	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, tmp, dMatrixE, 0, EMtrans);
	}
	else{
    	gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1, dMatrixG, EMtrans, 0, tmp);
    	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, tmp, dMatrixG, 0, EMtrans);
	}*/
	for(ii = 0; ii < 10; ii++){
		vl[ii] = 0.; ind[ii] = 0.; changed[ii] = 0.;
		eigenValues[ii] = 0.;
		eigenValuesSorted[ii] = 0.;
		for(jj = 0; jj < 10; jj++){
			vc[ii][jj] = Ddelta(ii,jj); eigenVectors[ii][jj] = Ddelta(ii,jj);
			EM[ii][jj] = 0.;
		}
	}
	for(int F = FMU; F <= FMax; F++){
		for(int FF = FMU; FF <= FMax; FF++){
			if(F == FF){
				C = p1(F) - p1(ik) - p1(JJ[level]);
				EM[F - FMU][FF - FMU] += aa[level] * C / 2;
				if(2 * JJ[level] - 1 > 0) EM[F - FMU][FF - FMU] += bb[level] * (1.5 * p1(C) - 2 * p1(ik) * p1(JJ[level])) / (4 * ik * (2 * ik -1) * JJ[level] * (2 * JJ[level] - 1));
			}
			EM[F - FMU][FF - FMU] += el_en_J(m, F, FF, level, B) + el_en_I(m, F, FF, level, B)/* + gsl_matrix_get(EMtrans, F - FMU, FF - FMU)*/;
		}
	}
	/*for(int F = FMU; F <= FMax; F++){
		for(int FF = FMU; FF <= FMax; FF++){
			cout << "m: " << m << "; F: " << F << "; FF: " << FF << "; en: " << EM[F - FMU][FF - FMU] << endl;
		}
	}
	cout << endl;*/
	/*if(dim == 4){
		EM[0][0] = 4;
		EM[0][1] = -30;
		EM[0][2] = 60;
		EM[0][3] = -35;

		EM[1][0] = -30;
		EM[1][1] = 300;
		EM[1][2] = -675;
		EM[1][3] = 420;

		EM[2][0] = 60;
		EM[2][1] = -675;
		EM[2][2] = 1620;
		EM[2][3] = -1050;

		EM[3][0] = -35;
		EM[3][1] = 420;
		EM[3][2] = -1050;
		EM[3][3] = 700;
	}*/
	for(k = 0; k < n; k++){
		ind[k] = maxind(k);
	}
	if(fabs(B) > 0){

		do{
			//cout << state << endl;
			mx = 0;
			for(k = 1; k < n - 1; k++){
				if(fabs(EM[k][ind[k]]) > fabs(EM[mx][ind[mx]])) mx = k;
			}
			l = ind[mx];
			fill_RM(mx,l,dim);
			rotate(dim);
			sum = 0;
			for(i = 0; i < dim; i++){
				for(j = i + 1; j < dim; j++){
					sum += EM[i][j];
				}
			}
			ind[mx] = maxind(mx);
			ind[l] = maxind(l);
		} while(fabs(sum) > prec);
	}
	for(ii = 0; ii < dim; ii++){
		eigenValues[ii] = EM[ii][ii];
		eigenValuesSorted[ii] = EM[ii][ii];
		/*n = 0;
		for(i = 0; i < dim; i++){ //kuru iipashvektoru paarbaudaam
			state = 1;
			for(j = 0; j < dim; j++){ //iipashvektora poziicija
				s = vc[j][i] * eigenValues[ii];
				c = 0;
				for(k = 0; k < dim; k++){
					c += EMs[j][k] * vc[k][i];
				}
				if(fabs(c - s) > prec * 100) state = 0;
				//cout << c - s << " " << state << endl;
			}
			//cout << endl;
			if(state == 1){ ind[ii] = i; n += 1; }
		}
		if(n < 1) cout << "Nav atrasts iipashveertiibai atbilstosh vektors !!!" << endl;
		if(n > 1) cout << "Atrasts vairaak kaa viens vektros iipashveertiibai !!!" << endl;*/
	}
	//if(dim == 4){
	for(ii = 0; ii < dim; ii++){
		//cout << "eigen: " << eigenValues[ii] << endl;
		for(jj = 0; jj < dim; jj++){
			eigenVectors[ii][jj] = vc[ii][jj];
			eigenVectorsSorted[ii][jj] = vc[ii][jj];
			//cout << vc[jj][ii] << endl;
		}
	}
	for(ii = 0; ii < dim; ii++){
		for(jj = dim - 1; jj > ii; jj--){
			if(eigenValues[jj-1] > eigenValues[jj]){
				swap = eigenValues[jj-1];
				eigenValues[jj-1] = eigenValues[jj];
				eigenValues[jj] = swap;
				for(kk = 0; kk < dim; kk++){
					swap = eigenVectors[kk][jj-1];
					eigenVectors[kk][jj-1] = eigenVectors[kk][jj];
					eigenVectors[kk][jj] = swap;
				}
				/*for(kk = 0; kk < dim; kk++){
					swap = eigenVectors[kk][jj-1];
					eigenVectors[kk][jj-1] = eigenVectors[kk][jj];
					eigenVectors[kk][jj] = swap;
				}*/
			}
		}
	}
}

double el_en_J_probe(int m, int f, int ff, int level, double B){
	double retVal;
	retVal = mi_Bora_MHz * gj_probe[level] * B * pow(-1., JJ_probe[level] + ik_probe + f + ff - m + 1) * sqrt((2*f + 1)*(2*ff + 1) * JJ_probe[level] * (JJ_probe[level] + 1) * (2*JJ_probe[level] + 1)) * three_j_wrap(2*f, 2, 2*ff, -2*m, 0, 2*m) * six_j_wrap(2*JJ_probe[level], 2*f, 2*ik_probe, 2*ff, 2*JJ_probe[level], 2);
	return retVal;
}

double el_en_I_probe(int m, int f, int ff, int level, double B){
	double retVal;
	retVal = mi_Bora_MHz * gi_probe * B * pow(-1., JJ_probe[level] + ik_probe + f + ff - m + 1) * sqrt((2*f + 1)*(2*ff + 1) * ik_probe * (ik_probe + 1) * (2*ik_probe + 1)) * three_j_wrap(2*f, 2, 2*ff, -2*m, 0, 2*m) * six_j_wrap(2*ik_probe, 2*f, LL_probe[level], 2*ff, 2*ik_probe, 2);
	return retVal;
}

void energySplitsProbe(int m, int FMin, int FMax, int level, double B){
    //if(level == 0) B = 0;
	int FMU;
	if(abs(m) > FMin) FMU = abs(m);
	else FMU = FMin;
	int diff = FMU - FMin;
	int dim = FMax - FMU + 1;
	n = dim;
	//cout << " m: " << m << "; dim: " << dim << endl;
	state = dim;
	int Fp, ii, jj, kk;
	double C;
	double energySelf;
	double energySide;
	double norm;
	double swap;
	double sides[10];
	/*gsl_matrix * EMtrans = gsl_matrix_calloc((int)(2*JJ[level]+1), (int)(2*JJ[level]+1));
	gsl_matrix * tmp = gsl_matrix_calloc((int)(2*JJ[level]+1), (int)(2*JJ[level]+1));
	for(int F = FMU; F <= FMax; F++){
		for(int FF = FMU; FF <= FMax; FF++){
			gsl_matrix_set(EMtrans, F - FMU, FF - FMU, el_en_J(m, F, FF, level, Btrans) + el_en_I(m, F, FF, level, Btrans));
		}
	}
	if(JJ[level] > 1){
    	gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1, dMatrixE, EMtrans, 0, tmp);
    	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, tmp, dMatrixE, 0, EMtrans);
	}
	else{
    	gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1, dMatrixG, EMtrans, 0, tmp);
    	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, tmp, dMatrixG, 0, EMtrans);
	}*/
	for(ii = 0; ii < 10; ii++){
		vl[ii] = 0.; ind[ii] = 0.; changed[ii] = 0.;
		eigenValuesProbe[ii] = 0.;
		eigenValuesSorted[ii] = 0.;
		for(jj = 0; jj < 10; jj++){
			vc[ii][jj] = Ddelta(ii,jj); eigenVectorsProbe[ii][jj] = Ddelta(ii,jj);
			EM[ii][jj] = 0.;
		}
	}
	for(int F = FMU; F <= FMax; F++){
		for(int FF = FMU; FF <= FMax; FF++){
			if(F == FF){
				C = p1(F) - p1(ik_probe) - p1(JJ_probe[level]);
				EM[F - FMU][FF - FMU] += aa_probe[level] * C / 2;
				if(2 * JJ_probe[level] - 1 > 0) EM[F - FMU][FF - FMU] += bb_probe[level] * (1.5 * p1(C) - 2 * p1(ik_probe) * p1(JJ_probe[level])) / (4 * ik_probe * (2 * ik_probe -1) * JJ_probe[level] * (2 * JJ_probe[level] - 1));
			}
			EM[F - FMU][FF - FMU] += el_en_J_probe(m, F, FF, level, B) + el_en_I_probe(m, F, FF, level, B)/* + gsl_matrix_get(EMtrans, F - FMU, FF - FMU)*/;
		}
	}
	/*for(int F = FMU; F <= FMax; F++){
		for(int FF = FMU; FF <= FMax; FF++){
			cout << "m: " << m << "; F: " << F << "; FF: " << FF << "; en: " << EM[F - FMU][FF - FMU] << endl;
		}
	}
	cout << endl;*/
	// if(dim == 4){
	// 	EM[0][0] = 4;
	// 	EM[0][1] = -30;
	// 	EM[0][2] = 60;
	// 	EM[0][3] = -35;

	// 	EM[1][0] = -30;
	// 	EM[1][1] = 300;
	// 	EM[1][2] = -675;
	// 	EM[1][3] = 420;

	// 	EM[2][0] = 60;
	// 	EM[2][1] = -675;
	// 	EM[2][2] = 1620;
	// 	EM[2][3] = -1050;

	// 	EM[3][0] = -35;
	// 	EM[3][1] = 420;
	// 	EM[3][2] = -1050;
	// 	EM[3][3] = 700;
	// }
	for(k = 0; k < n; k++){
		ind[k] = maxind(k);
	}
	if(fabs(B) > 0){

		do{
			//cout << state << endl;
			mx = 0;
			for(k = 1; k < n - 1; k++){
				if(fabs(EM[k][ind[k]]) > fabs(EM[mx][ind[mx]])) mx = k;
			}
			l = ind[mx];
			fill_RM(mx,l,dim);
			rotate(dim);
			sum = 0;
			for(i = 0; i < dim; i++){
				for(j = i + 1; j < dim; j++){
					sum += EM[i][j];
				}
			}
			ind[mx] = maxind(mx);
			ind[l] = maxind(l);
		} while(fabs(sum) > prec);
	}
	for(ii = 0; ii < dim; ii++){
		eigenValuesProbe[ii] = EM[ii][ii];
		eigenValuesSorted[ii] = EM[ii][ii];
		/*n = 0;
		for(i = 0; i < dim; i++){ //kuru iipashvektoru paarbaudaam
			state = 1;
			for(j = 0; j < dim; j++){ //iipashvektora poziicija
				s = vc[j][i] * eigenValues[ii];
				c = 0;
				for(k = 0; k < dim; k++){
					c += EMs[j][k] * vc[k][i];
				}
				if(fabs(c - s) > prec * 100) state = 0;
				//cout << c - s << " " << state << endl;
			}
			//cout << endl;
			if(state == 1){ ind[ii] = i; n += 1; }
		}
		if(n < 1) cout << "Nav atrasts iipashveertiibai atbilstosh vektors !!!" << endl;
		if(n > 1) cout << "Atrasts vairaak kaa viens vektros iipashveertiibai !!!" << endl;*/
	}
	//if(dim == 4){
	for(ii = 0; ii < dim; ii++){
		//cout << "eigen: " << eigenValues[ii] << endl;
		for(jj = 0; jj < dim; jj++){
			eigenVectorsProbe[ii][jj] = vc[ii][jj];
			eigenVectorsSorted[ii][jj] = vc[ii][jj];
			//cout << vc[jj][ii] << endl;
		}
	}
	for(ii = 0; ii < dim; ii++){
		for(jj = dim - 1; jj > ii; jj--){
			if(eigenValuesProbe[jj-1] > eigenValuesProbe[jj]){
				swap = eigenValuesProbe[jj-1];
				eigenValuesProbe[jj-1] = eigenValuesProbe[jj];
				eigenValuesProbe[jj] = swap;
				for(kk = 0; kk < dim; kk++){
					swap = eigenVectorsProbe[kk][jj-1];
					eigenVectorsProbe[kk][jj-1] = eigenVectorsProbe[kk][jj];
					eigenVectorsProbe[kk][jj] = swap;
				}
				/*for(kk = 0; kk < dim; kk++){
					swap = eigenVectors[kk][jj-1];
					eigenVectors[kk][jj-1] = eigenVectors[kk][jj];
					eigenVectors[kk][jj] = swap;
				}*/
			}
		}
	}
}
