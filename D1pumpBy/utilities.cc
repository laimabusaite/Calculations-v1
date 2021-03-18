//
// funkcija, kas izdrukaa liimenu energijas
void print_energy_splits(int Bi){
	for(int Fg = FStartGr; Fg <= FEndGr; Fg++){
		for(int m = -Fg; m <= Fg; m++){
			cout << "Fg: " << Fg << "; m: " << m << "; energy " << energySplitsMatrix[0][Fg - FStartGr][m + FEndGr][Bi] << endl;
		}
	}
	cout << endl;
	for(int Fe = FStartEx; Fe <= FEndEx; Fe++){
		for(int M = -Fe; M <= Fe; M++){
			cout << "Fe: " << Fe << "; M: " << M << "; energy " << energySplitsMatrix[1][Fe - FStartEx][M + FEndEx][Bi] << endl;
		}
	}
}
//
// funkcija, kas saglabaa liimenu energijas
void save_level_splits(){
	dataFileGr << "Mag "; int add = 0;
	for(int Fg = FStartGr; Fg <= FEndGr; Fg++){
		for(int m = -Fg; m <= Fg; m++){
			dataFileGr << Fg << "_" << m << " ";
		}
	}
	dataFileGr << "\n";
	dataFileEx << "Mag ";
	for(int Fe = FStartEx; Fe <= FEndEx; Fe++){
		for(int M = -Fe; M <= Fe; M++){
			dataFileEx << Fe << "_" << M << " ";
		}
	}
	dataFileEx << "\n";
	for(int i = 0; i <= steps; i++){
		//if(BSign > 0) add = steps + 1;
		//else add = 0;
		dataFileGr << i * step * BSign << " ";
		for(int Fg = FStartGr; Fg <= FEndGr; Fg++){
			for(int m = -Fg; m <= Fg; m++){
				dataFileGr << energySplitsMatrix[0][Fg - FStartGr][m + FEndGr][i + add] << " ";
			}
		}
		dataFileGr << "\n";
		dataFileEx << i * step * BSign << " ";
		for(int Fe = FStartEx; Fe <= FEndEx; Fe++){
			for(int M = -Fe; M <= Fe; M++){
				dataFileEx << energySplitsMatrix[1][Fe - FStartEx][M + FEndEx][i + add] << " ";
			}
		}
		dataFileEx << "\n";
		//cout << BSign << endl;
	}
}
//Probe funkcija, kas saglabaa liimenu energijas
void save_level_splits_probe(){
	dataFileGrProbe << "Mag "; int add = 0;
	for(int Fg = FStartGr_probe; Fg <= FEndGr_probe; Fg++){
		for(int m = -Fg; m <= Fg; m++){
			dataFileGrProbe << Fg << "_" << m << " ";
		}
	}
	dataFileGrProbe << "\n";
	dataFileExProbe << "Mag ";
	for(int Fe = FStartEx_probe; Fe <= FEndEx_probe; Fe++){
		for(int M = -Fe; M <= Fe; M++){
			dataFileExProbe << Fe << "_" << M << " ";
		}
	}
	dataFileExProbe << "\n";
	for(int i = 0; i <= steps; i++){
		//if(BSign > 0) add = steps + 1;
		//else add = 0;
		dataFileGrProbe << i * step * BSign << " ";
		for(int Fg = FStartGr_probe; Fg <= FEndGr_probe; Fg++){
			for(int m = -Fg; m <= Fg; m++){
				dataFileGrProbe << energySplitsMatrixProbe[0][Fg - FStartGr_probe][m + FEndGr_probe][i + add] << " ";
			}
		}
		dataFileGrProbe << "\n";
		dataFileExProbe << i * step * BSign << " ";
		for(int Fe = FStartEx_probe; Fe <= FEndEx_probe; Fe++){
			for(int M = -Fe; M <= Fe; M++){
				dataFileExProbe << energySplitsMatrixProbe[1][Fe - FStartEx_probe][M + FEndEx_probe][i + add] << " ";
			}
		}
		dataFileExProbe << "\n";
		//cout << BSign << endl;
	}
}
//
// funkcija, kas saglabaa liimenu sajaukshanaas sheemu
void save_level_mixing(){
	dataFileGr << "Mag "; int add = 0;
	for(int Fg = FStartGr; Fg <= FEndGr; Fg++){
		for(int m = -Fg; m <= Fg; m++){
			for(int FgMix = FStartGr; FgMix <= FEndGr; FgMix++){
				dataFileGr << Fg << "_" << FgMix  << "m=" << m << " ";
				//cout << "Fg: " << Fg << "; m: " << m << "; FgMix" << FgMix << "; value = " << levelMixingMatrix[0][Fg - FStartGr][FgMix - FStartGr][m + FEndGr][i] << endl;
			}
		}
	}
	dataFileEx << "Mag ";
	for(int Fe = FStartEx; Fe <= FEndEx; Fe++){
		for(int m = -Fe; m <= Fe; m++){
			for(int FeMix = FStartEx; FeMix <= FEndEx; FeMix++){
				dataFileEx << Fe << "_" << FeMix << "m=" << m << " ";
			}
		}
	}
	dataFileEx << "\n";
	dataFileGr << "\n";
	for(int i = 0; i <= steps; i++){
		dataFileGr << i * step * BSign << " ";
		//if(BSign > 0) add = steps + 1;
		//else add = 0;
		for(int Fg = FStartGr; Fg <= FEndGr; Fg++){
			for(int m = -Fg; m <= Fg; m++){
				for(int FgMix = FStartGr; FgMix <= FEndGr; FgMix++){
					dataFileGr << levelMixingMatrix[0][Fg - FStartGr][FgMix - FStartGr][m + FEndGr][i + add] << " ";
				}
			}
		}
		dataFileEx << i * step * BSign << " ";
		for(int Fe = FStartEx; Fe <= FEndEx; Fe++){
			for(int m = -Fe; m <= Fe; m++){
				for(int FeMix = FStartEx; FeMix <= FEndEx; FeMix++){
					dataFileEx << levelMixingMatrix[1][Fe - FStartEx][FeMix - FStartEx][m + FEndEx][i + add] << " ";
				}
			}
		}
		dataFileEx << "\n";
		dataFileGr << "\n";
	}
}
//
// funkcija, kas saglabaa mag liimenu paareju varbuutiibas
// argumenti:
// 1: pol - paarejas polarizaacija 0, +-1
// 2: Fg - pamatstaavokla F
void save_trans_prob(int pol, int Fg){
	int m;
	string Sign;
	string sign;
	double temp, ttemp;
	dataFilePr << "Mag ";
	for(int Fe = FStartEx; Fe <= FEndEx; Fe++){
		for(int M = -Fe; M <= Fe; M++){
			m = M - pol;
			if(abs(m) <= Fg){
				if(M < 0) Sign = "m";
				else if(M > 0) Sign = "p";
				else Sign = "";
				if(m < 0) sign = "m";
				else if(m > 0) sign = "p";
				else sign = "";
				dataFilePr << Fe << Sign << M << "G" << Fg << sign << m << " ";
			}
		}
	}
	dataFilePr << "\n";
	for(int i = 0; i <= steps; i++){
		dataFilePr << i * step * BSign << " ";
		for(int Fe = FStartEx; Fe <= FEndEx; Fe++){
			for(int M = -Fe; M <= Fe; M++){
				m = M - pol;
				if(abs(m) <= Fg){
					temp = 0.0;
					for(int F1mix = FStartEx; F1mix <= FEndEx; F1mix++){ 
						for(int F2mix = FStartGr; F2mix <= FEndGr; F2mix++){
							ttemp = levelMixingMatrix[1][Fe - FStartEx][F1mix - FStartEx][M + FEndEx][i] * levelMixingMatrix[0][Fg - FStartGr][F2mix - FStartGr][m + FEndGr][i] * three_j_matrix[F1mix][-M][F2mix][m] * reduced_matrix[F1mix][F2mix];
							temp += abs(ttemp);
						}
					}
					dataFilePr << temp * temp << " ";
				}
			}
		}
		dataFilePr << "\n";
	}
}
//
// funkcija, kas izdrukaa apdziivotiibas
//
double print_apdz(int type = 0){
	double apdzGrCnt = 0;
	double apdzExCnt = 0;
	gsl_complex tmp;
	for(int Fg = FStartGr; Fg <= FEndGr; Fg++){
		for(int m = -Fg; m <= Fg; m++){
			for(int Fgp = FStartGr; Fgp <= FEndGr; Fgp++){
				for(int mp = -Fgp; mp <= Fgp; mp++){
					tmp = gsl_vector_complex_get(resVector, vietaGr[Fg][m][Fgp][mp]);
					if(Fg == Fgp && m == mp) apdzGrCnt += GSL_REAL(tmp);
					if(Fg == Fgp && m == mp && type == 1) cout << "Fg -> " << Fg << "; m -> " << m << "; mp -> " << mp << "; tmp -> " << GSL_REAL(tmp) << " + " << GSL_IMAG(tmp) << " * i " << endl;
					else if(Fg == Fgp && (abs(m - mp) == 2 || abs(m - mp) == 4) &&  type == 2)  cout << "Fg -> " << Fg << "; m -> " << m << "; mp -> " << mp << "; tmp -> " << GSL_REAL(tmp) << " + " << GSL_IMAG(tmp) << " * i " << endl;
					else if(Fg == Fgp &&  type == 3)  cout << "Fg -> " << Fg << "; m -> " << m << "; mp -> " << mp << "; tmp -> " << GSL_REAL(tmp) << " + " << GSL_IMAG(tmp) << " * i " << endl;
				}
			}
		}
	}
	cout << endl;
	for(int Fe = FStartEx; Fe <= FEndEx; Fe++){
		for(int M = -Fe; M <= Fe; M++){
			for(int Fep = FStartEx; Fep <= FEndEx; Fep++){
				for(int Mp = -Fep; Mp <= Fep; Mp++){
					tmp = gsl_vector_complex_get(resVector, vietaEx[Fe][M][Fep][Mp]);
					if(Fe == Fep && M == Mp) apdzExCnt += GSL_REAL(tmp);
					if(Fe == Fep && M == Mp && type == 1) cout << "Fe -> " << Fe << "; M -> " << M <<  "; Mp -> " << Mp << "; tmp -> " << GSL_REAL(tmp) << " + " << GSL_IMAG(tmp) << " * i " << endl;
					//if(Fe == Fep && type == 2) cout << "Fe -> " << Fe << "; M -> " << M <<  "; Mp -> " << Mp << "; tmp -> " << GSL_REAL(tmp) << " + " << GSL_IMAG(tmp) << " * i " << endl;
				}
			}
		}
	}
	cout << endl << "ground -> " << apdzGrCnt << "; excited -> " << apdzExCnt << "; sum -> " << apdzGrCnt + apdzExCnt << endl;
}
//
// funkcija, kas izdrukaa apdziivotiibas
//
double print_apdz_sum(int type = 0){
	double apdzGrCnt = 0;
	double apdzExCnt = 0;
	gsl_complex tmp;
	if(type == 0){
	for(int Fg = FStartGr; Fg <= FEndGr; Fg++){
		for(int m = -Fg; m <= Fg; m++){
			tmp = gsl_vector_complex_get(resVector, vietaGr[Fg][m][Fg][m]);
			apdzGrCnt += GSL_REAL(tmp);
		}
		cout << "Fg = " << Fg << "; apdz = " << apdzGrCnt << endl;
		apdzGrCnt = 0;
	}
	cout << endl;}
	if(type == 1){
	for(int Fe = FStartEx; Fe <= FEndEx; Fe++){
		for(int M = -Fe; M <= Fe; M++){
			tmp = gsl_vector_complex_get(resVector, vietaEx[Fe][M][Fe][M]);
			apdzExCnt += GSL_REAL(tmp);
		}
		cout << "Fe = " << Fe << "; apdz = " << apdzExCnt << endl;
		apdzExCnt = 0;
	}
	cout << endl;}
}
//
// fluorescences intensitaate
double intensity(int step, int pol_u = pol, double teta_u = teta_i, double fii_u = fii_i){
	gsl_complex temp, tempElement, FF;
	double intensityR = 0.0, intensityI = 0.0;
	int q1, q2;
	for(int Fg = FStartGr; Fg <= FEndGr; Fg++){
		for(int m = -Fg; m <= Fg; m++){
			for(int Fe = FStartEx; Fe <= FEndEx; Fe++){
				for(int M = -Fe; M <= Fe; M++){
					for(int Fep = FStartEx; Fep <= FEndEx; Fep++){
						for(int Mp = -Fep; Mp <= Fep; Mp++){
							q1 = Mp - m;
							q2 = M - m;
							if(abs(q1) < 2 && abs(q2) < 2){
								/*E_mul = gsl_complex_mul(gsl_complex_conjugate(E(q1, pol_u, teta_u, fii_u)), E(q2, pol_u, teta_u, fii_u));
								for(int F1mix = FStartGr; F1mix <= FEndGr; F1mix++){
									for(int F2mix = FStartEx; F2mix <= FEndEx; F2mix++){
										for(int F3mix = FStartEx; F3mix <= FEndEx; F3mix++){
											tempElement = gsl_vector_complex_get(resVector, vietaEx[F2mix][M][F3mix][Mp]);
											FF = three_j_wrap(2 * F3mix, 2 * 1, 2 * F1mix, -(2 * Mp), 2 * q1, 2 * m) * three_j_wrap(2 * F2mix, 2 * 1, 2 * F1mix, -(2 * M), 2 * q2, 2 * m) * reducedJ(F3mix, F1mix, JJ[1], JJ[0]) * reducedJ(F2mix, F1mix, JJ[1], JJ[0]) * levelMixingMatrix[0][Fg - FStartGr][F1mix - FStartGr][m + FEndGr][step] * levelMixingMatrix[1][Fe - FStartEx][F2mix - FStartEx][M + FEndEx][step] * levelMixingMatrix[1][Fep - FStartEx][F2mix - FStartEx][Mp + FEndEx][step];
											temp = gsl_complex_mul_real(gsl_complex_mul(E_mul, tempElement), FF);
											intensityR += GSL_REAL(temp);
											intensityI += GSL_IMAG(temp);
										}
									}
								}
								*/
								FF = gsl_complex_mul(baseElementSimpleInd(Fg, Fep, m, Mp, step, pol_u, teta_u, fii_u), baseElementSimple(Fe, Fg, M, m, step, pol_u, teta_u, fii_u));
								tempElement = gsl_complex_rect(0,0);
								/*for(int F1mix = FStartEx; F1mix <= FEndEx; F1mix++){
									for(int F2mix = FStartEx; F2mix <= FEndEx; F2mix++){
										tempElement = gsl_complex_add(tempElement, gsl_complex_mul_real(gsl_vector_complex_get(resVector, vietaEx[F1mix][M][F2mix][Mp]),levelMixingMatrix[1][F1mix - FStartEx][Fe - FStartEx][M + FEndEx][step] * levelMixingMatrix[1][F2mix - FStartEx][Fe - FStartEx][Mp + FEndEx][step]));
									}
								}*/
								tempElement = gsl_vector_complex_get(resVector, vietaEx[Fep][Mp][Fe][M]);
								temp = gsl_complex_mul(tempElement, FF);
								intensityR += GSL_REAL(temp);
								intensityI += GSL_IMAG(temp);
							}
						}
					}
				}
			}
		}
	}
	//cout << "intensityR: " << intensityR << "; intensityI: " << intensityI << endl;
	return intensityR;
}
//
// fluorescences intensitaate atsevishkjai paarejai
double intensity_sel(int step, int Fe, int Fg, int pol_u = pol, double teta_u = teta_i, double fii_u = fii_i){
	gsl_complex temp, tempElement, E_mul;
	double intensityR = 0.0, intensityI = 0.0, FF;
	int q1, q2;
	for(int m = -Fg; m <= Fg; m++){
		for(int M = -Fe; M <= Fe; M++){
			for(int Mp = -Fe; Mp <= Fe; Mp++){
				q1 = Mp - m;
				q2 = M - m;
				if(abs(q1) < 2 && abs(q2) < 2){
					E_mul = gsl_complex_mul(gsl_complex_conjugate(E(q1, pol_u, teta_u, fii_u)), E(q2, pol_u, teta_u, fii_u));
					for(int F1mix = FStartGr; F1mix <= FEndGr; F1mix++){
						for(int F2mix = FStartEx; F2mix <= FEndEx; F2mix++){
							for(int F3mix = FStartEx; F3mix <= FEndEx; F3mix++){
								tempElement = gsl_vector_complex_get(resVector, vietaEx[F2mix][M][F3mix][Mp]);
								FF = three_j_wrap(2 * F3mix, 2 * 1, 2 * F1mix, -(2 * Mp), 2 * q1, 2 * m) * three_j_wrap(2 * F2mix, 2 * 1, 2 * F1mix, -(2 * M), 2 * q2, 2 * m) * reducedJ(F3mix, F1mix, JJ[1], JJ[0]) * reducedJ(F2mix, F1mix, JJ[1], JJ[0]) * levelMixingMatrix[0][Fg - FStartGr][F1mix - FStartGr][m + FEndGr][step] * levelMixingMatrix[1][Fe - FStartEx][F2mix - FStartEx][M + FEndEx][step] * levelMixingMatrix[1][Fe - FStartEx][F3mix - FStartEx][Mp + FEndEx][step];
								temp = gsl_complex_mul_real(gsl_complex_mul(E_mul, tempElement), FF);
								intensityR += GSL_REAL(temp);
								intensityI += GSL_IMAG(temp);
							}
						}
					}
				}
			}
		}
	}
	//cout << "intensityR: " << intensityR << "; intensityI: " << intensityI << endl;
	return intensityR;
}

// absorbcijas intensitaate
double absorption_intensity(int step, int pol_u = pol, double teta_u = teta_i, double fii_u = fii_i){
	gsl_complex temp0, temp, tempElement, FF, G_probe;
	double intensityR = 0.0, intensityI = 0.0;
	int q1, q2;
	for(int Fe = FStartEx_probe; Fe <= FEndEx_probe; Fe++){
		for(int M = -Fe; M <= Fe; M++){
			for(int Fg = FStartGr_probe; Fg <= FEndGr_probe; Fg++){
				for(int m = -Fg; m <= Fg; m++){
					for(int Fgp = FStartGr_probe; Fgp <= FEndGr_probe; Fgp++){
						for(int mp = -Fgp; mp <= Fgp; mp++){
							q1 = M - mp;
							q2 = M - m;
							if(abs(q1) < 2 && abs(q2) < 2){
								/*E_mul = gsl_complex_mul(gsl_complex_conjugate(E(q1, pol_u, teta_u, fii_u)), E(q2, pol_u, teta_u, fii_u));
								for(int F1mix = FStartGr; F1mix <= FEndGr; F1mix++){
									for(int F2mix = FStartEx; F2mix <= FEndEx; F2mix++){
										for(int F3mix = FStartEx; F3mix <= FEndEx; F3mix++){
											tempElement = gsl_vector_complex_get(resVector, vietaEx[F2mix][M][F3mix][Mp]);
											FF = three_j_wrap(2 * F3mix, 2 * 1, 2 * F1mix, -(2 * Mp), 2 * q1, 2 * m) * three_j_wrap(2 * F2mix, 2 * 1, 2 * F1mix, -(2 * M), 2 * q2, 2 * m) * reducedJ(F3mix, F1mix, JJ[1], JJ[0]) * reducedJ(F2mix, F1mix, JJ[1], JJ[0]) * levelMixingMatrix[0][Fg - FStartGr][F1mix - FStartGr][m + FEndGr][step] * levelMixingMatrix[1][Fe - FStartEx][F2mix - FStartEx][M + FEndEx][step] * levelMixingMatrix[1][Fep - FStartEx][F2mix - FStartEx][Mp + FEndEx][step];
											temp = gsl_complex_mul_real(gsl_complex_mul(E_mul, tempElement), FF);
											intensityR += GSL_REAL(temp);
											intensityI += GSL_IMAG(temp);
										}
									}
								}
								*/
								FF = gsl_complex_mul(baseElementSimpleProbe(Fe, Fg, M, m, step, pol_u, teta_u, fii_u), baseElementSimpleIndProbe(Fgp, Fe, mp, M, step, pol_u, teta_u, fii_u)); //gsl_complex baseElementSimpleInd(int Fg, int Fe, int m, int M, int step, int pol_u = pol, double teta_u = teta_i, double fii_u = fii_i)
								tempElement = gsl_complex_rect(0,0);
								/*for(int F1mix = FStartEx; F1mix <= FEndEx; F1mix++){
									for(int F2mix = FStartEx; F2mix <= FEndEx; F2mix++){
										tempElement = gsl_complex_add(tempElement, gsl_complex_mul_real(gsl_vector_complex_get(resVector, vietaEx[F1mix][M][F2mix][Mp]),levelMixingMatrix[1][F1mix - FStartEx][Fe - FStartEx][M + FEndEx][step] * levelMixingMatrix[1][F2mix - FStartEx][Fe - FStartEx][Mp + FEndEx][step]));
									}
								}*/
								tempElement = gsl_vector_complex_get(resVector, vietaGr[Fg][m][Fgp][mp]);
								temp0 = gsl_complex_mul(tempElement, FF);
								// G_probe = GammaPCalcMinusProbe(Fe, M, Fgp, mp, step);
								// G_probe = GammaPCalcProbe(Fe, M, Fgp, mp, step);

								G_probe = gsl_complex_add(GammaPCalcPlusProbe(Fe, M, Fgp, mp, step), GammaPCalcMinusProbe(Fe, M, Fg, m, step));
								// cout << Fe << " " << M << " " << Fgp << " " << mp << " " << GSL_REAL(G_probe) << " " << detuning_probe << " " << splitProbe(Fe, M, Fgp, mp, step) << endl;
								temp = gsl_complex_mul(G_probe, temp0); //temp0; //
								intensityR += GSL_REAL(temp);
								intensityI += GSL_IMAG(temp);
							}
						}
					}
				}
			}
		}
	}
	// cout << "intensityR: " << intensityR << "; intensityI: " << intensityI << endl;
	return intensityR;
	// return intensityI;
}

//
//funkcija, kas saglabaa failaa bliivumu matricas
void save_density_matrix(double B){
	//matrixFile << "LA LA LA 22";
	//cout << "LA LA AL";
	//ofstream matrixFile;

	// char filename [ FILENAME_MAX ];
	// sprintf(filename, "%s%.2f", "dms/densMatrix-det=m5-B=", B);

	// cout << "dm filename -> " << filename  << endl;
	// matrixFile.open(filename, ios::out | ios::trunc);
	gsl_complex tmp;
	double real, imag;

	for(int Fg = FStartGr; Fg <= FEndGr; Fg++){ //for(int Fg = Fg_d; Fg <= Fg_d; Fg++){ //
		for(int m = -Fg; m <= Fg; m++){

			for(int Fgp = FStartGr; Fgp <= FEndGr; Fgp++){ //for(int Fgp = Fg_d; Fgp <= Fg_d; Fgp++){ //
				for(int mp = -Fgp; mp <= Fgp; mp++){
					tmp = gsl_vector_complex_get(resVector, vietaGr[Fg][m][Fgp][mp]);
					real = GSL_REAL(tmp);
					imag = GSL_IMAG(tmp);
					if(abs(real) < prec) real = 0;
					if(abs(imag) < prec) imag = 0;
					// matrixFile << "(" << real << "," << imag << ") ";
					matrixFile << "" << std::fixed << setprecision(16) << real << "+I*" << imag << "\t";
					matrixFileGr << "" << std::fixed << setprecision(16) << real << "+I*" << imag << "\t";
				}
			}
			matrixFile << endl;
			matrixFileGr << endl;
		}
	}

	matrixFile << endl << endl;
	for(int Fe = FStartEx; Fe <= FEndEx; Fe++){ //for(int Fe = Fe_d; Fe <= Fe_d; Fe++){ //
		for(int M = Fe; M >= -Fe; M--){
			for(int Fep = FStartEx; Fep <= FEndEx; Fep++){ //for(int Fep = Fe_d; Fep <= Fe_d; Fep++){//
				for(int Mp = Fep; Mp >= -Fep; Mp--){
					tmp = gsl_vector_complex_get(resVector, vietaEx[Fe][M][Fep][Mp]);
					real = GSL_REAL(tmp);
					imag = GSL_IMAG(tmp);
					if(abs(real) < prec) real = 0;
					if(abs(imag) < prec) imag = 0;
					matrixFile << "" << std::fixed << setprecision(16) << real << "+I*" << imag << "\t";
					matrixFileEx << "" << std::fixed << setprecision(16) << real << "+I*" << imag << "\t";
				}
			}
			matrixFile << endl;
			matrixFileEx << endl;
		}
	}
	// matrixFile.close();
}

int Fact(int n){
	if(0>n) return -1;
	if(0 == n) return 1;
	else{
		return ( n* Fact(n-1));
	}
}

int lklimit(int J, int m, int mp){
	int k;
	k = max(m-mp,0);
	//cout << "kl_(" << m << "," << mp << ") = " << k;
	//cout << "; (J+m-k) = " << J+m-k << "; (J-k-mp) = " << J-k-mp << "; (k-m+mp) = " << k-m+mp << endl;
	return k;
	}

int uklimit(int J, int m, int mp){
	int k;
	k = min(J+m,J-mp);// retVal = min(retVal,J-mp); retVal = min(retVal,J+mp);
	//cout << "ku_(" << m << "," << mp << ") = " << k;
	//cout << "; (J+m-k) = " << J+m-k << "; (J-k-mp) = " << J-k-mp << "; (k-m+mp) = " << k-m+mp << endl;
	return k;
	}

void fill_rot_matrices(double theta = Pi/2){
	int m, mp, k, llim, ulim, vieta, J = 3;
	double de;
	for(m = -J; m <= J; m += 2){
		for(mp = -J; mp <= J; mp += 2){
			de = 0;
			llim = lklimit(J, m, mp); ulim = uklimit(J, m, mp);
			for(k = llim; k <= ulim; k++){
				de += pow(-1.0, (k-m+mp)/2) * sqrt(Fact((J+m)/2)*Fact((J-m)/2)*Fact((J+mp)/2)*Fact((J-mp)/2)) / (Fact((J+m-k)/2)*Fact(k/2)*Fact((J-k-mp)/2)*Fact((k-m+mp)/2)) * pow(cos(theta/2),J-k+(m-mp)/2) * pow(sin(theta/2),2*k-(m+mp)/2);
			}
			if(abs(de) < prec) de = 0;
			gsl_matrix_set(dMatrixE,(mp+J)/2,(m+J)/2,de);
		}
	}
	J = 1;
	for(m = -J; m <= J; m += 2){
		for(mp = -J; mp <= J; mp += 2){
			de = 0;
			llim = lklimit(J, m, mp); ulim = uklimit(J, m, mp);
			for(k = llim; k <= ulim; k++){
				de += pow(-1.0, (k-m+mp)/2) * sqrt(Fact((J+m)/2)*Fact((J-m)/2)*Fact((J+mp)/2)*Fact((J-mp)/2)) / (Fact((J+m-k)/2)*Fact(k/2)*Fact((J-k-mp)/2)*Fact((k-m+mp)/2)) * pow(cos(theta/2),J-k+(m-mp)/2) * pow(sin(theta/2),2*k-(m+mp)/2);
			}
			if(abs(de) < prec) de = 0;
			gsl_matrix_set(dMatrixG,(mp+J)/2,(m+J)/2,de);
		}
	}

}

void rotate_dens_matrix(int type = 0, int J = Fg_d, double theta = Pi/2){
	int m, mp, k, llim, ulim, vieta;
	double de; gsl_complex tv;
	gsl_matrix_complex * dMatrix = gsl_matrix_complex_calloc(2*J+1, 2*J+1);
	gsl_matrix_complex * tmp = gsl_matrix_complex_calloc(2*J+1, 2*J+1);
	rotated = gsl_matrix_complex_calloc(2*J+1, 2*J+1);
	for(m = -J; m <= J; m++){
		for(mp = -J; mp <= J; mp++){
			de = 0;
			llim = lklimit(J, m, mp); ulim = uklimit(J, m, mp);
			for(k = llim; k <= ulim; k++){
				de += pow(-1.0, k-m+mp) * sqrt(Fact(J+m)*Fact(J-m)*Fact(J+mp)*Fact(J-mp)) / (Fact(J+m-k)*Fact(k)*Fact(J-k-mp)*Fact(k-m+mp)) * pow(cos(theta/2),2*J-2*k+m-mp) * pow(sin(theta/2),2*k-m+mp);
			}
			if(abs(de) < prec) de = 0;
			gsl_matrix_complex_set(dMatrix,mp+J,m+J,gsl_complex_rect(de,0));
			//cout << "de_(" << mp << "," << m << ") = " << de << endl;
			if(type == 0) vieta = vietaGr[J][m][J][mp];
			else vieta = vietaEx[J][m][J][mp];
			//if(abs(m-mp) == 0 && abs(m) == 2)
			 gsl_matrix_complex_set(rotated,mp+J,m+J,gsl_vector_complex_get(resVector, vieta));
			//cout << " rho(" << m << "," << mp << ") = " << GSL_REAL(gsl_vector_complex_get(resVector, vieta)) << " + " << GSL_IMAG(gsl_vector_complex_get(resVector, vieta)) << " * i" << endl;
		}
		//cout << endl;
	}
	//cout << endl;
	gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), dMatrix, rotated, gsl_complex_rect(0,0), tmp);
	gsl_blas_zgemm (CblasNoTrans, CblasTrans, gsl_complex_rect(1,0), tmp, dMatrix, gsl_complex_rect(0,0), rotated);
	/*for(m = -J; m <= J; m++){
		for(mp = -J; mp <= J; mp++){
			tv = gsl_matrix_complex_get(rotated, m+J,mp+J);
			if(gsl_complex_abs(tv) < prec) tv = gsl_complex_rect(0,0);
			if(m == mp)
			 cout << " rho(" << m << "," << mp << ") = " << GSL_REAL(tv) << " + " << GSL_IMAG(tv) << " * i" << endl;;
		}
		cout << endl;
	}
	cout << endl;*/
	
}

// funkcija savāc un vidējo Doplera profilā blīvuma matricas
void collect_dens_matrix(int Bi, double weight){
    //cout << "svars = " << weight << endl;
    int Fg, Fgp, Fe, Fep, m, mp;
    dcomplex el; gsl_complex tv;

    for(Fg = FStartGr; Fg <= FEndGr; Fg++){
        for(m = -Fg; m <= Fg; m++){
            for(Fgp = FStartGr; Fgp <= FEndGr; Fgp++){
                for(mp = -Fgp; mp <= Fgp; mp++){
                    tv = gsl_vector_complex_get(resVector, vietaGr[Fg][m][Fgp][mp]);
                    rhoG[vietaMxGr[Fg][m]][vietaMxGr[Fgp][mp]] += GSL_REAL(tv)  + dc_imagunit * GSL_IMAG(tv) * weight;
                }
            }
        }
    }
    for(Fe = FStartEx; Fe <= FEndEx; Fe++){
        for(m = -Fe; m <= Fe; m++){
            for(Fep = FStartEx; Fep <= FEndEx; Fep++){
                for(mp = -Fep; mp <= Fep; mp++){
                    tv = gsl_vector_complex_get(resVector, vietaEx[Fe][m][Fep][mp]);
                    rhoE[vietaMxEx[Fe][m]][vietaMxEx[Fep][mp]] += GSL_REAL(tv)  + dc_imagunit * GSL_IMAG(tv) * weight;
                }
            }
        }
    }
}


// funkcija, kas pārveido blīvuma matricu |J,mj> bāzē
void base_f_to_j(int level = 0, int Bi = 0){
    int FStart, FEnd, F, Fp, M, Mp, J, m, mp;
    if(level == 0){
        FStart = FStartGr; FEnd = FEndGr;
    }
    else{
        FStart = FStartEx; FEnd = FEndEx;
    }
    for(m = -JJ[level]; m <= JJ[level]; m++){
        for(mp = - JJ[level]; mp <= JJ[level]; mp){
            for(F = FStart; F <= FEnd; F++){
                for(M = - F; M <= F; M++){
                    for(Fp = FStart; F <= FEnd; F++){
                        for(Mp = - Fp; Mp <= Fp; Mp++){
                        }
                    }
                }
            }
        }
    }
}

// funkcija izdrukā / saglabā vidējotās blīvuma matricas
void treat_collected(double B){
   	gsl_complex tmp; double sumEx = 0, sumGr = 0; dcomplex SM, moment, tv;
   	int Fg, Fgp, Fe, Fep, m, mp;
   	momFile << B << " ";
   	ofstream mxFile;
    char filename [ FILENAME_MAX ];
    for(Fe = FStartEx; Fe <= FEndEx; Fe++){
    	sprintf(filename, "%s%i%s%.0f", "pol/mxE/Fe=", Fe, "/matrix-B=", B);
    	mxFile.open(filename);
        moment = 0;
        for(m = -Fe; m <= Fe; m++){
            for(mp = -Fe; mp <= Fe; mp++){
                //cout << "rho(Fe = " << Fe << ", m = " << m << ", mp = " << mp << ") = " << rhoE[vietaMxEx[Fe][m]][vietaMxEx[Fe][mp]] << endl;
                tv = rhoE[vietaMxEx[Fe][m]][vietaMxEx[Fe][m]];
                if(m == mp) sumEx += real(tv);
                mxFile << real(tv) << " + I * " << imag(tv) << " ";
            }
            mxFile << endl;
        }
        mxFile.close();
        for(m = -Fe; m <= Fe; m++){
            moment += (double)m * rhoE[vietaMxEx[Fe][m]][vietaMxEx[Fe][m]];
        }
        cout << "Fe = " << Fe << "; moments = " << moment/sumEx << endl;
        momFile << real(moment) << " ";
        SM += moment;
    }
    cout << "Exc moment = " << SM << endl;
    momFile << real(SM) << " " << sumEx << " ";
    SM = 0;

    moment = 0;
    for(Fg = FStartGr; Fg <= FEndGr; Fg++){
        dcomplex moment = 0;
        for(m = -Fg; m <= Fg; m++){
            sumGr += real(rhoG[vietaMxGr[Fg][m]][vietaMxGr[Fg][m]]);
        }
        for(m = -Fg; m <= Fg; m++){
            moment += (double)m * rhoG[vietaMxGr[Fg][m]][vietaMxGr[Fg][m]];
        }
        cout << "Fg = " << Fg << "; moments = " << moment/sumGr << endl;
        momFile << real(moment) << " ";
        SM += moment;
    }
    cout << "Gr moment = " << SM << endl;
    momFile << real(SM) << " " << sumGr << endl;
    SM = 0;
    cout << "sumGr = " << sumGr << "; sumEx = " << sumEx << endl;

    int i, j, k;
    for(i = 0; i <= Edim; i++){
        for(j = 0; j<= Edim; j++){
            rhoE[i][j] = 0;
        }
    }
    for(i = 0; i <= Gdim; i++){
        for(j = 0; j<= Gdim; j++){
            rhoG[i][j] = 0;
        }
    }
    //set_dms_zero();

}

double count_dark_states(int m, int mp){
	gsl_complex tmp;
	double abso = 0;
	tmp = gsl_complex_div(gsl_vector_complex_get(resVector, vietaGr[Fg_d][m][Fg_d][mp]), gsl_complex_sqrt(gsl_complex_mul(gsl_vector_complex_get(resVector, vietaGr[Fg_d][m][Fg_d][m]), gsl_vector_complex_get(resVector, vietaGr[Fg_d][mp][Fg_d][mp]))));
	abso += gsl_complex_abs(tmp);
	return abso;
}

double count_bright_states(int mc = 0){
	gsl_complex tmp;
	double abso = 0;
	for(int Fg = FStartGr; Fg <= FEndGr; Fg++){
		for(int m = -Fg; m <= Fg; m++){
			for(int Fgp = FStartGr; Fgp <= FEndGr; Fgp++){
				for(int mp = -Fgp; mp <= Fgp; mp++){
					if(abs(m-mp) == 0 && m == mc && Fg == Fgp && Fg == Fg_d){
						tmp = gsl_vector_complex_get(resVector, vietaGr[Fg][m][Fgp][mp]);
						abso += gsl_complex_abs(tmp);
					}
				}
			}
		}
	}
	return abso;
}

double count_bright_parasites(){
	gsl_complex tmp;
	double abso = 0;
	for(int Fg = FStartGr; Fg <= FEndGr; Fg++){
		for(int m = -Fg; m <= Fg; m++){
			for(int Fgp = FStartGr; Fgp <= FEndGr; Fgp++){
				for(int mp = -Fgp; mp <= Fgp; mp++){
					if(abs(m-mp) == 2 && (m == 0 || mp == 0) && (Fgp == Fg_d || Fg == Fg_d)){
						tmp = gsl_vector_complex_get(resVector, vietaGr[Fg][m][Fgp][mp]);
						abso += gsl_complex_abs(tmp);
					}
				}
			}
		}
	}
	return abso;
}

bool FileExists(string strFilename) { 
  struct stat stFileInfo; 
  bool blnReturn; 
  int intStat; 

  // Attempt to get the file attributes 
  intStat = stat(strFilename.c_str(),&stFileInfo); 
  if(intStat == 0) { 
    // We were able to get the file attributes 
    // so the file obviously exists. 
    blnReturn = true; 
  } else { 
    // We were not able to get the file attributes. 
    // This may mean that we don't have permission to 
    // access the folder which contains this file. If you 
    // need to do that level of checking, lookup the 
    // return values of stat which will give you 
    // more details on why stat failed. 
    blnReturn = false; 
  } 
   
  return(blnReturn); 
}

