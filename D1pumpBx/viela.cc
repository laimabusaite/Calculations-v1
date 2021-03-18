/* liimenu parametri
int FStartGr, FEndGr, FDeltaGr, FStartEx, FEndEx, FDeltaEx;
*/
/*vielas elementa raksturlielumi
double ik; //kodola spins;
double S;//elektrona spins
double gi;//kodola spins*/

/*liimenju raksturlielumi
double LL[4];//liimenja elektonu orbitaalais moments
double JJ[4];//liimenja pilnais lenkiskais moments
double aa[4];//liimenja siikstruktuuras koeficients A;
double bb[4];//liimenja siikstruktuuras koeficients B;
double gj[4];//Landee faktors liimenim*/

/*//133Rb 6(2D3/2)//vielas raksturlielumi;/
	    double ik=3.5;
    	double gi=-0.00039885395;
  //pamatliimenja 5(2S1/2) koeficienti un ierosinaata liimenja 5(2P3/2) koeficienti
    	double LL[2] = { 0.0, 1 };
    	double JJ[2] = { 0.5, 1.5 };
	    double gj[2] = { 2.00233113, 1.3341 };//2.00233113;//=( 1 + ((Jg*(Jg+1)+S*(S+1)-Lgr*(Lgr+1))/(2*Jg*(Jg+1))) 
    	double aa[2] = { 2298.1579425, 16.605 };//1011.910813; //(MHz)
    	//double aa[2] = { 100000000011.9, 25000000000.009 };//1011.910813; //(MHz)
    	double bb[2] = { 0.0, -0.15 };//(MHz)
		double masa = 0.133/(6.02 * pow(10,23.0));
        double lambda = 894.59295986 * pow(10.0, -9.0);
		double Gamma_n = 4.575; //MHz
/*   
//85Rb D1;liimenji:2,3->2,3; 5(2S1/2)->5(2P1/2)//vielas raksturlielumi;/
	    double ik=2.5;
    	double gi=-0.00029364;
  //pamatliimenja 5(2S1/2) koeficienti un ierosinaata liimenja 5(2P3/2) koeficienti
    	double LL[2] = { 0.0, 1 };
    	double JJ[2] = { 0.5, 0.5 };
	    double gj[2] = { 2.00233113, 0.666 };//2.00233113;//=( 1 + ((Jg*(Jg+1)+S*(S+1)-Lgr*(Lgr+1))/(2*Jg*(Jg+1))) 
    	double aa[2] = { 1011.910813, 120.527 };//1011.910813; //(MHz)
    	//double aa[2] = { 100000000011.9, 25000000000.009 };//1011.910813; //(MHz)
    	double bb[2] = { 0.0, 0.0 };//(MHz)
		double masa = 1.409993199 * pow(10,-25.0);
        double lambda = 794.767282 * pow(10.0, -9.0);
		double Gamma_n = 5.75; //MHz
		// const char baseFolder[] = "D1-Rb85-noDoppler-HanleGeometry-LinObs/";
    // const char baseFolder[] = "D1-Rb85-fine-doppler/";
    // const char baseFolder[] = "D1-Rb85-no-doppler/";
    const char baseFolder[] = "Pump_D1-Rb85_";
    // const char baseFolder[] = "D1-Rb85-Doppler2.75Sigma/";
    // const char baseFolder[] = "D1-Rb85-DopplerSigmaB1420-1560/";
    // const char baseFolder[] = "D1-Rb85-LaserDetuning w Doppler/";
    // const char baseFolder[] = "D1-Rb85-DopplerSigma1/";

		const char baseFilename[] = "LIF_Rb85-D1-Rabi=";
    // const char dmfldname[] = "dms";
/*
///85Rb D2;liimenji:2,3->1,2,3,4; 5(2S1/2)->5(2P3/2)//vielas raksturlielumi;/
	    double ik=2.5;
    	double gi=-0.00029364;
  //pamatliimenja 5(2S1/2) koeficienti un ierosinaata liimenja 5(2P3/2) koeficienti
    	double LL[2] = { 0.0, 1 };
    	double JJ[2] = { 0.5, 1.5 };
	    double gj[2] = { 2.00233113, 1.3362 };//2.00233113;//=( 1 + ((Jg*(Jg+1)+S*(S+1)-Lgr*(Lgr+1))/(2*Jg*(Jg+1))) 
    	double aa[2] = { 1011.910813, 25.002 };//1011.910813; //(MHz)
    	//double aa[2] = { 100000000011.9, 25000000000.009 };//1011.910813; //(MHz)
    	double bb[2] = { 0.0, 25.79 };//(MHz)
		double masa = 1.409993199 * pow(10,-25.0);
        double lambda = 780.241368271 * pow(10.0, -9.0);
		double Gamma_n = 6.0666; //MHz
		const char baseFolder[] = "Pump_D2-Rb85_";
		const char baseFilename[] = "LIF_Rb85-D2-Rabi=";
   
      /*      
85 Rb	      ------------ Fe=4
ierosinaatie  ------------ Fe=3
limeji        ------------ Fe=2
 5(2P3/2)     ------------ Fe=1

5(2P1/2)      ------------ Fg=3
              ------------ Fg=2


pamatlimenji  ------------ Fg=3
              ------------ Fg=2

      *//*
//87Rb D2;liimenji:1,2->0,1,2,3; 5(2S1/2)->3(2P3/2)///vielas raksturlielumi;/
	    double ik=1.5;
    	double gi=-0.0009951414; 
  //pamatliimenja 5(2S1/2) koeficienti un ierosinaata liimenja 5(2P3/2) koeficienti
    	double LL[2] = { 0.0, 1 };
    	double JJ[2] = { 0.5, 1.5 };
	    double gj[2] = { 2.00233113, 1.3362 };//2.00233113;//=( 1 + ((Jg*(Jg+1)+S*(S+1)-Lgr*(Lgr+1))/(2*Jg*(Jg+1))) 
    	double aa[2] = { 3417.341305452145, 84.7185 };//1011.910813; //(MHz)
    	double bb[2] = { 0.0, 12.4965 };//(MHz)
		double masa = 1.443160648 * pow(10,-25.0);
        double lambda = 780.241209686 * pow(10.0, -9.0); //sho izmanto Doplera platuma aprekinahanai (mainot pa 10nm aprekinos nevajadzetu but butiskam izmainam - ta saka Linards)
		double Gamma_n = 6.0666; //MHz
		const char baseFolder[] = "D2-Rb87-detuning-off-ground-on-excited-off-2013/"; //folderis kura saglabas datus
		const char baseFilename[] = "Rb87-D2-Rabi="; //failu nosaukums ar kadu tiks saglabats
      
/*      
87 Rb	      ------------ Fe=3
ierosinaatie  ------------ Fe=2
limeji        ------------ Fe=1
              ------------ Fe=0


pamatlimenji  ------------ Fg=2
              ------------ Fg=1

//133Cs D2;liimenji:3,4->2,3,4,5; 6(2S1/2)->6(2P3/2)///vielas raksturlielumi/
	    double ik=3.5;
    	double gi=-0.00039885395;
  //pamatliimenja 6(2S1/2) koeficienti un ierosinaata liimenja 6(2P3/2) koeficienti
    	double LL[2] = { 0.0, 1 };
    	double JJ[2] = { 0.5, 1.5 };
	    double gj[2] = { 2.00254032, 1.334 };//2.00233113;//=( 1 + ((Jg*(Jg+1)+S*(S+1)-Lgr*(Lgr+1))/(2*Jg*(Jg+1))) 
    	double aa[2] = { 2298.1579425, 50.288827 };//1011.910813; //(MHz)
    	double bb[2] = { 0.0, -0.4934 };//(MHz)
		double masa = 2.206211106 * pow(10,-25.0);
        double lambda = 852.34727582 * pow(10.0, -9.0);
		double Gamma_n = 5.234; //MHz
        const char baseFolder[] = "Pump_D2-Cs133_";
        const char baseFilename[] = "LIF_Cs133-D2-Rabi=";
      
      /*
133 Cs	      ------------ Fe=5
ierosinaatie  ------------ Fe=4
limeji        ------------ Fe=3
              ------------ Fe=2


pamatlimenji  ------------ Fg=4
              ------------ Fg=3

//87Rb D1;liimenji:1,2->1,2 5(2S1/2)->3(2P1/2)///vielas raksturlielumi/
	    double ik=1.5;
    	double gi=-0.0009951414; 
  //pamatliimenja 5(2S1/2) koeficienti un ierosinaata liimenja 5(2P3/2) koeficienti
    	double LL[2] = { 0.0, 1 };
    	double JJ[2] = { 0.5, 0.5 };
	    double gj[2] = { 2.00233113, 0.666 };//2.00233113;//=( 1 + ((Jg*(Jg+1)+S*(S+1)-Lgr*(Lgr+1))/(2*Jg*(Jg+1))) 
    	double aa[2] = { 3417.341305452145, 407.24 };//1011.910813; //(MHz)
    	//double aa[2] = { 34170000.34130642, 4080000.328 };//1011.910813; //(MHz)
    	double bb[2] = { 0.0, 0.0 };//(MHz)
		double masa = 1.443160648 * pow(10,-25.0);
        double lambda = 794.978851016 * pow(10.0, -9.0);
		double Gamma_n = 5.75; //MHz
		// const char baseFolder[] = "D1-Rb87-noDoppler-HanleGeometry-LinObs/";
    // const char baseFolder[] = "D1-Rb87/";
    const char baseFolder[] = "D1-Rb87-DopplerSigma/";
    // const char baseFolder[] = "D1-Rb87-no-doppler/";
		const char baseFilename[] = "Rb87-D1-Rabi=";

/*
		 
      
//133Cs; D1 liimenji:3,4->3,4; 6(2S1/2)->6(2P1/2)///vielas raksturlielumi/*/
	    double ik=3.5;
    	double gi=-0.00039885395;
    //pamatliimenja koeficienti
    	double LL[2]= { 0.0, 1 };
    	double JJ[2]= { 0.5, 0.5 };
    	double gj[2]= { ( 1 + ((JJ[0]*(JJ[0]+1)+0.5*(0.5+1)-LL[0]*(LL[0]+1))/(2*JJ[1]*(JJ[1]+1))) ), ( 1 + ((JJ[1]*(JJ[1]+1)+0.5*(0.5+1)-LL[1]*(LL[1]+1))/(2*JJ[1]*(JJ[1]+1))) ) };//=( 1 + ((Jg*(Jg+1)+S*(S+1)-Lg*(Lg+1))/(2*Jg*(Jg+1))=2; 
    	double aa[2]= { 2298, 292 }; //(MHz)
    	double bb[2]= { 0.0, 0.0 };//(MHz)
		double masa = 0.133/(6.02 * pow(10,23.0));
        double lambda = 894.59295986 * pow(10.0, -9.0);
		double Gamma_n = 4.575; //MHz
		// const char baseFolder[] = "Pump_D1-Cs133-DopplerON-B300/";
        const char baseFolder[] = "Pump_D1-Cs133_";
		const char baseFilename[] = "LIF_Cs133-D1-Rabi=";
     
    //ierosinaata liimenja koeficienti
    	/*LL[1]=1;
    	JJ[1]=0.5;//=1/2
    	gj[1]=( 1 + ((JJ[1]*(JJ[1]+1)+0.5*(0.5+1)-LL[1]*(LL[1]+1))/(2*JJ[1]*(JJ[1]+1))) );
    	aa[1]=292;//(MHz)
    	bb[1]=0;//-0.38;//(MHz)*/

/*   
//39K D1;liimenji:1,2->1,2; 4(2S1/2)->4(2P1/2)//vielas raksturlielumi;/
        double ik=1.5;//Tiecke(2011) nuclear spin = 3/2
        double gi=-0.00014193489;//vecais->Rb85D1=0.00029364;// nuclear Lande factor - paliek nemainīgs? precizet vajag
  //pamatliimenja 4(2S1/2) koeficienti un ierosinaata liimenja 4(2P1/2) koeficienti
        double LL[2] = { 0.0, 1 };//pamata L vertiba un ierosinata stavokla L vertiba
        double JJ[2] = { 0.5, 0.5 };//pamata J vertiba un ierosin J vert
        double gj[2] = { 2.00229421, 0.666 };//2.00233113;//=( 1 + ((Jg*(Jg+1)+S*(S+1)-Lgr*(Lgr+1))/(2*Jg*(Jg+1))) (nevajadzetu mainities)pamat un ierosin parbaudit
        double aa[2] = { 230.8598601, 27.775 };//1011.910813; //(MHz) 
        //double aa[2] = { 100000000011.9, 25000000000.009 };//1011.910813; //(MHz)
        double bb[2] = { 0.0, 0.0 };//(MHz)
        double masa = 6.49242489 * pow(10,-26.0);//kg????
        double lambda = 770.108385049 * pow(10.0, -9.0);//m???
        double Gamma_n = 5.956; //MHz (Natural linewidth)
        const char baseFolder[] = "Pump_D1-K39_";
        const char baseFilename[] = "LIF_K39-D1-Rabi=";
/*

/*   
//39K D2;liimenji:1,2->0,1,2,3; 4(2S1/2)->4(2P1/2)//vielas raksturlielumi;/
        double ik=1.5;//Tiecke(2011) nuclear spin = 3/2
        double gi=-0.00014193489;//vecais->Rb85D1=0.00029364;// nuclear Lande factor - paliek nemainīgs? precizet vajag
  //pamatliimenja 4(2S1/2) koeficienti un ierosinaata liimenja 4(2P1/2) koeficienti
        double LL[2] = { 0.0, 1 };//pamata L vertiba un ierosinata stavokla L vertiba
        double JJ[2] = { 0.5, 1.5 };//pamata J vertiba un ierosin J vert
        double gj[2] = { 2.00229421, 1.333 };//2.00233113;//=( 1 + ((Jg*(Jg+1)+S*(S+1)-Lgr*(Lgr+1))/(2*Jg*(Jg+1))) (nevajadzetu mainities)pamat un ierosin parbaudit
        double aa[2] = { 230.8598601, 6.093 };//1011.910813; //(MHz) 
        //double aa[2] = { 100000000011.9, 25000000000.009 };//1011.910813; //(MHz)
        double bb[2] = { 0.0, 2.786 };//(MHz)
        double masa = 6.49242489 * pow(10,-26.0);//kg????
        double lambda = 766.700921822 * pow(10.0, -9.0);//m???
        double Gamma_n = 6.035; //MHz (Natural linewidth)
        const char baseFolder[] = "Pump_D2-K39_/";
        const char baseFilename[] = "LIF_K39-D2-Rabi=";
/**/
		
int FStartGr = int(abs(ik - JJ[0]));
int FEndGr = int(abs(ik + JJ[0]));
int FDeltaGr = int(FEndGr - FStartGr + 1);

int FStartEx = int(abs(ik - JJ[1]));
int FEndEx = int(abs(ik + JJ[1]));
int FDeltaEx = int(FEndEx - FStartEx + 1);

double splitsGr[10];
double splitsEx[10];

double GammaP[10][10];
