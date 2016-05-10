#include "MC_VEGA_SER_001.hh"

// Parallel code to simulate neutron scattering with U235 and Gd152
// by using Monte Carlo code
// developed by Abhijit Bhattacharyya
// Nuclear Physics Division, Barc, Mumbai - 400 085, INDIA

// Find maximum of three double numbers
//
double max3(double a, double b, double c)
{
  double temp = max2(a, c);
  double maximum = (temp > b) ? temp : b;
  return maximum;
}

// Find minimum of three double numbers
//
double min3(double a, double b, double c)
{
  double temp = min2(a, c);
  double minimum = (temp < b) ? temp : b;
  return minimum;
}

// returns the gradient for linear Interpolation routine
//
double funcinterp(double x0, double x1, double y0, double y1)
{
  return (y1 - y0) / (x1 - x0);
}

//  Linear Interpolation routine computed upto second order of accuracy
//
double interpol(double eval, double eng[maxf], double xs[maxf], long last)
{
  double x0, x1, x2, x3;
  double y0, y1, y2, y3;
  double term1, term2, term3, term4, xsecinterp;
  int i, j, keel = 0;

  if (eval <= eng[0]) 		return xs[0];
  if (eval >= eng[last]) 	return xs[last];
  if (eval < eng[last]){
    for (i = 1; i < last; i++){
      if (abs(eval - eng[i]) <= 1.0e-5) return xs[i];
      if (eval > eng[i]) {
	if (eval <= eng[i+1]){
	  keel = i;
	}
      }
    }
    x0 = eng[keel - 1];    	y0 = xs[keel - 1];
    x1 = eng[keel]; 		y1 = xs[keel];
    x2 = eng[keel + 1];		y2 = xs[keel + 1];
    x3 = eng[keel + 2];		y3 = xs[keel + 2];
  }

  term1 = y0;
  term2 = (eval - x0) * funcinterp(x0, x1, y0, y1);
  term3 = (eval - x0) * (eval - x1) * (funcinterp(x1, x2, y1, y2) - funcinterp(x0, x1, y0, y1)) / (x2 - x0);
  term4 = (eval - x0) * (eval - x1) * (eval - x2) * (1.0 / (x3 - x0)) * ((funcinterp(x2, x3, y2, y3) - funcinterp(x1, x2, y1, y2)) / (x3 - x1) - (funcinterp(x1, x2, y1, y2) - funcinterp(x0, x1, y0, y1)) / (x2 - x0));

  xsecinterp = term1 + term2 + term3 + term4;

  return xsecinterp;
}



// The code for reflection
//
void reflection(double x0, double y0, double z0, double xc, double yc, double zc, double d2walk)
{
  double ntersecx, ntersecy, ntersecz, newx, newy, newz, dnew, Rf2, dx, dy, dz;
  double dircosl, dircosm, dircosn, xlast, ylast, zlast;
  double tparam, t1, t2, lineA, lineB, lineC;

  lineA = x0 * x0 + y0 * y0 + z0 * z0 - Rreflector * Rreflector;
  lineB = 2.0 * (x0 * (xc - 1.0) + y0 * (yc - 1.0) + z0 * (zc - 1.0));
  lineC = ((x0 - xc) * (x0 - xc)) + ((y0 - yc) * (y0 - yc)) + ((z0 - zc) * (z0 - zc));
  t1 = (-lineA + pow((lineB * lineB - 4.0 * lineA * lineC), 0.5)) / (2.0 * lineC);
  t2 = (-lineA - pow((lineB * lineB - 4.0 * lineA * lineC), 0.5)) / (2.0 * lineC);
  if (t1 >= 0 && t1 <= 1) tparam = t1;
  if (t2 >= 0 && t2 <= 1) tparam = t2;
  ntersecx = x0 * (1.0 - tparam) + tparam * xc;
  ntersecy = y0 * (1.0 - tparam) + tparam * yc;
  ntersecz = z0 * (1.0 - tparam) + tparam * zc;
  dx = ntersecx - x0; dy = ntersecy - y0;  dz = ntersecz - z0;
  Rf2 = ntersecx*ntersecx + ntersecy*ntersecy + ntersecz*ntersecz;
  newx = ntersecx - x0 - (2.0 * ntersecx/Rf2) * (ntersecx * dx + ntersecy * dy + ntersecz * dz);
  newy = ntersecy - y0 - (2.0 * ntersecy/Rf2) * (ntersecx * dx + ntersecy * dy + ntersecz * dz);
  newz = ntersecz - z0 - (2.0 * ntersecz/Rf2) * (ntersecx * dx + ntersecy * dy + ntersecz * dz);
  dnew = pow(((newx - ntersecx) * (newx * ntersecx) + (newy - ntersecy) * (newy - ntersecy) + (newz - ntersecz) * (newz - ntersecz)), 0.5);
  dircosl = (newx - ntersecx) / dnew;
  dircosm = (newy - ntersecy) / dnew;
  dircosn = (newz - ntersecz) / dnew;
  xc = ntersecx + dircosl * d2walk;
  yc = ntersecy + dircosm * d2walk;
  zc = ntersecz + dircosn * d2walk;
  return;
}


int whichnucleus(double whichmod, double sigD2O, double sigC12, double sigBnat, double sigCH2)
{
  //whichmod = 15:D2O  16:Graphite   19:B  21:Plastic
  //whichmat = 17:Gd  18:U
  double tempo, isD = 0, isO = 0,isC = 0, isH = 0;
  int whatATOM = 0;

  if (whichmod == 16){	// GRAPHITE
    return whichmod;
  } else if (whichmod == 19){	// Natural Boron
    return whichmod;
  } else if (whichmod == 15){	// D2O
    isD = nD_D2O/sigD2O;      isO = nO_D2O/sigD2O;     // isO is around 0.01, 0.02,.... isD is around 0.3, 0.2,...
    tempo =static_cast<double>(rand()) / static_cast<double>(RAND_MAX + 1.0);
    while (tempo > isD && tempo > isO) {
	srand(time(NULL));
	tempo = static_cast<double>(rand()) / static_cast<double>(RAND_MAX + 1.0);
    }
    if (isD >= tempo && isO < tempo) whatATOM = 0;	// collides with D
    if (isO >= tempo && isD < tempo) whatATOM = 11;	// collides with O
    do {
	tempo = 0.0;
	tempo = static_cast<double>(rand()) / static_cast<double>(RAND_MAX + 1.0);
	whatATOM = (tempo > 0.33333) ? 0 : 11;
    }while ((isO > tempo && isD > tempo) || (isO < tempo && isD < tempo));
    return whatATOM;
  } else if (whichmod == 21){	// Plastic (C2H4)n
    isH = nH1Pl / sigCH2;     isC = nC12Pl / sigCH2;
    tempo =static_cast<double>(rand()) / static_cast<double>(RAND_MAX + 1.0);
    while (tempo > isC && tempo > isH) {
	srand(time(NULL));
	tempo = static_cast<double>(rand()) / static_cast<double>(RAND_MAX + 1.0);
    }
    if (isH >= tempo && isC < tempo) whatATOM = 20;	// collides with H
    if (isC >= tempo && isH < tempo) whatATOM = 4;	// collides with C
    do {
	tempo = 0.0;
	tempo = static_cast<double>(rand()) / static_cast<double>(RAND_MAX + 1.0);
	whatATOM = (tempo >= 0.5) ? 20 : 4;
    }while ((isO > tempo && isD > tempo) || (isO < tempo && isD < tempo));
    return whatATOM;
  }
  // return list: 0:D  4:C  11:O  16:G  19:B  20:H
}


// Scatter each neutron till (a) it is lost out of reflector (b) absorbed (c0 thermal
//
void scatter(struct particles *neu, int whichmod, int whichmat, int long id, bool rpt){
    double ener, sigE, sigI, sigR, sigF, sign2n;
    double sigH1, sigH2, sigO16, sigD2O;
    double sigBe9, sigC12, sigB10, sigB11;
    double sigU235, sigU238, sU235f, sU235r, sU238f, sU238r, sigUnat, sigGdnat, sigBnat;
    double sigGd154, sigGd155, sigGd156, sigGd157, sigGd158, sigGd160, sigCH2;
    double sigmod, sigmat, sigreflector;
    double modmass, efactorLAB;
    double signow, mfpnow, d = 500.0, displace, isD, isO;
    double xc, yc, zc, rc, delE, tmpU, tmpV, tmpW;
    double coslamCM;
    double efactorLABD, efactorLABO, efactorLABU, efactorLABGd;
    double csU235, csGd, csD, csO, CSMIN, CSMAX, RANGE, TOTALCS;
    double tempo, dpassed, dtowalk;
    double x0, y0, z0, r0;
    int i, j, whatATOM = 0, ibin;

    srand(time(NULL));
    ener = neu->En;	// Take the current energy of the present neutron

    //******************************* 1H1 **************************************
    sigE = interpol(ener, enH1E, xH1E, maxH1E);
    sigI = interpol(ener, enH1I, xH1I, maxH1I);
    sigH1 = sigE + sigI;

    //******************************* 1H2 **************************************
    sigE = interpol(ener, enH2E, xH2E, maxH2E);
    sigI = interpol(ener, enH2I, xH2I, maxH2I);
    sigH2 = sigE + sigI;

    //******************************* 8O16 **************************************
    sigE = interpol(ener, enO16E, xO16E, maxO16E);
    sigI = interpol(ener, enO16I, xO16I, maxO16I);
    sigR = interpol(ener, enO16R, xO16R, maxO16R);
    sigO16 = sigE + sigI + sigR;

    //******************************** D2O **************************************
    sigD2O = (nD_D2O * sigH2 + nO_D2O * sigO16);			// total SIGMA
    //isD = nD_D2O/sigD2O;      isO = nO_D2O/sigD2O;

    //********************************** B10 ************************************
    sigE = interpol(ener, enB10E, xB10E, maxB10E);
    sigI = interpol(ener, enB10I, xB10I, maxB10I);
    sigR = interpol(ener, enB10R, xB10R, maxB10R);
    sigB10 = sigE + sigI + sigR;

    //********************************** B11 ************************************
    sigE = interpol(ener, enB11E, xB11E, maxB11E);
    sigI = interpol(ener, enB11I, xB11I, maxB11I);
    sigR = interpol(ener, enB11R, xB11R, maxB11R);
    sigB11 = sigE + sigI + sigR;

    //******************************** Be9 **************************************
    sigE = interpol(ener, enBe9E, xBe9E, maxBe9E);
    sigI = interpol(ener, enBe9I, xBe9I, maxBe9I);
    sigR = interpol(ener, enBe9R, xBe9R, maxBe9R);
    sign2n = interpol(ener, enBe9N2N, xBe9N2N, maxBe9N2N);
    sigBe9 = sigE + sigI + sigR + sign2n;

    //***************************** GRAPHITE ***********************************
    sigE = interpol(ener, enC12E, xC12E, maxC12E);
    sigI = interpol(ener, enC12I, xC12I, maxC12I);
    sigR = interpol(ener, enC12R, xC12R, maxC12R);
    sigC12 = sigE + sigI + sigR;

    //***************************** U235 *************************************
    sigE = interpol(ener, enU235E, xU235E, maxU235E);
    sigF = interpol(ener, enU235F, xU235F, maxU235F);
    sigI = interpol(ener, enU235I, xU235I, maxU235I);
    sigR = interpol(ener, enU235R, xU235R, maxU235R);
    sU235f = sigF; sU235r = sigR;
    sigU235 = sigE + sigF + sigI + sigR;

    //***************************** U238 *************************************
    sigE = interpol(ener, enU238E, xU238E, maxU238E);
    sigF = interpol(ener, enU238F, xU238F, maxU238F);
    sigI = interpol(ener, enU238I, xU238I, maxU238I);
    sigR = interpol(ener, enU238R, xU238R, maxU238R);
    sU238f = sigF; sU238r = sigR;
    sigU238 = sigE + sigF + sigI + sigR;

    //***************************** NATURAL U *********************************
    sigUnat = (nU235 * sigU235 + nU238 * sigU238);

    //***************************** Gd154 **************************************
    sigE = interpol(ener, enGd154E, xGd154E, maxGd154E);
    sigI = interpol(ener, enGd154I, xGd154I, maxGd154I);
    sigR = interpol(ener, enGd154R, xGd154R, maxGd154R);
    sigGd154 = sigE + sigI + sigR;

    //***************************** Gd155 ***************************************
    sigE = interpol(ener, enGd155E, xGd155E, maxGd155E);
    sigI = interpol(ener, enGd155I, xGd155I, maxGd155I);
    sigR = interpol(ener, enGd155R, xGd155R, maxGd155R);
    sigGd155 = sigE + sigI + sigR;

    //***************************** Gd156 ***************************************
    sigE = interpol(ener, enGd156E, xGd156E, maxGd156E);
    sigI = interpol(ener, enGd156I, xGd156I, maxGd156I);
    sigR = interpol(ener, enGd156R, xGd156R, maxGd156R);
    sigGd156 = sigE + sigI + sigR;

    //***************************** Gd157 **************************************
    sigE = interpol(ener, enGd157E, xGd157E, maxGd157E);
    sigI = interpol(ener, enGd157I, xGd157I, maxGd157I);
    sigR = interpol(ener, enGd157R, xGd157R, maxGd157R);
    sigGd157 = sigE + sigI + sigR;

    //****************************** Gd158 **************************************
    sigE = interpol(ener, enGd158E, xGd158E, maxGd158E);
    sigI = interpol(ener, enGd158I, xGd158I, maxGd158I);
    sigR = interpol(ener, enGd158R, xGd158R, maxGd158R);
    sigGd158 = sigE + sigI + sigR;

    //****************************** Gd160 **************************************
    sigE = interpol(ener, enGd160E, xGd160E, maxGd160E);
    sigI = interpol(ener, enGd160I, xGd160I, maxGd160I);
    sigR = interpol(ener, enGd160R, xGd160R, maxGd160R);
    sigGd160 = sigE + sigI + sigR;

    //******************************* NATURAL Gd *********************************
    sigGdnat = (nGd154 * sigGd154 + nGd155 * sigGd155 + nGd156 * sigGd156 + nGd157 * sigGd157 + nGd158 * sigGd158 + nGd160 * sigGd160);

    //******************************* NATURAL B **********************************
    sigBnat = (nB10 * sigB10 + nB11 * sigB11);

    //******************************* PLASTICS (CH2)n ****************************
    sigCH2 = (nC12Pl * sigC12 + nH1Pl * sigH1);

    //***************************************************************************

    //******************************** CHOOSE MODERATOR SIGMA *******************
    if (whichmod == 15){		// if moderator is D2O
      sigmod = sigD2O * barn2cm;
    } else if (whichmod == 16){		// if moderator is Graphite
      sigmod = nC12gr * sigC12 * barn2cm;
    } else if (whichmod == 19){		// if moderator is Natural Boron
      sigmod = sigBnat * barn2cm;
    } else if (whichmod == 21){		// if moderator is Plastic (C2H4)n
      sigmod = sigCH2 * barn2cm;
    }

    //******************************** CHOOSE MATERIAL SIGMA **********************
    if (whichmat == 17){		// if material is Natural Gadolinium
      sigmat = sigGdnat * barn2cm;
    } else if (whichmat == 18){		// if material is natural uranium
      sigmat = sigUnat * barn2cm;
    }

    //********************************* CHOOSE REFLECTOR SIGMA ************************
    sigreflector = sigBe9 * barn2cm;

    // 	Angles for the CM frame
    //
    modmass = MatList[whichmod].amu;
    coslamCM = 2.0 * (static_cast<double>(rand()) / static_cast<double>(RAND_MAX + 1.0)) - 1.0;
    efactorLAB = (1.0 + modmass * modmass + 2.0 * modmass * coslamCM) / pow((1.0 + modmass), 2.0);


    // Now compute macroscopic CS depending on the type of nuclides
    rc = neu->r;

    if ((rc < Rscatter) || (rc > Rscatter + Delrscatter && rc < Rreflector)){
      signow = sigmod;
      mfpnow = 1.0/signow;	// Mean free path
    } else if ((rc >= Rreflector) && (rc < Rreflector+Treflector)){
      signow = sigreflector;
      mfpnow = 1.0 / signow;
    } else if ((rc >= Rscatter) && (rc <= Rscatter + Delrscatter)){
      signow = sigmat;
      mfpnow = 1.0/signow;	// Mean free path
    }

    // Compute random displacement d that a particle may travel
    do {
      tempo = static_cast<double>(rand()) / static_cast<double>(RAND_MAX + 1.0);
      d = -log(tempo) * mfpnow;
    } while (d > mfpnow);

    // Allow the neutron to move

    xc = neu->x0 + d * neu->U;
    yc = neu->y0 + d * neu->V;
    zc = neu->z0 + d * neu->W;
    rc = pow(((xc*xc) + (yc*yc) + (zc*zc)), 0.5);

    tmpU = neu->U; tmpV = neu->V;  tmpW = neu->W;
    // If particle returns back to source again then re-initialize to current energy
    rpt = (rc < Rsource) ? true : false;
    if (rpt) initneutron(neu, rpt, false); // false for makenue as it is re-init

    displace = abs(rc - neu->r);
    ibin = (int)(neu->En /(double)binstep + 0.5);
    //cout << neu->En << "   " << ibin << endl;

    if ((rc < Rscatter) || (rc > Rscatter + Delrscatter)) {
      if (rc >= Rreflector){
	x0 = neu->x; y0 = neu->y; z0 = neu->z; r0 = neu->r;
	tempo = static_cast<double>(rand()) / static_cast<double>(RAND_MAX + 1.0);
	d = -log(tempo) * mfpnow;
	if (d <= Treflector) {
	  tempo = static_cast<double>(rand())/static_cast<double>(RAND_MAX + 1.0);
	  if (signow <= tempo) {
	    neu->lost = true;	// absorb the praticle
	    return;
	  } else {	// refletion by boundary wall
	    dtowalk = pow((((x0 - xc) * (x0 - xc)) + ((y0 - yc) * (y0 - yc)) + ((z0 - zc) * (z0 - zc))), 0.5);	// now particle is at reflector
	    reflection(x0, y0, z0, xc, yc, zc, dtowalk);
	    delE = exp(-1.0 * nBe9 * signow * barn2cm) * efactorLAB;
	    GetDC(false, tmpU, tmpV, tmpW, MatList[1].amu);	// come to new point and collide
	    neu->E0 = neu->En;   	neu->En = neu->En * (delE);
	    neu->x0 = neu->x;  	neu->y0 = neu->y; 	neu->z0 = neu->z;
	    neu->r0 = neu->r;
	    neu->x = xc;  		neu->y = yc;  		neu->z = zc;
	    neu->r = rc;
	    return;	// return 2 main and scatter again
	  }
	} else if (d > Treflector){
	  neu->lost = true;  // neutron lost in open space outside the system
	  return;
	}
      }
      // Select the nucleus to collide with
      //
      whatATOM = whichnucleus(whichmod, sigD2O, sigC12, sigBnat, sigCH2);

      if (whatATOM == 0){ 					// D
	delE = exp(-1.0 * nD_D2O * sigH2 * barn2cm) * efactorLAB; 	// energy lost
	GetDC(false, tmpU, tmpV, tmpW, MatList[0].amu);			// Dir cosine
	neu->coll++; neu->colD++; totcol++;
      } else if (whatATOM == 4) {				// C // Plastic
	delE = exp(-1.0 * nC12Pl * sigC12 * barn2cm) * efactorLAB;
	GetDC(false, tmpU, tmpV, tmpW, MatList[4].amu);
	neu->coll++; totcol++;
      } else if (whatATOM == 11){				// O
	delE = exp(-1.0 * nO_D2O * sigO16 * barn2cm) * efactorLAB; 	// energy lost
	GetDC(false, tmpU, tmpV, tmpW, MatList[11].amu);  		// Dir cosine
	neu->coll++;  neu->colO++;  totcol++;
      } else if (whatATOM == 16){				// C // Graphite
	delE = exp(-1.0 * nC12gr * sigC12 * barn2cm) * efactorLAB;
	GetDC(false, tmpU, tmpV, tmpW, MatList[16].amu);
	neu->coll++;  totcol++;
      } else if ( whatATOM == 19){				// B
	delE = exp(-1.0 * (nB10nat * sigB10 + nB11nat * sigB11) * barn2cm) * efactorLAB;
	GetDC(false, tmpU, tmpV, tmpW, MatList[19].amu);
	neu->coll++; totcol++;
      } else if (whatATOM == 20){				// H
	delE = exp(-1.0 * nH1Pl * sigH1 * barn2cm) * efactorLAB;
	GetDC(false, tmpU, tmpV, tmpW, MatList[20].amu);
	neu->coll++;  totcol++;
      }

      // 	Update energy and update position
      //
      neu->E0 = neu->En;   	neu->En = neu->En * (delE);
      neu->x0 = neu->x;  	neu->y0 = neu->y; 	neu->z0 = neu->z;
      neu->r0 = neu->r;
      neu->x = xc;  		neu->y = yc;  		neu->z = zc;
      neu->r = rc;
      nbin[ibin][0]++;           nbin[ibin][1] = neu->En;
      totbin = max2(totbin, ibin);
      //cout << "NBIN (" << ibin << ",0)= " << nbin[ibin][0] << "  totbin = " << totbin << endl;
      if (neu->En <= Ethermal) {
	neu->nthermal = true;
	nthermal++;
	return;
      }
    }else if ((rc >= Rscatter) && (rc <= Rscatter + Delrscatter)) {
      if (itype == 1){		// fission
	tempo = static_cast<double>(rand())/static_cast<double>(RAND_MAX + 1.0);
	if ((sigmat) > tempo){
	  neu->fissloss = true;  nfiss++; nlost++; totcol++;
	  neu->lost = true;
	  return;
	}
      }
      if (itype == 2){	// absorption by Gd152
	    csGd 	= sigGdnat * barn2cm;
	    csD		= nD_D2O * sigH2 * barn2cm;
	    csO		= nO_D2O * sigO16 * barn2cm;
	    CSMIN 	= min3(csGd, csD, csO);
	    CSMAX	= max3(csGd, csD, csO);
	    RANGE	= CSMAX - CSMIN;
	    TOTALCS	= csGd + csD + csO;

	    do {
	      tempo = static_cast<double>(rand())/static_cast<double>(RAND_MAX + 1.0);
	    }while ((tempo > (csD + csO)/TOTALCS) && (tempo > csGd/TOTALCS));

	    if (tempo < (csD + csO)/TOTALCS) {
	      if (csD/TOTALCS>tempo && csO/TOTALCS<tempo) whatATOM = 1;	// H2
	      if (csD/TOTALCS<tempo && csO/TOTALCS>tempo) whatATOM = 2;	// O16
	    } else {
	      whatATOM = 0;						// Gd152
	    }
	    if (whatATOM == 0){
	      neu->lost = true; neu->Gdabs = true;
	      if (whichGd == 1) nabsGd++;
	      neu->coll++; totcol++;
	      nlost++;
	      return;
	    }
	    if (whatATOM == 1) {
	      delE = exp(-1.0 * nD_D2O * sigH2 * barn2cm) * efactorLAB; // energy lost
	      GetDC(false, tmpU, tmpV, tmpW, MD);			  // Dir cosine
	      neu->coll++; neu->colD++; totcol++;
	    }
	    if (whatATOM == 2) {
	      delE = exp(-1.0 * nO_D2O * sigO16 * barn2cm) * efactorLAB; // energy lost
	      GetDC(false, tmpU, tmpV, tmpW, MO);  			  // Dir cosine
	      neu->coll++;  neu->colO++;  totcol++;
	    }
	    // Update energy and position
	    //
	    neu->E0 = neu->En;  neu->En = neu->En * (delE);
	    neu->x0 = neu->x;  	neu->y0 = neu->y;   	neu->z0 = neu->z;
	    neu->r0 = neu->r;
	    neu->x = xc; neu->y = yc; neu->z = zc; neu->r = rc;
	    if (neu->En <= Ethermal) {
	      neu->nthermal = true;
	      nthermal++;
	      return;
	    }
      }
    }
    return;
}


//
// Get the Direction Cosines of the particles
void GetDC(bool makeneu, double &U, double &V, double &W, double Atwt)
{
  double rnd1, rnd2, rnd3, s = 10.0;
  double psilab, filab;
  double cosfilab, sinfilab, cospsilab, sinpsilab;

  while (s > 1.0) {
    srand(time(NULL));
    rnd1 = static_cast<double>(rand()) / static_cast<double>(RAND_MAX + 1.0);
    rnd2 = static_cast<double>(rand()) / static_cast<double>(RAND_MAX + 1.0);
    s = (rnd1 * rnd1) + (rnd2 * rnd2);
  }

  cosfilab = (rnd1*rnd1 - rnd2*rnd2) / s;
  sinfilab = (2.0 * rnd1 * rnd2) / s;

  rnd3 = 2.0 * static_cast<double>(rand()) / static_cast<double>(RAND_MAX + 1.0) - 1.0;
  if (makeneu == true) {
  cospsilab = rnd3;
  } else {
    cospsilab = (1.0 + Atwt * rnd3) / (1.0 + (Atwt * Atwt) + 2.0 * Atwt * rnd3);
  }
  sinpsilab = 1.0 - (cospsilab * cospsilab);

  if ((1.0 - W * W) < 1.0e-4){
    U = sinpsilab * cosfilab;
    V = sinpsilab * sinfilab;
    W = cospsilab;
  } else {
    U = U * cospsilab + (U * W * sinpsilab * cosfilab - V * sinpsilab * sinfilab) / (pow((1.0 - (W * W)), 0.5));
    V = V * cosfilab + (U * sinpsilab * sinfilab + V * W * sinpsilab * cosfilab) / (pow((1.0 - (W * W)), 0.5));
    W = W * cospsilab - (pow((1.0 - (W * W)), 0.5)) * sinpsilab * cosfilab;
  }
  return;
}

//
// Initialize neutrons
void initneutron(struct particles *neu, bool rpt, bool makeneu)
{
  double costheta, sintheta, cosphi, sinphi, phi;
  double tmpU, tmpV, tmpW;
  double ri, ro, rfactor, reff;


  // generate random polar angle theta and azimuthal angle phi
  srand(time(NULL));
  costheta = 2.0 * static_cast<double>(rand()) / static_cast<double>(RAND_MAX + 1.0) - 1.0;
  sintheta = pow((1.0 - costheta * costheta), 0.5);
  phi = 2.0 * PI() * static_cast<double>(rand()) / static_cast<double>(RAND_MAX + 1.0);
  cosphi = cos(phi);    sinphi = sin(phi);

  // generate random source point between ri and ro
  ri = 1.0e-3;    ro = Rsource;
  //rfactor = static_cast<double>(rand()) / static_cast<double>(RAND_MAX + 1.0);
  //reff = pow(((ri * ri * ri) + rfactor *((ro * ro * ro) - (ri * ri * ri))), (1.0 / 3.0));
  reff = ro;
  // generate x0, y0, z0 and r0
  neu->x0 = reff * sintheta * cosphi;
  neu->y0 = reff * sintheta * sinphi;
  neu->z0 = reff * costheta;
  neu->r0 = pow(((neu->x0 * neu->x0) + (neu->y0 * neu->y0) + (neu->z0 * neu->z0)), 0.50);



  // if not repeat particle, energy is initial energy
  if (!rpt) neu->E0 = En0;

  // In initialization, xo =x
  neu->x = neu->x0;
  neu->y = neu->y0;
  neu->z = neu->z0;
  neu->r = neu->r0;

  // Random number for initial direction cosines
  tmpU = static_cast<double>(rand()) / static_cast<double>(RAND_MAX + 1.0);
  tmpV = static_cast<double>(rand()) / static_cast<double>(RAND_MAX + 1.0);
  tmpW = static_cast<double>(rand()) / static_cast<double>(RAND_MAX + 1.0);

  // Enter the random numbers to get direction cosines function to get
  // present direction cosines and allot them to neutron
  GetDC(makeneu, tmpU, tmpV, tmpW, 0.0);	// the last entry is 0.0 as n is created and no collision is there
  neu->U = tmpU;  neu->V = tmpV; neu->W = tmpW;
  makeneu = false;

  if (!rpt) neu->En = neu->E0;

  // If it is first initialization initialize counters to zero
  if (!rpt){
    neu->coll = 0; neu->colD = 0; neu->colO = 0; neu->nthermal = 0;
    neu->lost = false; neu->Gdabs = false; neu->fissloss = false;
  }
  return;
}


//
// Initialize scatterers
void initscat(void)
{
  int i;

  //H2
  MatList[0].name	= "H2";
  MatList[0].isElement	= 'Y';
  MatList[0].compnum	= 0;
  for (i = 0; i < maxcomp; i++) MatList[0].comp[i] = 0;
  MatList[0].molefrac	= 100.0;		// only isotope for D2O
  MatList[0].amu	= MD;			// amu
  MatList[0].density	= 0.088;		// g/cc density for solid!!!!!!!

  //Be9
  MatList[1].name	= "Be9";
  MatList[1].isElement	= 'Y';
  MatList[1].compnum	= 0;
  for (i = 0; i < maxcomp; i++) MatList[1].comp[i] = 0;
  MatList[1].molefrac	= 100.0;		// only isotope for natural Be
  MatList[1].amu	= 9.012182;  		// amu
  MatList[1].density	= 1.848;		// g/cc
  nBe9			= (MatList[1].density * nAvo)/MatList[1].amu;

  //B10
  MatList[2].name	= "B10";
  MatList[2].isElement	= 'Y';
  MatList[2].compnum	= 0;
  for (i = 0; i < maxcomp; i++) MatList[2].comp[i] = 0;
  MatList[2].molefrac	= 19.9;			// mole frac for natural Boron
  MatList[2].amu	= 10.811;		// amu
  MatList[2].density	= 2.460;		// g/cc
  nB10			= (MatList[2].density * nAvo)/MatList[2].amu;

  //B11
  MatList[3].name	= "B11";
  MatList[3].isElement	= 'Y';
  MatList[3].compnum	= 0;
  for (i = 0; i < maxcomp; i++) MatList[3].comp[i] = 0;
  MatList[3].molefrac	= 80.1;			// mole frac for natural Boron
  MatList[3].amu	= 10.811;		// amu
  MatList[3].density	= 2.460;		// g/cc
  nB11			= (MatList[3].density * nAvo) / MatList[3].amu;

  //C12
  MatList[4].name	= "C12";
  MatList[4].isElement	= 'Y';
  MatList[4].compnum	= 0;
  for (i = 0; i < maxcomp; i++) MatList[4].comp[i] = 0;
  MatList[4].molefrac	= 100.0;		// mole frac for GRAPHITE
  MatList[4].amu	= 12.0107;		// amu
  MatList[4].density	= 2.267;		// g/cc

  //Gd154
  MatList[5].name	= "Gd154";
  MatList[5].isElement	= 'Y';
  MatList[5].compnum	= 0;
  for (i = 0; i < maxcomp; i++) MatList[5].comp[i] = 0;
  MatList[5].molefrac	= 2.18;
  MatList[5].amu	= MGd;			// amu
  MatList[5].density	= 7.901;		// g/cc

  //Gd155
  MatList[6].name	= "Gd155";
  MatList[6].isElement	= 'Y';
  MatList[6].compnum	= 0;
  for (i = 0; i < maxcomp; i++) MatList[6].comp[i] = 0;
  MatList[6].molefrac	= 14.8;
  MatList[6].amu	= MGd;			// amu
  MatList[6].density	= 7.901;		// g/cc

  //Gd156
  MatList[7].name	= "Gd156";
  MatList[7].isElement	= 'Y';
  MatList[7].compnum	= 0;
  for (i = 0; i < maxcomp; i++) MatList[7].comp[i] = 0;
  MatList[7].molefrac	= 20.47;
  MatList[7].amu	= MGd;			// amu
  MatList[7].density	= 7.901;		// g/cc

  //Gd157
  MatList[8].name	= "Gd157";
  MatList[8].isElement	= 'Y';
  MatList[8].compnum	= 0;
  for (i = 0; i < maxcomp; i++) MatList[8].comp[i] = 0;
  MatList[8].molefrac	= 15.65;
  MatList[8].amu	= MGd;			// amu
  MatList[8].density	= 7.901;		// g/cc

  //Gd158
  MatList[9].name	= "Gd158";
  MatList[9].isElement	= 'Y';
  MatList[9].compnum	= 0;
  for (i = 0; i < maxcomp; i++) MatList[9].comp[i] = 0;
  MatList[9].molefrac	= 24.84;
  MatList[9].amu	= MGd;			// amu
  MatList[9].density	= 7.901;		// g/cc

  //Gd160
  MatList[10].name	= "Gd160";
  MatList[10].isElement	= 'Y';
  MatList[10].compnum	= 0;
  for (i = 0; i < maxcomp; i++) MatList[10].comp[i] = 0;
  MatList[10].molefrac	= 21.86;
  MatList[10].amu	= MGd;			// amu
  MatList[10].density	= 7.901;		// g/cc

  //O16
  MatList[11].name	= "O16";
  MatList[11].isElement	= 'Y';
  MatList[11].compnum	= 0;
  for (i = 0; i < maxcomp; i++) MatList[11].comp[i] = 0;
  MatList[11].molefrac	= 100.0;
  MatList[11].amu	= MO;			// amu
  MatList[11].density	= 1.495;		// g/cc *** density of solid !!!!

  //U234
  MatList[12].name	= "U234";
  MatList[12].isElement	= 'Y';
  MatList[12].compnum	= 0;
  for (i = 0; i < maxcomp; i++) MatList[12].comp[i] = 0;
  MatList[12].molefrac	= 0.0056;
  MatList[12].amu	= 234.0409521;		// amu
  MatList[12].density	= 19.05;		// g/cc

  //U235
  MatList[13].name	= "U235";
  MatList[13].isElement	= 'Y';
  MatList[13].compnum	= 0;
  for (i = 0; i < maxcomp; i++) MatList[13].comp[i] = 0;
  MatList[13].molefrac	= 0.7205;
  MatList[13].amu	= 235.0439299;		// amu
  MatList[13].density	= 19.05;		// g/cc

  //U238
  MatList[14].name	= "U238";
  MatList[14].isElement	= 'Y';
  MatList[14].compnum	= 0;
  for (i = 0; i < maxcomp; i++) MatList[14].comp[i] = 0;
  MatList[14].molefrac	= 99.2739;
  MatList[14].amu	= 238.0507882;		// amu
  MatList[14].density	= 19.05;		// g/cc


  //D2O
  MatList[15].name	= "D2O";
  MatList[15].isElement	= 'N';		// Allotrope
  MatList[15].compnum	= 2;
  for (i = 0; i < maxcomp; i++) MatList[15].comp[i] = 0;
  MatList[15].comp[0] 	= 0;
  MatList[15].comp[1]	= 11;
  MatList[15].amu	= (2.0 * MatList[0].amu + 1.0 * MatList[11].amu);
  MatList[15].density 	= 1.111;		// g/cc
  nD_D2O 		= (MatList[15].density * nAvo * 2.0) / MD;
  nO_D2O 		= (MatList[15].density * nAvo * 1.0) / MO;


  //Graphite
  MatList[16].name	= "GRAPHITE";
  MatList[16].isElement	= 'N';
  MatList[16].compnum	= 0;
  for (i = 0; i < maxcomp; i++) MatList[16].comp[i] = 0;
  MatList[16].amu	= MatList[4].amu;	// C12
  MatList[16].density	= MatList[4].density;
  nC12gr		= (MatList[16].density *nAvo) / MatList[16].amu;

  //Natural Gd
  MatList[17].name	= "Natural Gd";
  MatList[17].isElement	= 'Y';
  MatList[17].compnum	= 0;
  for (i = 0; i < maxcomp; i++) MatList[17].comp[i] = 0;
  MatList[17].amu	= (MatList[5].molefrac/100.0*MatList[5].amu)+(MatList[6].molefrac/100.0*MatList[6].amu)
  +(MatList[7].molefrac/100.0*MatList[7].amu)+(MatList[8].molefrac/100.0*MatList[8].amu)+
  (MatList[9].molefrac/100.0*MatList[9].amu)+(MatList[10].molefrac/100.0*MatList[10].amu);
  MatList[17].density	= 7.901;		// g/cc
  nGd154		= ((Gdcon/100.0)*TotalGd_slab_wt)*(MatList[5].molefrac/100.0)*(nAvo/MatList[17].amu);
  nGd155		= ((Gdcon/100.0)*TotalGd_slab_wt)*(MatList[6].molefrac/100.0)*(nAvo/MatList[17].amu);
  nGd156		= ((Gdcon/100.0)*TotalGd_slab_wt)*(MatList[7].molefrac/100.0)*(nAvo/MatList[17].amu);
  nGd157		= ((Gdcon/100.0)*TotalGd_slab_wt)*(MatList[8].molefrac/100.0)*(nAvo/MatList[17].amu);
  nGd158		= ((Gdcon/100.0)*TotalGd_slab_wt)*(MatList[9].molefrac/100.0)*(nAvo/MatList[17].amu);
  nGd160		= ((Gdcon/100.0)*TotalGd_slab_wt)*(MatList[10].molefrac/100.0)*(nAvo/MatList[17].amu);

  //Natural Uranium
  MatList[18].name	= "Natural Uranium";
  MatList[18].isElement	= 'Y';   // Mixture of isotopes
  MatList[18].compnum	= 0;
  for (i = 0; i < maxcomp; i++) MatList[18].comp[i] = 0;
  MatList[18].amu	= (!enrich) ? (MatList[13].molefrac/100.0*MatList[13].amu)+(MatList[14].molefrac/100.0*MatList[14].amu) : (enrich/100.0*MatList[13].amu)+((1.0-(enrich/100.0))*MatList[14].amu);
  MatList[18].density	= 19.05;		// g/cc
  nU235			= (!enrich) ? (MatList[13].density * nAvo)/MatList[18].amu : (enrich/100.0)*MatList[18].density*nAvo/MatList[13].amu;
  nU238			= (!enrich) ? (MatList[14].density * nAvo)/MatList[18].amu : (1.0 - (enrich/100.0))*MatList[18].density*nAvo/MatList[18].amu;


  // Natural Boron
  MatList[19].name	= "Natural Boron";
  MatList[19].isElement = 'Y';
  MatList[19].compnum	= 0;
  for (i = 0; i < maxcomp; i++) MatList[19].comp[i] = 0;
  MatList[19].density	= 2.45;	// g/cc
  MatList[19].amu	= MatList[2].molefrac / 100.0 * MatList[2].amu + MatList[3].molefrac / 100.0 * MatList[3].amu;
  nB10nat		= MatList[19].density * (nAvo/MatList[19].amu) * MatList[2].molefrac/100.0;
  nB11nat		= MatList[19].density * (nAvo/MatList[19].amu) * MatList[3].molefrac/100.0;

  // Hydrogen
  MatList[20].name	= "Hydrogen";
  MatList[20].isElement = 'Y';
  MatList[20].compnum	= 0;
  for (i = 0; i < maxcomp; i++) MatList[20].comp[i] = 0;
  MatList[20].molefrac	= 100.0;
  MatList[20].density	= 0.0899e-3;			// g/cc
  MatList[20].amu	= 1.001;			// amu
  nH1			= (MatList[20].density*nAvo) / MatList[20].amu;


  // Plastics (C2H4)n
  MatList[21].name	= "Plastic C2H4";
  MatList[21].isElement	= 'N';
  MatList[21].compnum	= 2;
  for (i = 0; i < maxcomp; i++) MatList[21].comp[i] = 0;
  MatList[21].comp[0]	= 4;
  MatList[21].comp[1]	= 16;
  MatList[21].amu	= 2.0 * MatList[4].amu + 4.0 * MatList[20].amu;
  MatList[21].density	= 1260.0;			// g/cc
  nH1Pl			= (MatList[21].density * nAvo * 2.0) / MatList[20].amu;
  nC12Pl		= (MatList[21].density * nAvo * 4.0) / MatList[4].amu;


  for (i = 0; i < maxmat; i++){
    if (MatList[i].isElement == 'N') MatList[i].molefrac =  100.0;
  }
  cout << " All materials have been constructed" << endl;

  return;
}


void pause(void)
{
    std::cin.clear(); //clear errors
    std::cin.sync(); //clear the buffer
    std::cin.get();
    return;
}


//
// Input required data
void inputdata(void)
{
  int nchr = 100, expo1, expo2;
  char line[nchr], *thestrptr;
  string str1;
  double  xtra;

  iFile.open(infile00.c_str(), ifstream::in);
  if (iFile){

    iFile.getline(line, nchr);
    itype		= (int)strtod(line, &thestrptr);

    iFile.getline(line, nchr);
    Rsource             = strtod(line, &thestrptr);

    iFile.getline(line, nchr);
    Rscatter            = strtod(line, &thestrptr);

    iFile.getline(line, nchr);
    Delrscatter         = strtod(line, &thestrptr);

    iFile.getline(line, nchr);
    Rreflector          = strtod(line, &thestrptr);

    iFile.getline(line, nchr);
    Treflector          = strtod(line, &thestrptr);

    iFile.getline(line, nchr);
    TotalGd_slab_wt     = strtod(line, &thestrptr);

    iFile.getline(line, nchr);
    enrich              = strtod(line, &thestrptr);

    iFile.getline(line, nchr);
    Unatwt              = strtod(line, &thestrptr);

    iFile.getline(line, nchr);
    Gdcon               = strtod(line, &thestrptr);

    iFile.getline(line, nchr);
    whichmoderator      = (int)strtod(line, &thestrptr); // 15 = D2O 16 = Graphite, 21 Plastics

    iFile.getline(line, nchr);
    whichmat            = (int)strtod(line, &thestrptr); // 18 = U 17 = Gadolinium

    iFile.getline(line, nchr);
    nparticles          = (long)strtod(line, &thestrptr);

    iFile.getline(line, nchr);
    En0                 = strtod(line, &thestrptr);
  }
  iFile.close();

  expo1 = (int)log10(1.0*nparticles);
  expo2 = (int)log10(En0);
  xtra	= (expo2 > 4) ? (expo2/1.0e+6) : (expo2/1.0e+3);
  if (expo2 <= 4) str1 = " keV ";
  if (expo2 >= 5) str1 = " MeV ";

  cout << endl << " PARAMETERS READ BY THE SYSTEM.........\n";
  if (itype == 1) cout << " Type of reaction = FISSION ";
  if (itype == 1) cout << " Material used natural U" << endl;
  if (itype == 2) cout << " Type of reaction = Absorption ";
  cout << " by the material used :  ";
  if (itype == 2 && whichmat == 17) cout << " Natural Gd\n ";
  if (itype == 2) cout << " Moderator used in the system ---> ";
  if (itype == 2 && whichmoderator == 15) cout << "  D2O " << endl;
  if (itype == 2 && whichmoderator == 16) cout << "  Graphite " << endl;
  if (itype == 2 && whichmoderator == 21) cout << "  Plastics (C2H4)n " << endl;

  cout << " Source radius = " << Rsource << " cms \n";
  cout << " Radius of Scatterer (ID) = " << Rscatter << " cms \n";
  cout << " Thickness of the Scatterer = " << Delrscatter << " cms\n";
  cout << " Radius of the Reflector (ID) = " << Rreflector << " cms\n";
  cout << " Thickness of the Reflector = " << Treflector << " cms\n";
  cout << " Total number of neutrons = " << "10^" << expo1 << "  with energy " << En0 << " eV " << "  \n\n\n";
  return;
}


// 	Output some data on the screen
//
void outputdata(long int tcol, long int ncol2, long int coxy, long int cdeu, long int cfiss, long int nGd, long int nth, double avcol)
{
  double stddev = 0.0, dev = 0.0;
  string str1;

  //ostringstream str11, str12;
  //char tempf[50];

  //time_t secs=time(0);
  //tm *t;=localtime(&secs);
  //This prints the date in ISO format.
  //cout << t->tm_mon+1 << "  " << t->tm_mday << "  " << t->tm_year+1900 << endl;
  //printf("%04d-%02d-%02d\n",t->tm_year+1900,t->tm_mon+1,t->tm_mday);
  //This prints the date and time in ISO format.
  //printf("%04d-%02d-%02d %02d:%02d:%02d\n",t->tm_year+1900,t->tm_mon+1,t->tm_mday,t->tm_hour,t->tm_min,t->tm_sec);
  //pause();

  expo1 = (int)log10(1.0*nparticles);
  expo2 = (int)log10(En0);
  if (expo2 == 3) str1 = " keV ";
  if (expo2 == 6) str1 = " MeV ";
  str11 << expo1;  str12 << En0;
  outfile01 = outfile01+"_10^"+str11.str()+"_"+str12.str()
  +"REACTION"+strtype.str()+"_"+strd.str()+strm.str()+stry.str()+".dat";
  strcpy(tempf, outfile01.c_str());
  oFile1.open(tempf);

  avcol = tcol/(double)nparticles;
  dev = (ncol2/(double)nparticles) - (avcol * avcol);
  stddev = pow(dev, 0.5) / (double)nparticles;
  cout << endl << endl << "********************  RESULTS ********************" << endl << endl;
  cout << "Neutron Numbers = " << "10^" << expo1 << "  with energy = " << (En0) << "  eV " << endl;
  cout << "Moderator used --->  " << MatList[whichmoderator].name<< endl;
  cout << " Total numbers of collisions = " << tcol << "  MEAN = " << avcol << endl;
  //cout << " D/total collisions = " << cdeu/(double)tcol*100.0 << "  %" << endl;
  //cout << " O/total collisions = " << coxy/(double)tcol*100.0 << "  %" << endl;
  if (cfiss) cout << " fission/N = " << cfiss/(double)nparticles*100.0 << "  %"  << endl;
  if (nGd) cout << " Gd Absorption/N = " << nGd/(double)nparticles*100.0 << "  %" << endl;
  if (nth) cout << " Thermalized/N = " << nth/(double)nparticles*100.0 << "  %" << endl;
  cout << "  The Standard deviation = " << stddev << endl << endl;

  oFile1 << endl << endl << "********************  RESULTS ********************" << endl << endl;
  oFile1 << endl << " PARAMETERS READ BY THE SYSTEM.........\n";
  if (itype == 1) oFile1 << " Type of reaction = FISSION ";
  if (itype == 1) oFile1 << " Material used natural U" << endl;
  if (itype == 2) oFile1 << " Type of reaction = Absorption " << endl;
  if (itype == 2) oFile1 << " Moderator used in the system ---> ";
  if (itype == 2) if (whichmoderator == 15) oFile1 << "  D2O " << endl;
  if (itype == 2) if (whichmoderator == 16) oFile1 << "  Graphite " << endl;
  if (itype == 2) if (whichmoderator == 21) oFile1 << "  Plastics (C2H4)n " << endl;
  oFile1 << " Source radius = " << Rsource << " cms \n";
  oFile1 << " Radius of Scatterer (ID) = " << Rscatter << " cms \n";
  oFile1 << " Thickness of the Scatterer = " << Delrscatter << " cms\n";
  oFile1 << " Radius of the Reflector (ID) = " << Rreflector << " cms\n";
  oFile1 << " Thickness of the Reflector = " << Treflector << " cms\n";
  oFile1 << " Total number of neutrons = " << "10^" << expo1 << "  with energy " << En0 << " eV " << "  \n\n\n";
  oFile1 << " ----------------------  DATA -----------------------\n";
  oFile1 << "Neutron Numbers = " << "10^" << expo1 << "  with energy = " << (En0) << "  eV " << endl;
  oFile1 << "Moderator used --->  " << MatList[whichmoderator].name << endl;
  oFile1 << " Total numbers of collisions = " << tcol << "  MEAN = " << avcol << endl;
  //cout << " D/total collisions = " << cdeu/(double)tcol*100.0 << "  %" << endl;
  //cout << " O/total collisions = " << coxy/(double)tcol*100.0 << "  %" << endl;
  if (cfiss) oFile1 << " fission/N = " << cfiss/(double)nparticles*100.0 << "  %"  << endl;
  if (nGd) oFile1 << " Gd Absorption/N = " << nGd/(double)nparticles*100.0 << "  %" << endl;
  if (nth) oFile1 << " Thermalized/N = " << nth/(double)nparticles*100.0 << "  %" << endl;
  oFile1 << "  The Standard deviation = " << stddev << endl << endl;
  oFile1.close();
  return;
}


// Read data file for scatterers cross section data
// data downloaded from IAEA-NDS data explorer center
void readscatdata(void)
{
  int i;
  float k1, k2;
  iFile.open(infile01.c_str());				// H2 E
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enH2E[i] = k1; xH2E[i++] = k2;
  }
  maxH2E = i-1;
  iFile.close();
 // cout << " D-E data loaded... " << endl;


  iFile.open(infile02.c_str());				// H2 I
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enH2I[i] = k1; xH2I[i++] = k2;
  }
  maxH2I = i-1;
  iFile.close();
  // cout << " D-I data loaded...." << endl;

  iFile.open(infile03.c_str());				// Be9 N2N
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enBe9N2N[i] = k1; xBe9N2N[i++] = k2;
  }
  maxBe9N2N = i-1;
  iFile.close();
  // cout << " Be9-N2N data loaded...." << endl;

  iFile.open(infile04.c_str());				// Be9 E
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enBe9E[i] = k1; xBe9E[i++] = k2;
  }
  maxBe9E = i - 1;
  iFile.close();
  //cout << " Be9-E data loaded...." << endl;

  iFile.open(infile05.c_str());				// Be9 I
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enBe9I[i] = k1; xBe9I[i++] = k2;
  }
  maxBe9I = i - 1;
  iFile.close();
  //cout << " Be9-I data loaded...." << endl;

  iFile.open(infile06.c_str());				// Be9 R
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enBe9R[i] = k1; xBe9R[i++] = k2;
  }
  maxBe9R = i - 1;
  iFile.close();
  //cout << " Be9-R data loaded...." << endl;

  iFile.open(infile07.c_str());				// B10 E
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enB10E[i] = k1; xB10E[i++] = k2;
  }
  maxB10E = i - 1;
  iFile.close();
  //cout << " B10-E data loaded...." << endl;

  iFile.open(infile08.c_str());				// B10 I
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enB10I[i] = k1; xB10I[i++] = k2;
  }
  maxB10I = i - 1;
  iFile.close();
  //cout << " B10-I data loaded...." << endl;

  iFile.open(infile09.c_str());				// B10 R
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enB10R[i] = k1; xB10R[i++] = k2;
  }
  maxB10R = i - 1;
  iFile.close();
  //cout << " B10-R data loaded...." << endl;

  iFile.open(infile10.c_str());				// B11 E
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enB11E[i] = k1; xB11E[i++] = k2;
  }
  maxB11E = i - 1;
  iFile.close();
  //cout << " B11-E data loaded...." << endl;

  iFile.open(infile11.c_str());				// B11 I
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enB11I[i] = k1; xB11I[i++] = k2;
  }
  maxB11I = i - 1;
  iFile.close();
  //cout << " B11-I data loaded...." << endl;

  iFile.open(infile12.c_str());				// B11 R
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enB11R[i] = k1; xB11R[i++] = k2;
  }
  maxB11R = i - 1;
  iFile.close();
  //cout << " B11-R data loaded...." << endl;

  iFile.open(infile13.c_str());				// Gd154 E
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enGd154E[i] = k1; xGd154E[i++] = k2;
  }
  maxGd154E = i - 1;
  iFile.close();
  //cout << " Gd154-E data loaded...." << endl;

  iFile.open(infile14.c_str());				// Gd154 I
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enGd154I[i] = k1; xGd154I[i++] = k2;
  }
  maxGd154I = i - 1;
  iFile.close();
  //cout << " Gd154-I data loaded...." << endl;

  iFile.open(infile15.c_str());				// Gd154 R
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enGd154R[i] = k1; xGd154R[i++] = k2;
  }
  maxGd154R = i - 1;
  iFile.close();
  //cout << " Gd154-R data loaded...." << endl;

  iFile.open(infile16.c_str());				// Gd155 E
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enGd155E[i] = k1; xGd155E[i++] = k2;
  }
  maxGd155E = i - 1;
  iFile.close();
  //cout << " Gd155-E data loaded...." << endl;

  iFile.open(infile17.c_str());				// Gd155 I
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enGd155I[i] = k1; xGd155I[i++] = k2;
  }
  maxGd155I = i - 1;
  iFile.close();
  //cout << " Gd155-I data loaded...." << endl;

  iFile.open(infile18.c_str());				// Gd155 R
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enGd155R[i] = k1; xGd155R[i++] = k2;
  }
  maxGd155R = i - 1;
  iFile.close();
  //cout << " Gd155-R data loaded...." << endl;

  iFile.open(infile19.c_str());				// Gd156 E
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enGd156E[i] = k1; xGd156E[i++] = k2;
  }
  maxGd156E = i - 1;
  iFile.close();
  //cout << " Gd156-E data loaded...." << endl;

  iFile.open(infile20.c_str());				// Gd156I
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enGd156I[i] = k1; xGd156I[i++] = k2;
  }
  maxGd156I = i - 1;
  iFile.close();
  //cout << " Gd156-I data loaded...." << endl;

  iFile.open(infile21.c_str());				// Gd156 R
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enGd156R[i] = k1; xGd156R[i++] = k2;
  }
  maxGd156R = i - 1;
  iFile.close();
  //cout << " Gd156-R data loaded...." << endl;

  iFile.open(infile22.c_str());				// Gd157 E
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enGd157E[i] = k1; xGd157E[i++] = k2;
  }
  maxGd157E = i - 1;
  iFile.close();
  //cout << " Gd157-E data loaded...." << endl;

  iFile.open(infile23.c_str());				// Gd157 I
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enGd157I[i] = k1; xGd157I[i++] = k2;
  }
  maxGd157I = i - 1;
  iFile.close();
  //cout << " Gd157-I data loaded...." << endl;

  iFile.open(infile24.c_str());				// Gd157 R
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enGd157R[i] = k1; xGd157R[i++] = k2;
  }
  maxGd157R = i - 1;
  iFile.close();
  //cout << " Gd157-R data loaded...." << endl;

  iFile.open(infile25.c_str());				// Gd158 E
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enGd158E[i] = k1; xGd158E[i++] = k2;
  }
  maxGd158E = i - 1;
  iFile.close();
  //cout << " Gd158-E data loaded...." << endl;

  iFile.open(infile26.c_str());				// Gd158 I
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enGd158I[i] = k1; xGd158I[i++] = k2;
  }
  maxGd158I = i - 1;
  iFile.close();
  //cout << " Gd158-I data loaded...." << endl;

  iFile.open(infile27.c_str());				// Gd158 R
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enGd158R[i] = k1; xGd158R[i++] = k2;
  }
  maxGd158R = i - 1;
  iFile.close();
  //cout << " Gd158-R data loaded...." << endl;

  iFile.open(infile28.c_str());				// Gd160 E
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enGd160E[i] = k1; xGd160E[i++] = k2;
  }
  maxGd160E = i - 1;
  iFile.close();
  //cout << " Gd160-E data loaded...." << endl;

  iFile.open(infile29.c_str());				// Gd160 I
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enGd160I[i] = k1; xGd160I[i++] = k2;
  }
  maxGd160I = i - 1;
  iFile.close();
  //cout << " Gd160-I data loaded...." << endl;

  iFile.open(infile30.c_str());				// Gd160 R
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enGd160R[i] = k1; xGd160R[i++] = k2;
  }
  maxGd160R = i - 1;
  iFile.close();
  //cout << " Gd160-R data loaded...." << endl;

  iFile.open(infile31.c_str());				// C12 E
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enC12E[i] = k1; xC12E[i++] = k2;
  }
  maxC12E = i - 1;
  iFile.close();
  //cout << " C12-E data loaded...." << endl;

  iFile.open(infile32.c_str());				// C12 I
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enC12I[i] = k1; xC12I[i++] = k2;
  }
  maxC12I = i - 1;
  iFile.close();
  //cout << " C12-I data loaded...." << endl;

  iFile.open(infile33.c_str());				// C12 R
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enC12R[i] = k1; xC12R[i++] = k2;
  }
  maxC12R = i - 1;
  iFile.close();
  //cout << " C12-R data loaded...." << endl;

  iFile.open(infile34.c_str());				// O16 E
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enO16E[i] = k1; xO16E[i++] = k2;
  }
  maxO16E = i - 1;
  iFile.close();
  //cout << " O16-E data loaded...." << endl;

  iFile.open(infile35.c_str());				// O16I
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enO16I[i] = k1; xO16I[i++] = k2;
  }
  maxO16I = i - 1;
  iFile.close();
  //cout << " O16-I data loaded...." << endl;

  iFile.open(infile36.c_str());				// O16 R
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enO16R[i] = k1; xO16R[i++] = k2;
  }
  maxO16R = i - 1;
  iFile.close();
  //cout << " O16-R data loaded...." << endl;

  iFile.open(infile37.c_str());				// U234 E
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enU234E[i] = k1; xU234E[i++] = k2;
  }
  maxU234E = i - 1;
  iFile.close();
  //cout << " U234-E data loaded...." << endl;

  iFile.open(infile38.c_str());				// U234 F
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enU234F[i] = k1; xU234F[i++] = k2;
  }
  maxU234F = i - 1;
  iFile.close();
  //cout << " U234-F data loaded...." << endl;

  iFile.open(infile39.c_str());				// U234 I
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enU234I[i] = k1; xU234I[i++] = k2;
  }
  maxU234I = i - 1;
  iFile.close();
  //cout << " U234-I data loaded...." << endl;

  iFile.open(infile40.c_str());				// U234 R
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enU234R[i] = k1; xU234R[i++] = k2;
  }
  maxU234R = i - 1;
  iFile.close();
  //cout << " U234-R data loaded...." << endl;

  iFile.open(infile41.c_str());				// U235 E
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enU235E[i] = k1; xU235E[i++] = k2;
  }
  maxU235E = i - 1;
  iFile.close();
  //cout << " U235-E data loaded...." << endl;

  iFile.open(infile42.c_str());				// U235 F
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enU235F[i] = k1; xU235F[i++] = k2;
  }
  maxU235F = i - 1;
  iFile.close();
  //cout << " U235-F data loaded...." << endl;

  iFile.open(infile43.c_str());				// U235 I
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enU235I[i] = k1; xU235I[i++] = k2;
  }
  maxU235I = i - 1;
  iFile.close();
  //cout << " U235-I data loaded...." << endl;

  iFile.open(infile44.c_str());				// U235 R
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enU235R[i] = k1; xU235R[i++] = k2;
  }
  maxU235R = i - 1;
  iFile.close();
  //cout << " U235-R data loaded...." << endl;

  iFile.open(infile45.c_str());				// U238 E
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enU238E[i] = k1; xU238E[i++] = k2;
  }
  maxU238E = i - 1;
  iFile.close();
  //cout << " U238-E data loaded...." << endl;

  iFile.open(infile46.c_str());				// U238 F
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enU238F[i] = k1; xU238F[i++] = k2;
  }
  maxU238F = i - 1;
  iFile.close();
  //cout << " U238-F data loaded...." << endl;

  iFile.open(infile47.c_str());				// U238 I
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enU238I[i] = k1; xU238I[i++] = k2;
  }
  maxU238I = i - 1;
  iFile.close();
  //cout << " U238-I data loaded...." << endl;

  iFile.open(infile48.c_str());				// U238 R
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enU238R[i] = k1; xU238R[i++] = k2;
  }
  maxU238R = i - 1;
  iFile.close();
  //cout << " U238-R data loaded...." << endl;

  iFile.open(infile49.c_str());				// H1 E
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enH1E[i] = k1; xH1E[i++] = k2;
  }
  maxH1E = i - 1;
  iFile.close();
  //cout << " H1-E data loaded...." << endl;


  iFile.open(infile50.c_str());				// H1 I
  i = 0;
  while (true){
    iFile >> k1 >> k2;
    if (iFile.eof()) break;
    enH1I[i] = k1; xH1I[i++] = k2;
  }
  maxH1I = i - 1;
  iFile.close();
  //cout << " H1-I data loaded....." << endl;


  cout << " \nData files for scattering data have been read\n \n \n";
}


//
// This is main code
int main(int argc, char **argv){
  clock_t t1, t2;
  int i;
  long int ncol2 = 0, noxy=0, ndeu=0, nth=0;
  long int nGd = 0, coxy=0, cdeu=0, cfiss=0, tcol=0, cnfis=0;
  double nk[maxbin][2];
  double enrich, Gdwt;
  double avcol;
  bool rpt = false, makeneu = true;
  int day, mon, yr;

  time_t secs=time(0);
  tm *t=localtime(&secs);
  struct particles neutron;

    t1 = clock();


    inputdata();
    readscatdata();
    initscat();
    binstep = (int)((En0) / (double)maxbin + 0.5);
    day = t->tm_mday;
    mon = t->tm_mon+1;
    yr = t->tm_year+1900;
    strd << day;   strm << mon;   stry << yr;
    strtype << itype;

  for (i = 0; i <  nparticles; i++){
    neutron.tag = i;
    initneutron(&neutron, rpt, makeneu); 			// initialize this neutron
    while (!neutron.lost && neutron.En > Ethermal){		// if not lost continue for this neutron
      scatter(&neutron, whichmoderator, whichmat, i, rpt);	// if !lost + !thermal re-scatter
      if (neutron.lost>0 || neutron.nthermal>0) break;
    }
    noxy += neutron.colO;    ndeu += neutron.colD;    //totcol += neutron.coll;
    col2 += neutron.coll * neutron.coll;
  }
    outputdata(totcol, col2, noxy, ndeu, nfiss, nabsGd, nthermal, avcol);
    outfile00 = outfile00+"_10^"+str11.str()+"_"+str12.str()+
    "REACTION"+strtype.str()+"_"+strd.str()+strm.str()+stry.str()+".dat";
    strcpy(tempf, outfile00.c_str());
    oFile2.open(tempf);
    //oFile2.open(outfile00);
    for (i = 0; i < totbin; i++){
      oFile2 << i << "    " << nbin[i][0] << "    " << nbin[i][1] << endl;
      //cout << nk[i][0] << "   " << nk[i][1] << endl;
    }
    oFile2.close();
    //outputdata(totcol, col2, noxy, ndeu, nfiss, nabsGd, nthermal, avcol);
    t2 = clock();
    int timesec = (int)(((double)t2 - (double)t1)/CLOCKS_PER_SEC+0.5);
    cout << "  Time taken = " << timesec/60 << "  mins " << timesec % 60 << " secs " << endl << endl;

  return 0;
}
