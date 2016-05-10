#include <iostream>
#include <cstdio>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cmath>
#include <ctime>
#include <climits>
#include <string>

using namespace std;

// declare some small parameters and functions
#define PI() (2.0 * asin(1.0))
#define max2(a, b) ((a > b) ? a : b)
#define min2(a, b) ((a < b) ? a : b)
#define fnamlen 201


// declare structures of particles and nuclides
struct particles{
  unsigned long tag;
  double x0, y0, z0, r0, E0;
  double x, y, z, r, En;
  double U, V, W;
  unsigned long coll, colD, colO;
  bool nthermal, lost, Gdabs, fissloss;
};

struct nuclides{
  string name;
  char isElement;
  int compnum;
  int comp[5];
  double molefrac, amu, density;
};



// declare file related parameters
//	filenames
string pathData		= "/home/vega/WORK/Fission_Monte_Carlo/NUCLEAR_DATA/";
string pathInp          = "/home/vega/WORK/Fission_Monte_Carlo/Serial/";
string infile00   	= pathInp+"Nuclear_VEGA_SER_U_01.inp";		// _01 for 64Gd    _02 for 92U
string infile01  	= pathData+"all_FINAL/1H2_N_EL.dat";
string infile02  	= pathData+"all_FINAL/1H2_N_INEL.dat";
string infile03  	= pathData+"all_FINAL/4Be9_N2N.dat";
string infile04  	= pathData+"all_FINAL/4Be9_N_EL.dat";
string infile05  	= pathData+"all_FINAL/4Be9_N_INEL.dat";
string infile06  	= pathData+"all_FINAL/4Be9_N_RAD.dat";
string infile07  	= pathData+"all_FINAL/5B10_N_EL.dat";
string infile08  	= pathData+"all_FINAL/5B10_N_INEL.dat";
string infile09  	= pathData+"all_FINAL/5B10_N_RAD.dat";
string infile10  	= pathData+"all_FINAL/5B11_N_EL.dat";
string infile11  	= pathData+"all_FINAL/5B11_N_INEL.dat";
string infile12  	= pathData+"all_FINAL/5B11_N_RAD.dat";
string infile13  	= pathData+"all_FINAL/64Gd154_N_EL.dat";
string infile14  	= pathData+"all_FINAL/64Gd154_N_INEL.dat";
string infile15  	= pathData+"all_FINAL/64Gd154_N_RAD.dat";
string infile16  	= pathData+"all_FINAL/64Gd155_N_EL.dat";
string infile17  	= pathData+"all_FINAL/64Gd155_N_INEL.dat";
string infile18  	= pathData+"all_FINAL/64Gd155_N_RAD.dat";
string infile19  	= pathData+"all_FINAL/64Gd156_N_EL.dat";
string infile20  	= pathData+"all_FINAL/64Gd156_N_INEL.dat";
string infile21  	= pathData+"all_FINAL/64Gd156_N_RAD.dat";
string infile22  	= pathData+"all_FINAL/64Gd157_N_EL.dat";
string infile23  	= pathData+"all_FINAL/64Gd157_N_INEL.dat";
string infile24  	= pathData+"all_FINAL/64Gd157_N_RAD.dat";
string infile25  	= pathData+"all_FINAL/64Gd158_N_EL.dat";
string infile26  	= pathData+"all_FINAL/64Gd158_N_INEL.dat";
string infile27  	= pathData+"all_FINAL/64Gd158_N_RAD.dat";
string infile28  	= pathData+"all_FINAL/64Gd160_N_EL.dat";
string infile29  	= pathData+"all_FINAL/64Gd160_N_INEL.dat";
string infile30  	= pathData+"all_FINAL/64Gd160_N_RAD.dat";
string infile31  	= pathData+"all_FINAL/6C12_N_EL.dat";
string infile32  	= pathData+"all_FINAL/6C12_N_INEL.dat";
string infile33  	= pathData+"all_FINAL/6C12_N_RAD.dat";
string infile34  	= pathData+"all_FINAL/8O16_N_EL.dat";
string infile35  	= pathData+"all_FINAL/8O16_N_INEL.dat";
string infile36  	= pathData+"all_FINAL/8O16_N_RAD.dat";
string infile37  	= pathData+"all_FINAL/92U234_N_EL.dat";
string infile38  	= pathData+"all_FINAL/92U234_N_F.dat";
string infile39  	= pathData+"all_FINAL/92U234_N_INEL.dat";
string infile40  	= pathData+"all_FINAL/92U234_N_RAD.dat";
string infile41  	= pathData+"all_FINAL/92U235_N_EL.dat";
string infile42  	= pathData+"all_FINAL/92U235_N_F.dat";
string infile43  	= pathData+"all_FINAL/92U235_N_INEL.dat";
string infile44  	= pathData+"all_FINAL/92U235_N_RAD.dat";
string infile45  	= pathData+"all_FINAL/92U238_N_EL.dat";
string infile46  	= pathData+"all_FINAL/92U238_N_F.dat";
string infile47  	= pathData+"all_FINAL/92U238_N_INEL.dat";
string infile48  	= pathData+"all_FINAL/92U238_N_RAD.dat";
string infile49  	= pathData+"all_FINAL/1H1_N_EL.dat";
string infile50  	= pathData+"all_FINAL/1H1_N_INEL.dat";

string outfile00 	= "EBIN_count";
string outfile01	= "output";


//	file pointers
ifstream iFile;
ofstream oFile1, oFile2;


// declare variable
//	CONSTANTS
const double kg2gm = 1000.0;
const double atomicMassU = 1.660538921e-27 * kg2gm;   // convert kg to gm
const double AvoBarn = 0.6022;	// Avogadro no. multiplied by Barn.cm
const double nAvo = 6.022e+23;	// Avogadro no.
const double Massneu = 1.674927351e-27 * kg2gm; 	// mass of n in gms
const double MD = 2.01410177784;			// mass of D g/mol
const double MO = 15.99943;				// mass of O g/mol
const double MGd = 157.25;				// mass of Gd152
const double MU = 238.02891;				// mass of U235
const double Ethermal = 25.0e-3;		// energy of thermal neutron
const double barn2m = 1.0e-28;			// 1b = 1e-28 m^2
const double barn2cm = barn2m * 1.0e+4; 	// 1b = 1e-24 cm^2
const long maxf =10005;
const int maxmat = 40;
const int maxcomp = 5;
const long maxbin = 50;

//	List of Structures
nuclides MatList[maxmat];


//	integers
long int totcol = 0, col2 = 0, nthermal = 0, nlost = 0, nabsGd= 0, nfiss = 0;
long int nparticles;
int itype, whichGd; // 1 = fission, 2 = absorption by 64Gd152
int whichmoderator, whichmat; // 15 = D2O, 16 = graphite, 17 = Natural gd, 18 = natural U, 19 = Natural Boron;
long maxH1E=0, maxH1I=0,maxH2E=0, maxH2I=0,maxBe9N2N=0, maxBe9E=0, maxBe9I=0, maxBe9R=0;
long maxB11R=0, maxB11I=0, maxB11E=0, maxB10R=0, maxB10I=0, maxB10E=0;
long maxGd155R=0, maxGd155I=0, maxGd155E=0, maxGd154R=0, maxGd154I=0, maxGd154E=0;
long maxGd157R=0, maxGd157I=0, maxGd157E=0, maxGd156R=0, maxGd156I=0, maxGd156E=0;
long maxGd160R=0, maxGd160I=0, maxGd160E=0, maxGd158R=0, maxGd158I=0, maxGd158E=0;
long maxO16R=0, maxO16I=0, maxO16E=0, maxC12R=0, maxC12I=0, maxC12E=0;
long maxU234R=0, maxU234I=0, maxU234F=0, maxU234E=0;
long maxU235R=0, maxU235I=0, maxU235F=0, maxU235E=0;
long maxU238R=0, maxU238I=0, maxU238F=0, maxU238E=0;
int binstep, totbin, expo1, expo2;

//	stringS & STRINGS
ostringstream str11, str12, strd, strm, stry, strtype;
char tempf[100];

//	Doubles
double nH1, nBe9, nD_D2O, nO_D2O, nU234, nU235, nU238, nGd154, nGd155, nH1Pl, nC12Pl;
double nGd156, nGd157, nGd158, nGd160, nc12, nC12gr, nB10nat, nB11nat, nB10, nB11;
double TotalGd_slab_wt;
double enrich, Gdcon, Unatwt;
double En0;
double Rsource, Rscatter, Delrscatter, Rreflector, Treflector;
double nbin[maxbin][2];
double enH1E[maxf], xH1E[maxf], enH1I[maxf], xH1I[maxf];
double enH2E[maxf], xH2E[maxf], enH2I[maxf], xH2I[maxf];
double enBe9N2N[maxf], xBe9N2N[maxf], enBe9E[maxf], xBe9E[maxf];
double enBe9I[maxf], xBe9I[maxf], enBe9R[maxf], xBe9R[maxf];
double enB10E[maxf], xB10E[maxf], enB10I[maxf], xB10I[maxf], enB10R[maxf], xB10R[maxf];
double enB11E[maxf], xB11E[maxf], enB11I[maxf], xB11I[maxf], enB11R[maxf], xB11R[maxf];
double enGd154E[maxf], xGd154E[maxf], enGd154I[maxf], xGd154I[maxf], enGd154R[maxf], xGd154R[maxf];
double enGd155E[maxf], xGd155E[maxf], enGd155I[maxf], xGd155I[maxf], enGd155R[maxf], xGd155R[maxf];
double enGd156E[maxf], xGd156E[maxf], enGd156I[maxf], xGd156I[maxf], enGd156R[maxf], xGd156R[maxf];
double enGd157E[maxf], xGd157E[maxf], enGd157I[maxf], xGd157I[maxf], enGd157R[maxf], xGd157R[maxf];
double enGd158E[maxf], xGd158E[maxf], enGd158I[maxf], xGd158I[maxf], enGd158R[maxf], xGd158R[maxf];
double enGd160E[maxf], xGd160E[maxf], enGd160I[maxf], xGd160I[maxf], enGd160R[maxf], xGd160R[maxf];
double enC12E[maxf], xC12E[maxf], enC12I[maxf], xC12I[maxf], enC12R[maxf], xC12R[maxf];
double enO16E[maxf], xO16E[maxf], enO16I[maxf], xO16I[maxf], enO16R[maxf], xO16R[maxf];
double enU234E[maxf], xU234E[maxf], enU234F[maxf], xU234F[maxf], enU234I[maxf], xU234I[maxf], enU234R[maxf], xU234R[maxf];
double enU235E[maxf], xU235E[maxf], enU235F[maxf], xU235F[maxf], enU235I[maxf], xU235I[maxf], enU235R[maxf], xU235R[maxf];
double enU238E[maxf], xU238E[maxf], enU238F[maxf], xU238F[maxf], enU238I[maxf], xU238I[maxf], enU238R[maxf], xU238R[maxf];

//	Declare function prototypes
double max3(double, double, double);
double min3(double, double, double);
double funcinterp(double, double, double, double);
double interpol(double, double [maxf], double [maxf], int);
void reflection(double, double, double, double, double, double, double);
void scatter(struct particles *, struct nuclides *, struct nuclides *, struct nuclides *, int, long, bool);
void GetDC(bool, double &, double &, double &, double);
void initneutron(struct particles *, bool, bool);
void initscat(struct nuclides *, struct nuclides *, struct nuclides *,struct nuclides *, double, double);
void inputdata(void);
void outputdata(long int, long int, long int, long int, long int, long int, long int, long int, double);
void readscatdata(void);
void pause(void);
