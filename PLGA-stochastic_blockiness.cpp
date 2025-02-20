/*****************************************************************************************
	
	Written by Nate Lynd
	This program simulates a copolymerization with reversibility
	No cycles here

*****************************************************************************************/
#include <cmath>
#include <string>
#include <random>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <Eigen/Core>
#include <stdexcept>
#undef eigen_assert
#define eigen_assert(x) \
  if (!(x)) { throw (std::runtime_error("Put your message here")); }

#include "polymer.h"

using namespace std;
using namespace Eigen;

/*****************************************************************************************

	
	
*****************************************************************************************/

int		selectNewChain(int selected_chain, int chains);
void	MWD_report(vector<polymer> polymers,string distribution);
void	compositionDistribution_report(vector<polymer> polymers,string composition);
void    matrix_report(vector<polymer> polymers, string matrix, int chains, double volume);
void	maldi_model(vector<polymer> polymers, string maldi);
double	conversion(int A,int nA0,int B,int nB0);

int		totalDyads(vector<polymer> polymers);
int		totalTriads(vector<polymer> polymers);
int		totalTetrads(vector<polymer> polymers);
int		totalPentads(vector<polymer> polymers);


int main(int argc, char* argv[])
{
 	string		kinetics_file;
	string      sequence_file;
 	string		distribution;
 	string		composition;
 	string		maldi;
	string      matrix;
	string      chaintest;
	ofstream	kinetics[2];
		
	const double	avogadros = 6.022149e+23; // no. per mole.
	
	random_device rd;
	double      rand = rd();
	mt19937		rando(rd()); //30817122
	double		rando_maximo = 4294967295;
	double		a,b,c;
	int			selected_chain;
	int			selected_ester=0;
	int			ester_count=0;
	int			ester_id = 0;
	int			chain1=0;
	int			chain2=0;
	int			monomer1=0;
	int			monomer2=0;
	vector<int> chain1sequence;
	vector<int> chain2sequence;
	int			chain1length;
	int			chain2length;
	
	double		I0,A0,B0;
	double		cI,cA,cB,cP,cPA,cPB,cPAA,cPAB,cPBA,cPBB,cPAAA,cPAAB,cPABA,cPABB,cPBBB,cPBBA,cPBAB,cPBAA,cE;
	int			nP=0,nPA=0,nPB=0,nPAA=0,nPAB=0,nPBA=0,nPBB=0,nPAAA=0,nPABB=0,nPABA=0,nPAAB=0,nPBAA=0,nPBBB=0,nPBBA=0,nPBAB=0,nE=0;
	int			nAA, nAB, nBA, nBB, nIM;//nP is total number of polymers. nIM is total number of linear chains that have initiated
	double		kA,kB;
	double		kAA,kAB,kBA,kBB,kAmA,kAmB,kBmA,kBmB,kT; //regular rate constants (units depend on rate equation)
	double		kA_d, kB_d, kAA_d, kAB_d, kBA_d, kBB_d, kAmA_d, kAmB_d, kBmA_d, kBmB_d, kT_d; //discretized rate constants (units are all in seconds)
	double		termination_time;
	double		cat_ratio; //molar ratio of catalyst to initiator ([cat]/[in])
	double		volume;
	double		time = 0.;
	double		time_report = 0.;
	double		time_report2 = 0.;
	double		time_report_int=1.;
	double		time_report_int2 = 15.;
	int			report2_count = 0;
	
	int			counter = 0;
	int			chains; //This is the number of initiators / number of linear chains
	int			I,A,B;
	int			nI0,nA0,nB0;
	int			Atype = 1, Btype = 2;
	
	vector<int>			degreeOfPolymerization;
	
	vector<int>			initiated(0);
	
	vector<double>		r(15);
	vector<double>		lower(15),upper(15);
	double				r_tot;
	int					n_tot;
	
	vector<double>	dyad_fractions(4);
	vector<double>	triad_fractions(8);
	vector<double>	tetrad_fractions(16);
	vector<double>  pentad_fractions(32);
	nAA = 0; nAB = 0; nBA = 0; nBB = 0; nIM = 0;
	
	const int	iI=0,iA=1,iB=2,iXA=3,iXB=4; // enumeration

	// read input from command line:
	// kinetics file for structured data output:
	kinetics_file				=	string(argv[1]);
	// initiator molecules = number of chains
	I0							=	atof(argv[2]);
	// A monomers
	A0							=	atof(argv[3]);
	// B monomer
	B0							=	atof(argv[4]);
	// forward rate constants:
	kAA							=	atof(argv[5]);
	kAB							=	atof(argv[6]);
	kBA							=	atof(argv[7]);
	kBB							=	atof(argv[8]);
	// backward rate constants
	kAmA						=	atof(argv[9]);
	kAmB						=	atof(argv[10]);
	kBmA						=	atof(argv[11]);
	kBmB						=	atof(argv[12]);
	//transesterification rate constants
	kT							=	atof(argv[13]);
	// other stuff
	termination_time			=	atof(argv[14]);
	chains						=	atof(argv[15]);
	cat_ratio					=	atof(argv[16]);

	vector<polymer>		polymers(chains);//linear polymer chains

	time_report					=	0.;
	time_report2				=	0.;
	time_report_int				=	0.5;
	time_report_int2			=	15.;

	cout << showpoint;
	// for simplicity, we set kA = kAA and kB = kBB
	kA = kAA;
	kB = kBB;
	cout << "rA*rB = " << kAA/kAB*kBB/kBA << endl;
	cout << "kA = " << kA << endl;
	cout << "kB = " << kB << endl;
	cout << "kAA = " << kAA << endl;
	cout << "kAB = " << kAB << endl;
	cout << "rA = " << kAA/kAB << endl;
	cout << "kBB = " << kBB << endl;
	cout << "kBA = " << kBA << endl;
	cout << "rB = " << kBB/kBA << endl;
	cout << "kA-A = " << kAmA << endl;
	cout << "kA-B = " << kAmB << endl;
	cout << "r-A = " << kAmA/kAmB << endl;
	cout << "kB-B = " << kBmB << endl;
	cout << "kB-A = " << kBmA << endl;
	cout << "kT = " << kT << endl;
	cout << "I0 = " << I0 << endl;
	cout << "A0 = " << A0 << endl;
	cout << "B0 = " << B0 << endl;
	cout << "DP = " << (A0+B0)/I0 << endl;
	cout << fixed << "rand= " << rand << endl;


	kinetics_file = kinetics_file + ".dat";
	sequence_file = string(argv[1]) + ".seq";

	
	// The system volume will be based on the initiator concentration (number of chains)
	volume = double(chains)/avogadros/I0;
	// set up initial conditions
	nI0 = I = chains; // I(0)
	nA0 = A = int(A0*volume*avogadros); // A(0)
	nB0 = B = int(B0*volume*avogadros); // B(0)
	
	cout << "System volume is " << volume << " L" << endl;
	// this will set the number of monomers
	cout << "Number of I initiators: " << I << endl;
	cout << "Number of A monomers: " << A << endl;
	cout << "Number of B monomers: " << B << endl;
	cout << "Number of potential polymers: " << polymers.size() << endl;

	//Calculate discretized rate constants
	kA_d	=		(kA/volume)/avogadros;
	kB_d	=		(kB/volume)/avogadros;
	kAA_d	=		(kAA/volume)/avogadros;
	kAB_d	=		(kAB/volume)/avogadros;
	kBA_d	=		(kBA/volume)/avogadros;
	kBB_d	=		(kBB/volume)/avogadros;
	kAmA_d	=		kAmA;
	kAmB_d	=		kAmB;
	kBmA_d	=		kBmA;
	kBmB_d	=		kBmB;
	kT_d	=		(kT/volume)/avogadros;
	
	kinetics[0].open(kinetics_file.c_str(),ios::out);
	kinetics[0].precision(14);
	kinetics[0] << "System volume is " << volume << " L" << endl;
	kinetics[0] << "Number of I initiators: " << I << endl;
	kinetics[0] << "Number of A monomers: " << A << endl;
	kinetics[0] << "Number of B monomers: " << B << endl;
	kinetics[0] << "Number of potential polymers: " << polymers.size() << endl;
	kinetics[0] << "rA*rB = " << kAA/kAB*kBB/kBA << endl;
	kinetics[0] << "kA = " << kA << endl;
	kinetics[0] << "kB = " << kB << endl;
	kinetics[0] << "kAA = " << kAA << endl;
	kinetics[0] << "kAB = " << kAB << endl;
	kinetics[0] << "rA = " << kAA/kAB << endl;
	kinetics[0] << "kBB = " << kBB << endl;
	kinetics[0] << "kBA = " << kBA << endl;
	kinetics[0] << "rB = " << kBB/kBA << endl;
	kinetics[0] << "kA-A = " << kAmA << endl;
	kinetics[0] << "kA-B = " << kAmB << endl;
	kinetics[0] << "r-A = " << kAmA/kAmB << endl;
	kinetics[0] << "kB-B = " << kBmB << endl;
	kinetics[0] << "kB-A = " << kBmA << endl;
	kinetics[0] << "kT = " << kT << endl;
	kinetics[0] << "I0 = " << I0 << endl;
	kinetics[0] << "A0 = " << A0 << endl;
	kinetics[0] << "B0 = " << B0 << endl;
	kinetics[0] << "DP = " << (A0+B0)/I0 << endl;
	kinetics[0] << fixed << "rand= " << rand << endl;
	kinetics[0] << "volume =" << volume << endl;
	kinetics[0] << "catalyst:initiator ratio =" << cat_ratio << endl;
	kinetics[0] << "time\tI\tA\tB\tPA\tPB\tconv."<< endl;
	cout << "time\tI\tA\tB\tPA\tPB\tconv."<< endl;

	kinetics[1].open(sequence_file.c_str(),ios::out);
	kinetics[1].precision(14);
	kinetics[1] << "time"<<"\t"<< "AA" << "\t" << "AB" << "\t" << "BA" << "\t" << "BB"<<"\t"
		<<"AAA"<<"\t"<<"AAB"<<"\t"<<"ABA"<<"\t"<<"ABB"<<"\t"<<"BAA"<<"\t"<<"BAB"<<"\t"<<"BBA"<<"\t"<<"BBB"<<"\t"
		<<"AAAA"<<"\t"<<"AAAB"<<"\t"<<"AABA"<<"\t"<<"AABB"<<"\t"<<"ABAA"<<"\t"<<"ABAB"<<"\t"<<"ABBA"<<"\t"<<"ABBB"<<"\t"<<"BAAA"<<"\t"<<"BAAB"<<"\t"<<"BABA"<<"\t"<<"BABB"<<"\t"<<"BBAA"<<"\t"<<"BBAB"<<"\t"<<"BBBA"<<"\t"<<"BBBB"<<"\t"
		<<"AAAAA"<<"\t"<<"AAAAB"<<"\t"<<"AAABA"<<"\t"<<"AAABB"<<"\t"<<"AABAA"<<"\t"<<"AABAB"<<"\t"<<"AABBA"<<"\t"<<"AABBB"<<"\t"<<"ABAAA"<<"\t"<<"ABAAB"<<"\t"<<"ABABA"<<"\t"<<"ABABB"<<"\t"<<"ABBAA"<<"\t"<<"ABBAB"<<"\t"<<"ABBBA"<<"\t"
		<<"ABBBB"<<"\t"<<"BAAAA"<<"\t"<<"BAAAB"<<"\t"<<"BAABA"<<"\t"<<"BAABB"<<"\t"<<"BABAA"<<"\t"<<"BABAB"<<"\t"<<"BABBA"<<"\t"<<"BABBB"<<"\t"<<"BBAAA"<<"\t"<<"BBAAB"<<"\t"<<"BBABA"<<"\t"<<"BBABB"<<"\t"<<"BBBAA"<<"\t"<<"BBBAB"<<"\t"
		<<"BBBBA"<<"\t"<<"BBBBB"<<"\t"<<"repeat units" << "\t" <<"dyads" <<"\t"<<"triads"<<"\t"<<"tetrads"<<"\t"<<"pentads"<<"\t""n G's" <<"\t"<<"n GG's" << "\t"<<"n GGG's" << "\t" << "n GGGG's" << endl;

	// calculate concentrations for reporting:
	cI	=	double(I)/avogadros/volume;
	cA	=	double(A)/avogadros/volume;
	cB	=	double(B)/avogadros/volume;
	
	while (time < termination_time)
	{
		nP = 0; nPA = 0; nPB = 0; nPAA = 0; nPAB = 0; nPBA = 0; nPBB = 0; nPAAA = 0, nPABB = 0, nPABA = 0, nPAAB = 0, nPBAA = 0, nPBBB = 0, nPBBA = 0, nPBAB = 0; nE = 0; nIM = 0;
		//This starts over counting end group dyads and triads every time
		for (int i = 0; i <= polymers.size() - 1; i++)
		{
			if (polymers[i].getDegreeOfPolymerization() >= 1)
			{
				nIM += 1;
				nP += 1;
				nE += polymers[i].getDegreeOfPolymerization();

				if (polymers[i].endGroup() == 1)
				{
					nPA += 1;
					if (polymers[i].endGroupDyad() == 11)
					{
						nPAA += 1;
						if (polymers[i].endGroupTriad() == 111)
						{
							nPAAA += 1;
						}
						else if (polymers[i].endGroupTriad() == 211)
						{
							nPBAA += 1;
						}
					}
					else if (polymers[i].endGroupDyad() == 21)
					{
						nPBA += 1;
						if (polymers[i].endGroupTriad() == 121)
						{
							nPABA += 1;
						}
						else if (polymers[i].endGroupTriad() == 221)
						{
							nPBBA += 1;
						}
					}
				}
				else if (polymers[i].endGroup() == 2)
				{
					nPB += 1;
					if (polymers[i].endGroupDyad() == 22)
					{
						nPBB += 1;
						if (polymers[i].endGroupTriad() == 122)
						{
							nPABB += 1;
						}
						else if (polymers[i].endGroupTriad() == 222)
						{
							nPBBB += 1;
						}
					}
					else if (polymers[i].endGroupDyad() == 12)
					{
						nPAB += 1;
						if (polymers[i].endGroupTriad() == 112)
						{
							nPAAB += 1;
						}
						else if (polymers[i].endGroupTriad() == 212)
						{
							nPBAB += 1;
						}
					}
				}
			}
		}

		//cout << "nE= " << nE << endl;
		//cout << "nAA, nAB, nBA, nBB, nIM are " << nAA << " " << nAB << " " << nBA << " " << nBB << " " << nIM << endl;

		//update nAA, nAB, nBA, nBB with each reaction instead of counting constantly
		// reaction rates- calculated based on numbers of molecules
		// units of molecules/time
		r[0] = kA_d * I *cat_ratio * A;
		r[1] = kB_d * I * cat_ratio * B;
		r[2] = kAA_d * nPA * cat_ratio * A;
		r[3] = kAB_d * nPA * cat_ratio * B;
		r[4] = kBA_d * nPB * cat_ratio * A;
		r[5] = kBB_d * nPB * cat_ratio * B;
		r[6] = kAmA_d * nPAAA * cat_ratio;
		r[7] = kAmB_d * nPABB * cat_ratio;
		r[8] = kBmA_d * nPBAA * cat_ratio;
		r[9] = kBmB_d * nPBBB * cat_ratio;
		r[10] = kT_d * nE * nI0 * cat_ratio; //number esters * number ends
		r_tot = r[0] + r[1] + r[2] + r[3] + r[4] + r[5] + r[6] + r[7] + r[8] + r[9] + r[10];;

		// update time
		a = double(rando()) / rando_maximo; // time
		b = double(rando()) / rando_maximo; // reaction
		c = double(rando()) / rando_maximo; // chain proposal
		selected_chain = int(double(chains) * c);
		//cout << "selected_chain initial is " <<selected_chain << endl;
		

		time = time - log(a) /r_tot;  

		// calculate limits for reaction decision:
		lower[0] = 0.;
		for (int i = 0; i <= 9; i++)
		{
			lower[i + 1] = lower[i] + r[i] / r_tot;
			upper[i] = lower[i + 1];
		}
		upper[10] = 1.;

		// make decision
		if (b >= lower[0] && b <= upper[0])
		{
			//cout << "reaction 1!" << endl;
			// reaction 1 initiation with A monomer
			// which chain (that has DP = 0)
			if (polymers[selected_chain].hasInitiated() == false)
			{
				// do reaction 
				I = I - 1;
				A = A - 1;
				nAA = nAA + 1;
				//Do we also need to do nIM? I don't think so because we count it elsewhere
				polymers[selected_chain].addMonomerToEnd(Atype);
				polymers[selected_chain].addMonomerToEnd(Atype); //twice because each A monomer is actually a dimer
				initiated.push_back(selected_chain); // this keeps track of initiated chains
				//polymers[selected_chain].printPolymer();
			}
			else // find new chain
			{
				while (polymers[selected_chain].hasInitiated() == true)
				{
					selected_chain = int(double(chains) * double(rando()) / rando_maximo);
					//polymers[selected_chain].printPolymer();
					//cout << "\t1 selected_chain = " << selected_chain << endl;					
					if (polymers[selected_chain].hasInitiated() == false)
					{
						// do reaction
						I = I - 1;
						A = A - 1;
						nAA = nAA + 1;
						polymers[selected_chain].addMonomerToEnd(Atype);
						polymers[selected_chain].addMonomerToEnd(Atype);
						initiated.push_back(selected_chain);
						break;
					}
				}
			}
		}
		else if (b > lower[1] && b <= upper[1])
		{
			//cout << "reaction 2!" << endl;
			//reaction 2- Initiation with B
			if (polymers[selected_chain].hasInitiated() == false)
			{
				// do reaction
				I = I - 1;
				B = B - 1;
				nBB = nBB + 1;
				polymers[selected_chain].addMonomerToEnd(Btype);
				polymers[selected_chain].addMonomerToEnd(Btype);
				initiated.push_back(selected_chain);
			}
			else // find new chain
			{
				while (polymers[selected_chain].hasInitiated() == true)
				{
					selected_chain = int(double(chains) * double(rando()) / rando_maximo);
					//polymers[selected_chain].printPolymer();
					//cout << "\t2 selected_chain = " << selected_chain << endl;
					if (polymers[selected_chain].hasInitiated() == false)
					{
						// do reaction
						I = I - 1;
						B = B - 1;
						nBB = nBB + 1;
						polymers[selected_chain].addMonomerToEnd(Btype);
						polymers[selected_chain].addMonomerToEnd(Btype);
						initiated.push_back(selected_chain);
						break;
					}
				}
			}
		}

		// from here only from previously initiated chains
		else if (b > lower[2] && b <= upper[2])
		{
			// reaction 3  kAA*cPA*cA;
			//cout << "reaction 3!" << endl;
			if (polymers[selected_chain].endGroup() == Atype)
			{
				// do reaction
				A = A - 1;
				nAA = nAA + 2;
				polymers[selected_chain].addMonomerToEnd(Atype);
				polymers[selected_chain].addMonomerToEnd(Atype);
			}
			else // find new chain
			{
				while (polymers[selected_chain].endGroup() != Atype)
				{
					selected_chain = int(double(chains) * double(rando()) / rando_maximo);
					//polymers[selected_chain].printPolymer();
					//cout << "\t3 selected_chain = " << selected_chain << endl;
					if (polymers[selected_chain].endGroup() == Atype)
					{
						// do reaction
						A = A - 1;
						nAA = nAA + 2;
						polymers[selected_chain].addMonomerToEnd(Atype);
						polymers[selected_chain].addMonomerToEnd(Atype);
						//cout << "\t\t reaction executed!" << endl;
						//polymers[selected_chain].printPolymer();
						break;
					}
				}
			}
		}
		else if (b > lower[3] && b <= upper[3])
		{
			//cout << "reaction 4!" << endl;
			//reaction 4	kAB*cPA*cB;
			if (polymers[selected_chain].endGroup() == Atype)
			{
				// do reaction
				B = B - 1;
				nAB = nAB + 1; nBB = nBB + 1;
				polymers[selected_chain].addMonomerToEnd(Btype);
				polymers[selected_chain].addMonomerToEnd(Btype);
			}
			else // find new chain
			{
				while (polymers[selected_chain].endGroup() != Atype)
				{
					selected_chain = int(double(chains) * double(rando()) / rando_maximo);
					//polymers[selected_chain].printPolymer();
					//cout << "\t4 selected_chain = " << selected_chain << endl;
					if (polymers[selected_chain].endGroup() == Atype)
					{
						// do reaction
						B = B - 1;
						nAB = nAB + 1; nBB = nBB + 1;
						polymers[selected_chain].addMonomerToEnd(Btype);
						polymers[selected_chain].addMonomerToEnd(Btype);
						break;
					}
				}
			}
		}
		else if (b > lower[4] && b <= upper[4])
		{
			// reaction 5	kBA*cPB*cA;
			//cout << "reaction 5!" << endl;
			if (polymers[selected_chain].endGroup() == Btype)
			{
				// do reaction
				A = A - 1;
				nBA = nBA + 1; nAA = nAA + 1;
				polymers[selected_chain].addMonomerToEnd(Atype);
				polymers[selected_chain].addMonomerToEnd(Atype);
			}
			else // find new chain
			{
				while (polymers[selected_chain].endGroup() != Btype)
				{
					selected_chain = int(double(chains) * double(rando()) / rando_maximo);
					//polymers[selected_chain].printPolymer();
					//cout << "\t5 selected_chain = " << selected_chain << endl;
					if (polymers[selected_chain].endGroup() == Btype)
					{
						// do reaction
						A = A - 1;
						nBA = nBA + 1; nAA = nAA + 1;
						polymers[selected_chain].addMonomerToEnd(Atype);
						polymers[selected_chain].addMonomerToEnd(Atype);
						break;
					}
				}
			}
		}
		else if (b > lower[5] && b <= upper[5])
		{
			// reaction 6	kBB*cPB*cB
			//cout << "reaction 6!" << endl;
			if (polymers[selected_chain].endGroup() == Btype)
			{
				// do reaction
				B = B - 1;
				nBB = nBB + 2;
				polymers[selected_chain].addMonomerToEnd(Btype);
				polymers[selected_chain].addMonomerToEnd(Btype);
			}
			else if (polymers[selected_chain].hasInitiated()) // find new chain
			{
				while (polymers[selected_chain].endGroup() != Btype)
				{
					selected_chain = int(double(chains) * double(rando()) / rando_maximo);
					//polymers[selected_chain].printPolymer();
					//cout << "\t6 selected_chain = " << selected_chain << endl;
					if (polymers[selected_chain].endGroup() == Btype)
					{
						// do reaction
						B = B - 1;
						nBB = nBB + 2;
						polymers[selected_chain].addMonomerToEnd(Btype);
						polymers[selected_chain].addMonomerToEnd(Btype);
						break;
					}
				}
			}
		}
		else if (b > lower[6] && b <= upper[6])
		{
			// reaction 7	kAmA*cPAAA
			//cout << "reaction 7!" << endl;
			if (polymers[selected_chain].endGroupTriad() == 111)
			{
				// do reaction
				A = A + 1;
				nAA = nAA - 2;
				polymers[selected_chain].removeMonomerFromEnd();
				polymers[selected_chain].removeMonomerFromEnd();
			}
			else // find new chain
			{
				while (polymers[selected_chain].endGroupTriad() != 111) // this could be an issue
				{
					selected_chain = int(double(chains) * double(rando()) / rando_maximo);
					//polymers[selected_chain].printPolymer();
					//cout << "7 selected_chain = " << selected_chain << endl;
					if (polymers[selected_chain].endGroupTriad() == 111)
					{
						// do reaction
						A = A + 1;
						nAA = nAA - 2;
						polymers[selected_chain].removeMonomerFromEnd();
						polymers[selected_chain].removeMonomerFromEnd();
						break;
					}
				}
			}
		}
		else if (b > lower[7] && b <= upper[7])
		{
			// reaction 8	kAmB*cPABB
			//cout << "reaction 8!" << endl;
			if (polymers[selected_chain].endGroupTriad() == 122)
			{
				// do reaction
				B = B + 1;
				nAB = nAB - 1; nBB = nBB - 1;
				polymers[selected_chain].removeMonomerFromEnd();
				polymers[selected_chain].removeMonomerFromEnd();
			}
			else // find new chain
			{
				while (polymers[selected_chain].endGroupTriad() != 122)
				{
					selected_chain = int(double(chains) * double(rando()) / rando_maximo);
					//polymers[selected_chain].printPolymer();
					//cout << "8 selected_chain = " << selected_chain << endl;
					if (polymers[selected_chain].endGroupTriad() == 122)
					{
						// do reaction
						B = B + 1;
						nAB = nAB - 1; nBB = nBB - 1;
						polymers[selected_chain].removeMonomerFromEnd();
						polymers[selected_chain].removeMonomerFromEnd();
						break;
					}
				}
			}
		}
		else if (b > lower[8] && b <= upper[8])
		{
			// reaction 9	kBmA*cPBAA
			//cout << "reaction 9!" << endl;
			if (polymers[selected_chain].endGroupTriad() == 211)
			{
				// do reaction
				A = A + 1;
				nBA = nBA - 1; nAA = nAA - 1;
				polymers[selected_chain].removeMonomerFromEnd();
				polymers[selected_chain].removeMonomerFromEnd();
			}
			else // find new chain
			{
				while (polymers[selected_chain].endGroupTriad() != 211)
				{
					selected_chain = int(double(chains) * double(rando()) / rando_maximo);
					//polymers[selected_chain].printPolymer();
					//cout << "9 selected_chain = " << selected_chain << endl;
					if (polymers[selected_chain].endGroupTriad() == 211)
					{
						// do reaction
						A = A + 1;
						nBA = nBA - 1; nAA = nAA - 1;
						polymers[selected_chain].removeMonomerFromEnd();
						polymers[selected_chain].removeMonomerFromEnd();
						break;
					}
				}
			}
		}
		else if (b > lower[9] && b <= upper[9])
		{
			// reaction 10	kBmB*cPBBB
			//cout << "reaction 10!" << endl;
			if (polymers[selected_chain].endGroupTriad() == 222)
			{
				// do reaction
				B = B + 1;
				nBB = nBB - 2;
				polymers[selected_chain].removeMonomerFromEnd();
				polymers[selected_chain].removeMonomerFromEnd();
			}
			else // find new chain
			{
				while (polymers[selected_chain].endGroupTriad() != 222)
				{
					selected_chain = int(double(chains) * double(rando()) / rando_maximo);
					//polymers[selected_chain].printPolymer();
					//cout << "10 selected_chain = " << selected_chain << endl;
					if (polymers[selected_chain].endGroupTriad() == 222)
					{
						// do reaction
						B = B + 1;
						nBB = nBB - 2;
						polymers[selected_chain].removeMonomerFromEnd();
						polymers[selected_chain].removeMonomerFromEnd();
						break;
					}
				}
			}
		}

		//Transesterification: This probably should be more complicated than it is
		//Random transesterification for now
		else if (b > lower[10] && b <= upper[10])
		{
			//cout << "reaction 11, AA transesterification!" << endl;
			//Chain 1 is the chain which is being attacked and chain 2 is the attacking chain
			chain2 = selected_chain;
			selected_ester = ceil(double(nE) * double(rando()) / rando_maximo); //only check AA esters
			ester_count = 0;
			for (int i = 0; i <= polymers.size() - 1; i++) //go through polymers
			{
				for (int j = 0; j <= polymers[i].getDegreeOfPolymerization() - 2; j++)
				{

					ester_count += 1; // add one to ester count if correct ester
					//check if we've reached selected ester
					if (ester_count == selected_ester)
					{
						//store selected ester
						ester_id = j;
						chain1 = i;
						//cout << "before break, ester_count and selected ester are" << ester_count << " " << selected_ester << endl;
						break;
					}
				}
			}

			//Double check that we won't make a cycle- too complex
			//keep in mind dyad counts here are not going to be accurate
			if (chain2 == chain1)
			{
				while (chain2 == chain1)
				{
					chain2= int(double(chains) * double(rando() / rando_maximo));
					if (chain1 != chain2)
					{
						break;
					}
				}
			}
			//store original sequences
			chain1sequence = polymers[chain1].sequence;
			chain1length = polymers[chain1].getDegreeOfPolymerization();
			chain2sequence = polymers[chain2].sequence;
			chain2length = polymers[chain2].getDegreeOfPolymerization();
			//now we have about 1 billion ifs to go through to actually execute the reaction -- Not anymore!
					if (polymers[chain2].hasInitiated() == false)
					{
						//initiator attacking linear chain
						//cout << "initiator attacking linear chain" << endl;
						//One initiator is consumed
						I = I - 1;
						//Don't fix dyads- counts are wrong until recounted for sequence distribution
						//Remove monomers from end of chain1 and add to chain2
						for (int ii = ester_id+1; ii <= chain1length-1; ii++)
						{
							polymers[chain1].removeMonomerFromEnd();
							if (chain1sequence[ii] == 1)
							{
								polymers[chain2].addMonomerToEnd(Atype);
							}
							else if (chain1sequence[ii] == 2)
							{
								polymers[chain2].addMonomerToEnd(Btype);
							}
						}
						/*	cout << "chain 1 after is " << chain1 << endl;
						polymers[chain1].printPolymer();
						cout << "chain 2 after is " << chain2 << endl;
						polymers[chain2].printPolymer();*/
					}
					else if (polymers[chain2].hasInitiated() == true)
					{
						//linear chain attacking another linear chain
						//cout << "linear chain attacking another linear chain" << endl;
						//Remove monomers from end of chain1 and add to chain2
						for (int ii = ester_id + 1; ii <= chain1length - 1; ii++)
						{
							polymers[chain1].removeMonomerFromEnd();
							if (chain1sequence[ii] == 1)
							{
								polymers[chain2].addMonomerToEnd(Atype);
							}
							else if (chain1sequence[ii] == 2)
							{
								polymers[chain2].addMonomerToEnd(Btype);
							}
						}
						/*cout << "chain 1 after is " << chain1 << endl;
						polymers[chain1].printPolymer();
						cout << "chain 2 after is " << chain2 << endl;
						polymers[chain2].printPolymer();*/
					}
				
				
			
		}
		else
		{
			// throw error
			cout << " oh nonononononon shit shit shit!" << endl;
			cout << a << "\t" << b << "\t" << c << endl;
			cout << "selected chain = " << selected_chain << endl;
			cout << polymers[selected_chain].getDegreeOfPolymerization() << endl;
			cout << lower[0] << lower[1] << lower[2]<<lower[3]<<lower[4]<<lower[5]<<lower[6]<<lower[7]<<lower[8]<<lower[9]<<lower[10]<<upper[10]<<endl;
		}
//Report sequence motifs
		if (time >= time_report2)
		{
			// set all fractions to zero to start
			for (int i = 0; i <= 3; i++)
			{
				dyad_fractions[i] = 0.;
			}
			for (int i = 0; i <= 7; i++)
			{
				triad_fractions[i]=0.;
			}
			for (int i = 0; i <= 15; i++)
			{
				tetrad_fractions[i]=0.;
			}
			for (int i = 0; i <= 31; i++)
			{
				pentad_fractions[i]=0.;
			}
			int k = 0.; //counter for sequence motifs

			//calculate dyad populations
			for (int i = 0; i <= 1; i++)
			{
				for (int ii = 0; ii <= 1; ii++)
				{
					for (int j = 0; j <= polymers.size() - 1; j++)
					{
						dyad_fractions[k]+= polymers[j].numberDyadsOfType(i+1,ii+1);
					}
					k+=1;
				}
			}

			//calculate triad populations
			k=0.;
			for (int i = 0; i <= 1; i++)
			{
				for (int ii = 0; ii <= 1; ii++)
				{
					for (int iii = 0; iii <= 1; iii++)
					{
						for (int j = 0; j <= polymers.size() - 1; j++)
						{
							triad_fractions[k]+=polymers[j].numberTriadsOfType(i+1,ii+1,iii+1);
						}
						k+=1;
					}
				}
			}

			//calculate tetrad populations
			k=0.;
			for (int i = 0; i <= 1; i++)
			{
				for (int ii = 0; ii <= 1; ii++)
				{
					for (int iii = 0; iii <= 1; iii++)
					{
						for (int iiii = 0; iiii <= 1; iiii++)
						{
							for(int j=0;j<=polymers.size()-1;j++)
							{
								tetrad_fractions[k] += polymers[j].numberTetradsOfType(i + 1, ii + 1, iii + 1,iiii+1);
							}
							k += 1;
						}
					}
				}
			}

			//calculated pentad populations
			k = 0.;
			for (int i = 0; i <= 1; i++)
			{
				for (int ii = 0; ii <= 1; ii++)
				{
					for (int iii = 0; iii <= 1; iii++)
					{
						for (int iiii = 0; iiii <= 1; iiii++)
						{
							for (int iiiii = 0; iiiii <= 1; iiiii++)
							{
								for(int j=0; j<=polymers.size()-1;j++)
								{
									pentad_fractions[k] += polymers[j].numberPentadsOfType(i + 1, ii + 1, iii + 1, iiii + 1,iiiii+1);
								}
								k += 1;
							}
						}
					}
				}
			}
			//divide populations by total numbers and report
			kinetics[1]<<time<<"\t";
			for (int i = 0; i <= 3; i++)
			{
				if (totalDyads(polymers) > 0)
				{
					dyad_fractions[i] = dyad_fractions[i] /totalDyads(polymers);
					kinetics[1] << dyad_fractions[i] << "\t";
				}
				else if (totalDyads(polymers) == 0) //need this to avoid nAn's in first line
				{
					dyad_fractions[i] = 0;
					kinetics[1] << dyad_fractions[i] << "\t";
				}
			}
			for (int i = 0; i <= 7; i++)
			{
				if (totalTriads(polymers) > 0)
				{
					triad_fractions[i] = triad_fractions[i] / totalTriads(polymers);
					kinetics[1] << triad_fractions[i] << "\t";
				}
				else if (totalTriads(polymers) == 0)
				{
					triad_fractions[i] = 0;
					kinetics[1] << triad_fractions[i] << "\t";
				}
			}
			for (int i = 0; i <= 15; i++)
			{
				if (totalTetrads(polymers) > 0)
				{
					tetrad_fractions[i] = tetrad_fractions[i] / totalTetrads(polymers);
					kinetics[1] << tetrad_fractions[i] << "\t";
				}
				else if (totalTetrads(polymers) == 0)
				{
					tetrad_fractions[i] = 0;
					kinetics[1] << tetrad_fractions[i] << "\t";
				}
			}
			
			for (int i = 0; i <= 31; i++)
			{
				if (totalPentads(polymers) > 0)
				{
					pentad_fractions[i] = pentad_fractions[i] / totalPentads(polymers);
					kinetics[1] << pentad_fractions[i] << "\t";
				}
				else if (totalPentads(polymers) == 0)
				{
					pentad_fractions[i] = 0;
					kinetics[1] << pentad_fractions[i] << "\t";
				}
			}

			//kinetics[1] << totalDyads(polymers) << "\t" << totalTriads(polymers) << "\t" << totalTetrads(polymers) << "\t" << totalPentads(polymers) << "\t";
			// Also want to report number of G's, number of GG's, GGG's, GGGG's.
			int nrepeats=0, gone = 0, gtwo = 0, gthree = 0, gfour = 0;
			for (int i = 0; i <= polymers.size()-1; i++)
			{
				for (int ii = 0; ii <= polymers[i].getDegreeOfPolymerization()-1; ii++)
				{
					nrepeats += 1;
					if (polymers[i].repeatUnitAtPosition(ii) == 1)
					{
						gone += 1;
					}
				}
				gtwo += polymers[i].numberDyadsOfType(1, 1);
				gthree += polymers[i].numberTriadsOfType(1, 1, 1);
				gfour += polymers[i].numberTetradsOfType(1, 1, 1, 1);
			}
			kinetics[1] << nrepeats<<"\t"<<totalDyads(polymers) << "\t" << totalTriads(polymers) << "\t" << totalTetrads(polymers) << "\t" << totalPentads(polymers) << "\t" << gone << "\t" << gtwo << "\t" << gthree << "\t" << gfour;
			kinetics[1]<<endl;
			//kinetics[1] << nAA << "\t" << nAB << "\t" << nBA << "\t" << nBB << endl;

			if (time >= 60 && time<360)
			{
				time_report_int2 = 60;
			}
			else if (time >= 360)
			{
				time_report_int2 = 360;
			}
			time_report2 = time_report2 + time_report_int2;
		}

	
			cI = double(I) / avogadros / volume;
			cA = double(A) / avogadros / volume;
			cB = double(B) / avogadros / volume;
			cPA = double(nPA) / avogadros / volume;
			cPB = double(nPB) / avogadros / volume;

			if (time >= time_report)
			{
				kinetics[0] << time << "\t" << cI << "\t" << cA << "\t" << cB << "\t" << cPA << "\t" << cPB << "\t" << conversion(A, nA0, B, nB0) << endl;
				cout << time << "\t" << cI << "\t" << cA << "\t" << cB << "\t" << cPA << "\t" << cPB << "\t" << conversion(A, nA0, B, nB0) << endl;

				time_report = time_report + time_report_int;
			}
	
		//Still calculate concentrations to report
					
		if (time > termination_time)
		{
			cout << "time up" << endl;
			break;
		}	
		

		counter += 1;
	}
	cI = double(I) / avogadros / volume;
	cA = double(A) / avogadros / volume;
	cB = double(B) / avogadros / volume;
	cPA = double(nPA) / avogadros / volume;
	cPB = double(nPB) / avogadros / volume;
	kinetics[0] << time << "\t" << cI << "\t" << cA << "\t" << cB << "\t" << cPA << "\t" << cPB << "\t" << conversion(A, nA0, B, nB0) << endl;
	cout << time << "\t" << cI << "\t" << cA << "\t" << cB << "\t" << cPA << "\t" << cPB << "\t" << conversion(A, nA0, B, nB0) << endl;

	kinetics[0].close();
	kinetics[1].close();
	// report out - configure all for Mathematica:
	distribution	= string(argv[1]) + ".mwd";
	composition		= string(argv[1]) + ".composition";
	matrix			= string(argv[1]) + ".matrix";
	maldi			= string(argv[1]) + ".maldi";
	
	//MWD_report(polymers,distribution);
	//compositionDistribution_report(polymers,composition);
	//matrix_report(polymers, matrix, chains, volume);
	//maldi_model(polymers,maldi);
	
	return 0;
}
/****************************************************************************************/
int	selectNewChain(int selected_chain, int chains)
{
	if (selected_chain < (chains - 1))
	{
		return (selected_chain + 1);
	}
	else if (selected_chain == (chains - 1))
	{
		return 0;
	}
	else
	{
		cout << "selectNewChain error at " << selected_chain << endl;
		return 0;
	}
}
/****************************************************************************************/
// molecular weight distribution
void MWD_report(vector<polymer> polymers,string distribution)
{
	// returns M[1] = Mn, and M[2] = Mw
	int total_polymers = polymers.size();
	int	ones, twos;
	int max_degree_of_polymerization = 0;
	ofstream kinetics;
	double Mn = 0.;
	double Mw = 0.; double Nw = 0.;
	double PDI = 0.;
	// need to know the maximum degree of polymerization so that even if there are some 
	// chains that are only 1 or 2, they will not create a segmentation fault
	for (int i=0;i<=total_polymers-1;i++)
	{
		if (max_degree_of_polymerization < polymers[i].getDegreeOfPolymerization())
		{
			max_degree_of_polymerization = polymers[i].getDegreeOfPolymerization();
		}
	}

	MatrixXd Mat = MatrixXd::Zero(max_degree_of_polymerization + 1, max_degree_of_polymerization + 1);
	VectorXd Lin = VectorXd::Zero(max_degree_of_polymerization);

	for (int i = 0; i <= total_polymers - 1; i++)
	{
		ones=0., twos=0.;
		for (int j = 0; j <= polymers[i].getDegreeOfPolymerization() - 1; j++)
		{
			// count all polymers with length j
			if (polymers[i].repeatUnitAtPosition(j) == 1)
			{
				ones += 1;
			}
			else if (polymers[i].repeatUnitAtPosition(j) == 2)
			{
				twos += 1;
			}
		}
		// count chain in right place based on ones and twos
		Mat(ones, twos) += 1.; // chain
	}
	// now that we have this whole matrix counting the number of chains for each i,j
	// we can calculate molecular weight averages
	// Not convinced any of this is right...

	kinetics.open(distribution.c_str(), ios::out);
	kinetics.precision(14);
	// here

	for (int N = 0; N <= max_degree_of_polymerization - 1; N++)
	{
		for (int i = 0; i <= N-1; i++)
		{
			Mn = Mn + double(Mat(i, N - i)) * double(N + 1);
		}
	}
	Mn /= double(total_polymers);

	for (int N = 0; N <= max_degree_of_polymerization - 1; N++)
	{
		for (int i = 0; i <=N-1; i++)
		{
			Mw = Mw + double(Mat(i, N - i)) * double(N + 1) * double(N + 1); 
		}
	}
	for (int N = 0; N <= max_degree_of_polymerization - 1; N++)
	{
		for (int i = 0; i <= N - 1; i++)
		{
			Nw = Nw + double(Mat(i, N - i)) * double(N + 1);
		}
	}
	Mw = Mw / Nw;
	PDI = Mw / Mn;

	kinetics << "***** Molecular Weight Distribution Report *****" << endl;
	kinetics << "Nn = " << Mn << "\tNw = " << Mw << "\tPDI = " << PDI << endl;
	kinetics << "gridView = {";
	for (int i = 0; i <= Mat.rows() - 1; i++)
	{
		kinetics << "{";
		for (int j = 0; j <= Mat.cols() - 1; j++)
		{
			kinetics << Mat(i, j) << ",";
		}
		kinetics << Mat(i, Mat.cols() - 1) << "},";
	}
	kinetics << "};" << endl;

	for (int N = 0; N <= max_degree_of_polymerization - 1; N++)
	{
		for (int i = 0; i <= total_polymers - 1; i++)
		{
			if (polymers[i].getDegreeOfPolymerization() == (N + 1))
			{
				Lin(N) += 1.;
			}
		}
	}
	Lin /= double(total_polymers);

	kinetics << "distribution={";
	for (int i = 0; i <= Lin.size() - 2; i++)
	{
		kinetics << "{" << i + 1 << "," << Lin(i) << "},";
	}
	kinetics << "{" << Lin.size() << "," << Lin(Lin.size() - 1) << "}};";
	kinetics << endl;

	kinetics.close();
	
	return;
}
/****************************************************************************************/
// composition distribution
void compositionDistribution_report(vector<polymer> polymers,string composition_file)
{
	//next version should bin and histogram.
	ofstream kinetics;
	vector<double>	composition(polymers.size());
	// for each chain
	for (int i=0;i<=polymers.size()-1;i++)
	{
		composition[i] = polymers[i].compositionOf(1);
	}
	kinetics.open(composition_file.c_str(),ios::out);
	kinetics.precision(14);
	kinetics << "composition = {";
	for (int i=0;i<=polymers.size()-2;i++)
	{
		kinetics << "{" << i+1 << "," << composition[i] << "},";
	}
	kinetics << "{" << polymers.size() << "," << composition[polymers.size()-1] << "}};" << endl;
	kinetics << "compositionHistogram = {";
	for (int i=0;i<=polymers.size()-2;i++)
	{
		kinetics << "{" << composition[i] << "},";
	}
	kinetics << "{" << composition[polymers.size()-1] << "}};" << endl;
	kinetics.close();
	return;
}
/****************************************************************************************/
// output matrix of polymer chains at end of polymerization, in easy to read format

void	matrix_report(vector<polymer> polymers, string matrix, int chains, double volume)
{
	// this function literally just outputs the matrix describing all the chains
	ofstream kinetics;

	kinetics.open(matrix.c_str(), ios::out);
	kinetics.precision(14);

	//kinetics<<"volume = "<< volume << endl;
	//kinetics<<"chains = "<<chains << endl;
	//kinetics<<"matrix = {" <<endl;
	for (int i = 0; i <= polymers.size()-1; i++)
	{
		for (int j = 0; j <= polymers[i].getDegreeOfPolymerization() - 1; j++)
		{
			kinetics << polymers[i].sequence[j] << ",";
		}
		kinetics << endl;
	}

	kinetics << "end" << endl;
	kinetics.close();
	return;
}

/****************************************************************************************/
// output for chain test

void chaintest_report(vector<polymer> polymers, string chaintest, int A, int nA0, int B, int nB0)
{
	//Just output conversion at time t=10 minutes for testing effect of number of chains
	ofstream kinetics;

	kinetics.open(chaintest.c_str(), ios::out);
	kinetics.precision(14);

	kinetics << "conversion at t=10 min = " << conversion(A, nA0, B, nB0) << endl;

	kinetics.close();
	return;
}
/****************************************************************************************/
double	conversion(int A,int nA0,int B,int nB0)
{
	return ((double(nA0-A)+double(nB0-B))/double(nA0+nB0));
}
/****************************************************************************************/
int totalDyads(vector<polymer> polymers)
{
	int total_dyads=0;
	for (int i = 0; i <= polymers.size()-1; i++)
	{
		if (polymers[i].getDegreeOfPolymerization() >= 2)
		{
			total_dyads+=(polymers[i].getDegreeOfPolymerization()-1);
		}
	}
	return total_dyads;
}
/****************************************************************************************/
int totalTriads(vector<polymer> polymers)
{
	int total_triads = 0;
	for (int i = 0; i <= polymers.size() - 1; i++)
	{
		if (polymers[i].getDegreeOfPolymerization() >= 3)
		{
			total_triads += (polymers[i].getDegreeOfPolymerization() - 2);
		}
	}
	return total_triads;
}
/****************************************************************************************/
int totalTetrads(vector<polymer> polymers)
{
	int total_tetrads = 0;
	for (int i = 0; i <= polymers.size() - 1; i++)
	{
		if (polymers[i].getDegreeOfPolymerization() >= 4)
		{
			total_tetrads += (polymers[i].getDegreeOfPolymerization() - 3);
		}
	}
	return total_tetrads;
}
/****************************************************************************************/
int totalPentads(vector<polymer> polymers)
{
	int total_pentads = 0;
	for (int i = 0; i <= polymers.size() - 1; i++)
	{
		if (polymers[i].getDegreeOfPolymerization() >= 3)
		{
			total_pentads += (polymers[i].getDegreeOfPolymerization() - 4);
		}
	}
	return total_pentads;
}

/*****************************************************************************************/
void maldi_model(vector<polymer> polymers, string maldi)
{
	int total_polymers= polymers.size();
	int ones, twos;
	int max_degree_of_polymerization = 0;
	int nlines = 0;
	double Mn = 0;
	double Mw = 0;
	ofstream kinetics;

	for (int i = 0; i <= total_polymers - 1; i++)
	{
		if (max_degree_of_polymerization < polymers[i].getDegreeOfPolymerization())
		{
			max_degree_of_polymerization = polymers[i].getDegreeOfPolymerization();
		}
	}
	MatrixXd Mat = MatrixXd::Zero(max_degree_of_polymerization + 1, max_degree_of_polymerization + 1);
	for (int i = 0; i <= total_polymers - 1; i++)
	{
		ones = 0;	twos = 0;
		for (int j = 0; j <= polymers[i].getDegreeOfPolymerization() - 1; j++)
		{
			// count all polymers with length j
			if (polymers[i].repeatUnitAtPosition(j) == 1)
			{
				ones += 1;
			}
			else if (polymers[i].repeatUnitAtPosition(j) == 2)
			{
				twos += 1;
			}
		}
		// count chain in right place based on ones and twos
		Mat(ones, twos) += 1.; // chain
	}
	kinetics.open(maldi.c_str(), ios::out);

	//Probably not the most efficient to do this whole loop twice but want to have the number of lines before everything else in the file
	for (int i = 0; i <= Mat.rows() - 1; i++)
	{
		for (int j = 0; j <= Mat.cols() - 1; j++)
		{
			if (Mat(i, j) > 0)
			{
				nlines += 1;
			}
		}
	}
	kinetics << "number lines= " << nlines << endl;
	kinetics << "***** MALDI model data, linear chains *****" << endl;
	kinetics << "G" << "\t" << "L" << "\t" << "population" << endl;
	for (int i = 0; i <= Mat.rows() - 1; i++)
	{
		for (int j = 0; j <= Mat.cols() - 1; j++)
		{
			if (Mat(i, j) > 0)
			{
				kinetics << i << "\t" << j << "\t" << Mat(i, j) << endl;
			}
		}
	}

	kinetics << "end" << endl;

	/*kinetics << "*****Even L and G chains*******" << endl;
	kinetics << "mass" << "\t" << "number" << endl;
	for (int i = 0; i <= Mat.rows() - 1; i++)
	{
		for (int j = 0; j <= Mat.cols() - 1; j++)
		{
			if (i % 2 == 0 && j % 2 == 0)
			{
				kinetics << (i * 58.005480) + (j * 72.021130) + 209.188135 << "\t" << Mat(i, j) << endl;
				//cout << "I%2" << i%2 << "j%2" << j%2 << endl;
			}
		}

	}
	kinetics << "*****Odd L, Even G chains*******" << endl;
	kinetics << "mass" << "\t" << "number" << endl;
	for (int i = 0; i <= Mat.rows() - 1; i++)
	{
		for (int j = 0; j <= Mat.cols() - 1; j++)
		{
			if (i % 2 == 0 && j % 2 == 1)
			{
				kinetics << (i * 58.005480) + (j * 72.021130) + 209.188135 << "\t" << Mat(i, j) << endl;
			}
		}

	}
	kinetics << "*****Even L, Odd G chains*******" << endl;
	kinetics << "mass" << "\t" << "number" << endl;
	for (int i = 0; i <= Mat.rows() - 1; i++)
	{
		for (int j = 0; j <= Mat.cols() - 1; j++)
		{
			if (i % 2 == 1 && j % 2 == 0)
			{
				kinetics << (i * 58.005480) + (j * 72.021130) + 209.188135 << "\t" << Mat(i, j) << endl;
			}
		}

	}
	kinetics << "*****Odd L and G chains*******" << endl;
	kinetics << "mass" << "\t" << "number" << endl;
	for (int i = 0; i <= Mat.rows() - 1; i++)
	{
		for (int j = 0; j <= Mat.cols() - 1; j++)
		{
			if (i % 2 == 1 && j % 2 == 1)
			{
				kinetics << (i * 58.005480) + (j * 72.021130) + 209.188135 << "\t" << Mat(i, j) << endl;
			}
		}

	}*/
	kinetics.close();
	return;
}
/******************************************************************************************************/
void maldi_cycles(vector<polymer> cycles, string cycle) //Separate file for macrocycles for easily putting into isospec
{
	int total_polymers = cycles.size();
	int ones, twos;
	int max_degree_of_polymerization = 0;
	int nlines = 0;
	ofstream kinetics;

	for (int i = 0; i <= total_polymers - 1; i++)
	{
		if (max_degree_of_polymerization < cycles[i].getDegreeOfPolymerization())
		{
			max_degree_of_polymerization = cycles[i].getDegreeOfPolymerization();
		}
	}
	MatrixXd Mat = MatrixXd::Zero(max_degree_of_polymerization + 1, max_degree_of_polymerization + 1);
	for (int i = 0; i <= total_polymers - 1; i++)
	{
		ones = 0;	twos = 0;
		for (int j = 0; j <= cycles[i].getDegreeOfPolymerization() - 1; j++)
		{
			// count all polymers with length j
			if (cycles[i].repeatUnitAtPosition(j) == 1)
			{
				ones += 1;
			}
			else if (cycles[i].repeatUnitAtPosition(j) == 2)
			{
				twos += 1;
			}
		}
		// count chain in right place based on ones and twos
		if (ones > 0 || twos > 0)
		{
			Mat(ones, twos) += 1.;
		}
	}
	kinetics.open(cycle.c_str(), ios::out);

	//Probably not the most efficient to do this whole loop twice but want to have the number of lines before everything else in the file
	for (int i = 0; i <= Mat.rows() - 1; i++)
	{
		for (int j = 0; j <= Mat.cols() - 1; j++)
		{
			if (Mat(i, j) > 0)
			{
				nlines += 1;
			}
		}
	}
	kinetics << "number lines= " << nlines << endl;
	kinetics << "***** MALDI model data, cyclic chains *****" << endl;
	kinetics << "G" << "\t" << "L" << "\t" << "population" << endl;
	for (int i = 0; i <= Mat.rows() - 1; i++)
	{
		for (int j = 0; j <= Mat.cols() - 1; j++)
		{
			if (Mat(i, j) > 0)
			{
				kinetics << i << "\t" << j << "\t" << Mat(i, j) << endl;
			}
		}
	}

	kinetics << "end" << endl;

	kinetics.close();
	return;
}