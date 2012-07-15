//move.cpp
//by Jeff Sung
//Last updated 2011-09-29 at 06:16
//move.cpp takes a PDB, reads in selections of atom numbers to move by  
//specific vectors, and moves those atoms, writing everything to a new pdb.
//The code below contains a simplistic (read: lazy) method of parsing PDBs.
//
//Syntax for file 'input'
//Odd Numbered lines: list of atom numbers to be moved (indexing starts from 1)
//Even Numbered lines: vector (x y z) by which to move the atoms listed in the 
//previous line.
//Will move arbitrary numbers of atom selections
//
//Notes: no error checking for existence or syntax of file 'input'
//       input and output are PDB files; no error checking for correct syntax
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

using namespace std;

struct atom{
/*
	string record;
	string serial;
	string name;
	char altLoc;
	string resName;
	char chainID;
	string resSeq;
	char iCode;
	string x, y, z;
	string occupancy;
	string tempFactor;
	string element;
	string charge;
*/
	int serial;
	int resSeq;
	double x, y, z;
	string before;
	string after;
	
	friend ostream &operator<<(ostream &os, const atom &at);
	friend bool operator==(atom &a1, atom &a2);
};

struct residue{
	vector <atom> atoms;
	double x, y, z; //movement vector
};

double str2double(const string &str);
int str2int(const string & str);
atom readatom(const string & str);
string trim(const string & str);
vector <atom> split(const string & str);


double str2double(const string &str){
	stringstream ss(str);
	double d;
	if (!(ss>>d)) d=0;
	return d;
}
int str2int(const string &str){
	stringstream ss(str);
	int i;
	if (!(ss>>i)) i=0;
	return i;
}
atom readatom(const string &str){
	atom at;
/*
	at.record=trim(str.substr(0,6));
	at.serial=trim(str.substr(6,5));
	at.name=trim(str.substr(12,4));
	at.altLoc=str[16];
	at.resName=trim(str.substr(17,3));
	at.chainID=str[21];
	at.resSeq=trim(str.substr(22,4));
	at.iCode=str[26];
	at.x=trim(str.substr(30,8));
	at.y=trim(str.substr(38,8));
	at.z=trim(str.substr(46,8));
	at.occupancy=trim(str.substr(54,6));
	at.tempFactor=trim(str.substr(60,6));
	at.element=trim(str.substr(76,2));
	at.charge=trim(str.substr(78,2));
*/
	at.before=str.substr(0,30);
	at.after=str.substr(54,26);
	at.serial=str2int(str.substr(6,5));
	at.resSeq=str2int(str.substr(22,4));
	at.x=str2double(str.substr(30,8));
	at.y=str2double(str.substr(38,8));
	at.z=str2double(str.substr(46,8));
	return at;
}

string trim(const string & strorig){
/*
	string str=strorig;
	string::size_type pos = str.find_last_not_of(' ');
	if(pos != string::npos) {
		str.erase(pos + 1);
		pos = str.find_first_not_of(' ');
		if(pos != string::npos) str.erase(0, pos);
	}
	else str.erase(str.begin(), str.end());
	return str;
*/
	return strorig;
}

ostream &operator<<(ostream &os, const atom &at){
/*
	os<<setiosflags(ios_base::right|ios_base::fixed);
	os<<setw(6)<<at.record;
	os<<setw(5)<<at.serial;
	os<<setw(4)<<at.name;
	os<<setw(1)<<at.altLoc;
	os<<setw(3)<<at.resName;	
	os<<setw(1)<<at.chainID;
	os<<setw(4)<<at.resSeq;
	os<<setw(1)<<at.iCode;
	os<<setw(8)<<at.x;
	os<<setw(8)<<at.y;
	os<<setw(8)<<at.z;
	os<<setw(6)<<at.occupancy;
	os<<setw(6)<<at.tempFactor;
	os<<setw(2)<<at.element;
	os<<setw(2)<<at.charge;
*/
	os<<at.before;
	os<<setiosflags(ios_base::right|ios_base::fixed);
	os<<setw(8)<<setprecision(3)<<at.x;
	os<<setw(8)<<setprecision(3)<<at.y;
	os<<setw(8)<<setprecision(3)<<at.z;
	os<<at.after;
	return os;
}

bool operator==(atom &a1, atom &a2){
	return a1.serial==a2.serial;
}

vector <atom> split(const string & str){
	vector <atom> ret(1);
	atom at;
	stringstream ss(str);
	string s2;
	int i=0;
	for (; ss>>s2;){
		at.serial=str2int(s2);
		ret[i++]=at;
		if (i==ret.size()) ret.resize(2*i);
	}
	ret.resize(i);
	return ret;
}

int main(int argc,char *argv[]){
	if (argc < 3){
		cout<<"Syntax: "<<argv[0]<<" pdb_in pdb_out"<<endl;
	}
	else{
		int which=0, i, j=0;
		atom curat;
		residue curesid;
		string line, s2;
		double d;

		//Read in the atoms that need to be moved and the requisite vectors
		ifstream infile("input");
		vector <residue> movelist(1);
		while (getline(infile,line)){
			movelist[which].atoms=split(line);/*
			for (j=0; j<movelist[which].atoms.size(); j++)
				cout<<"movelist["<<which<<"].atoms["<<j<<"].serial="<<movelist[which].atoms[j].serial<<endl;*/
			infile>>movelist[which].x>>movelist[which].y>>movelist[which].z;
			which++;
			infile.ignore(100,'\n');
			if (movelist.size()==which) movelist.resize(2*which);
		}	
		movelist.resize(which);
		infile.close();
		
		//Read in the given PDB, move any necessary atoms, and then write the resultant PDB
		infile.open(argv[1]);
		ofstream outfile(argv[2]);
		while (getline(infile,line)){
			if (line.find("ATOM")==string::npos&&line.find("HETATM")==string::npos){
				outfile<<line<<endl;
				continue;
			}
			curat=readatom(line);
			i=0;
			for (i=0; i<movelist.size(); i++){
				for (j=0; j<movelist[i].atoms.size(); j++){
					if (movelist[i].atoms[j]==curat){
						//Move atom
						curat.x+=movelist[i].x;
						curat.y+=movelist[i].y;
						curat.z+=movelist[i].z;
					}
				}
			}
			//Write
			outfile<<curat<<endl;	
		}
		infile.close();
		outfile.close();

	}
	return 0;
}
