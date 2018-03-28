#include<stdio.h>
#include<stdlib.h>
#include<string>
#include<string.h>
#include<math.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<map>
#include<algorithm>

using namespace std;

bool ANG_FLAG = false;
bool DIH_FLAG = false;

typedef struct vec {
    double x;
    double y;
    double z;

    vec(){ x=0; y=0; z=0; }

    vec(double * arr){
        x = arr[0];
        y = arr[1];
        z = arr[2];
    }
    double * array(int d, int r = 1){
        double * arr = (double*) malloc(d * sizeof(double));
        if (d>0) arr[0] = x;
        if (d>1) arr[1] = y;
        if (d>2) arr[2] = z;
        if (d>3)
            for (int i=3; i<d; i++)
                arr[i] = r;
        return arr;
    }
    vec& operator-(const vec& v){
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }
    vec& operator+(const vec& v){
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }
    vec& operator*(const vec& v){
        x *= v.x;
        y *= v.y;
        z *= v.z;
        return *this;
    }
} vec;

// Atom structure
typedef struct atom {
    string name;    //H11, C3, etc
    string res;     //Residue name
    int n;          //atom number
    double x;       
    double y;
    double z;
    vec vector(){
        vec v;
        v.x = x;
        v.y = y;
        v.z = z;
        return v;
    }
} atom;

typedef struct angle {
    atom * a;
    angle() {
        a = (atom*) malloc(3 * sizeof(atom));
    };
    angle(vector<atom*> b){
        angle();
        a[0] = *(b[0]);
        a[1] = *(b[1]);
        a[2] = *(b[2]);
    };
} angle;

typedef struct dihedral {
    atom * a;
    dihedral() {
        a = (atom*) malloc(4 * sizeof(atom));
    };
    dihedral(vector<atom*> b){
        dihedral();
        a[0] = *(b[0]);
        a[1] = *(b[1]);
        a[2] = *(b[2]);
        a[3] = *(b[3]);
    };
} dihedral;

// Statistic data structure
typedef struct stat {
    int n;          //atom identifier ?
    double d_sum;   //Overal sum   
    double d_avg;   //Average over all data
    double sd;      //Standard Deviation
    string res;
    vector<double> data;
    stat() {    //constructor
        n=0;
        d_sum=0;
        d_avg=0;
        sd=0;
        data = vector<double>();
    }
} stat;

// Trim word tail
static inline std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
        std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

// Shortest distance between two points on 3-dimensional space
double dist(atom A, atom B){
    double d = sqrt((A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y)+(A.z-B.z)*(A.z-B.z));
    return d;
}

double modulo(vec A){
    double m = sqrt(A.x*A.x+A.y*A.y+A.z*A.z);
    return m;
}

double dotproduct(vec A, vec B){
    double dp = A.x*B.x + A.y*B.y + A.z*B.z;
    return dp;
}

double calculate_angle(angle ang){
    atom A = ang.a[0];
    atom B = ang.a[1];
    atom C = ang.a[2];
    vec vec1 = A.vector()-B.vector();
    vec vec2 = C.vector()-B.vector();
    
    double r = acos ((dotproduct(vec1,vec2))/(modulo(vec1)*modulo(vec2)));
    return r;
}

double calculate_angle(vector< atom * > ang){
    atom A = *ang[0];
    atom B = *ang[1];
    atom C = *ang[2];
    vec vec1 = A.vector()-B.vector();
    vec vec2 = C.vector()-B.vector();
    
    double r = acos ((dotproduct(vec1,vec2))/(modulo(vec1)*modulo(vec2)));
    return r;
}

double calculate_angle(vec atomsvec[2]){
    vec vec1 = atomsvec[0];
    vec vec2 = atomsvec[1];
    
    double r = acos ((dotproduct(vec1,vec2))/(modulo(vec1)*modulo(vec2)));
    return r;
}

double calculate_angle(vec vec1, vec vec2){
    double r = acos ((dotproduct(vec1,vec2))/(modulo(vec1)*modulo(vec2)));
    return r;
}

double * matprodarray(double ** M, double * v1){
    int m, n, mn; // m linhas, n colunas
    //mn = (int)sizeof(M)/(int)sizeof(double*);
    //n = (int)sizeof(M[0])/sizeof(double);
    //m = m/n;

    //int vsize = (int)sizeof(v1)/(int)sizeof(double);

    //if (vsize != n){
    //    printf("Produto incompativel: matriz %d x %d por vetor %d x 1", m, n, vsize);
    //    return NULL;
    //}
 
    double * v2 = (double*) malloc(3 * sizeof(double));
    v2[0]=0; v2[1]=0; v2[2]=0;

    for (int i=0; i<2; i++)
        for (int j=0; j<4; j++)
            v2[i]+=M[i][j]*v1[j];

    return v2;
}

double calculate_dihedral(vector <atom*> dih){
    atom A = *dih[0];
    atom B = *dih[1]; //P1
    atom C = *dih[2]; //P2
    atom D = *dih[3];
    
    double x1, x2, y1, y2, z1, z2;
    x1 = B.x;
    x2 = C.x;
    y1 = B.y;
    y2 = C.y;
    z1 = B.z;
    z2 = C.z;

    double Dy, Dx, Dz, D1, D2;
    Dx = x2-x1;
    Dy = y2-y1;
    Dz = z2-z1;
    
    D1 = sqrt(Dx*Dx + Dz*Dz);
    D2 = sqrt(Dx*Dx + Dy*Dy + Dz*Dz);

    double ** M;
    double m03, m13, m23;
    M = (double**) malloc(2*sizeof(double*));
    M[0] = (double*) malloc(4*sizeof(double));
    M[1] = (double*) malloc(4*sizeof(double));

    m03 = (-x1*Dz + z1*Dx)/D1;
    m13 = Dy*(-x1*Dx - z1*Dz)/D1/D2 + y1*D1/D2;
    m23 = (x1*Dx + y1*Dy + z1*Dz)/D2;

    M[0][0] = Dz/D1;       M[0][1] = 0;      M[0][2] = -Dx/D1;      M[0][3] = m03;
    M[1][0] = Dx*Dy/D1/D2; M[1][1] = -D1/D2; M[1][2] = Dy*Dz/D1/D2; M[1][3] = m13;
    //M[2][0] = -Dx/D1;      M[2][1] = -Dy/D2; M[2][2] = -Dz/D2;      M[2][3] = m23;
    //M[3][0] = 0;           M[3][1] = 0;      M[3][2] = 0;           M[3][3] = 1;
    
    vec vec1 = vec(matprodarray(M, (A.vector()).array(4)));
    vec vec2 = vec(matprodarray(M, (D.vector()).array(4)));
    return calculate_angle(vec1, vec2);
}

//Print basic help with usage samples
void print_help(){
    printf("USAGE\n    hydrogen-dist file.pdb [options]\n\n");
    printf("EXAMPLE\n    hydrogen-dist file.pdb             Calcula distancias de todos os frames \n");
    printf("    hydrogen-dist file.pdb -b 1 -e 10  Apenas do frame 1 ao 10 \n");
    printf("    hydrogen-dist file.pdb -d 5        Salva apenas distancias menores que 5 \n");
    printf("    hydrogen-dist file.pdb --all       Salva apenas um arquivo com as médias sobre todos\n");
    printf("                                       os frames com os respectivos desvios padrão.\n");
}

int isnumber(string str) {
    for (int i=0; i<str.size(); i++)
        if (str[i] < 48 || str[i] > 57)
            return 0;
    return 1;
}

int main(int argc, char *argv[]){
    bool first_frame=true, all=false;
    int reference_count=0, frames=0;
    double t, max_dist=-1, begin=-1, end=-1;
    char * input_name, * fname_groups;
    string fname_prefix, line_string="", gline_string, reference_name;
    ifstream input, groupsf; 
    vector< int > other_count;  //Count each residues atoms
    vector< atom > reference;   //Atoms of reference molecule
    vector< vector< atom > > atoms;     
    vector< vector< vector< vector< double > > > > distances;
    vector< angle > AngleAtomsID;
    vector< vector < atom* > > AngleAtoms;
    vector< dihedral > DihedralAtomsID;
    vector< vector < atom* > > DihedralAtoms;
    map< string, int > rn;
    map< string, map< string, int > > types;
    map< string, map< string, string > > groups;
    map< string, map< string, map < string, stat > > > groups_stats;

    if (argc<2){
        print_help();
        return 1;
    }
    
    //-----------------------------------------------------------------------------
    // Arguments read
    //-----------------------------------------------------------------------------
    for (int i=1; i<argc; i++){
        if (argv[i][0]=='-') {
            // -d 5.0 or --max-dist 5.0 -> maximum distance to consider
            if (!strcmp(argv[i],"--max-dist") || !strcmp(argv[i], "-d")) {
                if (argc>i+1){
                    stringstream ss(string(argv[++i]));
                    ss >> max_dist;
                }
                else {
                    printf("Expected a distance after '%s'.\nExiting.\n", argv[i]);
                    return 1;
                }
            } else
            if (!strcmp(argv[i],"--all")) {
                all=true;
            } else
            // -b 2 -> first frame
            if (!strcmp(argv[i],"-b")) {
                if (argc>i+1){
                    stringstream ss(string(argv[++i]));
                    ss >> begin;
                }
            } else
            // -e 100 -> last frame
            if (!strcmp(argv[i],"-e")) {
                if (argc>i+1) {
                    stringstream ss(string(argv[++i]));
                    ss >> end;
                }
            } else
            // -o output -> output prefix
            if (!strcmp(argv[i],"-o")) {
                fname_prefix = string(argv[++i]);
            } else
            // -g groupsfile -> groups definition input file
            if (!strcmp(argv[i],"-g")) {
                fname_groups = argv[++i];
                groupsf.open(fname_groups, ifstream::in);
                if (!groupsf){
                    cout<<"Couldn't open groups file "<<string(argv[i])<<endl;
                    return 1;
                }
                printf("Using groups file: %s\n", fname_groups);
            } else
            // -r resname or --reference resname -> reference molecule name
            if (!strcmp(argv[i],"--reference") || !strcmp(argv[i], "-r")) {
                reference_name = argv[++i];
            } else {
                printf("Unrecognized option: %s\nExiting.\n", argv[i]);
                return 1;
            }
        } else {
            // Input PDB file
            input_name = argv[i];
            input.open(input_name, ifstream::in);
            if (!input){
                cout<<"Couldn't open file "<<string(argv[i])<<endl;
                return 1;
            }
            else printf("Input file: %s\n", input_name);
        }
    }

    if (!input.is_open()){
        printf("No input file specified!\n");
        print_help();
        return 1;
    }

    if (reference_name == ""){
        printf("Reference molecule: ");
        cin >> reference_name;
    }

    string filename_str(input_name);
    if (fname_prefix == "")
        fname_prefix = string(filename_str.begin(), filename_str.end()-4);

    if (max_dist>0) printf("Output distances < %.3lf ", max_dist);
    else printf("Output all distances ");

    if (begin>0 || end>=0){
        if (end>=begin)
            printf("from frames %.1f to %.1f\n", begin, end);
        else if (end>0) {
            printf("- ERROR\nBegining frame (-b) must be greater then or equal end \
frame (-e).\n");
            printf("Exiting.\n");
            return(1);
        } else 
            printf("starting from frame %.1f\n", begin);
    } else
        printf("from all frames.\n");


    //-----------------------------------------------------------------------------
    //-----------------------------------------------------------------------------


    //-----------------------------------------------------------------------------
    // Groups file read
    //-----------------------------------------------------------------------------
    FILE * output_angle, * output_dihedral;
    map < string, FILE * > output_stat;
    map < string, FILE * > output_grp_stat;
    stringstream gline_stream;
    if (groupsf){
        printf("Reading groups file...\n");
        string word, rname, gname, aname, section, token;
	int gnumber=0;
        while(getline(groupsf,gline_string)){
            //clearing the stream
            gline_stream.str(string());
            gline_stream.clear();
            gline_stream << gline_string;
            if (gline_string[0] == '[') {
                gline_stream >> word;     //garbage ('[')
                gline_stream >> section;    //residue name
                
                if (section == "angle"){
                    ANG_FLAG = true;
                    string output_angle_name = fname_prefix+"_"+section+"s.dat";
                    output_angle=fopen(output_angle_name.c_str(), "w");
                    printf("Opening %s for writing angles values\n", output_angle_name.c_str());

                    string output_stat_angle_name = fname_prefix+"_"+section+"s_AVG.dat";
                    output_stat[section]=fopen(output_stat_angle_name.c_str(), "w");
                    printf("Opening %s for writing angles statistics\n", output_stat_angle_name.c_str());
                }
                else if (section == "dihedral"){
                    DIH_FLAG = true;
                    string output_dihedral_name = fname_prefix+"_"+section+"s.dat";
                    output_dihedral=fopen(output_dihedral_name.c_str(), "w");
                    printf("Opening %s for writing dihedrals values\n", output_dihedral_name.c_str());

                    string output_stat_dihedral_name = fname_prefix+"_"+section+"s_AVG.dat";
                    output_stat[section]=fopen(output_stat_dihedral_name.c_str(), "w");
                    printf("Opening %s for writing dihedrals statistics\n", output_stat_dihedral_name.c_str());
                }
                else if (section == "energy"){
                    string output_stat_energy_name = fname_prefix+"_"+section+"AVG.dat";
                    output_stat[section]=fopen(output_stat_energy_name.c_str(), "w");
                    printf("Opening %s for writing energy statistics\n", output_stat_energy_name.c_str());
                }
                else {
                    rname = section;
                    section = "group";
                    string output_stat_grp_name = fname_prefix+"_H-dist_"+rname+"_GROUP.dat";
                    output_grp_stat[rname]=fopen(output_stat_grp_name.c_str(), "w");
                    printf("Opening %s for writing groups statistics\n", output_stat_grp_name.c_str());

                    groups[rname] = map< string, string >();

                    types[rname] = map<string,int>();
                    gnumber=0;
                }
                gline_stream >> word;     //garbage
                //clearing the stream
                gline_stream.str(string());
                gline_stream.clear();
                getline(groupsf,gline_string);
                gline_stream << gline_string;
            }
            if (section == "angle"){
                struct angle Angle;
                for (int a=0; a<3; a++){
                    gline_stream >> token; //atom number or residue name
                    if (!isnumber(token)) { //token é residuo
                        Angle.a[a].res=rtrim(token);
                        gline_stream >> token;
                        Angle.a[a].name=rtrim(token);
                        Angle.a[a].n=-1;
                    } else 
                        Angle.a[a].n = atoi(rtrim(token).c_str());
                }
                AngleAtomsID.push_back(Angle);
            } 
            else if (section == "dihedral"){
                struct dihedral Dihedral;
                for (int a=0; a<4; a++){
                    gline_stream >> token; //atom number or residue name
                    if (!isnumber(token)) { //token é residuo
                        Dihedral.a[a].res=rtrim(token);
                        gline_stream >> token;
                        Dihedral.a[a].name=rtrim(token);
                        Dihedral.a[a].n=-1;
                    } else 
                        Dihedral.a[a].n = atoi(rtrim(token).c_str());
                }
                DihedralAtomsID.push_back(Dihedral);
            }
            else if (section == "energy"){

            }
            else if (section == "group"){
                gline_stream >> gname;    //group label

                types[rname][gname]=++gnumber;

                while (gline_stream >> aname) {//atom name
                    groups[rname][aname] = gname;
                }
            }
            else cout << "Something is wrong: section isn't defined properly\n";
        }
    }
    vector< vector<double> > angles(AngleAtomsID.size());
    vector< vector<double> > dihedrals(DihedralAtomsID.size());
    //-----------------------------------------------------------------------------
    //-----------------------------------------------------------------------------

    vector < FILE * > output, output_dist_stat;

    bool skip=false, finish=false;

    //----------------------------------------------------------------------------
    // File read section - Specific for PDB format, with strict text formatation.
    //----------------------------------------------------------------------------
    line_string="";
    while(getline(input,line_string)){


        stringstream line_stream(line_string);
        string var, word;
        line_stream >> var;
            
        if(var!="TITLE") continue;  //Searches line contaning "t = x"
        
        //Decides whether to consider the frame or not, or terminate
        //the read
        while(line_stream >> word){
            if (word=="t="){
                line_stream >> t;
                if (end>=0 && t>end) finish=true;
                if (t<begin) skip=true;
                else skip=false;
                break;
            }
        }
        if (skip) continue;
        if (finish) break;
        
        frames++;

        while(var!="ATOM") {
            //Clear stream (fixed an undefined behaviour)---
            line_string="";
            getline(input,line_string,'\n');
            line_stream.str("");
            line_stream.clear();
            //----------------------------------------------
            line_stream << line_string;
            line_stream >> var;
        }
        
        reference.clear();
        reference_count=0;
        for(int j=0; j<atoms.size(); j++){
            atoms[j].clear();
        }

        if (first_frame) printf("Output files:\n");

        while(var=="ATOM"){
            int n;
            string number(line_string.begin()+7, line_string.begin()+11);
            stringstream integer(number);
            integer >> n;
            integer.str(string());
            integer.clear();

            string name(line_string.begin()+13, line_string.begin()+17);
            string res(line_string.begin()+17, line_string.begin()+21);

            double x, y, z;
            string real_str(line_string.begin()+32, line_string.begin()+39);
            stringstream real(real_str);
            real >> x;
            real.str(string());
            real.clear();

            real_str = string(line_string.begin()+40, line_string.begin()+47);
            real << real_str;
            real >> y;
            real.str(string());
            real.clear();

            real_str = string(line_string.begin()+48, line_string.begin()+55);
            real << real_str;
            real >> z;
            real.str(string());
            real.clear();
            
            struct atom atm;
            atm.n = n;
            atm.name = rtrim(name);
            atm.res = rtrim(res);
            atm.x = x;
            atm.y = y;
            atm.z = z;

            if (res==reference_name) {
                reference.push_back(atm);
                if (first_frame) reference_count++;

            } else {
                if (rn.find(res)==rn.end()){    //se foi encontrado novo residuo
                    distances.resize(distances.size()+1);

                    rn[res]=atoms.size();
                    if (first_frame) other_count.resize(rn[res]+1,0);
                    atoms.resize(rn[res]+1);
                    if(all){
                        string output_name = fname_prefix+"_dist_"+res+".dat";
                        output.push_back(fopen(output_name.c_str(), "w"));
                        
                        printf("\t%s\n", output_name.c_str());

                        fprintf(output[rn[res]], "Time   %s   atom  %s  atom     Dist\n", \
                            reference_name.c_str(), res.c_str());
                    }
                    string output_stat_name=fname_prefix+"_dist_"+res+"_AVG.dat";
                    output_dist_stat.push_back(fopen(output_stat_name.c_str(), "w"));

                    printf("\t%s\n", output_stat_name.c_str());

                    fprintf(output_dist_stat[rn[res]], \
                        " Id1   Label1  Id2   Label2     Dist     Sigma  Frequency\n");

                }
                atoms[rn[res]].push_back(atm);
                if (first_frame) other_count[rn[res]]++;
            }
            line_string="";
            getline(input,line_string);
            line_stream.str("");
            line_stream.clear();
            line_stream << line_string;
            line_stream >> var;
        }
        if (first_frame){
            printf("%d atoms on %s\n", reference_count, reference_name.c_str());
            for(map<string, int>::iterator it = rn.begin(); it != rn.end(); it++){
                printf("%d atoms on %s\n", other_count[it->second], it->first.c_str());
                distances[it->second].resize(reference_count);
                for(int i=0; i<reference_count; i++){
                    distances[it->second][i].resize(other_count[it->second]);
                }
            }
            // Angles
            if ( ANG_FLAG ) {
                cout << "Checking for atoms defined in angles section... ";
                for(int ang=0; ang<AngleAtomsID.size(); ang++){
                    AngleAtoms.push_back( vector <atom*>(3) );
                    AngleAtoms[ang][0] = NULL;
                    AngleAtoms[ang][1] = NULL;
                    AngleAtoms[ang][2] = NULL;
                    for(int at=0; at<3; at++)
                        // atom defined by its number
                        if(AngleAtomsID[ang].a[at].n != -1){
                            // atom belongs to reference residue
                            for(int i=0; i<reference.size(); i++) //atomos do residuo de referencia
                                if(reference[i].n == AngleAtomsID[ang].a[at].n){
                                    AngleAtoms[ang][at] = &reference[i];
                                    break;
                                }
                            if (AngleAtoms[ang][at] != NULL) continue;
                            // atom belongs to reference residue
                            for(int i=0; i<rn.size(); i++){ //residuos
                                for(int j=0; j<atoms[i].size(); j++) //atomos do residuo
                                    if(atoms[i][j].n == AngleAtomsID[ang].a[at].n){
                                        AngleAtoms[ang][at] = &atoms[i][j];
                                        break;
                                    }
                                if (AngleAtoms[ang][at] != NULL) break;
                            }
                            if (AngleAtoms[ang][at] == NULL){
                                printf("There is no atom of number n. %d on the system!\n",\
                                    AngleAtomsID[ang].a[at].n);
                                    return(1);
                            }
                        // atom defined by its residue and name
                        } else {
                            // atom belongs to reference residue
                            if (AngleAtomsID[ang].a[at].res == reference_name) {
                                for(int i=0; i<reference.size(); i++){ //atomos do residuo de referencia
                                    if(reference[i].name == AngleAtomsID[ang].a[at].name){
                                        AngleAtoms[ang][at] = &reference[i];
                                        break;
                                    }
                                }
                                if (AngleAtoms[ang][at] == NULL){
                                    printf("Couldn't find atom %s in reference residue!\n",\
                                        AngleAtomsID[ang].a[at].name.c_str());
                                    return(1);
                                }
                            }
                            // atom belongs to other residue
                            else {
                                for(map<string, int>::iterator it = rn.begin(); it != rn.end(); it++){
                                    if(it->first == AngleAtomsID[ang].a[at].res){
                                        int i=it->second;
                                        for(int j=0; j<atoms[i].size(); j++) { //atomos do residuo
                                            //cout << "comparing names '" << atoms[i][j].name << "' and '" << \
                                                AngleAtomsID[ang].a[at].name << "'" << endl;
                                            if(atoms[i][j].name == AngleAtomsID[ang].a[at].name){
                                                AngleAtoms[ang][at] = &atoms[i][j];
                                                break;
                                            }
                                        }
                                        if (AngleAtoms[ang][at] != NULL) break;
                                        printf("Couldn't find atom %s in residue %s!\n",\
                                            AngleAtoms[ang][at]->name.c_str(), it->first.c_str());
                                        return(1);
                                    }
                                }
                                if (AngleAtoms[ang][at] == NULL){
                                    printf("Couldn't find atom %s in any residues!\n",\
                                        AngleAtoms[ang][at]->name.c_str());
                                    return(1);
                                }
                            }
                        }
                }
                fprintf(output_angle, "T         N.1   Res1 Name1  N.2   Res2 Name2  N.3   Res3 Name3   Angle\n");
                cout << " ok." << endl;
            }
            // Dihedrals
            if ( DIH_FLAG ) {
                cout << "Checking for atoms defined in dihedrals section... ";
                for(int dih=0; dih<DihedralAtomsID.size(); dih++){
                    DihedralAtoms.push_back( vector <atom*>(4) );
                    DihedralAtoms[dih][0] = NULL;
                    DihedralAtoms[dih][1] = NULL;
                    DihedralAtoms[dih][2] = NULL;
                    DihedralAtoms[dih][3] = NULL;
                    for(int at=0; at<4; at++)
                        // atom defined by its number
                        if(DihedralAtomsID[dih].a[at].n != -1){
                            // atom belongs to reference residue
                            for(int i=0; i<reference.size(); i++) //atomos do residuo de referencia
                                if(reference[i].n == DihedralAtomsID[dih].a[at].n){
                                    DihedralAtoms[dih][at] = &reference[i];
                                    break;
                                }
                            if (DihedralAtoms[dih][at] != NULL) continue;
                            // atom belongs to reference residue
                            for(int i=0; i<rn.size(); i++){ //residuos
                                for(int j=0; j<atoms[i].size(); j++) //atomos do residuo
                                    if(atoms[i][j].n == DihedralAtomsID[dih].a[at].n){
                                        DihedralAtoms[dih][at] = &atoms[i][j];
                                        break;
                                    }
                                if (DihedralAtoms[dih][at] != NULL) break;
                            }
                            if (DihedralAtoms[dih][at] == NULL){
                                printf("There is no atom of number n. %d on the system!\n",\
                                    DihedralAtomsID[dih].a[at].n);
                                    return(1);
                            }
                        // atom defined by its residue and name
                        } else {
                            // atom belongs to reference residue
                            if (DihedralAtomsID[dih].a[at].res == reference_name) {
                                for(int i=0; i<reference.size(); i++){ //atomos do residuo de referencia
                                    if(reference[i].name == DihedralAtomsID[dih].a[at].name){
                                        DihedralAtoms[dih][at] = &reference[i];
                                        break;
                                    }
                                }
                                if (DihedralAtoms[dih][at] == NULL){
                                    printf("Couldn't find atom %s in reference residue!\n",\
                                        DihedralAtomsID[dih].a[at].name.c_str());
                                    return(1);
                                }
                            }
                            // atom belongs to other residue
                            else {
                                for(map<string, int>::iterator it = rn.begin(); it != rn.end(); it++){
                                    if(it->first == DihedralAtomsID[dih].a[at].res){
                                        int i=it->second;
                                        for(int j=0; j<atoms[i].size(); j++) { //atomos do residuo
                                            //cout << "comparing names '" << atoms[i][j].name << "' and '" << \
                                                DihedralAtomsID[dih].a[at].name << "'" << endl;
                                            if(atoms[i][j].name == DihedralAtomsID[dih].a[at].name){
                                                DihedralAtoms[dih][at] = &atoms[i][j];
                                                break;
                                            }
                                        }
                                        if (DihedralAtoms[dih][at] != NULL) break;
                                        printf("Couldn't find atom %s in residue %s!\n",\
                                            DihedralAtoms[dih][at]->name.c_str(), it->first.c_str());
                                        return(1);
                                    }
                                }
                                if (DihedralAtoms[dih][at] == NULL){
                                    printf("Couldn't find atom %s in any residues!\n",\
                                        DihedralAtoms[dih][at]->name.c_str());
                                    return(1);
                                }
                            }
                        }
                }
                fprintf(output_dihedral, "T         N.1   Res1 Name1  N.2   Res2 Name2  N.3   Res3 Name3   N.4   Res4 Name4   Dihedral\n");
                cout << " ok." << endl;
            }

            cout << "Primeiro frame analisado." << endl;
        }
        first_frame=false;

        // Calculo das distancias
        for(int j=0; j<rn.size(); j++){ //residuos
            for(int i=0; i<reference.size(); i++){ //atomos do residuo de referencia
                for(int k=0; k<atoms[j].size(); k++){ //atomos do residuo
                    atom A = reference[i], B = atoms[j][k];
                    double D = dist(A,B);
                    if (max_dist<0 || D<max_dist){
                        distances[j][i][k].push_back(D);
                        if (all) 
                            fprintf(output[j], "%-10.1f%4d%5d%5s%5s%8.3lf\n", \
                            t, i, k, (A.name).c_str(), (B.name).c_str(), D);
                    }
                }       
            }
        }

        // Calculo dos angulos
        for(int ang=0; ang<AngleAtoms.size(); ang++){
            
            atom A = *(AngleAtoms[ang][0]);
            atom B = *(AngleAtoms[ang][1]);
            atom C = *(AngleAtoms[ang][2]);
            double D1 = dist(A,B); 
            double D2 = dist(B,C);
            //if (max_dist<0 || (D1<max_dist && D2<max_dist)){
                double Angle = calculate_angle(AngleAtoms[ang]);
                angles[ang].push_back(Angle);
            //}
            //if (all) 
                fprintf(output_angle, "%-10.1f%-4d%6s%6s  %-4d%6s%6s  %-4d%6s%6s%8.3f\n", t, \
                AngleAtoms[ang][0]->n, AngleAtoms[ang][0]->res.c_str(), AngleAtoms[ang][0]->name.c_str(),\
                AngleAtoms[ang][1]->n, AngleAtoms[ang][1]->res.c_str(), AngleAtoms[ang][1]->name.c_str(),\
                AngleAtoms[ang][2]->n, AngleAtoms[ang][2]->res.c_str(), AngleAtoms[ang][2]->name.c_str(),\
                Angle);
        }
        // Calculo dos diedros
        for(int dih=0; dih<DihedralAtoms.size(); dih++){
            //atom A = *(DihedralAtoms[dih][0]);
            //atom B = *(DihedralAtoms[dih][1]);
            //atom C = *(DihedralAtoms[dih][2]);
            //double D1 = dist(A,B); 
            //double D2 = dist(B,C);
            //if (max_dist<0 || (D1<max_dist && D2<max_dist)){
                double Dihedral = calculate_dihedral(DihedralAtoms[dih]);
                dihedrals[dih].push_back(Dihedral);
            //}
            //if (all) 
                fprintf(output_dihedral, "%-10.1f%-4d%6s%6s  %-4d%6s%6s  %-4d%6s%6s  %-4d%6s%6s%8.3f\n", t, \
                DihedralAtoms[dih][0]->n, DihedralAtoms[dih][0]->res.c_str(), DihedralAtoms[dih][0]->name.c_str(),\
                DihedralAtoms[dih][1]->n, DihedralAtoms[dih][1]->res.c_str(), DihedralAtoms[dih][1]->name.c_str(),\
                DihedralAtoms[dih][2]->n, DihedralAtoms[dih][2]->res.c_str(), DihedralAtoms[dih][2]->name.c_str(),\
                DihedralAtoms[dih][3]->n, DihedralAtoms[dih][3]->res.c_str(), DihedralAtoms[dih][3]->name.c_str(),\
                Dihedral);
        }
    }
    //--------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------

    //----------------------------------------------------------------------------
    // Distances and statistics calculations, and output.
    //----------------------------------------------------------------------------
    printf("%d frames analysed\n", frames);
    printf("Calculating averages and distances\n");
    string reference_grp, res_grp;
    stat data;
    int m=0;
    for(int i=0; i<distances.size();i++){
        for(int j=0; j<distances[i].size();j++){
            for(int k=0; k<distances[i][j].size();k++){
		string rname=atoms[i][k].res;
		string aname=atoms[i][k].name;
                //Group section
                bool count_this=false;
                if (groupsf.is_open()){
                    if(groups.find(rname)!=groups.end() && \
                        groups[reference_name].find(reference[j].name)!=groups[reference_name].end()){
                        if(groups[rname].find(aname)!=groups[rname].end()){
                            reference_grp = groups[reference_name][reference[j].name];
                            res_grp = groups[rname][aname];
                            if(groups_stats.find(reference_grp)==groups_stats.end()){
                                groups_stats[reference_grp] = map< string, map< string, stat > >();
                                groups_stats[reference_grp][rname] = map< string, stat >();
                                groups_stats[reference_grp][rname][res_grp] = stat();
                            }
                            else if (groups_stats[reference_grp].find(rname)==groups_stats[reference_grp].end())
                                groups_stats[reference_grp][rname] = map< string, stat >();

                            else if (groups_stats[reference_grp][rname].find(res_grp)==\
                                groups_stats[reference_grp][rname].end())
                                groups_stats[reference_grp][rname][res_grp] = stat();

                            data = groups_stats[reference_grp][rname][res_grp];
                            data.res = rname;
                            count_this=true;
                        }
                    }
                }
                double avg=0, sum=0;
                int N = distances[i][j][k].size();
                //if (N) {
                    ///////////////
                    double d, stddev=0;
                    if (count_this){
                        for(int l=0; l<N;l++){
                            d = distances[i][j][k][l];
                            sum += d;

                            data.data.push_back(d);
                        }
                        data.d_sum += sum;
                        data.n += N;
                        groups_stats[reference_grp][rname][res_grp] = data;

                    } else
                        for(int l=0; l<N;l++){
                            d = distances[i][j][k][l];
                            sum += d;
                        }

                if (N){
                    avg = sum/N;
                    for(int l=0; l<N; l++){
                        d = distances[i][j][k][l];
                        stddev += (d-avg)*(d-avg);
                    }
                    stddev = sqrt(stddev/(double)N);
                    fprintf(output_dist_stat[i], "%4d %8s %4d %8s %8.3lf %9.3lf %10d\n",\
                    j, reference[j].name.c_str(), k, atoms[i][k].name.c_str(), avg, stddev, N);

                }
            }
        }
    }
    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------

    for(map<string, FILE*>::iterator it = output_grp_stat.begin(); \
        it!=output_grp_stat.end(); it++){
        fprintf(it->second, " Id1   Label1  Id2   Label2     Dist     Sigma  Frequency\n");
    }

    int ni=0, nj=0;
    if (groupsf.is_open()){
        for(map< string, map< string, map< string, stat > > >::iterator reference_it = groups_stats.begin();\
            reference_it!=groups_stats.end(); reference_it++){
	    ni=types[reference_name][reference_it->first];
            for(map< string, map< string, stat > >::iterator res_n = (reference_it->second).begin(); \
                res_n!=(reference_it->second).end(); res_n++){
                for(map< string, stat >::iterator res_it = (res_n->second).begin(); \
                    res_it!=(res_n->second).end(); res_it++){
                    
                    int n = (res_it->second).n;
                    double d_sum = (res_it->second).d_sum;
                    double avg = (n?(double)d_sum/n:0);
                    string rname = (res_it->second).res;
                    (res_it->second).d_avg = avg;
                    vector<double> data = (res_it->second).data;
                    double d, stddev=0;
                    if (n){
                        for(int i=0; i<data.size(); i++){
                            d = data[i];
                            stddev += (d-avg)*(d-avg);
                        }
                        stddev = sqrt(stddev/(double)n);
                    }
                    nj=types[rname][res_it->first];
                    fprintf(output_grp_stat[rname], "%4d %8s %4d %8s %8.3lf %9.3lf %10d\n", \
                        ni, reference_it->first.c_str(), nj, res_it->first.c_str(), avg, stddev, n);
                }
            }
        }
    }
    cout << "Done with distances." << endl;

    //----------------------------------------------------------------------------
    // Angles and statistics calculation, and output.
    //----------------------------------------------------------------------------

    if ( ANG_FLAG ) {
        fprintf(output_stat["angle"], " N.1  Res1  Name1  N.2  Res2  Name2  N.3  Res3  Name3  Avg. angle  Sigma  Frequency\n");
        stat anglesdata;
        double avg=0;
        double stddev=0;
        for (int ang=0; ang<angles.size(); ang++){
            int n = angles[ang].size();
            for (int f=0; f<n; f++){
                avg+=angles[ang][f];
            }
            avg /= n;
            for (int f=0; f<n; f++){
                double anglevalue = angles[ang][f];
                stddev += (avg-anglevalue)*(avg-anglevalue);
            }
            stddev = sqrt(stddev/(double)n);
            fprintf(output_stat["angle"], "%4d %6s %6s %4d %6s %6s %4d %6s %6s %8.3lf %9.3lf %7d\n", \
                AngleAtoms[ang][0]->n, AngleAtoms[ang][0]->res.c_str(), AngleAtoms[ang][0]->name.c_str(),\
                AngleAtoms[ang][1]->n, AngleAtoms[ang][1]->res.c_str(), AngleAtoms[ang][1]->name.c_str(),\
                AngleAtoms[ang][2]->n, AngleAtoms[ang][2]->res.c_str(), AngleAtoms[ang][2]->name.c_str(),\
                avg, stddev, n);
        }
        cout << "Done with angles." << endl;
        fclose(output_angle);
    }

    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------

    //----------------------------------------------------------------------------
    // Dihedrals and statistics calculation, and output.
    //----------------------------------------------------------------------------

    if ( DIH_FLAG ) {
        fprintf(output_stat["dihedral"], " N.1  Res1  Name1  N.2  Res2  Name2  N.3  Res3  Name3  N.4  Res4  Name4  Avg. dihedral  Sigma  Frequency\n");
        stat dihedralsdata;
        double avg=0;
        double stddev=0;
        for (int dih=0; dih<dihedrals.size(); dih++){
            int n = dihedrals[dih].size();
            for (int f=0; f<n; f++){
                avg+=dihedrals[dih][f];
            }
            avg /= n;
            for (int f=0; f<n; f++){
                double dihedralvalue = dihedrals[dih][f];
                stddev += (avg-dihedralvalue)*(avg-dihedralvalue);
            }
            stddev = sqrt(stddev/(double)n);
            fprintf(output_stat["dihedral"], "%4d %6s %6s %4d %6s %6s %4d %6s %6s %4d %6s %6s %8.3lf %9.3lf %7d\n", \
                DihedralAtoms[dih][0]->n, DihedralAtoms[dih][0]->res.c_str(), DihedralAtoms[dih][0]->name.c_str(),\
                DihedralAtoms[dih][1]->n, DihedralAtoms[dih][1]->res.c_str(), DihedralAtoms[dih][1]->name.c_str(),\
                DihedralAtoms[dih][2]->n, DihedralAtoms[dih][2]->res.c_str(), DihedralAtoms[dih][2]->name.c_str(),\
                DihedralAtoms[dih][3]->n, DihedralAtoms[dih][3]->res.c_str(), DihedralAtoms[dih][3]->name.c_str(),\
                avg, stddev, n);
        }
        cout << "Done with dihedrals." << endl;
        fclose(output_dihedral);
    }
        
    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------


    for(int i=0; i<output_dist_stat.size(); i++){
        if (all) fclose(output[i]);
        fclose(output_dist_stat[i]);
    }
    
    if (groupsf.is_open()){ 
        for(map< string, FILE * >::iterator it=output_grp_stat.begin(); it!=output_grp_stat.end(); it++)
            fclose(it->second);
        for(map< string, FILE * >::iterator it=output_stat.begin(); it!=output_stat.end(); it++)
            fclose(it->second);
    }
    printf("Done.\n");

return 0;
}
