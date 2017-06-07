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
#include <algorithm>

using namespace std;

// Atom structure
typedef struct atom {
    string name;    //H11, C3, etc
    string res;     //Residue name
    double x;       
    double y;
    double z;
} atom;

// Statistic data structure
typedef struct stat {
    int n;          
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

//Print basic help with usage samples
void print_help(){
    printf("USAGE\n    hydrogen-dist file.pdb [options]\n\n");
    printf("EXAMPLE\n    hydrogen-dist file.pdb             Calcula distancias de todos os frames \n");
    printf("    hydrogen-dist file.pdb -b 1 -e 10  Apenas do frame 1 ao 10 \n");
    printf("    hydrogen-dist file.pdb -d 5        Salva apenas distancias menores que 5 \n");
    printf("    hydrogen-dist file.pdb --all       Salva apenas um arquivo com as médias sobre todos\n");
    printf("                                       os frames com os respectivos desvios padrão.\n");
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
            printf("- ERROR\nBegining frame (-b) must be greater or equal then end \
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
    map < string, FILE * > output_grp_stat;
    stringstream gline_stream;
    if (groupsf){
        printf("Reading groups file...\n");
        string word, rname, gname, aname;
	int gnumber=0;
        while(getline(groupsf,gline_string)){
            //clearing the stream
            gline_stream.str(string());
            gline_stream.clear();
            gline_stream << gline_string;
            if (gline_string[0] == '[') {
                gline_stream >> word;     //garbage (']')
                gline_stream >> rname;    //residue name
                
                string output_stat_grp_name = fname_prefix+"_H-dist_"+rname+"_GROUP.dat";
                output_grp_stat[rname]=fopen(output_stat_grp_name.c_str(), "w");
                printf("Opening %s for writing groups statistics...\n", output_stat_grp_name.c_str());

                gline_stream >> word;     //garbage
                groups[rname] = map< string, string >();
                
                //clearing the stream
                gline_stream.str(string());
                gline_stream.clear();

                getline(groupsf,gline_string);
                gline_stream << gline_string;
                types[rname] = map<string,int>();
                gnumber=0;
            }
            gline_stream >> gname;    //group label

            types[rname][gname]=++gnumber;

            while (gline_stream >> aname) {//atom name
                groups[rname][aname] = gname;
            }
        }
    }
    //-----------------------------------------------------------------------------
    //-----------------------------------------------------------------------------

    vector < FILE * > output, output_stat;

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
                    output_stat.push_back(fopen(output_stat_name.c_str(), "w"));

                    printf("\t%s\n", output_stat_name.c_str());

                    fprintf(output_stat[rn[res]], \
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
        }

        for(int j=0; j<rn.size(); j++){ //residuos
            for(int i=0; i<reference.size(); i++){ //atomos da BCD
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
        first_frame=false;
    }
    //--------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------

    //----------------------------------------------------------------------------
    // Distances and statistics calculation, and output.
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
                    fprintf(output_stat[i], "%4d %8s %4d %8s %8.3lf %9.3lf %10d\n",\
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


    for(int i=0; i<output_stat.size(); i++){
        if (all) fclose(output[i]);
        fclose(output_stat[i]);
    }
    for(map< string, FILE * >::iterator it=output_grp_stat.begin(); it!=output_grp_stat.end(); it++){
        if (groupsf.is_open()) fclose(it->second);
    }
    printf("Done.\n");

return 0;
}
