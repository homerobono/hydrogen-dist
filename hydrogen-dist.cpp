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

typedef struct atom {
    string name;
    string res;
    double x;
    double y;
    double z;
} atom;

typedef struct stat {
    int n;
    double d_sum;
    double d_avg;
    double sd;
    string res;
    vector<double> data;
    stat() {
        n=0;
        d_sum=0;
        d_avg=0;
        sd=0;
        data = vector<double>();
    }
} stat;

// trim from end
static inline std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(),
                    std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
            return s;
}


double dist(atom A, atom B){
    double d = sqrt((A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y)+(A.z-B.z)*(A.z-B.z));
    return d;
}

void print_help(){
    printf("USAGE\n    ./hydrogen-dist file.pdb [options]\n\n");
    printf("EXAMPLE\n    ./hydrogen-dist file.pdb             Calcula distancias de todos os frames \n");
    printf("    ./hydrogen-dist file.pdb -b 1 -e 10  Apenas do frame 1 ao 10 \n");
    printf("    ./hydrogen-dist file.pdb -d 5        Salva apenas distancias menores que 5 \n");
    printf("    ./hydrogen-dist file.pdb --average-only   Salva apenas um arquivo com as médias sobre todos\n");
    printf("                                              os frames com os respectivos desvios padrão.\n");
}

int main(int argc, char *argv[]){
    vector< vector< atom > > atoms;
    vector< atom > bcd;
    vector< vector< vector< vector< double > > > > statistics;
    map< string, int > rn;
    map< string, map< string, int > > tipos;
    map< string, map< string, string > > groups;
    map< string, map< string, map < string, stat > > > groups_stats;
    double t, max_dist=-1, begin=-1, end=-1;
    int bcd_count=0, frames=0;
    vector<int> other_count;
    bool first_run=true, avg_only=false;
    ifstream input, groupsf; 
    char * input_name, * fname_groups;
    string fname_prefix, line_string="", gline_string;

    //for (int i=0; i<argc; i++) printf("%s\n", argv[i]);

    if (argc<2){
        print_help();
        return 1;
    }

    //Arguments read
    for (int i=1; i<argc; i++){
        if (argv[i][0]=='-') {
            // -d 5.0 ou --max-dist 5.0
            if (!strcmp(argv[i],"--max-dist") || !strcmp(argv[i], "-d")) {
                if (argc>i+1){
                    //max_dist=stof(string(argv[++i]));
                    stringstream ss(string(argv[++i]));
                    ss >> max_dist;
                }
                else {
                    printf("Expected a distance after '%s'.\nExiting.\n", argv[i]);
                    return 1;
                }
            } else
            if (!strcmp(argv[i],"--average-only")) {
                avg_only=true;
            } else
            // -b 2
            if (!strcmp(argv[i],"-b")) {
                if (argc>i+1){
                    //begin=stof(string(argv[++i]));
                    stringstream ss(string(argv[++i]));
                    ss >> begin;
                }
            } else
            // -e 100
            if (!strcmp(argv[i],"-e")) {
                if (argc>i+1) {
                    //end=stof(string(argv[++i]));
                    stringstream ss(string(argv[++i]));
                    ss >> end;
                }
            } else
            if (!strcmp(argv[i],"-o")) {
                fname_prefix = string(argv[++i]);
            } else
            if (!strcmp(argv[i],"-g")) {
                fname_groups = argv[++i];
                groupsf.open(fname_groups, ifstream::in);
                if (!groupsf){
                    cout<<"Couldn't open groups file "<<string(argv[i])<<endl;
                    return 1;
                }
                printf("Using groups file: %s\n", fname_groups);
            } else {
                printf("Unrecognized option: %s\nExiting.\n", argv[i]);
                return 1;
            }
        } else {
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
    
    string filename_str(input_name);
    if (fname_prefix=="")
        fname_prefix = string(filename_str.begin(), filename_str.end()-4);
    
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
                gline_stream >> word;     //garbage
                gline_stream >> rname;    //residue name
                
                string output_stat_grp_name = fname_prefix+"_H-dist_"+rname+"_AVG_GRP.dat";
                output_grp_stat[rname]=fopen(output_stat_grp_name.c_str(), "w");
                printf("Opening %s for writing groups statistics...\n", output_stat_grp_name.c_str());

                gline_stream >> word;     //garbage
                groups[rname] = map< string, string >();
                
                //clearing the stream
                gline_stream.str(string());
                gline_stream.clear();

                getline(groupsf,gline_string);
                gline_stream << gline_string;
                tipos[rname] = map<string,int>();
                gnumber=0;
            }
            //gline_stream >> gnumber;    //group label
            gline_stream >> gname;    //group label

            tipos[rname][gname]=++gnumber;

            while (gline_stream >> aname) {//atom name
                groups[rname][aname] = gname;
            }
        }/*
        for(map< string, map< string, int > >::iterator it = tipos.begin(); it!=tipos.end(); it++){
            cout << it->first << endl;
            for(map< string, int >::iterator it2 = it->second.begin(); it2!=it->second.end(); it2++){
                cout << it2->first << " " << it2->second << endl;
            }
        }*/
    }
    /*
    for(map< string, map< string, pair< string, int > > >::iterator it1 = groups.begin(); it1!=groups.end(); it1++){
        cout << it1->first << endl;
        for(map< string, pair< string,int > >::iterator it2 = (it1->second).begin(); it2!=(it1->second).end(); it2++){
            cout << it2->first << ": " << (it2->second).first << endl;
        }
        cout<< endl;
    }
    */

    if (max_dist>0) printf("Output distances < %.3lf ", max_dist);
    else printf("Output all distances ");

    if (begin>0 || end>=0){
        if (end>=begin)
            printf("from frames %.1f to %.1f\n", begin, end);
        else if (end>0) {
            printf("- ERROR\nBegining frame (-b) must be greater or equal then end frame (-e).\n");
            printf("Exiting.\n");
            return(1);
        } else 
            printf("starting from frame %.1f\n", begin);
    } else
        printf("from all frames.\n");

    vector < FILE * > output;
    vector < FILE * > output_stat;

    int n_of_args = 7;

    bool skip=false, finish=false;

    line_string="";
    while(getline(input,line_string)){
        stringstream line_stream(line_string);
        string var, word;
        line_stream >> var;
            
       // cout << "reading input" << endl;
       // cout << "var " << var << endl;
        
        if(var!="TITLE") continue;
        
        while(line_stream >> word){
            if (word=="t="){
                line_stream >> t;
               // cout << "reading t = " << t << endl;
                if (end>=0 && t>end) finish=true;
                if (t<begin) skip=true;
                else skip=false;
                break;
            }
        }
        if (skip) continue;
        if (finish) break;
        
        frames++;

        //cout << "reading atom" << endl;
        while(var!="ATOM") {
            line_string="";
            getline(input,line_string,'\n');
            //cout<< line_string << endl;
            line_stream.str("");
            line_stream.clear();
            line_stream << line_string;
            line_stream >> var;
            //cout << "var " << var << endl;
        }
        //cout << "readed" << endl;
        
        bcd.clear();
        bcd_count=0;
        for(int j=0; j<atoms.size(); j++){
            atoms[j].clear();
            //other_count[j]=0;
        }

        if (first_run) printf("Output files:\n");

        while(var=="ATOM"){
            string name(line_string.begin()+13, line_string.begin()+17);
            //cout << line_string << endl;
            if (name[0]=='H'){
                string res(line_string.begin()+17, line_string.begin()+21);

                //remove(name.begin(), name.end(), ' ');
                //remove(res.begin(), res.end(), ' ');
                
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

                //cout<<"res "<<res<<endl;
                if (res=="F6NH") {
                    bcd.push_back(atm);

                    if (first_run) bcd_count++;
                    //cout << "bcd atoms " << bcd_count << endl;
                } else {
                    if (rn.find(res)==rn.end()){    //se foi encontrado novo residuo
                        statistics.resize(statistics.size()+1);

                        rn[res]=atoms.size();
                        if (first_run) other_count.resize(rn[res]+1,0);
                        //cout << "rn[" << res << "] = " << rn[res] << endl;
                        //cout << "other count size " << other_count.size() <<endl;
                        //cout << "other count[i] = " << other_count.size() <<endl;
                        atoms.resize(rn[res]+1);
                        if(!avg_only){
                            string output_name = fname_prefix+"_H-dist_"+res+".dat";
                            output.push_back(fopen(output_name.c_str(), "w"));
                            
                            printf("\t%s\n", output_name.c_str());

                            fprintf(output[rn[res]], "T     F6NH atom %s atom     Distance\n", res.c_str());
                        }
                        string output_stat_name = fname_prefix+"_H-dist_"+res+"_AVG.dat";
                        output_stat.push_back(fopen(output_stat_name.c_str(), "w"));

                        printf("\t%s\n", output_stat_name.c_str());

                        fprintf(output_stat[rn[res]], " Id1   Label1  Id2   Label2     Dist     Sigma  Frequency\n");

                    }
                    atoms[rn[res]].push_back(atm);
                    if (first_run) other_count[rn[res]]++;
                    //cout << "res " << res << " other count " << other_count[rn[res]] << endl;
                }
                //printf("%.1f   ", t);
            }
            line_string="";
            getline(input,line_string);
            line_stream.str("");
            line_stream.clear();
            line_stream << line_string;
            line_stream >> var;
        }
        if (first_run){
            printf("%d hydrogens on BCD\n", bcd_count);
            for(map<string, int>::iterator it = rn.begin(); it != rn.end(); it++){
                printf("%d hydrogens on %s\n", other_count[it->second], it->first.c_str());
                statistics[it->second].resize(bcd_count);
                //cout << "itsec " << it->second <<endl;
                for(int i=0; i<bcd_count; i++){
                    statistics[it->second][i].resize(other_count[it->second]);
                    //cout << "othercount " << other_count[it->second] <<endl;
                }
            }
        }

        for(int j=0; j<rn.size(); j++){ //residuos
            for(int i=0; i<bcd.size(); i++){ //atomos da BCD
                for(int k=0; k<atoms[j].size(); k++){ //atomos do residuo
                    atom A = bcd[i], B = atoms[j][k];
                    double D = dist(A,B);
                    if (max_dist<0 || D<max_dist){
                        statistics[j][i][k].push_back(D);
                        if (!avg_only) {
                            fprintf(output[j], "%-10.1f%4d%5d%5s%5s%8.3lf\n", t, i, k, (bcd[i].name).c_str(), (atoms[j][k].name).c_str(), D);
                        }
                    }

                }       
            }
        }
        first_run=false;
    }
 

    printf("%d frames analysed\n", frames);
    printf("Calculating averages and statistics\n");
    string bcd_grp, res_grp;
    stat data;
    int m=0;
    for(int i=0; i<statistics.size();i++){
        for(int j=0; j<statistics[i].size();j++){
            for(int k=0; k<statistics[i][j].size();k++){
		string rname=atoms[i][k].res;
		string aname=atoms[i][k].name;
                //Group section
                bool count_this=false;
                if (groupsf.is_open()){
                    if(groups.find(rname)!=groups.end() && groups["F6NH"].find(bcd[j].name)!=groups["F6NH"].end()){
                        if(groups[rname].find(aname)!=groups[rname].end()){
                            //cout<< "found " << rname << " and " << aname <<endl;
                            bcd_grp = groups["F6NH"][bcd[j].name];
                            res_grp = groups[rname][aname];
                            //cout<< "resname " << rname << " resgrp " << res_grp <<endl;
                            if(groups_stats.find(bcd_grp)==groups_stats.end()){
                                groups_stats[bcd_grp] = map< string, map< string, stat > >();
                                groups_stats[bcd_grp][rname] = map< string, stat >();
                                groups_stats[bcd_grp][rname][res_grp] = stat();
                            }
                            else if (groups_stats[bcd_grp].find(rname)==groups_stats[bcd_grp].end())
                                groups_stats[bcd_grp][rname] = map< string, stat >();

                            else if (groups_stats[bcd_grp][rname].find(res_grp)==groups_stats[bcd_grp][rname].end())
                                groups_stats[bcd_grp][rname][res_grp] = stat();

                            data = groups_stats[bcd_grp][rname][res_grp];
                            data.res = rname;
                            count_this=true;
                        }
                    }
                }
                double avg=0, sum=0;
                int N = statistics[i][j][k].size();
                //if (N) {
                    ///////////////
                    double d, stddev=0;
                    if (count_this){
                        for(int l=0; l<N;l++){
                            d = statistics[i][j][k][l];
                            sum += d;

                            data.data.push_back(d);
                        }
                        data.d_sum += sum;
                        data.n += N;
                        groups_stats[bcd_grp][rname][res_grp] = data;
                        //cout << " res name " << rname << " res group " << res_grp << endl;

                    } else
                        for(int l=0; l<N;l++){
                            d = statistics[i][j][k][l];
                            sum += d;
                        }

                if (N){
                    avg = sum/N;
                    for(int l=0; l<N; l++){
                        d = statistics[i][j][k][l];
                        stddev += (d-avg)*(d-avg);
                    }
                    stddev = sqrt(stddev/(double)N);
                    fprintf(output_stat[i], "%4d %8s %4d %8s %8.3lf %9.3lf %10d\n", j, bcd[j].name.c_str(), k, atoms[i][k].name.c_str(), avg, stddev, N);

                }
            }
        }
    }

    for(map<string, FILE*>::iterator it = output_grp_stat.begin(); it!=output_grp_stat.end(); it++){
        fprintf(it->second, " Id1   Label1  Id2   Label2     Dist     Sigma  Frequency\n");
    }

    int ni=0, nj=0;
    if (groupsf.is_open()){
        for(map< string, map< string, map< string, stat > > >::iterator bcd_it = groups_stats.begin(); bcd_it!=groups_stats.end(); bcd_it++){
	    ni=tipos["F6NH"][bcd_it->first];
            for(map< string, map< string, stat > >::iterator res_n = (bcd_it->second).begin(); res_n!=(bcd_it->second).end(); res_n++){
                for(map< string, stat >::iterator res_it = (res_n->second).begin(); res_it!=(res_n->second).end(); res_it++){
                    //cout << "res first " << res_it->first << "res sec " << res_it->second  << endl;
                    
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
                    nj=tipos[rname][res_it->first];
                    fprintf(output_grp_stat[rname], "%4d %8s %4d %8s %8.3lf %9.3lf %10d\n", ni, bcd_it->first.c_str(), nj, res_it->first.c_str(), avg, stddev, n);
                }
            }
        }
    }


    for(int i=0; i<output_stat.size(); i++){
        if (!avg_only) fclose(output[i]);
        fclose(output_stat[i]);
    }
    for(map< string, FILE * >::iterator it=output_grp_stat.begin(); it!=output_grp_stat.end(); it++){
        if (groupsf.is_open()) fclose(it->second);
    }
    printf("Done.\n");

return 0;
}
