#include<iostream>
#include<vector>
#include <random>
#include <time.h>
#include <fstream>
#include <sstream>
#include <list>
#include <string>

using namespace std;
#include "parameters.h"

void write_pop(list<vector<int>>& cells, int g, int wmax, string fname){
    string gens = to_string(g);
    gens.insert(gens.begin(), 5 - gens.size(), '0');
    std::ofstream cellfile;
    std::string fnamecell=".csv";
    fnamecell = fname+"gen_" + gens +fnamecell;
    cellfile.open(fnamecell);
    int counter = 0;
    for (const auto& i:cells) {
            if(counter<wmax){
                for(const auto& j : i) cellfile << j << ",";
                cellfile << endl;
            }
        counter++;
    }
    cellfile.close();
}

float selection(vector<int> cell){
    float p_death_1 = 0.;
    float p_death_2 = 0.;
    float p_death = 0;
    for(const auto& chr:cell){
        p_death_1 +=abs(4.0-float(chr))/4.0;
        p_death_2 +=abs(2.0-float(chr))/2.0;
    }

    p_death_2 += 0.43*cell.size();

    p_death = min(p_death_1,p_death_2);

    p_death/=(float)cell.size();
    return p_death;
}

list<vector<int>>::iterator divide(list<vector<int>>& cells, list<vector<int>>::iterator it, int nchrom, float pmis, int maxchrom, mt19937& gen, float p_death){
    uniform_real_distribution<> dis(0.0, 1.0);
    if(p_death>dis(gen)){
        it = cells.erase(it);
        return it;
    }
    bool keep_0=true;
        bool del_i=false;
        vector<int> c0(*it);
        for(auto i=0; i<nchrom; i++){
            if(dis(gen)<pmis){
                int ctot = 2*c0[i];
                int csplit = rand()%ctot;
                c0[i]= csplit;
                if(csplit<1 | csplit>maxchrom){keep_0=false;}
                (*it)[i] = ctot-csplit;
                if(ctot-csplit<1 | ctot-csplit>maxchrom) {del_i = true;}
            }
        }
        if(keep_0) cells.push_front(c0);

        if(del_i) {it = cells.erase(it);} else{++it;}

        return it;
}
int main(int argc , char *argv[])
{
    mt19937 gen(time(NULL));

    parameters p(argv[1]);

    //vector<int> cell(nchrom,2); / diploid
    //vector<int> cell = {4,4,3,3,3,4,6,3,3,3,4,4,4,6,4,3,4,4,4,4,3,4}; //SNU-16
    list<vector<int>> cells;
    int nchrom=p.cell.size();



    for(auto i =0; i<p.init_size; i++){
        cells.push_back(p.cell);
    }
    list<vector<int>>::iterator it;
    uniform_real_distribution<> dis(0.0, 1.0);
    for(auto g = 0; g<p.G; g++){

        it= cells.begin();
        while(it!=cells.end()){
            it = divide(cells, it, nchrom, p.pmis, p.maxchrom, gen,p.p_death);
        }
        if(g%p.output_gens==0){
            cout << "gen: " << g << "; pop. size: " << cells.size() << endl;
            write_pop(cells,g,p.wmax,p.fname);
        }

        if(cells.size()>p.max_size){

            it = cells.begin();
            while(it!=cells.end()){
                if(dis(gen)>p.bf){
                    it = cells.erase(it);
                }else{++it;}
            }
        }
    }

    while(cells.size()<p.max_size){
        it = cells.begin();
        while(it!=cells.end()){
            it = divide(cells, it, nchrom, p.pmis, p.maxchrom, gen,p.p_death);
        }
    }

    //cout << cells.size() << endl;
    write_pop(cells,p.G,p.wmax,p.fname);









    return 0;
}
