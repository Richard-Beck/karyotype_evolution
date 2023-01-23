#include <iostream>
#include<vector>
#include<map>
#include<random>
#include<list>
#include<time.h>
#include<fstream>
#include<algorithm>

using namespace std;

#include "fitness_landscapes.h"
#include "clone_logic.h"
#include "write.h"


map<vector<int>,karyotype> m;
list<karyotype> mutants;

int main(int argc, char *argv[])
{

    float p, sd_mutation,mean_mutation;
    float dt = 0.1;
    if(argc<4) {
            cout << "Using default arguments" << endl;
            p = 0.001;
            sd_mutation=0.1;
            mean_mutation=-0.14;
    }else{
            p = stof(argv[1]);
            sd_mutation=stof(argv[2]);
            mean_mutation=stof(argv[3]);
    }

    string write_dir = make_subdir("output");
    write_log(p,sd_mutation,mean_mutation,dt,write_dir);
    int target_output_size=1000;
    mt19937 gen(time(NULL));
    srand(time(NULL));

    // instantiate cell population
    vector<int> k1(10,2);

    // instantiate fitness peaks
    int npeaks = 3;
    int ndiff = 3;
    list<vector<int>> pks;
    vector<float> heights;
    vector<float> sigmas;
    uniform_real_distribution<> dis(0.0, 1.0);
    vector<int> new_pk = k1;
    uniform_int_distribution<> choose_chrom(0, k1.size()-1);
    float height = 1.0;
    for(int i = 0; i<npeaks; i++){

        for(int j = 0; j<ndiff; j++){
            int ms_index = choose_chrom(gen);
            uniform_int_distribution<> choose_cn(1, 2*new_pk[ms_index]);
            new_pk[ms_index] = choose_cn(gen);
        }
        pks.push_back(new_pk);
        for(const auto xx:new_pk) cout << xx << ",";
        cout << endl;
        float sigma =  dis(gen)*5;
        height += dis(gen)*0.5;
        heights.push_back(height);
        sigmas.push_back(sigma);
    }
    // instantiate fitness landscape
    //hoc_landscape f(sd_mutation,mean_mutation);
    gaussian_landscape f(pks,heights,sigmas);

    write_landscape(f,write_dir);

    float f0 = f.get_fitness(k1);
    karyotype c1(k1,10000,f0);
    m[k1]=c1;


    cout << "starting sim..." << endl;
    // iterate over generations
    for(int i=0; i<4000; i++){
        float mean_fitness=0;

        int ncells = 0;
        auto it = m.begin();
        float min_fitness = it->second.fitness;
        float max_fitness = min_fitness;

        // perform cell divisions
        for (auto& [key, value] : m){
            value.divide(p,gen,mutants,dt);
            mean_fitness=(value.n*value.fitness+ncells*mean_fitness)/(ncells+value.n);
            min_fitness = min(min_fitness,value.fitness);
            max_fitness = max(max_fitness,value.fitness);
            ncells+=value.n;
        }
        // assign any mis-segregated cells to their clone

        for(auto cl:mutants){
            // if the clone already exists, just add one cell
            if (auto search = m.find(cl.cn); search != m.end()){
                (search->second).n+=1;
            }else{ //if this is a new clone then generate a new entry in the map
                // the following call to get fitness needs rethought so we can
                // call the same way regardless of the fitness landscape
                //float fitness = f.get_fitness(cl.cn,cl.fitness,gen);
                float fitness = f.get_fitness(cl.cn);
                karyotype c2(cl.cn,1,fitness);

                m[cl.cn]=c2;
            }
        }
        mutants.clear();
        //rank_selection(m,gen);
        //roulette_selection(m,gen,min_fitness,max_fitness);

        cout << "Ncells: " << ncells << "; mean fitness: " << mean_fitness<< "; min fitness: " << min_fitness<< "; max fitness: " << max_fitness << "; Nclones: " << m.size() << endl;

        if(ncells>2000000){
            float p_output = (float)target_output_size/(float)ncells;
            write_pop(m, i, p_output, write_dir, gen);
            for (auto it = m.begin(); it != m.end();){
                    poisson_distribution<> dpois((it->second).n*0.01);
                    (it->second).n = dpois(gen);
                    if((it->second).n==0){
                        it = m.erase(it);
                    }else{++it;}
            }
        }


    }

    return 0;
}
