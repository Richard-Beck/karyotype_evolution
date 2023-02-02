#include <iostream>
#include<vector>
#include <cctype>
#include<map>
#include<random>
#include<list>
#include<time.h>
#include<fstream>
#include <sstream>
#include<algorithm>

using namespace std;

#include "setup.h"
#include "fitness_landscapes.h"
#include "clone_logic.h"
#include "write.h"
#include "polyh.h"



map<vector<int>,karyotype> m;
list<karyotype> mutants;

int main(int argc, char *argv[])
{

    float sd_mutation,mean_mutation;


    string config_file_path = "config/test_A.txt";
    //config_file_path = "output/secondTest_rep_001/00000/proc_data/config.txt";

    if(argc<2) {
        cout << "using default config file (hopefully one is there...)" << endl;
    }else{
            config_file_path = argv[1];
    }
    parameters par(config_file_path);

    string write_dir = make_subdir(par.output_dir);
    write_log(par.p,0,0,par.dt,write_dir);

    mt19937 gen(time(NULL));
    srand(time(NULL));

    // instantiate cell population
    vector<int> k1=par.init_kary;

    // instantiate fitness peaks

    cout << "initializing landscape...";
    // instantiate fitness landscape
    //hoc_landscape f(par.sd_mutation,par.mean_mutation);
    fitness_landscape *f;
    list<vector<int>> peaks;
    gaussian_landscape gl(peaks);
    polyharmonic_landscape pl(peaks);
    cout << par.fitness_landscape_file << endl;
    if(par.fitness_landscape_file=="not supplied"){
        peaks.push_back(k1);
        vector<float> heights = {1.0};
        vector<float> sigmas = {5.0};
        f = &gl;
        f->init(peaks, heights, sigmas);
    }else{
        par.read_landscape_file();
        if(par.fitness_landscape_type=="gaussian"){
            f = &gl;
            f->init(par.peaks, par.heights, par.sigma);
        }
        if(par.fitness_landscape_type=="polyh"){
            f = &pl;
            f->init(par.peaks, par.f, par.w,par.v,2);
        }
    }
    cout << f->get_fitness(k1) << endl;
    cout << "complete" << endl;

    //if(par.fitness_landscape_type=="gaussian") write_landscape(f,write_dir);

    float f0 = f->get_fitness(k1);
    karyotype c1(k1,par.init_size,f0);
    m[k1]=c1;


    cout << "starting sim..." << endl;
    // iterate over generations
    for(int i=0; i<par.Nsteps; i++){
        float mean_fitness=0;

        int ncells = 0;
        auto it = m.begin();
        float min_fitness = it->second.fitness;
        float max_fitness = min_fitness;

        // perform cell divisions
        for (auto& [key, value] : m){
            value.divide(par.p,gen,mutants,par.dt);
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
                float fitness = f->get_fitness(cl.cn);
                karyotype c2(cl.cn,1,fitness);

                m[cl.cn]=c2;
            }
        }
        mutants.clear();
        //rank_selection(m,gen);
        //roulette_selection(m,gen,min_fitness,max_fitness);

        cout << "Ncells: " << ncells << "; mean fitness: " << mean_fitness<< "; min fitness: " << min_fitness<< "; max fitness: " << max_fitness << "; Nclones: " << m.size() << endl;

        if(ncells>par.max_size){
            float p_output = (float)par.target_output_size/(float)ncells;
            write_pop(m, i, p_output, write_dir, gen);
            for (auto it = m.begin(); it != m.end();){
                    poisson_distribution<> dpois((it->second).n*par.bf);
                    (it->second).n = dpois(gen);
                    if((it->second).n==0){
                        it = m.erase(it);
                    }else{++it;}
            }
        }


    }

    return 0;
}
