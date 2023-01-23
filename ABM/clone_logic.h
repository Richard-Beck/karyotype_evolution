

struct karyotype{
    int n=0;
    vector<int> cn;
    float fitness=1;
    karyotype()=default;
    karyotype(vector<int>&,int,float);
    void divide(float,mt19937&,list<karyotype>&);
    void divide(float,mt19937&,list<karyotype>&,float);
};

karyotype::karyotype(vector<int> &founder, int N, float f){

    cn=founder;
    fitness = f;
    n=N;

}

// all cells divide
void karyotype::divide(float p, mt19937& gen, list<karyotype>& mutants){

    binomial_distribution<> d1(cn.size(), p);

    for(int i=0;i<n;i++){
        int n_mis = d1(gen);
        if(n_mis>0){
            //if misseg happens, decrement cell count by 1
            n--;
            // vectors to store new daughters:
            vector<int> d1 = cn;
            vector<int> d2 = cn;
            // if any chromosome copy numbers go to zero, daughters not valid.
            bool d1_valid=true;
            bool d2_valid=true;
            // following is not strictly accurate since missegregations can be overwritten.
            // may be necessary to fix this.
            for(int j=0;j<n_mis;j++){
                int m = rand()%d1.size(); // choose a chromosome at random
                // following is not accurate since if cn_m1==cn[m], no mis-segregation has occurred.
                int cn_m1 = rand()%(2*cn[m]); // choose a mis-segregated copy number
                int cn_m2 = 2*cn[m]-cn_m1; // second daughter gets the inverse

                // set invalid flags if invalid cn generated:
                if(cn_m1==0) d1_valid = false;
                if(cn_m2==0) d2_valid = false;

                d1[m] = cn_m1;
                d2[m] = cn_m2;
            }
            if(d1_valid){
                    karyotype tmp(d1,1,fitness);
                    mutants.push_back(tmp);
            }
            if(d2_valid){
                    karyotype tmp(d2,1,fitness);
                    mutants.push_back(tmp);
            }
        };
    }
    n*=2; // all the cells that didn't missegregate, divide.


}

//in overloaded version, karyotype fitness is interpreted as net growth rate and
// only fitness*dt*n get to divide. This does have problems though because using the
// net growth rate means we will underestimate the number of divisions.
void karyotype::divide(float p, mt19937& gen, list<karyotype>& mutants, float dt){

    float p_div = dt*abs(fitness);
    binomial_distribution<> d0(n, p_div);
    int n_divs = d0(gen);
    if(fitness<0){
        n-=n_divs;
        return;
    }

    binomial_distribution<> d1(cn.size(), p);

    for(int i=0;i<n_divs;i++){
        int n_mis = d1(gen);
        if(n_mis>0){
            //if misseg happens, decrement cell count by 1
            n--;
            // vectors to store new daughters:
            vector<int> d1 = cn;
            vector<int> d2 = cn;
            // if any chromosome copy numbers go to zero, daughters not valid.
            bool d1_valid=true;
            bool d2_valid=true;
            // following is not strictly accurate since missegregations can be overwritten.
            // may be necessary to fix this.
            for(int j=0;j<n_mis;j++){
                int m = rand()%d1.size(); // choose a chromosome at random
                // following is not accurate since if cn_m1==cn[m], no mis-segregation has occurred.
                int cn_m1 = rand()%(2*cn[m]); // choose a mis-segregated copy number
                int cn_m2 = 2*cn[m]-cn_m1; // second daughter gets the inverse

                // set invalid flags if invalid cn generated:
                if(cn_m1==0) d1_valid = false;
                if(cn_m2==0) d2_valid = false;

                d1[m] = cn_m1;
                d2[m] = cn_m2;
            }
            if(d1_valid){
                    karyotype tmp(d1,1,fitness);
                    mutants.push_back(tmp);
            }
            if(d2_valid){
                    karyotype tmp(d2,1,fitness);
                    mutants.push_back(tmp);
            }
        }
        else{
                n++;
        }
    }



}



void rank_selection(std::map<vector<int>,karyotype> &m, mt19937& gen){

        // example:
        // say there are 5 clones in m and their fitness in order of appearance in m is:
        // 4,1,2,9,3
        // we make a vector of pairs with the first item in the pair fitness, second item their order in m:
        // {{4,1}, {1,2}, {2,3}, {9,4}, {3,5}}
        // now sort them:
        // x_sort = {{1,2}, {2,3},{3,5},{4,1},{9,4}}
        // if we iterate along the i'th element of x_sort using:
        // x_rank[x_sort[i]] = i;
        // we get x_rank = {4,1,2,5,3};

        // pair vector stores fitness and the order of each clone in m
        vector<pair<float,int>> fitness_rank;
        // ranks holds the indexes of m ranked in fitness orders
        vector<int> ranks(m.size(),0);

        // populate fitness_rank vector
        int counter = 0;
        for (auto it = m.begin(); it != m.end();){
            fitness_rank.push_back(make_pair((it->second).fitness, counter));
            it++;
            counter++;
        }

        // sort by fitness
        sort(fitness_rank.begin(), fitness_rank.end());

        // populate rank vector
        for(int a = 0; a<fitness_rank.size(); a++){
            ranks[fitness_rank[a].second] = a;
        }

        // kill unfit cells based on rank, since we know
        // m.begin()+counter has rank ranks[counter]
        counter = 0;
        for (auto it = m.begin(); it != m.end();){
            float p_survival = (float)ranks[counter]/(float)m.size();
            binomial_distribution<> d1((it->second).n, p_survival);
            (it->second).n = d1(gen);
            if((it->second).n<1){
                it = m.erase(it);
            }else{
                ++it;
            }
            counter++;
        }

}

void roulette_selection(std::map<vector<int>,karyotype> &m, mt19937& gen, float min_fitness, float max_fitness){
        for (auto it = m.begin(); it != m.end();){
            float p_survival = ((it->second).fitness-min_fitness)/(max_fitness-min_fitness);
            p_survival = min(p_survival,(float)0.9);
            p_survival = max((float)0.1,p_survival);
            binomial_distribution<> d1((it->second).n, p_survival);
            (it->second).n = d1(gen);
            if((it->second).n<1){
                it = m.erase(it);
            }else{
                ++it;
            }
        }

}

