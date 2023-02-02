// generic fitness landscape that should be able to hold gaussian/polyh/random_field data
// and access the appropriate methods.
class fitness_landscape{
    public:
        //used in all landscapes
        list<vector<int>> peaks;
        fitness_landscape(list<vector<int>>& p){
            peaks = p;
        }
        virtual float get_fitness(vector<int>&)=0;

        //overloaded initializer signatures for various landscapes
        virtual void init(list<vector<int>>&,vector<float>&,vector<float>&){};
        virtual void init(list<vector<int>>&,vector<float>&, vector<float>&,vector<float>&, int){};
};

class gaussian_landscape: public fitness_landscape{
    public:
        vector<float> heights;
        vector<float> sigmas;
        gaussian_landscape(list<vector<int>>& p):fitness_landscape(p) {}
        void init(list<vector<int>>& p, vector<float>& h,vector<float>& s){
            peaks=p;
            heights=h;
            sigmas=s;
        };
        float get_fitness(vector<int>&);

};

float gaussian_landscape::get_fitness(vector<int>& cn){
    int i = 0;
    float ff = 0.0;
    for(const auto pk:peaks){
        float ffi = 0.0;
        for(int j =0; j<cn.size(); j++){
            ffi+=pow((cn[j]-pk[j])/sigmas[i],2);
        }
        ffi = heights[i]*exp(-ffi);
        ff = max(ff,ffi);
        i++;
    }
    return(ff);
}

class polyharmonic_landscape: public fitness_landscape{
    public:
        int k;
        vector<float> f;
        vector<float> w;
        vector<float> v;
        polyharmonic_landscape(list<vector<int>>& p):fitness_landscape(p) {

        }
        void init(list<vector<int>>& p, vector<float>& ff, vector<float>& ww,
                               vector<float>& vv, int kk){
            peaks=p;
            k=kk;
            f=ff;
            v=vv;
            w=ww;
        }
        float rbf(float r);
        float get_fitness(vector<int>&);

};

float polyharmonic_landscape::rbf(float r){
    float val;
    if(r>=1){
        if(k%2==0){
            val=pow(r,k)*log(r);
        }else{
            val = pow(r,k);
        }
    }else{
        if(k%2==0){
            val=pow(r,k-1)*log(pow(r,r));
        }else{
            val = pow(r,k);
        }
    }
    return val;
}

float polyharmonic_landscape::get_fitness(vector<int>& pt){
    vector<int> p2;
    p2.push_back(1);
    for(int i=0; i<pt.size(); i++) p2.push_back(pt[i]);

    float fitness = 0;
    vector<float> r; // consider pre-allocating
    for(const auto knot_i:peaks){
        float d_i = 0;
        for(int j = 0; j<pt.size(); j++){
            d_i+=pow(knot_i[j]-pt[j],2);
        }
        r.push_back(rbf(sqrt((float)d_i)));
    }
    for(int i = 0; i<r.size(); i++){
        fitness+=w[i]*r[i];
    }
    for(int i = 0; i<p2.size(); i++){
        fitness+=v[i]*p2[i];
    }
    return fitness;
}

struct hoc_landscape{

    float sigma =1;
    float mean =1.;
    float f0=1.0;
    map<vector<int>,float> visited;

    hoc_landscape() = default;
    hoc_landscape(float s, float m){
        sigma=s;
        mean=m;
        }
    float get_fitness(vector<int>&);
    float get_fitness(vector<int>&, float, mt19937&);
};

float hoc_landscape::get_fitness(vector<int>& cn){
    return f0;
}


float hoc_landscape::get_fitness(vector<int>& cn, float f, mt19937& gen){
    if (auto search = visited.find(cn); search != visited.end()){
            return(search->second);
            }else{ //if this is a new clone then generate a new entry in the map
                normal_distribution<> d{f+mean,sigma};
                f=d(gen);
                visited[cn]=f;
            }
    return f;
}



