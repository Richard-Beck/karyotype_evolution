// generic fitness landscape that should be able to hold gaussian/polyh/random_field data
// and access the appropriate methods.
struct fitness_landscape{

    string landscape_type = "gaussian";

    //used in all landscapes
    list<vector<int>> peaks;

    //gaussian landscape
    vector<float> heights;
    vector<float> sigmas;

    //polyh
    int k;
    vector<float> f;
    vector<float> w;
    vector<float> v;

    fitness_landscape() = default;

    // initializers for gaussian landscape
    void Init(string,list<vector<int>>, vector<float>,vector<float>);
    void Init(string,vector<int>& , int , int , mt19937&);

    //initializers for polyharmonic interpolator
    void Init(string,list<vector<int>>&, vector<float>&, vector<float>&, vector<float>&, int);


    float rbf(float);

    // fitness function
    float get_gaussian_fitness(vector<int>&);
    float get_poly_fitness(vector<int>&);
    float get_fitness(vector<int>&);
};


void fitness_landscape::Init(string type,list<vector<int>> p0, vector<float> h0, vector<float> s0){

    landscape_type = type;
    peaks=p0;
    heights=h0;
    sigmas = s0;
}

void fitness_landscape::Init(string type,vector<int>& k1, int npeaks, int ndiff, mt19937& gen){


    landscape_type = type;
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
        peaks.push_back(new_pk);
        for(const auto xx:new_pk) cout << xx << ",";
        cout << endl;
        float sigmax =  dis(gen)*5;
        height += dis(gen)*0.5;
        heights.push_back(height);
        sigmas.push_back(sigmax);
    }
}

void fitness_landscape::Init(string type,list<vector<int>>& knots, vector<float>& fvals,
                                                     vector<float>& wvals,vector<float>& vvals, int kk){

    landscape_type = type;
    peaks=knots;
    k=kk;
    f=fvals;
    w=wvals;
    v=vvals;
}

float fitness_landscape::rbf(float r){
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

float fitness_landscape::get_gaussian_fitness(vector<int>& cn){
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

float fitness_landscape::get_poly_fitness(vector<int>& pt){
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

float fitness_landscape::get_fitness(vector<int>& cn){

    if(landscape_type=="gaussian") return get_gaussian_fitness(cn);
    if(landscape_type=="polyh") return get_poly_fitness(cn);
}

struct gaussian_landscape{

    list<vector<int>> peaks;
    vector<float> heights;
    vector<float> sigma;
    gaussian_landscape() = default;
    void Init(list<vector<int>>, vector<float>,vector<float>);
    void Init(vector<int>& , int , int , mt19937&);
    float get_fitness(vector<int>&);
};

void gaussian_landscape::Init(list<vector<int>> p0, vector<float> h0, vector<float> s0){

    peaks=p0;
    heights=h0;
    sigma = s0;
}

void gaussian_landscape::Init(vector<int>& k1, int npeaks, int ndiff, mt19937& gen){


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
        peaks.push_back(new_pk);
        for(const auto xx:new_pk) cout << xx << ",";
        cout << endl;
        float sigmax =  dis(gen)*5;
        height += dis(gen)*0.5;
        heights.push_back(height);
        sigma.push_back(sigmax);
    }
}


float gaussian_landscape::get_fitness(vector<int>& cn){
    int i = 0;
    float ff = 0.0;
    for(const auto pk:peaks){
        float ffi = 0.0;
        for(int j =0; j<cn.size(); j++){
            ffi+=pow((cn[j]-pk[j])/sigma[i],2);
        }
        ffi = heights[i]*exp(-ffi);
        ff = max(ff,ffi);
        i++;
    }
    return(ff);
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



