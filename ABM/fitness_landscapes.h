

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

