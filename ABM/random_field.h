struct random_field{

    vector<float> centers;
    float wavelength;
    random_field() = default;
    void Init(list<vector<float>>, float);
    void Init(vector<int>& , int , float , mt19937&);
    float get_fitness(vector<int>&);

};

random_field::Init(list<vector<float>> c, float w){
        centers = c;
        wavelength=w;
}

random_field::Init(vector<int>& k0, int npeaks, float w, mt19937&){
        wavelength=w;
        uniform_real_distribution<> dis(0.0, 1.0);
        for(int i = 0; i<npeaks; i++){
            vector<float> p;
            for(int j = 0; j<k0.size(); j++){
                p.push_back(dis(gen));
            }
            centers.push_back(p);
        }
}

random_field::get_fitness(vector<int>& cn){
}
