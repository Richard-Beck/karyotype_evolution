struct parameters{

    float p=0.001;
    float wavelength=1.;
    float centre = 0.;
    float scale = 1.;
    float sd_mutation,mean_mutation;
    float dt = 0.1;
    int Nsteps = 1000;
    int target_output_size=1000;
    string fitness_landscape_type;
    string fitness_landscape_file="not_supplied";
    string output_dir = "output";
    int output_gens=1000;
    int init_size=10000;
    int max_size=2000000;
    float bf=0.01;//,0.05,0.1,0.5};
    int maxchrom=8; // copy number max limit
    int G = 500;
    float p_death=0.0;
    vector<int> init_kary;

    list<vector<int>> peaks;
    vector<float> heights;
    vector<float> sigma;

    //for the polyh landscape
    vector<float> f;
    vector<float> w;
    vector<float> v;

    parameters(string path);
    void read_landscape_file();
    void read_gaussian_file();
    void read_polyh_file();
    void read_random_file();

};

parameters::parameters(string path){

   fstream fin;
   fin.open(path, ios::in);
   std::string tmp, row;
   vector<string> words;
   char delim = ',';

   while(std::getline(fin, row)){
        words.clear();
        stringstream tmp2(row);
        while (std::getline(tmp2, tmp, delim)) {
        words.push_back(tmp);

        }

        if(words[0]=="init_size") init_size=stoi(words[1]);
        if(words[0]=="fitness_landscape_type") fitness_landscape_type=words[1];
        if(words[0]=="fitness_landscape_file") {
            // remove leading whitespace from string filepath:
            string s = words[1];
            s.erase(std::remove_if(s.begin(), s.end(), ::isspace),s.end());
            fitness_landscape_file=s;
        }
        if(words[0]=="output_dir") {
            // remove leading whitespace from string filepath:
            string s = words[1];
            s.erase(std::remove_if(s.begin(), s.end(), ::isspace),s.end());
            output_dir=s;
        }
        if(words[0]=="dt") dt=stof(words[1]);
        if(words[0]=="p") p=stof(words[1]);
        if(words[0]=="wavelength") wavelength=stof(words[1]);
        if(words[0]=="scale") scale=stof(words[1]);
        if(words[0]=="centre") centre=stof(words[1]);
        if(words[0]=="Nsteps") Nsteps=stoi(words[1]);
        if(words[0]=="init_size") init_size=stof(words[1]);
        if(words[0]=="init_kary"){
            for(int i=1; i<words.size(); i++){
                init_kary.push_back(stoi(words[i]));
            }
        }
    }

}

void parameters::read_landscape_file(){
    if(fitness_landscape_type=="gaussian") read_gaussian_file();
    if(fitness_landscape_type=="polyh") read_polyh_file();
    if(fitness_landscape_type=="random") read_random_file();
}
void parameters::read_gaussian_file(){

   fstream fin;
   fin.open(fitness_landscape_file, ios::in);
   std::string tmp, row;
   vector<string> words;
   char delim = ',';
   while(std::getline(fin, row)){

        words.clear();
        stringstream tmp2(row);
        while (std::getline(tmp2, tmp, delim)) {
            words.push_back(tmp);
        }
        int rowsize = words.size();
        vector<int> peak;
        for(int i = 0; i< (rowsize-2); i++){
            peak.push_back(stoi(words[i]));
        }
        peaks.push_back(peak);
        heights.push_back(stof(words[rowsize-2]));
        sigma.push_back(stof(words[rowsize-1]));
    }

}

void parameters::read_random_file(){

   fstream fin;
   fin.open(fitness_landscape_file, ios::in);
   std::string tmp, row;
   vector<string> words;
   char delim = ',';
   while(std::getline(fin, row)){

        words.clear();
        stringstream tmp2(row);
        while (std::getline(tmp2, tmp, delim)) {
            words.push_back(tmp);
        }
        int rowsize = words.size();
        vector<int> peak;
        for(int i = 0; i< rowsize; i++){
            peak.push_back(stoi(words[i]));
        }
        peaks.push_back(peak);
    }

}

void parameters::read_polyh_file(){

   fstream fin;
   fin.open(fitness_landscape_file, ios::in);
   std::string tmp, row;
   vector<string> words;
   char delim = ',';
   while(std::getline(fin, row)){
        // the last row of the input is v. Approach is to clear v and reset it every time through the loop,
        // so at the end of the loop v is correct
        v.clear();
        words.clear();
        stringstream tmp2(row);
        while (std::getline(tmp2, tmp, delim)) {
            words.push_back(tmp);
        }
        int rowsize = words.size();
        vector<int> peak;
        for(int i = 0; i< (rowsize-1); i++){
            peak.push_back(stoi(words[i]));
            v.push_back(stof(words[i]));
        }
        peaks.push_back(peak);
        w.push_back(stof(words[rowsize-1]));
        v.push_back(stof(words[rowsize-1]));
    }

   // remove the last elements (which are v)
   w.pop_back();
   peaks.pop_back();

}


