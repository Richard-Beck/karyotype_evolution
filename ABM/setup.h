struct parameters{

    float p=0.001;
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

    parameters(string path);
    void read_landscape_file();

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


