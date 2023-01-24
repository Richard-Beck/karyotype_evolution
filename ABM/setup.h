struct parameters{

    float p=0.001;
    float sd_mutation,mean_mutation;
    float dt = 0.1;
    int Nsteps = 1000;
    int target_output_size=1000;
    string fitness_landscape_type, fitness_landscape_file;
    int output_gens=1000;
    int init_size=10000;
    int max_size=100000;
    float bf=0.01;//,0.05,0.1,0.5};
    int maxchrom=8; // copy number max limit
    int G = 500;
    float p_death=0.0;
    vector<int> init_kary;
    parameters(string path);

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
        if(words[0]=="fitness_landscape_file") fitness_landscape_file=words[1];
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
