#ifndef PARAMETERS_H_INCLUDED
#define PARAMETERS_H_INCLUDED



struct parameters{
    int output_gens=1000;
    int init_size=1;
    int max_size=100000;
    float bf=0.01;//,0.05,0.1,0.5};
    int maxchrom=8; // copy number max limit
    int G = 200;
    float pmis=0.001;//,0.0001,0.00001};//{0.01,0.001,0.0001};
    int wmax = 1000; // max no of cells to write to file
    float p_death=0.0;
    string fname="set_fname";
    vector<int> cell;
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
        if(words[0]=="max_size") max_size=stoi(words[1]);
        if(words[0]=="maxchrom") maxchrom=stoi(words[1]);
        if(words[0]=="bf") bf=stof(words[1]);
        if(words[0]=="G") G=stoi(words[1]);
        if(words[0]=="pmis") pmis=stof(words[1]);
        if(words[0]=="wmax") wmax=stoi(words[1]);
        if(words[0]=="p_death") p_death=stof(words[1]);
        if(words[0]=="fname") fname=words[1];
        if(words[0]=="output_gens") output_gens=stoi(words[1]);
        if(words[0]=="cell"){
            for(int i=1; i<words.size(); i++){
                cell.push_back(stoi(words[i]));
            }
        }
    }
}

#endif // PARAMETERS_H_INCLUDED
