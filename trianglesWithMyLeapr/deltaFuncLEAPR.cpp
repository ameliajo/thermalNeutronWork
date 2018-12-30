
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>



int main() {
  std::string line;
  int i = 0;
  std::vector<std::string> v;
  std::vector<std::vector<double>> v2;
  std::ifstream fin;
  double val;
  fin.open( ( "inputVals.txt" ) );
  if ( fin.is_open()) {
    while ( getline ( fin, line )) {
      std::vector<double> thisVec {};
      std::stringstream ss ( line );
      while (ss >> val){
        thisVec.push_back(val);
      }
      //for ( auto x : thisVec ){ std::cout << x << std::endl; }
      if (thisVec.size() > 0){ v2.push_back(thisVec);}
      //if ( getline ( ss, line, ',')) {
      //  v.push_back( line );
     // }
    }
    for ( auto x : v2[0] ){ std::cout << x << std::endl; }
    for ( auto x : v2[1] ){ std::cout << x << std::endl; }
    //for ( auto x : v2[2] ){ std::cout << x << std::endl; }
  }
  return 0;
}

  /*
std::vector<double> vec; 
std::vector<std::string> y;
std::ifstream file("inputVals.txt");
if (file.is_open()) {
    std::string line;
    while (getline(file, line)) {
      
        std::cout << line << std::endl;
        printf("%s", line.c_str());
        std::cout << '\n' << std::endl;
        y.push_back(line);
    }
    file.close();
}
    //std::vector<double> alphaVals;
    //for (std::string x : y){
    //  std::cout << x << std::endl;
    //  alphaVals.push_back(stod(x));
   // }

return 0;
}
*/
  /*
  std::vector<double> rainfall;    // a vector to hold rainfall data
  std::ifstream inputFile("inputVals.txt");
  if (inputFile) {
      double val1 = 0; 
      while ( inputFile >> val1) {
        //std::cout << val1 << std::endl;
        rainfall.push_back(val1);
  }
  for ( auto x : rainfall){ std::cout << x << std::endl; }
  }
  
  return 0;
}
*/


     // close the file


      /*
  std::vector<std::string> my_arr;
  std::ifstream dict_file("inputVals.txt");
  std::string line;

  while(std::getline(dict_file, line)) {
     new_line;
    new_line = line + "\n";
    my_arr.push_back(new_line);
  }

  //std::vector<double> alphaVals = my_arr[0];
  //std::vector<double> betaVals = my_arr[1];
  //std::vector<double> rhoVals = my_arr[2];
  //for ( auto x : alphaVals ){ std::cout << x << std::endl; }
  return 0;
}



*/

/*
//#include "leapr.cpp"
#include <iostream>
#include "/Users/amelia/Dropbox (MIT)/Senior Fall/Masters/LEAPR/leapr/src/contin/contin.h"
#include <fstream>


int main() {
  std::cout << "Hello world!" << std::endl;


  int ntempr, iprint, nphon, mat, npr, iel, ncold, nss, nalpha, nbeta, 
      lat, ni, nd, nka, mss;
  double za, awr, spr, aws, sps, delta, twt, c, tbeta, dka, b7;
  std::vector<double> alpha, beta, temp_vec, rho, oscE, oscW, 
    kappa;


  int itemp = 0;
  //nout   = 24;                                                 // Card 1
  //title  = "h in h20, shortened endf model";                   // Card 2
  ntempr = 1;       iprint = 1;      nphon = 100;              // Card 3
  mat    = 101;     za     = 1001.0;                           // Card 4
  awr    = 0.99917; spr    = 20.449; npr   = 2;   iel = 0;   ncold = 0; // Card 5
  nss    = 1;       b7     = 1.;     aws   = 1.1; sps = 3.8883; mss = 1; 
                                                                     // Card 6
  nalpha = 5;       nbeta  = 7;      lat   = 1;                // Card 7
  alpha  = { 0.01008, 0.015, 0.0252, 0.033, 0.050406 };        // Card 8
  beta   = { 0.00000, 0.006375, 0.012750, 0.025500, 0.038250, 0.05100, 0.065750 }; // Card 9
  temp_vec   = { 296.0 };                                          // Card 10 
  delta  = 0.00255;     ni     = 67;                           // Card 11
  rho    = { 0.00000, 0.00050, 0.00100, 0.00200, 0.00350, 0.00500, 
    0.00750, 0.01000, 0.01300, 0.01650, 0.02000, 0.02450, 0.02900, 
    0.03400, 0.03950, 0.04500, 0.05060, 0.05620, 0.06220, 0.06860, 
    0.07500, 0.08300, 0.09100, 0.09900, 0.10700, 0.11500, 0.11970, 
    0.12140, 0.12180, 0.11950, 0.11250, 0.10650, 0.10050, 0.09542, 
    0.09126, 0.08710, 0.08390, 0.08070, 0.07798, 0.07574, 0.07350, 
    0.07162, 0.06974, 0.06804, 0.06652, 0.06500, 0.06340, 0.06180, 
    0.06022, 0.05866, 0.05710, 0.05586, 0.05462, 0.05350, 0.05250, 
    0.05150, 0.05042, 0.04934, 0.04822, 0.04706, 0.04590, 0.04478, 
    0.04366, 0.04288, 0.04244, 0.04200, 0. };
                                                                   // Card 12
  twt    = 0.055556;     c      = 0.0;    tbeta = 0.444444;    // Card 13
  nd     = 2;                                                  // Card 14
  //oscE  = { 205.0,    0.48};                                  // Card 15
  oscE   = { 35.8,    0.48};                                  // Card 15
  oscW   = { 0.166667, 0.333333 };                             // Card 16
  nka    = 4;       dka    = 0.01;                             // Card 17
  kappa  = { 0.1, 0.2, 0.4, 0.7 };                             // Card 18

  double bk = 8.617385e-5;
  double therm = 0.0253;

  // Loop over scatterers and temperatures
  int isecs = 0;
  double arat  = 1.0; // This is for scaling alpha and beta values later
  bool done = false;


  if ( isecs == 0 ){ std::cout << "Principal scatterer" << std::endl; }
  if ( isecs >  0 ){ 
    std::cout << "Secondary scatterer" << std::endl;
    arat = aws / awr; 
  }
 
  double temp = temp_vec[itemp];
  double tev = bk * temp;
  double sc = 1.0;
  if ( lat == 1 ){ sc = therm/tev; }
  double scaling = sc/arat;

 
  std::vector<std::vector<std::vector<double>>> sym_sab( alpha.size(),
    std::vector<std::vector<double>> (beta.size(), 
    std::vector<double> ( ntempr, 0.0 ) ) );



  auto lambda_s_t_eff = contin( itemp, nphon, delta, tbeta, scaling, tev,
                                sc, rho, alpha, beta, sym_sab );
  double lambda_s = std::get<0>(lambda_s_t_eff);
  std::cout << lambda_s << std::endl;




  return 0;
}
*/
