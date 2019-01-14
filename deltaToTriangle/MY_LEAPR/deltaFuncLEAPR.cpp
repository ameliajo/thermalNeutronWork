#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "/Users/amelia/Dropbox (MIT)/Senior Fall/Masters/LEAPR/leapr/src/contin/contin.h"
#include "/Users/amelia/Dropbox (MIT)/Senior Fall/Masters/LEAPR/leapr/src/discre/discre.h"



template <typename T>
auto runLEAPR( const T& alpha, const T& beta, const T& rho, int nd, 
  T oscE = T(), T oscW = T() ){

  int ntempr, iprint, nphon, mat, npr, iel, ncold, nss, nalpha, nbeta, 
      lat, ni, nka, mss;
  double za, awr, spr, aws, sps, delta, twt, c, tbeta, dka, b7;
  std::vector<double> temp_vec,kappa;

  /*
  std::cout << "Running Amelia's LEAPR" << std::endl;
  std::cout << "  Alpha ";
  for ( auto& a : alpha ) { std::cout << a << " ";}
  std::cout << "\n  Beta  ";
  for ( auto& b : beta ) { std::cout << b << " ";}
  std::cout << "\n  ND    " << nd;
  std::cout << "\n  Osc E ";
  for ( auto& b : oscE ) { std::cout << b << " ";}
  std::cout << "\n  Osc W ";
  for ( auto& b : oscW ) { std::cout << b << " ";}
  std::cout << std::endl;
  std::cout << "\n  Rho   ";
  for ( auto& b : rho ) { std::cout << b << " ";}
  std::cout << std::endl;
  */



  int itemp = 0;
  //nout   = 24;                                                        // Card 1
  //title  = "h in h20, shortened endf model";                          // Card 2
  ntempr = 1;       iprint = 1;      nphon = 100;                       // Card 3
  mat    = 101;     za     = 1001.0;                                    // Card 4
  awr    = 0.99917; spr    = 20.449; npr   = 2;   iel = 0;   ncold = 0; // Card 5
  nss    = 1;       b7     = 1.;     aws   = 1.1; sps = 3.8883; mss = 1; 
                                                                        // Card 6
  nalpha = 5;       nbeta  = 7;      lat   = 0;                         // Card 7
  temp_vec   = { 296.0 };                                               // Card 10 
  delta  = 0.00255;     ni     = 67;                                    // Card 11
                                                                        // Card 12
  twt    = 0.0;     c      = 0.0;    tbeta = 0.5;             // Card 13
  //nd     = 2;                                                           // Card 14
  //oscE  = { 205.0,    0.48};                                          // Card 15
  //oscW   = { 0.166667, 0.333333 };                                    // Card 16
  nka    = 4;       dka    = 0.01;                                      // Card 17
  kappa  = { 0.1, 0.2, 0.4, 0.7 };                                      // Card 18



  // If we are running a case with no delta functions, then we want to make sure
  // that our continuous piece integrates to 1
  if (nd == 0){ tbeta = 1.0; }


  std::cout << "\n  tbeta " << tbeta;
  std::cout << "\n\n" << std::endl;


  double bk = 8.617385e-5;
  double therm = 0.0253;

  // Loop over scatterers and temperatures
  int isecs = 0;
  double arat  = 1.0; // This is for scaling alpha and beta values later
  bool done = false;


  //if ( isecs == 0 ){ std::cout << "Principal scatterer" << std::endl; }
  if ( isecs >  0 ){ 
    //std::cout << "Secondary scatterer" << std::endl;
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



  std::vector<double> t_eff_vec ( temp_vec.size(), 0.0 );
  auto lambda_s_t_eff = contin( itemp, nphon, delta, tbeta, scaling, tev,
                                sc, rho, alpha, beta, sym_sab );
  double lambda_s = std::get<0>(lambda_s_t_eff);
  t_eff_vec[itemp] = std::get<1>(lambda_s_t_eff) * temp;

  if ( oscE.size() > 0 ){
    discre( itemp, sc, scaling, tev, lambda_s, twt, tbeta, alpha,
    beta, temp_vec, oscE, oscW, t_eff_vec, sym_sab );
  }


  return sym_sab;
}


int main() {
  std::string line;
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
      if (thisVec.size() > 0){ v2.push_back(thisVec);}
    }
  }
  std::vector<double> alphaVals = v2[0];
  std::vector<double> betaVals  = v2[1];
  std::vector<double> rhoVals   = v2[2];
  std::vector<double> oscE, oscW;
  if (v2.size() > 3){
    oscE = v2[3];
    oscW = v2[4];
  }
  int nd = oscE.size();
  auto sym_sab = runLEAPR(alphaVals,betaVals,rhoVals,nd,oscE,oscW);
  std::ofstream outFile("outputSAB.txt");
  for (const auto& a : sym_sab){
    for (const auto& b : a){
      for (const auto& c : b){
        outFile << c << " ";
      } 
    }
  }
  return 0;
}


