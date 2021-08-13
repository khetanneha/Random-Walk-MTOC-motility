// For pseudo-random number generators and distributions
#include <random> 
#include <iostream> 
using namespace std;



double getUniRand(){

// Use random_device to generate a seed for Mersenne twister engine.
std::random_device rd;    

// Use Mersenne twister engine to generate pseudo-random numbers.
std::mt19937 engine(rd());

// "Filter" MT engine's output to generate pseudo-random double values,
// **uniformly distributed** on the closed interval [0, 1].
// (Note that the range is [inclusive, inclusive].)

std::uniform_real_distribution<double> dist(0.0 , 1.0);

return dist(engine);
};

/*
double getUniRand(){

// Use random_device to generate a seed for Mersenne twister engine.
std::random_device rd;    

// Use Mersenne twister engine to generate pseudo-random numbers.
std::mt19937 engine(rd());

// "Filter" MT engine's output to generate pseudo-random double values,
// **uniformly distributed** on the closed interval [0, 1].
// (Note that the range is [inclusive, inclusive].)

std::uniform_real_distribution<double> dist(0.0 , 1.0);

// Generate pseudo-random number.
//double x = dist(engine);


double x;
for( int i = 0 ; i < 100000 ; i++){

//cout << x << endl ;
cout << dist(engine) << endl ;
};
return dist(engine);
};
*/

