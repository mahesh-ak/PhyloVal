#include "Util.hpp"

#include <iostream>
#include <iomanip> //for output formatting
#include <fstream>
#include <cmath>
#include <dirent.h> 
#include <cstdlib>


/*
 * Print a vector and its name.
 */
void Util::PrintVector(std::string name, std::vector<int> vector) {
    //print name
    std::cout << name << " (vector)" << std::endl;

    //print vector
    for (unsigned i = 0; i < vector.size(); i++)
        std::cout << std::setw(3) << vector[i];
    std::cout << std::endl;
}

/*
 * Print a two-dimensional vector.
 */
void Util::PrintMatrix(std::string name, std::vector<std::vector<int> >vector) {
    //print name
    std::cout << name << " (two-dimensional vector)" << std::endl;

    //print matrix
    for (unsigned i = 0; i < vector.size(); i++) {
        for (unsigned j = 0; j < vector.size(); j++)
            std::cout << std::setw(3) << vector[i][j];
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

/*
 * Calculate the binomial coefficient, binom(n,2)
 */
long Util::Choose2(int n) {
    long ln = long(n);
    return (ln - 1 + (ln&1)) * (ln >> 1);
}

/*
 * Calculate the binomial coefficient, binom(n,k)
 */
long Util::Choose(int n, int k) {
    long result = 1;
    if(k > n/2)
        k = n-k;

    for(int i = 1; i <= k; ++i)
    {
        long numerator = n - k + i;
        if(numerator % i == 0)
        {
            result *= numerator/i;
        }
        else if(result % i == 0)
        {
            result /= i;
            result *= numerator;
        }
        else
        {
            result *= numerator;
            result /= i;
        }
    }

    return result;
}

/*
 * Load contents of a file into a one-line string
 */
std::string Util::LoadFileToString(std::string filename) {
    std::ifstream infile(filename.c_str());

    if (!infile) {
        std::cerr << "Could not open file: " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

   std::string lines;
   std::string line;
   while (getline(infile, line)) {
      lines += line;
   }

   infile.close();

   return lines;
}
