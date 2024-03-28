#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <string>

/*
 * Various utility routines
 */

namespace Util {
    void PrintVector(std::string name, std::vector<int> vector);
    void PrintMatrix(std::string name, std::vector<std::vector<int> > vector);
    long Choose2(int n);
    long Choose(int n, int k);
    std::string LoadFileToString(std::string filename);

}

#endif
