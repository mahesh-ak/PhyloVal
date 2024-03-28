#ifndef QDIST_H
#define QDIST_H

#include "Tree.hpp"

long SubCubicQDist(Tree* t1, Tree* t2, 
                   long &b1, long &b2,
                   long &shared_butterflies,
                   long &diff_butterflies);

#endif