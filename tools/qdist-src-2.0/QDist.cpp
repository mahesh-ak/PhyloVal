#include "QDist.hpp"
#include <iostream>

#include "TreeUtil.hpp"
#include "Util.hpp"
#include "Matrix.hpp"


static long CountButterflies(Tree *t);
static void Count(Tree* t1, Tree* t2, long &shared, long &diff);


////////////////////////////////////////////////////////////////////////////////////////////
// The sub-cubic algorithm for quartet distance
//
// described in:
//   A sub-cubic time algorithm for computing the quartet distance between two general trees
//   by Thomas Mailund, Jesper Nielsen and Christian N.S. Pedersen
////////////////////////////////////////////////////////////////////////////////////////////
long SubCubicQDist(Tree* t1, Tree* t2, long &b1, long &b2, long &shared, long &diff) {

    // 1. B
    b1 = CountButterflies(t1);

    // 2. B'
    b2 = CountButterflies(t2);

    // 3. shared_B(T,T') and 
    // 4. diff_B(T,T')
    Count(t1, t2, shared, diff);

    // qdist(T,T') = B + B' - 2*shared_B(T,T') - diff_B(T,T')
    long qdist = b1 + b2 - 2*shared - diff;

    return qdist;
}



/*
 * Calculates the total number of butterflies in tree.
 */
static long CountButterflies(Tree *t) {

    //find leaf set sizes
    std::vector< int > leafSetSizes = TreeUtil::SubtreeLeafSetSizes(t);

    long butterflies = 0;



    //count for every inner node
    for(int ni = 0; ni < t->NumInternalNodes(); ni++) {
        InternalNode* iNode = t->GetInternalNode(ni);
        const std::vector<DirectedEdge*> &edges = iNode->GetEdges();
        const int numSubtrees = edges.size();
        //make sure that the degree of the inner node is at least three
        if(numSubtrees < 3)
            continue;

        //Sum of subtree leaves.
        long S = 0;
        //Sum of squared subtree leaves.
        long S2 = 0;

        for(int i = 0; i < numSubtrees; ++i)
        {
            long subtreeLeaves = leafSetSizes[edges[i]->GetEdgeId()];
            S += subtreeLeaves;
            S2 += subtreeLeaves * subtreeLeaves;
        }

        for(int i = 0; i < numSubtrees; ++i)
        {
            long subtreeLeaves = leafSetSizes[edges[i]->GetEdgeId()];
            butterflies += Util::Choose2(subtreeLeaves)
                * (S*S - S2 - 2*S*subtreeLeaves + 2*subtreeLeaves*subtreeLeaves);
        }
    }

    //divide shared butterflies by four because of symmetry.
    return butterflies / 4;
}



/*
 * Calculates either shared butterflies or both shared and different butterflies.
 */
static void Count(Tree* t1, Tree* t2, long &shared, long &diff) {

    //find shared leaf set sizes
    std::vector< std::vector<int> > sharedLeafSetSizes = TreeUtil::CalcSharedLeafSetSizes(t1, t2);

    //shared_B(T,T')
    double sharedButterflies = 0;
    //diff_B(T,T')
    double differentButterflies = 0;



    Matrix<double> I;
    Matrix<long> Imark;
    Matrix<double> I1markmark;
    Matrix<double> I1markmarkmark;
    Matrix<double> I2markmark;
    Matrix<double> I2markmarkmark;
    std::vector<long> R;
    std::vector<long> C;
    std::vector<long> Rmark;
    std::vector<long> Cmark;
    std::vector<long> Rmarkmark;
    std::vector<long> Cmarkmark;
    std::vector<long> Cmarkmarkmark;
    std::vector<long> Rmarkmarkmark;



    //count for every pair of inner nodes
    for (int n1i = 0; n1i < t1->NumInternalNodes(); n1i++) {
        InternalNode* iNode1 = t1->GetInternalNode(n1i);
        const std::vector<DirectedEdge*> &edges1 = iNode1->GetEdges();
        const int numSubtrees1 = edges1.size();
        //make sure that the degree of the inner node is at least three
        if (numSubtrees1 < 3)
            continue;

        const std::vector<unsigned> &iEdgesIdxs1 = iNode1->GetInternalEdgesIdxs();



        for (int n2i = 0; n2i < t2->NumInternalNodes(); n2i++) {
            InternalNode* iNode2 = t2->GetInternalNode(n2i);
            const std::vector<DirectedEdge*> &edges2 = iNode2->GetEdges();
            const int numSubtrees2 = edges2.size();
            //make sure that the degree of the inner node is at least three
            if (numSubtrees2 < 3)
                continue;

            const std::vector<unsigned> &iEdgesIdxs2 = iNode2->GetInternalEdgesIdxs();


            //matrix containing shared leaf set sizes for subtrees associated with the two inner nodes
            I.resize(numSubtrees1, numSubtrees2);
            //vector containing row sums of I
            R.clear();
            R.resize(numSubtrees1, 0);
            //vector containing column sums of I
            C.clear();
            C.resize(numSubtrees2, 0);

            //sum of all entries
            long M = 0;

            int i;
            int j;
            for (i = 0; i < numSubtrees1; i++)
                for (j = 0; j < numSubtrees2; j++) {
                    int numSharedLeaves = sharedLeafSetSizes[edges1[i]->GetEdgeId()][edges2[j]->GetEdgeId()];
                    I(i, j) = numSharedLeaves;
                    R[i] += numSharedLeaves;
                    C[j] += numSharedLeaves;
                    M += numSharedLeaves;
                }

            //I'
            Imark.resize(numSubtrees1, numSubtrees2);
            //R'
            Rmark.clear();
            Rmark.resize(numSubtrees1, 0);
            //C'
            Cmark.clear();
            Cmark.resize(numSubtrees2, 0);
            long Mmark = 0;

            for (i = 0; i < numSubtrees1; i++)
                for (j = 0; j < numSubtrees2; j++) {
                    long tmp = long(I(i,j)) * (M - R[i] - C[j] + long(I(i,j)));
                    Imark(i, j) = tmp;
                    Rmark[i] += tmp;
                    Cmark[j] += tmp;
                    Mmark += tmp;
                }

            //R''
            Rmarkmark.resize(numSubtrees1);
            for (i = 0; i < numSubtrees1; i++) {
                long sum = 0;
                for (j = 0; j < numSubtrees2; j++)
                    sum += long(I(i,j)) * (C[j] - long(I(i,j)));
                Rmarkmark[i] = sum;
            }
            
            //C''
            Cmarkmark.resize(numSubtrees2);
            for (j = 0; j < numSubtrees2; j++) {
                long sum = 0;
                for (i = 0; i < numSubtrees1; i++)
                    sum += long(I(i,j)) * (R[i] - long(I(i,j)));
                Cmarkmark[j] = sum;
            }

            //count shared butterflies for this pair of inner nodes
            long tmpShared = 0;

            //that means for all pairs of edges going to the inner nodes
            for (unsigned ti = 0; ti < iEdgesIdxs1.size(); ti++) {
                i = iEdgesIdxs1[ti];

                for (unsigned tj = 0; tj < iEdgesIdxs2.size(); tj++) {
                    j = iEdgesIdxs2[tj];

                    if (I(i,j) >= 2) {
                        long tmp = Util::Choose2(long(I(i,j))) *
                            (Mmark - Rmark[i] - Cmark[j] + Imark(i,j)
                             + (long(I(i,j)) - R[i] - C[j]) * (M - R[i] - C[j] + long(I(i,j)))
                             + Rmarkmark[i] - long(I(i,j)) * (C[j] - long(I(i,j))) 
                             + Cmarkmark[j] - long(I(i,j)) * (R[i] - long(I(i,j))));
                        tmpShared += tmp;
                    }
                }
            }

            //add contribution to overall count
            sharedButterflies += tmpShared;


            ////////////////////////////
            //THIS PART FOR DIFF BUTTS
            ////////////////////////////
                
            //number of rows and columns in I
            int rows = numSubtrees1;
            int cols = numSubtrees2;

            //count different butterflies for this pair of inner nodes
            long tmpDiff = 0;


            //there are two ways to calculate and it depends on the shape of I
            //SOLUTION 1: if I is wide and low
            if (rows < cols) {

                //C'''
                Cmarkmarkmark.resize(numSubtrees2);
                for (j = 0; j < numSubtrees2; j++) {
                    long sum = 0;
                    for (i = 0; i < numSubtrees1; i++)
                        sum += long(I(i,j) * I(i,j));
                    Cmarkmarkmark[j] = sum;
                }

                //I1'''
                I1markmark.resize(numSubtrees1, numSubtrees1);
                Matrix<double>::Mult(I, Matrix<double>::NO_TRANSPOSE,
                                     I, Matrix<double>::TRANSPOSE,
                                     I1markmark);
                I1markmarkmark.resize(numSubtrees1, numSubtrees2);
                Matrix<double>::Mult(I1markmark, I, I1markmarkmark);


                //count different butterflies for this pair of inner nodes
                //that means for all pairs of edges going to the inner nodes
                for (unsigned ti = 0; ti < iEdgesIdxs1.size(); ti++) {
                    i = iEdgesIdxs1[ti];

                    for (unsigned tj = 0; tj < iEdgesIdxs2.size(); tj++) {
                        j = iEdgesIdxs2[tj];

                        long tmp = long(I(i,j)) * ((M - R[i] - C[j] + long(I(i,j))) * (R[i] - long(I(i,j))) * (C[j] - long(I(i,j))) 
                                                   + (R[i] - long(I(i,j))) * (long(I(i,j)) * (R[i] - long(I(i,j))) - Cmarkmark[j])
                                                   + (C[j] - long(I(i,j))) * (long(I(i,j)) * (C[j] - long(I(i,j))) - Rmarkmark[i])
                                                   + long(I1markmarkmark(i,j)) - long(I(i,j)) * long(I1markmark(i,i)) - long(I(i,j)) * (Cmarkmarkmark[j] - long(I(i,j) * I(i,j))));

                        tmpDiff += tmp;

                    }
                }
            }

            //SOLUTION 1: if I is wide and low
            else {

                //R'''
                Rmarkmarkmark.resize(numSubtrees1);
                for (i = 0; i < numSubtrees1; i++) {
                    long sum = 0;
                    for (j = 0; j < numSubtrees2; j++)
                        sum += long(I(i,j) * I(i,j));
                    Rmarkmarkmark[i] = sum;
                }

                //I2'''
                I2markmark.resize(numSubtrees2, numSubtrees2);
                Matrix<double>::Mult(I, Matrix<double>::TRANSPOSE,
                                     I, Matrix<double>::NO_TRANSPOSE,
                                     I2markmark);
                I2markmarkmark.resize(numSubtrees1, numSubtrees2);
                Matrix<double>::Mult(I, I2markmark, I2markmarkmark);

                //count different butterflies for this pair of inner nodes
                //that means for all pairs of edges going to the inner nodes
                for (unsigned ti = 0; ti < iEdgesIdxs1.size(); ti++) {
                    i = iEdgesIdxs1[ti];

                    for (unsigned tj = 0; tj < iEdgesIdxs2.size(); tj++) {
                        j = iEdgesIdxs2[tj];

                        long tmp = long(I(i,j)) * ((M - R[i] - C[j] + long(I(i,j))) * (R[i] - long(I(i,j))) * (C[j] - long(I(i,j)))
                                                   + (R[i] - long(I(i,j))) * (long(I(i,j)) * (R[i] - long(I(i,j))) - Cmarkmark[j])
                                                   + (C[j] - long(I(i,j))) * (long(I(i,j)) * (C[j] - long(I(i,j))) - Rmarkmark[i])
                                                   + long(I2markmarkmark(i,j)) - long(I(i,j)) * long(I2markmark(j,j)) - long(I(i,j)) * (Rmarkmarkmark[i] - long(I(i,j) * I(i,j))));

                        tmpDiff += tmp;

                    }
                }
            }

            //add contribution to overall count
            differentButterflies += tmpDiff;
        }
    }

    //make the result permanent
    //divide shared butterflies by four because of symmetry.
    shared = long(sharedButterflies / 4.0);
    //divide different butterflies by four because of multiple pairs of edges
    diff = long(differentButterflies / 4.0);
}

