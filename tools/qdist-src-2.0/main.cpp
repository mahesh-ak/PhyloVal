
#include <iostream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <set>

#include "Util.hpp"
#include "NewickParser.hpp"
#include "Util.hpp"
#include "Tree.hpp"
#include "TreeUtil.hpp"
#include "QDist.hpp"



static const std::string &ExtractLeafLabel(const LeafNode *n) {
    return n->GetLabel();
}

int main(int argc, char** argv) {

    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << " tree1 tree2" << std::endl;
        std::cout << "  Where:" << std::endl;
        std::cout << "    tree1 and tree2 are files each containing one tree in newic" << std::endl;
        std::cout << "    format. All leaves in the two trees should be labeled, and" << std::endl;
        std::cout << "    the two trees should have the same set of leaves." << std::endl;
        std::cout << std::endl;
        std::cout << "Prints the quartet-distance between tree1 and tree2 and various" << std::endl;
        std::cout << "summary statistics:" << std::endl;
        std::cout << "    N      - The number of leaves in the trees (should be the same for both)." << std::endl;
        std::cout << "    B1     - The number of butterfly quartets in the first tree." << std::endl;
        std::cout << "    B2     - The number of butterfly quartets in the second tree." << std::endl;
        std::cout << "    S      - The number of shared butterfly quartets." << std::endl;
        std::cout << "    D      - The number of different butterfly quartets." << std::endl;
        std::cout << "    Norm B - The normalized shared butterflies, i.e. S / min(B1,B2)." << std::endl;
        std::cout << "    Q      - The quartet distance between the two trees." << std::endl;
        std::cout << "    Norm Q - The normalized quartet distance, i.e. Q / (N choose 4)." << std::endl;
        std::cout << std::endl;
        return 1;
    }

    Tree* tree1;
    Tree* tree2;

    //load newick strings from files
    std::string filename1 = argv[1];
    std::string input1 = Util::LoadFileToString(filename1);

    std::string filename2 = argv[2];
    std::string input2 = Util::LoadFileToString(filename2);

    //parse strings
    NewickParser* parser = new NewickParser();
    tree1 = parser->Parse(input1);
    tree2 = parser->Parse(input2);

    // Check that the leaf lists match
    std::set<std::string> leaves1, leaves2;
    
    std::transform(tree1->GetLeafNodes().begin(), tree1->GetLeafNodes().end(), 
                   std::inserter(leaves1, leaves1.begin()), ExtractLeafLabel);
    std::transform(tree2->GetLeafNodes().begin(), tree2->GetLeafNodes().end(), 
                   std::inserter(leaves2, leaves2.begin()), ExtractLeafLabel);
    
    if (leaves1 != leaves2) {
        std::cout << "The two trees do not have the same leaf sets!" << std::endl;
        return 1;
    }
    
    
    //make identical numbering of the leaves
    TreeUtil::RenumberTreeAccordingToOther(tree2, tree1);    

    TreeUtil::CheckTree(tree1);
    TreeUtil::CheckTree(tree2);

    long n = leaves1.size();
    long max_qdist = Util::Choose(n, 4);
    
    long qdist, b1, b2, shared_b, diff_b;
    qdist = SubCubicQDist(tree1, tree2, b1, b2, shared_b, diff_b);
    
    long min_b = std::min(b1, b2);
    
    std::cout << "N\tB1\tB2\tS\tD\tNorm B\tQ\tNorm Q" << std::endl;
    std::cout << n << '\t' << b1 << '\t' << b2 << '\t' << shared_b << '\t' << diff_b << '\t' << (double(shared_b) / min_b)  << '\t' << qdist << '\t' << (double(qdist) / max_qdist) << std::endl;


	return 0;
}

