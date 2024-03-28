#include "NewickParser.hpp"
#include "Util.hpp"
#include "TreeUtil.hpp"
#include "QDist.hpp"

#include <cstdlib>
#include <iostream>
#include <sstream>



template<typename T>
std::string toString(const T& x)
{
    std::ostringstream oss;
    oss << x;
    return oss.str();
}



void testTrees(Tree* tree1, Tree* tree2)
{
#if 0 // FIXME: left over from when we implemented more of the algorithms.
    long result1 = QDist::SubCubicQDist(tree1, tree2);
    long result2 = QDist::CubicQDist(tree1, tree2);
    long result3 = QDist::QuarticQDist(tree1, tree2);

    bool fail = false;

    if(result1 != result2)
    {
        std::cout << "Sub-cubic and cubic qdists disagree." << std::endl;
        fail = true;
    }

    if(result1 != result3)
    {
        std::cout << "Sub-cubic and quartic qdists disagree." << std::endl;
        fail = true;
    }

    if(result2 != result3)
    {
        std::cout << "Cubic and quartic qdists disagree." << std::endl;
        fail = true;
    }

    if(fail)
        exit(-1);
#endif
}



int main(int argc, char** argv) {

    Tree* tree1;
    Tree* tree2;

    const std::string FILE_PREFIX = "testdata/small";
    const std::string FILE_SUFFIX = ".tree";
    const unsigned N_FILES = 5;

    NewickParser* parser = new NewickParser();

    for(unsigned i = 1; i <= N_FILES; ++i)
    {
        std::string filename1 = FILE_PREFIX + toString(i) + FILE_SUFFIX;
        std::string input1 = Util::LoadFileToString(filename1);
        tree1 = parser->Parse(input1);

        for(unsigned j = 1; j <= N_FILES; ++j)
        {
            std::string filename2 = FILE_PREFIX + toString(j) + FILE_SUFFIX;
            std::string input2 = Util::LoadFileToString(filename2);
            tree2 = parser->Parse(input2);

            TreeUtil::RenumberTreeAccordingToOther(tree2, tree1);    
    
            TreeUtil::CheckTree(tree1);
            TreeUtil::CheckTree(tree2);

            testTrees(tree1, tree2);
        }
    }

	return 0;
}
