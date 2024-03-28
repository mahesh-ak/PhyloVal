#include "NewickParser.hpp"
#include <iostream>
#include <cstdlib>

/*
 * A recursive descend parser for the Newick format. It builds a tree data
 * structure, as specified in Tree.h
 * 
 * It parses strings built from the grammar specified below. Each of the
 * six first rules are handled by a separate function. The two last rules
 * are handled implicitly. The parser handles each of the examples given
 * below, which are found at the Wikipedia page on Newick format:
 *
 *    http://en.wikipedia.org/wiki/Newick_format
 * 
 *
 * WARNING:
 *   to work with the qdist algorithms, leaves need to named
 *
 *
 * GRAMMAR:
 *
 *   Tree      --> Subtree ";" | Branch ";"
 *   Subtree   --> Leaf | Internal
 *   Leaf      --> Name
 *   Internal  --> "(" BranchSet ")" Name
 *   BranchSet --> Branch | BranchSet "," Branch
 *   Branch    --> Subtree Length
 *   Name      --> empty | string
 *   Length    --> empty | ":" number
 *
 *
 * EXAMPLES:
 *
 *   (,,(,));                               no nodes are named
 *   (A,B,(C,D));                           leaf nodes are named
 *   (A,B,(C,D)E)F;                         all nodes are named
 *   (:0.1,:0.2,(:0.3,:0.4):0.5);           all but root node have a distance to parent
 *   (:0.1,:0.2,(:0.3,:0.4):0.5):0.0;       all have a distance to parent
 *   (A:0.1,B:0.2,(C:0.3,D:0.4):0.5);       distances and leaf names (popular)
 *   (A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;     distances and all names
 *   ((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;    a tree rooted on a leaf node (rare)
 */

NewickParser::NewickParser() {

}

NewickParser::~NewickParser() {
}

/*
 * Helper function. Removes all white-space from a string
 */
std::string trim(std::string string) {
    int whitespaceSize = 3;
    char whitespace[] = {' ', '\t', '\n'};

    for (int i = 0; i < whitespaceSize; i++) {
        char c = whitespace[i];
        std::string::size_type loc = 0;
        loc = string.find(c, 0);
        while (loc != std::string::npos) {
            string = string.erase(loc,1);
            loc = string.find(c, 0);
        }
    }

    return string;
}

/*
 * Helper function.
 * Find the index of the right parenthesis that matches the left parenthesis that is the first char of this string
 * Return the index or -1 if no parenthesis matches
 */
int findMatchingRightParenthesis(std::string string) {
    if (string[0] != '(') {
        std::cerr << "NewickParser ERROR: expected string with left parenthesis as first char." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string::size_type index = 0;
    int parenthesesDepth = 0;
    while (index < string.size()) {
        if (string[index] == '(')
            parenthesesDepth++;
        if (string[index] == ')') {
            if (parenthesesDepth == 1)
                return index;
            else
                parenthesesDepth--;
        }
        index++;
    }

    //no matching parenthesis found
    return -1;
}

/*
 * First rule of the grammar.
 * Tree --> Subtree ";" rule
 */
Tree* NewickParser::Parse(std::string string) {
    //reset state
    ResetParser();

    //remove whitespace
    string = trim(string);
    //make sure the string ends with semicolon
    int lastIndex = string.size() - 1;
    if (string[lastIndex] == ';')
        string.erase(lastIndex);

    Node* root = ParseSubtree(string);

    Tree* tree = new Tree();
    tree->SetRoot(root);
    tree->SetInternalNodeList(GetInternalNodeList());
    tree->SetLeafNodeList(GetLeafNodeList());
    tree->SetEdgeList(GetDirectedEdgeList());

    return tree;
}

/*
 * Subtree --> Leaf | Internal
 */
Node* NewickParser::ParseSubtree(std::string string) {
    //std::cout << "ENTER subtree < " << string << " >" << std::endl;
    Node* subtree;

    if (string.empty())
        subtree = ParseLeafNode(string);
    else if (string[0] == '(')
        subtree = ParseInternalNode(string);
    else
        subtree = ParseLeafNode(string);

    return subtree;
}

/*
 * Internal --> "(" BranchSet ")" Name
 */
Node* NewickParser::ParseInternalNode(std::string string) {
    //std::cout << "ENTER internal node < " << string << " >"  << std::endl;

    if (string[0] != '(') {
        std::cerr << "NewickParser ERROR: Expected ( as first char in internal node." << std::endl;
        exit(EXIT_FAILURE);
    }

    //find last right parenthesis
    std::string::size_type indexOfLastRightParen = string.find_last_of(')');
    if (indexOfLastRightParen == std::string::npos) {
        std::cerr << "NewickParser ERROR: Internal node does not have )." << std::endl;
        exit(EXIT_FAILURE);
    }

    //find branch set
    int branchSetStringLength = indexOfLastRightParen - 1;
    std::vector<DirectedEdge*> branchSet = ParseBranchSet(string.substr(1,branchSetStringLength));


    //remaining number of chars after the right parenthesis
    int remainingChars = string.size()-1  - indexOfLastRightParen;
    //find name
    std::string name;
    if (remainingChars == 0)
        name = "Internal NONAME";
    else
        name = string.substr(indexOfLastRightParen + 1, remainingChars);
    
    //construct the internal node
    InternalNode* internalNode = new InternalNode(name, GetNextInternalId());
    //and add all the branches
    std::string::size_type i;
    for (i = 0; i < branchSet.size(); i++) {
        DirectedEdge* edge = branchSet[i];

        edge->SetFromNode(internalNode);

        DirectedEdge* backEdge = new DirectedEdge(GetNextEdgeId());
        //add edge to parser list
        AddDirectedEdge(backEdge);
        backEdge->SetFromNode(edge->GetToNode());
        backEdge->SetToNode(internalNode);
        internalNode->AddEdge(edge);
        edge->GetToNode()->AddEdge(backEdge);

        edge->SetBackEdge(backEdge);
        backEdge->SetBackEdge(edge);
    }

    //add node to parser list
    AddInternalNode(internalNode);
    return internalNode;
}

/*
 * Leaf --> Name
 */
Node* NewickParser::ParseLeafNode(std::string string) {
    //std::cout << "ENTER leaf node < " << string << " >"  << std::endl;

    std::string name = string;
    if (string.empty())
        name = "Leaf NONAME";

    LeafNode* leafNode = new LeafNode(name, GetNextLeafId());
    //add node to parser list
    AddLeafNode(leafNode);
    return leafNode;
}

/*
 * BranchSet --> Branch | BranchSet "," Branch
 */
std::vector<DirectedEdge*> NewickParser::ParseBranchSet(std::string string) {
    //std::cout << "ENTER branch set < " << string << " >"  << std::endl;

    std::vector<DirectedEdge*> branchSet;

    bool lookingForLastBranch = true;
    //separate string into comma separated branches
    while (lookingForLastBranch) {
        std::string nextBranchString = "";

        //check if the remains of the branch set is only the empty string
        if (string.empty()) {
            //nextBranchString remains empty
            lookingForLastBranch = false;
        }
        else if (string[0] == '(') { //next branch is an internal node
            int indexOfLastRightParen = findMatchingRightParenthesis(string);
            if (indexOfLastRightParen == -1) {
                std::cerr << "NewickParser ERROR: Could not find end last parenthesis of internal node." << std::endl;
                exit(EXIT_FAILURE);
            }

            //find first comma after the right parenthesis
            std::string::size_type indexOfFirstCommaAfterParen = string.find_first_of(',', indexOfLastRightParen);

            if (indexOfFirstCommaAfterParen != std::string::npos) { //there is more to come
                //make a branch from beginning until comma
                nextBranchString = string.substr(0,indexOfFirstCommaAfterParen);
                //continue working on everything after the comma
                string = string.substr(indexOfFirstCommaAfterParen + 1);
            }
            else { //this is the last branch
                nextBranchString = string;
                lookingForLastBranch = false;
            }

        }
        else { //next branch is a leaf
            std::string::size_type indexOfFirstComma = string.find_first_of(',');

            if (indexOfFirstComma != std::string::npos) { //there is more to come
            
                //make a branch from beginning until comma
                nextBranchString = string.substr(0,indexOfFirstComma);
                //continue working on everything after the comma
                string = string.substr(indexOfFirstComma + 1);


            }
            else { //this is the last branch
                nextBranchString = string;
                lookingForLastBranch = false;
            }
        }
        DirectedEdge* branch = ParseBranch(nextBranchString);
        branchSet.push_back(branch);
    }
    

    //make a branch from remaining string
    //DirectedEdge* lastBranch = ParseBranch(string);
    //branchSet.push_back(lastBranch);

    return branchSet;
}

/*
 * Branch --> Subtree Length
 */
DirectedEdge* NewickParser::ParseBranch(std::string string) {
    //std::cout << "ENTER branch < " << string << " >"  << std::endl;

    Node* node;
    if (string.empty()) { //no reason to look for length
        node = ParseSubtree(string);
    }
    else {

        //look for colon from this index (0 if leaf)
        int indexToStartColonSearch = 0;

        if (string[0] == '(') { //this branch is an internal node
            int indexOfLastRightParen = findMatchingRightParenthesis(string);
            if (indexOfLastRightParen == -1) {
                std::cerr << "NewickParser ERROR: Could not find end last parenthesis of internal node." << std::endl;
                exit(EXIT_FAILURE);
            }
            indexToStartColonSearch = indexOfLastRightParen;
        }

        //look for colon which indicates length specification
        std::string::size_type indexOfFirstColon = string.find_first_of(':', indexToStartColonSearch);


        if (indexOfFirstColon == std::string::npos) //no length
            node = ParseSubtree(string);
        else {
            node = ParseSubtree(string.substr(0, indexOfFirstColon));
            //TODO use length for something
            //length is everything after the colon. not used at the moment
            std::string length = string.substr(indexOfFirstColon + 1);
        }
    }

    DirectedEdge* branch = new DirectedEdge(GetNextEdgeId());
    //add edge to parser list
    AddDirectedEdge(branch);
    branch->SetToNode(node);
    return branch;
}

