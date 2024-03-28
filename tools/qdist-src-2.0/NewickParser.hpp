#ifndef NEWICK_PARSER_H
#define NEWICK_PARSER_H

#include <string>
#include <vector>

#include "Tree.hpp"
#include "InternalNode.hpp"
#include "LeafNode.hpp"
#include "DirectedEdge.hpp"

/*
 * A parser for the Newick format
 * url: http://en.wikipedia.org/wiki/Newick_format
 */

class Node;

class NewickParser {
public:
    NewickParser();
    ~NewickParser();

    Tree* Parse(std::string string);

private:
    int internalIdCount;
    int leafIdCount;
    int edgeIdCount;
    std::vector<InternalNode*> internalNodes;
    std::vector<LeafNode*> leafNodes;
    std::vector<DirectedEdge*> directedEdges;

    void ResetParser() {
        internalIdCount = 0;
        leafIdCount = 0;
        edgeIdCount = 0;
        internalNodes = std::vector<InternalNode*>();
        leafNodes = std::vector<LeafNode*>();
        directedEdges = std::vector<DirectedEdge*>();
    }
    int GetNextInternalId() { return internalIdCount++; }
    int GetNextLeafId()     { return leafIdCount++; }
    int GetNextEdgeId()     { return edgeIdCount++; }
    
    void AddInternalNode(InternalNode* internalNode) { internalNodes.push_back(internalNode); }
    void AddLeafNode(LeafNode* leafNode)             { leafNodes.push_back(leafNode); }
    void AddDirectedEdge(DirectedEdge* directedEdge) { directedEdges.push_back(directedEdge); }
    
    std::vector<InternalNode*> GetInternalNodeList() { return internalNodes; }
    std::vector<LeafNode*> GetLeafNodeList()         { return leafNodes; }
    std::vector<DirectedEdge*> GetDirectedEdgeList() { return directedEdges; }    

    //parsing functions
    Node* ParseSubtree(std::string string);
    Node* ParseInternalNode(std::string string);
    Node* ParseLeafNode(std::string string);
    std::vector<DirectedEdge*> ParseBranchSet(std::string string);
    DirectedEdge* ParseBranch(std::string string);

};

#endif
