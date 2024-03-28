#ifndef TREE_H
#define TREE_H

#include <vector>
#include <cstddef>

class Node;
class InternalNode;
class LeafNode;
class DirectedEdge;

/*
 * The Tree class is a container for all the parts of a tree, including:
 *  - an entry point, called root
 *  - a list containing all internal nodes
 *  - a list containing all leaf nodes
 *  - a list containing all edges
 */

//INVARIANT: for each InternalNode/LeafNode/DirectedEdge in the three lists, list index
//           and internalId/leafId/edgeId is the same.

class Tree {
public:
    Tree()
        : root(NULL),
          internalNodes(),
          leafNodes(),
          edges()
    {}

    ~Tree()
    {}

    void SetRoot(Node* root) { this->root = root; }
    Node* GetRoot()          { return root; }

    InternalNode* GetInternalNode(int internalId)
    { return internalNodes[internalId]; }
    LeafNode* GetLeafNode(int leafId)
    { return leafNodes[leafId]; }
    DirectedEdge* GetEdge(int edgeId)
    { return edges[edgeId]; }

    int NumInternalNodes() const { return internalNodes.size(); }
    int NumLeafNodes()     const { return leafNodes.size(); }
    int NumEdges()         const { return edges.size(); }

    void SetInternalNodeList(const std::vector<InternalNode*> &internalNodes)
    { this->internalNodes = internalNodes; }
    void SetLeafNodeList(const std::vector<LeafNode*> &leafNodes)
    { this->leafNodes = leafNodes; }
    void SetEdgeList(const std::vector<DirectedEdge*> &edges)
    { this->edges = edges; }

    const std::vector<InternalNode*> &GetInternalNodes()
    { return internalNodes; }
    const std::vector<LeafNode*> &GetLeafNodes()
    { return leafNodes; }
    const std::vector<DirectedEdge*> &GetEdges()
    { return edges; }

private:
    Node* root;
    std::vector<InternalNode*> internalNodes;
    std::vector<LeafNode*> leafNodes;
    std::vector<DirectedEdge*> edges;
};

#endif
