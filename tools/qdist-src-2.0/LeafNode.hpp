#ifndef LEAF_NODE_H
#define LEAF_NODE_H

#include "Node.hpp"
#include <string>

/*
 * The LeafNode includes
 *  - label
 *  - id
 *  - a single directed edge pointing to an internal node
 */

class LeafNode : public Node {
public: 
    LeafNode(const std::string &label, int leafId)
        : Node(label),
          leafId(leafId),
          edge(NULL)
    {}

    ~LeafNode()
    {}

    void SetLeafId(int leafId)       { this->leafId = leafId; }
    int GetLeafId()                  { return leafId; }

    void AddEdge(DirectedEdge* edge) { this->edge = edge; }
    DirectedEdge* GetEdge()          { return edge; }

    bool isLeaf() const { return true; }

private:
    int leafId;
    DirectedEdge* edge;
};

#endif
