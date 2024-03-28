#ifndef INTERNAL_NODE_H
#define INTERNAL_NODE_H

#include "Node.hpp"
#include "DirectedEdge.hpp"

#include <string>
#include <vector>

/*
 * The InternalNode includes
 *  - label
 *  - id
 *  - a list containing all directed edges
 */

class InternalNode : public Node {
public: 
    InternalNode(const std::string &label, int internalId)
        : Node(label),
          internalId(internalId),
          edges(), internalEdgesIdxs()
    {}

    ~InternalNode()
    {}

    void SetInternalId(int internalId)
    { this->internalId = internalId; }

    int GetInternalId() { return internalId; }

    void AddEdge(DirectedEdge* edge)
    {
        if(edge->GetToNode()->isInternal())
            internalEdgesIdxs.push_back(edges.size());
        edges.push_back(edge);
    }
    const std::vector<DirectedEdge*> &GetEdges()             const { return edges; }
    const std::vector<unsigned>      &GetInternalEdgesIdxs() const { return internalEdgesIdxs; }

    bool isLeaf() const { return false; }

private:
    int internalId;
    std::vector<DirectedEdge*> edges;
    std::vector<unsigned> internalEdgesIdxs;
};

#endif
