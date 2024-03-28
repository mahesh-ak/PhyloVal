#ifndef DIRECTED_EDGE_H
#define DIRECTED_EDGE_H

/*
 * The DirectedEdge includes
 *  - id
 *  - an origin node
 *  - an end node
 *  - a pointer to the oppositely directed edge
 */

//forward declaration
class Node;

class DirectedEdge {
public: 
    DirectedEdge(int edgeId)
        : edgeId(edgeId),
          fromNode(NULL),
          toNode(NULL),
          backEdge(NULL)
    {}

    ~DirectedEdge()
    {}

    void SetEdgeId(int edgeId)              { this->edgeId = edgeId; }
    int GetEdgeId()                         { return edgeId; }
    void SetToNode(Node* toNode)            { this->toNode = toNode; }
    Node* GetToNode()                       { return toNode; }
    const Node* GetToNode()           const { return toNode; }
    void SetFromNode(Node* fromNode)        { this->fromNode = fromNode; }
    Node* GetFromNode()                     { return fromNode; }
    const Node* GetFromNode()         const { return fromNode; }
    void SetBackEdge(DirectedEdge* backEdge)
                                            { this->backEdge = backEdge; }
    DirectedEdge* GetBackEdge()             { return backEdge; }
    const DirectedEdge* GetBackEdge() const { return backEdge; }

private:
    int edgeId;
    Node* fromNode;
    Node* toNode;

    DirectedEdge* backEdge;
};

#endif
