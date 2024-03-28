#ifndef NODE_H
#define NODE_H
#include <vector>
#include <string>

/*
 * The Node is an abstract class that serves as super class for LeafNode and InternalNode
 * - dictates the need for the function AddEdge in all sub classes
 */

//forward declaration
class DirectedEdge;

class Node {
public:

    Node(const std::string &label)
        : label(label)
    {}

    virtual ~Node()
    {}

    void SetLabel(const std::string &label) { this->label = label; }
    const std::string &GetLabel() const     { return label; }

    virtual void AddEdge(DirectedEdge* edge) = 0;

    virtual bool isLeaf() const = 0;
    bool isInternal() const { return !isLeaf(); }

private:
    std::string label;
    int nodeId;

};

#endif
