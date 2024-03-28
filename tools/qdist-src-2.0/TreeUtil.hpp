#ifndef TREE_UTIL_H
#define TREE_UTIL_H

#include "Tree.hpp"
#include "InternalNode.hpp"
#include "LeafNode.hpp"
#include "DirectedEdge.hpp"

/*
 * Various utility routines that work on trees
 */

class Node;
class Path;
class Center;

class TreeUtil {
public:
    TreeUtil();
    ~TreeUtil();

    static void PrintTree(Tree* tree);
    static void PrintSubtree(Node* node, Node* fromNode, int indent);
    static void CheckTree(Tree* tree);
    static void CheckSubtree(Node* node, Node* fromNode);
    static void RenumberTreeAccordingToOther(Tree* tree, Tree* other);

    static std::vector<int> SubtreeLeafSetSizes(Tree* tree);
    static std::vector< std::vector<int> > CalcSharedLeafSetSizes(Tree* t1, Tree* t2);

    static Path* FindPath(LeafNode* fromNode, LeafNode* toNode);
    static Center FindCenter(Tree* tree, LeafNode* a, LeafNode* b, LeafNode* c);
    static std::vector<Center> MakeCenterArrayFromPath(Tree* tree, Path* path);

    static std::vector<LeafNode*> CollectLeavesInSubtree(DirectedEdge* subtreeEdge);

private:
    static int CountLeavesDownwards(Node* node, Node* fromNode, std::vector<int>* subtreeLeafSetSizes);
    static void CalcLeavesUpwards(Tree* tree, std::vector<int>* subtreeLeafSetSizes);

    static void CalcSharedLeafSetSizesDownDown(DirectedEdge* e1, DirectedEdge* e2, std::vector< std::vector<int> >* sharedLeafSetSizes);

    static bool FindPathRecursive(DirectedEdge* edge, LeafNode* endNode, std::vector<DirectedEdge*>*);

    static std::vector<DirectedEdge*> CollectEdgesPointingAwayFromRoot(Tree* tree);
    static void CollectEdgesRecursive(Node* node, Node* fromNode, std::vector<DirectedEdge*>* downEdges);

    static void CollectLeavesRecursive(DirectedEdge* subtreeEdge, std::vector<LeafNode*>* leaves);
};

/*
 * A class for representing a path between two nodes
 * includes:
 *  origin node
 *  end node
 *  a list of all directed edges on the path
 */
class Path {
public:
    Path(Node* fromNode, Node* toNode, const std::vector<DirectedEdge*> &edges) {
        this->fromNode = fromNode;
        this->toNode = toNode;
        this->edges = edges;
    }
    ~Path() {
        this->fromNode = 0;
        this->toNode = 0;
    }

    Node* GetFromNode() {return fromNode;}
    Node* GetToNode() {return toNode;}
    const std::vector<DirectedEdge*> &GetEdges() {return edges;}

private:
    Node* fromNode;
    Node* toNode;
    std::vector<DirectedEdge*> edges;
};

/*
 * A class for representing a center of a triplet (a,b,c)
 * includes:
 *  the center node
 *  three edges pointing from the center towards a,b and c respectively
 */
class Center {
public:
    Center()
        : center(NULL), aEdge(NULL), bEdge(NULL), cEdge(NULL)
    {}

    Center(Node* center, DirectedEdge* aEdge, DirectedEdge* bEdge, DirectedEdge* cEdge)
        : center(center), aEdge(aEdge), bEdge(bEdge), cEdge(cEdge)
    {}

    Node* GetCenterNode() {return center;}
    DirectedEdge* GetAEdge() {return aEdge;}
    DirectedEdge* GetBEdge() {return bEdge;}
    DirectedEdge* GetCEdge() {return cEdge;}

    bool isNull() const { return center == NULL; }

private:
    Node* center;
    DirectedEdge* aEdge;
    DirectedEdge* bEdge;
    DirectedEdge* cEdge;
};


#endif
