#include "TreeUtil.hpp"
#include <iostream>
#include <string>
#include <assert.h>
#include <utility>
#include <algorithm>

TreeUtil::TreeUtil() {
}

TreeUtil::~TreeUtil() {
}

/*
 * A function for pretty-printing a tree
 */
void TreeUtil::PrintTree(Tree* tree) {
    std::cout << "----------------------------------------------------" << std::endl;
    std::cout << "Printing tree with ";
    std::cout << tree->NumInternalNodes() << " inner nodes";
    std::cout << ", " << tree->NumLeafNodes() << " leaf nodes";
    std::cout << " and " << tree->NumEdges() << " edges" << std::endl;

    TreeUtil::PrintSubtree(tree->GetRoot(), NULL, 0);

    std::cout << "----------------------------------------------------" << std::endl;
}

// helper function - prints a (sub)tree
void TreeUtil::PrintSubtree(Node* node, Node* fromNode=NULL, int indent=0) {
    for (int i = 0; i < indent; i++)
        std::cout << " ";

    std::cout << node->GetLabel();

    if (node->isLeaf()) {
        LeafNode* leaf = (LeafNode*)node;
        std::cout << "leaf_id=" << leaf->GetLeafId();
    }

    std::cout << std::endl;

    if (node->isInternal()) {
        InternalNode* internal = (InternalNode*)node;
        std::string::size_type i;
        for (i = 0; i < internal->GetEdges().size(); i++) {
            Node* neighbor = internal->GetEdges()[i]->GetToNode();
            if (neighbor != fromNode)
                TreeUtil::PrintSubtree(neighbor, internal, indent+1);
        }
    }

}

/*
 * A function verifying that a tree is properly connected
 * - only used for debugging purposes
 */
void TreeUtil::CheckTree(Tree* tree) {
    //std::cout << "CHECKING INNER ID's..." << std::endl;
    for (int i = 0; i < tree->NumInternalNodes(); i++) {
        assert(i == tree->GetInternalNode(i)->GetInternalId());
    }

    //std::cout << "CHECKING LEAF ID's..." << std::endl;
    for (int i = 0; i < tree->NumLeafNodes(); i++) {
        assert(i == tree->GetLeafNode(i)->GetLeafId());
    }

    //std::cout << "CHECKING EDGE ID's..." << std::endl;
    for (int i = 0; i < tree->NumEdges(); i++) {
        assert(i == tree->GetEdge(i)->GetEdgeId());
    }

    //check that references in the tree are set
    TreeUtil::CheckSubtree(tree->GetRoot(), NULL);
}

// helper function for the CheckTree function - checks a subtree
void TreeUtil::CheckSubtree(Node* node, Node* fromNode=NULL) {
    if (node->isInternal()) {
        InternalNode* internal = (InternalNode*)node;
        std::string::size_type i;
        for (i = 0; i < internal->GetEdges().size(); i++) {
            DirectedEdge* edge = internal->GetEdges()[i];
            Node* thisNode = edge->GetFromNode();
            Node* neighbor = edge->GetToNode();

            //check stuff
            assert(thisNode != NULL);
            assert(thisNode == internal);
            assert(neighbor != NULL);

            //continue checking
            if (neighbor != fromNode)
                TreeUtil::CheckSubtree(neighbor, internal);
        }
        //std::cout << "CHECKED InternalNode: " << internal->GetLabel() << std::endl;
    }
    else if (node->isLeaf()) {
        LeafNode* leaf = (LeafNode*)node;
        DirectedEdge* edge = leaf->GetEdge();

        //check stuff
        assert(edge != NULL);

        Node* thisNode = edge->GetFromNode();
        Node* neighbor = edge->GetToNode();

        //check stuff
        assert(thisNode != NULL);
        assert(thisNode == leaf);
        assert(neighbor != NULL);
        assert(neighbor == fromNode);

        //std::cout << "CHECKED LeafNode: " << leaf->GetLabel() << std::endl;
    }
}


/*
 * Renumber the leaves in one tree such that both trees have the same leaf-label-leaf-id correspondance
 */
void TreeUtil::RenumberTreeAccordingToOther(Tree* tree, Tree* other) {
    //make a new vector to store the leaf nodes
    std::vector<LeafNode*> newOrderLeaves(tree->NumLeafNodes());

    //renumber labels in 'tree' so they correspond to numbers from 'other'
    const std::vector<LeafNode*> &leaves = tree->GetLeafNodes();
    const std::vector<LeafNode*> &otherLeaves = other->GetLeafNodes();
    for (int i = 0; i < tree->NumLeafNodes(); i++) {
        LeafNode* leaf = leaves[i];
        for (int j = 0; j < other->NumLeafNodes(); j++) {
            LeafNode* otherLeaf = otherLeaves[j];
            if (leaf->GetLabel() == otherLeaf->GetLabel()) {
                leaf->SetLeafId(otherLeaf->GetLeafId());
                newOrderLeaves[leaf->GetLeafId()] = leaf;
            }
        }
    }

    tree->SetLeafNodeList(newOrderLeaves);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Subtree Leaf Set Size
////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * Calculate for each directed edge, the number of leaves in the subtree
 * identified by the edge.
 */
std::vector<int> TreeUtil::SubtreeLeafSetSizes(Tree* tree) {
    std::vector<int> subtreeLeafSetSizes(tree->NumEdges(), -1);

    TreeUtil::CountLeavesDownwards(tree->GetRoot(), NULL, &subtreeLeafSetSizes);
    TreeUtil::CalcLeavesUpwards(tree, &subtreeLeafSetSizes);

    return subtreeLeafSetSizes;
}

/*
 * Helper function for SubtreeLeafSetSizes.
 * Recursively count the number of leaves in each subtree denoted by an
 * edge pointing away from the root.
 */
int TreeUtil::CountLeavesDownwards(Node* node, Node* fromNode, std::vector<int>* subtreeLeafSetSizes) {
    int count = 0;
    if (node->isInternal()) {
        InternalNode* internal = (InternalNode*)node;

        const std::vector<DirectedEdge*> &edges = internal->GetEdges();
        for (std::vector<DirectedEdge*>::size_type i = 0; i < edges.size(); i++) {
            DirectedEdge* edge = edges[i];
            Node* toNode = edge->GetToNode();

            if (toNode != fromNode) {
                int tmpCount = TreeUtil::CountLeavesDownwards(toNode, internal, subtreeLeafSetSizes);
                //store the count for the subtree identified by edge
                (*subtreeLeafSetSizes)[edge->GetEdgeId()] = tmpCount;
                //add to the overall count
                count += tmpCount;
            }

            //clean up
            edge = 0;
            toNode = 0;
        }

        internal = 0;
    }
    else if (node->isLeaf()) {
        count = 1;
    }
    
    return count;
}

/*
 * Helper function for SubtreeLeafSetSizes.
 * Calculate the number of leaves in subtrees pointing towards the root.
 */
void TreeUtil::CalcLeavesUpwards(Tree* tree, std::vector<int>* subtreeLeafSetSizes) {
    for (int i = 0; i < tree->NumEdges(); i++) {
        if ((*subtreeLeafSetSizes)[i] == -1) {
            int backEdgeId = tree->GetEdge(i)->GetBackEdge()->GetEdgeId();
            (*subtreeLeafSetSizes)[i] = tree->NumLeafNodes() - (*subtreeLeafSetSizes)[backEdgeId];
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Shared Leaf Set Size
////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * Calculate, for each pair of directed edges, the number of leaves common to the two subtrees
 * identified by the two edges.
 */
std::vector< std::vector<int> > TreeUtil::CalcSharedLeafSetSizes(Tree* t1, Tree* t2) {
    //calculate the sizes of each subtree in the two trees
    std::vector<int> t1LeafSetSizes = TreeUtil::SubtreeLeafSetSizes(t1);
    std::vector<int> t2LeafSetSizes = TreeUtil::SubtreeLeafSetSizes(t2);

    //collect all edges pointing away from the root in both trees
    std::vector<DirectedEdge*> t1DownEdges = TreeUtil::CollectEdgesPointingAwayFromRoot(t1);
    std::vector<DirectedEdge*> t2DownEdges = TreeUtil::CollectEdgesPointingAwayFromRoot(t2);

    //a matrix to hold the result. initially filled with -1 in each entry
    std::vector< std::vector<int> > sharedLeafSetSizes(t1->NumEdges(), std::vector<int>(t2->NumEdges(), -1));


    //calculate shared leaf set sizes for all down-down pairs of edges
    for (std::vector<DirectedEdge*>::size_type i = 0; i < t1DownEdges.size(); i++) {
        for (std::vector<DirectedEdge*>::size_type j = 0; j < t2DownEdges.size(); j++) {
            TreeUtil::CalcSharedLeafSetSizesDownDown(t1DownEdges[i], t2DownEdges[j], &sharedLeafSetSizes);
        }
    }

    //calculate shared leaf set sizes for remaining pairs of edges
    for (std::vector<DirectedEdge*>::size_type i = 0; i < t1DownEdges.size(); i++) {
        for (std::vector<DirectedEdge*>::size_type j = 0; j < t2DownEdges.size(); j++) {
            int t1DownEdgeId = t1DownEdges[i]->GetEdgeId();
            int t1UpEdgeId = t1DownEdges[i]->GetBackEdge()->GetEdgeId();
            int t2DownEdgeId = t2DownEdges[j]->GetEdgeId();
            int t2UpEdgeId = t2DownEdges[j]->GetBackEdge()->GetEdgeId();

            //up-down
            sharedLeafSetSizes[t1UpEdgeId][t2DownEdgeId] = t2LeafSetSizes[t2DownEdgeId] - sharedLeafSetSizes[t1DownEdgeId][t2DownEdgeId];
            //down-up
            sharedLeafSetSizes[t1DownEdgeId][t2UpEdgeId] = t1LeafSetSizes[t1DownEdgeId] - sharedLeafSetSizes[t1DownEdgeId][t2DownEdgeId];
            //up-up
            sharedLeafSetSizes[t1UpEdgeId][t2UpEdgeId] = t1->NumLeafNodes() - (t1LeafSetSizes[t1DownEdgeId] + t2LeafSetSizes[t2DownEdgeId] - sharedLeafSetSizes[t1DownEdgeId][t2DownEdgeId]);

        }
    }

    return sharedLeafSetSizes;
}

/*
 * Helper function for CalcSharedLeafSetSizes.
 * Recursively calculate the shared leaf set sizes for all pairs of subtrees in the two subtrees given
 */
void TreeUtil::CalcSharedLeafSetSizesDownDown(DirectedEdge* e1, DirectedEdge* e2, std::vector< std::vector<int> >* sharedLeafSetSizes) {

    Node* n1 = e1->GetToNode();
    Node* n2 = e2->GetToNode();

    int e1_id = e1->GetEdgeId();
    int e2_id = e2->GetEdgeId();

    //shared leaf set size not yet calculated
    if ((*sharedLeafSetSizes)[e1_id][e2_id] == -1) {

        //first is a leaf
        if (n1->isLeaf()) {
            LeafNode* leaf1 = (LeafNode*)n1;

            //second is also leaf
            if (n2->isLeaf()) {
                LeafNode* leaf2 = (LeafNode*)n2;

                //compare them
                if (leaf1->GetLeafId() == leaf2->GetLeafId())
                    (*sharedLeafSetSizes)[e1_id][e2_id] = 1; //same leaves
                else
                    (*sharedLeafSetSizes)[e1_id][e2_id] = 0; //diff leaves
            }
            //second is internal. recurse
            else {
                InternalNode* internal2 = (InternalNode*)n2;
                int sum = 0;
                for (std::vector<DirectedEdge*>::size_type i = 0; i < internal2->GetEdges().size(); i++) {
                    DirectedEdge* e2sub = internal2->GetEdges()[i];
                    if (e2sub != e2->GetBackEdge()) { //ensure that the edge points downwards
                        CalcSharedLeafSetSizesDownDown(e1, e2sub, sharedLeafSetSizes);
                        sum += (*sharedLeafSetSizes)[e1_id][e2sub->GetEdgeId()];
                    }
                }
                (*sharedLeafSetSizes)[e1_id][e2_id] = sum;
            }
        }
        //first is internal. recurse
        else {
            InternalNode* internal1 = (InternalNode*)n1;
            int sum = 0;
            for (std::vector<DirectedEdge*>::size_type i = 0; i < internal1->GetEdges().size(); i++) {
                DirectedEdge* e1sub = internal1->GetEdges()[i];
                if (e1sub != e1->GetBackEdge()) { //ensure that the edge points downwards
                    CalcSharedLeafSetSizesDownDown(e1sub, e2, sharedLeafSetSizes);
                    sum += (*sharedLeafSetSizes)[e1sub->GetEdgeId()][e2_id];
                }
            }
            (*sharedLeafSetSizes)[e1_id][e2_id] = sum;
        }
    }

}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Finding paths
////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * Find a path between two leaves by a full traversal of the tree
 */
Path* TreeUtil::FindPath(LeafNode* fromNode, LeafNode* toNode) {
    std::vector<DirectedEdge*>* edges = new std::vector<DirectedEdge*>();
    edges->reserve(100);
    bool success = TreeUtil::FindPathRecursive(fromNode->GetEdge(), toNode, edges);
    assert(success);
    std::reverse(edges->begin(), edges->end());

    //create the path
    Path* path = new Path(fromNode, toNode, *edges);
    return path;
}

/*
 * Helper function for the FindPath function
 */
bool TreeUtil::FindPathRecursive(DirectedEdge* currentEdge, LeafNode* endNode, std::vector<DirectedEdge*>* path) {
    Node* toNode = currentEdge->GetToNode();
    if (toNode->isLeaf()) {
        LeafNode* leaf = (LeafNode*) toNode;
        if (leaf == endNode) {
            path->push_back(currentEdge);
            return true;
        }
        else
            return false;
    }    
    else {
        InternalNode* internal = (InternalNode*) toNode;
        const std::vector<DirectedEdge*> &edges = internal->GetEdges();
        for (unsigned i = 0; i < edges.size(); i++) {
            DirectedEdge* tmpEdge = edges[i];
            //don't go backwards
            if (tmpEdge == currentEdge->GetBackEdge())
                continue;
            bool success = TreeUtil::FindPathRecursive(tmpEdge, endNode, path);
            if (success) {
                path->push_back(currentEdge);
                return true;
            }
        }
        return false;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Finding centers
////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * Find the center of the triplet (a,b,c)
 */
Center TreeUtil::FindCenter(Tree* tree, LeafNode* a, LeafNode* b, LeafNode* c) {
    //find each path between pairs of leaves
    Path* pathAB = TreeUtil::FindPath(a, b);
    Path* pathBC = TreeUtil::FindPath(b, c);
    Path* pathCA = TreeUtil::FindPath(c, a);

    // find the center from the three paths
    std::vector<int> internalsFound = std::vector<int>(tree->NumInternalNodes());

    int centerId = -1;

    for (unsigned i = 0; i < pathAB->GetEdges().size(); i++) {
        Node* node = pathAB->GetEdges()[i]->GetToNode();
        if (node->isInternal()) {
            InternalNode* internal = (InternalNode*)node;
            internalsFound[internal->GetInternalId()] += 1;

            //clean up
            internal = 0;
        }

        //clean up
        node = 0;
    }

    for (unsigned i = 0; i < pathBC->GetEdges().size(); i++) {
        Node* node = pathBC->GetEdges()[i]->GetToNode();
        if (node->isInternal()) {
            InternalNode* internal = (InternalNode*)node;
            internalsFound[internal->GetInternalId()] += 1;

            //clean up
            internal = 0;
        }

        //clean up
        node = 0;
    }

    for (unsigned i = 0; i < pathCA->GetEdges().size(); i++) {
        Node* node = pathCA->GetEdges()[i]->GetToNode();
        if (node->isInternal()) {
            InternalNode* internal = (InternalNode*)node;
            internalsFound[internal->GetInternalId()] += 1;

            //check if the third occurrence was found
            if (internalsFound[internal->GetInternalId()] == 3) {
                centerId = internal->GetInternalId();
                break;
            }

            //clean up
            internal = 0;
        }

        //clean up
        node = 0;
    }

    assert(centerId != -1);

    Node* centerNode = tree->GetInternalNode(centerId);

    DirectedEdge* aEdge = NULL;
    for (unsigned i = 0; i < pathAB->GetEdges().size(); i++) {
        DirectedEdge* edge = pathAB->GetEdges()[i];
        if (edge->GetToNode() == centerNode) {
            aEdge = edge->GetBackEdge();
            break;
        }
        //clean up
        edge = 0;
    }

    DirectedEdge* bEdge = NULL;
    for (unsigned i = 0; i < pathBC->GetEdges().size(); i++) {
        DirectedEdge* edge = pathBC->GetEdges()[i];
        if (edge->GetToNode() == centerNode) {
            bEdge = edge->GetBackEdge();
            break;
        }
        //clean up
        edge = 0;
    }

    DirectedEdge* cEdge = NULL;
    for (unsigned i = 0; i < pathCA->GetEdges().size(); i++) {
        DirectedEdge* edge = pathCA->GetEdges()[i];
        if (edge->GetToNode() == centerNode) {
            cEdge = edge->GetBackEdge();
            break;
        }
        //clean up
        edge = 0;
    }

    //clean up paths
    delete pathAB;
    delete pathBC;
    delete pathCA;

    assert(aEdge != NULL);
    assert(bEdge != NULL);
    assert(cEdge != NULL);

    return Center(centerNode, aEdge, bEdge, cEdge);
}

/*
 * Construct an array holding every center found on the given path
 */
std::vector<Center> TreeUtil::MakeCenterArrayFromPath(Tree* tree, Path* path) {
    std::vector<Center> centers(tree->NumLeafNodes());

    //traverse the path
    for (unsigned i = 0; i < path->GetEdges().size(); i++) { //don't consider last edge as it points to end leaf
        Node* node = path->GetEdges()[i]->GetToNode();

        //for all internal nodes found
        if (node->isInternal()) {
            //the center node
            InternalNode* internal = (InternalNode*) node;
            //Ta path
            DirectedEdge* backEdge = path->GetEdges()[i]->GetBackEdge();
            //Tb path
            DirectedEdge* forwardEdge = path->GetEdges()[i+1];
            
            //for all subtrees of the center node
            for (unsigned j = 0; j < internal->GetEdges().size(); j++) {
                DirectedEdge* edge = internal->GetEdges()[j];
                //not including Ta and Tb
                if (edge != backEdge && edge != forwardEdge) {
                    //find all leaves in the subtree
                    std::vector<LeafNode*> leaves = TreeUtil::CollectLeavesInSubtree(edge);
                    
                    //construct centers for a, b and every leaf in the subtree
                    for (unsigned k = 0; k < leaves.size(); k++) {
                        LeafNode* leaf = leaves[k];
                        centers[leaf->GetLeafId()] = Center(internal, backEdge, forwardEdge, edge);
                    }
                }
            }
        }        
    }
    
    return centers;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// General tree utility routines
////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * Collect a list of all leaves in the subtree identified by a directed edge
 */
std::vector<LeafNode*> TreeUtil::CollectLeavesInSubtree(DirectedEdge* subtreeEdge) {
    std::vector<LeafNode*> leaves;
    leaves.reserve(100);
    TreeUtil::CollectLeavesRecursive(subtreeEdge, &leaves);
    return leaves;
}

/*
 * Helper routine for CollectLeavesInSubtree
 */
void TreeUtil::CollectLeavesRecursive(DirectedEdge* subtreeEdge, std::vector<LeafNode*>* leaves) {
    Node* toNode = subtreeEdge->GetToNode();
    if (toNode->isLeaf()) {
        leaves->push_back((LeafNode*) toNode);
    }    
    else {
        InternalNode* internal = (InternalNode*) toNode;
        const std::vector<DirectedEdge*> &edges = internal->GetEdges();
        for (unsigned i = 0; i < edges.size(); i++) {
            DirectedEdge* tmpEdge = edges[i];
            //don't go backwards
            if (tmpEdge == subtreeEdge->GetBackEdge())
                continue;
            else
                TreeUtil::CollectLeavesRecursive(tmpEdge, leaves);
        }
    }
}

/*
 * Collect all directed edges pointing downwards from the root in the given tree
 */
std::vector<DirectedEdge*> TreeUtil::CollectEdgesPointingAwayFromRoot(Tree* tree) {
    std::vector<DirectedEdge*> downEdges;
    TreeUtil::CollectEdgesRecursive(tree->GetRoot(), NULL, &downEdges);

    return downEdges;
}

/*
 * Helper routine for CollectEdgesPointingAwayFromRoot
 */
void TreeUtil::CollectEdgesRecursive(Node* node, Node* fromNode, std::vector<DirectedEdge*>* downEdges) {
    if (node->isInternal()) {
        InternalNode* internal = (InternalNode*)node;

        const std::vector<DirectedEdge*> &edges = internal->GetEdges();
        for (std::vector<DirectedEdge*>::size_type i = 0; i < edges.size(); i++) {
            DirectedEdge* edge = edges[i];
            Node* toNode = edge->GetToNode();

            if (toNode != fromNode) {
                downEdges->push_back(edge);
                TreeUtil::CollectEdgesRecursive(toNode, internal, downEdges);
            }
        }
    }
}
