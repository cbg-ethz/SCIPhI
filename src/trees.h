/**
 * SCIPhI: Single-cell mutation identification via phylogenetic inference
 * <p>
 * Copyright (C) 2018 ETH Zurich, Jochen Singer
 * <p>
 * This file is part of SCIPhI.
 * <p>
 * SCIPhI is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * <p>
 * SCIPhI is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * <p>
 * You should have received a copy of the GNU General Public License
 * along with SCIPhI. If not, see <http://www.gnu.org/licenses/>.
 *
 * @author: Jochen Singer
 */
#ifndef TREES_H_
#define TREES_H_


#include <stdbool.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <queue>

#include "boost/multi_array.hpp"

#include "sciphi_config.h"
#include "scoreTree.h"
#include "output.h"

using namespace std;

// This function returns the node ID of the sibling of node "vertex"
unsigned getSibling(typename Config<SampleTree>::TGraph const & g, unsigned vertex)
{
    unsigned parentNode = source(*in_edges(vertex, g).first, g);
    if (target(*out_edges(parentNode,  g).first, g) == vertex)
    {
        return target(*(out_edges(parentNode, g).first + 1), g);
    }
    return target(*out_edges(parentNode, g).first, g);
}

//TODO write specialization for sampletree
// This function draws a new sibling node for a given vertex
template <typename TTreeType>
unsigned getNewSibling(Config<TTreeType> const & config, unsigned vertex)
{
    boost::dynamic_bitset<> bitSet(num_vertices(config.getTree()), 0); 
    unsigned numElemInSubtree = 0;

    ExtractNodesBFSVisitor vis(bitSet, numElemInSubtree);
    breadth_first_search(config.getTree(), vertex, visitor(vis));

    bitSet[source(*in_edges(vertex, config.getTree()).first, config.getTree())] = true;
    ++numElemInSubtree;

    unsigned nonDescendantsRank = rand() % (num_vertices(config.getTree()) - 1 - numElemInSubtree);

	std::vector<unsigned> ancestors;
    unsigned i = 0;
	for(unsigned counter = 0; ; ++i)
    {
		if(bitSet[i]==false)
        {
            if (counter == nonDescendantsRank)
            {
                break;
            }
			++counter;
        }
    }

	return i;
}

// This function prunes and then re-attaches nodes
inline
void
pruneAndReAttach(Config<SampleTree> & config)
{
    // pick a node to move with its subtree (artifical root is excluded)
    unsigned nodeToMove = rand() % (num_vertices(config.getTree()) - 1);
    unsigned parentNodeToMove = source(*in_edges(nodeToMove, config.getTree()).first, config.getTree());
    while (parentNodeToMove == num_vertices(config.getTree()) - 1)
    {
        nodeToMove = rand() % (num_vertices(config.getTree()) - 1);
        parentNodeToMove = source(*in_edges(nodeToMove, config.getTree()).first, config.getTree());
    }
    unsigned sibling = target(*out_edges(parentNodeToMove, config.getTree()).first, config.getTree());
    if (sibling == nodeToMove)
        sibling = target(*(++out_edges(parentNodeToMove, config.getTree()).first), config.getTree());
    unsigned grandParentNodeToMove = source(*in_edges(parentNodeToMove, config.getTree()).first, config.getTree());

    add_edge(grandParentNodeToMove, sibling, config.getTree());
    remove_edge(grandParentNodeToMove, parentNodeToMove, config.getTree());
    remove_edge(parentNodeToMove, sibling, config.getTree());
    
    unsigned newSibling = getNewSibling(config, nodeToMove);      

    unsigned newGrandParent = source(*in_edges(newSibling, config.getTree()).first, config.getTree());
    remove_edge(newGrandParent, newSibling, config.getTree());
    add_edge(newGrandParent, parentNodeToMove, config.getTree());
    add_edge(parentNodeToMove, newSibling, config.getTree());
}

// This function swaps the labels of ndoes
inline
void
swapNodeLabels(Config<MutationTree> & config)
{

    unsigned firstNode = rand() % config.getNumMutations() + 1; 
    unsigned secondNode = rand() % config.getNumMutations() + 1; 
    while (firstNode == secondNode)
        secondNode = rand() % config.getNumMutations() + 1; 

    std::swap(config.getTree()[firstNode].mutation, config.getTree()[secondNode].mutation);
}

// This function swaps the labels of ndoes
inline
void
swapNodeLabels(Config<SampleTree> & config)
{

    unsigned firstNode = rand() % config.getNumSamples(); 
    unsigned secondNode = rand() % config.getNumSamples(); 
    while (firstNode == secondNode)
        secondNode = rand() % config.getNumSamples(); 

    unsigned numInternalNodes = config.getNumSamples() - 1;
    // swap the node labels
    firstNode += numInternalNodes;
    secondNode += numInternalNodes;

    std::swap(config.getTree()[firstNode].sample, config.getTree()[secondNode].sample);
}

// This function swaps to subtrees
inline
void
swapSubtrees(Config<MutationTree> & config)
{
    unsigned firstNode = rand() % config.getNumMutations() + 1; // exlcude root
    unsigned secondNode = rand() % config.getNumMutations() + 1; // exlcude root


    while (firstNode == secondNode)
        secondNode = rand() % config.getNumMutations() + 1; // exlcude root

    unsigned parentFirstNode = source(*in_edges(firstNode,config.getTree()).first, config.getTree());
    unsigned parentSecondNode = source(*in_edges(secondNode,config.getTree()).first, config.getTree());
    

    // remove and add connection of fist node
    remove_edge(parentFirstNode, firstNode, config.getTree());
    add_edge(parentSecondNode, firstNode, config.getTree());
    
    // remove and add connection of second node
    remove_edge(parentSecondNode, secondNode, config.getTree());
    add_edge(parentFirstNode, secondNode, config.getTree());
}

inline
void
swapSubtrees(Config<SampleTree> & /*config*/)
{}


inline
void createInitialTree(Config<SampleTree> & config)
{
    unsigned i = 1;
    for (; i < config.getNumSamples() - 1; ++i)
    {
        add_edge((i-1)/2, i, config.getTree());
        config.getTree()[i].sample = -1;
    }
    for (; i < 2 * config.getNumSamples() - 1; ++i)
    {
        add_edge((i-1)/2, i, config.getTree());
        config.getTree()[i].sample = i - (config.getNumSamples() - 1);
    }   

    // add artificial node
    add_edge(2 * config.getNumSamples() - 1, 0, config.getTree());

    config.setTmpTree(config.getTree());
}

inline
void createRandomTree(Config<MutationTree> & config)
{
    std::vector<unsigned> mutToVertices;
    mutToVertices.resize(config.getNumMutations());
    for (unsigned i = 0; i < config.getNumMutations(); ++i)
        mutToVertices[i] = i;

    for (unsigned i = 0; i < config.getNumMutations() * 100; ++i)
        std::swap(mutToVertices[rand() % config.getNumMutations()], mutToVertices[rand() % config.getNumMutations()]);

    for (unsigned i = 1; i < config.getNumMutations() + 1; ++i) // exclude root
    {
        add_edge(rand() % (i), i, config.getTree());
        config.getTree()[i].mutation = mutToVertices[i-1];
    }

    config.setTmpTree(config.getTree());
}

inline
void createInitialTree(Config<MutationTree> & config)
{
    createRandomTree(config);
}

// Forward declaration
template <typename TTreeType>
unsigned getBestAttachmentPosition(Config<TTreeType> & config,
							  unsigned int attachment);
// find the node a mutation attaches to best
inline
void getMutation2NodeAssignment(std::vector<unsigned> & mut2Node,
                                Config<SampleTree> config)
{
    mut2Node.resize(config.getNumMutations());
    for (unsigned i = 0; i < config.getNumMutations(); ++i)
    {
        mut2Node[i] = getBestAttachmentPosition(config, i);
    }
}

// This funtion fills a two dimensional table with the mutation to sample assignment
inline
void getMut2SampleAssignment(boost::multi_array<bool, 2> & mut2Sample,
                            Config<SampleTree> & config,
                            std::vector<unsigned> const & mut2Node)
{
    typedef boost::multi_array<bool, 2> TArray;
    typedef TArray::index           TArrayIndex; 

    mut2Sample.resize(boost::extents[config.getNumSamples()][config.getNumMutations()]);
    for (TArrayIndex i = 0; i < config.getNumSamples(); ++i)
    {
        for (TArrayIndex j = 0; j < config.getNumMutations(); ++j)
        {
            mut2Sample[i][j] = 0;
        }
    }
    
    std::vector<unsigned> sampleNode2SampleID;
    sampleNode2SampleID.resize(config.getNumSamples());
	for (unsigned int i = config.getNumSamples() - 1; i < 2 * config.getNumSamples() - 1; ++i)
    {
        sampleNode2SampleID[i - (config.getNumSamples() - 1)] = config.getTree()[i].sample;
    }

    for (unsigned i = 0; i < config.getNumSamples(); ++i)
    {
        unsigned posInTree = config.getNumSamples() - 1 + i;
        while (posInTree != num_vertices(config.getTree()) - 1)
        {
            for (unsigned j = 0; j < config.getNumMutations(); ++j)
            {
                if (mut2Node[j] == posInTree)
                {
                    mut2Sample[sampleNode2SampleID[i]][j] = true;
                }
            }
            posInTree = source(*in_edges(posInTree, config.getTree()).first, config.getTree());
        }
    }
}

// This function retrieves the mutations assigned to a node
template<typename TTreeType>
std::vector<std::vector<unsigned> >
getMutationsOfNode(Config<TTreeType> & config)
{
    std::vector<std::vector<unsigned> > mutationsOfNodes(num_vertices(config.getTree()));

    typename Config<TTreeType>::TAttachmentScores attachments;
    attachments.resize(num_vertices(config.getTree()) - 1);

    for (unsigned attach = 0; attach < config.getNumMutations(); ++attach)
    {
        unsigned node = getBestAttachmentPosition(config, attach);
        mutationsOfNodes[node].push_back(attach);
    }
    return mutationsOfNodes;
}

// This function simplifies a tree bu removing nodes without mutations assigned
boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, Vertex<SimpleTree>> 
simplifyTree(Config<SampleTree> & config)
{
    typedef boost::adjacency_list<boost::vecS, boost::vecS, 
                            boost::bidirectionalS, Vertex<SimpleTree>> TGraph;
    TGraph newTree;
    std::vector<std::vector<unsigned> > mutationsOfNodes = getMutationsOfNode(config);    
    std::stack<unsigned> processedVertices;

    SimplifyTreeDFSVisitor vis(newTree, mutationsOfNodes, processedVertices);
    unsigned rootVertex = target(*out_edges(num_vertices(config.getTree()) - 1, config.getTree()).first, config.getTree());
    depth_first_search(config.getTree(), visitor(vis).root_vertex(rootVertex));

    return newTree;
}

// This function prints the tree
inline
void printTree(Config<MutationTree>::TGraph const & g, std::string const & fileName)
{
    std::ofstream ofs(fileName.c_str());
    write_graphviz(ofs, g, boost::make_label_writer(boost::get(&Config<MutationTree>::TGraph::vertex_property_type::mutation, g)));
}

// This function prints the mutation to sample assignment to disc
inline
void printMutation2SampleAssignment(Config<SampleTree> & config,
        std::string const & outfileName)
{

    std::vector<unsigned> mut2Node;
    getMutation2NodeAssignment(mut2Node, config);

    boost::multi_array<bool, 2> mut2Sample;
    getMut2SampleAssignment(mut2Sample,
                            config,
                            mut2Node);

    std::ofstream outFile;
    outFile.open(outfileName);
    outFile << "chrom\tposition\t";
	for (unsigned int i = config.getNumSamples() - 1; i < 2 * config.getNumSamples() - 1; ++i)
        outFile << "cell" << i << "\t";
    outFile << "\n";
    for (unsigned int j = 0; j < config.getNumMutations(); ++j)
    {
        outFile << std::get<0>(config.indexToPosition[j]) << "\t" << std::get<1>(config.indexToPosition[j]) << "\t";
        for (unsigned int i = 0; i < config.getNumSamples(); ++i)
        {
			outFile << mut2Sample[i][j] << "\t";
		}
		outFile << std::endl;
	}
    outFile.close();
}

// This class is experimental and computes the number of possibilites a 
// chromosome could have been lost after a mutation occured
class ComputeNumPlacementsDFSVisitor : public boost::default_dfs_visitor
{
    Config<SampleTree> const & config;
    std::vector<unsigned> & numPlacements;

public:

    ComputeNumPlacementsDFSVisitor(Config<SampleTree> const & config_,
            std::vector<unsigned> & numPlacements_) :
        config(config_),
        numPlacements(numPlacements_)
    {}

    template <typename TVertex >
    void discover_vertex(TVertex v, boost::adjacency_list<boost::vecS, 
                                  boost::vecS, 
                                  boost::bidirectionalS, 
                                  Vertex<SampleTree>> const & g) const
    {
        if (v == num_vertices(g) - 1)
        {
            numPlacements[v] = 0;
            return;
        }
        return;

    }
    
    template <typename TVertex >
    void finish_vertex(TVertex v, boost::adjacency_list<boost::vecS, 
                                  boost::vecS, 
                                  boost::bidirectionalS, 
                                  Vertex<SampleTree>> const & g) const
    {
        if (v == num_vertices(g) - 1)
        {
            return;
        }

        if (g[v].sample != -1)
        {
            numPlacements[v] = 0;
            return ;
        }

        auto it = out_edges(v,g).first;
        numPlacements[v] = numPlacements[target(*it, g)] + numPlacements[target(*(it+1), g)];
        if (g[target(*it, g)].sample == -1)
        {
            ++numPlacements[v];
        }
        if (g[target(*(it+1), g)].sample == -1)
        {
            ++numPlacements[v];
        }
        numPlacements.back() += numPlacements[v];
        return ;
    }
};

unsigned getNumPlacements(Config<SampleTree> const & config)
{
    static std::vector<unsigned> numPlacements; // static because the vector should be re-used
    numPlacements.resize(num_vertices(config.getTree()), 0);
    ComputeNumPlacementsDFSVisitor vis(config, numPlacements);
    depth_first_search(config.getTree(), visitor(vis).root_vertex(num_vertices(config.getTree()) - 1));
   
    return numPlacements.back();
}

/*
template <typename TTreeType>
unsigned
getSimpleDistance(TTreeType const & trueTree, 
                  TTreeType const & newTree)
{
    unsigned result = 0;

    std::vector<unsigned> parentsTrueTree;
    parentsTrueTree.resize(num_vertices(trueTree));
    for (unsigned i = 1; i < num_vertices(trueTree); ++i)
        parentsTrueTree[i] = source(*in_edges(i, trueTree).first, trueTree);
    
    std::vector<unsigned> parentsNewTree;
    parentsNewTree.resize(num_vertices(newTree));
    for (unsigned i = 1; i < num_vertices(newTree); ++i)
        parentsNewTree[i] = source(*in_edges(i, newTree).first, newTree);

    for (unsigned i = 1; i < parentsTrueTree.size(); ++i)
        if (parentsTrueTree[i] != parentsNewTree[i])
            ++result;

    return result;
}

template <typename TTreeType>
bool
isDuplicateTreeFast(std::vector<TTreeType> const & trees, TTreeType const & newTree)
{
    for (unsigned i = 0; i < trees.size(); ++i)
        if (getSimpleDistance(trees[i], newTree) != 0)
            return true;

    return false;
}
*/
#endif /* TREES_H_ */
