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

#define SEQAN_HAS_ZLIB 1

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <time.h>
#include <map>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/basic.h>
#include <seqan/arg_parse.h>

#include <boost/random/uniform_01.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/negative_binomial.hpp>

using namespace seqan;

struct Options
{
    String<char> outFileMutations;
    String<char> outFileMutationAssignments;
    String<char> outFileMutationAssignmentsWithOutDropOuts;
    String<char> outFileTreeStructure;
    String<char> outFileParameters;
    String<char> pileupPrefix;
	unsigned numSamples;
	unsigned numMutations;
    unsigned numPositions;
	double dropOutRate;
	double missingInformation;
    double clbm; // chromosome loss before mutation happend
    double clam; // chromosome loss after mutation happend
    bool assignMutationsToLeafs;
	unsigned seed;

    Options() :
		numSamples(10),
		numMutations(50),
        numPositions(100000),
		dropOutRate(0.0),
		missingInformation(0.0),
        clbm(0.0),
        clam(0.0),
        assignMutationsToLeafs(true),
		seed(42)
	{};
};

ArgumentParser::ParseResult
parseCommandLine(Options & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("create_single_cell_seq_tree");
    setShortDescription(parser, "Subsample paired-end fastq files");
    setVersion(parser, "0.1");
    setDate(parser, "December 2015");
	
	addOption(parser, ArgParseOption("o", "OutputFilePrefix", "Prefix of output files", ArgParseArgument::OUTPUT_FILE, "OUTPREFIX"));
	addOption(parser, ArgParseOption("ns", "NumSamples", "Number of samples to generate.", ArgParseArgument::INTEGER, "NUMSAMPLES"));
	addOption(parser, ArgParseOption("nm", "NumMutations", "Number of mutations to generate.", ArgParseArgument::INTEGER, "NUMMUTATIONS"));
	addOption(parser, ArgParseOption("np", "NumPositions", "Number of positions to generate.", ArgParseArgument::INTEGER, "NUMMUTATIONS"));
	addOption(parser, ArgParseOption("dor", "DropOutRate", "The allelic drop out rate.", ArgParseArgument::DOUBLE, "DROPOUTRATE"));
	addOption(parser, ArgParseOption("clbm", "ChromLostBeforeMutation", "Fraction of homocygous mutations.", ArgParseArgument::DOUBLE, "cygocitycoeff"));
	addOption(parser, ArgParseOption("clam", "ChromLosAfterMut", "The allelic drop out rate.", ArgParseArgument::DOUBLE, "DROPOUTRATE"));
	addOption(parser, ArgParseOption("mi", "MissingInformation", "Fraction of sites with no sequencing information.", ArgParseArgument::DOUBLE, "MISSINGINFORMATION"));
	addOption(parser, ArgParseOption("s", "Seed", "Seed for random number generator.", ArgParseArgument::INTEGER, "SEED"));
	addOption(parser, ArgParseOption("nml", "NoMutationsInLeafs", "Do not assign mutations to leafes."));
    
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    String<char> outFilePrefix;
    if (isSet(parser, "o"))
        getOptionValue(outFilePrefix, parser, "o");
	else
        outFilePrefix = "tree";

    options.outFileMutations = outFilePrefix;
    append(options.outFileMutations, "_mutations.tsv");
    options.outFileMutationAssignments = outFilePrefix;
    append(options.outFileMutationAssignments, "_mut2sample.tsv");
    options.outFileMutationAssignmentsWithOutDropOuts = outFilePrefix;
    append(options.outFileMutationAssignmentsWithOutDropOuts, "_mut2sample_noDropOuts.tsv");
    options.outFileTreeStructure = outFilePrefix;
    append(options.outFileTreeStructure, "_structure.dot");
    options.pileupPrefix = outFilePrefix;
    options.outFileParameters = outFilePrefix;
    append(options.outFileParameters, "_params.txt");
    
    if (isSet(parser, "s"))
        getOptionValue(options.seed, parser, "s");
	else
		options.seed = time(NULL);
    if (isSet(parser, "ns"))
        getOptionValue(options.numSamples, parser, "ns");
    if (isSet(parser, "nm"))
        getOptionValue(options.numMutations, parser, "nm");
    if (isSet(parser, "np"))
        getOptionValue(options.numPositions, parser, "np");
    if (isSet(parser, "dor"))
        getOptionValue(options.dropOutRate, parser, "dor");
    if (isSet(parser, "clbm"))
        getOptionValue(options.clbm, parser, "clbm");
    if (isSet(parser, "clam"))
        getOptionValue(options.clam, parser, "clam");
    if (isSet(parser, "mi"))
        getOptionValue(options.missingInformation, parser, "mi");
    if (isSet(parser, "nml"))
        options.assignMutationsToLeafs = false;
    return ArgumentParser::PARSE_OK;
}

std::vector<std::pair<unsigned, unsigned>> treeStructureToChildVector(std::vector<unsigned> treeStructure)
{
    std::vector<std::pair<unsigned, unsigned>> childVector;
    childVector.resize(treeStructure.size(), std::pair<unsigned, unsigned>(UINT_MAX, UINT_MAX));

    unsigned parent = UINT_MAX;
    for (size_t i = 1; i < treeStructure.size(); ++i)
    {
        parent = treeStructure[i];
        if (childVector[parent].first == UINT_MAX)
        {
            childVector[parent].first = i;
        }
        else
        {
            childVector[parent].second = i;
        }
    }
    return childVector;
}

void createTreeStructure(std::vector<unsigned> & treeStructure, 
						 unsigned numMutations)
{
	treeStructure.resize(numMutations + 1, 0);		//One additional node for the root
	for (unsigned i = 1; i < treeStructure.size(); ++i)
		treeStructure[i] = rand() % i;
}

void createTreeStructureSampleTree(std::vector<unsigned> & treeStructure, 
						 unsigned numSamples)
{
    std::vector<unsigned> numChildren;
    numChildren.resize(numSamples - 1, 0);
	treeStructure.resize(2 * numSamples - 1, - 1);
    unsigned i = 1;
	for (; i < numSamples - 1; ++i)
    {
        unsigned id = rand() % i;
        while (numChildren[id] >= 2)
            id = rand() % i;

		treeStructure[i] = id;
        ++numChildren[id];
    }

    unsigned posLeaf = 0;
	for (; i < 2 * numSamples - 1; ++i)
    {
        while (numChildren[posLeaf] >= 2)
            ++posLeaf;
		
        treeStructure[i] = posLeaf;
        ++numChildren[posLeaf];
    }
  
    for (unsigned i = 0; i < numChildren.size(); ++i)
        assert(numChildren[i] == 2);
}

void assignSamplesToNodesMutTree(std::vector<unsigned> & sampleToNodeAssignment,
						  std::vector<unsigned> const & treeStructure, 
						  unsigned numSamples)
{
	sampleToNodeAssignment.resize(numSamples);
    unsigned i = 0;
	for (; i < treeStructure.size() - 1; ++i)
		sampleToNodeAssignment[i] = i + 1;
	for (; i < sampleToNodeAssignment.size(); ++i)
		sampleToNodeAssignment[i] = (rand() % (treeStructure.size() - 1)) + 1;
}


//unsigned getSubTreeSize(std::vector<unsigned> const & treeStructure,
//        std::vector<std::pair<unsigned, unsigned>> const & childVector,
//        unsigned node)
//{
//    unsigned size = 0;
//    while(mutationToNodeAssignment[node][mutation] != mutType)
//    {
//        //left node
//        if (childVector[node].first != -1 && mutationToNodeAssignment[childVector[node].first][mutation] == 0)
//        {
//            node = childVector[node].first;
//        }
//        else if (childVector[node].second != -1 && mutationToNodeAssignment[childVector[node].second][mutation] == 0)
//        {
//            node = childVector[node].second;
//        }
//        else
//        {
//            ++size;
//            node = treeStructure[node];
//        }
//    }
//    return size;
//}

struct GetSize
{
    unsigned size;
    unsigned startNode = UINT_MAX;

    GetSize() : size(0) {};
    GetSize(unsigned startNode_) : size(0), startNode(startNode_) {};

    bool operator()() 
    {
        ++this->size;
        return false;
    };
    bool operator()(unsigned currentNode, bool isLeaf)
    {
        if(!isLeaf && currentNode != this -> startNode)
        {
            ++this->size;
        }
        return false;
    };
};

struct GetHomoNodeID
{
    unsigned nodePos = UINT_MAX;
    unsigned rootNode = UINT_MAX;
    unsigned counter = 0;
    unsigned id = UINT_MAX;

    GetHomoNodeID(unsigned counter, unsigned rootNode_) :nodePos(counter), rootNode(rootNode_) {};

    bool operator()(unsigned nodeID, bool isLeaf) 
    {
        if (!isLeaf)
        {
            if (nodePos == counter)
            {
                id = nodeID;
                return true;
            }
        ++this->counter;
        }
        return false;
    };
};


template <typename TType>
void determineHomoNodeId(std::vector<unsigned> const & treeStructure,
        std::vector<std::pair<unsigned, unsigned>> const & childVector,
        unsigned node,
        TType & obj)
{
    std::vector<bool> nodeVisits;
    nodeVisits.resize(treeStructure.size(), 0);
    unsigned currentNode = node;
    unsigned test = 0;
    while(!nodeVisits[node])
    {
        if (childVector[currentNode].first != UINT_MAX && !nodeVisits[childVector[currentNode].first])
        {
            currentNode = childVector[currentNode].first;
        }
        else if (childVector[currentNode].second != UINT_MAX && !nodeVisits[childVector[currentNode].second])
        {
            currentNode = childVector[currentNode].second;
        }
        else
        {
            nodeVisits[currentNode] = 1;
            if (obj(currentNode, childVector[currentNode].first == UINT_MAX))
            {
                break;
            }
            currentNode = treeStructure[currentNode];
            ++test;
        }
    }
}

void passMutToSubTree(std::vector<std::vector<unsigned>> & mutationToNodeAssignment,
        std::vector<unsigned> const & treeStructure,
        std::vector<std::pair<unsigned, unsigned>> const & childVector,
        unsigned node,
        unsigned mutation,
        unsigned mutType)
{
    unsigned currentNode = node;
    std::cout << node << " " << mutation << std::endl;
    while(mutationToNodeAssignment[node][mutation] != mutType)
    {
        //left node
        if (childVector[currentNode].first != UINT_MAX && mutationToNodeAssignment[childVector[currentNode].first][mutation] != mutType)
        {
            currentNode = childVector[currentNode].first;
        }
        else if (childVector[currentNode].second != UINT_MAX && mutationToNodeAssignment[childVector[currentNode].second][mutation] != mutType)
        {
            currentNode = childVector[currentNode].second;
        }
        else
        {
            mutationToNodeAssignment[currentNode][mutation] = mutType;
            currentNode = treeStructure[currentNode];
        }
    }
}

void assignMutation(
        std::vector<std::vector<unsigned>> & mutationToNodeAssignment,
        std::vector<unsigned> const & treeStructure,
        std::vector<std::pair<unsigned, unsigned>> const & childVector,
        unsigned node,
        unsigned mutation,
        std::vector<std::vector<unsigned>> & heteroMutsPerNode,
        //std::vector<unsigned> const & numPlacements,
        //unsigned sumNumPlacemanets,
        Options const & options)
{
    unsigned mutType = 1;
    //std::cout <<  options.clbm << std::endl;
    if (static_cast<double>(rand())/static_cast<double>(RAND_MAX) < options.clbm)
    {
	    mutType = 3;
    }
    else
    {
        heteroMutsPerNode[node].push_back(mutation);
    }
    passMutToSubTree(mutationToNodeAssignment, treeStructure, childVector, node, mutation, mutType);

    //GetSize getSize(node);
    //dfs(treeStructure, childVector, node, getSize);


    
    //double insertHomo = static_cast<double>(rand()) / RAND_MAX ;
    //std::cout << (static_cast<double>(numPlacements[node])/static_cast<double>(sumNumPlacemanets)) << " " << insertHomo << std::endl;
    ////if (mutType == 1 && insertHomo < options.clam * (static_cast<double>(numPlacements[node])/static_cast<double>(sumNumPlacemanets))){
    //if (mutType == 1 && insertHomo < options.clam){ // * (static_cast<double>(numPlacements[node])/static_cast<double>(sumNumPlacemanets))){
    //    std::cout << getSize.size << " " << node << std::endl;
    //    unsigned homoNodeId = rand() % getSize.size;
    //    GetHomoNodeID getHomeNodeId(homoNodeId, node);
    //    dfs(treeStructure, childVector, node, getHomeNodeId);
    //    if (rand() % 2 == 0)
    //    {
    //        mutType = 2;
    //    }
    //    else
    //    {
    //        mutType = 3;
    //    }
    //    passMutToSubTree(mutationToNodeAssignment, treeStructure, childVector, getHomeNodeId.id, mutation, mutType);
    //}
    
}

unsigned getBin(std::vector<double> & probSubHetero,
        double prob)
{
    double sum = 0;
    for (unsigned i = 0; i < probSubHetero.size(); ++i)
    {
        sum += probSubHetero[i];
        if (prob <= sum)
        {
            return i;
        }
    }
    return UINT_MAX;
}

void makeSubTreeHomozygous(
        std::vector<std::vector<unsigned>> & mutationToNodeAssignment,
        std::vector<unsigned> const & treeStructure,
        std::vector<std::pair<unsigned, unsigned>> const & childVector,
        std::vector<std::vector<unsigned>> & heteroMutsPerNode,
        std::vector<unsigned> const & numPlacements,
        unsigned sumNumPlacemanets,
        Options const & options)
{
    double test = 0;
    std::vector<double> probSubHetero;
    probSubHetero.resize(numPlacements.size());
    for (unsigned node = 0; node < probSubHetero.size(); ++node)
    {
        probSubHetero[node] = static_cast<double>(numPlacements[node])/static_cast<double>(sumNumPlacemanets);
        test += probSubHetero[node];
    }

    unsigned numHeteroMuts = 0;
    for (unsigned i = 0; i < heteroMutsPerNode.size(); ++i)
    {
        if (probSubHetero[i] > 0)
        {
            numHeteroMuts += heteroMutsPerNode[i].size();
        }
    }

    std::vector<std::vector<unsigned>> homoPerNode;
    homoPerNode.resize(heteroMutsPerNode.size());
    unsigned numSampleMuts = std::ceil(static_cast<double>(numHeteroMuts) * options.clam);

    double probFullNodes = 0;
    for (unsigned node = 0; node < numSampleMuts; ++node)
    {
        double prob = static_cast<double>(rand())/static_cast<double>(RAND_MAX) * (1.0 - probFullNodes);
        unsigned nodeId = getBin(probSubHetero, prob);
        if (nodeId == UINT_MAX)
        {
            std::cerr << "Undifined behaviour!" << std::endl;
            std::exit(1);
        }
        unsigned counter = rand() % numPlacements[nodeId];
        GetHomoNodeID getHomeNodeId(counter, nodeId);
        determineHomoNodeId(treeStructure, childVector, nodeId, getHomeNodeId);

        unsigned mutType;
        if (rand() % 2 == 0)
        {
            mutType = 4; // chrom reference
        }
        else
        {
            mutType = 5; // chrom alternative
        }
        passMutToSubTree(mutationToNodeAssignment, treeStructure, childVector, getHomeNodeId.id, heteroMutsPerNode[nodeId][homoPerNode[nodeId].size()], mutType);
        homoPerNode[nodeId].push_back(getHomeNodeId.id);

        if (homoPerNode[nodeId].size() == heteroMutsPerNode[nodeId].size())
        {
            probFullNodes += probSubHetero[nodeId];
            probSubHetero[nodeId] = 0;
        }
    }
}

std::vector<std::vector<unsigned>> assignMutationToNodesSampleTree(std::vector<unsigned> const & treeStructure, 
        std::vector<std::pair<unsigned, unsigned>> const & childVector, 
        std::vector<unsigned> const & numPlacements,
        unsigned sumNumPlacemanets,
        Options const & options)
{
	std::vector<std::vector<unsigned>> mutationToNodeAssignment;
    mutationToNodeAssignment.resize(treeStructure.size());
    for (size_t i = 0; i < mutationToNodeAssignment.size(); ++i)
    {
        mutationToNodeAssignment[i].resize(options.numMutations, 0);
    }
    std::vector<std::vector<unsigned>> heteroMutsPerNode;
    heteroMutsPerNode.resize(mutationToNodeAssignment.size());

    unsigned i = 0;
    if (options.assignMutationsToLeafs)
    {
        for (; i < treeStructure.size() && i < options.numMutations; ++i)
        {
            //assignMutation(mutationToNodeAssignment, treeStructure, childVector, i, i, numPlacements, sumNumPlacemanets, options);
            assignMutation(mutationToNodeAssignment, treeStructure, childVector, i, i, heteroMutsPerNode, options);
        }
        for (; i < options.numMutations; ++i)
        {
            //assignMutation(mutationToNodeAssignment, treeStructure, childVector, rand() % (treeStructure.size()), i, numPlacements, sumNumPlacemanets, options);
            assignMutation(mutationToNodeAssignment, treeStructure, childVector, rand() % (treeStructure.size()), i, heteroMutsPerNode, options);
        }
    }
    else
    {
        for (; i < options.numSamples - 1 && i < options.numMutations; ++i)
        {
            //assignMutation(mutationToNodeAssignment, treeStructure, childVector, i, i, numPlacements, sumNumPlacemanets, options);
            assignMutation(mutationToNodeAssignment, treeStructure, childVector, i, i, heteroMutsPerNode,options);
        }
        for (; i < options.numMutations; ++i)
        {
            //assignMutation(mutationToNodeAssignment, treeStructure, childVector, rand() % (options.numSamples - 1), i, numPlacements, sumNumPlacemanets, options);
            assignMutation(mutationToNodeAssignment, treeStructure, childVector, rand() % (options.numSamples - 1), i, heteroMutsPerNode, options);
        }
    }

    makeSubTreeHomozygous(
        mutationToNodeAssignment,
        treeStructure,
        childVector,
        heteroMutsPerNode,
        numPlacements,
        sumNumPlacemanets,
        options);
    return mutationToNodeAssignment;
}

std::vector<std::vector<unsigned>> extractMutationToSampleAssigments(std::vector<std::vector<unsigned>> const & mutationToNodeAssignment, 
        std::vector<std::pair<unsigned, unsigned>> const & childVector, 
        Options options)
{
    std::vector<std::vector<unsigned>> result;
    result.resize(options.numSamples);
    for (size_t i = 0; i < result.size(); ++i)
        result[i].resize(options.numMutations);

    unsigned counter = 0;
    for (size_t i = 0; i < childVector.size(); ++i)
    {
        if(childVector[i].first == UINT_MAX)
        {
            for (size_t j = 0; j < options.numMutations; ++j)
            {
                result[counter][j] = mutationToNodeAssignment[i][j];
            }
            ++counter;
        }
    }
    return result;
}


unsigned getNumAllMutationsMutTree(std::vector<unsigned> const & sampleToNodeAssignment,
                            std::vector<unsigned> const & treeStructure)
{
    unsigned counter = 0;
	for (unsigned int i = 0; i < sampleToNodeAssignment.size(); ++i)
	{
        for (unsigned node = sampleToNodeAssignment[i]; node > 0; node = treeStructure[node])
		{
            ++counter;
        }
    }
    return counter;
}

//std::vector<std::vector<unsigned>> createMutationToSampleAssigments(std::vector<std::pair<unsigned, unsigned>> const & mutationToNodeAssignment, std::vector<unsigned> const & treeStructure, Options const & options)
//{
//    std::vector<std::vector<unsigned>> mutationToSample;
//    mutationToSample.resize(options.numSamples);
//    for (unsigned i = 0; i < mutationToSample.size(); ++i)
//        mutationToSample[i].resize(options.numMutations, 0);
//
//    // assign the mutations in the root to all samples
//    for (unsigned int i = 0; i < options.numMutations; ++i)
//    {
//        if (mutationToNodeAssignment[i].first == 0)
//        {
//            for (unsigned j = 0; j < options.numSamples; ++j)
//            {
//                if (mutationToNodeAssignment[i].second == 0)
//                {
//                    mutationToSample[j][i] = 1;
//                }
//                else if (mutationToNodeAssignment[i].second == 1)
//                {
//                    std::cout << "TEST" << std::endl;
//                    mutationToSample[j][i] = 3;
//                }
//                else
//                {
//                    std::cerr << "Wrong cygocity!" << std::endl;
//                    std::exit (EXIT_FAILURE);
//                }
//            }
//        }
//    }
//    for (unsigned i = 0; i < options.numSamples; ++i)
//    {
//        unsigned posInTree = options.numSamples - 1 + i;
//        while (posInTree != 0)
//        {
//            for (unsigned j = 0; j < options.numMutations; ++j)
//            {
//                if (mutationToNodeAssignment[j].first == posInTree)
//                {
//                    if (mutationToNodeAssignment[j].second == 0)
//                    {
//                        std::cout << "BEE: " << j << " " << i << std::endl;
//                        mutationToSample[i][j] = 1; // 2 chrom: 1 chrom ref, 1 chrom mut
//                    }
//                    else if (mutationToNodeAssignment[j].second == 1)
//                    {
//                        std::cout << "TEST" << std::endl;
//                        mutationToSample[i][j] = 3; // 1 chrom mut
//                    }
//                    else
//                    {
//                        std::cerr << "Wrong cygocity!" << std::endl;
//                        std::exit (EXIT_FAILURE);
//                    }
//                }
//            }
//            posInTree = treeStructure[posInTree];
//        }
//    }
//    return mutationToSample;
//}
                        
std::vector<std::vector<unsigned>> 
insertDropOuts(std::vector<std::vector<unsigned>> const & mutationToSample, Options const & options)
{
    std::vector<std::vector<unsigned>> mutationToSampleWithDropouts = mutationToSample;

    if (options.dropOutRate > 0)
    {
        for (size_t i = 0; i < mutationToSampleWithDropouts.size(); ++i)
        {
            for (size_t j = 0; j < mutationToSampleWithDropouts[i].size(); ++j)
            {
                if (mutationToSampleWithDropouts[i][j] == 1)
                {
                    if (static_cast<double>(rand())/static_cast<double>(RAND_MAX) < options.dropOutRate)
                    {
                        if(static_cast<double>(rand())/static_cast<double>(RAND_MAX) < 0.5)
                        {
                            mutationToSampleWithDropouts[i][j] = 2; // 1 chrom reference
                        }
                        else
                        {
                            mutationToSampleWithDropouts[i][j] = 3; // 1 chrom mut
                        }
                    }
                }
            }
        }
    }
    return mutationToSampleWithDropouts;
}

void writeMutationAssigmentsSampleTree(std::vector<std::vector<unsigned>> const & mutationToSample,
                        std::vector<unsigned> & positionVec,
                        std::string const & fileName,
						Options const & options)
{
    //std::vector<std::vector<bool>> mutationToSample;
    //createMutationToSampleAssigments(mutationToSample, mutationToNodeAssignment, treeStructure, options);

    std::ofstream outFile;
    outFile.open(toCString(fileName.c_str()));
    outFile << "#position\t";
	for (unsigned int i = options.numSamples - 1; i < 2 * options.numSamples - 1; ++i)
        outFile << "cell" << i << "\t";
    outFile << "\n";
    for (unsigned int j = 0; j < options.numMutations; ++j)
    {
        outFile << positionVec[j] << "\t";
        for (unsigned int i = 0; i < options.numSamples; ++i)
        {
			outFile << mutationToSample[i][j] << "\t";
		}
		outFile << std::endl;
	}
    outFile.close();
}

void writeTree(std::vector<unsigned> sampleToNodeAssignment,
			   std::vector<unsigned> & treeStructure,
			   Options const & options)
{
    std::ofstream outFile;
    outFile.open(toCString(options.outFileTreeStructure));
	outFile << "digraph mutationTree{" << std::endl;
    outFile << "\trankdir = BT;" << std::endl;
	for (unsigned int i = 1; i < treeStructure.size(); ++i)
	{
		outFile << "\t" << i << " -> " << treeStructure[i] << "[dir=back];" << std::endl;
	}
	for (unsigned int i = 0; i < sampleToNodeAssignment.size(); ++i)
	{
		outFile << "\tS" << i << " [shape=box];" << std::endl;
	}
	for (unsigned int i = 0; i < sampleToNodeAssignment.size(); ++i)
	{
		outFile << "\tS" << i << " -> " << sampleToNodeAssignment[i] << "[dir=none];" << std::endl;
	}
    outFile << "}" << std::endl;
    outFile.close();
}

void sampleMutationPositionsOnGenome(std::map<unsigned, unsigned> & positionMap,
                                    std::vector<unsigned> & positionVec,
                                    Options const & options)
{
    unsigned pos;
    auto it = positionMap.begin();
    positionVec.resize(options.numMutations);
    for (size_t i = 0; i < options.numMutations; ++i)
    {
        pos = rand() % options.numPositions + 1;
        it = positionMap.find(pos);
        if (it != positionMap.end())
        {
            --i;
        }
        else
        {
            positionMap.insert(std::pair<unsigned, unsigned>(pos, i));
            positionVec[i] = pos;
        }
    }
}

void writeTreeSampleTree(std::vector<std::vector<unsigned>> const & mutationToNodeAssignment,
			   std::vector<unsigned> const & treeStructure,
			   std::vector<unsigned> const & positionVec,
			   Options const & options)
{
    std::ofstream outFile;
    outFile.open(toCString(options.outFileTreeStructure));
	outFile << "digraph mutationTree{" << std::endl;
    outFile << "\trankdir = BT;" << std::endl;
    unsigned i = 1;
	for (; i < options.numSamples - 1; ++i)
	{
		outFile << "\t" << positionVec[i] << " -> " << positionVec[treeStructure[i]] << "[dir=back];" << std::endl;
	}
	for (; i < 2 * options.numSamples - 1; ++i)
	{
		outFile << "\tS" << i << " -> " << positionVec[treeStructure[i]] << "[dir=none];" << std::endl;
	}

	for (i = 0; i < 2 * options.numSamples - 1; ++i)
    {
        if (i < options.numSamples - 1)
        {
    		outFile << "\t" << positionVec[i] << " [label=\"M:" << i << " "; //[shape=box];" << std::endl;
        }
        else
        {
		    outFile << "\tS" << i << " [label=\"S" << i << " - M:"; // << i << " [shape=box];" << std::endl;
        }
        for (unsigned j = 0; j < options.numMutations; ++j)
        {
            if (mutationToNodeAssignment[i][j] == 1 || mutationToNodeAssignment[i][j] == 3) // hetero, dropout
            {
                outFile  << mutationToNodeAssignment[i][j] << "_" << positionVec[j] << ";";
            }
            else if (mutationToNodeAssignment[i][j] == 2) // homo
            {
                outFile << mutationToNodeAssignment[i][j]<< "_" << positionVec[j] << ";";
            }
        }
        if (i < options.numSamples - 1)
        {
		    outFile << "\"];" << std::endl;
        }
        else
        {
            outFile << "\",shape=box];" << std::endl;
        }
	}

    outFile << "}" << std::endl;
    outFile.close();
}

std::vector<unsigned> getNumPlacements(std::vector<unsigned> const & treeStructure,std::vector<std::pair<unsigned, unsigned>> const & childVector)
{
    std::vector<bool> nodeVisits;
    nodeVisits.resize(treeStructure.size(), 0);
    unsigned currentNode = 0;

    std::vector<unsigned> numPlacements;
    numPlacements.resize(treeStructure.size(), 0);

    
    while(!nodeVisits[0])
    {
        if (childVector[currentNode].first != UINT_MAX && !nodeVisits[childVector[currentNode].first])
        {
            currentNode = childVector[currentNode].first;
        }
        else if (childVector[currentNode].second != UINT_MAX && !nodeVisits[childVector[currentNode].second])
        {
            currentNode = childVector[currentNode].second;
        }
        else
        {
            if (childVector[currentNode].first != UINT_MAX && childVector[currentNode].second != UINT_MAX)
            {
                numPlacements[currentNode] = numPlacements[childVector[currentNode].first] + numPlacements[childVector[currentNode].second];
                if (childVector[childVector[currentNode].first].first != UINT_MAX)
                {
                    ++numPlacements[currentNode];
                }
                if (childVector[childVector[currentNode].second].first != UINT_MAX)
                {
                    ++numPlacements[currentNode];
                }
            }
            nodeVisits[currentNode] = 1;
            currentNode = treeStructure[currentNode];
        }
    }
    return numPlacements;
}


int main(int argc, const char* argv[]){

    //---------------------------------------------------------------
    // parse the command line.
    //---------------------------------------------------------------
    ArgumentParser parser;
    Options options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    //std::cout << "create_single_cell_seq_tree\n"
    //          << "================\n\n";

    srand(options.seed + 2); // srand(1) is reserved

	std::vector<unsigned> treeStructure;
	
    if (options.numSamples > options.numMutations)
    {
	    createTreeStructure(treeStructure, options.numMutations);

	    std::vector<unsigned> sampleToNodeAssignment;
        assignSamplesToNodesMutTree(sampleToNodeAssignment, treeStructure, options.numSamples);

        //std::vector<std::vector<std::tuple<unsigned, unsigned>>> mutationCalls;
        //createMutationCalls(mutationCalls,
        //                    sampleToNodeAssignment,
        //                    treeStructure,
        //                    options);

        //writeMutationCalls(mutationCalls,
        //                   sampleToNodeAssignment,
        //                   treeStructure,
        //                   options);

        writeTree(sampleToNodeAssignment,
                  treeStructure,
                  options);

        //writePileups(mutationCalls,
        //           sampleToNodeAssignment,
        //           treeStructure,
        //           options);
    }
    else
    {

	    createTreeStructureSampleTree(treeStructure, options.numSamples);

        std::vector<std::pair<unsigned, unsigned>> childVector = treeStructureToChildVector(treeStructure);
        
        std::vector<unsigned> numPlacements = getNumPlacements(treeStructure, childVector);
        unsigned sumNumPlacemanets = 0;
        for (unsigned i = 0; i < numPlacements.size(); ++i)
        {
            sumNumPlacemanets += numPlacements[i];
        }
        //std::cout << std::endl;

	    std::vector<std::vector<unsigned>> mutationToNodeAssignment = assignMutationToNodesSampleTree(treeStructure, childVector, numPlacements, sumNumPlacemanets, options);

        std::vector<std::vector<unsigned>> mutationToSample = extractMutationToSampleAssigments(mutationToNodeAssignment, childVector, options);

        std::vector<std::vector<unsigned>> mutationToSampleWithDropouts = insertDropOuts(mutationToSample, options);

        std::map<unsigned, unsigned> positionMap;
        std::vector<unsigned> positionVec;
        sampleMutationPositionsOnGenome(positionMap,
                                        positionVec,
                                        options);

        writeMutationAssigmentsSampleTree(mutationToSampleWithDropouts,
                        positionVec,
                        toCString(options.outFileMutationAssignments),
						options);

        writeTreeSampleTree(mutationToNodeAssignment,
                  treeStructure,
                  positionVec,
                  options);

    }

	return 0;	
}

//bool mutationAllowedMutTree(std::vector<std::vector<std::tuple<unsigned, unsigned>>> const & mutationCalls,
//        unsigned sample,
//        unsigned mutation)
//{
//    unsigned numMutations = 0;
//    for (unsigned i = 0; i < mutationCalls.size(); ++i)
//    {
//        double freq = static_cast<double>(std::get<1>(mutationCalls[i][mutation])) / static_cast<double>(std::get<0>(mutationCalls[i][mutation]));
//        if (i != sample && 0.3 < freq && freq < 0.7)
//            ++numMutations;
//    }
//    return numMutations >= 2;
//}



//void createMutationCalls(std::vector<std::vector<std::tuple<unsigned, unsigned>>> & mutationCalls,
//						 std::vector<unsigned> sampleToNodeAssignment,
//						 std::vector<unsigned> & treeStructure,
//						 Options const & options)
//{
//
//  	std::default_random_engine generator (options.seed);
//	
//    // create the coverage
//	mutationCalls.resize(sampleToNodeAssignment.size());
//
//    double prob = (options.covMean / options.covVar);
//    double r = options.covMean / (1.0 - options.covMean/options.covVar) - options.covMean;
//
//    double randFromUnif;
//	for (unsigned int i = 0; i < sampleToNodeAssignment.size(); ++i)
//	{
//		mutationCalls[i].resize(treeStructure.size());
//		for (unsigned int j = 0; j < treeStructure.size(); ++j)
//		{
//            boost::math::negative_binomial_distribution<> negBinCov(r, prob);
//            randFromUnif = static_cast<double>(rand())/static_cast<double>(RAND_MAX); 
//            double randFromNegBinCov = quantile(negBinCov, randFromUnif);
//
//            std::get<0>(mutationCalls[i][j]) = randFromNegBinCov;
//        }
//    }
//
//    //create the wild type counts
//    double wildAlpha = options.wildMean * options.wildVar;
//    double wildBeta = options.wildVar - options.wildMean * options.wildVar;
//	for (unsigned int i = 0; i < sampleToNodeAssignment.size(); ++i)
//	{
//		mutationCalls[i].resize(treeStructure.size());
//		for (unsigned int j = 0; j < treeStructure.size(); ++j)
//		{
//            randFromUnif = static_cast<double>(rand())/static_cast<double>(RAND_MAX); 
//            boost::math::beta_distribution<> betaWild(wildAlpha, wildBeta);
//            double randFromBetaWild = quantile(betaWild, randFromUnif);
//
//            boost::math::binomial_distribution<> betaBinWild(std::get<0>(mutationCalls[i][j]), randFromBetaWild);
//            randFromUnif = static_cast<double>(rand())/static_cast<double>(RAND_MAX); 
//            double randFromBetaBinWild = quantile(betaBinWild, randFromUnif);
//
//            std::get<1>(mutationCalls[i][j]) = randFromBetaBinWild;
//        }
//    }
//
//    //create the mutation counts
//    double altAlpha = options.altMean * options.altVar;
//    double altBeta = options.altVar - options.altMean * options.altVar;
//	for (unsigned int i = 0; i < sampleToNodeAssignment.size(); ++i)
//	{
//        for (unsigned node = sampleToNodeAssignment[i]; node > 0; node = treeStructure[node])
//		{
//            randFromUnif = static_cast<double>(rand())/static_cast<double>(RAND_MAX); 
//            boost::math::beta_distribution<> betaAlt(altAlpha, altBeta);
//            double randFromBetaAlt = quantile(betaAlt, randFromUnif);
//
//            boost::math::binomial_distribution<> betaBinAlt(std::get<0>(mutationCalls[i][node]), randFromBetaAlt);
//            randFromUnif = static_cast<double>(rand())/static_cast<double>(RAND_MAX); 
//            double randFromBetaBinAlt = quantile(betaBinAlt, randFromUnif);
//
//            std::get<1>(mutationCalls[i][node]) = randFromBetaBinAlt;
//        }
//    }
//
//    // create the missing information
//    unsigned numErrors = options.missingInformation * options.numSamples * options.numMutations;
//    std::cout << "num Missing: " << numErrors << std::endl;
//	for (unsigned int i = 0; i < numErrors; )
//	{
//        unsigned sample = rand() % options.numSamples;
//        unsigned mutation = rand() % options.numMutations;
//        if (mutationAllowedMutTree(mutationCalls, sample, mutation))
//        {
//            ++i;
//            std::get<0>(mutationCalls[sample][mutation]) = 0;
//            std::get<1>(mutationCalls[sample][mutation]) = 0;
//        }
//	}
//
//    //create drop out rate
//    numErrors = options.dropOutRate * getNumAllMutationsMutTree(sampleToNodeAssignment, treeStructure);
//    std::cout << "num drop out: " << numErrors << std::endl;
//	for (unsigned int i = 0; i < numErrors; )
//	{
//        unsigned sample = rand() % options.numSamples;
//        unsigned mutation = rand() % options.numMutations;
//        if (std::get<0>(mutationCalls[sample][mutation]) != 0)
//        {
//            double freq = static_cast<double>(std::get<1>(mutationCalls[sample][mutation])) / static_cast<double>(std::get<0>(mutationCalls[sample][mutation]));
//            if (0.3 < freq && freq < 0.7 && mutationAllowedMutTree(mutationCalls, sample, mutation))
//            {
//                if ((rand() % 2) == 0)
//                {
//                    std::binomial_distribution<int> binDist (std::get<0>(mutationCalls[sample][mutation]), 0.01);
//                    std::get<1>(mutationCalls[sample][mutation]) = binDist(generator);
//                }
//                else
//                {
//                    std::binomial_distribution<int> binDist (std::get<0>(mutationCalls[sample][mutation]), 0.9);
//                    std::get<1>(mutationCalls[sample][mutation]) = binDist(generator);
//                }
//                ++i;
//            }
//		}
//	}
//}
//
//void createMutationCallsSampleTree(std::vector<std::vector<std::tuple<unsigned, unsigned>>> & mutationCalls,
//                         std::vector<std::vector<bool>> const & mutationToSample,
//						 Options const & options)
//{
//
//  	std::default_random_engine generator (options.seed);
//	
//    // create the coverage
//	mutationCalls.resize(options.numSamples);
//
//    double prob = (options.covMean / options.covVar);
//    double r = options.covMean / (1.0 - options.covMean/options.covVar) - options.covMean;
//
//    double randFromUnif;
//	for (unsigned int i = 0; i < mutationCalls.size(); ++i)
//	{
//		mutationCalls[i].resize(options.numMutations);
//		for (unsigned int j = 0; j < mutationCalls[i].size(); ++j)
//		{
//            boost::math::negative_binomial_distribution<> negBinCov(r, prob);
//            randFromUnif = static_cast<double>(rand())/static_cast<double>(RAND_MAX); 
//            double randFromNegBinCov = quantile(negBinCov, randFromUnif);
//
//            std::get<0>(mutationCalls[i][j]) = randFromNegBinCov;
//        }
//    }
//
//    //create the wild type counts
//    double wildAlpha = options.wildMean * options.wildVar;
//    double wildBeta = options.wildVar - options.wildMean * options.wildVar;
//
//    //create the mutation counts
//    double altAlpha = options.altMean * options.altVar;
//    double altBeta = options.altVar - options.altMean * options.altVar;
//
//	for (unsigned int i = 0; i < mutationToSample.size(); ++i)
//	{
//		for (unsigned int j = 0; j < mutationToSample[i].size(); ++j)
//		{
//            // create the background noise
//            if (mutationToSample[i][j] == 0)
//            {
//                randFromUnif = static_cast<double>(rand())/static_cast<double>(RAND_MAX); 
//                boost::math::beta_distribution<> betaWild(wildAlpha, wildBeta);
//                double randFromBetaWild = quantile(betaWild, randFromUnif);
//
//                boost::math::binomial_distribution<> betaBinWild(std::get<0>(mutationCalls[i][j]), randFromBetaWild);
//                randFromUnif = static_cast<double>(rand())/static_cast<double>(RAND_MAX); 
//                double randFromBetaBinWild = quantile(betaBinWild, randFromUnif);
//
//                std::get<1>(mutationCalls[i][j]) = randFromBetaBinWild;
//            }
//            // create the mutation nucleotide counts
//            else
//            {
//                randFromUnif = static_cast<double>(rand())/static_cast<double>(RAND_MAX); 
//                boost::math::beta_distribution<> betaAlt(altAlpha, altBeta);
//                double randFromBetaAlt = quantile(betaAlt, randFromUnif);
//
//                boost::math::binomial_distribution<> betaBinAlt(std::get<0>(mutationCalls[i][j]), randFromBetaAlt);
//                randFromUnif = static_cast<double>(rand())/static_cast<double>(RAND_MAX); 
//                double randFromBetaBinAlt = quantile(betaBinAlt, randFromUnif);
//
//                std::get<1>(mutationCalls[i][j]) = randFromBetaBinAlt;
//
//            }
//        }
//    }
//}
//
//void writeMutationCalls(std::vector<std::vector<std::tuple<unsigned, unsigned>>> & mutationCalls,
//						std::vector<unsigned> sampleToNodeAssignment,
//						std::vector<unsigned> & treeStructure,
//						Options const & options)
//{
//    std::ofstream outFile;
//    outFile.open(toCString(options.outFileMutations));
//    outFile << "#";
//	for (unsigned int i = 0; i < sampleToNodeAssignment.size(); ++i)
//        outFile << "cell" << i << "\t";
//    outFile << "\n";
//    for (unsigned int j = 1; j < treeStructure.size(); ++j) //do not write out the root
//    {
//        for (unsigned int i = 0; i < sampleToNodeAssignment.size(); ++i)
//        {
//			outFile << std::get<0>(mutationCalls[i][j]) << "\t" << std::get<1>(mutationCalls[i][j]) << "\t";
//		}
//		outFile << std::endl;
//	}
//    outFile.close();
//}
//						 
//void writeMutationCallsSampleTree(std::vector<std::vector<std::tuple<unsigned, unsigned>>> & mutationCalls,
//						std::vector<unsigned> mutationToNodeAssignment,
//						Options const & options)
//{
//    std::ofstream outFile;
//    outFile.open(toCString(options.outFileMutations));
//    outFile << "#";
//	for (unsigned int i = options.numSamples - 1; i < 2 * options.numSamples - 1; ++i)
//        outFile << "cell" << i << "\t";
//    outFile << "\n";
//    for (unsigned int j = 0; j < options.numMutations; ++j) //do not write out the root
//    {
//        for (unsigned int i = 0; i < options.numSamples; ++i)
//        {
//			outFile << std::get<0>(mutationCalls[i][j]) << "\t" << std::get<1>(mutationCalls[i][j]) << "\t";
//		}
//		outFile << std::endl;
//	}
//    outFile.close();
//}
//
                           

//void writePileups(std::vector<std::vector<std::tuple<unsigned, unsigned>>> & mutationCalls,
//                  std::vector<unsigned> sampleToNodeAssignment,
//                  std::vector<unsigned> & treeStructure,
//                  Options const & options)
//{
//    for (unsigned int i = 0; i < sampleToNodeAssignment.size(); ++i)
//    {
//        String<char> outName = options.pileupPrefix;
//        append(outName, "_cell");
//        append(outName, std::to_string(i));
//        append(outName, ".pileup");
//
//        std::ofstream outFile;
//        outFile.open(toCString(outName));
//        for (unsigned int j = 1; j < treeStructure.size(); ++j) //do not write out the root
//        {
//            outFile << "chr" << j << "\t1000\tA\t" << std::get<0>(mutationCalls[i][j]) << "\t" ;
//            unsigned k = 0;
//            for (; k < std::get<1>(mutationCalls[i][j]); ++k)
//            {
//                outFile << (rand() % 2 == 1 ? "C" : "c");
//            }
//            for (; k < std::get<0>(mutationCalls[i][j]); ++k)
//            {
//                outFile << (rand() % 2 == 1 ? "." : ",");
//            }
//            outFile << "\t";
//            k = 0;
//            for (; k < std::get<0>(mutationCalls[i][j]); ++k)
//            {
//                outFile << "D";
//            }
//            outFile << "\n";
//
//        }
//        outFile.close();
//    }
//}
//
//
//void writePileupsSampleTree(std::vector<std::vector<std::tuple<unsigned, unsigned>>> const & mutationCalls,
//                  std::map<unsigned, unsigned> & positions,
//                  Options const & options)
//{
//    double prob = (options.covMean / options.covVar);
//    double r = options.covMean / (1.0 - options.covMean/options.covVar) - options.covMean;
//
//    // variables for the coverage
//    double randFromUnif;
//    double randFromNegBinCov;
//
//    // variables for the alternative counts
//    double wildAlpha = options.wildMean * options.wildVar;
//    double wildBeta = options.wildVar - options.wildMean * options.wildVar;
//    double randFromBetaWild;
//    double randFromBetaBinWild;
//
//    for (unsigned int i = 0; i < mutationCalls.size(); ++i)
//    {
//        String<char> outName = options.pileupPrefix;
//        append(outName, "_cell");
//        append(outName, std::to_string(i));
//        append(outName, ".pileup");
//
//        std::ofstream outFile;
//        outFile.open(toCString(outName));
//
//        auto it = positions.begin();
//        unsigned cov = 0;
//        unsigned alt = 0;
//        for (unsigned pos = 0; pos < options.numPositions; ++pos)
//        {
//            //create the coverage
//            boost::math::negative_binomial_distribution<> negBinCov(r, prob);
//            randFromUnif = static_cast<double>(rand())/static_cast<double>(RAND_MAX); 
//            randFromNegBinCov = quantile(negBinCov, randFromUnif);
//            cov = randFromNegBinCov;
//
//            //create the noise
//            randFromUnif = static_cast<double>(rand())/static_cast<double>(RAND_MAX); 
//            boost::math::beta_distribution<> betaWild(wildAlpha, wildBeta);
//            randFromBetaWild = quantile(betaWild, randFromUnif);
//
//            boost::math::binomial_distribution<> betaBinWild(cov, randFromBetaWild);
//            randFromUnif = static_cast<double>(rand())/static_cast<double>(RAND_MAX); 
//            randFromBetaBinWild = quantile(betaBinWild, randFromUnif);
//            alt = randFromBetaBinWild;
//
//            it = positions.find(pos);
//            if (it != positions.end())
//            {
//                cov = std::get<0>(mutationCalls[i][it->second]);
//                alt = std::get<1>(mutationCalls[i][it->second]);
//            }
//            
//            outFile << "chr1" << "\t" << pos << "\tA\t";
//            outFile << cov << "\t" ;
//            unsigned k = 0;
//            for (; k < alt; ++k)
//            {
//                outFile << (rand() % 2 == 1 ? "C" : "c");
//            }
//            for (; k < cov; ++k)
//            {
//                outFile << (rand() % 2 == 1 ? "." : ",");
//            }
//            outFile << "\t";
//            k = 0;
//            for (; k < cov; ++k)
//            {
//                outFile << "D";
//            }
//            outFile << "\n";
//            
//        }
//        outFile.close();
//    }
//}
//


//void writeParameters(Options const & options)
//{
//    double wildAlpha = options.wildMean * options.wildVar;
//    double wildBeta = options.wildVar - options.wildMean * options.wildVar;
//    double altAlpha = options.altMean * options.altVar;
//    double altBeta = options.altVar - options.altMean * options.altVar;
//
//    std::ofstream outFile;
//    outFile.open(toCString(options.outFileParameters));
//	outFile << "background frequency:\t" << options.wildMean << "\t" << options.wildVar << std::endl;
//	outFile << "wildAlpha:\t" << wildAlpha << std::endl;
//	outFile << "wildBeta:\t" << wildBeta << std::endl;
//	outFile << "average mutation frequency:\t" << options.altMean << "\t" << options.altVar << std::endl;
//	outFile << "altAlpha:\t" << altAlpha << std::endl;
//	outFile << "altBeta:\t" << altBeta << std::endl;
//}


