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
    double cpn; // fraction of mutations with copy number changes
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
        cpn(0.0),
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
	addOption(parser, ArgParseOption("cpn", "CopyNumber", "The copy number rate.", ArgParseArgument::DOUBLE, "COPYNUMBERRATE"));
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
    if (isSet(parser, "cpn"))
        getOptionValue(options.cpn, parser, "cpn");
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

//struct GetSize
//{
//    unsigned size;
//    unsigned startNode = UINT_MAX;
//
//    GetSize() : size(0) {};
//    GetSize(unsigned startNode_) : size(0), startNode(startNode_) {};
//
//    bool operator()() 
//    {
//        ++this->size;
//        return false;
//    };
//    bool operator()(unsigned currentNode, bool isLeaf)
//    {
//        if(!isLeaf && currentNode != this -> startNode)
//        {
//            ++this->size;
//        }
//        return false;
//    };
//};
//
//struct GetHomoNodeID
//{
//    unsigned nodePos = UINT_MAX;
//    unsigned rootNode = UINT_MAX;
//    unsigned counter = 0;
//    unsigned id = UINT_MAX;
//
//    GetHomoNodeID(unsigned counter, unsigned rootNode_) :nodePos(counter), rootNode(rootNode_) {};
//
//    bool operator()(unsigned nodeID, bool isLeaf) 
//    {
//        if (!isLeaf)
//        {
//            if (nodePos == counter)
//            {
//                id = nodeID;
//                return true;
//            }
//        ++this->counter;
//        }
//        return false;
//    };
//};
//
//
//template <typename TType>
//void determineHomoNodeId(std::vector<unsigned> const & treeStructure,
//        std::vector<std::pair<unsigned, unsigned>> const & childVector,
//        unsigned node,
//        TType & obj)
//{
//    std::vector<bool> nodeVisits;
//    nodeVisits.resize(treeStructure.size(), 0);
//    unsigned currentNode = node;
//    unsigned test = 0;
//    while(!nodeVisits[node])
//    {
//        if (childVector[currentNode].first != UINT_MAX && !nodeVisits[childVector[currentNode].first])
//        {
//            currentNode = childVector[currentNode].first;
//        }
//        else if (childVector[currentNode].second != UINT_MAX && !nodeVisits[childVector[currentNode].second])
//        {
//            currentNode = childVector[currentNode].second;
//        }
//        else
//        {
//            nodeVisits[currentNode] = 1;
//            if (obj(currentNode, childVector[currentNode].first == UINT_MAX))
//            {
//                break;
//            }
//            currentNode = treeStructure[currentNode];
//            ++test;
//        }
//    }
//}

void passMutToSubTree(std::vector<std::vector<int>> & mutationToNodeAssignment,
        std::vector<unsigned> const & treeStructure,
        std::vector<std::pair<unsigned, unsigned>> const & childVector,
        unsigned node,
        unsigned mutation,
        int mutType)
{
    unsigned currentNode = node;
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
        std::vector<std::vector<int>> & mutationToNodeAssignment,
        std::vector<unsigned> const & treeStructure,
        std::vector<std::pair<unsigned, unsigned>> const & childVector,
        unsigned node,
        unsigned mutation,
        std::vector<std::vector<unsigned>> & heteroMutsPerNode,
        Options const & options)
{
    int mutType = 1;
    if (static_cast<double>(rand())/static_cast<double>(RAND_MAX) < options.clbm)
    {
	    mutType = 3;
    }
    else if (static_cast<double>(rand())/static_cast<double>(RAND_MAX) < options.cpn)
    {
        double r = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
        int cp = 5 + std::ceil(std::log2(1.0/(1.0-r)));
        if (static_cast<double>(rand())/static_cast<double>(RAND_MAX) >= 0.0) // change this if the mutated allele should be duplicated
        {
            cp *= -1;
        }
        mutType = cp;
    }

    passMutToSubTree(mutationToNodeAssignment, treeStructure, childVector, node, mutation, mutType);

}

std::vector<std::vector<int>> assignMutationToNodesSampleTree(std::vector<unsigned> const & treeStructure, 
        std::vector<std::pair<unsigned, unsigned>> const & childVector, 
        std::vector<unsigned> const & numPlacements,
        unsigned sumNumPlacemanets,
        Options const & options)
{
	std::vector<std::vector<int>> mutationToNodeAssignment;
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

    //makeSubTreeHomozygous(
    //    mutationToNodeAssignment,
    //    treeStructure,
    //    childVector,
    //    heteroMutsPerNode,
    //    numPlacements,
    //    sumNumPlacemanets,
    //    options);
    return mutationToNodeAssignment;
}

std::vector<std::vector<int>> extractMutationToSampleAssigments(std::vector<std::vector<int>> const & mutationToNodeAssignment, 
        std::vector<std::pair<unsigned, unsigned>> const & childVector, 
        Options options)
{
    std::vector<std::vector<int>> result;
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

int copyNumberDropOut(int cp, Options const & options)
{
    unsigned absCp = std::abs(cp) - 3;
    std::vector<bool> allels(absCp, 1);
    unsigned sum = 0;

    while(sum == 0)
    {
        for(unsigned i = 0; i < absCp; ++i)
        {
            static_cast<double>(rand())/static_cast<double>(RAND_MAX) < options.dropOutRate/2.0 ? allels[i] = 0 : allels[i] = 1;
        }

        sum = 0;
        std::for_each(allels.begin(), allels.end(), [&] (bool n)
        {
            sum += n;
        });
    }

    unsigned ref = 0;
    unsigned mut = 0;
    if (cp > 0)
    {
        ref = allels[0];
        if (ref == 0)
        {
            return 3;
        }
        std::for_each(allels.begin() + 1, allels.end(), [&] (bool n)
        {
            mut += n;
        });
        if (mut == 0)
        {
            return 2;
        }
        if (mut == 1)
        {
            return 1;
        }
        return 4 + mut;
    }
    else
    {
        mut = allels.back();
        if (mut == 0)
        {
            return 2;
        }
        std::for_each(allels.begin(), allels.end() - 1, [&] (bool n)
        {
            ref += n;
        });
        if (ref == 0)
        {
            return 3;
        }
        if (ref == 1)
        {
            return 1;
        }
        return -1 * (4 + ref);
    }
}

std::vector<std::vector<int>> 
insertDropOuts(std::vector<std::vector<int>> const & mutationToSample, Options const & options)
{
    std::vector<std::vector<int>> mutationToSampleWithDropouts = mutationToSample;

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
                if (mutationToSampleWithDropouts[i][j] < 0 || mutationToSampleWithDropouts[i][j] >= 6)
                {
                    std::cout << "0: " << mutationToSampleWithDropouts[i][j] << std::endl;
                    mutationToSampleWithDropouts[i][j] = copyNumberDropOut(mutationToSampleWithDropouts[i][j], options);
                    std::cout << "1: " << mutationToSampleWithDropouts[i][j] << std::endl;
                }
            }
        }
    }
    return mutationToSampleWithDropouts;
}

void writeMutationAssigmentsSampleTree(std::vector<std::vector<int>> const & mutationToSample,
                        std::vector<unsigned> & positionVec,
                        std::string const & fileName,
						Options const & options)
{
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

void writeTreeSampleTree(std::vector<std::vector<int>> const & mutationToNodeAssignment,
			   std::vector<unsigned> const & treeStructure,
			   std::vector<unsigned> const & positionVec,
			   Options const & options)
{
    std::ofstream outFile;
    outFile.open(toCString(options.outFileTreeStructure));
	outFile << "digraph mutationTree{" << std::endl;
    outFile << "rankdir = BT;" << std::endl;
    unsigned i = 1;
	for (; i < options.numSamples - 1; ++i)
	{
		outFile << i << " -> " << treeStructure[i] << "[dir=back];" << std::endl;
	}
	for (; i < 2 * options.numSamples - 1; ++i)
	{
		outFile << i << " -> " << treeStructure[i] << "[dir=none];" << std::endl;
	}

	for (i = 0; i < 2 * options.numSamples - 1; ++i)
    {
        if (i < options.numSamples - 1)
        {
    		outFile << i << " [label=\"" << i << "\\n"; //[shape=box];" << std::endl;
        }
        else
        {
		    outFile << i << " [label=\"" << i << ".bam\\n"; // << i << " [shape=box];" << std::endl;
        }
        for (unsigned j = 0; j < options.numMutations; ++j)
        {
            outFile << mutationToNodeAssignment[i][j] << "_" << positionVec[j] << "\\n";
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
	
    createTreeStructureSampleTree(treeStructure, options.numSamples);

    std::vector<std::pair<unsigned, unsigned>> childVector = treeStructureToChildVector(treeStructure);
    
    std::vector<unsigned> numPlacements = getNumPlacements(treeStructure, childVector);
    unsigned sumNumPlacemanets = 0;
    for (unsigned i = 0; i < numPlacements.size(); ++i)
    {
        sumNumPlacemanets += numPlacements[i];
    }

    std::vector<std::vector<int>> mutationToNodeAssignment = assignMutationToNodesSampleTree(treeStructure, childVector, numPlacements, sumNumPlacemanets, options);

    std::vector<std::vector<int>> mutationToSample = extractMutationToSampleAssigments(mutationToNodeAssignment, childVector, options);

    std::vector<std::vector<int>> mutationToSampleWithDropouts = insertDropOuts(mutationToSample, options);

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

	return 0;	
}

