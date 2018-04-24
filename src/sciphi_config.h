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
#ifndef CONFIG_H
#define CONFIG_H

#include <boost/graph/graphviz.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/dynamic_bitset.hpp>

#include <array>
#include <vector>
#include <tuple>
#include <random>
#include <utility>
#include <cassert>
#include <unordered_map>

#include "noise_counts.h"
#include "logScores.h"
#include "attachmentScores.h"


struct MutationTree {};
struct SampleTree {};
struct SimpleTree {};

template <typename TTreeType>
struct Vertex{};

template <>
struct Vertex<MutationTree> {
    unsigned mutation;
};

template <>
struct Vertex<SampleTree> {
    int sample = -1;  // -1 indicates this is an inner node
};
template <>
struct Vertex<SimpleTree> {
    std::vector<unsigned> mutations;
    int sample = -1;
};

// This class is used to get statistics on the parameter learned
struct ParamsCounter
{
    std::vector<double> mu;
    std::vector<double> nu;
    std::vector<double> lambda; // experimental
    std::vector<double> wildAlpha;
    std::vector<double> wildBeta;
    std::vector<double> mutAlpha;
    std::vector<double> mutBeta;

    ParamsCounter(){};

    unsigned resize(unsigned newSize)
    {
        mu.resize(newSize, 0);
        nu.resize(newSize, 0);
        lambda.resize(newSize, 0);
        wildAlpha.resize(newSize, 0);
        wildBeta.resize(newSize, 0);
        mutAlpha.resize(newSize, 0);
        mutBeta.resize(newSize, 0);
        return newSize;
    }
};

// This class combines all data structures such that they dont have to be passed between functions
// TODO: create several sub-classes and use OOP!!!!!!!
template <typename TTreeType>
class Config{

	public:

    // Some type definitions
    typedef boost::adjacency_list<boost::vecS, 
                                  boost::vecS, 
                                  boost::bidirectionalS,   
                                  Vertex<TTreeType>>                    TGraph;                 
	typedef std::vector<bool>                                           TSampleNodes;           
	typedef LogScores                                                   TLogScores;
    
	typedef AttachmentScores<5>                                         TAttachmentScores;
	typedef std::vector<std::vector<std::tuple<unsigned, unsigned> > >  TData;
    typedef std::tuple<double, double>                                  TParamsTuple;
    typedef std::array<TParamsTuple, 6>                                 TParams;
    typedef std::tuple<double, unsigned, unsigned >                     TLearningParamsTuple;
    typedef std::vector<TLearningParamsTuple>                           TLearningParams;
    typedef boost::dynamic_bitset<>                                     TBitSet;

    // A practical enum for easy access to the parameters
    enum ParamType{
        wildOverDis = 0, 
        mutationOverDis = 1, 
        wildMean = 2, 
        mu = 3, 
        nu = 4, 
        lambda = 5, 
        mutationMean = 6};

    // Manu small helper functions, mostly getter and setter
    void updateParamsCounter();
    
    TGraph & getTree();
    TGraph const & getTree() const;
    void setTree(TGraph const & newTree);

    TGraph & getTmpTree();
    void setTmpTree(TGraph const & newTree);

	unsigned getNumSamples() const;
	void setNumSamples(unsigned newNumSamples);
	
	ParamType getParamToOptimize() const;
	void setParamToOptimize(ParamType newParamToOptimize);
	
    unsigned getNumMutations() const;
    unsigned getNumAttachments() const;

	TSampleNodes & getSampleNodes();
	TSampleNodes const & getSampleNodes() const;

	TLogScores & getLogScores();
	TLogScores const & getLogScores() const;
	void setLogScores(TLogScores newLogScores);
    
    TLogScores & getTmpLogScores();
	void setTmpLogScores(TLogScores & newTmpLogScores);

	TData & getData();
	TData const & getData() const;

	TData & getCompleteData();
	TData const & getCompleteData() const;

	TAttachmentScores & getTmpAttachmentScore();
	TAttachmentScores const & getTmpAttachmentScore() const;
	
    double getParam(ParamType param);
    double getParam(ParamType param) const;
	void setParam(ParamType param, double newParam);

    double getTmpParam(ParamType param);
	void setTmpParam(ParamType param, double newParam);
    
    double getSDParam(ParamType param);
	void setSDParam(ParamType param, double newParam);

    unsigned getSDCountParam(ParamType param);
	void setSDCountParam(ParamType param, unsigned newParam);
    
    unsigned getSDTrialsParam(ParamType param);
	void setSDTrialsParam(ParamType param, unsigned newParam);

    unsigned getMoveTyp();
    void setMoveTyp(unsigned newMoveType);

	int getOverDisInterval();

    void resetParameters();

    void printParameters();

    bool updateContainers(unsigned currentLoop);

    void initMutInSampleCounter();

	std::default_random_engine & getGenerator();

    // THe members of the config class
    // the main tree structure - we make use of a pair to restore 
    // a tree if the proposed one is discarded
    std::pair<TGraph, TGraph>                   tree;
    
    // the parameters of the tree
    TParams                                     params;                 // <overDisWild, overDisWildTmp>
                                                                        // <overDisMutation, overDisMutationTmp>
                                                                        // <wildMean, wildMeanTmp>
    
    // the parameters for learning the tree parameters               
    TLearningParams                             learningParams;         // <sdWildOverDis, numTrails, numSuc>
                                                                        // <sdMutationOverDis, numTrails, numSuc>
                                                                        // <sdWildMean, numTrails, numSuc>

    // the probabilities with which one of the moves 
    // (including parameter estimation) is performed
    std::array<double, 4>                       moveProbs;
	std::default_random_engine                  generator; //(seed);
    // the log scores for each mutation for each sample
    unsigned                                    loops;
    unsigned                                    reps;
    unsigned                                    fixedSeed;
    char                                        scoreType = 'm';
    double                                      paramsEstimateRate;
    double                                      priorMutationRate;
    unsigned                                    uniqTreshold;
    std::tuple<double, double>                  dataUsageRate; 
    unsigned                                    sampleLoops;
    unsigned                                    errorRateEstLoops;
    
    NoiseCounts                                 noiseCounts;
    std::pair<TLogScores, TLogScores>			logScores;
   
	TSampleNodes	                            sampleNodes;
	TData 				                        data;
	TData 				                        completeData;
    TAttachmentScores                           _tmpAttachmentScore;
	unsigned 			                        numSamples;
    unsigned                                    moveType;

    ParamType                                   paramToOptimize;
    std::string                                 outFilePrefix;
    std::string                                 refFileName;
    std::string                                 bamFileNames;
    std::string                                 exclusionFileName;
    std::string                                 mutationExclusionFileName;
    std::string                                 inFileName;
    std::vector<std::string>                    cellNames;
    std::vector<std::tuple<std::string, unsigned, char, char>> indexToPosition;
    std::vector<std::vector<TAttachmentScores::TAttachmentScore>> mutInSampleCounter;
    unsigned                                    numUniqMuts;
    unsigned                                    numCellWithMutationMin;
    unsigned                                    minCoverage;
    unsigned                                    minCoverageAcrossCells;
    unsigned                                    numMinCoverageAcrossCells;
    unsigned                                    minSupport;
    double                                      minFreq;
    unsigned                                    minCovInControlBulk;
    unsigned                                    maxSupInControlBulk;
    ParamsCounter                               paramsCounter;
    double                                      noiseScore;
    bool                                        learnZygocity;
    bool                                        computeMixScore;
    std::array<unsigned, 2>                     numMutPlacements;
    bool                                        estimateSeqErrorRate;
    std::array<double, 2>                       clamPrior;
    double                                      meanFilter;      
    unsigned                                    minCovNormalCell;
    unsigned                                    maxNumberNormalCellMutated;
    bool                                        useNormalCellInTree;
    //double mu;


	Config() :
        params{{TParamsTuple{100.0,100.0}, 
                TParamsTuple{2, 2}, 
                TParamsTuple{0.001, 0.001},
                TParamsTuple{0.9, 0.9},
                TParamsTuple{0, 0},
                TParamsTuple{0, 0}}},
        learningParams{{TLearningParamsTuple{5.0, 0, 0}, 
                        TLearningParamsTuple{0.1, 0, 0}, 
                        TLearningParamsTuple{0.01, 0, 0},
                        TLearningParamsTuple{0.01, 0, 0},
                        TLearningParamsTuple{0.01, 0, 0},
                        TLearningParamsTuple{0.01, 0, 0}}},
        moveProbs{{0.64, 0.16, 0 , 0.2}},
        //overDisInterval(5),
        generator(42),
        loops(1000000),
        reps(1),
        fixedSeed(42),
        scoreType('s'),
        paramsEstimateRate(0.2),
        priorMutationRate(0.0001),
        uniqTreshold(0),
        dataUsageRate{0, 0.1},
        sampleLoops(100000),
        errorRateEstLoops(100000),
        numUniqMuts(0),
        numCellWithMutationMin(1),
        minCoverage(1),
        minCoverageAcrossCells(0),
        numMinCoverageAcrossCells(0),
        minSupport(0),
        minFreq(0),
        minCovInControlBulk(6),
        maxSupInControlBulk(2),
        noiseScore(0),
        learnZygocity(false),
        computeMixScore(false),
        numMutPlacements({{0,0}}),
        estimateSeqErrorRate(true),
        clamPrior({{200,10000}}),
        meanFilter(0.25),
        minCovNormalCell(5),
        maxNumberNormalCellMutated(0),
        useNormalCellInTree(false)
    {};
};

template <typename TTreeType>
typename Config<TTreeType>::TGraph & 
Config<TTreeType>::getTree(){return this->tree.first;}

template<typename TTreeType>
typename Config<TTreeType>::TGraph const & 
Config<TTreeType>::getTree() const {return this->tree.first;}

template<typename TTreeType>
void 
Config<TTreeType>::setTree(Config::TGraph const & newTree){this->tree.first = newTree;}

template<typename TTreeType>
void 
Config<TTreeType>::updateParamsCounter()
{
    this->paramsCounter.wildAlpha.push_back(this->getParam(Config::wildOverDis) * this->getParam(Config::wildMean));
    this->paramsCounter.wildBeta.push_back(this->getParam(Config::wildOverDis) - this->getParam(Config::wildMean) * this->getParam(Config::wildOverDis));
    this->paramsCounter.mutAlpha.push_back(this->getParam(Config::mutationOverDis) * this->getParam(Config::mutationMean));
    this->paramsCounter.mutBeta.push_back(this->getParam(Config::mutationOverDis) - (0.5 -this->getParam(Config::wildMean) ) * this->getParam(Config::mutationOverDis));
    this->paramsCounter.mu.push_back(this->getParam(Config::mu));
    this->paramsCounter.nu.push_back(this->getParam(Config::nu));
    this->paramsCounter.lambda.push_back(this->getParam(Config::lambda));
}

template<typename TTreeType>
typename Config<TTreeType>::TGraph & 
Config<TTreeType>::getTmpTree(){return this->tree.second;}

template<typename TTreeType>
void 
Config<TTreeType>::setTmpTree(Config::TGraph const & newTree){this->tree.second = newTree;}

template<typename TTreeType>
unsigned 
Config<TTreeType>::getNumSamples() const {return this->numSamples;}

template<typename TTreeType>
void 
Config<TTreeType>::setNumSamples(unsigned newNumSamples){this->numSamples = newNumSamples;}

template<typename TTreeType>
typename Config<TTreeType>::ParamType 
Config<TTreeType>::getParamToOptimize() const {return this->paramToOptimize;}

template<typename TTreeType>
void 
Config<TTreeType>::setParamToOptimize(Config<TTreeType>::ParamType newParamToOptimize){this->paramToOptimize = newParamToOptimize;}

template<typename TTreeType>
unsigned 
Config<TTreeType>::getNumMutations() const 
{
    assert(this->getData().size() > 0);    
    return this->getData()[0].size();
}

template<typename TTreeType>
unsigned 
Config<TTreeType>::getNumAttachments() const 
{
    assert(this->getData().size() > 0);    
    return this->getData()[0].size();
}

template<typename TTreeType>
typename Config<TTreeType>::TSampleNodes & 
Config<TTreeType>::getSampleNodes(){return this->sampleNodes;}

template<typename TTreeType>
typename Config<TTreeType>::TSampleNodes const & 
Config<TTreeType>::getSampleNodes() const {return this->sampleNodes;}

template<typename TTreeType>
typename Config<TTreeType>::TLogScores & 
Config<TTreeType>::getLogScores(){return this->logScores.first;}

template<typename TTreeType>
typename Config<TTreeType>::TLogScores const & 
Config<TTreeType>::getLogScores() const {return this->logScores.first;}

template<typename TTreeType>
void 
Config<TTreeType>::setLogScores(Config<TTreeType>::TLogScores newLogScores){this->logScores.first = newLogScores;}

template<typename TTreeType>
typename Config<TTreeType>::TLogScores & 
Config<TTreeType>::getTmpLogScores(){return this->logScores.second;}

template<typename TTreeType>
void 
Config<TTreeType>::setTmpLogScores(Config<TTreeType>::TLogScores & newTmpLogScores){this->logScores.second = newTmpLogScores;}

template<typename TTreeType>
typename Config<TTreeType>::TData & 
Config<TTreeType>::getData(){return this->data;}

template<typename TTreeType>
typename Config<TTreeType>::TData const & 
Config<TTreeType>::getData() const {return this->data;}

template<typename TTreeType>
typename Config<TTreeType>::TData & 
Config<TTreeType>::getCompleteData(){return this->completeData;}

template<typename TTreeType>
typename Config<TTreeType>::TData const & 
Config<TTreeType>::getCompleteData() const {return this->completeData;}

template<typename TTreeType>
typename Config<TTreeType>::TAttachmentScores & 
Config<TTreeType>::getTmpAttachmentScore(){return this->_tmpAttachmentScore;}

template<typename TTreeType>
typename Config<TTreeType>::TAttachmentScores const & 
Config<TTreeType>::getTmpAttachmentScore() const {return this->_tmpAttachmentScore;}

template<typename TTreeType>
double 
Config<TTreeType>::getParam(Config<TTreeType>::ParamType param)
{
    if (param == this->mutationMean)
        return 0.5 - (2.0/3.0 * this->getParam(this->wildMean));
    return std::get<0>(this->params[param]);
}
template<typename TTreeType>
double 
Config<TTreeType>::getParam(Config<TTreeType>::ParamType param) const
{
    if (param == this->mutationMean)
        return 0.5 - (2.0/3.0 * this->getParam(this->wildMean));
    return std::get<0>(this->params[param]);
}

template<typename TTreeType>
void 
Config<TTreeType>::setParam(Config<TTreeType>::ParamType param, double newParam){std::get<0>(this->params[param]) = newParam;}

template<typename TTreeType>
double 
Config<TTreeType>::getTmpParam(Config<TTreeType>::ParamType param){return std::get<1>(this->params[param]);}

template<typename TTreeType>
void 
Config<TTreeType>::setTmpParam(Config<TTreeType>::ParamType param, double newParam){std::get<1>(this->params[param]) = newParam;}

template<typename TTreeType>
double 
Config<TTreeType>::getSDParam(Config<TTreeType>::ParamType param){return std::get<0>(this->learningParams[param]);}

template<typename TTreeType>
void 
Config<TTreeType>::setSDParam(Config<TTreeType>::ParamType param, double newParam){std::get<0>(this->learningParams[param]) = newParam;}

template<typename TTreeType>
unsigned 
Config<TTreeType>::getSDCountParam(Config<TTreeType>::ParamType param){return std::get<1>(this->learningParams[param]);}

template<typename TTreeType>
void 
Config<TTreeType>::setSDCountParam(Config<TTreeType>::ParamType param, unsigned newParam){std::get<1>(this->learningParams[param]) = newParam;}

template<typename TTreeType>
unsigned 
Config<TTreeType>::getSDTrialsParam(Config<TTreeType>::ParamType param){return std::get<2>(this->learningParams[param]);}

template<typename TTreeType>
void 
Config<TTreeType>::setSDTrialsParam(Config<TTreeType>::ParamType param, unsigned newParam){std::get<2>(this->learningParams[param]) = newParam;}

template<typename TTreeType>
unsigned 
Config<TTreeType>::getMoveTyp(){return this->moveType;}

template<typename TTreeType>
void 
Config<TTreeType>::setMoveTyp(unsigned newMoveType){this->moveType = newMoveType;}

template <typename TTreeType>
std::default_random_engine & 
Config<TTreeType>::getGenerator(){return this->generator;}

// This function is used to reset a newly drawn parameter if the proposed tree is rejected
template <typename TTreeType>
void 
Config<TTreeType>::resetParameters()
{
    switch (this->getParamToOptimize()) 
    {
        case(wildOverDis) :
        {
            this->setParam(wildOverDis, this->getTmpParam(wildOverDis));
            break;
        }
        case(mutationOverDis) :
        {
            this->setParam(mutationOverDis, this->getTmpParam(mutationOverDis));
            break;
        }
        case(wildMean) :
        {
            this->setParam(wildMean, this->getTmpParam(wildMean));
            break;
        }
        case(mu) :
        {
            this->setParam(mu, this->getTmpParam(mu));
            break;
        }
        case(nu) :
        {
            this->setParam(nu, this->getTmpParam(nu));
            break;
        }
        case(lambda) :
        {
            this->setParam(lambda, this->getTmpParam(lambda));
            break;
        }
        case(mutationMean) :
        {
            assert(false);
            break;
        }
    }
    this->setLogScores(this->getTmpLogScores());
};

template <typename TTreeType>
void 
Config<TTreeType>::initMutInSampleCounter()
{
    std::cout << this->getNumSamples() << " " << this->getCompleteData()[0].size() << std::endl;
    
    this->mutInSampleCounter.resize(this->getNumSamples());
    for (size_t i = 0; i < this->mutInSampleCounter.size(); ++i)
    {
        this->mutInSampleCounter[i].resize(this->getCompleteData()[0].size());
        for (unsigned j = 0; j < this->mutInSampleCounter[i].size(); ++j)
        {
            this->mutInSampleCounter[i][j].setMinusInfinity();
        }
    }
}

// This function is used to increase the number of mutations used for the tree inference over time
template <typename TTreeType>
bool
Config<TTreeType>::updateContainers(unsigned currentLoop)
{
    if ((currentLoop < this->loops && // the current iteration is still smaller than the total number of iterations
        static_cast<double>(currentLoop) / static_cast<double>(this->loops) > std::get<0>(this->dataUsageRate))) // ||
        //currentLoop == this->loops + this->sampleLoops - 1)
    {   unsigned newDataSize = -1;
        if (currentLoop == this->loops + this-> sampleLoops - 1)
        {
            newDataSize = this->getCompleteData()[0].size();
        }
        else
        {
            std::cout << "numUniqMuts: " << this->numUniqMuts << std::endl;
            std::get<0>(this->dataUsageRate) += std::get<1>(this->dataUsageRate); // increase the data usage rate by a pre defined fraction
            std::cout << "dataUsage<0>: " << std::get<0>(this->dataUsageRate) << " " << std::get<1>(this->dataUsageRate) << std::endl;
            std::cout << this->getCompleteData()[0].size() << " " <<  this->numUniqMuts << std::endl;
            newDataSize = std::round(std::get<0>(this->dataUsageRate) * (this->getCompleteData()[0].size() - this->numUniqMuts));
        }
        std::cout << "newDataSize: " << newDataSize << " " << std::get<0>(this->dataUsageRate) * (this->getCompleteData()[0].size() - this->numUniqMuts) << std::endl;
        size_t j = this->getNumMutations();
        for (; j < newDataSize; ++j)
        {
            for (size_t i = 0; i < this->getNumSamples(); ++i)
            {
                this->getData()[i].push_back(this->getCompleteData()[i][j]);
            }
        }


        // resize the log score vector
        this->getLogScores().resizeNumCells(this->getNumSamples());
        this->getTmpLogScores().resizeNumCells(this->getNumSamples());
        this->getLogScores().resizeNumMuts(this->getNumMutations());

        this->getTmpAttachmentScore().resize(2 * this->getNumSamples() - 1);
	
        computeLogScoresOP(*this);
        this->setTmpLogScores(this->getLogScores());


        return true;
    }
    return false;
};

template <typename TTreeType>
void
Config<TTreeType>::printParameters(){

    std::cout << "num Samples:\t" << this->getNumSamples() << std::endl;
    std::cout << "total # mut:\t" << this->getCompleteData()[0].size() << "\tcurrently used:\t" << this->getNumMutations() << std::endl;
    std::cout << "normal     - freq:    " << this->getParam(Config::wildMean) <<      " tmp:  " << this->getTmpParam(Config::wildMean) << " SD: " << this->getSDParam(Config::wildMean) << " count: " << this->getSDCountParam(Config::wildMean)  << " trails: " << this->getSDTrialsParam(Config::wildMean) << std::endl; 
    std::cout << "normal     - overDis: " << this->getParam(Config::wildOverDis) << " tmp:  " << this->getTmpParam(Config::wildOverDis)<< " SD: " << this->getSDParam(Config::wildOverDis) << " count: " << this->getSDCountParam(Config::wildOverDis) << " trails: " << this->getSDTrialsParam(Config::wildOverDis) << std::endl;    
    std::cout << "normla     - alpha:   " << this->getParam(Config::wildOverDis) * this->getParam(Config::wildMean) << " beta: " << this->getParam(Config::wildOverDis) - this->getParam(Config::wildMean) * this->getParam(Config::wildOverDis) << std::endl;
    std::cout << "mutation   - overDis: " << this->getParam(Config::mutationOverDis) << " tmp:  " << this->getTmpParam(Config::mutationOverDis) << " SD: " << this->getSDParam(Config::mutationOverDis) << " count: " << this->getSDCountParam(Config::mutationOverDis) << " trails: " << this->getSDTrialsParam(Config::mutationOverDis) << std::endl;
    std::cout << "mutation   - alpha: " << this->getParam(Config::mutationOverDis) * this->getParam(Config::mutationMean) << " beta: " << this->getParam(Config::mutationOverDis) - (0.5 -this->getParam(Config::wildMean) ) * this->getParam(Config::mutationOverDis) << std::endl;
    std::cout << "drop: " << this->getParam(Config::mu) << " SD: " << this->getSDParam(Config::mu) << " count: " << this->getSDCountParam(Config::mu) << " trails: " << this->getSDTrialsParam(Config::mu) << std::endl;
    std::cout << "lambda: " << this->getParam(Config::nu) << " SD: " << this->getSDParam(Config::nu) << " count: " << this->getSDCountParam(Config::nu) << " trails: " << this->getSDTrialsParam(Config::nu) << std::endl;
    //std::cout << "lambda: " << this->getParam(Config::lambda) << " SD: " << this->getSDParam(Config::lambda) << " count: " << this->getSDCountParam(Config::lambda) << " trails: " << this->getSDTrialsParam(Config::lambda) << std::endl;

};


class ExtractNodesBFSVisitor : public boost::default_bfs_visitor
{
    boost::dynamic_bitset<> & bitSet;
    unsigned & numElements;

public:

    ExtractNodesBFSVisitor(boost::dynamic_bitset<> & bitSet_, unsigned &numElements_) : 
        bitSet(bitSet_),
        numElements(numElements_)
    {}

    template < typename Vertex, typename Graph >
    void discover_vertex(Vertex v, const Graph & g)
    {
        (void)g;
        this->bitSet[v] = true;
        ++this->numElements;
    }
};


class SimplifyTreeDFSVisitor : public boost::default_dfs_visitor
{
    typedef boost::adjacency_list<boost::vecS, boost::vecS, 
                            boost::bidirectionalS, Vertex<SimpleTree>> TGraph;
    typedef std::vector<std::vector<unsigned> > TMutationsOfNodes;
    typedef std::stack<unsigned> TStack;

    TGraph & newGraph;
    TMutationsOfNodes const & mutationsOfNodes;
    TStack & processedVertices;

public:

    SimplifyTreeDFSVisitor(TGraph & graph_, 
            TMutationsOfNodes const & mutationsOfNodes_, 
            TStack & processedVertices_) :
        newGraph(graph_),
        mutationsOfNodes(mutationsOfNodes_),
        processedVertices(processedVertices_)
    {}
    
    template < typename Vertex, typename Graph >
    void discover_vertex(Vertex v, const Graph & g)
    {
        (void)g;

        if (v == (num_vertices(g) - 1))
        {
            return;
        }

        // Add vertex if
        //      vertex is root
        //      vertex contains mutation
        //      vertex is sample
        if (source(*in_edges(v, g).first, g) == num_vertices(g) - 1 ||
                this->mutationsOfNodes[v].size() > 0 || 
                g[v].sample != -1)
        {
            unsigned newVertex = add_vertex(newGraph);
            newGraph[newVertex].mutations = this->mutationsOfNodes[v];
            if (g[v].sample != -1)
            {
                newGraph[newVertex].sample = g[v].sample;
            }

            // add edge if current vertex is not the root
            if (source(*in_edges(v, g).first, g) != num_vertices(g) - 1)
            {
                add_edge(this->processedVertices.top(), newVertex, this->newGraph);
            }
            this->processedVertices.push(newVertex);

        }
    }

    template <typename TVertex >
    void finish_vertex(TVertex v, boost::adjacency_list<boost::vecS, 
                                  boost::vecS, 
                                  boost::bidirectionalS, 
                                  Vertex<SampleTree>> const & g) const
    {
        if (v == (num_vertices(g) - 1))
        {
            return;
        }

        // Remove added vertex if
        //      vertex is root
        //      vertex contains mutation
        //      vertex is sample
        if (source(*in_edges(v, g).first, g) == num_vertices(g) - 1 ||
                this->mutationsOfNodes[v].size() > 0 || 
                g[v].sample != -1)
        {
            this->processedVertices.pop();
        }
    }
};



class my_label_writer {
    public:

    my_label_writer(boost::adjacency_list<boost::vecS,boost::vecS, boost::bidirectionalS, Vertex<SimpleTree>> const & simpleTree_,
            std::vector<std::tuple<std::string, unsigned, char, char>> const & indexToPosition_,
            std::vector<std::string> const & cellNames_) : 
        simpleTree(simpleTree_),
        indexToPosition(indexToPosition_),
        cellNames(cellNames_)
    {}

    template <class VertexOrEdge>
    void operator()(std::ostream& out, const VertexOrEdge& v) const {

        if (simpleTree[v].sample == -1)
        {
            out << "[label=\"";
        }
        else
        {
            out << "[shape=box,label=\"" << cellNames[simpleTree[v].sample] << " ";
        }
        for (unsigned i = 0; i < simpleTree[v].mutations.size(); ++i)
        {
            out << std::get<0>(this->indexToPosition[simpleTree[v].mutations[i]]) << "_" << std::to_string(std::get<1>(this->indexToPosition[simpleTree[v].mutations[i]])) <<  ":";
        }
        out << "\"]";
    }

    private:
    boost::adjacency_list<boost::vecS,boost::vecS, boost::bidirectionalS, Vertex<SimpleTree>> const & simpleTree;
    std::vector<std::tuple<std::string, unsigned, char, char>> const & indexToPosition;
    std::vector<std::string> const & cellNames;
};

// this prints the graph as it is currently used
class my_label_writer_debug {
    public:

    my_label_writer_debug(boost::adjacency_list<boost::vecS,boost::vecS, boost::bidirectionalS, Vertex<SampleTree>> const & sampleTree_,
            std::vector<std::tuple<std::string, unsigned, char, char>> const & indexToPosition_) : 
        sampleTree(sampleTree_)
    {}

    template <class VertexOrEdge>
    void operator()(std::ostream& out, const VertexOrEdge& v) const {

        if (sampleTree[v].sample == -1)
        {
            out << "[label=\"" << v;
        }
        else
        {
            out << "[shape=box,label=\"" << v;
        }
        out << "\"]";
    }

    private:
    boost::adjacency_list<boost::vecS,boost::vecS, boost::bidirectionalS, Vertex<SampleTree>> const & sampleTree;
};

template <typename TTreeType>
double
computeWildLogScore(Config<TTreeType> & config, double altCount, double coverage);
template <typename TTreeType>
double
computeMutLogScore(Config<TTreeType> & config, double altCount, double coverage);

std::array<double, 3> getStats(std::vector<double> & values)
{
    std::array<double, 3> result;
    std::sort(values.begin(), values.end());
    result[0] = values[values.size()/2];

    double sum = 0;
    for (unsigned i = 0; i < values.size(); ++i)
    {
        sum += values[i];
    }
    result[1] = sum / static_cast<double>(values.size());

    double variance = 0;
    for (unsigned i = 0; i < values.size(); ++i)
    {
        variance += std::pow(values[i]-result[1], 2);
    }
    result[2] = std::sqrt(variance/static_cast<double>(values.size()));

    return result;
}

template<typename TTreeType>
void
writeParameters(Config<TTreeType> & config, std::string const & fileName)
{
    std::ofstream outFile;
    outFile.open(fileName);
	outFile << "background frequency:\t" << config.getParam(Config<TTreeType>::wildMean) << std::endl;
    std::array<double, 3> stats = getStats(config.paramsCounter.wildAlpha);
	outFile << "wildAlpha:\t" << stats[0]<< "\t" << stats[1] << "\t" << stats[2] << std::endl;
    stats = getStats(config.paramsCounter.wildBeta);
	outFile << "wildBeta:\t" << stats[0]<< "\t" << stats[1] << "\t" << stats[2] << std::endl;
	outFile << "average mutation frequency:\t" << config.getParam(Config<TTreeType>::mutationMean) << std::endl;
    stats = getStats(config.paramsCounter.mutAlpha);
	outFile << "altAlpha:\t" << stats[0]<< "\t" << stats[1] << "\t" << stats[2] << std::endl;
    stats = getStats(config.paramsCounter.mutBeta);
	outFile << "altBeta:\t" << stats[0]<< "\t" << stats[1] << "\t" << stats[2] << std::endl;
    stats = getStats(config.paramsCounter.mu);
	outFile << "mu:\t" << stats[0]<< "\t" << stats[1] << "\t" << stats[2] << std::endl;
    stats = getStats(config.paramsCounter.nu);
	outFile << "nu:\t" << stats[0]<< "\t" << stats[1] << "\t" << stats[2] << std::endl;
    stats = getStats(config.paramsCounter.lambda);
	//outFile << "lambda:\t" << stats[0]<< "\t" << stats[1] << "\t" << stats[2] << std::endl;
}

#endif
