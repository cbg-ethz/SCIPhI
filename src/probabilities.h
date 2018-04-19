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
#ifndef PROBABILITIES_H
#define PROBABILITIES_H

#include <vector>
#include <stdlib.h>
#include <cmath>
#include <string.h>
#include <float.h>
#include <iostream>
#include <random>
#include <fstream>
#include <sstream>
#include <limits.h>
#include "rand.h"
#include "sciphi_config.h"
#include "trees.h"

#include <dlib/global_optimization.h>

#include <boost/math/special_functions/digamma.hpp>

inline 
double 
logBetaBinCountsTerm(double sup, double cov)
{
    return std::lgamma(cov + 1.0) -  std::lgamma(sup + 1.0) - std::lgamma(cov - sup + 1.0);
}

inline 
double 
logBetaBinMixedTerm(double sup, double cov, double mean, double overDis)
{
    return std::lgamma(sup + mean * overDis) + std::lgamma(cov - sup + overDis * (1.0 - mean)) - std::lgamma(cov + overDis);
}

inline 
double 
logBetaBinParamsTerm(double mean, double overDis)
{
    return std::lgamma(overDis) - std::lgamma(mean * overDis) - std::lgamma(overDis * (1.0 - mean));
}

inline 
double 
logBetaBinPDF(double sup, double cov, double mean, double overDis)
{
    if (cov == 0)
        return 0;

    return logBetaBinCountsTerm(sup, cov) +
        logBetaBinMixedTerm(sup, cov, mean, overDis) + 
        logBetaBinParamsTerm(mean, overDis);
}

inline 
double 
logBetaBinPDFOP(double sup, double cov, double mean, double overDis)
{
    if (cov == 0)
        return 0;

    return //logBetaBinCountsTerm(sup, cov) + // we do not need the counts term as we are only interested in the score between iterations
        logBetaBinMixedTerm(sup, cov, mean, overDis) + 
        logBetaBinParamsTerm(mean, overDis);
}

inline 
double 
betaBinPDF(double sup, double cov, double mean, double overDis)
{
	return std::exp(logBetaBinPDF(sup, cov, mean, overDis));
}

// Add two values in real space by first exponentiating
inline
double
addLogProb(double x, double y)
{
    double maxScore;
    double minScore;
    if (x > y)
    {
        maxScore = x;
        minScore = y;
    }
    else
    {
        maxScore = y;
        minScore = x;
    }
    return std::log(1.0 + std::exp(minScore-maxScore)) + maxScore;
}

// Subtract two values in real space by first exponentiating
inline
double
subLogProb(double x, double y)
{
    double maxScore;
    double minScore;
    if (x > y)
    {
        maxScore = x;
        minScore = y;
    }
    else
    {
        maxScore = y;
        minScore = x;
    }
    return log(1.0 - exp(minScore-maxScore)) + maxScore;
}

// Add two weighted values in real space by first exponentiating
inline
double
addLogProbWeight(double x, double y, double nu) // = (1.0 - nu) * x + nu * y
{
    if (nu == 0.0)
    {
        return x;
    }
    else if (nu == 1.0)
    {
        return y;
    }

    if (x > y)
    {
        return std::log((1.0 - nu) + nu * std::exp(y-x)) + x;
    }
    return std::log(nu + (1.0 - nu) * std::exp(x-y)) + y;
}

template <typename TTreeType>
double
computeWildLogScore(Config<TTreeType> const & config, double altCount, double coverage)
{
    return logBetaBinPDF(altCount, 
        coverage, 
        config.getParam(Config<TTreeType>::wildMean), 
        config.getParam(Config<TTreeType>::wildOverDis));
}

template <typename TTreeType>
double
computeWildLogScoreOP(Config<TTreeType> const & config, double altCount, double coverage)
{
    return logBetaBinPDFOP(altCount, 
        coverage, 
        config.getParam(Config<TTreeType>::wildMean), 
        config.getParam(Config<TTreeType>::wildOverDis));
}

template <typename TTreeType>
double
computeRawMutLogScore(Config<TTreeType> const & config, double altCount, double coverage)
{
    return logBetaBinPDF(altCount, 
                coverage, 
                config.getParam(Config<TTreeType>::mutationMean), 
                config.getParam(Config<TTreeType>::mutationOverDis));
}

template <typename TTreeType>
double
computeRawMutLogScoreOP(Config<TTreeType> & config, double altCount, double coverage)
{
    return logBetaBinPDFOP(altCount, 
                coverage, 
                config.getParam(Config<TTreeType>::mutationMean), 
                config.getParam(Config<TTreeType>::mutationOverDis));
}


template <typename TTreeType>
double
computeHomoLogScoreOP(Config<TTreeType> & config, double altCount, double coverage)
{
    return logBetaBinPDFOP(coverage-altCount, 
        coverage, 
        config.getParam(Config<TTreeType>::wildMean), 
        config.getParam(Config<TTreeType>::wildOverDis));
}

// This function combines the wild type beta-bin, the hetero beta-bin and the homo mutand beta-bin
template <typename TTreeType>
double
computeHeteroLogScoreOP(Config<TTreeType> & config, double altCount, double coverage)
{
    double mu = config.getParam(Config<TTreeType>::mu);
    double oneMinusMu = 1.0 - mu;
    double logHomo = addLogProbWeight(computeWildLogScoreOP(config, altCount, coverage), computeHomoLogScoreOP(config, altCount,  coverage), 0.5);
    return addLogProbWeight(computeRawMutLogScoreOP(config, altCount, coverage), logHomo, oneMinusMu);
}

template <typename TTreeType>
void computeWildLogScoresOP(Config<TTreeType> & config)
{
	for (unsigned int i = 0; i < config.getLogScores().size(); ++i)
	{
		for (unsigned int j = 0; j < config.getLogScores()[0].size(); ++j)
		{
            std::get<0>(config.getLogScores()[i][j]) = computeWildLogScoreOP(config, std::get<1>(config.getData()[i][j]), std::get<0>(config.getData()[i][j]));
		}
	}
}

template <typename TTreeType>
void computeLogScoresOP(Config<TTreeType> & config)
{
	for (unsigned int i = 0; i < config.getLogScores().numCells(); ++i)
	{
		for (unsigned int j = 0; j < config.getLogScores().numMuts(); ++j)
		{
            config.getLogScores().wtScore(i, j) = computeWildLogScoreOP(config, std::get<1>(config.getData()[i][j]), std::get<0>(config.getData()[i][j]));
            config.getLogScores().hetScore(i, j) = computeHeteroLogScoreOP(config, std::get<1>(config.getData()[i][j]), std::get<0>(config.getData()[i][j]));
            config.getLogScores().homScore(i, j) = computeHomoLogScoreOP(config, std::get<1>(config.getData()[i][j]), std::get<0>(config.getData()[i][j]));
		}
	}
}

inline
void
computeNoiseScore(Config<SampleTree> & config)
{

    long double mean = config.getParam(Config<SampleTree>::wildMean);
    long double overDis = config.getParam(Config<SampleTree>::wildOverDis);
    long double alpha = mean * overDis;
    long double beta = overDis * (1.0 - mean);

    config.noiseScore = config.noiseCounts.numPos * logBetaBinParamsTerm(mean, overDis);
    for (unsigned i = 0; i < config.noiseCounts.sup.size(); ++i)
    {
        config.noiseScore += config.noiseCounts.sup[i].second * std::lgamma(config.noiseCounts.sup[i].first + alpha);
    }
    for (unsigned i = 0; i < config.noiseCounts.covMinusSup.size(); ++i)
    {
        config.noiseScore += config.noiseCounts.covMinusSup[i].second * std::lgamma(config.noiseCounts.covMinusSup[i].first + beta); 
    }   
    for (unsigned i = 0; i < config.noiseCounts.cov.size(); ++i)
    {
        config.noiseScore -= config.noiseCounts.cov[i].second * std::lgamma(config.noiseCounts.cov[i].first + overDis); 
    }
}

// This class is used to pass the score from the root towards the leaves
class PassScoreToChildrenBFSVisitor : public boost::default_bfs_visitor
{
    Config<SampleTree> & config;
    Config<SampleTree>::TAttachmentScores & attachmentScores;
    Config<SampleTree>::TAttachmentScores::TAttachmentScore & sumScore;

public:

    PassScoreToChildrenBFSVisitor(Config<SampleTree> & config_, 
            Config<SampleTree>::TAttachmentScores & attachmentScores_, 
            Config<SampleTree>::TAttachmentScores::TAttachmentScore & sumScore_, 
            unsigned attachment_) : 
        config(config_),
        attachmentScores(attachmentScores_),
        sumScore(sumScore_)
    {}

    template <typename TVertex >
    void discover_vertex(TVertex v, boost::adjacency_list<boost::vecS, 
                                  boost::vecS, 
                                  boost::bidirectionalS, 
                                  Vertex<SampleTree>> const & g) const
    {
        (void)g;
        if (v == num_vertices(config.getTree()) - 1)
        {
            return;
        }

        unsigned parentNode = source(*in_edges(v, g).first, g);
        if (parentNode == num_vertices(config.getTree()) - 1)
        {
            sumScore.addInRealSpace(attachmentScores[v]);
            return;
        }

        if (g[v].sample == -1)
        {
            sumScore.addInRealSpace(attachmentScores[v]);
            attachmentScores[v].homScore() = addLogProb(attachmentScores[v].homScore(), attachmentScores[parentNode].homScore());
            attachmentScores[v].hetScore() = addLogProb(attachmentScores[v].hetScore(), attachmentScores[parentNode].hetScore());
            if(config.computeMixScore)
            {
                if (isnan(attachmentScores[v].mixWildScore()))
                {
                    attachmentScores[v].mixWildScore() = attachmentScores[parentNode].mixWildScore();
                    attachmentScores[v].mixHomScore() = attachmentScores[parentNode].mixHomScore();
                }
                else
                {
                    attachmentScores[v].mixWildScore() = addLogProb(attachmentScores[v].mixWildScore(), attachmentScores[parentNode].mixWildScore());
                    attachmentScores[v].mixHomScore() = addLogProb(attachmentScores[v].mixHomScore(), attachmentScores[parentNode].mixHomScore());
                }
            }
        }
        else
        {
            sumScore.hetScore() = addLogProb(sumScore.hetScore(), attachmentScores[v].hetScore());
            attachmentScores[v].homScore() = attachmentScores[parentNode].homScore();
            attachmentScores[v].hetScore() = addLogProb(attachmentScores[v].hetScore(), attachmentScores[parentNode].hetScore());
            
            if(config.computeMixScore)
            {
                unsigned grandParent = source(*in_edges(parentNode, g).first, g);
                if (grandParent == num_vertices(config.getTree()) - 1)
                {
                    attachmentScores[v].mixWildScore() = -INFINITY;
                    attachmentScores[v].mixHomScore() = -INFINITY;
                }
                else
                {
                    attachmentScores[v].mixWildScore() = attachmentScores[parentNode].mixWildScore();
                    attachmentScores[v].mixHomScore() = attachmentScores[parentNode].mixHomScore();
                }
            }
        }
    }
};


// experimental phase
class ComputeLostScoreDFSVisitor : public boost::default_dfs_visitor
{
    Config<SampleTree> & config;
    Config<SampleTree>::TAttachmentScores const & attachmentScores;
    std::vector<double> & lostScores;
    double & leafScore;

public:

    ComputeLostScoreDFSVisitor(Config<SampleTree> & config_, 
            Config<SampleTree>::TAttachmentScores const & attachmentScores_,
            std::vector<double> & lostScores_,
            double & leafScore_,
            unsigned attachment_) : 
        config(config_),
        attachmentScores(attachmentScores_),
        lostScores(lostScores_),
        leafScore(leafScore_)
    {}

    template <typename TVertex >
    void discover_vertex(TVertex v, boost::adjacency_list<boost::vecS, 
                                  boost::vecS, 
                                  boost::bidirectionalS, 
                                  Vertex<SampleTree>> const & g) //const
    {
        // this is the artificial node which needs to be excludede
        if (v == num_vertices(config.getTree()) - 1)
        {
            return;
        }

        unsigned parentNode = source(*in_edges(v, g).first, g);

        // this is the root node which needs to be excludede
        if (parentNode == num_vertices(config.getTree()) - 1)
        {
            return;
        }
        
        if (g[v].sample == -1) // inner node
        {
            // if the parent node is the root then there should be no summation of the nodes
            if (source(*in_edges(parentNode, g).first, g) == num_vertices(config.getTree()) - 1)
            {
                //std::cout << v << " " << lostScores[v] << " " << attachmentScores[parentNode].hetScore() - attachmentScores[v].hetScore() << std::endl;
                lostScores[v] = attachmentScores[parentNode].hetScore() - attachmentScores[v].hetScore();
                leafScore = addLogProb(leafScore, lostScores[v]);
            }
            else
            {
                //std::cout << v << " " << parentNode << " "  << getSibling(g, v) << " " << std::exp(lostScores[v]) << " " << std::exp(attachmentScores[parentNode].hetScore() - attachmentScores[v].hetScore()) << " " << std::exp(lostScores[parentNode]) << " " << std::exp(attachmentScores[getSibling(g, v)].hetScore()) << " " << std::exp(lostScores[parentNode] + attachmentScores[getSibling(g, v)].hetScore());
                lostScores[v] = attachmentScores[parentNode].hetScore() - attachmentScores[v].hetScore();
                lostScores[v] = addLogProb(lostScores[v], lostScores[parentNode] + attachmentScores[getSibling(g, v)].hetScore());
                //std::cout << " " <<  std::exp(lostScores[v]) << std::endl;
                leafScore = addLogProb(leafScore, lostScores[v]);
            }
        }
        else // leaf
        {
            if (source(*in_edges(parentNode, g).first, g) == num_vertices(config.getTree()) - 1)
            {
                //std::cout << v << " " << lostScores[v] << " " << attachmentScores[parentNode].hetScore() - attachmentScores[v].hetScore() << std::endl;
                lostScores[v] = -INFINITY;
            }
            else
            {
                lostScores[v] = leafScore;
            }
        }
    }

    template <typename TVertex >
    void finish_vertex(TVertex v, boost::adjacency_list<boost::vecS, 
                                  boost::vecS, 
                                  boost::bidirectionalS, 
                                  Vertex<SampleTree>> const & g) 
    {
        if (v == num_vertices(config.getTree()) - 1)
        {
            return;
        }

        unsigned parentNode = source(*in_edges(v, g).first, g);
        if (parentNode == num_vertices(config.getTree()) - 1)
        {
            return;
        }
        if (g[v].sample == -1) // inner node
        {
            // if the parent node is the root then there should be no summation of the nodes
            if (source(*in_edges(parentNode, g).first, g) == num_vertices(config.getTree()) - 1)
            {
                leafScore = subLogProb(leafScore, lostScores[v]);
            }
            else
            {
                leafScore = subLogProb(leafScore, lostScores[v]);
            }
        }

    }
};

/*
 * This function is used to optimize mean and overdispersion
 * for a given loci
 */
struct OptimizeBetaBinMeanOverDis
{
    std::vector<std::pair<unsigned, unsigned>> const & counts;

    OptimizeBetaBinMeanOverDis(
            std::vector<std::pair<unsigned, unsigned>> const & counts_) :
        counts(counts_)
    {};

    double operator()(const dlib::matrix<double,0,1> & x) const
    {
        double result = 0;
        for (size_t cell = 0; cell < this->counts.size(); ++cell)
        {
            result += logBetaBinPDF(this->counts[cell].first, this->counts[cell].second, x(0), x(1));
        }
        return result;
    };
};

struct OptimizeBetaBinMeanOverDisDerivates
{
    std::vector<std::pair<unsigned, unsigned>> const & counts;

    OptimizeBetaBinMeanOverDisDerivates(
            std::vector<std::pair<unsigned, unsigned>> const & counts_) :
        counts(counts_)
    {}

    dlib::matrix<double> operator()(const dlib::matrix<double,0,1> & x) const
    {
        double mean = x(0);
        double overDis = x(1);
        dlib::matrix<double,0,1> res = {0,0};

        double temp = 0;
        unsigned counter = 0;
        for (size_t cell = 0; cell < this->counts.size(); ++cell)
        {
            unsigned k = this->counts[cell].first;
            unsigned n = this->counts[cell].second;
            temp += overDis *  boost::math::digamma(k + mean * overDis)
                - overDis *  boost::math::digamma(n - k + overDis - overDis * mean);
            ++counter;
        }
        res(0) = counter * (-overDis *  boost::math::digamma(mean * overDis) + overDis *  boost::math::digamma(overDis - overDis * mean)) + temp;

        temp = 0;
        for (size_t cell = 0; cell < this->counts.size(); ++cell)
        {
            unsigned k = this->counts[cell].first;
            unsigned n = this->counts[cell].second;
            temp += mean * boost::math::digamma(k + mean * overDis) +
                (1.0 - mean) * boost::math::digamma(n - k + overDis - overDis * mean) -
                 boost::math::digamma(n + overDis);
        }
        res(1) = counter * (boost::math::digamma(overDis) - mean * boost::math::digamma(mean * overDis) - (1.0 - mean) * boost::math::digamma(overDis - overDis * mean)) + temp;

        return res;
    }
};

struct OptimizeBetaBinOverDis
{
    std::vector<std::pair<unsigned, unsigned>> const & counts;
    double meanFilter;

    OptimizeBetaBinOverDis(
            std::vector<std::pair<unsigned, unsigned>> const & counts_, 
            double meanFilter_) : 
        counts(counts_),
        meanFilter(meanFilter_)
    {};

    double operator()(const dlib::matrix<double,0,1> & x) const
    {
        double result = 0;
        for (size_t cell = 0; cell < this->counts.size(); ++cell)
        {
             result += logBetaBinPDF(this->counts[cell].first, this->counts[cell].second, meanFilter, x(0));
        }
        return result;
    };
};

struct OptimizeBetaBinOverDisDerivates
{
    std::vector<std::pair<unsigned, unsigned>> const & counts;
    double meanFilter;

    OptimizeBetaBinOverDisDerivates(
            std::vector<std::pair<unsigned, unsigned>> const & counts_, 
            double meanFilter_) : 
        counts(counts_),
        meanFilter(meanFilter_)
    {};

    dlib::matrix<double> operator()(const dlib::matrix<double,0,1> & x) const
    {
        double mean = this->meanFilter;
        double overDis = x(0);
        dlib::matrix<double,0,1> res = {0};

        double temp = 0;
        unsigned counter = 0;
        for (size_t cell = 0; cell < this->counts.size(); ++cell)
        {
            unsigned k = this->counts[cell].first;
            unsigned n = this->counts[cell].second;
            temp += mean * boost::math::digamma(k + mean * overDis) +
                (1.0 - mean) * boost::math::digamma(n - k + overDis - overDis * mean) -
                 boost::math::digamma(n + overDis);
            ++counter;
        }
        res(0) = counter * (boost::math::digamma(overDis) - mean * boost::math::digamma(mean * overDis) - (1.0 - mean) * boost::math::digamma(overDis - overDis * mean)) + temp;

        return res;
    }
};



#endif
