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
#ifndef SCORETREE_H_
#define SCORETREE_H_


#include <stdbool.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <sstream>
#include <float.h>
#include <math.h>
#include <cmath>
#include <queue>
#include <tuple>
#include "sciphi_config.h"
#include "trees.h"
#include "probabilities.h"
#include "attachmentScores.h"

#include <boost/math/distributions/beta.hpp>

using namespace std;

//extern double epsilon;
//std::vector<std::string> getGeneNames(string fileName, unsigned numSamples);

template <typename TTreeType>
unsigned getBestAttachmentPosition(Config<TTreeType> & config,
							  unsigned int attachment);

/*
 * This class is used to compute the attachment scores.
 * In order to do so, the scores of the child nodes are used
 * to compute the score of the current node.
 */
class ComputeScoreDFSVisitor : public boost::default_dfs_visitor
{

public:
    Config<SampleTree> & config;
    Config<SampleTree>::TAttachmentScores & attachmentScores;
    unsigned gene;
    Config<SampleTree>::TAttachmentScores::TAttachmentScore & scoreSum;

    ComputeScoreDFSVisitor(Config<SampleTree> & config_, 
            Config<SampleTree>::TAttachmentScores & attachmentScores_, 
            unsigned gene_,
            Config<SampleTree>::TAttachmentScores::TAttachmentScore & scoreSum_) :
        config(config_),
        attachmentScores(attachmentScores_),
        gene(gene_),
        scoreSum(scoreSum_)
    {}
    
    template <typename TVertex >
    void finish_vertex(TVertex v, boost::adjacency_list<boost::vecS, 
                                  boost::vecS, 
                                  boost::bidirectionalS, 
                                  Vertex<SampleTree>> const & g) const
    {
        // ignore artificial root
        if (v == num_vertices(config.getTree()) - 1)
        {
            return;
        }

        // resue container
        Config<SampleTree>::TAttachmentScores & attachmentScores = this->attachmentScores;

        // compute scores for leaf nodes
        if (g[v].sample != -1)
        {
            attachmentScores.computeLogHetScoreLeaf(v, this->config.getLogScores().wtScore(g[v].sample, gene), this->config.getLogScores().hetScore(g[v].sample, gene));
            attachmentScores.computeLogHomScoreLeaf(v, this->config.getLogScores().wtScore(g[v].sample, gene), this->config.getLogScores().homScore(g[v].sample, gene));
            if (config.computeMixScore) // experimental
            {
                attachmentScores.computeLogMixWildScoreLeaf(v);
                attachmentScores.computeLogMixHomScoreLeaf(v);
            }
            // add hetero score to overall tree score
            scoreSum.hetScore() = addLogProb(scoreSum.hetScore(), attachmentScores[v].hetScore());

            return ;
        }

        // compute score of inner node and add it to the overall tree socre
        auto it = out_edges(v,g).first;
        attachmentScores.computeLogHetScoreInnerNode(v, attachmentScores.hetScore(target(*it, g)), attachmentScores.hetScore(target(*(it+1), g)));
        attachmentScores.computeLogHomScoreInnerNode(v, attachmentScores.homScore(target(*it, g)), attachmentScores.homScore(target(*(it+1), g)));
        scoreSum.hetScore() = addLogProb(scoreSum.hetScore(), attachmentScores[v].hetScore());
        scoreSum.homScore() = addLogProb(scoreSum.homScore(), attachmentScores[v].homScore());

        if (config.computeMixScore) // experimental
        {

            unsigned leftNodeId = target(*it, g);
            unsigned rightNodeId = target(*(it + 1), g);
            bool innerNodeLeft = g[leftNodeId].sample == -1;
            bool innerNodeRight = g[rightNodeId].sample == -1;
            attachmentScores.computeLogMixWildScoreInnerNode(
                    v,
                    innerNodeLeft, 
                    attachmentScores.mixWildScore(leftNodeId),
                    attachmentScores.hetScore(leftNodeId),
                    innerNodeRight, 
                    attachmentScores.mixWildScore(rightNodeId),
                    attachmentScores.hetScore(rightNodeId));
            attachmentScores.computeLogMixHomScoreInnerNode(
                    v,
                    innerNodeLeft, 
                    attachmentScores.homScore(leftNodeId), 
                    attachmentScores.mixHomScore(leftNodeId),
                    attachmentScores.hetScore(leftNodeId),
                    innerNodeRight, 
                    attachmentScores.homScore(rightNodeId), 
                    attachmentScores.mixHomScore(rightNodeId),
                    attachmentScores.hetScore(rightNodeId));

            if(innerNodeLeft || innerNodeRight)
            {
                scoreSum.mixWildScore() = addLogProb(scoreSum.mixWildScore(), attachmentScores[v].mixWildScore());
                scoreSum.mixHomScore() = addLogProb(scoreSum.mixHomScore(), attachmentScores[v].mixHomScore());
            }
        }
        return ;
    }
};

// This function invokes the computation of the tree score
void getAllAttachmentScores(typename Config<SampleTree>::TAttachmentScores & attachmentScores,
                            Config<SampleTree>::TAttachmentScores::TAttachmentScore & scoreSum,
                            Config<SampleTree> & config,
							unsigned int attachment)
{
    scoreSum.setMinusInfinity();
    ComputeScoreDFSVisitor vis(config, attachmentScores, attachment, scoreSum);
    unsigned rootVertex = target(*out_edges(num_vertices(config.getTree()) - 1, config.getTree()).first, config.getTree());

    depth_first_search(config.getTree(), visitor(vis).root_vertex(rootVertex));
}

// this function computes the maximum likelihood score of a tree
template <typename TTreeType>
double getBestAttachmentScore(Config<TTreeType> & config,
							  unsigned int attachment)
{
    Config<SampleTree>::TAttachmentScores::TAttachmentScore scoreSum;
    scoreSum.setMinusInfinity();
	
    typename Config<TTreeType>::TAttachmentScores & attachmentScores = config.getTmpAttachmentScore();
    getAllAttachmentScores(attachmentScores, scoreSum, config, attachment);

	attachmentScores[0].computeFinalScore(
            config.getParam(Config<TTreeType>::nu),
            config.getParam(Config<TTreeType>::lambda),
            config.getNumSamples() -1,
            config.numMutPlacements[0],
            config.getTree()[0].sample != -1,
            config.computeMixScore);
	double bestScore = attachmentScores[0].finalScore();
    for (unsigned i = 1; i < attachmentScores.size(); ++i)
    {
        attachmentScores[i].computeFinalScore(
            config.getParam(Config<TTreeType>::nu),
            config.getParam(Config<TTreeType>::lambda),
            config.getNumSamples() -1,
            config.numMutPlacements[0],
            config.getTree()[i].sample != -1,
            config.computeMixScore);
        if (attachmentScores[i].finalScore() > bestScore)
        {
            bestScore = attachmentScores[i].finalScore();
        }
    }

    return bestScore;
}

// This function computes the best attachment node for a given mutation
template <typename TTreeType>
unsigned getBestAttachmentPosition(Config<TTreeType> & config,
							  unsigned int attachment){
	
    typename Config<TTreeType>::TAttachmentScores::TAttachmentScore scoreSum;
    scoreSum.setMinusInfinity();
    typename Config<TTreeType>::TAttachmentScores & attachmentScores = config.getTmpAttachmentScore();
    getAllAttachmentScores(attachmentScores, scoreSum, config, attachment);
    

    unsigned bestPos = 0;
	attachmentScores[0].computeFinalScore(
            config.getParam(Config<TTreeType>::nu),
            config.getParam(Config<TTreeType>::lambda),
            config.getNumSamples() -1,
            config.numMutPlacements[0],
            config.getTree()[0].sample != -1,
            config.computeMixScore);
	double bestScore = attachmentScores[0].finalScore();
    for (unsigned i = 1; i < attachmentScores.size(); ++i)
    {
        attachmentScores[i].computeFinalScore(
            config.getParam(Config<TTreeType>::nu),
            config.getParam(Config<TTreeType>::lambda),
            config.getNumSamples() -1,
            config.numMutPlacements[0],
            config.getTree()[i].sample != -1,
            config.computeMixScore);
        if (attachmentScores[i].finalScore() > bestScore)
        {
            bestScore = attachmentScores[i].finalScore();
            bestPos = i;
        }
    }

    return bestPos;
}

// This function puts a prior on lambda and is experimental
template <typename TTreeType>
double clamPriorLog(Config<TTreeType> & config)
{
    boost::math::beta_distribution<> betaPDF(config.clamPrior[0],config.clamPrior[1]);
    return std::log(boost::math::pdf(betaPDF, config.getParam(Config<TTreeType>::lambda)));
}

// This function compute the maximum likelihood score of a tree
inline
double maxScoreTree(Config<SampleTree> & config)
{
    double noAttachmentScore = 0;
    for (unsigned attachment = 0; attachment < config.getNumMutations(); ++ attachment)
    {
        for (unsigned i = 0; i < config.getLogScores().numCells(); ++i)
        {
            noAttachmentScore += config.getLogScores().wtScore(i, attachment);
        }
    }

  	double treeScore = 0.0;
  	for(unsigned mut = 0; mut < config.getNumMutations(); ++mut){
  		treeScore += getBestAttachmentScore(config, mut);
    }
    
    if (config.computeMixScore) // experimental
    {
	    return treeScore + noAttachmentScore + config.noiseScore + clamPriorLog(config);
    }
	return treeScore + noAttachmentScore + config.noiseScore;
}

// This function computes the final score of the tree
template <typename TTreeType>
double getSumAttachmentScore(Config<TTreeType> & config,
                             typename Config<TTreeType>::TAttachmentScores & attachmentScores,
							 unsigned int attachment)
{
    Config<SampleTree>::TAttachmentScores::TAttachmentScore scoreSum;
    getAllAttachmentScores(attachmentScores, scoreSum, config, attachment);

    scoreSum.computeFinalScore(
            config.getParam(Config<TTreeType>::nu), 
            config.getParam(Config<TTreeType>::lambda), 
            config.getNumSamples() - 1,
            config.numMutPlacements[0], 
            false,
            config.computeMixScore);

    return scoreSum.finalScore();
}

template <typename TTreeType>
double getSumAttachmentScore(Config<TTreeType> & config,
							 unsigned int attachment)
{

    typename Config<TTreeType>::TAttachmentScores & attachmentScores = config.getTmpAttachmentScore();
    return getSumAttachmentScore(config, attachmentScores, attachment);

}

/* computes the log score for the complete tree using the sumScore scheme, where likelihoods of all attachment points of a sample are added */
template <typename TTreeType>
double sumScoreTree(Config<TTreeType> & config)
{
    double noAttachmentScore = 0;
    for (unsigned attachment = 0; attachment < config.getNumMutations(); ++ attachment)
    {
        for (unsigned i = 0; i < config.getLogScores().numCells(); ++i)
        {
            noAttachmentScore += config.getLogScores().wtScore(i, attachment);
        }
    }

	double sumTreeScore = 0.0;

	for(std::size_t attachment=0; attachment < config.getNumAttachments(); attachment++)
    {
		sumTreeScore += getSumAttachmentScore(config, attachment);
	}

    if (config.computeMixScore)
    {
	    return sumTreeScore + noAttachmentScore + config.noiseScore + clamPriorLog(config);
    }

	return sumTreeScore + noAttachmentScore + config.noiseScore;
}

template <typename TTreeType>
double scoreTree(Config<TTreeType> & config)
{
	double result = -DBL_MAX;

	if(config.scoreType=='m')
    {
		result = maxScoreTree(config);
    }
	else
    {
		result = sumScoreTree(config);
    }
    
	return result;
}
#endif /* SCORETREE_H_ */
