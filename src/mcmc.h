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
#ifndef MCMC_H
#define MCMC_H

#include "sciphi_config.h"
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
#include "trees.h"
#include "rand.h"
#include "scoreTree.h"
#include "probabilities.h"

template <typename TTreeType>
void 
printLogScores(Config<TTreeType> & config)
{
	for (unsigned int i = 0; i < config.getLogScores().size(); ++i)
	{
		for (unsigned int j = 0; j < config.getLogScores()[0].size(); ++j)
		{
            std::cout << std::get<0>(config.getLogScores()[i][j]) << " | " << std::get<1>(config.getLogScores()[i][j]) << "\t";
		}
        std::cout << std::endl;
	}
}

// This function keeps the coeffients between 0 and 1
// by mirroring at 0 or 1
double adjusteCoefficient(double mixtureCoEff)
{
    while(true)
    {
        if(mixtureCoEff < 0.0)
        {
            mixtureCoEff = std::abs(mixtureCoEff);
        }
        else if (mixtureCoEff > 1.0)
        {
            mixtureCoEff = 1.0 - std::abs(1 - mixtureCoEff);
        }
        else
        {
            break;
        }
    }
    return mixtureCoEff;
}

// This function assures that the overdispersion of the heterozygous mutation
// cannot take a U form (alpha and beta are very close or bigger than 1)
double adjustMutOverDis(double mutationOverDis)
{
    if (mutationOverDis < 2.0)
    {
        return 2.0 + (2.0 - mutationOverDis);
    }
    return mutationOverDis;
}

template <typename TTreeType>
void
changeParameters(Config<TTreeType> & config)
{
    // store the old scores, in case the new parameter will be rejected
    config.setTmpLogScores(config.getLogScores());

    //choose which parameter to optimize
    if (config.computeMixScore && config.learnZygocity) 
    {
        config.setParamToOptimize((typename Config<TTreeType>::ParamType)(rand() % 6));
    }
    else if (config.learnZygocity || config.computeMixScore)
    {
        config.setParamToOptimize((typename Config<TTreeType>::ParamType)(rand() % 5));
        if (config.getParamToOptimize() == Config<TTreeType>::nu && config.computeMixScore)
        {
            config.setParamToOptimize(Config<TTreeType>::lambda);
        }
    }
    else
    {
        config.setParamToOptimize((typename Config<TTreeType>::ParamType)(rand() % 4));
    }
    typename Config<TTreeType>::ParamType const & param = config.getParamToOptimize();

    // update the standard deviaten of the chosen parameter every 1000 steps
    if (config.getSDTrialsParam(param) == 1000)
    {
        double sucRate = static_cast<double>(config.getSDCountParam(param))/static_cast<double>(config.getSDTrialsParam(param));
        config.setSDParam(param, config.getSDParam(param) * (std::exp((sucRate - 0.5)/10.0)));
        config.setSDCountParam(param, 0);
        config.setSDTrialsParam(param, 0);
    }
    config.setSDTrialsParam(param, config.getSDTrialsParam(param) + 1);

    // store old param value in case the new one will be rejected
    config.setTmpParam(param, config.getParam(param));
    double oldParamValue = config.getParam(param);
    
    // draw new paramter value
    std::normal_distribution<double> distribution(oldParamValue, config.getSDParam(param));
    double newParamValue = distribution(config.getGenerator());

    if (param == Config<TTreeType>::ParamType::mu || 
            param == Config<TTreeType>::ParamType::nu ||
            param == Config<TTreeType>::ParamType::lambda) 
    {
        newParamValue = adjusteCoefficient(newParamValue);
    }
    else if (param == Config<TTreeType>::ParamType::mutationOverDis)
    {
        newParamValue = adjustMutOverDis(newParamValue);
    }

    config.setParam(param, std::abs(newParamValue));
    if (param == Config<TTreeType>::ParamType::wildOverDis || param == Config<TTreeType>::ParamType::wildMean)
    {
        computeNoiseScore(config);
    }
    computeLogScoresOP(config);
}

template <typename TTreeType>
void 
proposeNextConfiguration(Config<TTreeType> & config)
{
	config.setMoveTyp(sampleRandomMove(config.moveProbs));      // pick the move type

	if(config.getMoveTyp() == 4) /* change one of the parameters */
    {
        changeParameters(config);
    }
    else
    {
        config.setTmpTree(config.getTree());
        config.numMutPlacements[1] = config.numMutPlacements[0];
        if(config.getMoveTyp()==1){ /* prune and re-attach */
            pruneAndReAttach(config);
        }
        else if(config.getMoveTyp()==2){ /* swap two node labels  */
            swapNodeLabels(config);
        }
        else if(config.getMoveTyp()==3){  /*  swap two subtrees  */
            swapSubtrees(config);    
        }
        else
        {
            std::cout << "This should not have happened!" << std::endl;
        }
        config.numMutPlacements[0] = getNumPlacements(config);
    }
}

template <typename TTreeType>
void manageBestTrees(Config<TTreeType> & /*config*/,
                     double & /*bestTreeLogScore*/,
                     double /*currTreeLogScore*/,
                     TTreeType & /*bestTrees*/,
                     double /*minDist*/)
{}

inline
void
manageBestTrees(Config<SampleTree> & config,
                     double & bestTreeLogScore,
                     double currTreeLogScore,
                     std::vector<Config<SampleTree>::TGraph> & bestTrees,
                     std::array<std::tuple<double, double>, 6> & bestParams,
                     double /*tag*/)
{
    if(currTreeLogScore > bestTreeLogScore)
    {   
        std::cout << "The new best score is: " <<  currTreeLogScore << std::endl;
        //save the current state to disk
        writeIndex(config);
        bestTreeLogScore = currTreeLogScore;
        bestTrees[0] = config.getTree();
        bestParams = config.params;

        //save the current best tree to disk
        boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, Vertex<SimpleTree>> newTreeBest = simplifyTree(config);
        std::ofstream ofs(config.outFilePrefix + ".gv");
        write_graphviz(ofs, newTreeBest, my_label_writer(newTreeBest, config.indexToPosition, config.cellNames, config.cellColours, config.cellClusters));
        ofs.close();
    }
}


// This function is actually doing the sampling of the mutations from the posterior distribution
inline
void 
updateMutInSampleCounts(Config<SampleTree> & config)
{
    Config<SampleTree>::TAttachmentScores & attachmentScores = config.getTmpAttachmentScore();
    Config<SampleTree>::TAttachmentScores::TAttachmentScore scoreSum;

    // for all mutations
    for (size_t attachment = 0; attachment < config.getNumAttachments(); ++attachment)
    {
        // get the attachment scores
        scoreSum.setMinusInfinity();
        getAllAttachmentScores(attachmentScores, scoreSum, config, attachment);

        // traverse the attachment scores to the leaves
        scoreSum.setMinusInfinity();
        PassScoreToChildrenBFSVisitor visBFS(config, attachmentScores, scoreSum, attachment);
        breadth_first_search(config.getTree(), num_vertices(config.getTree()) - 1, visitor(visBFS));

        // compute the probability of the tree to be in hetero- or homozygous state
        double logPHet = scoreSum.hetScore() + std::log((1.0 - config.getParam(Config<SampleTree>::nu))) + std::log(config.getNumSamples() * 2 - 1);
        double logPHom = scoreSum.homScore() + std::log( config.getParam(Config<SampleTree>::nu)) + std::log(config.getNumSamples() - 1);
        double logPD = addLogProb(logPHet, logPHom);
        logPHet -= logPD;
        logPHom -= logPD;
        
        // for all attachment positions
        for (unsigned i = 0; i < attachmentScores.size(); ++i)
        {
            if (config.getTree()[i].sample != -1)
            {
                // normalize probabilities
                attachmentScores[i] -= scoreSum;

                // compute final score
                attachmentScores[i].finalScore() = addLogProb(attachmentScores[i].hetScore() + logPHet, attachmentScores[i].homScore() + logPHom);

                // update the counter
                config.mutInSampleCounter[config.getTree()[i].sample][attachment].addInRealSpace(attachmentScores[i]);
            }
        }
    }
}


// this function is in experimental phase
inline
void 
updateMutInSampleCountsMix(Config<SampleTree> & config)
{
    Config<SampleTree>::TAttachmentScores & attachmentScores = config.getTmpAttachmentScore();
    Config<SampleTree>::TAttachmentScores::TAttachmentScore scoreSum;
    Config<SampleTree>::TAttachmentScores::TAttachmentScore scoreSum2;
    for (size_t attachment = 0; attachment < config.getNumAttachments(); ++attachment)
    {
        scoreSum2.setMinusInfinity();
        getAllAttachmentScores(attachmentScores, scoreSum, config, attachment);
        
        std::vector<double> lostScores;
        lostScores.resize(attachmentScores.size(), -INFINITY);
        double sumScoreReal = -INFINITY;
        ComputeLostScoreDFSVisitor vis(config, config._tmpAttachmentScore, lostScores, sumScoreReal, attachment);
        depth_first_search(config.getTree(), visitor(vis).root_vertex(num_vertices(config.getTree()) - 1));

        scoreSum.setMinusInfinity();
        PassScoreToChildrenBFSVisitor visBFS(config, attachmentScores, scoreSum, attachment);
        breadth_first_search(config.getTree(), num_vertices(config.getTree()) - 1, visitor(visBFS));
        
        double nu = config.getParam(Config<SampleTree>::nu);
        if (nu == 0)
        {
            nu += 1e-12;
        }

        double lambda = config.getParam(Config<SampleTree>::lambda);
        if (lambda == 0)
        {
            lambda += 1e-12;
        }

        double logPHet = scoreSum.hetScore() + std::log(1.0 - nu) + std::log(1.0 - lambda) - std::log(config.getNumSamples() * 2 - 1);
        double logPHom = scoreSum.homScore() + std::log(nu) - std::log(config.getNumSamples() -1);
        double logMixWild = scoreSum.mixWildScore() + std::log(1.0 - nu) + std::log(lambda) + std::log(0.5) - std::log(config.numMutPlacements[0]);
        double logMixHom = scoreSum.mixHomScore() + std::log(1.0 - nu) + std::log(lambda) + std::log(0.5) - std::log(config.numMutPlacements[0]);
        double logPD = addLogProb(addLogProb(logPHet, logPHom), addLogProb(logMixWild, logMixHom));
        logPHet -= logPD;
        logPHom -= logPD;
        logMixWild -= logPD;
        logMixHom -= logPD;

        for (unsigned i = 0; i < attachmentScores.size(); ++i)
        {
            if (config.getTree()[i].sample != -1)
            {
                // normalize probabilities
                attachmentScores[i] -= scoreSum;
                lostScores[i] -= scoreSum.mixWildScore();

                // if a node is attached to the root directly there is no lost-score
                if (lostScores[i] != -INFINITY)
                {
                    attachmentScores[i].mixWildScore() = subLogProb(attachmentScores[i].mixWildScore(), lostScores[i]);
                    attachmentScores[i].finalScore() = addLogProb(
                            addLogProb(attachmentScores[i].hetScore() + logPHet, attachmentScores[i].homScore() + logPHom), 
                            addLogProb(attachmentScores[i].mixWildScore() + logMixWild, attachmentScores[i].mixHomScore() + logMixHom)
                            );
                }
                else
                {
                    attachmentScores[i].finalScore() = addLogProb(attachmentScores[i].hetScore() + logPHet, attachmentScores[i].homScore() + logPHom);
                }

                // recompute the scores for storing them
                attachmentScores[i].hetScore() = attachmentScores[i].hetScore() + logPHet;
                attachmentScores[i].homScore() = attachmentScores[i].homScore() + logPHom;
                attachmentScores[i].mixWildScore() = attachmentScores[i].mixWildScore() + logMixWild;
                attachmentScores[i].mixHomScore() = attachmentScores[i].mixHomScore() + logMixHom;

                config.mutInSampleCounter[config.getTree()[i].sample][attachment].hetScore() = addLogProb(config.mutInSampleCounter[config.getTree()[i].sample][attachment].hetScore(), attachmentScores[i].hetScore());
                config.mutInSampleCounter[config.getTree()[i].sample][attachment].homScore() = addLogProb(config.mutInSampleCounter[config.getTree()[i].sample][attachment].homScore(), attachmentScores[i].homScore());
                if (attachmentScores[i].mixWildScore() != -INFINITY)
                {
                    config.mutInSampleCounter[config.getTree()[i].sample][attachment].mixWildScore() = addLogProb(config.mutInSampleCounter[config.getTree()[i].sample][attachment].mixWildScore(), attachmentScores[i].mixWildScore());
                }
                if (attachmentScores[i].mixHomScore() != -INFINITY)
                {
                    config.mutInSampleCounter[config.getTree()[i].sample][attachment].mixHomScore() = addLogProb(config.mutInSampleCounter[config.getTree()[i].sample][attachment].mixHomScore(), attachmentScores[i].mixHomScore());
                }
                if (attachmentScores[i].finalScore() != -INFINITY)
                {
                    config.mutInSampleCounter[config.getTree()[i].sample][attachment].finalScore() = addLogProb(config.mutInSampleCounter[config.getTree()[i].sample][attachment].finalScore(), attachmentScores[i].finalScore());
                }
            }
        }
    }
}

template <typename TTreeType>
double
runMCMC(std::vector<typename Config<TTreeType>::TGraph> & bestTrees, 
               std::array<std::tuple<double, double>, 6> & bestParams,
			   Config<TTreeType> & config,
			   std::vector<std::vector<unsigned>> & sampleTrees)
{


	int minDist = INT_MAX;                       // initialize distance to true tree if given
	double bestTreeLogScore = -DBL_MAX;          // initialize best tree score

	for(unsigned r=0; r<config.reps; r++){   // repeat the MCMC, start over with random tree each time, only best score and list of best trees is kept between repetitions
        
        config.updateContainers(0);
        computeNoiseScore(config);
        
        bestTrees.resize(1);
        bestTrees[0] = config.getTree();
        bestParams = config.params;
        config.numMutPlacements[0] = getNumPlacements(config);
        config.numMutPlacements[1] = config.numMutPlacements[0];

		double currTreeLogScore = scoreTree(config);       // get score of random tree
        manageBestTrees(config,
                bestTreeLogScore, 
                currTreeLogScore,
                bestTrees,
                bestParams,
                minDist);


        double propTreeLogScore;
        double random;

        for(unsigned it=0; it<config.loops + config.sampleLoops; it++){        // run the iterations of the MCMC

            // Mutations are added in batches of 10% of all mutations 10 times
            // This function handels the resizing of the corresponding containers
            if (config.updateContainers(it))
            {
                minDist = INT_MAX;                       // initialize distance to true tree if given
                bestTreeLogScore = -DBL_MAX;          // initialize best tree score
                config.getTree() = bestTrees[0];
		        currTreeLogScore = scoreTree(config);       // get score of random tree
                manageBestTrees(config,
                        bestTreeLogScore, 
                        currTreeLogScore,
                        bestTrees,
                        bestParams,
                        minDist);
            }
        	
            if(it % 10000 == 0)
			{ 
                config.printParameters();
                std::cout.precision(15);
                std::cout << "iterations: " << it << std::endl;
                std::cout << "score: " << currTreeLogScore << std::endl;
                std::cout << "bestScore: " << bestTreeLogScore << std::endl;
                std::cout << "num trees: " << bestTrees.size() << std::endl;
			}

            // sample the next tree configuration
            proposeNextConfiguration(config); 
            propTreeLogScore = scoreTree(config);


            random = (double) rand() / RAND_MAX;
            
            // the proposed tree is accepted
            if (random < std::exp((propTreeLogScore-currTreeLogScore)))
            {
                // update counter of parameter estimation
                if (config.getMoveTyp() == 4)
                    config.setSDCountParam(config.getParamToOptimize(), config.getSDCountParam(config.getParamToOptimize()) + 1);
                currTreeLogScore  = propTreeLogScore;

                manageBestTrees(config,
                        bestTreeLogScore, 
                        currTreeLogScore,
                        bestTrees,
                        bestParams,
                        minDist);
            }
            else // the proposed tree is rejected
            {
                if (config.getMoveTyp() == 4)
                {
                    config.resetParameters();
                }
                else
                {
                    config.getTree().swap(config.getTmpTree());
                    config.numMutPlacements[0] = config.numMutPlacements[1]; 
                }
            }
           
            // sample from the posterior distribution
            if (it >= config.loops)
            {
                if (config.computeMixScore) // experimental
                {
                    updateMutInSampleCountsMix(config);
                }
                else
                {
                    updateMutInSampleCounts(config);
                }
                config.updateParamsCounter();
            }
        }
	} 
	
 	return bestTreeLogScore;
}

inline
void 
normalizeMutationCounts(Config<SampleTree> & config)
{
    for (unsigned int j = 0; j < config.getNumMutations(); ++j)
    {
        for (unsigned int i = 0; i < config.getNumSamples(); ++i)
        {
            config.mutInSampleCounter[i][j].exp();
			if (j < config.getNumMutations() - config.numUniqMuts)
            {
                config.mutInSampleCounter[i][j] /= static_cast<double>(config.sampleLoops);
            }
		}
	}
}

#endif
