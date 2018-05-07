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
#define _SECURE_SCL 1

#include <stdbool.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <iostream>
#include <sstream>
#include <time.h>
#include <tuple>
#include <cmath>
#include <assert.h>

#include "sciphi_config.h"
#include "mcmc.h"
#include "readData.h"
#include "output.h"

#include <boost/multi_array.hpp>
#include <boost/program_options.hpp>


template <typename TTreeType>
int readParameters(Config<TTreeType> & config, int argc, char* argv[]);

template <typename TTreeType>
std::vector<double> setMoveProbs(Config<TTreeType> const & config);

// This function resets weights for the different move operations
// The weights determine how much time is spend in the different 
// move/parameter estimation steps.
template <typename TTreeType>
void setMoveProbs(Config<TTreeType> & config, double newParamsEstimateRate){

	config.moveProbs[0] = config.moveProbs[0] / (1.0 - config.moveProbs[3]);
	config.moveProbs[0] = config.moveProbs[0] * (1.0 - newParamsEstimateRate);
	config.moveProbs[1] = config.moveProbs[1] / (1.0 - config.moveProbs[3]);
	config.moveProbs[1] = config.moveProbs[1] * (1.0 - newParamsEstimateRate);
	config.moveProbs[2] = config.moveProbs[2] / (1.0 - config.moveProbs[3]);
	config.moveProbs[2] = config.moveProbs[2] * (1.0 - newParamsEstimateRate);
    config.moveProbs[3] = 1 - config.moveProbs[0] - config.moveProbs[1] - config.moveProbs[2];

	return;
}


// Function to read the parameters from the command line.
template <typename TTreeType>
int readParameters(Config<TTreeType> & config, int argc, char* argv[]){

    boost::program_options::options_description generic("Generic options");
	generic.add_options()("help,h", "Print this help");

	// configuration options
	boost::program_options::options_description parseConfig("Configuration");
	parseConfig.add_options()
		(",i", boost::program_options::value<decltype(config.refFileName)>(&config.refFileName), "Reference file.")
		("il", boost::program_options::value<decltype(config.loadName)>(&config.loadName), "Directory from which to read intermediate results.")
		("ol", boost::program_options::value<decltype(config.saveName)>(&config.saveName), "Directory to store intermediate results.")
		("ex", boost::program_options::value<decltype(config.exclusionFileName)>(&config.exclusionFileName), "Filename of exclusion list (VCF format).")
		("me", boost::program_options::value<decltype(config.mutationExclusionFileName)>(&config.mutationExclusionFileName), "Filename of mutations to exclude during the sequencing error rate estimation (VCF format).")
		("in", boost::program_options::value<decltype(config.bamFileNames)>(&config.bamFileNames), "Name of the BAM files used to create the mpileup.")
		(",o", boost::program_options::value<decltype(config.outFilePrefix)>(&config.outFilePrefix), "Prefix of output files.")
		(",r", boost::program_options::value<decltype(config.reps)>(&config.reps), "Number of repetitions. [1]")
		(",l", boost::program_options::value<decltype(config.loops)>(&config.loops), "Maximal number of iterations per repetition. [1000000]")
		("lz", boost::program_options::value<decltype(config.learnZygocity)>(&config.learnZygocity), "Set to 1 if zygocity should be learned. [0]")
		("ls", boost::program_options::value<decltype(config.sampleLoops)>(&config.sampleLoops), "Number of sample iterations. [100000]")
		("pr", boost::program_options::value<decltype(config.priorMutationRate)>(&config.priorMutationRate), "Prior mutation rate [0.0001].")
		("uniq,u", boost::program_options::value<decltype(config.uniqTreshold)>(&config.uniqTreshold), "Filter mutations showing up to this number of cells showing the mutations. [0]")
		(",e", boost::program_options::value<decltype(config.paramsEstimateRate)>(&config.paramsEstimateRate), "Paramter estimation rate, i.e. the frection of loops used to estimate the different parameters. [0.2]")
		("seed", boost::program_options::value<decltype(config.fixedSeed)>(&config.fixedSeed), "Seed for the random number generator. [42]")
		(",t", boost::program_options::value<decltype(config.scoreType)>(&config.scoreType), "Tree score type [m (max), s (sum)]. [s]")
		("wildOverDis", boost::program_options::value<double>(&std::get<0>(config.params[0])), "Overdispersion for wild type. [100]")
		("mutationOverDis", boost::program_options::value<double>(&std::get<0>(config.params[1])), "Overdispersion for mutant type. [2]")
		("wildMean", boost::program_options::value<double>(&std::get<0>(config.params[2])), "Mean error rate. [0.001]")
		("cwm", boost::program_options::value<decltype(config.numCellWithMutationMin)>(&config.numCellWithMutationMin), "Number of cells requiered to have a mutation in order to be called. [1]")
		("mc", boost::program_options::value<decltype(config.minCoverage)>(&config.minCoverage), "Minimum coverage required per cell. [1]")
		("mcc", boost::program_options::value<decltype(config.minCoverageAcrossCells)>(&config.minCoverageAcrossCells), "Minimum coverage required by at least x cells which support the alternative. [0]")
		("nmc", boost::program_options::value<decltype(config.numMinCoverageAcrossCells)>(&config.numMinCoverageAcrossCells), "Number of cells which need a minimum coverage which support the alternative.")
		("ms", boost::program_options::value<decltype(config.minSupport)>(&config.minSupport), "Minimum number of reads required to support the alternative. [0]")
		("mf", boost::program_options::value<decltype(config.minFreq)>(&config.minFreq), "Minimum required frequency of reads supporting the alternative per cell. [0]")
		("mff", boost::program_options::value<decltype(config.meanFilter)>(&config.meanFilter), "Mean of acceptable variant allele frequency across all cells for a specific locus. [0.25]")
		("bns", boost::program_options::value<decltype(config.maxSupInControlBulk)>(&config.maxSupInControlBulk), "Loci with up to this number of alternative supporting reads in the bulk control sample will be skiped. [2]")
		("bnc", boost::program_options::value<decltype(config.minCovInControlBulk)>(&config.minCovInControlBulk), "Minimum required coverage of reads in the bulk control sample. [6]")
		("mnc", boost::program_options::value<decltype(config.maxNumberNormalCellMutated)>(&config.maxNumberNormalCellMutated), "Maximum number of control cells allowed to be mutated. [0]")
		("unc", boost::program_options::value<decltype(config.useNormalCellsInTree)>(&config.useNormalCellsInTree), "Use normal cells for tree reconstruction. [false]")
		//("cms", boost::program_options::value<decltype(config.computeMixScore)>(&config.computeMixScore), "Compute the mix score = allow a mutation to be homozygous or lost in a subtree.") - feature still in experimental phase
		//("mspa", boost::program_options::value<double>(&config.clamPrior[0]), "The alpha parameter of the beta function used as prior on the rate of loosing a chromosome after some mutation occurred.") - feature still in experimental phase
		//("mspb", boost::program_options::value<double>(&config.clamPrior[1]), "The beta parameter of the beta function used as prior on the rate of loosing a chromosome after some mutation occurred.") - feature still in experimental phase
		("ese", boost::program_options::value<decltype(config.estimateSeqErrorRate)>(&config.estimateSeqErrorRate), "Estimate the sequencing error rate. [1]");

   // hidden options, i.e., input files
	boost::program_options::options_description hidden("Hidden options");
	hidden.add_options()
		("input-file", boost::program_options::value<std::string>()->required(), "input files");

	boost::program_options::options_description cmdline_options;
	cmdline_options.add(generic).add(parseConfig).add(hidden);
	boost::program_options::options_description visible("Allowed options");
	visible.add(generic).add(parseConfig);
	boost::program_options::positional_options_description p;
	p.add("input-file", -1);
	boost::program_options::variables_map global_options;

	/* 0.0) parse program options */
	try
	{
		boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), global_options);

		// show help options
		if (global_options.count("help"))
		{
			std::cout << visible << '\n';
			exit(EXIT_SUCCESS);
		}

		boost::program_options::notify(global_options);
	}
	catch (boost::program_options::required_option& e)
	{
		if (e.get_option_name() == "--input-file")
		{
			std::cerr << "ERROR: You have provided no input files. sciphi requires at least one input file.\n";
		}
		else
		{
			std::cerr << "ERROR: " << e.what() << '\n';
		}
		exit(EXIT_FAILURE);
	}
	catch (boost::program_options::error& e)
	{
		std::cerr << "ERROR: " << e.what() << '\n';
		exit(EXIT_FAILURE);
	} 

    if (global_options.count("-e"))
    {
        setMoveProbs(config, config.paramsEstimateRate);
    }
        
    config.inFileName = global_options["input-file"].as<std::string>();

    std::get<1>(config.params[0]) = std::get<0>(config.params[0]);
    std::get<1>(config.params[1]) = std::get<0>(config.params[1]);
    std::get<1>(config.params[2]) = std::get<0>(config.params[2]);
	
    return 0;
}

int main(int argc, char* argv[])
{
    typedef Config<SampleTree> TConfig;
	TConfig config{};

    //readNucInfo(config);
    //readGraph(config);
    //return 0;

    // read the command line arguments
    std::cout << "Reading the config file: ... " << std::flush;
	readParameters(config, argc, argv);
    std::cout << "done!" << std::endl;

    // extract the mutation data from the files
    if (config.loadName == "" )
    {
        std::cout << "Reading the mpileup file: " << std::flush;
        getData(config);
        createInitialTree(config);
        std::cout << "done!" << std::endl;
    }
    else
    {
        std::cout << "Reading the stored results file: " << std::flush;
        readCellNames(config);
        readGraph(config);
        readNucInfo(config);
        std::cout << "done!" << std::endl;

    }

    // After the data is read the log scores need to be computed.
	computeLogScoresOP(config);
    config.paramsCounter.resize(0);

	std::vector<typename TConfig::TGraph> optimalTrees;         // list of all optimal trees (as ancestor matrices)
    std::array<std::tuple<double, double>, 6> optimalParams;    // The optimal parameters
	std::vector<std::vector<unsigned>> sampleTrees;             // list where tree samples are stored, if sampling based on posterior distribution is needed

	/* initialize the random number generator, either with a user defined seed, or a random number */
	if(config.fixedSeed == static_cast<unsigned>(-1)){
		initRand();                                  // initialize random number generator
	}
	else{
		srand(config.fixedSeed);
	}

	// Find best scoring trees by MCMC
    config.initMutInSampleCounter(); 
    runMCMC(optimalTrees, optimalParams, config, sampleTrees);

    // Write the results to disk
    printMutation2SampleAssignment(config, config.outFilePrefix + "_mut2Sample.tsv");
    normalizeMutationCounts(config);
    printMutationProbability(config, config.outFilePrefix + ".probs");
    writeParameters(config, config.outFilePrefix + ".params.txt");
    writeVCF(config, config.outFilePrefix + ".vcf");
    
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, Vertex<SimpleTree>> TGraph;
    // the following can be included for debugging purposes
    //TGraph newTree = simplifyTree(config);
    //std::ofstream ofs(config.outFilePrefix + "_last.gv");
    //write_graphviz(ofs, newTree, my_label_writer(newTree, config.indexToPosition, config.cellNames));

    config.getTree() = optimalTrees[0];
    writeIndex(config);

    config.params = optimalParams;
    computeLogScoresOP(config);
    TGraph newTreeBest = simplifyTree(config);
    std::ofstream ofs2(config.outFilePrefix + ".gv");
    write_graphviz(ofs2, newTreeBest, my_label_writer(newTreeBest, config.indexToPosition, config.cellNames));
}
