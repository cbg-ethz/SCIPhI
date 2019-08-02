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
#include "version.h"

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

    std::string lastName; // Dummy for a deprecated option value
	// configuration options
	boost::program_options::options_description parseConfig("Configuration");
	parseConfig.add_options()
		("in", boost::program_options::value<decltype(config.bamFileNames)>(&config.bamFileNames), "Name of the BAM files used to create the mpileup.")
		("il", boost::program_options::value<decltype(config.loadName)>(&config.loadName), "Directory from which to read intermediate results.")
		("ex", boost::program_options::value<decltype(config.exclusionFileName)>(&config.exclusionFileName), "File name of exclusion list (VCF format), containing loci which should be ignored.")
        ("me", boost::program_options::value<decltype(config.mutationExclusionFileName)>(&config.mutationExclusionFileName), "File name of mutations to exclude during the sequencing error rate estimation (VCF format).")
        ("inc", boost::program_options::value<decltype(config.variantInclusionFileName)>(&config.variantInclusionFileName), "File name of inclusion list (VCF format) containing Variants (CHROM, POS, REF, ALT) that should be included.")
		(",o", boost::program_options::value<decltype(config.outFilePrefix)>(&config.outFilePrefix), "Prefix of output files.")
		("ol", boost::program_options::value<decltype(config.lastName)>(&lastName), "This option is deprecated. The index will be saved in a folder specified with -o in \"last_index\". If desired, one can continue from here to explore more search space.")
		("sa", boost::program_options::value<decltype(config.sampling)>(&config.sampling), "Sampling step. If a value x different from 0 is specified, every x iteration, after the burn in phase, an index will be writen to disk to provide a posterior sampling. [0]")
		//(",r", boost::program_options::value<decltype(config.reps)>(&config.reps), "Number of repetitions. [1]")
		(",l", boost::program_options::value<decltype(config.loops)>(&config.loops), "Maximal number of iterations before sampling form the posterior distribution of the mutation to cell assignment. [1000000]")
		("ls", boost::program_options::value<decltype(config.sampleLoops)>(&config.sampleLoops), "Number of iterations in which the mutation to cell assignment is sampled. [100000]")
		("lz", boost::program_options::value<decltype(config.learnZygocity)>(&config.learnZygocity), "Set to 1 if zygosity should be learned. The zygosity rate is the fraction of mutations which are homozygous in all cells. [0]")
		("zyg", boost::program_options::value<double>(&std::get<0>(config.params[Config<TTreeType>::nu])), "Zygosity rate. [0]")
		("pr", boost::program_options::value<decltype(config.priorMutationRate)>(&config.priorMutationRate), "Prior mutation rate [0.0001].")
		("uniq,u", boost::program_options::value<decltype(config.uniqTreshold)>(&config.uniqTreshold), "Mark loci with this number of cells being mutated as \"PASS\". [1]")
		(",e", boost::program_options::value<decltype(config.paramsEstimateRate)>(&config.paramsEstimateRate), "Parameter estimation rate, i.e. the fraction of loops used to estimate the different parameters. [0.2]")
        ("ur",boost::program_options::value<double>(&std::get<1>(config.dataUsageRate)), "Data usage rate increment steps. In order to speed up the algorithm one can choose to iteratively add more and more of the candidate loci while learning the tree model. For example, 0.1 would mean to start with 10% of the data, then use 20%, and so on. Per default, all data is used from the start. [1]")
		("seed", boost::program_options::value<decltype(config.fixedSeed)>(&config.fixedSeed), "Seed for the random number generator. [42]")
		(",t", boost::program_options::value<decltype(config.scoreType)>(&config.scoreType), "Tree score type [m (max), s (sum)]. [s]")
		("wildMean", boost::program_options::value<double>(&std::get<0>(config.params[2])), "Mean error rate. If the sequencing error rate should not be learned \"--ese 0\" one can specify it. [0.001]")
		("wildOverDis", boost::program_options::value<double>(&std::get<0>(config.params[0])), "Initial overdispersion for wild type. The overdispersion is learned during the tree traversal. [100]")
		("mutationOverDis", boost::program_options::value<double>(&std::get<0>(config.params[1])), "Initial overdispersion for mutant type. The overdispersion is learned during the tree traversal.[2]")
		("ese", boost::program_options::value<decltype(config.estimateSeqErrorRate)>(&config.estimateSeqErrorRate), "Estimate the sequencing error rate. [1]")
		("sub", boost::program_options::value<decltype(config.sub)>(&config.sub), "PCR substitution rate. An error early during the PCR amplification can result in high allele fractions. These are typically filtered out by requiring at least two cells to be mutated. However, sometimes the PCR error in one cell is a real mutation in others. If in such a scenario the coverage for the cell with the PCR error is high, SCIPhI will introduce false positives. This can be avoided by specifying the PCR error rate. As these events are rather unlikely very small values, such as 0.0000001 should be chosen. In practice 0 works very well.  [0]")
		("md", boost::program_options::value<decltype(config.minDist)>(&config.minDist), "Window size in which only x (see option \"--mmw\") mutations are allowed to be present. Some sequences are prone to produce sequencing errors, which can be avoided by specifying a range in which only x mutations are allowed. [1]")
		("mmw", boost::program_options::value<decltype(config.maxMutPerWindow)>(&config.maxMutPerWindow), "Maximum number of mutations allowed per window. [1]")
		("cwm", boost::program_options::value<decltype(config.numCellWithMutationMin)>(&config.numCellWithMutationMin), "Number of tumor cells required to have a mutation in order to be called. [2]")
		("mnp", boost::program_options::value<decltype(config.minNumCellsPassFilter)>(&config.minNumCellsPassFilter), "Number of cells which need to pass the filters described below. [2]")
		("mc", boost::program_options::value<decltype(config.minCoverage)>(&config.minCoverage), "Minimum coverage required per cell. [1]")
		("ms", boost::program_options::value<decltype(config.minSupport)>(&config.minSupport), "Minimum number of reads required to support the alternative. [3]")
		("mf", boost::program_options::value<decltype(config.minFreq)>(&config.minFreq), "Minimum required frequency of reads supporting the alternative per cell. [0]")
		("mff", boost::program_options::value<decltype(config.meanFilter)>(&config.meanFilter), "Mean of acceptable variant allele frequency across all cells for a specific locus. Mapping artifacts may result in low allele frequencies across cells. In order to filter these events out we apply a log-likelihood ratio test where the sequencing error model has a mean of this value. [0.25]")
		("bns", boost::program_options::value<decltype(config.maxSupInControlBulk)>(&config.maxSupInControlBulk), "Loci with up to this number of alternative supporting reads in the bulk control sample will be skipped as germline. [2]")
		("bnc", boost::program_options::value<decltype(config.minCovInControlBulk)>(&config.minCovInControlBulk), "Minimum required coverage of reads in the bulk control sample. [6]")
		("ncf", boost::program_options::value<decltype(config.normalCellFilter)>(&config.normalCellFilter), "Normal cell filter. Currently there are three options: (0) Do not use the normal cells for filtering; (1) use a simple filtering scheme excluding mutations if the probability of being mutated is higher than not being mutated for any cell independently; (2) filter mutations where the probability that at least one cell is mutated is higher than no cell is mutated. Note that in contrast to (1) the cells are not independent and cells with no alternative support need to be explained via dropout events. [1]")
		("mnc", boost::program_options::value<decltype(config.maxNumberNormalCellMutated)>(&config.maxNumberNormalCellMutated), "Maximum number of control cells allowed to be mutated. [0]")
		("unc", boost::program_options::value<decltype(config.useNormalCellsInTree)>(&config.useNormalCellsInTree), "Use normal cells for tree reconstruction. [false]");
		//("cms", boost::program_options::value<decltype(config.computeMixScore)>(&config.computeMixScore), "Compute the mix score = allow a mutation to be homozygous or lost in a subtree.") - feature still in experimental phase
		//("mspa", boost::program_options::value<double>(&config.clamPrior[0]), "The alpha parameter of the beta function used as prior on the rate of loosing a chromosome after some mutation occurred.") - feature still in experimental phase
		//("mspb", boost::program_options::value<double>(&config.clamPrior[1]), "The beta parameter of the beta function used as prior on the rate of loosing a chromosome after some mutation occurred.") - feature still in experimental phase

   // hidden options, i.e., input files
	boost::program_options::options_description hidden("Hidden options");
	hidden.add_options()
		("input-file", boost::program_options::value<std::string>()->required(), "input files")
		("smt", boost::program_options::value<decltype(config.mutToMaxName)>(&config.mutToMaxName), "Store to save mutations distribution of MAP tree.");

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
    config.bestName = config.outFilePrefix + "/best_index/";
    config.lastName = config.outFilePrefix + "/last_index/";
    config.samplingName = config.outFilePrefix + "/sample_index";

    std::get<1>(config.params[0]) = std::get<0>(config.params[0]);
    std::get<1>(config.params[1]) = std::get<0>(config.params[1]);
    std::get<1>(config.params[2]) = std::get<0>(config.params[2]);
    std::get<1>(config.params[3]) = std::get<0>(config.params[3]);
    std::get<1>(config.params[4]) = std::get<0>(config.params[4]);
    std::get<1>(config.params[5]) = std::get<0>(config.params[5]);
	
    return 0;
}

int main(int argc, char* argv[])
{
    typedef Config<SampleTree> TConfig;
	TConfig config{};
	
    // read the command line arguments
    std::cout << "SCIPhI v" << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << "-" << VERSION_STRING << std::endl;
    std::cout << "Reading the config file: ... " << std::flush;
	readParameters(config, argc, argv);
    std::cout << "done!" << std::endl;

    /* initialize the random number generator, either with a user defined seed, or a random number */
	if(config.fixedSeed == static_cast<unsigned>(-1)){
		initRand();                                  // initialize random number generator
	}
	else{
		srand(config.fixedSeed);
	}


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

    if (config.mutToMaxName != "")
    {
        config.updateContainers(0);
        writeMutToMaxTree(config);
        return 0;
    }
    
	std::vector<typename TConfig::TGraph> optimalTrees;         // list of all optimal trees (as ancestor matrices)
    std::array<std::tuple<double, double>, 6> optimalParams;    // The optimal parameters
	std::vector<std::vector<unsigned>> sampleTrees;             // list where tree samples are stored, if sampling based on posterior distribution is needed

	// Find best scoring trees by MCMC
    config.initMutInSampleCounter(); 
    runMCMC(optimalTrees, optimalParams, config, sampleTrees);

    // Write the results to disk
    printMutation2SampleAssignment(config, config.outFilePrefix + "_mut2Sample.tsv");
    normalizeMutationCounts(config);
    printMutationProbability(config, config.outFilePrefix + ".probs");
    writeParameters(config, config.outFilePrefix + ".params.txt");
    writeVCF(config, config.outFilePrefix + ".vcf");
    
    writeFinalIndex(config);

    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, Vertex<SimpleTree>> TGraph;
    config.getTree() = optimalTrees[0];
    config.params = optimalParams;
    computeLogScoresOP(config);
    TGraph newTreeBest = simplifyTree(config);
    std::ofstream ofs2(config.outFilePrefix + ".gv");
    write_graphviz(ofs2, newTreeBest, my_label_writer(newTreeBest, config.indexToPosition, config.cellNames, config.cellColours, config.cellClusters));
}
