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

#include "sciphi_config.h"

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <iomanip>
#include <ctime>
#include <math.h>

void writeToFile(std::string content, std::string fileName);
std::string getGraphVizFileContentNumbers(int* parents, int n);
void getGraphVizFileContentNames(std::string & output, 
								 Config<MutationTree> & config,
								 std::vector<unsigned> const & parents, 
								 std::vector<std::string> const & geneNames, 
								 bool attachSamples, 
								 std::vector<std::vector<bool>> const & ancMatrix);
std::string getBestAttachmentString(std::vector<std::vector<bool>> const & ancMatrix, 
									Config<MutationTree> & config,
									std::vector<std::string> geneNames);
bool** attachmentPoints(std::vector<std::vector<bool>> const & ancMatrix, 
						Config<MutationTree> const & config);
void printSampleTrees(std::vector<int*> list, int n, std::string fileName);
void printGraphVizFile(int* parents, int n);
void printScoreKimSimonTree(Config<MutationTree> & config,
							char scoreType);

inline
void 
printMutationProbability(Config<SampleTree> const & config,
        std::string const & outfileName)
{

    std::ofstream outFile;
    outFile.open(outfileName);
    outFile << "chrom\tposition\t";
	for (unsigned int i = 0; i < config.getNumSamples(); ++i)
        outFile << config.cellNames[i] << "\t";
    outFile << "\n";
    for (unsigned int j = 0; j < config.getNumMutations(); ++j)
    {
        outFile << std::get<0>(config.indexToPosition[j]) << "\t" << std::get<1>(config.indexToPosition[j]) << "\t";
        for (unsigned int i = 0; i < config.getNumSamples(); ++i)
        {
			//outFile << config.mutInSampleCounter[i][j] << "\t";
			outFile << config.mutInSampleCounter[i][j].hetScore() << "|" << 
                       config.mutInSampleCounter[i][j].homScore() << "|" <<
                       config.mutInSampleCounter[i][j].finalScore() << "\t";
		}
		outFile << std::endl;
	}
    outFile.close();
}

inline
void 
writeVCFHeader(Config<SampleTree> const & config,
               std::ofstream & outFile)
{
    outFile << "##fileformat=VCFv4.1\n";
   
    // get current date
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    outFile << "##fileDate=" << std::put_time(&tm, "%d-%m-%Y") <<"\n";
    
    outFile << "##source=SCATE\n";
    outFile << "##FILTER=<ID=LowQual,Description=\"Low quality\">\n";
    outFile << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">\n";
    outFile << "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for alt alleles\">\n";
    outFile << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">\n";
    outFile << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";
    outFile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    outFile << "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n";
    
    outFile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    if (config.bamFileNames != "")
    {
        for (unsigned i = 0; i < config.getNumSamples(); ++i)
        {
            outFile << "\t" << config.cellNames[i];
        }
    }
    else
    {
        for (unsigned i = 0; i < config.getNumSamples(); ++i)
        {
            outFile << "\tcell" << i;
        }
    }
    outFile << "\n";
    outFile.close();
}


inline
void 
writeVCFSingleEntry(std::ofstream & outFile,
        Config<SampleTree> const & config,
         std::vector<unsigned> & entries)
{
    // chrom
    outFile << std::get<0>(config.indexToPosition[entries[0]]) << "\t";
    // pos
    outFile << std::get<1>(config.indexToPosition[entries[0]]) << "\t";
    // id
    outFile << "*\t";
    // ref
    outFile << std::get<2>(config.indexToPosition[entries[0]]) << "\t";
    // alt
    outFile << std::get<3>(config.indexToPosition[entries[0]]) << "\t";
    // qual
    outFile << "*\t";
    // filter
    outFile << "PASS\t";
    // info
    outFile << "DP:";
    unsigned dpAll = 0;
    for (unsigned sample = 0; sample < config.getNumSamples(); ++sample)
        dpAll += std::get<0>(config.getData()[entries[0]][sample]);
    outFile << dpAll;
    outFile << "\n";
    
}

inline
bool
checkIfMutationIsUniq(Config<SampleTree> const & config,
         std::vector<unsigned> & entries,
         unsigned mutationPos)
{
    unsigned numMutatedCells = 0;
    for (unsigned sample = 0; sample < config.getNumSamples(); ++sample)
    {
        if (config.mutInSampleCounter[sample][mutationPos].finalScore() >= 0.95)
        {
            ++numMutatedCells;
        }
    }
    return numMutatedCells <= config.uniqTreshold;
}

inline
void 
writeVCFMultiEntry(std::ofstream & outFile,
        Config<SampleTree> const & config,
         std::vector<unsigned> & entries)
{
    for (unsigned alt = 0; alt < entries.size(); ++alt)
    {
        // chrom
        outFile << std::get<0>(config.indexToPosition[entries[alt]]) << "\t";
        // pos
        outFile << std::get<1>(config.indexToPosition[entries[alt]]) << "\t";
        // id
        outFile << "*\t";
        // ref
        outFile << std::get<2>(config.indexToPosition[entries[alt]]) << "\t";
        // alt
        outFile << std::get<3>(config.indexToPosition[entries[alt]]) << "\t";
        // qual
        outFile << "*\t";
        // filter
        if (checkIfMutationIsUniq(config, entries, entries[alt]))
        {
            outFile << ".\t";
        }
        else
        {
            outFile << "PASS\t";
        }
        // info
        outFile << "DP:";
        unsigned dpAll = 0;
        for (unsigned sample = 0; sample < config.getNumSamples(); ++sample)
            dpAll += std::get<0>(config.getData()[sample][entries[alt]]);
        outFile << dpAll << "\t";
        // format
        outFile << "GT:AD:DP:GQ:PL\t";
        for (unsigned sample = 0; sample < config.getNumSamples(); ++sample)
        {
            if (config.mutInSampleCounter[sample][entries[alt]].finalScore() > 0.95)
            {
                outFile << "0/1:";
            }
            else
            {
                outFile << "0/0:";
            }
            outFile << std::get<1>(config.getData()[sample][entries[alt]]) << ":";
            outFile << std::get<0>(config.getData()[sample][entries[alt]]) << ":";
            outFile << static_cast<unsigned>(-10 * std::log10(std::max(0.00000001, 1 - config.mutInSampleCounter[sample][entries[alt]].finalScore()))) << "\t";
        }
        outFile << "\n";
    }
}

inline
void 
writeVCFEntry(std::ofstream	& outFile,
        Config<SampleTree> const & config,
         std::vector<unsigned> & entries)
{
        writeVCFMultiEntry(outFile, config, entries);
}


inline
void 
writeVCF(Config<SampleTree> const & config,
        std::string const & outfileName)
{
    // Create the output file
    std::ofstream outFile;
    outFile.open(outfileName);

    // print the header information
    writeVCFHeader(config, outFile);

    // create a map to sort the mutations according to chromosome and position
    std::map<std::pair<std::string, unsigned>, std::vector<unsigned>> posMap;
    for (unsigned i = 0; i < config.indexToPosition.size(); ++i)
    {
        std::pair<std::string, unsigned> newPos (std::get<0>(config.indexToPosition[i]),std::get<1>(config.indexToPosition[i]));
        auto it = posMap.find(newPos);
        if (it != posMap.end())
        {
            it->second.push_back(i);
        }
        else
        {
            posMap[newPos] = {i};
        }
    }

    for (auto it = posMap.begin(); it != posMap.end(); ++it)
    {
        writeVCFEntry(outFile, config, it->second);
    }
}

void writeTree(Config<SampleTree> const & config)
{
    std::cout << config.saveName + "/tree.gv" << std::endl;
    std::ofstream ofs(config.saveName + "/tree.gv");
    write_graphviz(ofs, config.getTree(), my_label_writer_complete(config.getTree()));

}

void writeNucInfo(Config<SampleTree> const & config)
{
    string makeDir = "mkdir -p " + config.saveName;
    std::system(makeDir.c_str());
    std::ofstream outFile;
    outFile.open(config.saveName + "/nuc.tsv");
    
    outFile << "=numSamples=" << "\n";
    outFile << config.getNumSamples() << "\n";
    
    outFile << "=params=" << "\n";
    outFile << config.getParam(Config<SampleTree>::wildMean) << "\t"<< config.getSDParam(Config<SampleTree>::wildMean) << "\t" << config.getSDCountParam(Config<SampleTree>::wildMean) << "\t" << config.getSDTrialsParam(Config<SampleTree>::wildMean) << "\n";
    outFile << config.getParam(Config<SampleTree>::wildOverDis) << "\t" << config.getSDParam(Config<SampleTree>::wildOverDis) << "\t" << config.getSDCountParam(Config<SampleTree>::wildOverDis) << "\t" << config.getSDTrialsParam(Config<SampleTree>::wildOverDis) << "\n";
    outFile << config.getParam(Config<SampleTree>::mutationOverDis) << "\t"<< config.getSDParam(Config<SampleTree>::mutationOverDis) << "\t" << config.getSDCountParam(Config<SampleTree>::mutationOverDis) << "\t" << config.getSDTrialsParam(Config<SampleTree>::mutationOverDis) << "\n";
    outFile << config.getParam(Config<SampleTree>::mu) << "\t"<< config.getSDParam(Config<SampleTree>::mu) << "\t" << config.getSDCountParam(Config<SampleTree>::mu) << "\t" << config.getSDTrialsParam(Config<SampleTree>::mu) << "\n";
    outFile << config.getParam(Config<SampleTree>::nu) << "\t"<< config.getSDParam(Config<SampleTree>::nu) << "\t" << config.getSDCountParam(Config<SampleTree>::nu) << "\t" << config.getSDTrialsParam(Config<SampleTree>::nu) << "\n";

    std::cout << config.indexToPosition.size() << " " << config.getData().size() << std::endl;

    outFile << "=mutations=" << "\n";
    for (unsigned i = 0; i < config.indexToPosition.size(); ++i)
    {
        // chrom
        outFile << std::get<0>(config.indexToPosition[i]) << "\t";
        // pos
        outFile << std::get<1>(config.indexToPosition[i]) << "\t";
        // ref
        outFile << std::get<2>(config.indexToPosition[i]) << "\t";
        // alt
        outFile << std::get<3>(config.indexToPosition[i]);

        for (unsigned j = 0; j < config.getData().size(); ++j)
        {
            outFile << "\t" << std::get<0>(config.getData()[j][i]) << "\t" << std::get<1>(config.getData()[j][i]);
        }
        outFile << "\n";
    }

    outFile << "=background=" << "\n";
    if (config.noiseCounts.cov.size() > 0)
    {
        outFile << config.noiseCounts.cov[0].first << "\t" << config.noiseCounts.cov[0].second;
        for (unsigned i = 1; i < config.noiseCounts.cov.size(); ++i)
        {
            outFile << "\t" << config.noiseCounts.cov[i].first << "\t" << config.noiseCounts.cov[i].second;
        }
        outFile << "\n";
    }

    if (config.noiseCounts.sup.size() > 0)
    {
        outFile << config.noiseCounts.sup[0].first << "\t" << config.noiseCounts.sup[0].second;
        for (unsigned i = 1; i < config.noiseCounts.sup.size(); ++i)
        {
            outFile << "\t" << config.noiseCounts.sup[i].first << "\t" << config.noiseCounts.sup[i].second;
        }
        outFile << "\n";
    }

    if (config.noiseCounts.covMinusSup.size() > 0)
    {
        outFile << config.noiseCounts.covMinusSup[0].first << "\t" << config.noiseCounts.covMinusSup[0].second;
        for (unsigned i = 1; i < config.noiseCounts.covMinusSup.size(); ++i)
        {
            outFile << "\t" << config.noiseCounts.covMinusSup[i].first << "\t" << config.noiseCounts.covMinusSup[i].second;
        }
        outFile << "\n";
    }
}

void writeIndex(Config<SampleTree> const & config)
{
    if (config.saveName != "")
    {
        writeTree(config);
        writeNucInfo(config);
    }
}


#endif /* OUTPUT_H_ */
