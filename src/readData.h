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
#ifndef READDATA_H
#define READDATA_H

#include <boost/algorithm/string.hpp>
#include <boost/math/distributions/chi_squared.hpp>

#include <cctype>
#include <set>

#include "sciphi_config.h"
#include "probabilities.h"
#include "noise_counts.h"

#include <dlib/global_optimization.h>

// This functions is used to extract information about 
// sequencing errors. Only if a single cell shows a mutation
// the error rates will be collected.
inline
void updateSeqErrorStats(unsigned  & seqErrors,
    unsigned & seqErrorsCombCov,
    std::vector<std::array<unsigned, 5>> const & counts)
{

    for (unsigned short j = 0; j < 4; ++j) // loop over the nucleotides
    {
        for (unsigned i = 0; i < counts.size(); ++i) // loop over all cells
        {
            seqErrors += counts[i][j]; 
        }
    }
    for (unsigned i = 0; i < counts.size(); ++i) // update the coverage
    {
        seqErrorsCombCov += counts[i][4];
    }
}

inline
unsigned charToIndex(char c)
{
    switch (std::toupper(c))
    {
        case('A'):
            return 0;
        case('C'):
            return 1;
        case('G'):
            return 2;
        case('T'):
            return 3;
        case('-'):
        case('+'):
            return 4;
        case('^'):
            return 5;
        case('$'):
            return 6;
        case('.'):
        case(','):
            return 7;
        case('*'):
            return 8;
    }

    if (('A' <= c && c <= 'Z') || ('a' <= c && c <= 'z'))
    {
        return 9;
    }

    std::cerr << "Unknown character \"" << c << "\" in pileup sequence!" << std::endl;
    std::exit(EXIT_FAILURE);

    return 10;
}

inline
char indexToChar(unsigned index)
{
    switch (index)
    {
        case(0):
            return 'A';
        case(1):
            return 'C';
        case(2):
            return 'G';
        case(3):
            return 'T';
    }
    return 'N';
}

// This function is used to skip indels in the pileup
unsigned skipIndels(std::string const & nucleotides,
                    unsigned currentPos)
{
    if(nucleotides[currentPos] != '-' && nucleotides[currentPos] != '+')
    {
        return currentPos;
    }

    unsigned numIndels = 0;
    unsigned i = currentPos + 1; // skip the '-' or '+'
    for (; i < nucleotides.size(); ++i)
    {
        if (nucleotides[i] >= '0' && nucleotides[i] <= '9')
        {
            numIndels *= 10;
            numIndels += static_cast<int>(nucleotides[i]) - 48;
        }
        else
        {
            break;
        }
    }
    return i + (numIndels - 1);
}

inline 
void extractCellNucCountInformation(std::array<unsigned, 5> & counts,
        std::string const & nucleotides)
{
    //counts = {{0, 0, 0, 0, 0}}; // (a,c,g,t, coverage)
    for (size_t j = 0; j < nucleotides.size(); ++j) // loop over the nucleotides of a cell in the pileup
    {
        unsigned index = charToIndex(nucleotides[j]);
        if (index < 4)
        {
            ++counts[index];
        }
        else if (index == 9)
        {
            --counts[4];
        }
        else if (index == 4)
        {
            j = skipIndels(nucleotides, j);
        }
        else if (index == 5)
        {
            ++j;
        }
        else if (index == 6)
        {
            continue;
        }

    }
}

// This function is used to collect the sequencing information of all cells
inline 
void extractSeqInformation(std::array<unsigned, 5> & count,
        std::vector<std::string> const & splitVec,
        unsigned position)
{
        count = {{0,0,0,0,0}}; // (a,c,g,t, coverage)
        count[4] = std::stoi(splitVec[position * 3 + 3]);
        extractCellNucCountInformation(count, splitVec[position * 3 + 4]);
}

// This function is used to collect the sequencing information of all cells
inline 
void extractSeqInformation(std::vector<std::array<unsigned, 5>> & counts,
        std::vector<std::string> const & splitVec,
        std::vector<unsigned> const & positions)
{
    for (unsigned i = 0; i < positions.size(); ++i) // loop over the cells
    {
        extractSeqInformation(counts[i], splitVec, positions[i]);
    }
}

// This function updates an array with sequencing error/noies information.
// The array is later used when computing the tree score to account for all
// the non mutated positions.
inline 
void addNoiseCounts(std::vector<std::array<unsigned, 5>> & counts,
        GappedNoiseCounts & gappedNoiseCounts)
{
    for (unsigned cellId = 0; cellId < counts.size(); ++cellId) // loop over the cells
    {
        unsigned numErrors = 0;
        for (unsigned nucId = 0; nucId < 4; ++nucId) // check all 4 different nucleotides
        {
            numErrors += counts[cellId][nucId];
        }
        gappedNoiseCounts.add(numErrors, counts[cellId][4]);
    }
}

template <typename TTreeType> 
void applyCoverageFilterPerCell(std::vector<std::array<unsigned, 5>> & counts,
        Config<TTreeType> const & config)
{
    for (unsigned i = 0; i < counts.size(); ++i)
    {
        if(counts[i][4] < config.minCoverage)
        {
            for (unsigned j = 0; j < 5; ++j)
            {
                counts[i][j] = 0;
            }
        }
    }
}

bool passSuppFilter(unsigned altCount, 
                    unsigned minSupport)
{
    if (altCount >= minSupport)
    {
        return true;
    }
    return false;
}

bool passFreqFilter(double altCount, 
                    double coverage, 
                    double minFreq)
{
    if (altCount/coverage >= minFreq)
    {
        return true;
    }
    return false;
}

bool passCovFilter(unsigned coverage, 
                    unsigned minCov)
{
    if (coverage >= minCov)
    {
        return true;
    }
    return false;
}

template <typename TTreeType> 
bool applyFilterAcrossCells(std::vector<std::array<unsigned, 5>> & counts,
        Config<TTreeType> const & config,
        unsigned nucId)
{
    unsigned numCellsAboveThresh = 0;
    for (unsigned i = 0; i < counts.size(); ++i)
    {
        if (passCovFilter(counts[i][4], config.minCoverage) &&
                passFreqFilter(counts[i][nucId],counts[i][4], config.minFreq) &&
                passSuppFilter(counts[i][nucId], config.minSupport))
        {
            ++numCellsAboveThresh;
            if(numCellsAboveThresh >= config.numMinCoverageAcrossCells)
            {
                return true;
            }
        }
    }
    return false;
}


inline
void 
readExclusionList(std::set<std::tuple<std::string, std::string>> & exMap,
                  std::string const & fileName)
{
    std::ifstream inputStream(fileName);
    std::string currLine;                    // line currently processed
    std::vector<std::string> splitVec;  // vector storing the words in currLine dived by a tab
    while (std::getline(inputStream, currLine))
    {
        if (currLine[0] != '#')
        {
            boost::split(splitVec, currLine, boost::is_any_of("\t"));
            exMap.insert(std::tuple<std::string, std::string>(splitVec[0], std::to_string(std::stoi(splitVec[1]))));
        }
    }
}

template <typename TTreeType>
void
readCellInformation(std::vector<unsigned> & tumorCellPos, std::vector<unsigned> & normalCellPos, unsigned & tumorBulkPos, unsigned & normalBulkPos, Config<TTreeType> & config)
{

    tumorCellPos.resize(0);
    normalCellPos.resize(0);
    tumorBulkPos = UINT_MAX;
    normalBulkPos = UINT_MAX;

    std::ifstream inputStream(config.bamFileNames);

    std::vector<std::string> splitVec;
    std::vector<std::string> splitVecEntry;

    unsigned counter = 0;
    std::string currLine;

    while(getline(inputStream, currLine))
    {
        if (currLine != "")
        {
            boost::split(splitVec, currLine, boost::is_any_of("\t"));
            if (splitVec.back() == "BN")
            {
                if (normalBulkPos != UINT_MAX)
                {
                    std::cout << "WARNING: Multiple bulk control files provided. Only using the first one!";
                }
                else
                {
                    normalBulkPos = counter;
                }
            }
            if (splitVec.back() == "BT")
            {
            }
            if (splitVec.back() == "CN")
            {
                normalCellPos.push_back(counter);
            }
            if (splitVec.back() == "CT" || (config.useNormalCellInTree && splitVec.back() == "CN"))
            {
                tumorCellPos.push_back(counter);
                boost::split(splitVecEntry, splitVec[0], boost::is_any_of("/"));
                config.cellNames.push_back(splitVecEntry.back());
            }
            ++counter;
        }
    }
}


inline 
double logNChoose2(unsigned numMut)
{
    return log(static_cast<double>(numMut)) + log((static_cast<double>(numMut) - 1.0) / 2.0);
}

inline 
double logNChooseK(unsigned n, unsigned k, double logNChoosekMinusOne)
{
    if (k == 0)
        return 0;
    return logNChoosekMinusOne + log(static_cast<double>(n + 1 - k) / static_cast<double>(k));
}

inline 
double logNChooseK(unsigned n, unsigned k)
{
    if (k == 0)
        return 0;
    return std::lgamma(n+1) - std::lgamma(k+1) - std::lgamma(n - k + 1);
}

inline 
double sumValuesInLogSpace(std::vector<double>::const_iterator itBegin, std::vector<double>::const_iterator itEnd)
{
    std::vector<double>::const_iterator it = itBegin;
    double maxLogValue = *it;
    ++it;
    for (;it != itEnd; ++it)
    {
        if (maxLogValue < *it)
        {
            maxLogValue = *it;
        }
    }

    it = itBegin;
    double h1 = 0;
    for (;it != itEnd; ++it)
    {
        h1 += exp(*it - maxLogValue);
    }

    return log(h1) + maxLogValue;

}

template<typename TTreeType>
void estimateSeqErrorRate(Config<TTreeType> & config,
    std::set<std::tuple<std::string, std::string>> const & exMap,
    std::set<std::tuple<std::string, std::string>> const & errExMap,
    std::vector<unsigned> const & tumorCellPos,
    std::vector<unsigned> const & normalCellPos)
{
    //estimate the error rate
    std::ifstream inputStream(config.inFileName);
    if (inputStream.fail())
    {
        std::cerr << "The pileup file provided does not exist or you lack permission to access it." << std::endl;
        return;
    }
    std::string currLine;
    std::vector<std::string> splitVec;
    std::vector<std::array<unsigned, 5>> tumor_counts(tumorCellPos.size(), {{0,0,0,0,0}});
    std::vector<std::array<unsigned, 5>> normal_counts(normalCellPos.size(), {{0,0,0,0,0}});
    unsigned seqErrors = 0;
    unsigned seqErrorsCombCov = 0;
    for (size_t lineNumber = 0; lineNumber < config.errorRateEstLoops && getline(inputStream, currLine); ++lineNumber)
    {
        boost::split(splitVec, currLine, boost::is_any_of("\t"));
        auto itEx = exMap.find(std::tuple<std::string, std::string>(splitVec[0], splitVec[1]));
        auto itErrEx = errExMap.find(std::tuple<std::string, std::string>(splitVec[0], splitVec[1]));
        if (itEx == exMap.end() && itErrEx == errExMap.end())
        {
            extractSeqInformation(tumor_counts, splitVec, tumorCellPos);
            updateSeqErrorStats(seqErrors, seqErrorsCombCov, tumor_counts);
            extractSeqInformation(normal_counts, splitVec, normalCellPos);
            updateSeqErrorStats(seqErrors, seqErrorsCombCov, normal_counts);
        }
        else
        {
            --lineNumber;
        }
    }
    ++seqErrors; // add pseudo counts
    inputStream.close();
    
    // set the new estimated error rate

    if(seqErrorsCombCov > 0)
    {
        config.setParam(Config<TTreeType>::wildMean, static_cast<double>(seqErrors)/static_cast<double>(seqErrorsCombCov));
        config.setTmpParam(Config<TTreeType>::wildMean, config.getParam(Config<TTreeType>::wildMean));
    }
}   

void printCurrentChrom(std::string & currentChrom, std::string const & testChrom)
{
    if (testChrom != currentChrom)
    {
        currentChrom = testChrom;
        std::cout << currentChrom << std::endl << std::flush;
    }
}

double updateLogH1Temp(double logH1Temp, unsigned numCells, unsigned numMuts)
{
    return logH1Temp + 2 * logNChooseK(numCells, numMuts) - std::log(2*numMuts - 1) - logNChooseK(2*numCells, 2*numMuts);
}

bool mustH0Win(double & logH1Max, double logH1Temp, double logNumCells, double logH0)
{
    if (logH1Temp >= logH1Max)
    {
        logH1Max = logH1Temp;
        return false;
    }
    else
    {
        if (logH1Max + logNumCells < logH0)
        {
            return true;
        }
    }
    return false;
}

template <typename TTreeType, typename TData>
void insertData(Config<TTreeType> & config, TData const & data)
{
    for (size_t i = 0; i < data.size(); ++i)
    {
        config.indexToPosition.push_back(std::get<0>(data[i]));
        for (size_t j = 0; j < std::get<1>(data[i]).size(); ++j)
        {
            config.getCompleteData()[j].push_back(std::get<1>(data[i])[j]);
        }
    }
}

template <typename TTreeType>
bool
passNormalFilter(std::array<unsigned, 5> const & normalCounts, Config<TTreeType> const & config)
{
    if (normalCounts[4] >= config.minCovInControlBulk)
    {
        for (unsigned i = 0; i < 4; ++i)
        {
            if (normalCounts[i] >= config.maxSupInControlBulk)
            {
                return false;
            }
        }
        return true;
    }
    return false;
}

template <typename TTreeType>
bool
passNormalCellFilter(std::vector<std::array<unsigned, 5>> const & normalCellCounts, Config<TTreeType> const & config)
{
    unsigned maxCov = 0;
    unsigned numMutated = 0;
    for (size_t i = 0; i < normalCellCounts.size(); ++i)
    {
        if (normalCellCounts[i][4] > maxCov)
        {
            maxCov = normalCellCounts[i][4];
        }
        for (size_t j = 0; j < 4; ++j)
        {
            if (computeWildLogScore(config, normalCellCounts[i][j], normalCellCounts[i][4]) < computeRawMutLogScore(config, normalCellCounts[i][j], normalCellCounts[i][4]))
            {
                ++numMutated;
            }
        }
    }
    
    if (maxCov >= config.minCovNormalCell && numMutated <= config.maxNumberNormalCellMutated)
    {
        return true;
    }
    return false;
}

template <typename TTreeType>
bool computeProbCellsAreMutated(
    Config<TTreeType> const & config,
    std::vector<long double> & logProbs,
    std::vector<long double> & tempLogProbs,
    std::vector<double> & logProbTempValues,
    std::vector<std::array<unsigned, 5>> & filteredCounts,
    std::vector<double> & cellsNotMutated,
    std::vector<double> & cellsMutated,
    unsigned currentChar)
{

    unsigned numCells = config.getNumSamples();
    double logNumCells = log(numCells);

    tempLogProbs[0] = 0.0;
    for (size_t i = 1; i < filteredCounts.size() + 1; ++i)
    {
        cellsNotMutated[i-1] = computeWildLogScore(config, filteredCounts[i - 1][currentChar], filteredCounts[i - 1][4]);
        cellsMutated[i-1] = computeRawMutLogScore(config, filteredCounts[i - 1][currentChar], filteredCounts[i - 1][4]);
        tempLogProbs[i] = tempLogProbs[i - 1] + cellsNotMutated[i - 1];

    }
    swap(logProbs, tempLogProbs);
    double logH0 = logProbs.back() + log(1.0 - config.priorMutationRate); 
    logProbTempValues[0] = logH0;

    // compute the probabilitues of observing 1, 2, 3, ... mutations
    double logH1Max = -DBL_MAX; // the current best alternative score
    long double logNOverK = 0;  // helper to efficiently compute nChooseK
    size_t numMut = 1;          // number of mutations currently computet
    
    for (; numMut <= numCells; ++numMut)
    {
        double logProbAllPrevCellsMutated = logProbs[numMut - 1];
        double currentCellMutated = cellsMutated[numMut - 1];
        tempLogProbs[numMut] = logProbAllPrevCellsMutated + currentCellMutated;
        for (size_t i = numMut + 1; i < filteredCounts.size() + 1; ++i)
        {
            double previousCellNotMutated = logProbs[i - 1];
            currentCellMutated = cellsMutated[i - 1]; 
            double previousCellMutated = tempLogProbs[i -1];
            double currentCellNotMutated = cellsNotMutated[i - 1]; 
            tempLogProbs[i] = addLogProb(previousCellNotMutated + currentCellMutated,  
                        previousCellMutated + currentCellNotMutated);

        }
        swap(logProbs, tempLogProbs);
        logNOverK = logNChooseK(numCells, numMut, logNOverK);
        double logH1Temp = logProbs.back() + log(config.priorMutationRate) - logNOverK;
        logH1Temp = updateLogH1Temp(logH1Temp, numCells, numMut);

        // check whether the alternative hyothesis can win
        bool h0Wins = mustH0Win(logH1Max, logH1Temp, logNumCells, logH0);
        logProbTempValues[numMut] = logH1Temp;
        if (h0Wins)
        { 
            return true;
        }
    }
    return false;
}


template <typename TTreeType>
bool readMpileupFile(Config<TTreeType> & config)
{
    typedef std::tuple<unsigned, unsigned>                                  TCountType;
    typedef std::tuple<std::string, unsigned, char, char>                   TPositionInfo;
    typedef std::tuple<TPositionInfo, std::vector<TCountType>, long double> TDataEntry;
    typedef std::vector<TDataEntry>                                         TData;
    
    // read the cellular information
    std::vector<unsigned> tumorCellPos;
    std::vector<unsigned> normalCellPos;
    unsigned tumorBulkPos; 
    unsigned normalBulkPos;
    readCellInformation(tumorCellPos, normalCellPos, tumorBulkPos, normalBulkPos, config);

    // the variables for number of samples
    unsigned numCells = tumorCellPos.size();
    config.setNumSamples(numCells);
    
    // resize the container holding the nucleotide counts of likely mutations
    config.getCompleteData().resize(numCells);
    config.getData().resize(numCells);
    double logNumCells = log(numCells);
    
    // determine which positions to skip
    std::set<std::tuple<std::string, std::string>> exMap;
    readExclusionList(exMap, config.exclusionFileName);
    
    // determine which positions to skip during the error rate estimation
    std::set<std::tuple<std::string, std::string>> errExMap;
    readExclusionList(errExMap, config.mutationExclusionFileName);
    
    TData data;

    std::array<unsigned, 5> normalBulkCounts;
    std::vector<std::array<unsigned, 5>> counts(numCells, {{0,0,0,0,0}});           // vector to hold the nucleotide information {a,c,g,t,coverage}
    std::vector<std::array<unsigned, 5>> countsNormal(normalCellPos.size(), {{0,0,0,0,0}});           // vector to hold the nucleotide information {a,c,g,t,coverage}
    std::vector<std::array<unsigned, 5>> filteredCounts(numCells, {{0,0,0,0,0}});   // vector to hold the filtered nucleotide information {a,c,g,t,coverage}
    std::vector<TCountType> tempCounts(numCells);                                   // helper vector to store nucleotide information
    std::vector<long double> logProbs(numCells + 1, 0);                             // probabilities of observing 0, 1, 2, 3, 4 ... mutations
    std::vector<long double> tempLogProbs(numCells + 1, 0);                         // helper array for probabilities of observing 0, 1, 2, 3, 4 ... mutations
    std::vector<double> logProbTempValues(numCells + 1);

    std::vector<double> cellsNotMutated(numCells);
    std::vector<double> cellsMutated(numCells);


    if (config.estimateSeqErrorRate)
    {
        estimateSeqErrorRate(config, exMap, errExMap, tumorCellPos, normalCellPos);
    }

    std::ofstream ofs(config.outFilePrefix + "_posMutated.tsv");

    std::ifstream inputStream(config.inFileName, std::ifstream::in);
    std::string currentChrom = "";
    std::string currLine;
    std::vector<std::string> splitVec;

    GappedNoiseCounts gappedNoiseCounts;
    unsigned counter = 0;

    while (getline(inputStream, currLine))
    {
        // solit the current line into easily accessible chunks
        boost::split(splitVec, currLine, boost::is_any_of("\t"));
        // print the progress
        printCurrentChrom(currentChrom, splitVec[0]);

        // check if the current pos is to be excluded
        auto it = exMap.find(std::tuple<std::string, std::string>(splitVec[0], splitVec[1]));
        if (it == exMap.end())
        {
            if (normalBulkPos != UINT_MAX)
            {
                extractSeqInformation(normalBulkCounts, splitVec, normalBulkPos);
                if (!passNormalFilter(normalBulkCounts, config))
                {
                    continue;
                }
            }

            extractSeqInformation(countsNormal, splitVec, normalCellPos);
            if (!passNormalCellFilter(countsNormal, config))
            {
                continue;
            }

            extractSeqInformation(counts, splitVec, tumorCellPos);

            filteredCounts = counts;
            applyCoverageFilterPerCell(filteredCounts, config);
            //applySuppFreqFilterPerCell(filteredCounts, config);

            bool positionMutated = false;
            for (unsigned short j = 0; j < 4; ++j)
            {
                if(j == charToIndex(splitVec[2][0]))
                {
                    continue;
                }
                bool h0Wins = !applyFilterAcrossCells(filteredCounts, config, j);

                double logH0 = -DBL_MAX;
                double logH1 = -DBL_MAX;
                if(!h0Wins)
                {
                    h0Wins = computeProbCellsAreMutated(config, 
                            logProbs, 
                            tempLogProbs, 
                            logProbTempValues, 
                            filteredCounts,
                            cellsNotMutated,
                            cellsMutated,
                            j);
                }

                if (h0Wins)
                {
                    logH1 = -DBL_MAX;
                    logH0 = DBL_MAX;
                }
                else
                {
                    logH0 = sumValuesInLogSpace(logProbTempValues.begin(),logProbTempValues.begin() + config.numCellWithMutationMin);
                    logH1 = sumValuesInLogSpace(logProbTempValues.begin() + config.numCellWithMutationMin, logProbTempValues.end());
                }
                

                if (logH1 > logH0)
                {
                    positionMutated = true;

                    static std::vector<std::pair<unsigned, unsigned>> testCounts;
                    testCounts.resize(0);
                    for (size_t cell = 0; cell < counts.size(); ++cell)
                    {
                        if(cellsNotMutated[cell] < cellsMutated[cell])
                        {
                            testCounts.push_back({counts[cell][j], counts[cell][4]});
                        }
                    }

                    dlib::matrix<double,0,1> startingPointMeanOverDis = {0.05,5.0};
                    dlib::matrix<double,0,1> dLibMinMeanOverDis = {0.001,0.1};
                    dlib::matrix<double,0,1> dLibMaxMeanOverDis = {0.999,10000.0};
                    OptimizeBetaBinMeanOverDis optBetaBinMeanOverDis(testCounts);
                    OptimizeBetaBinMeanOverDisDerivates optBetaBinDerMeanOverDis(testCounts);
                    double resultMeanOverDis = dlib::find_max_box_constrained(dlib::bfgs_search_strategy(),  // Use BFGS search algorithm
                        dlib::objective_delta_stop_strategy(1e-7), // Stop when the change in rosen() is less than 1e-7
                        optBetaBinMeanOverDis, 
                        optBetaBinDerMeanOverDis, 
                        startingPointMeanOverDis, 
                        dLibMinMeanOverDis,
                        dLibMaxMeanOverDis);
                    
                    dlib::matrix<double,0,1> startingPointOverDis = {2.0};
                    dlib::matrix<double,0,1> dLibMinOverDis = {0.1};
                    dlib::matrix<double,0,1> dLibMaxOverDis = {10000.0};
                    OptimizeBetaBinOverDis optBetaBinOverDis(testCounts, config.meanFilter);
                    OptimizeBetaBinOverDisDerivates optBetaBinDerOverDis(testCounts, config.meanFilter);
                    double resultOverDis = dlib::find_max_box_constrained(dlib::bfgs_search_strategy(),  // Use BFGS search algorithm
                        dlib::objective_delta_stop_strategy(1e-7), // Stop when the change in rosen() is less than 1e-7
                        optBetaBinOverDis, 
                        optBetaBinDerOverDis, 
                        startingPointOverDis, 
                        dLibMinOverDis,
                        dLibMaxOverDis);


                    double pValue;
                    if (resultOverDis > resultMeanOverDis) // if true the results are within the optimization limit
                    {
                        pValue = 1;
                    }
                    else
                    {
                        double testStat = -2 * (resultOverDis - resultMeanOverDis);
                        boost::math::chi_squared mydist(1);
                        pValue = 1 - boost::math::cdf(mydist, testStat);
                    }

                    if(pValue > 0.05 || startingPointMeanOverDis >= config.meanFilter)
                    {
                        //std::cout << "HIT" << std::endl;
                        unsigned numAffectetCells = 0;
                        for (size_t cell = 0; cell < counts.size(); ++cell)
                        {
                            std::get<0>(tempCounts[cell]) = counts[cell][4];
                            std::get<1>(tempCounts[cell]) = counts[cell][j];
                            if (counts[cell][j] > 0)
                            {
                                ++numAffectetCells;
                            }
                        }
                
                        if (numAffectetCells > 1)
                        {
                          data.push_back(TDataEntry(TPositionInfo(splitVec[0], std::stoi(splitVec[1]), splitVec[2][0], indexToChar(j)), tempCounts, logH1/logH0));
                          ofs << splitVec[0] << "\t" << splitVec[1] << "\t" <<  splitVec[2][0] << "\t" <<  indexToChar(j) << "\t";
                          for (auto it = tempCounts.begin(); it != tempCounts.end(); ++it)
                              ofs << std::get<0>(*it) << "\t" << std::get<1>(*it) << "\t";
                          ofs<< "\n";
                        }
                        else
                        {
                            data.push_back(TDataEntry(TPositionInfo(splitVec[0], std::stoi(splitVec[1]), splitVec[2][0], indexToChar(j)), tempCounts, -DBL_MAX));
                            ++config.numUniqMuts;
                        }
                    }
                }
            }
            if(!positionMutated)
            {
                addNoiseCounts(counts, gappedNoiseCounts);
            }
        }
        ++counter;
    }

    config.noiseCounts = NoiseCounts(gappedNoiseCounts);
    computeNoiseScore(config);

    //std::sort(begin(data), end(data), 
    //    [](TDataEntry const &t1, TDataEntry const &t2) {return get<2>(t1) > get<2>(t2);}
    //);


    insertData(config, data);
    return 0;
}


template <typename TTreeType>
void getData(Config<TTreeType> & config)
{
    readMpileupFile(config);
}
#endif

