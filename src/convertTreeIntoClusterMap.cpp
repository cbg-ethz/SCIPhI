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


#include <boost/multi_array.hpp>
#include <boost/program_options.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/graph/adjacency_list.hpp>

struct Files
{
    std::string inFile;
    std::string inFileProbs;
    std::string outFile;
};

struct DotVertex {
    std::string name;
    std::string label;
    int peripheries;
};

struct DotEdge {
    std::string label;
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
        DotVertex, DotEdge> graph_t;

class PrintTreeAsMap : public boost::default_dfs_visitor
{
    boost::dynamic_properties & dp;
    std::vector<std::string> & muts;
    std::set<std::string> & mutsSet;
    std::vector<std::string> & cells;
    std::set<std::string> * cellsSet;

public:

    PrintTreeAsMap(boost::dynamic_properties & dp_, 
            std::vector<std::string> & muts_, 
            std::set<std::string> & mutsSet_, 
            std::vector<std::string> & cells_) : 
        dp(dp_),
        muts(muts_),
        mutsSet(mutsSet_),
        cells(cells_) 
    {}
    
    template <typename TVertex >
    void discover_vertex(TVertex v, graph_t const & g) const
    {
        std::string label = boost::get("label", dp, v);
        std::vector<std::string> splitVec;
        boost::split(splitVec, label, boost::is_any_of(":\n "), boost::token_compress_on);
        if(out_degree(v, g) != 0)
        {
            for (unsigned i = 0; i < splitVec.size() - 1; ++i)
            {
                auto it = mutsSet.find(splitVec[i]);
                if (it == mutsSet.end())
                {
                    muts.push_back(splitVec[i]);
                    mutsSet.insert(splitVec[i]);
                }
            }
        }
        else
        {
            for (unsigned i = 1; i < splitVec.size() - 1; ++i)
            {
                auto it = mutsSet.find(splitVec[i]);
                if (it == mutsSet.end())
                {
                    muts.push_back(splitVec[i]);
                    mutsSet.insert(splitVec[i]);
                }
            }
        }
        return ;
    }

    template <typename TVertex >
    void finish_vertex(TVertex v, graph_t const & g) const
    {
        std::string label = boost::get("label", dp, v);
        if(out_degree(v, g) == 0)
        {
            std::vector<std::string> splitVec;
            boost::split(splitVec, label, boost::is_any_of(" :\n"));
            cells.push_back(splitVec[0]);
        }
        return ;
    }
};


int readParameters(Files & files, int argc, char* argv[]){

    boost::program_options::options_description generic("Generic options");
	generic.add_options()("help,h", "Print this help");

	// configuration options
	boost::program_options::options_description parseConfig("Configuration");
	parseConfig.add_options()
		(",i", boost::program_options::value<decltype(files.inFile)>(&files.inFile), "The tree to convert.")
		(",p", boost::program_options::value<decltype(files.inFileProbs)>(&files.inFileProbs), "The tree to convert.")
		(",o", boost::program_options::value<decltype(files.outFile)>(&files.outFile), "The map file to write.");

	boost::program_options::options_description cmdline_options;
	cmdline_options.add(generic).add(parseConfig);
	boost::program_options::options_description visible("Allowed options");
	visible.add(generic).add(parseConfig);
	boost::program_options::variables_map global_options;

	/* 0.0) parse program options */
	try
	{
		boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).run(), global_options);

		// show help options
		if (global_options.count("help"))
		{
			std::cout << visible << '\n';
			exit(EXIT_SUCCESS);
		}

		boost::program_options::notify(global_options);
	}
	catch (boost::program_options::error& e)
	{
		std::cerr << "ERROR: " << e.what() << '\n';
		exit(EXIT_FAILURE);
	} 

    return 0;
}

int main(int argc, char* argv[])
{

    Files files;
    readParameters(files, argc, argv);

    graph_t tree;
    boost::dynamic_properties dp(boost::ignore_other_properties);

    dp.property("node_id",     boost::get(&DotVertex::name,        tree));
    dp.property("label",       boost::get(&DotVertex::label,       tree));
    dp.property("peripheries", boost::get(&DotVertex::peripheries, tree));
    dp.property("label",       boost::get(&DotEdge::label,         tree));

    std::ifstream ifs(files.inFile);
    bool status = boost::read_graphviz(ifs, tree, dp);

    std::vector<std::string> muts;
    std::set<std::string> mutsSet;
    std::vector<std::string> cells;
    PrintTreeAsMap vis(dp, muts, mutsSet, cells);
    
    depth_first_search(tree, visitor(vis).root_vertex(0));
    mutsSet.clear();

    std::unordered_map<std::string,unsigned> mutMap;
    for (unsigned i = 0; i < muts.size(); ++i)
    {
        mutMap[muts[i]] = i;
    }
    
    std::unordered_map<std::string,unsigned> cellsMap;
    for (unsigned i = 0; i < cells.size(); ++i)
    {
        cellsMap[cells[i]] = i;
    }

    std::vector<std::pair<std::string, unsigned>> cellOrder(cellsMap.begin(), cellsMap.end());
    std::sort(cellOrder.begin(), cellOrder.end(), [](auto &left, auto &right) {return left.second < right.second;});
    std::vector<std::pair<std::string, unsigned>> mutOrder(mutMap.begin(), mutMap.end());
    std::sort(mutOrder.begin(), mutOrder.end(), [](auto &left, auto &right) {return left.second < right.second;});

    std::vector<std::vector<long double>> structuredProbs;
    structuredProbs.resize(cells.size());
    for (unsigned i = 0; i < structuredProbs.size(); ++i)
    {
        structuredProbs[i].resize(muts.size(), 0);
    }

    std::ifstream file(files.inFileProbs);
    std::string   line;
    std::getline(file, line);  // read up-to the first tab (discard tab).
    std::vector<std::string> cellNames;
    boost::split(cellNames, line, boost::is_any_of("\t"), boost::token_compress_on);
    while(std::getline(file, line))
    {
        std::vector<std::string> splitVec;
        boost::split(splitVec, line, boost::is_any_of("\t"), boost::token_compress_on);

        std::string id = splitVec[0] + "_" + splitVec[1];

        auto it = mutsSet.find(id);
        if (it == mutsSet.end())
        {
            for (unsigned i = 2; i < splitVec.size()-1; ++i)
            {
                std::vector<std::string> splitEntry;
                boost::split(splitEntry, splitVec[i], boost::is_any_of("|"), boost::token_compress_on);
                structuredProbs[cellsMap[cellNames[i]]][mutMap[id]] = std::stod(splitEntry.back());
                //std::cout << "i: " << i << " cellNames[" << i << "] " << cellNames[i] << " id: " << id << " mutMap[" << id << "] " << mutMap[id] << " "  << structuredProbs[cellsMap[cellNames[i]]][mutMap[id]] << std::endl;
                mutsSet.insert(id);
            }
        }
    }

    std::ofstream ofFile(files.outFile);
    ofFile << "cellName\t";
    for (unsigned i = 0; i < cellOrder.size(); ++i)
    {
        ofFile << cellOrder[i].first << "\t";
    }
     ofFile << std::endl;
    for (unsigned i = 0; i < muts.size(); ++i)
    {
        ofFile << mutOrder[i].first << "\t";
        for (unsigned j = 0; j < cells.size(); ++j)
        {
            //std::cout << j << " " << structuredProbs.size() <<  " " << i << " " << structuredProbs[j].size() << std::endl;
            ofFile << structuredProbs[j][i] << "\t";
        }
        ofFile << std::endl;
    }
    ofFile.close();
}
