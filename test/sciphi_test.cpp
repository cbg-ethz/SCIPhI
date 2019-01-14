#define BOOST_TEST_MODULE ScateTest
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "src/sciphi_config.h"
#include "src/mcmc.h"
#include "src/scoreTree.h"
#include "src/readData.h"

BOOST_AUTO_TEST_CASE( sample_tree_init )
{
    Config<SampleTree> config;
    config.setNumSamples(3);

    createInitialTree(config);

    BOOST_CHECK_MESSAGE(num_vertices(config.getTree()) == 6, num_vertices(config.getTree()) << " != 5");

    BOOST_CHECK(source(*in_edges(1, config.getTree()).first, config.getTree()) == 0);
    BOOST_CHECK(source(*in_edges(2, config.getTree()).first, config.getTree()) == 0);
    BOOST_CHECK(source(*in_edges(3, config.getTree()).first, config.getTree()) == 1);
    BOOST_CHECK(source(*in_edges(4, config.getTree()).first, config.getTree()) == 1);

    BOOST_CHECK(config.getTree()[0].sample == -1);
    BOOST_CHECK(config.getTree()[1].sample == -1);
    BOOST_CHECK(config.getTree()[2].sample == 0);
    BOOST_CHECK(config.getTree()[3].sample == 1);
    BOOST_CHECK(config.getTree()[4].sample == 2);
    
    // artificial root
    BOOST_CHECK_MESSAGE(source(*in_edges(0, config.getTree()).first, config.getTree()) == 5, source(*in_edges(0, config.getTree()).first, config.getTree()) << " != 5");
}

BOOST_AUTO_TEST_CASE( get_sibling )
{

    for (unsigned i = 0; i < 10000; ++i)
    {
        Config<SampleTree> config;
        config.setNumSamples(7);

        createInitialTree(config);

        unsigned newNonDescendant = getNewSibling(config, 1);

        BOOST_CHECK_MESSAGE(newNonDescendant != 0, newNonDescendant << " == 0");
        BOOST_CHECK_MESSAGE(newNonDescendant != 1, newNonDescendant << " == 1");
        BOOST_CHECK_MESSAGE(newNonDescendant != 3, newNonDescendant  << " == 3");
        BOOST_CHECK_MESSAGE(newNonDescendant != 4, newNonDescendant  << " == 4");
        BOOST_CHECK_MESSAGE(newNonDescendant != 7, newNonDescendant  << " == 7");
        BOOST_CHECK_MESSAGE(newNonDescendant != 8, newNonDescendant  << " == 8");
        BOOST_CHECK_MESSAGE(newNonDescendant != 9, newNonDescendant  << " == 9");
        BOOST_CHECK_MESSAGE(newNonDescendant != 10, newNonDescendant  << " == 10");
        BOOST_CHECK_MESSAGE(newNonDescendant != 13, newNonDescendant  << " == 13");
    }
}

BOOST_AUTO_TEST_CASE( sample_tree_prune_and_reattach )
{

    std::vector<unsigned> counter = {0, 0, 0, 0, 0};
    for (unsigned i = 0; i < 10000; ++i)
    {
        bool hit = false;
        Config<SampleTree> config;
        config.setNumSamples(3);
        createInitialTree(config);

        pruneAndReAttach(config);

        if ( source(*in_edges(0, config.getTree()).first, config.getTree()) == 1 &&
             source(*in_edges(1, config.getTree()).first, config.getTree()) == 5 &&
             source(*in_edges(2, config.getTree()).first, config.getTree()) == 0 &&
             source(*in_edges(3, config.getTree()).first, config.getTree()) == 0 &&
             source(*in_edges(4, config.getTree()).first, config.getTree()) == 1)
        {
            hit =true;
            ++counter[0];
        }
        else if ( source(*in_edges(0, config.getTree()).first, config.getTree()) == 5 &&
             source(*in_edges(1, config.getTree()).first, config.getTree()) == 0 &&
             source(*in_edges(2, config.getTree()).first, config.getTree()) == 1 &&
             source(*in_edges(3, config.getTree()).first, config.getTree()) == 0 &&
             source(*in_edges(4, config.getTree()).first, config.getTree()) == 1)
        {
            hit =true;
            ++counter[1];
        }
        else if ( source(*in_edges(0, config.getTree()).first, config.getTree()) == 5 &&
             source(*in_edges(1, config.getTree()).first, config.getTree()) == 0 &&
             source(*in_edges(2, config.getTree()).first, config.getTree()) == 0 &&
             source(*in_edges(3, config.getTree()).first, config.getTree()) == 1 &&
             source(*in_edges(4, config.getTree()).first, config.getTree()) == 1)
        {
            hit =true;
            ++counter[2];
        }
        else if ( source(*in_edges(0, config.getTree()).first, config.getTree()) == 5 &&
             source(*in_edges(1, config.getTree()).first, config.getTree()) == 0 &&
             source(*in_edges(2, config.getTree()).first, config.getTree()) == 1 &&
             source(*in_edges(3, config.getTree()).first, config.getTree()) == 1 &&
             source(*in_edges(4, config.getTree()).first, config.getTree()) == 0)
        {
            hit =true;
            ++counter[3];
        }
        else if ( source(*in_edges(0, config.getTree()).first, config.getTree()) == 1 &&
             source(*in_edges(1, config.getTree()).first, config.getTree()) == 5 &&
             source(*in_edges(2, config.getTree()).first, config.getTree()) == 0 &&
             source(*in_edges(3, config.getTree()).first, config.getTree()) == 1 &&
             source(*in_edges(4, config.getTree()).first, config.getTree()) == 0)
        {
            hit =true;
            ++counter[4];
        }

        BOOST_CHECK_MESSAGE(hit, "Invalid tree configuration.");
        if (!hit)
            break;
    }
    BOOST_CHECK_MESSAGE(counter[0] > 0, "1 configuration never reached");
    BOOST_CHECK_MESSAGE(counter[1] > 0, "2 configuration never reached");
    BOOST_CHECK_MESSAGE(counter[2] > 0, "3 configuration never reached");
    BOOST_CHECK_MESSAGE(counter[3] > 0, "4 configuration never reached");
    BOOST_CHECK_MESSAGE(counter[4] > 0, "5 configuration never reached");

}

BOOST_AUTO_TEST_CASE( sample_tree_swap_node_labels )
{
    std::vector<unsigned> counter = {0, 0, 0};
    for (unsigned i = 0; i < 10000; ++i)
    {
        bool hit = false;
        Config<SampleTree> config;
        config.setNumSamples(3);
        createInitialTree(config);

        BOOST_CHECK(source(*in_edges(0, config.getTree()).first, config.getTree()) == 5);
        BOOST_CHECK(source(*in_edges(1, config.getTree()).first, config.getTree()) == 0);
        BOOST_CHECK(source(*in_edges(2, config.getTree()).first, config.getTree()) == 0);
        BOOST_CHECK(source(*in_edges(3, config.getTree()).first, config.getTree()) == 1);
        BOOST_CHECK(source(*in_edges(4, config.getTree()).first, config.getTree()) == 1);

        swapNodeLabels(config);

        if ( config.getTree()[0].sample == -1 &&
             config.getTree()[1].sample == -1 &&
             config.getTree()[2].sample == 1 &&
             config.getTree()[3].sample == 0 &&
             config.getTree()[4].sample == 2)
        {
            hit =true;
            ++counter[0];
        }
        else if ( config.getTree()[0].sample == -1 &&
             config.getTree()[1].sample == -1 &&
             config.getTree()[2].sample == 0 &&
             config.getTree()[3].sample == 2 &&
             config.getTree()[4].sample == 1)
        {
            hit =true;
            ++counter[1];
        }
        else if ( config.getTree()[0].sample == -1 &&
             config.getTree()[1].sample == -1 &&
             config.getTree()[2].sample == 2 &&
             config.getTree()[3].sample == 1 &&
             config.getTree()[4].sample == 0)
        {
            hit =true;
            ++counter[2];
        }

        BOOST_CHECK_MESSAGE(hit, "Invalid tree configuration.");
        if (!hit)
        {
            for (unsigned i = 0; i < num_vertices(config.getTree()); ++i)
               std::cout << "Sample of node " << i << " = " << config.getTree()[i].sample << std::endl;
            break;
        }
    }
    for (unsigned i = 0; i < counter.size(); ++i)
        BOOST_CHECK_MESSAGE(counter[i] > 0, i << ". configuration never reached");

}

BOOST_AUTO_TEST_CASE( sample_tree_get_mutation_of_nodes )
{
    Config<SampleTree> config;
    config.setNumSamples(3);
    createInitialTree(config);

    BOOST_CHECK(source(*in_edges(0, config.getTree()).first, config.getTree()) == 5);
    BOOST_CHECK(source(*in_edges(1, config.getTree()).first, config.getTree()) == 0);
    BOOST_CHECK(source(*in_edges(2, config.getTree()).first, config.getTree()) == 0);
    BOOST_CHECK(source(*in_edges(3, config.getTree()).first, config.getTree()) == 1);
    BOOST_CHECK(source(*in_edges(4, config.getTree()).first, config.getTree()) == 1);

    config.data.resize(3);
    config.data[0].resize(3, std::tuple<unsigned, unsigned>{10, 0});
    config.data[1].resize(3, std::tuple<unsigned, unsigned>{10, 0});
    config.data[2].resize(3, std::tuple<unsigned, unsigned>{10, 0});

    // cell 2 and 3 have mutation 1
    std::get<1>(config.data[1][0]) = 10;
    std::get<1>(config.data[2][0]) = 10;

    // cell 2 has mutation 2
    std::get<1>(config.data[1][1]) = 10;
    
    // cell 2 has mutation 3
    std::get<1>(config.data[1][2]) = 10;

    config.getLogScores().resizeNumCells(config.getNumSamples());
    config.getTmpLogScores().resizeNumCells(config.getNumSamples());
    config.getLogScores().resizeNumMuts(config.getNumMutations());
    computeLogScoresOP(config);
    config._tmpAttachmentScore.resize(5);

    std::vector<std::vector<unsigned> > mutationsOfNodes = getMutationsOfNode(config);


    BOOST_CHECK_MESSAGE(mutationsOfNodes.size() == 6, "mutationsOfNodes.size() = " << mutationsOfNodes.size() << " != 6");
    BOOST_CHECK_MESSAGE(mutationsOfNodes[0].size() == 0, "mutationsOfNodes[0].size() = " << mutationsOfNodes[0].size() << " != 0");
    BOOST_CHECK_MESSAGE(mutationsOfNodes[1][0] == 0, "mutationsOfNodes[1][0] = " << mutationsOfNodes[1][0] << " != 0");
    BOOST_CHECK_MESSAGE(mutationsOfNodes[2].size() == 0, "mutationsOfNodes[2].size() = " << mutationsOfNodes[2].size() << " != 0");
    BOOST_CHECK_MESSAGE(mutationsOfNodes[3][0] == 1, "mutationsOfNodes[3][0] = " << mutationsOfNodes[3][0] << " != 1");
    BOOST_CHECK_MESSAGE(mutationsOfNodes[3][1] == 2, "mutationsOfNodes[3][1] = " << mutationsOfNodes[3][1] << " != 2");
    BOOST_CHECK_MESSAGE(mutationsOfNodes[4].size() == 0, "mutationsOfNodes[4].size() = " << mutationsOfNodes[4].size() << " != 0");
    BOOST_CHECK_MESSAGE(mutationsOfNodes[5].size() == 0, "mutationsOfNodes[5].size() = " << mutationsOfNodes[5].size() << " != 0");
}

BOOST_AUTO_TEST_CASE( sample_tree_get_num_placements )
{
    {
        Config<SampleTree> config;
        config.setNumSamples(4);

        add_edge(0, 1, config.getTree());
        config.getTree()[0].sample = -1;
        config.getTree()[1].sample = -1;
        add_edge(0, 2, config.getTree());
        config.getTree()[2].sample = -1;
        add_edge(1, 3, config.getTree());
        config.getTree()[3].sample = 0;
        add_edge(1, 4, config.getTree());
        config.getTree()[4].sample = 1;
        add_edge(2, 5, config.getTree());
        config.getTree()[5].sample = 2;
        add_edge(2, 6, config.getTree());
        config.getTree()[6].sample = 3;
        add_edge(7, 0, config.getTree());
        config.getTree()[7].sample = -1;

        unsigned numPlacements = getNumPlacements(config);

        BOOST_CHECK_MESSAGE(numPlacements == 2, "numPlacements = " << numPlacements << " != 2");
    }
    {
        Config<SampleTree> config;
        config.setNumSamples(4);

        add_edge(0, 1, config.getTree());
        config.getTree()[0].sample = -1;
        config.getTree()[1].sample = -1;
        add_edge(1, 2, config.getTree());
        config.getTree()[2].sample = -1;
        add_edge(2, 3, config.getTree());
        config.getTree()[3].sample = 0;
        add_edge(2, 4, config.getTree());
        config.getTree()[4].sample = 1;
        add_edge(1, 5, config.getTree());
        config.getTree()[5].sample = 2;
        add_edge(0, 6, config.getTree());
        config.getTree()[6].sample = 3;
        add_edge(7, 0, config.getTree());
        config.getTree()[7].sample = -1;

        unsigned numPlacements = getNumPlacements(config);

        BOOST_CHECK_MESSAGE(numPlacements == 3, "numPlacements = " << numPlacements << " != 3");
    }
}

BOOST_AUTO_TEST_CASE( sample_tree_compute_score )
{
    {
        {
            Config<SampleTree> config;
            config.setNumSamples(8);

            add_edge(0, 1, config.getTree());
            config.getTree()[0].sample = -1;
            config.getTree()[1].sample = -1;
            add_edge(0, 2, config.getTree());
            config.getTree()[2].sample = -1;
            add_edge(1, 3, config.getTree());
            config.getTree()[3].sample = -1;
            add_edge(1, 4, config.getTree());
            config.getTree()[4].sample = -1;
            add_edge(3, 5, config.getTree());
            config.getTree()[5].sample = -1;
            add_edge(3, 6, config.getTree());
            config.getTree()[6].sample = -1;

            add_edge(5, 7, config.getTree());
            config.getTree()[7].sample = 0;
            add_edge(5, 8, config.getTree());
            config.getTree()[8].sample = 1;

            add_edge(6, 9, config.getTree());
            config.getTree()[9].sample = 2;
            add_edge(6, 10, config.getTree());
            config.getTree()[10].sample = 3;

            add_edge(4, 11, config.getTree());
            config.getTree()[11].sample = 4;
            add_edge(4, 12, config.getTree());
            config.getTree()[12].sample = 5;

            add_edge(2, 13, config.getTree());
            config.getTree()[13].sample = 6;
            add_edge(2, 14, config.getTree());
            config.getTree()[14].sample = 7;
            
            add_edge(15, 0, config.getTree());

            std::get<0>(config.logScores).resizeNumCells(8);
            std::get<0>(config.logScores).resizeNumMuts(1);
            std::get<0>(config.logScores).wtScore(0,0) = std::log(0.01);
            std::get<0>(config.logScores).hetScore(0,0) = std::log(0.09);
            std::get<0>(config.logScores).homScore(0,0) = std::log(0.05);

            std::get<0>(config.logScores).wtScore(1,0) = std::log(0.02);
            std::get<0>(config.logScores).hetScore(1,0) = std::log(0.08);
            std::get<0>(config.logScores).homScore(1,0) = std::log(0.04);

            std::get<0>(config.logScores).wtScore(2,0) = std::log(0.03);
            std::get<0>(config.logScores).hetScore(2,0) = std::log(0.07);
            std::get<0>(config.logScores).homScore(2,0) = std::log(0.03);

            std::get<0>(config.logScores).wtScore(3,0) = std::log(0.04);
            std::get<0>(config.logScores).hetScore(3,0) = std::log(0.05);
            std::get<0>(config.logScores).homScore(3,0) = std::log(0.02);

            std::get<0>(config.logScores).wtScore(4,0) = std::log(0.05);
            std::get<0>(config.logScores).hetScore(4,0) = std::log(0.04);
            std::get<0>(config.logScores).homScore(4,0) = std::log(0.08);

            std::get<0>(config.logScores).wtScore(5,0) = std::log(0.06);
            std::get<0>(config.logScores).hetScore(5,0) = std::log(0.04);
            std::get<0>(config.logScores).homScore(5,0) = std::log(0.07);

            std::get<0>(config.logScores).wtScore(6,0) = std::log(0.07);
            std::get<0>(config.logScores).hetScore(6,0) = std::log(0.03);
            std::get<0>(config.logScores).homScore(6,0) = std::log(0.06);

            std::get<0>(config.logScores).wtScore(7,0) = std::log(0.08);
            std::get<0>(config.logScores).hetScore(7,0) = std::log(0.02);
            std::get<0>(config.logScores).homScore(7,0) = std::log(0.05);

            config.noiseScore = -1000;
            config.data.resize(8);
            config.data[0].resize(1);
            config._tmpAttachmentScore.resize(15);

            scoreTree(config);

            double result = std::log(0.09) + std::log(0.08) + std::log(0.07) + std::log(0.05) + std::log(0.04) + std::log(0.04) - std::log(0.01) - std::log(0.02) - std::log(0.03) - std::log(0.04) - std::log(0.05) - std::log(0.06);
            BOOST_CHECK_MESSAGE(std::abs(config._tmpAttachmentScore.hetScore(1)-result) <= std::numeric_limits<double>::epsilon(), std::exp(config._tmpAttachmentScore.hetScore(1)) << " != " << result);
            result = std::log(0.05) + std::log(0.04) + std::log(0.03) + std::log(0.02) + std::log(0.08) + std::log(0.07) - std::log(0.01) - std::log(0.02) - std::log(0.03) - std::log(0.04) - std::log(0.05) - std::log(0.06);
            BOOST_CHECK_MESSAGE(std::abs(config._tmpAttachmentScore.homScore(1)-result) <= 10 * std::numeric_limits<double>::epsilon(), config._tmpAttachmentScore.homScore(1) << " != " << result);

            result = 1.0/15.0 * 0.8 * exp(config._tmpAttachmentScore.hetScore(1)) + 1.0/7.0 * 0.2 * exp(config._tmpAttachmentScore.homScore(1));
            result = std::log(result);
            config._tmpAttachmentScore[1].computeFinalScore(0.2, 0.3, 7, getNumPlacements(config), false, config.computeMixScore);
            BOOST_CHECK_MESSAGE(std::abs(result - config._tmpAttachmentScore[1].finalScore()) <=  10.0 * std::numeric_limits<double>::epsilon(), result << " != " << config._tmpAttachmentScore[1].finalScore());
        }
    }
}

BOOST_AUTO_TEST_CASE( sample_tree_simplify_tree )
{
    Config<SampleTree> config;
    config.setNumSamples(3);
    createInitialTree(config);

    config.data.resize(3);
    config.data[0].resize(3, std::tuple<unsigned, unsigned>{10, 0});
    config.data[1].resize(3, std::tuple<unsigned, unsigned>{10, 0});
    config.data[2].resize(3, std::tuple<unsigned, unsigned>{10, 0});

    // cell  1, 2 and 3 have mutation 1
    std::get<1>(config.data[0][0]) = 10;
    std::get<1>(config.data[1][0]) = 10;
    std::get<1>(config.data[2][0]) = 10;

    // cell 2 has mutation 2
    std::get<1>(config.data[1][1]) = 10;
    
    // cell 2 has mutation 3
    std::get<1>(config.data[1][2]) = 10;

    config.getLogScores().resizeNumCells(config.getNumSamples());
    config.getTmpLogScores().resizeNumCells(config.getNumSamples());
    config.getLogScores().resizeNumMuts(config.getNumMutations());
    config._tmpAttachmentScore.resize(5);
    computeLogScoresOP(config);


    typedef boost::adjacency_list<boost::vecS, boost::vecS, 
                            boost::bidirectionalS, Vertex<SimpleTree>> TGraph;
    TGraph newTree = simplifyTree(config);
    

    BOOST_CHECK_MESSAGE(num_vertices(newTree) == 4, "num_vertices(newTree) = " << num_vertices(newTree) << " != 4");
    BOOST_CHECK_MESSAGE(newTree[0].mutations[0] == 0, "newTree[0].mutations[0] = " << newTree[0].mutations[0] << " != 0");
    BOOST_CHECK_MESSAGE(newTree[1].mutations[0] == 1, "newTree[1].mutations[0] = " << newTree[1].mutations[0] << " != 1");
    BOOST_CHECK_MESSAGE(newTree[1].mutations[1] == 2, "newTree[1].mutations[1] = " << newTree[1].mutations[1] << " != 2");
}

BOOST_AUTO_TEST_CASE( read_data_skip_indels )
{
    string testNucs = ".,..,";
    BOOST_CHECK_MESSAGE(skipIndels(testNucs, 3) == 3, "skipIndels = " << skipIndels(testNucs, 3) << " != 3");
    testNucs = ".,.-3ACG.,";
    BOOST_CHECK_MESSAGE(skipIndels(testNucs, 3) == 7, "skipIndels = " << skipIndels(testNucs, 3) << " != 7");
    testNucs = ".,.+13ACGGGGGGGGGGG.,";
    BOOST_CHECK_MESSAGE(skipIndels(testNucs, 3) == 18, "skipIndels = " << skipIndels(testNucs, 3) << " != 18");
}

BOOST_AUTO_TEST_CASE( read_data_extract_seq_information )
{
    std::vector<string> testLine= {"chr", "13", "REF", "5", "...,,^A", "IIIII","5", "AC$GTN", "IIIII","5", ".-3AACA.T.", "IIIII"};
    std::vector<std::array<unsigned, 5>> counts;
    counts.resize(3);
    counts[0]= {{0, 0, 0, 0, 0}};
    extractSeqInformation(counts[0], testLine, 0);
    BOOST_CHECK_MESSAGE(counts[0][0] == 0, " = counts[1][0] " << counts[0][0] << " != 0");
    BOOST_CHECK_MESSAGE(counts[0][1] == 0, " = counts[1][0] " << counts[0][1] << " != 0");
    BOOST_CHECK_MESSAGE(counts[0][2] == 0, " = counts[1][0] " << counts[0][2] << " != 0");
    BOOST_CHECK_MESSAGE(counts[0][3] == 0, " = counts[1][0] " << counts[0][3] << " != 0");
    BOOST_CHECK_MESSAGE(counts[0][4] == 5, " = counts[1][0] " << counts[0][4] << " != 5");

    counts[1]= {{0, 0, 0, 0, 0}};
    extractSeqInformation(counts[1], testLine, 1);
    BOOST_CHECK_MESSAGE(counts[1][0] == 1, " = counts[1][0] " << counts[1][0] << " != 1");
    BOOST_CHECK_MESSAGE(counts[1][1] == 1, " = counts[1][0] " << counts[1][1] << " != 1");
    BOOST_CHECK_MESSAGE(counts[1][2] == 1, " = counts[1][0] " << counts[1][2] << " != 1");
    BOOST_CHECK_MESSAGE(counts[1][3] == 1, " = counts[1][0] " << counts[1][3] << " != 1");
    BOOST_CHECK_MESSAGE(counts[1][4] == 4, " = counts[1][0] " << counts[1][4] << " != 4");
    
    counts[2]= {{0, 0, 0, 0, 0}};
    extractSeqInformation(counts[2], testLine, 2);
    BOOST_CHECK_MESSAGE(counts[2][0] == 1, " = counts[1][0] " << counts[2][0] << " != 1");
    BOOST_CHECK_MESSAGE(counts[2][1] == 0, " = counts[1][0] " << counts[2][1] << " != 1");
    BOOST_CHECK_MESSAGE(counts[2][2] == 0, " = counts[1][0] " << counts[2][2] << " != 1");
    BOOST_CHECK_MESSAGE(counts[2][3] == 1, " = counts[1][0] " << counts[2][3] << " != 1");
    BOOST_CHECK_MESSAGE(counts[2][4] == 5, " = counts[1][0] " << counts[2][4] << " != 5");
}

BOOST_AUTO_TEST_CASE( read_data_seq_error_stats )
{
    std::vector<std::array<unsigned, 5>> counts;
    counts.resize(3);
    counts[0]= {{1, 0, 1, 0, 5}};
    counts[1]= {{1, 1, 0, 0, 4}};
    counts[2]= {{0, 0, 0, 3, 6}};

    unsigned errors = 0;
    unsigned cov = 0;
    updateSeqErrorStats(errors, cov, counts);
    
    BOOST_CHECK_MESSAGE(errors == 7, "errors: " << errors << " != 7");
    BOOST_CHECK_MESSAGE(cov == 15, "cov: " << cov << " != 15");
}

BOOST_AUTO_TEST_CASE( read_data_log_n_choose_k )
{

    std::cout.precision(15);

    BOOST_CHECK_MESSAGE(std::round(exp(logNChooseK(10, 0, 0))) == 1, exp(logNChooseK(10, 0, 1)) << " != " << 1);   
    BOOST_CHECK_MESSAGE(std::round(exp(logNChooseK(10, 0, log(2)))) == 1, exp(logNChooseK(10, 0, log(2))) << " != " << 1);   
    BOOST_CHECK_MESSAGE(std::round(exp(logNChooseK(10, 10, log(10)))) == 1, exp(logNChooseK(10, 10, log(10))) << " != " << 1);   
    BOOST_CHECK_MESSAGE(std::round(exp(logNChooseK(10, 6, log(252)))) == 210, exp(logNChooseK(10, 6, log(252))) << " != " << 210);   

    BOOST_CHECK_MESSAGE(std::round(exp(logNChoose2(2))) == 1, exp(logNChoose2(2)) << " != " << 1);   
    BOOST_CHECK_MESSAGE(std::round(exp(logNChoose2(3))) == 3, exp(logNChoose2(3)) << " != " << 3);   
    BOOST_CHECK_MESSAGE(std::round(exp(logNChoose2(10))) == 45, exp(logNChoose2(10)) << " != " << 45);   
}

BOOST_AUTO_TEST_CASE( read_data_pass_cov_filter )
{
    std::vector<string> testLine= {"chr", "13", "REF", "2", "..", "II","2", "AC", "II","5", ".-3AACA.T.", "IIIII"};
    std::vector<std::array<unsigned, 5>> counts;
    counts.resize(3);
    std::vector<unsigned> pos = {0, 1, 2};

    extractSeqInformation(counts, testLine, pos);

    BOOST_CHECK_MESSAGE(passCovFilter(counts[0][4], 3) == false, "passCovFilter(counts[0][4], 3) = " << passCovFilter(counts[0][4], 3));
    BOOST_CHECK_MESSAGE(passCovFilter(counts[1][4], 3) == false, "passCovFilter(counts[1][4], 3) = " << passCovFilter(counts[1][4], 3));
    BOOST_CHECK_MESSAGE(passCovFilter(counts[2][4], 3) == true, "passCovFilter(counts[2][4], 3) = " << passCovFilter(counts[2][4], 3));
}

BOOST_AUTO_TEST_CASE( read_data_pass_supp_filter )
{
    std::vector<string> testLine= {"chr", "13", "REF", "6", "..ACGG", "IIIIII","2", "AC", "II","5", ".-3AACAAT.", "IIIII"};
    std::vector<std::array<unsigned, 5>> counts;
    counts.resize(3);

    std::vector<unsigned> pos = {0, 1, 2};

    extractSeqInformation(counts, testLine, pos);

    BOOST_CHECK_MESSAGE(passSuppFilter(counts[0][0], 2) == false, "passSuppFilter(counts[0][0], 2) = " << passSuppFilter(counts[0][0], 2));
    BOOST_CHECK_MESSAGE(passSuppFilter(counts[0][1], 1) == true, "passSuppFilter(counts[0][1], 1) = " << passSuppFilter(counts[0][1], 1));
    BOOST_CHECK_MESSAGE(passSuppFilter(counts[0][2], 2) == true, "passSuppFilter(counts[0][2], 2) = " << passSuppFilter(counts[0][2], 2));
}

BOOST_AUTO_TEST_CASE( read_data_pass_freq_filter )
{
    std::vector<string> testLine= {"chr", "13", "REF", "6", "..ACGG", "IIIIII","2", "AC", "II","5", ".-3AACAAT.", "IIIII"};
    std::vector<std::array<unsigned, 5>> counts;
    counts.resize(3);

    std::vector<unsigned> pos = {0, 1, 2};

    extractSeqInformation(counts, testLine, pos);

    BOOST_CHECK_MESSAGE(passFreqFilter(counts[0][2], counts[0][4], 0.4) == false, "passFreqFilter(counts[0][2], counts[0][4], 0.4) = " << passFreqFilter(counts[0][2], counts[0][4], 0.4));
    BOOST_CHECK_MESSAGE(passFreqFilter(counts[0][2], counts[0][4], 0.33) == true, "passFreqFilter(counts[0][2], counts[0][4], 0.4) = " << passFreqFilter(counts[0][2], counts[0][4], 0.4));
}

BOOST_AUTO_TEST_CASE( read_data_apply_filter_across_cells )
{
    Config<SampleTree> config;
    config.setNumSamples(3);
    config.minCoverageAcrossCells = 5;
    config.minCoverage = 5;
    config.minFreq = 0.25;
    config.minSupport = 2;
    config.minNumCellsPassFilter = 2;

    std::vector<string> testLine= {"chr", "13", "REF", "4", "..AA", "IIII","2", "AC", "II","5", ".-3AACAAT.", "IIIII"};
    std::vector<std::array<unsigned, 5>> counts;
    counts.resize(3);
    std::vector<unsigned> pos = {0, 1, 2};

    extractSeqInformation(counts, testLine, pos);

    BOOST_CHECK_MESSAGE(applyFilterAcrossCells(counts, config, 0) == false, "applyFilterAcrossCells(counts, config, 0) = " << applyFilterAcrossCells(counts, config, 0));

    testLine= {"chr", "13", "REF", "5", "...AA", "IIIII","2", "AC", "II","5", ".-3AACAAT.", "IIIII"};
    extractSeqInformation(counts, testLine, pos);

    BOOST_CHECK_MESSAGE(applyFilterAcrossCells(counts, config, 0) == true, "applyFilterAcrossCells(counts, config, 0) = " << applyFilterAcrossCells(counts, config, 0));

    testLine= {"chr", "13", "REF", "10", ",,,,,....A", "IIIIIIIIII","2", "AC", "II","5", ".-3AACAAT.", "IIIII"};
    extractSeqInformation(counts, testLine, pos);

    BOOST_CHECK_MESSAGE(applyFilterAcrossCells(counts, config, 0) == false, "applyFilterAcrossCells(counts, config, 0) = " << applyFilterAcrossCells(counts, config, 0));

    config.minCoverageAcrossCells = 3;
    testLine= {"chr", "13", "REF", "5", "..A", "IIIII","2", "AC", "II","5", ".-3AACAAT.", "IIIII"};
    extractSeqInformation(counts, testLine, pos);

    BOOST_CHECK_MESSAGE(applyFilterAcrossCells(counts, config, 0) == false, "applyFilterAcrossCells(counts, config, 0) = " << applyFilterAcrossCells(counts, config, 0));

}
//BOOST_AUTO_TEST_CASE( compute_log_scores )
//{
//    Config<SampleTree> config;
//    config.getData().resize(1);
//    config.getData()[0].push_back(std::pair<unsigned, unsigned>(50, 25));
//    config.getData()[0].push_back(std::pair<unsigned, unsigned>(50, 2));
//    config.getData()[0].push_back(std::pair<unsigned, unsigned>(50, 45));
//
//    config.getLogScores().resize(1);
//    config.getLogScores()[0].resize(3);
//
//    computeLogScores(config);
//    
//    BOOST_CHECK_MESSAGE((int)(std::get<0>(config.getLogScores()[0][0]) * 1000.0) == -10759, "logScores[0][0] = " << std::get<0>(config.getLogScores()[0][0]) << " != -10.75");
//    BOOST_CHECK_MESSAGE((int)(std::get<0>(config.getLogScores()[0][1]) * 1000.0) == -3415, "logScores[0][1] " << std::get<0>(config.getLogScores()[0][1]) << " != -3.415");
//    BOOST_CHECK_MESSAGE((int)(std::get<0>(config.getLogScores()[0][2]) * 1000.0) == -21374, "logScores[0][2] " << std::get<0>(config.getLogScores()[0][2]) << " != -21.37");
//
//    BOOST_CHECK_MESSAGE((int)(std::get<1>(config.getLogScores()[0][0]) * 1000.0) == -3103, "logScores[0][0] " << std::get<1>(config.getLogScores()[0][0]) << " != -3.103");
//    BOOST_CHECK_MESSAGE((int)(std::get<1>(config.getLogScores()[0][1]) * 1000.0) == -8037, "logScores[0][1] " << std::get<1>(config.getLogScores()[0][1]) << " != -8.037");
//    BOOST_CHECK_MESSAGE((int)(std::get<1>(config.getLogScores()[0][2]) * 1000.0) == -6154, "logScores[0][2] " << std::get<1>(config.getLogScores()[0][2]) << " != -6.154");
//    
//    config.getLogScores()[0].clear();
//    config.getLogScores()[0].resize(3);
//    
//    computeWildLogScores(config);
//    BOOST_CHECK_MESSAGE((int)(std::get<0>(config.getLogScores()[0][0]) * 1000.0) == -10759, "logScores[0][0] = " << std::get<0>(config.getLogScores()[0][0]) << " != -10.75");
//    BOOST_CHECK_MESSAGE((int)(std::get<0>(config.getLogScores()[0][1]) * 1000.0) == -3415, "logScores[0][1] " << std::get<0>(config.getLogScores()[0][1]) << " != -3.415");
//    BOOST_CHECK_MESSAGE((int)(std::get<0>(config.getLogScores()[0][2]) * 1000.0) == -21374, "logScores[0][2] " << std::get<0>(config.getLogScores()[0][2]) << " != -21.37");
//
//    computeMutLogScores(config);
//    BOOST_CHECK_MESSAGE((int)(std::get<1>(config.getLogScores()[0][0]) * 1000.0) == -3103, "logScores[0][0] " << std::get<1>(config.getLogScores()[0][0]) << " != -3.103");
//    BOOST_CHECK_MESSAGE((int)(std::get<1>(config.getLogScores()[0][1]) * 1000.0) == -8037, "logScores[0][1] " << std::get<1>(config.getLogScores()[0][1]) << " != -8.037");
//    BOOST_CHECK_MESSAGE((int)(std::get<1>(config.getLogScores()[0][2]) * 1000.0) == -6154, "logScores[0][2] " << std::get<1>(config.getLogScores()[0][2]) << " != -6.154");
//}
//
//BOOST_AUTO_TEST_CASE( get_best_attachment_score )
//{
//    Config<SampleTree> config;
//    config.setNumSamples(3);
//    createInitialTree(config);
//
//
//    config.getData().resize(3);
//    config.getData()[0].push_back(std::pair<unsigned, unsigned>(50, 25));
//    config.getData()[1].push_back(std::pair<unsigned, unsigned>(50, 25));
//    config.getData()[2].push_back(std::pair<unsigned, unsigned>(50, 25));
//
//    config.getData()[0].push_back(std::pair<unsigned, unsigned>(50, 2));
//    config.getData()[1].push_back(std::pair<unsigned, unsigned>(50, 25));
//    config.getData()[2].push_back(std::pair<unsigned, unsigned>(50, 25));
//    
//    config.getData()[0].push_back(std::pair<unsigned, unsigned>(50, 2));
//    config.getData()[1].push_back(std::pair<unsigned, unsigned>(50, 25));
//    config.getData()[2].push_back(std::pair<unsigned, unsigned>(50, 2));
//
//    config.getData()[0].push_back(std::pair<unsigned, unsigned>(50, 2));
//    config.getData()[1].push_back(std::pair<unsigned, unsigned>(50, 2));
//    config.getData()[2].push_back(std::pair<unsigned, unsigned>(50, 25));
//
//    config.getData()[0].push_back(std::pair<unsigned, unsigned>(50, 25));
//    config.getData()[1].push_back(std::pair<unsigned, unsigned>(50, 2));
//    config.getData()[2].push_back(std::pair<unsigned, unsigned>(50, 2));
//
//    config.getTmpAttachmentScore().resize(5);
//    config.getLogScores().resize(3);
//    for (unsigned i = 0; i < config.getLogScores().size(); ++i)
//        config.getLogScores()[i].resize(5);
//
//    computeLogScores(config);
//
//    BOOST_CHECK_MESSAGE((int)(getBestAttachmentScore(config, 0) * 1000) == -9310, "The best attachment score should be " << -9.31041 << " but is " << getBestAttachmentScore(config, 0) << ".");
//    BOOST_CHECK_MESSAGE((int)(getBestAttachmentScore(config, 1) * 1000) == -9621, "The best attachment score should be " << -9.62194 << " but is " << getBestAttachmentScore(config, 1) << ".");
//    BOOST_CHECK_MESSAGE((int)(getBestAttachmentScore(config, 2) * 1000) == -9933, "The best attachment score should be " << -9.93347 << " but is " << getBestAttachmentScore(config, 2) << ".");
//    BOOST_CHECK_MESSAGE((int)(getBestAttachmentScore(config, 3) * 1000) == -9933, "The best attachment score should be " << -9.93347 << " but is " << getBestAttachmentScore(config, 3) << ".");
//    BOOST_CHECK_MESSAGE((int)(getBestAttachmentScore(config, 4) * 1000) == -9933, "The best attachment score should be " << -9.93347 << " but is " << getBestAttachmentScore(config, 4) << ".");
//}
