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

#ifndef ATTACHMENTSCORES_H
#define ATTACHMENTSCORES_H

#include "sciphi_config.h"

// Forward declaration
double addLogProbWeight(double x, double y, double nu);
double addLogProb(double x, double y);

/*
 * This class stores values associated with the attachment of mutatations
 * to nodes in the phylogenetic tree.
 * Here N is the number of different attachment scores (currently only 5 is
 * used).
 */
template <unsigned N>
struct AttachmentScore{

    typedef std::array<double, N> TAttachmentScore;

    TAttachmentScore attachmentScore;

    // init everything to -INFINITY because we are doing the computation in 
    // log space and exp(-INFINITY) = 0
    AttachmentScore() :
        attachmentScore({{-INFINITY,-INFINITY,-INFINITY,-INFINITY,-INFINITY}}) {};
    
    std::array<double, N> makeArray()
    {
        std::array<double, N> result;
        result[0]=this->hetScore();
        result[1]=this->homScore();
        result[N-1]=this->finalScore();

        if (N == 3)
        {
            return result;
        }

        result[2]=this->mixWildScore();
        result[3]=this->mixHomScore();
        return result;
    }
    
    // transform attachment scores into log
    void log()
    {
        for (unsigned i = 0; i < N; ++i)
        {
            this->attachmentScore[i] = std::log(this->attachmentScore[i]);
        }
    }

    // exponentiate attachment scores
    void exp()
    {
        for (unsigned i = 0; i < N; ++i)
        {
            this->attachmentScore[i] = std::exp(this->attachmentScore[i]);
        }
    }

    double & hetScore()
    {
        return this->attachmentScore[0];
    }
    double const & hetScore() const
    {
        return this->attachmentScore[0];
    }

    double & homScore()
    {
        return this->attachmentScore[1];
    }
    double const & homScore() const
    {
        return this->attachmentScore[1];
    }
    
    double & mixWildScore()
    {
        return this->attachmentScore[2];
    }
    double const & mixWildScore() const
    {
        return this->attachmentScore[2];
    }
    
    double & mixHomScore()
    {
        return this->attachmentScore[3];
    }
    double const & mixHomScore() const
    {
        return this->attachmentScore[3];
    }

    double & finalScore()
    {
        return this->attachmentScore[N-1];
    }
    double const & finalScore() const
    {
        return this->attachmentScore[N-1];
    }

    AttachmentScore & operator+=(AttachmentScore const & rightSide)
    {
        this->hetScore() += rightSide.hetScore();
        this->homScore() += rightSide.homScore();
        if (N > 3)
        {
            // this check is necessary as there is not always a mix score available to be added
            if(!isnan(rightSide.mixWildScore())) 
            {
                this->mixWildScore() += rightSide.mixWildScore();
            }
            if(!isnan(rightSide.mixHomScore()))
            {
                this->mixHomScore() += rightSide.mixHomScore();
            }
        }
        this->finalScore() += rightSide.finalScore();

        return *this;
    }
    AttachmentScore & operator+=(TAttachmentScore & rightSide)
    {
        this->hetScore() += rightSide.hetScore();
        this->homScore() += rightSide.homScore();
        if (N > 3)
        {
            // this check is necessary as there is not always a mix score available to be added
            if(!isnan(rightSide.mixWildScore()))
            {
                this->mixWildScore() += rightSide.mixWildScore();
            }
            if(!isnan(rightSide.mixHomScore()))
            {
                this->mixHomScore() += rightSide.mixHomScore();
            }
        }
        this->finalScore() += rightSide.finalScore();

        return *this;
    }

    AttachmentScore & operator-=(AttachmentScore & rightSide)
    {
        this->hetScore() -= rightSide.hetScore();
        this->homScore() -= rightSide.homScore();
        if (N > 3)
        {
            // this check is necessary as there is not always a mix score available to be added
            if(!isnan(rightSide.mixWildScore()))
            {
                this->mixWildScore() -= rightSide.mixWildScore();
            }
            if(!isnan(rightSide.mixHomScore()))
            {
                this->mixHomScore() -= rightSide.mixHomScore();
            }
        }
        this->finalScore() -= rightSide.finalScore();

        return *this;
    }
    
    bool operator!=(AttachmentScore & rightSide)
    {
        if (this->hetScore() != rightSide.hetScore())
            return false;
        if (this->homScore() != rightSide.homScore())
            return false;
        if (this->mixWildScore() != rightSide.mixWildScore())
            return false;
        if (this->mixHomScore() != rightSide.mixHomScore())
            return false;
        if (this->finalScore() != rightSide.finalScore())
            return false;
        return true;
    }

    // this function takes two scores in log space, exponentiates them and 
    // returns the log of the sum
    void addInRealSpace(AttachmentScore const & rightSide)
    {
        this->hetScore() = addLogProb(this->hetScore(), rightSide.hetScore());
        this->homScore() = addLogProb(this->homScore(), rightSide.homScore());
        if (N > 3)
        {
            if(!isnan(rightSide.mixWildScore()))
            {
                this->mixWildScore() = addLogProb(this->mixWildScore(), rightSide.mixWildScore());
            }
            if(!isnan(rightSide.mixHomScore()))
            {
                this->mixHomScore() = addLogProb(this->mixHomScore(), rightSide.mixHomScore());
            }
        }
        this->finalScore() = addLogProb(this->finalScore(), rightSide.finalScore());
    }

    AttachmentScore & operator/=(double rightSide)
    {
        for (unsigned i = 0; i < N; ++i)
        {
            this->attachmentScore[i] /= rightSide;
        }
        return *this;
    }
   
    // this functions computes the final score of an attachment point
    void computeFinalScore(
            double nu, 
            double lambda,
            unsigned numNodes,
            unsigned numMutPlacements,
            bool isLeaf,
            bool useMixture)
    {
        // if the node is a leaf just return the weighted hetero score
        if (isLeaf)
        {
            this->finalScore() = this->hetScore() - std::log(numNodes * 2 + 1);
            return;
        }

        // this is currently alwas true as the mixture is still experimental
        // combine the hetero and the homozygous score
        if (!useMixture || 
                isnan(this->mixWildScore()))
        {
            this->finalScore() = addLogProbWeight(this->hetScore() - std::log(numNodes * 2 + 1), this->homScore() - std::log(numNodes), nu);
            return;
        }

        double mixHomo = addLogProbWeight(this->mixWildScore(), this->mixHomScore(), 0.5) - std::log(numMutPlacements);
        double mix = addLogProbWeight(this->hetScore() - std::log(numNodes * 2 + 1), mixHomo, lambda);
        double score = addLogProbWeight(mix, this->homScore()  - std::log(numNodes), nu);

        this->finalScore() = score;
    }
   
    AttachmentScore cellProbReal(AttachmentScore const & sumReal)
    {
        AttachmentScore result = {};
        result.hetScore()=this->hetScore()/sumReal.hetScore();
        result.homScore()=this->homScore()/sumReal.homScore();
        result.mixWildScore()=this->mixWildScore()/sumReal.mixWildScore();
        result.mixHomScore()=this->mixHomScore()/sumReal.mixHomScore();
        result.finalScore()=this->finalScore()/sumReal.finalScore();

        return result;
    }

    void setMinusInfinity()
    {
        this->hetScore() = -INFINITY;
        this->homScore()= -INFINITY;
        this->mixWildScore() = -INFINITY;
        this->mixHomScore() = -INFINITY;
        this->finalScore() = -INFINITY;
    }
};

template <unsigned N>
std::ostream& operator<<(std::ostream& os, AttachmentScore<N> const & obj)
{
    os << obj.hetScore() << "|" << obj.homScore(); 
    if (N == 3)
    {
        return os << "|" << obj.finalScore();
    }
    os << "|" << obj.mixWildScore() << "|" << obj.mixHomScore() << "|" << obj.finalScore();
    return os;
}

std::ostream& operator<<(std::ostream& os, AttachmentScore<3> const & obj)
{
    os << obj.hetScore() << "|" << obj.homScore() << "|" << obj.finalScore();
    return os;
}

std::ostream& operator<<(std::ostream& os, AttachmentScore<5> const & obj)
{
    os << obj.hetScore() << "|" << obj.homScore() << "|" << obj.mixWildScore() << "|" << obj.mixHomScore() << "|" << obj.finalScore();
    return os;
}

template <unsigned N>
AttachmentScore<N> operator+(AttachmentScore<N> & leftSide, AttachmentScore<N> & rightSide)
{
    AttachmentScore<N> result;
    for (unsigned i = 0; i < N; ++i)
    {
        result.attachmentScore[i] =leftSide.attachmentScore[i] + rightSide.attachmentScore[i];
    }
    return result;
}

/*
 * This class combines many AttachmentScores into one object
 */
template <unsigned N>
struct AttachmentScores
{
    typedef AttachmentScore<N> TAttachmentScore;
    typedef std::vector<AttachmentScore<N>> TAttachmentScores;

    TAttachmentScores attachmentScores;

    AttachmentScore<N> & operator[](unsigned pos)
    {
        return this->attachmentScores[pos];
    }
    
    AttachmentScore<N> const & operator[](unsigned pos) const
    {
        return this->attachmentScores[pos];
    }

    double & hetScore(unsigned attachPoint)
    {
        return this->attachmentScores[attachPoint].hetScore();
    }
    double const & hetScore(unsigned attachPoint) const
    {
        return this->attachmentScores[attachPoint].hetScore();
    }

    double & homScore(unsigned attachPoint)
    {
        return this->attachmentScores[attachPoint].homScore();
    }
    double const & homScore(unsigned attachPoint) const
    {
        return this->attachmentScores[attachPoint].homScore();
    }
    
    double & mixWildScore(unsigned attachPoint)
    {
        return this->attachmentScores[attachPoint].mixWildScore();
    }
    double const & mixWildScore(unsigned attachPoint) const
    {
        return this->attachmentScores[attachPoint].mixWildScore();
    }
    
    double & mixHomScore(unsigned attachPoint)
    {
        return this->attachmentScores[attachPoint].mixHomScore();
    }
    double const & mixHomScore(unsigned attachPoint) const
    {
        return this->attachmentScores[attachPoint].mixHomScore();
    }
    
    unsigned size()
    {
        return this->attachmentScores.size();
    }
    
    void resize(unsigned newSize)
    {
        this->attachmentScores.resize(newSize);
    }

    AttachmentScore<N> max()
    {
        AttachmentScore<N> result = this->attachmentScores[0];

        for (unsigned i = 1; i < this->size(); ++i)
        {
            if (result.hetScore() < this->attachmentScores[i].hetScore())
            {
                result.hetScore()=this->attachmentScores[i].hetScore();
            }
            if (result.homScore() < this->attachmentScores[i].homScore())
            {
                result.homScore() = this->attachmentScores[i].homScore();
            }
            if (result.mixWildScore() < this->attachmentScores[i].mixWildScore())
            {
                result.mixWildScore() = this->attachmentScores[i].mixWildScore();
            }
            if (result.mixHomScore() < this->attachmentScores[i].mixHomScore())
            {
                result.mixHomScore() = this->attachmentScores[i].mixHomScore();
            }
        }
        return result;
    }

    /*
    // this function sum the individual attachment scores provided a max value (in log space)
    AttachmentScore<N> sumInReal(AttachmentScore<N> const & max)
    {
        AttachmentScore<N> result;
        for (unsigned i = 1; i < this->size(); ++i)
        {
            result.hetScore() += std::exp(this->attachmentScores[i].hetScore() - max.hetScore());
            result.homScore() += std::exp(this->attachmentScores[i].homScore() - max.homScore());
            result.mixWildScore() += std::exp(this->attachmentScores[i].mixWildScore() - max.mixWildScore());
            result.mixHomScore() += std::exp(this->attachmentScores[i].mixHomScore() - max.mixHomScore());

        }
        result.hetScore() = std::log(result.hetScore()) + max.hetScore();
        result.homScore() = std::log(result.homScore()) + max.homScore();
        result.mixWildScore() = std::log(result.mixWildScore()) + max.mixWildScore();
        result.mixHomScore() = std::log(result.mixHomScore()) + max.mixHomScore();
        return result;
    }

    AttachmentScore<N> sumInReal()
    {
        AttachmentScore<N> max = this->max();
        return sumInReal(max);
    }

    template <typename TTree>
    AttachmentScore<N> cellSumReal(TTree const & tree, unsigned currentNode, AttachmentScore<N> const & max)
    {
        AttachmentScore<N> result = {};
        while(currentNode != num_vertices(tree) - 1)
        {
            result.hetScore() += std::exp(this->attachmentScores[currentNode].hetScore() - max.hetScore());
            result.homScore() += std::exp(this->attachmentScores[currentNode].homScore() - max.homScore());
            result.mixWildScore() += std::exp(this->attachmentScores[currentNode].mixWildScore() - max.mixWildScore());
            result.mixHomScore() += std::exp(this->attachmentScores[currentNode].mixHomScore() - max.mixHomScore());
            currentNode = source(*in_edges(currentNode, tree).first, tree);
        }
        return result;
    }

    template <typename TTree>
    AttachmentScore<N> cellSumReal(TTree const & tree, unsigned currentNode)
    {
        AttachmentScore<N> max = this->max();
        return cellSum(tree, currentNode, max);
    }

    template <typename TTree>
    AttachmentScore<N> cellProbReal(TTree const & tree, unsigned currentNode)
    {
        AttachmentScore<N> max = this->max();
        AttachmentScore<N> sumReal = this->sumReal(max);
        AttachmentScore<N> cellSumReal = this->cellSumReal(tree, currentNode, max);
        return cellSumReal.cellProbReal(sumReal);
    }
    */
    void computeLogHetScoreLeaf(unsigned attachPoint, double wtScore, double hetScore)
    {
        this->hetScore(attachPoint) = hetScore - wtScore;    
    }
    
    void computeLogHetScoreInnerNode(unsigned attachPoint, double hetScoreLeft, double hetScoreRigth)
    {
        this->hetScore(attachPoint) = hetScoreLeft + hetScoreRigth;
    }

    void computeLogHomScoreLeaf(unsigned attachPoint, double wtScore, double homScore)
    {
        this->homScore(attachPoint) = homScore - wtScore;    
    }
    
    void computeLogHomScoreInnerNode(unsigned attachPoint, double homScoreLeft, double homScoreRigth)
    {
        this->homScore(attachPoint) = homScoreLeft + homScoreRigth;
    }
    
    double computeLogMixScoreLeaf()
    {
        return std::nan("");
    }

    void computeLogMixWildScoreLeaf(unsigned attachPoint)
    {
        this->mixWildScore(attachPoint) = computeLogMixScoreLeaf();
    }
    
    void computeLogMixHomScoreLeaf(unsigned attachPoint)
    {
        this->mixHomScore(attachPoint) = computeLogMixScoreLeaf();
    }

    // This function is experimental
    // It is used to compute the mixture homozygous alternative score.
    // Here the chromosom without a mutation is lost after the mutatiin
    // manifested itself.
    void computeLogMixHomScoreInnerNode(
            unsigned attachPoint, 
            bool innerNodeLeft, 
            double homScoreLeft, 
            double mixHomScoreLeft,
            double hetScoreLeft,
            bool innerNodeRight, 
            double homScoreRight, 
            double mixHomScoreRight,
            double hetScoreRight)
    {
        // if both children are leaves no mix score can be computed
        if (!innerNodeLeft && !innerNodeRight)
        {
            this->mixHomScore(attachPoint) = std::nan("");
            return;
        }

        // if the left child is an inner node a mix score can be computed
        double result = std::nan("");
        if (innerNodeLeft && !innerNodeRight)
        {
            // compute the mixture score from the homozygous alternative
            // and hetero score of the children
            result = homScoreLeft + hetScoreRight;
            
            // if the is a mixture score in the inner node add the mixture score
            if (!isnan(mixHomScoreLeft))
            {
                result = addLogProb(result, mixHomScoreLeft + hetScoreRight);
            }
            this->mixHomScore(attachPoint) = result;
            return;
        }
        
        // if the right child is an inner node a mix score can be computed
        if (!innerNodeLeft && innerNodeRight)
        {
            // compute the mixture score from the homozygous alternative
            // and hetero score of the children
            result = homScoreRight + hetScoreLeft;
            
            // if the is a mixture score in the inner node add the mixture score
            if (!isnan(mixHomScoreRight))
            {
                result = addLogProb(result, mixHomScoreRight + hetScoreLeft);
            }
            this->mixHomScore(attachPoint) = result;
            return;
        }

        // if both children are inner nodes proceed similar to above
        result = homScoreLeft + hetScoreRight;
        result = addLogProb(result, homScoreRight + hetScoreLeft);
        if (!isnan(mixHomScoreLeft))
        {
            result = addLogProb(result, mixHomScoreLeft + hetScoreRight);
        }
        if (!isnan(mixHomScoreRight))
        {
            result = addLogProb(result, mixHomScoreRight + hetScoreLeft);
        }
        this->mixHomScore(attachPoint) = result;
    }


    // This function is experimental.
    // It is used to compute the mixture wildtype alternative score.
    // Here the chromosome with the mutation is lost
    void computeLogMixWildScoreInnerNode(
            unsigned attachPoint,
            bool innerNodeLeft, 
            double mixWildScoreLeft,
            double hetScoreLeft,
            bool innerNodeRight, 
            double mixWildScoreRight,
            double hetScoreRight)
    {
        // if both children are leaves no mix score can be computed
        if (!innerNodeLeft && !innerNodeRight)
        {
            this->mixWildScore(attachPoint) = std::nan("");
            return;
        }

        double result = std::nan("");
        if (innerNodeLeft && !innerNodeRight)
        {
            result = hetScoreRight;
            if (!isnan(mixWildScoreLeft))
            {
                result = addLogProb(result, mixWildScoreLeft + hetScoreRight);
            }
            this->mixWildScore(attachPoint) = result;
            return;
        }
        
        if (!innerNodeLeft && innerNodeRight)
        {
            result = hetScoreLeft;
            if (!isnan(mixWildScoreRight))
            {
                result = addLogProb(result, mixWildScoreRight + hetScoreLeft);
            }
            this->mixWildScore(attachPoint) = result;
            return;
        }

        result = addLogProb(hetScoreLeft, hetScoreRight);
        if (!isnan(mixWildScoreLeft))
        {
            result = addLogProb(result, mixWildScoreLeft + hetScoreRight);
        }
        if (!isnan(mixWildScoreRight))
        {
            result = addLogProb(result, mixWildScoreRight + hetScoreLeft);
        }
        this->mixWildScore(attachPoint) = result;
    }


    // Experimental function
    // This function computes the final score of an inner node using the 
    // appropriate weights.
    void computeLogFinalScoreInnerNode(
            unsigned attachPoint, 
            bool innerNodeLeft, 
            bool innerNodeRight, 
            double nu, 
            double lambda,
            unsigned numNodes,
            unsigned numMutPlacements)
    {
        if (!innerNodeLeft && !innerNodeRight)
        {
            this->finalScore(attachPoint) = std::log(1.0 / static_cast<double>(numNodes)) +
                    addLogProbWeight(// hetero or all homo
                            this->hetScore(attachPoint), 
                                    addLogProbWeight(
                                            this->wildScore(attachPoint),
                                            this->homScore(attachPoint),
                                            0.5),
                            nu);
            return;
        }
        else {
            this->finalScore(attachPoint) = 
                    addLogProbWeight(// hetero or homo
                            std::log(1.0 / static_cast<double>(numNodes)) + this->hetScore(attachPoint), 
                            addLogProbWeight(// all homozygous
                                    (1.0 / static_cast<double>(numMutPlacements)) + addLogProbWeight(
                                            this->wildScore(attachPoint),
                                            this->homScore(attachPoint),
                                            0.5),
                                    addLogProbWeight(
                                            this->mixWildScore(attachPoint),
                                            this->mixHomScore(attachPoint),
                                            0.5),
                                    lambda), 
                            nu);

        }
    }
};

#endif
