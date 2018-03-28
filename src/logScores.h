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
#ifndef LOGSCORES_H
#define LOGSCORES_H

struct LogScores
{
    typedef std::vector<std::vector<std::array<double, 4> > > TLogScores;

    TLogScores logScores;

    double & wtScore(unsigned attachPoint, unsigned geneID)
    {
        return this->logScores[attachPoint][geneID][0];
    }
    double const & wtScore(unsigned attachPoint, unsigned geneID) const
    {
        return this->logScores[attachPoint][geneID][0];
    }
    
    double & hetScore(unsigned attachPoint, unsigned geneID)
    {
        return this->logScores[attachPoint][geneID][1];
    }
    double const & hetScore(unsigned attachPoint, unsigned geneID) const
    {
        return this->logScores[attachPoint][geneID][1];
    }

    double & homScore(unsigned attachPoint, unsigned geneID)
    {
        return this->logScores[attachPoint][geneID][2];
    }
    double const & homScore(unsigned attachPoint, unsigned geneID) const
    {
        return this->logScores[attachPoint][geneID][2];
    }
    
    unsigned numCells()
    {
        return this->logScores.size();
    }

    void resizeNumCells(unsigned newSize)
    {
        return this->logScores.resize(newSize);
    }

    unsigned numMuts()
    {
        return this->logScores[0].size();
    }
    
    void resizeNumMuts(unsigned newSize)
    {
        for (unsigned i = 0; i < this->numCells(); ++i)
        {
            this->logScores[i].resize(newSize, {{-INFINITY,-INFINITY,-INFINITY,-INFINITY}});
        }
        return;
    }

};

#endif
