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
#ifndef NOISE_COUNTS_H
#define NOISE_COUNTS_H

struct GappedNoiseCounts
{
    std::vector<uint64_t> cov;
    std::vector<uint64_t> sup;
    std::vector<uint64_t> covMinusSup;


    void add(std::vector<uint64_t> & vec, unsigned data)
    {
        if (data >= vec.size())
        {
            vec.resize(data + 1, 0);
        }
        ++vec[data];
    }

    void add(unsigned newSup, unsigned newCov)
    {
        this->add(this->cov, newCov);
        this->add(this->sup, newSup);
        this->add(this->covMinusSup, newCov - newSup);
    }

};

std::ostream& operator<<(std::ostream & os, const GappedNoiseCounts & noiseCounts)  
{  
    os << "Coverage: " << noiseCounts.cov.size() << ":";
    for (unsigned i = 0; i < noiseCounts.cov.size(); ++i)
        os << "\t" << i<< ":" << noiseCounts.cov[i];
    os << "\nSupport:  " << noiseCounts.sup.size() << ":";
    for (unsigned i = 0; i < noiseCounts.sup.size(); ++i)
        os << "\t" << i<< ":" << noiseCounts.sup[i];
    os << "\nCov-Sup  : " << noiseCounts.covMinusSup.size() << ":";
    for (unsigned i = 0; i < noiseCounts.covMinusSup.size(); ++i)
        os << "\t" << i<< ":" << noiseCounts.covMinusSup[i];
    return os;  
}  

struct NoiseCounts
{
    std::vector<std::pair<uint32_t, uint64_t>> cov;
    std::vector<std::pair<uint32_t, uint64_t>> sup;
    std::vector<std::pair<uint32_t, uint64_t>> covMinusSup;
    uint64_t numPos;

    NoiseCounts(){}

    NoiseCounts(GappedNoiseCounts const & gappedNoiseCounts)
    {
        this->numPos = 0;

        this->cov.resize(0);
        for (unsigned i = 0; i < gappedNoiseCounts.cov.size(); ++i)
        {
            if (gappedNoiseCounts.cov[i] > 0)
            {
                this->cov.resize(this->cov.size() + 1, std::make_pair(i, gappedNoiseCounts.cov[i]));
                this->numPos += gappedNoiseCounts.cov[i];
            }
        }
        this->sup.resize(0);
        for (unsigned i = 0; i < gappedNoiseCounts.sup.size(); ++i)
        {
            if (gappedNoiseCounts.sup[i] > 0)
            {
                this->sup.resize(this->sup.size() + 1, std::make_pair(i, gappedNoiseCounts.sup[i]));
            }
        }
        this->covMinusSup.resize(0);
        for (unsigned i = 0; i < gappedNoiseCounts.covMinusSup.size(); ++i)
        {
            if (gappedNoiseCounts.covMinusSup[i] > 0)
            {
                this->covMinusSup.resize(this->covMinusSup.size() + 1, std::make_pair(i, gappedNoiseCounts.covMinusSup[i]));
            }
        }
    }
};

std::ostream& operator<<(std::ostream & os, const NoiseCounts & noiseCounts)  
{  
    os << "Coverage: " << noiseCounts.cov.size() << ":";
    for (unsigned i = 0; i < noiseCounts.cov.size(); ++i)
        os << "\t" << noiseCounts.cov[i].first << ":" << noiseCounts.cov[i].second;
    os << "\nSupport:  " << noiseCounts.sup.size() << ":";
    for (unsigned i = 0; i < noiseCounts.sup.size(); ++i)
        os << "\t" << noiseCounts.sup[i].first << ":" << noiseCounts.sup[i].second;
    os << "\nCov-Sup  : " << noiseCounts.covMinusSup.size() << ":";
    for (unsigned i = 0; i < noiseCounts.covMinusSup.size(); ++i)
        os << "\t" << noiseCounts.covMinusSup[i].first << ":" << noiseCounts.covMinusSup[i].second;
    return os;  
}  
#endif
