
#include "AbstractSearchAlgorithm.h"

double AbstractSearchAlgorithm::diff_clock(const clock_t &clock1, const clock_t &clock2) {
    double diff_ticks = double(clock1) - clock2;
    double diff_ms = (diff_ticks) / (double(CLOCKS_PER_SEC) / 1000);
    return diff_ms / 1000;

}

void AbstractSearchAlgorithm::generate_stats(const Node &current_node) {
    std::cout << "Problem : " << getProblemName() << std::endl;
    std::cout << "total nodes explored: " << getExpanded() << std::endl;
    if (getEndStatus()) {
        auto path_start = current_node.getPathTilNow().begin();
        auto path_end = current_node.getPathTilNow().end();
        for (auto item = path_start + 1; item != path_end; ++item) {
            if (item->first == (item - 1)->first - 1 && item->second == (item - 1)->second - 1)
                std::cout << "LU" << '-';
            if (item->first == (item - 1)->first - 1 && item->second == (item - 1)->second)
                std::cout << "U" << '-';
            if (item->first == (item - 1)->first + 1 && item->second == (item - 1)->second - 1)
                std::cout << "LD" << '-';
            if (item->first == (item - 1)->first + 1 && item->second == (item - 1)->second)
                std::cout << "D" << '-';
            if (item->first == (item - 1)->first && item->second == (item - 1)->second + 1)
                std::cout << "R" << '-';
            if (item->first == (item - 1)->first && item->second == (item - 1)->second - 1)
                std::cout << "L" << '-';
            if (item->first == (item - 1)->first + 1 && item->second == (item - 1)->second + 1)
                std::cout << "RD" << '-';
            if (item->first == (item - 1)->first - 1 && item->second == (item - 1)->second + 1)
                std::cout << "RU" << '-';
        }
        std::cout << " solution cost: " << current_node.getActualCost();
        calcEBF(current_node.getDepth());
        calcDN(current_node.getDepth());
    } else {
        std::cout << "Failed";
        calcEBF(getMax());
        calcDN(getMax());
    }
    std::cout << std::endl;
    std::cout << "d/N : " << getDN() << std::endl;
    std::cout << "time in seconds: " << diff_clock(_current_time, _start_time) << std::endl;
    std::cout << "EBF : " << getEBF() << std::endl;
    std::cout << "Min cutoff : " << getMin() << std::endl;
    std::cout << "Max cutoff : " << getMax() << std::endl;
    std::cout << "Avg cutoff : " << getAvg() << std::endl;
}

double AbstractSearchAlgorithm::getAvg() const {
    if (_no_of_cutoffs == 0)
        return 0;
    return double(_sum_of_cutoffs_depths) / _no_of_cutoffs;
}

void AbstractSearchAlgorithm::addCutoffToSum(int cut_off_depth) {
    _sum_of_cutoffs_depths += cut_off_depth;
    ++_no_of_cutoffs;
}

void AbstractSearchAlgorithm::calcEBF(int depth) {
    if (!depth)
        return;
    setEBF(pow(getExpanded(), pow(depth, -1)));
}

void AbstractSearchAlgorithm::update_cutoffs(int cutoff_depth) {
    if (_min > cutoff_depth)
        _min = cutoff_depth;
    if (_max < cutoff_depth)
        _max = cutoff_depth;
    addCutoffToSum(cutoff_depth);
}

void AbstractSearchAlgorithm::calcDN(int depth) {
    if (!getExpanded())
        return;
    setDN(double(depth) / getExpanded());
}