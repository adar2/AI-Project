#include <utility>

#ifndef AI_PROJECT_HEURISTICFUNCTIONS_H
#define AI_PROJECT_HEURISTICFUNCTIONS_H

double chebyshev_distance(const std::pair<int, int> &p1, const std::pair<int, int> &p2, double min_val);

double normalized_euclidean_distance(const std::pair<int, int> &p1, const std::pair<int, int> &p2, double min_val);

#endif //AI_PROJECT_HEURISTICFUNCTIONS_H
