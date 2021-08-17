

#include "BiDirectionalAStar.h"
#include "PairHashing.h"
#include<unordered_map>
#include <set>

using std::multiset;
using std::unordered_map;

int BiDirectionalAStar::run_algorithm(double **array, int dimension, int *source, int *goal, float time_limit) {
    //open list front search
    multiset<shared_ptr<Node>, lessCompNodePointers> frontier_front = multiset<shared_ptr<Node>, lessCompNodePointers>();
    //open list backward search
    multiset<shared_ptr<Node>, lessCompNodePointers> frontier_back = multiset<shared_ptr<Node>, lessCompNodePointers>();
    //closed list front search
    unordered_map<pair<int, int>, shared_ptr<Node>, pair_hash> visited_front = unordered_map<pair<int, int>, shared_ptr<Node>, pair_hash>();
    //closed list backward search
    unordered_map<pair<int, int>, shared_ptr<Node>, pair_hash> visited_back = unordered_map<pair<int, int>, shared_ptr<Node>, pair_hash>();
    // two nodes to hold the solution, one for each direction.
    shared_ptr<Node> current_node, sol1, sol2;
    shared_ptr<vector<shared_ptr<Node>>> successors;
    int expand_counter;
    double sum, min = std::numeric_limits<double>::max();
    // h_cost is the heuristic cost of the node to goal node
    double h_cost = _heuristic_function(pair<int, int>(source[0], source[1]), pair<int, int>(goal[0], goal[1]),
                                        getMinOfCostMatrix());
    // source node to start a search from it to goal node, and from goal node to source node.
    shared_ptr<Node> sourceNode(new Node(h_cost, 0, source[0], source[1], 0));
    // inserts source node coordinates to its path
    sourceNode->insertElementToPath(pair<int, int>(sourceNode->getRow(), sourceNode->getCol()));
    // the heuristic function is symmetric so source to goal is the same as goal to source
    // init goal node with its cost as we are actually searching for the goal node and we need to consider its weight.
    shared_ptr<Node> goalNode(new Node(h_cost, array[goal[0]][goal[1]], goal[0], goal[1], 0));
    // inserts goal node coordinates to its path
    goalNode->insertElementToPath(pair<int, int>(goalNode->getRow(), goalNode->getCol()));
    // init front and back priority queues with source and goal nodes.
    frontier_front.insert(sourceNode);
    frontier_back.insert(goalNode);
    while (!frontier_front.empty() && !frontier_back.empty()) {
        setCurrentTime(clock());
        //check whether node v is at the top of both frontier_front and frontier_back
        if (*(*frontier_front.begin()) == *(*frontier_back.begin())) {
            update_cutoffs((*frontier_front.begin())->getDepth());
            setEndStatus(true);
            break;
        }
        // time out
        if (diff_clock(getCurrentTime(), getStartTime()) >= time_limit) {
            generate_stats(*current_node, getAvgHeuristicValue());
            return 1;
        }
        if (!frontier_front.empty()) {
            // get the node with the smallest h_cost, tie breaking is done by minimal h(n)
            current_node = *frontier_front.begin();
            frontier_front.erase(frontier_front.begin());
            if (visited_front.find(pair<int, int>(current_node->getRow(), current_node->getCol())) !=
                visited_front.end())
                continue;
            ++getInstance()._expanded;
            // mark that node as visited in the front search
            visited_front[pair<int, int>(current_node->getRow(), current_node->getCol())] = current_node;
            // expand counter is used to detect cutoffs and record the for stats.
            expand_counter = 0;
            // loop over node actions
            successors = current_node->successors(array, dimension);
            for (const auto &successor : *successors) {
                expand_counter++;
                h_cost = _heuristic_function(pair<int, int>(successor->getRow(), successor->getCol()),
                                             pair<int, int>(goalNode->getRow(), goalNode->getCol()),
                                             getMinOfCostMatrix());
                successor->setEvaluationCost(successor->getActualCost() + h_cost);
                sumNodeHeuristic(h_cost);
                auto visited_iterator = visited_front.find(pair<int, int>(successor->getRow(), successor->getCol()));
                // if the node is not visited and not in the open list, add it to open list.
                if (visited_iterator == visited_front.end())
                    frontier_front.insert(successor);
                    // if the node has been visited with higher h_cost, remove it from visited and insert it to open list with better h_cost
                else if (visited_iterator != visited_front.end() &&
                         visited_iterator->second->getEvaluationCost() > successor->getEvaluationCost()) {
                    visited_front.erase(visited_iterator);
                    frontier_front.insert(successor);
                    // if we didnt get into one of the above cond, we didnt expended this node so decrement expend_counter.
                } else expand_counter--;
            }
            if (!expand_counter) {
                // cut off occurrence
                update_cutoffs(current_node->getDepth());
            }
        }
        //check again as it may happen between the forward step to the backward step.
        //check whether node v is at the top of both frontier_front and frontier_back
        if (!frontier_front.empty() && !frontier_back.empty() &&
            *(*frontier_front.begin()) == *(*frontier_back.begin())) {
            update_cutoffs((*frontier_front.begin())->getDepth());
            setEndStatus(true);
            break;
        }
        // time out
        if (diff_clock(getCurrentTime(), getStartTime()) >= time_limit) {
            generate_stats(*current_node, getAvgHeuristicValue());
            return 1;
        }
        if (!frontier_back.empty()) {
            // get the node with the smallest h_cost, tie breaking is done by minimal h(n)
            current_node = *frontier_back.begin();
            frontier_back.erase(frontier_back.begin());
            if (visited_back.find(pair<int, int>(current_node->getRow(), current_node->getCol())) != visited_back.end())
                continue;
            ++getInstance()._expanded;
            // mark that node as visited in the front search
            visited_back[pair<int, int>(current_node->getRow(), current_node->getCol())] = current_node;
            // expand counter is used to detect cutoffs and record the for stats.
            expand_counter = 0;
            // loop over node actions
            successors = current_node->successors(array, dimension);
            for (const auto &successor: *successors) {
                expand_counter++;
                h_cost = _heuristic_function(pair<int, int>(successor->getRow(), successor->getCol()),
                                             pair<int, int>(sourceNode->getRow(), sourceNode->getCol()),
                                             getMinOfCostMatrix());
                successor->setEvaluationCost(successor->getActualCost() + h_cost);
                sumNodeHeuristic(h_cost);
                // use std::find_if with lambada function as node pointers comparator.
                auto visited_iterator = visited_back.find(pair<int, int>(successor->getRow(), successor->getCol()));
                // if the node is not visited and not in the open list, add it to open list.
                if (visited_iterator == visited_back.end())
                    frontier_back.insert(successor);
                    // if the node has been visited with higher h_cost, remove it from visited and insert it to open list with better h_cost
                else if (visited_iterator != visited_back.end() &&
                         visited_iterator->second->getEvaluationCost() > successor->getEvaluationCost()) {
                    visited_back.erase(visited_iterator);
                    frontier_back.insert(successor);

                }
                    // if we didnt get into one of the above cond, we didnt expended this node so decrement expend_counter.
                else expand_counter--;
            }
            if (!expand_counter) {
                // cut off occurrence
                update_cutoffs(current_node->getDepth());
            }
        }
    }
    // post search phase find the optimal cost, loop through open and closed lists and look for the path with minimum cost which is the optimal path
    if (getEndStatus()) {
        // multiset that holds the nodes in the front open list that are exists also in the backward open list.
        multiset<shared_ptr<Node>, lessCompNodePointers> intersection_front = multiset<shared_ptr<Node>, lessCompNodePointers>();
        // multiset that holds the nodes in the backward open list that are exists also in the front open list.
        multiset<shared_ptr<Node>, lessCompNodePointers> intersection_back = multiset<shared_ptr<Node>, lessCompNodePointers>();
        auto cmp = [](const shared_ptr<Node> &p1, const shared_ptr<Node> &p2) { return *p1 < *p2 && *p1 == *p2; };
        // insert the intersecting nodes to the multiset respectively
        std::set_intersection(frontier_front.begin(), frontier_front.end(), frontier_back.begin(), frontier_back.end(),
                              std::inserter(intersection_front, intersection_front.begin()), cmp);
        std::set_intersection(frontier_back.begin(), frontier_back.end(), frontier_front.begin(), frontier_front.end(),
                              std::inserter(intersection_back, intersection_back.begin()), cmp);
        const auto intersection_back_begin = intersection_back.begin();
        const auto intersection_back_end = intersection_back.end();
        for (const auto &element : intersection_front) {
            setCurrentTime(clock());
            if (diff_clock(getCurrentTime(), getStartTime()) >= time_limit) {
                setEndStatus(false);
                generate_stats(*current_node, getAvgHeuristicValue());
                return 1;
            }
            auto found = std::find_if(intersection_back_begin, intersection_back_end,
                                      [&element](const shared_ptr<Node> &node) { return *node == *element; });
            if (found != intersection_back_end) {
                sum = element->getActualCost() + (*found)->getActualCost();
                if (sum < min) {
                    min = sum;
                    sol1 = element;
                    sol2 = (*found);
                }
            }

        }
        const auto visited_back_end = visited_back.end();
        for (const auto &element: visited_front) {
            setCurrentTime(clock());
            if (diff_clock(getCurrentTime(), getStartTime()) >= time_limit) {
                setEndStatus(false);
                generate_stats(*current_node, getAvgHeuristicValue());
                return 1;
            }
            auto found = visited_back.find(pair<int, int>(element.second->getRow(), element.second->getCol()));
            if (found != visited_back_end) {
                sum = element.second->getActualCost() + found->second->getActualCost();
                if (sum < min) {
                    min = sum;
                    sol1 = element.second;
                    sol2 = (*found).second;
                }
            }

        }
        for (auto node = sol2->getPathTilNow().rbegin() + 1; node != sol2->getPathTilNow().rend(); ++node) {
            sol1->insertElementToPath(*node);
            sol1->setDepth(sol1->getDepth() + 1);
        }
        sol1->setActualCost(sol1->getActualCost() + sol2->getActualCost() - array[sol1->getRow()][sol1->getCol()]);
        generate_stats(*sol1, getAvgHeuristicValue());
        return 0;
    }
    generate_stats(*sourceNode, getAvgHeuristicValue());

    return 1;
}

BiDirectionalAStar &BiDirectionalAStar::getInstance() {
    static BiDirectionalAStar instance;
    return instance;
}



