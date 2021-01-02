//
// Created by r00t on 12/13/20.
//

#include "IterativeDeepeningSearch.h"


using std::unordered_map;
using std::stack;
using std::cout;
using std::endl;

IterativeDeepeningSearch &IterativeDeepeningSearch::getInstance() {
    static IterativeDeepeningSearch instance;
    return instance;
}


pair<shared_ptr<Node>, bool>
IterativeDeepeningSearch::DLS(int **array, int dimension,const shared_ptr<Node>& root,const shared_ptr<Node>& goal, int limit, float time_limit) {
    bool any_remaining = false, time_out;
    int solution_depth;
    shared_ptr<vector<shared_ptr<Node>>> successors;
    stack<shared_ptr<Node>> frontier = stack<shared_ptr<Node>>();
    unordered_map<pair<int, int>, bool, pair_hash> explored = unordered_map<pair<int, int>, bool, pair_hash>();
    shared_ptr<Node> current_node, node;
    frontier.push(root);
    while (!frontier.empty()) {
        setCurrentTime(clock());
        time_out = diff_clock(getCurrentTime(), getStartTime()) >= time_limit;
        current_node = frontier.top();
        frontier.pop();
        if (*current_node == *goal || time_out) {
            if (time_out) {
                return {nullptr, false};
            }
            setEndStatus(true);
            solution_depth = current_node->getDepth();
            setDN(double(solution_depth) / getExplored());
            calcEbf(solution_depth);
            return {current_node, true};
        }
        if (current_node->getDepth() < limit) {
            //increase the explored counter by one for the current node been expanded
            getInstance()._explored++;
            explored[pair<int, int>(current_node->getRow(), current_node->getCol())] = true;
            successors = current_node->successors(array,dimension);
            // because we use stack , we push successors in reversed order.
            for(auto i = successors->rbegin(); i!=successors->rend();++i ){
                    if(explored[pair<int,int>((*i)->getRow(),(*i)->getCol())]){
                        continue;
                    }
                    frontier.push((*i));
            }
        } else {
            any_remaining = true;
            update_cutoffs(current_node->getDepth());
        }

    }

    return {nullptr, any_remaining};
}

int IterativeDeepeningSearch::run_algorithm(int **array, int dimension, int *source, int *goal, float time_limit) {
    shared_ptr<Node> found;
    bool any_remaining;
    shared_ptr<Node> root ( new Node(0, 0, source[0], source[1], 0));
    root->insertElementToPath(pair<int, int>(source[0], source[1]));
    shared_ptr<Node>target( new Node(0, 0, goal[0], goal[1], 0));
    int max_depth = std::numeric_limits<int>::max();
    for (int i = 0; i < max_depth; ++i) {
        auto result = DLS(array, dimension, root, target, i, time_limit);
        found = result.first;
        any_remaining = result.second;
        if (found != nullptr) {
            // found the node
//            print_path(array, dimension, *found);
            generate_stats(*found);
            return 0;// return success.
        }
        if (diff_clock(getCurrentTime(), getStartTime()) >= time_limit || !any_remaining){
            //time out or end of search tree reached
            generate_stats(*root);
            break;

        }
    }
    return 1;//return failure.
}




