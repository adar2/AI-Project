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


pair<Node *, bool> IterativeDeepeningSearch::DLS(int **array, int dimension, Node *root, Node *goal, int limit) {
    bool any_remaining = false;
    int expand_counter, solution_depth;
    stack<Node *> frontier = stack<Node *>();
    unordered_map<pair<int, int>, bool, pair_hash> explored = unordered_map<pair<int, int>, bool, pair_hash>();
    Node *current_node, *node;
    int row = 0, col = 0;
    frontier.push(root);
    while (!frontier.empty()) {
        current_node = frontier.top();
        frontier.pop();
        if (*current_node == *goal) {
            // free all the nodes in the stack
            while (!frontier.empty()) {
                node = frontier.top();
                frontier.pop();
                delete node;
            }
            // update the number of nodes explored in the call goal found.
            IterativeDeepeningSearch &instance = getInstance();
            solution_depth = current_node->getDepth();
            instance._explored += explored.size();
            instance.setDN(double(solution_depth) / instance.getExplored());
            instance.setEbf(pow(instance.getExplored(), pow(solution_depth, -1)));
            return {current_node, true};
        }
        if (current_node->getDepth() < limit) {
            expand_counter = 0;
            explored[pair<int, int>(current_node->getRow(), current_node->getCol())] = true;
            for (int i = ACTIONS_SIZE - 1; i >= 0; --i) {
                switch (actions(i)) {
                    case U:
                        row = current_node->getRow() - 1;
                        col = current_node->getCol();
                        break;
                    case RU:
                        row = current_node->getRow() - 1;
                        col = current_node->getCol() + 1;
                        break;
                    case R:
                        row = current_node->getRow();
                        col = current_node->getCol() + 1;
                        break;
                    case RD:
                        row = current_node->getRow() + 1;
                        col = current_node->getCol() + 1;
                        break;
                    case D:
                        row = current_node->getRow() + 1;
                        col = current_node->getCol();
                        break;
                    case LD:
                        row = current_node->getRow() + 1;
                        col = current_node->getCol() - 1;
                        break;
                    case L:
                        row = current_node->getRow();
                        col = current_node->getCol() - 1;
                        break;
                    case LU:
                        row = current_node->getRow() - 1;
                        col = current_node->getCol() - 1;
                        break;
                }
                if (row < 0 || row >= dimension || col < 0 || col >= dimension || array[row][col] < 0)
                    continue;
                if (explored[pair<int, int>(row, col)])
                    continue;
                int cost = array[row][col] + current_node->getActualCost();
                node = new Node(cost, cost, row, col, current_node->getDepth() + 1);
                node->setPathTilNow(current_node->getPathTilNow());
                node->insertElementToPath(pair<int, int>(row, col));
                frontier.push(node);
                expand_counter++;
            }
            if (!expand_counter) {
                // cut off occurrence
                if (getInstance().getMin() == 0 || getInstance().getMin() > current_node->getDepth())
                    getInstance().setMin(current_node->getDepth());
                if (getInstance().getMax() == 0 || getInstance().getMax() < current_node->getDepth())
                    getInstance().setMax(current_node->getDepth());
                getInstance().addCutoffToSum(current_node->getDepth());
            }
        } else {
            any_remaining = true;
        }

        // compare pointers by addresses.
        if (current_node != root) {
            delete current_node;
        }
    }
    // maintain the number of nodes explored between function calls
    getInstance()._explored += explored.size();
    return {nullptr, any_remaining};
}

int IterativeDeepeningSearch::run_algorithm(int **array, int dimension, int *source, int *goal, float time_limit) {
    Node *found;
    bool any_remaining;
    Node *root = new Node(0, 0, source[0], source[1], 0);
    root->insertElementToPath(pair<int, int>(source[0], source[1]));
    Node *target = new Node(0, 0, goal[0], goal[1], 0);
    int max_depth = pow(dimension, 2);
    for (int i = 0; i < max_depth; ++i) {
        auto result = DLS(array, dimension, root, target, i);
        found = result.first;
        any_remaining = result.second;
        if (found != nullptr) {
            // found the node
            delete root;
            delete target;
            generate_stats(*found);
            delete found;
            return 0;// return success.
        }
        if (!any_remaining) {
            cout << "no nodes left" << endl;
            return 1;
        }
    }
    delete root;
    delete target;
    return 1;//return failure.
}



