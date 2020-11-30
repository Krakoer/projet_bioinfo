#include <bits/stdc++.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

typedef struct Node {
    bool paired;      // If the node is a base paired or unapired base
    bool pseudoknot;  // If the node is involved in a pseudoknot
    bool openBP;
    vector<Node *> children;
    Node *parent;
} Node;

Node *newNode(bool paired, Node *parent, bool pseudo = false, bool openBP = false) {
    Node *temp = new Node;
    temp->paired = paired;
    temp->parent = parent;
    temp->pseudoknot = pseudo;
    temp->openBP = openBP;
    return temp;
}

bool check_from_node_rec(Node *node, Node *motif) {
    if (node->parent == nullptr) {
        return false;
    } else if (!node->paired && !motif->paired) {
        return true;
    } else if (motif->openBP && node->openBP) {
        return true;
    } else if (motif->openBP && node->paired) {
        return true;
    } else if (motif->paired && node->paired) {
        if (motif->children.size() == node->children.size()) {
            for (int i = 0; i < node->children.size(); i++) {
                if (!check_from_node_rec(node->children[i], motif->children[i])) {
                    return false;
                }
            }
            return true;
        }
    }
    return false;
}

int count_pattern(Node *tree, Node *pattern) {
    /*
    Count the number of occurences of the given pattern in the given tree
    */

    queue<Node *> q;  // We use a queue system to go through the tree again

    q.push(tree);
    int count = 0;
    while (!q.empty()) {
        int n = q.size();

        while (n > 0) {
            Node *p = q.front();
            q.pop();
            if (p->paired && check_from_node_rec(p, pattern)) {  // We naivly perform a matching test with each base paired node being the root
                count++;
                //cout << "Pattern matches !!" << count << endl;
            }

            for (int i = 0; i < p->children.size(); i++) q.push(p->children[i]);
            n--;
        }
    }

    return count;
}

// https://stackoverflow.com/questions/5878775/how-to-find-and-replace-string
std::string ReplaceString(std::string subject, const std::string &search,
                          const std::string &replace) {
    size_t pos = 0;
    while ((pos = subject.find(search, pos)) != std::string::npos) {
        subject.replace(pos, search.length(), replace);
        pos += replace.length();
    }
    return subject;
}

Node *tree_from_string(string line, bool motif = false) {
    /*
    Given a dot-bracket string, the function returns the equivalent tree as a netsed Node structure
    */
    Node *root = newNode(true, nullptr);
    Node *cur_node = root;

    string pseudoknot_chars = "[]{}<>";

    for (int i = 0; i < line.length(); i++) {  // Go through the string
        char car = line[i];
        if (car != '&') {                                                                                            // & are ignored for now
            if (car == '.') {                                                                                        // If it is an unpaired base
                cur_node->children.push_back(newNode(false, cur_node, pseudoknot_chars.find(car) != string::npos));  // Push an unpaired base Node, and check if it is involved in a pseudoknot
            } else if (car == '(') {                                                                                 // If it is a new base paired create a new Node and push it
                Node *child = newNode(true, cur_node);                                                               //
                cur_node->children.push_back(child);                                                                 //
                cur_node = child;                                                                                    //
            } else if (car == ')') {                                                                                 // Finaly, if it is a base paired already created we go up in the tree
                cur_node = cur_node->parent;
            } else if (car == '*') {
                cur_node->openBP = true;
            }
        } else {
            std::cout << "Wrong character encountered : " << car << '\n';
            return nullptr;
        }
    }
    if (motif)
        return root->children[0];
    return root;
}

vector<int> findChaines(string line) {
    // Return posistions of &
    vector<int> pos;

    for (int i = 0; i < line.length(); i++) {
        char car = line[i];
        if (car == '&') {
            pos.push_back(i);
        }
    }

    return pos;
}

vector<Node *> load_from_file(ifstream *input_file, int line_offset) {
    /*
    Loads a set of tree from a given file. Each line is considered as a tree, and the optional
    line_offset parameter allows to skip lines at the begining of the file (header for example)
    */
    cout << line_offset << endl;
    vector<Node *> trees;
    string line;
    for (int i = 0; i < line_offset; i++) {
        getline(*input_file, line);
    }

    while (getline(*input_file, line)) {
        trees.push_back(tree_from_string(line));
    }

    return trees;
}

// from https://www.geeksforgeeks.org/generic-tree-level-order-traversal/
void print_tree(Node *root) {
    /*
    Print each layer of the tree : [-] is a base paired and [.] is unpaired
    */
    if (root == nullptr) return;

    queue<Node *> q;
    q.push(root);
    while (!q.empty()) {
        int n = q.size();

        while (n > 0) {
            Node *p = q.front();
            q.pop();
            cout << '[';
            if (p->paired) {
                cout << '-';
            } else {
                cout << '.';
            }
            if (p->openBP) {
                cout << '*';
            }
            cout << "] ";

            for (int i = 0; i < p->children.size(); i++) q.push(p->children[i]);
            n--;
        }

        cout << endl;
    }
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        std::cout << "Usage : ./dotbpattern input_tree input_patterns output_path\n./dotbpattern --help for help\n";

        if (argc >= 2 && argv[1] == (string) "--help") {
            cout << "Input tree file must be a single line file with a secondary structure in dot bracket notation" << endl
                 << "Input pattern file must contains one line per pattern, in a dot-brackets notation." << endl;
        }
        return 1;
    }

    char *input_path, *output_path, *pattern_path;
    ifstream input_file, pattern_file;
    ofstream output_file;

    input_path = argv[1];
    pattern_path = argv[2];
    output_path = argv[3];

    // Load the files
    input_file.open(input_path, ios::in);
    pattern_file.open(pattern_path, ios::in);
    output_file.open(output_path, ios::out);

    Node *tree = load_from_file(&input_file, 2).back();
    vector<Node *> patterns = load_from_file(&pattern_file, 0);

    // Close the files
    input_file.close();
    pattern_file.close();

    // We just loop through the patterns and count the number of occurence of each pattern
    // The we write this number in the output file
    if (tree != nullptr && patterns.size() > 0) {
        cout << "File loaded correcly\n";
        cout << "Loaded " << patterns.size() << " patterns\n";
        for (int i = 0; i < patterns.size(); i++) {
            int count = count_pattern(tree, patterns[i]);
            output_file << count << endl;
            cout << "Pattern " << i << " found " << count << " times.\n";
        }
        output_file.close();
        return 0;
    }

    else {
        output_file.close();
        cout << "An error occured\n";
        return 1;
    }
}