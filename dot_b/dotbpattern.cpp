#include <bits/stdc++.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

typedef struct Node {
    bool paired;  // If the node is a base paired or unapired base
    vector<Node *> children;
    Node *parent;
} Node;

Node *newNode(bool paired, Node *parent) {
    Node *temp = new Node;
    temp->paired = paired;
    temp->parent = parent;
    return temp;
}

bool check_from_root(Node *root, Node *pattern) {
    /*
    Takes a tree and a pattern, and checks if the pattern matches from the given root
    */
    queue<Node *> r, p;  // We use a queue system to go through the tree
    r.push(root);
    p.push(pattern);
    while (!r.empty() && !p.empty()) {
        int n = r.size();
        int m = p.size();

        while (n > 0 && m > 0) {
            Node *n1 = r.front();
            r.pop();
            Node *n2 = p.front();
            p.pop();

            if (n1->paired != n2->paired) {  // If the nodes are different, the pattern dosent match
                return false;
            }

            if (n2->children.size() == 0) {  // If the pattern has a base paired leaf, transform the equivalent base paired node of the tree to a leaf
                n1->children.clear();
            }

            for (int i = 0; i < n1->children.size(); i++) r.push(n1->children[i]);  // Push children of the nodes to the queue
            for (int i = 0; i < n2->children.size(); i++) p.push(n2->children[i]);
            n--;
            m--;
        }

        if (n != m) {
            // if n and m are not both 0, graphs dont match
            return false;
        }
    }
    // If p is not empty, it means that at a certain level, the pattern had less children than the given tree, so the graphs dont match
    return p.empty();
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
            if (p->paired && check_from_root(p, pattern)) {  // We naivly perform a matching test with each base paired node being the root
                count++;
                //cout << "Pattern matches !!" << count << endl;
            }

            for (int i = 0; i < p->children.size(); i++) q.push(p->children[i]);
            n--;
        }
    }

    return count;
}

Node *tree_from_string(string line) {
    /*
    Given a dot-bracket string, the function returns the equivalent tree as a netsed Node structure
    For now, only . and () are considered : & are ignored, and [, ] , { and  } are considered as .
    */
    Node *root = newNode(true, nullptr);
    Node *cur_node = root;

    char unpaired_chars[] = {'.', '[', ']', '{', '}'};

    for (int i = 0; i < line.length(); i++) {  // Go through the string
        char car = line[i];
        if (car != '&') {                                                                        // & are ignored for now
            if (find(begin(unpaired_chars), end(unpaired_chars), car) != end(unpaired_chars)) {  // If it is an unpaired base
                cur_node->children.push_back(newNode(false, cur_node));                          // Push an unpaired base Node
            } else if (car == '(') {                                                             // If it is a new base paired create a new Node and push it
                Node *child = newNode(true, cur_node);                                           //
                cur_node->children.push_back(child);                                             //
                cur_node = child;                                                                //
            } else if (car == ')') {                                                             // Finaly, if it is a base paired already created we go up in the tree
                cur_node = cur_node->parent;
            } else {
                std::cout << "Wrong character encountered : " << car << '\n';
                return nullptr;
            }
        }
    }

    return root;
}

vector<Node *> load_from_file(ifstream *input_file, int line_offset = 0) {
    /*
    Loads a set of tree from a given file. Each line is considered as a tree, and the optional
    line_offset parameter allows to skip lines at the begining of the file (header for example)
    */
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
    if (root == NULL) return;

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

    Node *tree = load_from_file(&input_file).back();
    vector<Node *> patterns = load_from_file(&pattern_file);

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