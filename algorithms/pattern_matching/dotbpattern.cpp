#include <bits/stdc++.h>

#include <algorithm>
#include <boost/filesystem.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

typedef struct Node {
    bool paired;      // If the node is a base paired or unapired base
    bool pseudoknot;  // If the node is involved in a pseudoknot
    bool openBP;      // If it is an open base, represented as a (*)
    int pos;          // Position in the string
    int posPaire;     // Position of the paire, if it is a base paired
    vector<Node *> children;
    Node *parent;
} Node;

Node *newNode(bool paired, Node *parent, int pos, bool pseudo = false, bool openBP = false) {
    Node *temp = new Node;
    temp->paired = paired;
    temp->parent = parent;
    temp->pseudoknot = pseudo;
    temp->openBP = openBP;
    temp->pos = pos;
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

vector<array<int, 2>> find_pattern(Node *tree, Node *pattern) {
    /*
    Count the number of occurences of the given pattern in the given tree and the positions
    */

    queue<Node *> q;  // We use a queue system to go through the tree again
    vector<array<int, 2>> pos;

    q.push(tree);

    while (!q.empty()) {
        int n = q.size();

        while (n > 0) {
            Node *p = q.front();
            q.pop();
            if (p->paired && check_from_node_rec(p, pattern)) {  // We naivly perform a matching test with each base paired node being the root
                array<int, 2> t = {p->pos, p->posPaire};
                pos.push_back(t);
            }

            for (int i = 0; i < p->children.size(); i++) q.push(p->children[i]);
            n--;
        }
    }

    return pos;
}

void print_output(char *name, vector<vector<int>> pos) {
    ofstream output_file;
    char *filename;
    strcpy(filename, name);
    strcat(filename, ".out");

    output_file.open(filename, ios::out);

    if (output_file.is_open()) {
        cout << "Printing output\n";

        for (int i = 0; i < pos.size(); i++) {
            // Go through the motifs
            vector<int> motif = pos[i];
            if (motif.size() > 0) {
                output_file << to_string(i) << '\t';
                for (int j = 0; j < motif.size(); j++) {
                    output_file << motif[j] << ';';
                }
                //output_file << '\t' << chaines.size() + 1 << '\n';
                output_file << '\n';
            }
        }
        output_file.close();
    }
}

vector<string> sep_chains(string line) {
    /*
    Breaks paired base which are not on the same chain (sparated by '&')
    */

    vector<tuple<int, int>> chains;

    int curChaine = 0, pos, chain;
    for (int i = 0; i < line.length(); i++) {
        char car = line[i];
        if (car == '(') {
            chains.push_back({i, curChaine});
        } else if (car == '&') {
            curChaine++;
        } else if (car == ')') {
            tie(pos, chain) = chains.back();
            chains.pop_back();
            if (chain != curChaine) {
                line[i] = '.';
                line[pos] = '.';
            }
        }
    }

    size_t p = 0;
    vector<string> res;
    while ((p = line.find('&')) != string::npos) {
        res.push_back(line.substr(0, p));
        line.erase(0, p + 1);
    }
    res.push_back(line);

    return res;
}

vector<Node *> tree_from_string(string line, bool motif = false) {
    /*
    Given a dot-bracket string, the function returns the equivalents trees as a netsed Node structure
    */
    vector<Node *> res;
    vector<string> lines = sep_chains(line);
    for (int l = 0; l < lines.size(); l++) {
        Node *root = newNode(true, nullptr, -1);
        Node *cur_node = root;

        string pseudoknot_chars = "[]{}<>";
        string str = lines[l];

        for (int i = 0; i < str.length(); i++) {  // Go through the string
            char car = str[i];
            if (car != '&') {                                                   // & are ignored for now
                if (car == '.') {                                               // If it is an unpaired base
                    cur_node->children.push_back(newNode(false, cur_node, i));  // Push an unpaired base Node
                } else if (car == '(') {                                        // If it is a new base paired create a new Node and push it
                    Node *child = newNode(true, cur_node, i);                   //
                    cur_node->children.push_back(child);                        //
                    cur_node = child;                                           //
                } else if (car == ')') {                                        //
                    cur_node->posPaire = i;                                     // Finaly, if it is a base paired already created we go up in the tree
                    cur_node = cur_node->parent;

                } else if (car == '*') {
                    cur_node->openBP = true;
                } else if (pseudoknot_chars.find(car) != string::npos) {
                    cur_node->children.push_back(newNode(false, cur_node, i, true));
                } else {
                    std::cout << "Wrong character encountered : " << car << '\n';
                }
            }
        }
        if (motif)
            res.push_back(root->children[0]);
        else {
            res.push_back(root);
        }
    }
    return res;
}

vector<vector<Node *>> load_from_file(ifstream *input_file, int line_offset = 0, bool motif = false, bool xdbn = false) {
    /*
    Loads a set of tree from a given file. Each line is considered as a tree, and the optional
    line_offset parameter allows to skip lines at the begining of the file (header for example)
    */
    vector<vector<Node *>> trees;
    string line;
    for (int i = 0; i < line_offset; i++) {
        getline(*input_file, line);
    }

    while (getline(*input_file, line)) {
        trees.push_back(tree_from_string(line, motif));
        if (xdbn)
            getline(*input_file, line);
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
            if (p->pseudoknot) {
                cout << '[';
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
        std::cout << "Usage : ./dotbpattern input_tree input_patterns output_path [options]\n./dotbpattern --help for help\n";

        if (argc >= 2 && argv[1] == (string) "--help") {
            cout << "Input tree file must be a single line file with a secondary structure in dot bracket notation" << endl
                 << "Input pattern file must contains one line per pattern, in a dot-brackets notation." << endl
                 << "Use --xdbn option to load xdbn instead of bdn" << endl;
        }
        return 1;
    }

    bool xdbn = false;
    for (int i = 0; i < argc; i++) {
        string arg = argv[i];
        if (arg.find("--xdbn") != string::npos) {
            xdbn = true;
        }
    }

    string input_path, output_path, pattern_path;
    ifstream input_file, pattern_file;
    ofstream output_file;

    input_path = argv[1];
    pattern_path = argv[2];
    output_path = argv[3];

    vector<string> input_files;

    pattern_file.open(pattern_path, ios::in);
    vector<vector<Node *>> patterns = load_from_file(&pattern_file, 0, true);
    pattern_file.close();

    boost::filesystem::path directory(argv[1]);
    if (boost::filesystem::exists(directory) && boost::filesystem::is_directory(directory)) {
        boost::filesystem::directory_iterator begin(directory);
        boost::filesystem::directory_iterator end;

        while (begin != end) {
            //cout << *begin << " ";
            input_files.push_back(begin->path().string());
            ++begin;
        }
        std::cout << "\n";
    }

    for (int f = 0; f < input_files.size(); f++) {
        cout << "Processing " << input_files[f] << endl;
        // Load the files
        input_file.open(input_files[f], ios::in);

        string name = input_files[f].substr(input_files[f].size() - 8, 4);
        vector<vector<Node *>> trees = load_from_file(&input_file, 2, false, xdbn);

        input_file.close();

        output_file.open(output_path + "/" + name + ".out", ios::out);

        for (int m = 0; m < trees.size(); m++) {  // For each model
            vector<Node *> chains = trees[m];
            for (int c = 0; c < chains.size(); c++) {  // And for each model
                Node *chain = chains[c];
                //cout << "Modèle N° " << m << "; Chaine n° " << c << endl;
                output_file << name << "-" << c << "-" << m << '\t';

                for (int p = 0; p < patterns.size(); p++) {
                    vector<array<int, 2>> positions = find_pattern(chain, patterns[p][0]);  // Considering pattenrs have 1 chaine only
                    if (positions.size() > 0) {
                        // ex : 99:12-19;58-75;
                        output_file << p << ':';
                        cout << "Pattern n° " << p << " found at pos ";
                        for (int pos = 0; pos < positions.size(); pos++) {
                            output_file << positions[pos][0] << '-' << positions[pos][1] << ';';
                            cout << positions[pos][0] << " ";
                        }
                        output_file << "\t";
                        cout << endl;
                    }
                }

                output_file << endl;
            }
        }

        output_file.close();
    }

    return 0;
}