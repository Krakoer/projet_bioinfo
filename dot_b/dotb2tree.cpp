#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <bits/stdc++.h>

using namespace std;

typedef struct Node
{
    int key1;
    int key2;
    vector<Node *> children;
    Node *parent;
} Node;

Node *newNode(int key1, int key2, Node *parent)
{
    Node *temp = new Node;
    temp->key1 = key1;
    temp->key2 = key2;
    temp->parent = parent;
    return temp;
}

Node *from_db(ifstream *input_file)
{
    Node *root = newNode(0, 0, nullptr);
    char name[4];

    //Read the file
    string line;
    int line_nb = 0;
    while (line_nb < 3)
    {
        getline(*input_file, line);
        line_nb++;
    }

    //std::cout << line;
    Node *cur_node = root;
    for (int i = 0; i < line.length(); i++)
    {
        char car = line[i];
        std::cout << "i = " << i << "   "
                  << "char = " << car << '\n';

        if (car == '.')
        {
            cur_node->children.push_back(newNode(i + 1, i + 1, cur_node));
        }
        else if (car == '(')
        {
            Node *child = newNode(i + 1, 0, cur_node);
            cur_node->children.push_back(child);
            cur_node = child;
        }
        else if (car == ')')
        {
            cur_node->key2 = i + 1;
            cur_node = cur_node->parent;
        }
        else
        {
            std::cout << "Wrong character encountered : " << car << '\n';
            return nullptr;
        }
    }

    return root;
}

void export_tree(ostream *output_file, Node *tree)
{
     if (tree == NULL)
        return;

    queue<Node *> q;
    q.push(tree);
    int i = 0;
    while (!q.empty())
    {
        int n = q.size();

        while (n > 0)
        {
            Node *p = q.front();
            q.pop();
            
            
            for (int j = 0; j < p->children.size(); j++)
                q.push(p->children[j]);
            n--;
        }

        cout << endl;
    }
}

// from https://www.geeksforgeeks.org/generic-tree-level-order-traversal/
void print_tree(Node *root)
{
    if (root == NULL)
        return;

    queue<Node *> q;
    q.push(root);
    while (!q.empty())
    {
        int n = q.size();

        while (n > 0)
        {
            Node *p = q.front();
            q.pop();
            cout << '[' << p->key1 << ';' << p->key2 << "] ";

            for (int i = 0; i < p->children.size(); i++)
                q.push(p->children[i]);
            n--;
        }

        cout << endl;
    }
}

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        std::cout << "Usage : ./dot2tree input_file";
        return 1;
    }

    char *input_path, *output_path;
    ifstream input_file;
    ofstream output_file;

    input_path = argv[1];
    output_path = argv[2];

    input_file.open(input_path, ios::in);
    output_file.open(output_path, ios::out);

    Node *tree = from_db(&input_file);
    if (tree != nullptr)
    {
        print_tree(tree);
        export_tree(&output_file, tree);

        std::cout << "File treated correcly\n";

        return 0;
    }

    else
    {
        std::cout << "An error occured\n";
        return 1;
    }
}