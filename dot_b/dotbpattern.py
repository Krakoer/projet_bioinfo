from anytree import Node, RenderTree
import sys


def tree_from_string(s):
    tree = Node("root")
    cur_node = tree
    nodes = []
    unpaired_cars = ['.', '[', ']', '{', '}', '<', '>']
    i = 0
    for car in s:
        if car != '&' and car != '*':
            if car in unpaired_cars:
                nodes.append(Node(i, parent=cur_node, paired=False))
            elif car == '(':
                child = Node(i, parent=cur_node, paired=True)
                cur_node = child
            elif car == ')':
                cur_node = cur_node.parent
            else:
                print(f"Wring charater encounterd : {car}")
        i += 1
    return tree


def load_from_file(path):
    trees = []
    f = open(path, 'r')
    ll = f.readline()
    f.close()
    for l in ll:
        trees.append(tree_from_string(l))
    return trees


if __name__ == "__main__":
    # if(len(sys.argv) < 4):
    #     print("Usage : python dotbpattern.py input_tree input_patterns output_file")
    input_path = sys.argv[1]
    trees = load_from_file(input_path)

    for pre, fill, node in RenderTree(trees[0]):
        print("%s%s" % (pre, node.name))
