import json
import sys

def main():
    if len(sys.argv) < 3:
        #print(f"Usage : python {sys.argv[0]} json_file output_path")
        return

    with open(sys.argv[1]) as f:
        data = json.load(f)

    print(data['num_isoCanonPairs'])

if __name__ == '__main__':
    main()