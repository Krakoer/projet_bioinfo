import os
import sys
from pathlib import Path

def main():
    if len(sys.argv) != 2:
        print(f"Usage : {sys.argv[0]} path_to_pdb_files")
        exit()

    path_pdb = sys.argv[1]
    pdb_files = [str(path) for path in Path(path_pdb).glob('**/*.pdb')]
    if not os.path.exists("./json_files"):
        os.popen("mkdir json_files").close()
    already_treated = [str(path).split('/')[-1].split('.')[0] for path in Path('./json_files').glob('**/*.json')]

    for f in pdb_files:
        name = f.split('/')[-1].split('.')[0]
        if name not in already_treated:
            cmd = f"./x3dna-dssr -i={f} -o=./json_files/{name}.json --json"
            os.popen(cmd).close()
            os.popen("rm dssr*").close()

if __name__ == '__main__':
    main()