import json
from tqdm import tqdm
import sys
from pathlib import Path


def fix_par(s):
        pstack = []
        fixed = [c for c in s]

        for i, c in enumerate(s):
            if c == '(':
                pstack.append(i)
            elif c == ')':
                if len(pstack) == 0:
                    fixed[i] = '.'
                else:
                    pstack.pop()

        while len(pstack) > 0:
            fixed[pstack.pop()] = '.'

        return ''.join(fixed)

def main():
    if len(sys.argv) < 3:
        #print(f"Usage : python {sys.argv[0]} json_file output_path")
        return

    json_pathlist = [str(path) for path in Path(sys.argv[1]).glob('**/*.json')]
    log_file = open('./log.txt', 'w')

    for json_path in tqdm(json_pathlist):
        name = json_path.split('/')[-1].split('.')[0]
        #print(name)
        with open(json_path) as f:
            try:
                data = json.load(f)
            except:
                log_file.write(name+"-JSON"'\n')
                continue

        f.close()
        out_file = open(sys.argv[2]+"/"+name+".xdbn", "w")

        #chains  = data['chains']
        if 'isoCanonPairs' in data:
            for isoPair in data['isoCanonPairs']:
                pos1 = isoPair['nt1'][3:]
                pos2 = isoPair['nt2'][3:]
                out_file.write(pos1+'-'+pos2+';')

        else:
            out_file.write(';')

        out_file.write('\n')

        if('dbn' in data):
            chains = list(data['dbn'].values())[1:]

            for chain_obj in chains:
                seq = chain_obj['bseq']
                db_string = chain_obj['sstr']

                out_file.write(seq)
                out_file.write('\n')
                out_file.write(fix_par(db_string))
                out_file.write('\n')

        else:
            log_file.write(name+'\n')

        out_file.close()
    log_file.close()


if __name__ == '__main__':
    main()
