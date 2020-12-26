import json
from tqdm import tqdm
import sys
from pathlib import Path
import re

canonPairs = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}


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


def isCanon(link):
    assert('bp' in link and 'LW' in link)
    bp1 = link['bp'][0].capitalize()
    bp2 = link['bp'][2].capitalize()
    try:
        answer = canonPairs[bp1] == bp2 and link["LW"] == 'cWW'
        return answer, True
    except:
        return False, False


def main():
    if len(sys.argv) < 3:
        #print(f"Usage : python {sys.argv[0]} json_files output_path")
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

        #out_file.write("HEAD_WILL_GO_HERE\n")

        try:
            pairs= data["pairs"]
        except:
            log_file.write(name+"-Pairs not in keys"'\n')
            continue

        nonCanonPairs = []
        header = ""
        for i, pair in enumerate(pairs):
            canon, success = isCanon(pair)
            if(not success):
                log_file.write(f"{name}-Pair n°{i}"'\n')
                continue
            if not canon:
                pos1=re.findall('[0-9]+', pair['nt1'])[-1]
                pos2=re.findall('[0-9]+', pair['nt2'])[-1]
                nonCanonPairs.append([int(pos1), int(pos2)])
                header += f"{pos1}-{pos2};"

        out_file.write(header+'\n')
                

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


def test():
    with open('../../RNA_files/json/5MQ0.json') as f:
        try:
            data = json.load(f)
        except:
            print("JSOn pb")
    pairs= data["pairs"]
    nonCanonPairs = []

    for i, pair in enumerate(pairs):
        canon, success = isCanon(pair)
        if(not success):
            #log_file.write(f"{name}-Pair n°{i}"'\n')
            print("Errooooor")
            continue
        if not canon:
            print(f"Non canon pair found index {i+1}")
            print(pair)
            nonCanonPairs.append([re.findall('[0-9]+', pair['nt1'])[-1], re.findall('[0-9]+', pair['nt2'])[-1]])
    print(nonCanonPairs)
    


if __name__ == '__main__':
    # main()
    test()